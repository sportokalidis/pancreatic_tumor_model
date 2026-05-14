// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey
// BioDynaMo collaboration. Apache-2.0 license.
//
// -----------------------------------------------------------------------------
#ifndef PANCREATIC_TUMOR_MODEL_H_
#define PANCREATIC_TUMOR_MODEL_H_

#include "biodynamo.h"
#include "core/environment/uniform_grid_environment.h"
#include "params/sim_param.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <mutex>

namespace bdm {
namespace pancreatic_tumor {

// Convenience accessor for SimParam — requires full BioDynaMo environment.
// Defined here (not in sim_param.h) to avoid pulling Param's full definition
// into the param-group header.
inline const SimParam* SP() {
  return Simulation::GetActive()->GetParam()->Get<SimParam>();
}

// ============================================================================
// Math helpers
// ============================================================================
inline real_t Clamp(real_t v, real_t lo, real_t hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}
inline Real3 ClampPoint(const Real3& pos, real_t lo, real_t hi) {
  return {Clamp(pos[0], lo, hi), Clamp(pos[1], lo, hi), Clamp(pos[2], lo, hi)};
}
// Hill saturation: x / (K + x)
inline real_t Sat(real_t x, real_t K) {
  return (x <= 0.0) ? 0.0 : x / (K + x);
}
// Convert a per-day rate to a per-step probability via the exact Poisson mapping.
inline real_t ProbFromRate(real_t rate_per_day, real_t dt_day) {
  if (rate_per_day <= 0.0) return 0.0;
  return 1.0 - std::exp(-rate_per_day * dt_day);
}

// ============================================================================
// Agent types
// ============================================================================

// Each agent stores a color_ field and implements Initialize so that daughter
// cells inherit it correctly after Divide().  The color_ field is the only
// custom field here — it is enough for visualization.

class TumorCell : public Cell {
  BDM_AGENT_HEADER(TumorCell, Cell, 1);

 public:
  TumorCell() = default;
  explicit TumorCell(const Real3& p) {
    SetPosition(p);
    SetDiameter(2.0 * SP()->cell_radius_um);
    color_ = SP()->color_tumor;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    color_ = bdm_static_cast<TumorCell*>(event.existing_agent)->color_;
  }

  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }

 private:
  int color_ = 0;
};

class StellateCell : public Cell {
  BDM_AGENT_HEADER(StellateCell, Cell, 1);

 public:
  StellateCell() = default;
  explicit StellateCell(const Real3& p) {
    SetPosition(p);
    SetDiameter(2.0 * SP()->cell_radius_um);
    color_ = SP()->color_psc;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    color_ = bdm_static_cast<StellateCell*>(event.existing_agent)->color_;
  }

  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }

 private:
  int color_ = 0;
};

class EffectorTCell : public Cell {
  BDM_AGENT_HEADER(EffectorTCell, Cell, 1);

 public:
  EffectorTCell() = default;
  explicit EffectorTCell(const Real3& p) {
    SetPosition(p);
    SetDiameter(2.0 * SP()->cell_radius_um);
    color_ = SP()->color_eff;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    color_ = bdm_static_cast<EffectorTCell*>(event.existing_agent)->color_;
  }

  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }

 private:
  int color_ = 0;
};

class NKCell : public Cell {
  BDM_AGENT_HEADER(NKCell, Cell, 1);

 public:
  NKCell() = default;
  explicit NKCell(const Real3& p) {
    SetPosition(p);
    SetDiameter(2.0 * SP()->cell_radius_um);
    color_ = SP()->color_nk;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    color_ = bdm_static_cast<NKCell*>(event.existing_agent)->color_;
  }

  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }

 private:
  int color_ = 0;
};

class HelperTCell : public Cell {
  BDM_AGENT_HEADER(HelperTCell, Cell, 1);

 public:
  HelperTCell() = default;
  explicit HelperTCell(const Real3& p) {
    SetPosition(p);
    SetDiameter(2.0 * SP()->cell_radius_um);
    color_ = SP()->color_helper;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    color_ = bdm_static_cast<HelperTCell*>(event.existing_agent)->color_;
  }

  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }

 private:
  int color_ = 0;
};

class TRegCell : public Cell {
  BDM_AGENT_HEADER(TRegCell, Cell, 1);

 public:
  TRegCell() = default;
  explicit TRegCell(const Real3& p) {
    SetPosition(p);
    SetDiameter(2.0 * SP()->cell_radius_um);
    color_ = SP()->color_treg;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    color_ = bdm_static_cast<TRegCell*>(event.existing_agent)->color_;
  }

  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }

 private:
  int color_ = 0;
};

// ============================================================================
// Population census
// ============================================================================
struct Counts { size_t C=0, P=0, E=0, N=0, H=0, R=0; };

// Global census: computed at most once per step via double-checked locking.
// Thread-safe: multiple behaviors call RefreshIfNeeded() concurrently but only
// one executes the ForEachAgent scan, the rest read the cached result.
struct GlobalCensus {
  std::atomic<size_t> step_cached{std::numeric_limits<size_t>::max()};
  Counts              cnt;
  std::mutex          mtx;

  static GlobalCensus& Instance() { static GlobalCensus gc; return gc; }

  void RefreshIfNeeded() {
    auto* sim   = Simulation::GetActive();
    size_t step = sim->GetScheduler()->GetSimulatedSteps();

    // Fast path: already cached for this step.
    if (step_cached.load(std::memory_order_acquire) == step) return;

    std::lock_guard<std::mutex> lock(mtx);
    // Re-check inside lock (double-checked locking).
    if (step_cached.load(std::memory_order_relaxed) == step) return;

    Counts c;
    sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
      if      (dynamic_cast<TumorCell*>(a))     ++c.C;
      else if (dynamic_cast<StellateCell*>(a))  ++c.P;
      else if (dynamic_cast<EffectorTCell*>(a)) ++c.E;
      else if (dynamic_cast<NKCell*>(a))        ++c.N;
      else if (dynamic_cast<HelperTCell*>(a))   ++c.H;
      else if (dynamic_cast<TRegCell*>(a))      ++c.R;
    });
    cnt = c;
    step_cached.store(step, std::memory_order_release);
  }

  const Counts& Get() const { return cnt; }
};

// Local neighborhood census — counts within a radius around a given agent.
struct LocalNeighborhoodCounter {
  struct Fun : public Functor<void, Agent*, real_t> {
    Counts* out;
    explicit Fun(Counts* o) : out(o) {}
    void operator()(Agent* nb, real_t /*sq*/) override {
      if      (dynamic_cast<TumorCell*>(nb))     ++out->C;
      else if (dynamic_cast<StellateCell*>(nb))  ++out->P;
      else if (dynamic_cast<EffectorTCell*>(nb)) ++out->E;
      else if (dynamic_cast<NKCell*>(nb))        ++out->N;
      else if (dynamic_cast<HelperTCell*>(nb))   ++out->H;
      else if (dynamic_cast<TRegCell*>(nb))      ++out->R;
    }
  };

  static Counts Around(const Agent& agent, real_t radius_um) {
    Counts out;
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    Fun   f(&out);
    ctxt->ForEachNeighbor(f, agent, radius_um * radius_um);
    return out;
  }
};

// Unified accessor: returns global or local counts depending on params.
inline Counts GetCounts(const Agent& self) {
  const auto* sp = SP();
  if (!sp->use_local_counts) {
    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();
    return gc.Get();
  }
  return LocalNeighborhoodCounter::Around(self, sp->local_radius_um);
}

// Always returns global counts — used by SourceBehavior because immune cell
// trafficking from the periphery is a systemic process, not a local one.
inline Counts GetGlobalCounts() {
  auto& gc = GlobalCensus::Instance();
  gc.RefreshIfNeeded();
  return gc.Get();
}

// ============================================================================
// Drug state (Paper Section 5)
// Holds M_gem, M_abr (drug concentrations) and the Anti-CD47 active flag.
// Updated at most once per step via double-checked locking (same pattern as
// GlobalCensus) so behaviors reading it always see the step-consistent value.
// ============================================================================
struct DrugState {
  std::atomic<size_t> step_cached{std::numeric_limits<size_t>::max()};
  std::mutex mtx;

  real_t M_gem = 0.0;
  real_t M_abr = 0.0;
  bool   acd47_active = false;

  static DrugState& Instance() { static DrugState ds; return ds; }

  void RefreshIfNeeded() {
    auto* sim   = Simulation::GetActive();
    size_t step = sim->GetScheduler()->GetSimulatedSteps();
    if (step_cached.load(std::memory_order_acquire) == step) return;

    std::lock_guard<std::mutex> lock(mtx);
    if (step_cached.load(std::memory_order_relaxed) == step) return;

    const auto* sp = SP();
    // Short-circuit: no treatment enabled
    if (!sp->treat_gem && !sp->treat_abr && !sp->treat_acd47) {
      step_cached.store(step, std::memory_order_release);
      return;
    }

    const real_t dt_day    = sp->dt_minutes / 1440.0;
    const real_t day_off   = 7.0;  // paper day 7 = simulation day 0
    const real_t spd_f     = 1440.0 / sp->dt_minutes;  // steps per day

    // Convert paper-day to step index
    auto s_of = [&](real_t pday) -> size_t {
      real_t sim_day = pday - day_off;
      return (sim_day <= 0.0) ? 0 : static_cast<size_t>(std::round(sim_day * spd_f));
    };
    // True if this step falls on an injection grid point
    auto is_inject = [&](size_t s0, size_t s1, size_t freq) -> bool {
      return step >= s0 && step <= s1 && ((step - s0) % freq == 0);
    };

    // Gemcitabine: exponential decay then pulse injection
    if (sp->treat_gem) {
      M_gem *= std::exp(-sp->gem_gamma * dt_day);
      size_t s0   = s_of(sp->treat_start_day);
      size_t s1   = s_of(sp->gem_end_day);
      size_t freq = static_cast<size_t>(std::round(sp->gem_freq_days * spd_f));
      if (freq > 0 && is_inject(s0, s1, freq)) M_gem += sp->gem_dose;
    }

    // Abraxane
    if (sp->treat_abr) {
      M_abr *= std::exp(-sp->abr_gamma * dt_day);
      size_t s0   = s_of(sp->treat_start_day);
      size_t s1   = s_of(sp->abr_end_day);
      size_t freq = static_cast<size_t>(std::round(sp->abr_freq_days * spd_f));
      if (freq > 0 && is_inject(s0, s1, freq)) M_abr += sp->abr_dose;
    }

    // Anti-CD47: binary active flag (no PK curve)
    acd47_active = sp->treat_acd47
                && step >= s_of(sp->treat_start_day)
                && step <= s_of(sp->acd47_end_day);

    step_cached.store(step, std::memory_order_release);
  }
};

// Volume ratio V_total / V_local.
// In global mode = 1.0 (no adjustment).
// In local mode, multiplying raw local counts by this factor gives their
// "global equivalent", so every rate formula stays numerically identical
// at uniform cell density while still capturing spatial heterogeneity.
inline real_t DensityCompensation() {
  const auto* sp = SP();
  if (!sp->use_local_counts) return 1.0;
  real_t side    = sp->max_bound - sp->min_bound;
  real_t V_total = side * side * side;
  real_t r       = sp->local_radius_um;
  real_t V_local = (4.0 / 3.0) * M_PI * r * r * r;
  return V_total / V_local;
}

// Density-compensated counts for use inside behaviors.
// In global mode dc=1: values equal raw integer counts (same as before).
// In local mode dc=V_total/V_local: local counts scaled to global equivalent.
struct EffCounts {
  real_t C = 0, P = 0, E = 0, N = 0, H = 0, R = 0;
  EffCounts(const Counts& c, real_t dc)
      : C(static_cast<real_t>(c.C) * dc),
        P(static_cast<real_t>(c.P) * dc),
        E(static_cast<real_t>(c.E) * dc),
        N(static_cast<real_t>(c.N) * dc),
        H(static_cast<real_t>(c.H) * dc),
        R(static_cast<real_t>(c.R) * dc) {}
};

// ============================================================================
// Population logger
// Owns the CSV file. Created once in Simulate(), stays alive for the run.
// Behaviors write via PopulationLogger::Instance().
// ============================================================================
struct PopulationLogger {
  std::ofstream csv;

  static PopulationLogger& Instance() {
    static PopulationLogger logger;
    return logger;
  }

  void Open(const std::string& dir) {
    std::filesystem::create_directories(dir);
    std::string path = dir + "/populations.csv";
    csv.open(path);
    csv << "step,days,C,P,E,N,H,R,total,M_gem,M_abr\n";
    csv.flush();
  }

  void Write(size_t step, real_t t_day,
             size_t C, size_t P, size_t E, size_t N, size_t H, size_t R,
             real_t M_gem = 0.0, real_t M_abr = 0.0) {
    csv << step << "," << t_day << ","
        << C << "," << P << "," << E << ","
        << N << "," << H << "," << R << ","
        << (C + P + E + N + H + R) << ","
        << M_gem << "," << M_abr << "\n";
    csv.flush();
  }
};

// ============================================================================
// Behaviors
//
// Design rules:
//   • AlwaysCopyToNew() in constructor  → behavior auto-copies to daughter
//   • No manual AddBehavior after Divide()
//   • No ClampPoint inside behaviors (BoundSpaceMode::kClosed handles walls)
//   • Death and division are independent events — both can fire in one step
//   • gate_C_K removed — not in reference ODE
// ============================================================================

class TumorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TumorBehavior, Behavior, 1);

 public:
  TumorBehavior() { AlwaysCopyToNew(); }

  void Run(Agent* a) override {
    auto* c    = bdm_static_cast<TumorCell*>(a);
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto* rng  = Simulation::GetActive()->GetRandom();
    const auto* sp = SP();

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const EffCounts cnt(GetCounts(*c), DensityCompensation());

    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();

    // --- Division (Eq. 2.1 growth terms) ---
    // (k_c + mu_c·P)·C·(1-C/K_C): mu_c·P is LINEAR in P (not Hill)
    // Crowding uses global C — K_C is a systemic resource/space constraint.
    // Local C would allow cells in Poisson-sparse pockets to divide past K_C
    // globally, causing 30-40% overshoot above the ODE steady state.
    const real_t C_global_crowd = sp->use_local_counts
        ? static_cast<real_t>(GetGlobalCounts().C) : cnt.C;
    real_t crowd  = 1.0 - Clamp(C_global_crowd / sp->K_C, 0.0, 1.0);
    real_t boostP = sp->c_boost_from_P * cnt.P;  // mu_c·P
    real_t div_rate = (sp->c_base_div + boostP) * crowd;

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      c->Divide();
      c->SetCellColor(sp->color_tumor_div);
    } else {
      c->SetCellColor(sp->color_tumor);
    }

    // --- Killing (Eq. 2.1 death terms) — independent of division ---
    // b_c·N·C (bilinear in N) and d_c·E·C/(1+r1·R) (bilinear in E)
    // Treg suppression of CTL killing is cytokine-mediated (TGF-β/IL-10) and
    // diffuses beyond the local radius — use global R to avoid singularity at
    // R_local=0, which causes certain death for ~9% of cells per step.
    const real_t R_global_c = sp->use_local_counts
        ? static_cast<real_t>(GetGlobalCounts().R) : cnt.R;
    real_t inhibit_R = 1.0 / (1.0 + sp->c_R_blocks_E * R_global_c);
    real_t killE = sp->c_kill_by_E * cnt.E * inhibit_R;
    real_t killN = sp->c_kill_by_N * cnt.N;

    // --- Drug kill (Eq. 5.1): c_c·(1−e^{−M})·C ---
    real_t drug_kill = sp->gem_c_c * (1.0 - std::exp(-ds.M_gem))
                     + sp->abr_c_c * (1.0 - std::exp(-ds.M_abr));

    if (rng->Uniform(0, 1) < ProbFromRate(killE + killN + drug_kill, dt_day)) {
      ctxt->RemoveAgent(c->GetUid());
    }
  }
};

class PSCBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(PSCBehavior, Behavior, 1);

 public:
  PSCBehavior() { AlwaysCopyToNew(); }

  void Run(Agent* a) override {
    auto* psc  = bdm_static_cast<StellateCell*>(a);
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto* rng  = Simulation::GetActive()->GetRandom();
    const auto* sp = SP();

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const EffCounts cnt(GetCounts(*psc), DensityCompensation());

    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();

    // --- Death (Eq. 5.2): lambda_p·P + c_p·(1−e^{−M})·P ---
    real_t drug_kill = sp->gem_c_p * (1.0 - std::exp(-ds.M_gem))
                     + sp->abr_c_p * (1.0 - std::exp(-ds.M_abr));
    if (rng->Uniform(0, 1) < ProbFromRate(sp->p_base_death + drug_kill, dt_day)) {
      ctxt->RemoveAgent(psc->GetUid());
      return;
    }

    // --- Division (Eq. 2.2) ---
    // (k_p + f_p*C/(mu_p+C))*P*(1-a_p*P)
    const real_t P_global_crowd = sp->use_local_counts
        ? static_cast<real_t>(GetGlobalCounts().P) : cnt.P;
    real_t crowd  = 1.0 - Clamp(P_global_crowd / sp->K_P, 0.0, 1.0);
    real_t boostC = sp->p_boost_from_C * Sat(cnt.C, sp->p_boost_from_C_K);
    real_t div_rate = (sp->p_base_div + boostC) * crowd;

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      psc->Divide();
    }
  }
};

class EffectorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(EffectorBehavior, Behavior, 1);

 public:
  EffectorBehavior() { AlwaysCopyToNew(); }

  void Run(Agent* a) override {
    auto* e    = bdm_static_cast<EffectorTCell*>(a);
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto* rng  = Simulation::GetActive()->GetRandom();
    const auto* sp = SP();

    // Random walk — enables immune infiltration of the tumor sphere
    if (sp->use_local_counts) {
      const Real3 pos = e->GetPosition();
      const real_t s  = sp->immune_step_um;
      e->SetPosition(ClampPoint({pos[0] + rng->Gaus(0.0, s),
                                 pos[1] + rng->Gaus(0.0, s),
                                 pos[2] + rng->Gaus(0.0, s)},
                                sp->min_bound, sp->max_bound));
    }

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const EffCounts cnt(GetCounts(*e), DensityCompensation());

    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();

    // --- Death (Eq. 5.3 death terms): b_e + c_e·C + δ_e·R + c_immune·(1−e^{−M}) ---
    const real_t R_global_e = sp->use_local_counts
        ? static_cast<real_t>(GetGlobalCounts().R) : cnt.R;
    real_t drug_kill = sp->gem_c_immune * (1.0 - std::exp(-ds.M_gem))
                     + sp->abr_c_immune * (1.0 - std::exp(-ds.M_abr));
    real_t die = sp->e_base_death
               + sp->e_inact_by_C * cnt.C
               + sp->e_suppr_by_R * R_global_e
               + drug_kill;

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(e->GetUid());
      return;
    }

    // --- Per-cell proliferation (Eq. 2.3 + Eq. 5.3 Anti-CD47 term) ---
    // p_e·H/(g_e+H) from the ODE; α·v_a boost when Anti-CD47 is active.
    real_t div_rate = sp->e_help_from_H * Sat(cnt.H, sp->e_help_from_H_K);
    if (ds.acd47_active) div_rate += sp->acd47_e_boost;

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      e->Divide();
    }
  }
};

class NKBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(NKBehavior, Behavior, 1);

 public:
  NKBehavior() { AlwaysCopyToNew(); }

  void Run(Agent* a) override {
    auto* n    = bdm_static_cast<NKCell*>(a);
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto* rng  = Simulation::GetActive()->GetRandom();
    const auto* sp = SP();

    if (sp->use_local_counts) {
      const Real3 pos = n->GetPosition();
      const real_t s  = sp->immune_step_um;
      n->SetPosition(ClampPoint({pos[0] + rng->Gaus(0.0, s),
                                 pos[1] + rng->Gaus(0.0, s),
                                 pos[2] + rng->Gaus(0.0, s)},
                                sp->min_bound, sp->max_bound));
    }

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const EffCounts cnt(GetCounts(*n), DensityCompensation());

    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();

    // --- Death (Eq. 5.4 death terms): b_n + c_n·C + δ_n·R + c_immune·(1−e^{−M}) ---
    const real_t R_global_n = sp->use_local_counts
        ? static_cast<real_t>(GetGlobalCounts().R) : cnt.R;
    real_t drug_kill = sp->gem_c_immune * (1.0 - std::exp(-ds.M_gem))
                     + sp->abr_c_immune * (1.0 - std::exp(-ds.M_abr));
    real_t die = sp->n_base_death
               + sp->n_inact_by_C * cnt.C
               + sp->n_suppr_by_R * R_global_n
               + drug_kill;

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(n->GetUid());
      return;
    }

    // --- Per-cell proliferation (Eq. 2.4, p_n·H·N/(g_n+H) term) ---
    // Constant influx a_n handled by SourceBehavior.
    real_t div_rate = sp->n_help_from_H * Sat(cnt.H, sp->n_help_from_H_K);

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      n->Divide();
    }
  }
};

class HelperBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(HelperBehavior, Behavior, 1);

 public:
  HelperBehavior() { AlwaysCopyToNew(); }

  void Run(Agent* a) override {
    auto* h    = bdm_static_cast<HelperTCell*>(a);
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto* rng  = Simulation::GetActive()->GetRandom();
    const auto* sp = SP();

    if (sp->use_local_counts) {
      const Real3 pos = h->GetPosition();
      const real_t s  = sp->immune_step_um;
      h->SetPosition(ClampPoint({pos[0] + rng->Gaus(0.0, s),
                                 pos[1] + rng->Gaus(0.0, s),
                                 pos[2] + rng->Gaus(0.0, s)},
                                sp->min_bound, sp->max_bound));
    }

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const EffCounts cnt(GetCounts(*h), DensityCompensation());

    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();

    // --- Death (Eq. 5.5 death terms): b_h + δ_h·R + c_immune·(1−e^{−M}) ---
    const real_t R_global_h = sp->use_local_counts
        ? static_cast<real_t>(GetGlobalCounts().R) : cnt.R;
    real_t drug_kill = sp->gem_c_immune * (1.0 - std::exp(-ds.M_gem))
                     + sp->abr_c_immune * (1.0 - std::exp(-ds.M_abr));
    real_t die = sp->h_base_death
               + sp->h_suppr_by_R * R_global_h
               + drug_kill;

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(h->GetUid());
      return;
    }

    // --- Per-cell self-activation (Eq. 2.5, p_h·H²/(g_h+H) = p_h·H/(g_h+H) per cell) ---
    // Constant influx a_h handled by SourceBehavior.
    real_t div_rate = sp->h_self_act * Sat(cnt.H, sp->h_self_act_K);

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      h->Divide();
    }
  }
};

class TRegBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TRegBehavior, Behavior, 1);

 public:
  TRegBehavior() { AlwaysCopyToNew(); }

  void Run(Agent* a) override {
    auto* r    = bdm_static_cast<TRegCell*>(a);
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto* rng  = Simulation::GetActive()->GetRandom();
    const auto* sp = SP();

    if (sp->use_local_counts) {
      const Real3 pos = r->GetPosition();
      const real_t s  = sp->immune_step_um;
      r->SetPosition(ClampPoint({pos[0] + rng->Gaus(0.0, s),
                                 pos[1] + rng->Gaus(0.0, s),
                                 pos[2] + rng->Gaus(0.0, s)},
                                sp->min_bound, sp->max_bound));
    }

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const EffCounts cnt(GetCounts(*r), DensityCompensation());

    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();

    // --- Death (Eq. 5.6 death terms): δ_r + r·N + c_immune·(1−e^{−M}) ---
    real_t drug_kill = sp->gem_c_immune * (1.0 - std::exp(-ds.M_gem))
                     + sp->abr_c_immune * (1.0 - std::exp(-ds.M_abr));
    real_t die = sp->r_decay
               + sp->r_cleared_by_N * cnt.N
               + drug_kill;

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(r->GetUid());
      return;
    }

    // --- Per-cell proliferation (Eq. 2.6, p_r·H·R/(g_r+H) term) ---
    // Absolute sources a, a_r·E, b_r·H handled by SourceBehavior.
    real_t div_rate = sp->r_prolif_by_H * Sat(cnt.H, sp->r_prolif_by_H_K);

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      r->Divide();
    }
  }
};

// ============================================================================
// SourceBehavior
// Implements density-independent (constant) immune recruitment terms from
// the ODE: a_e, a_n, a_h, a_r  (Eqs. 2.3–2.6).
//
// Each is a constant influx (cells/day) independent of current population.
// Attached to the reporter agent; runs every step.  Cells are added via
// ctxt->AddAgent so they appear in the next step.
//
// Why here and not in per-cell behaviors: a_e*E (wrong) vs a_e (correct).
// Dividing by current count would undercount when population is near zero.
// ============================================================================
class SourceBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(SourceBehavior, Behavior, 1);

 public:
  void Run(Agent* agent) override {
    auto* sim  = Simulation::GetActive();
    auto* ctxt = sim->GetExecutionContext();
    auto* rng  = sim->GetRandom();
    const auto* sp = SP();

    const real_t dt_day = sp->dt_minutes / 1440.0;
    // Immune trafficking is systemic — always use global counts regardless
    // of use_local_counts, since the reporter cell has no meaningful neighborhood.
    const Counts cnt    = GetGlobalCounts();
    const real_t lo     = sp->min_bound;
    const real_t hi     = sp->max_bound;

    auto rpos = [&]() -> Real3 {
      return {rng->Uniform(lo, hi), rng->Uniform(lo, hi), rng->Uniform(lo, hi)};
    };
    auto spawn_e = [&]() {
      auto* e = new EffectorTCell(rpos());
      e->AddBehavior(new EffectorBehavior());
      ctxt->AddAgent(e);
    };
    auto spawn_n = [&]() {
      auto* n = new NKCell(rpos());
      n->AddBehavior(new NKBehavior());
      ctxt->AddAgent(n);
    };
    auto spawn_h = [&]() {
      auto* h = new HelperTCell(rpos());
      h->AddBehavior(new HelperBehavior());
      ctxt->AddAgent(h);
    };
    auto spawn_r = [&]() {
      auto* r = new TRegCell(rpos());
      r->AddBehavior(new TRegBehavior());
      ctxt->AddAgent(r);
    };

    // --- a_e: constant E influx (Eq. 2.3 — unconditional, no crowding) ---
    if (rng->Uniform(0, 1) < ProbFromRate(sp->e_base_birth, dt_day)) spawn_e();

    // --- a_n: constant N influx (Eq. 2.4) ---
    if (rng->Uniform(0, 1) < ProbFromRate(sp->n_base_birth, dt_day)) spawn_n();

    // --- a_h: constant H influx (Eq. 2.5) ---
    if (rng->Uniform(0, 1) < ProbFromRate(sp->h_base_birth, dt_day)) spawn_h();

    // --- Eq. 2.6 three R source terms ---
    // (1) constant a
    if (rng->Uniform(0, 1) < ProbFromRate(sp->r_base_src, dt_day)) spawn_r();

    // (2) a_r·E — absolute induction by CTLs
    real_t rate_r_e = sp->r_induced_by_E * static_cast<real_t>(cnt.E);
    if (rng->Uniform(0, 1) < ProbFromRate(rate_r_e, dt_day)) spawn_r();

    // (3) b_r·H — absolute induction by Helper T cells
    real_t rate_r_h = sp->r_induced_by_H * static_cast<real_t>(cnt.H);
    if (rng->Uniform(0, 1) < ProbFromRate(rate_r_h, dt_day)) spawn_r();
  }
};

// ============================================================================
// Reporter
// A lightweight invisible agent that runs once per simulated day and writes
// population counts to the CSV via PopulationLogger.
// It reads counts from GlobalCensus (already computed this step by behaviors).
// ============================================================================
class ReporterCell : public Cell {
  BDM_AGENT_HEADER(ReporterCell, Cell, 1);

 public:
  ReporterCell() { SetDiameter(0.1); }
};

class ReportPopCounts : public Behavior {
  BDM_BEHAVIOR_HEADER(ReportPopCounts, Behavior, 1);

 public:
  void Run(Agent* /*unused*/) override {
    auto* sim = Simulation::GetActive();
    const auto* sp = SP();
    const size_t steps     = sim->GetScheduler()->GetSimulatedSteps();
    const real_t dt_day    = sp->dt_minutes / 1440.0;
    const real_t t_day     = dt_day * static_cast<real_t>(steps);
    const size_t steps_per_day = static_cast<size_t>(
        std::max(1.0, 1440.0 / sp->dt_minutes));

    if (steps % steps_per_day != 0) return;

    // Reuse the census and drug state already computed this step by behaviors.
    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();
    const Counts& c = gc.Get();

    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();

    const real_t day_since_start = t_day + 7.0;  // paper starts at day 7
    std::cout << "[day " << day_since_start << "] "
              << "C=" << c.C << " P=" << c.P
              << " E=" << c.E << " N=" << c.N
              << " H=" << c.H << " R=" << c.R;
    if (ds.M_gem > 1e-9 || ds.M_abr > 1e-9)
      std::cout << "  M_gem=" << ds.M_gem << " M_abr=" << ds.M_abr;
    std::cout << "\n";

    PopulationLogger::Instance().Write(
        steps, day_since_start, c.C, c.P, c.E, c.N, c.H, c.R,
        ds.M_gem, ds.M_abr);
  }
};

// ============================================================================
// Simulate
// ============================================================================
inline int Simulate(int argc, const char** argv) {
  // Register SimParam with BioDynaMo's param system BEFORE constructing
  // Simulation. Param::Param() copies all registered groups into the live
  // param store; Get<SimParam>() will segfault if this is skipped.
  Param::RegisterParamGroup(new SimParam());

  // Phase 1: read JSON into a temporary SimParam to extract bounds/dt for the
  // setup lambda. The actual params are loaded into sp after Simulation is up.
  // Optional --config <file> flag overrides the default params.json path.
  std::string config_file  = "params.json";
  std::string output_override;
  for (int i = 1; i + 1 < argc; ++i) {
    if (std::string(argv[i]) == "--config")  config_file     = argv[i + 1];
    if (std::string(argv[i]) == "--output")  output_override = argv[i + 1];
  }

  SimParam tmp;
  tmp.LoadParams(config_file);

  auto setp = [&tmp](Param* param) {
    param->bound_space          = Param::BoundSpaceMode::kClosed;
    param->min_bound            = tmp.min_bound;
    param->max_bound            = tmp.max_bound;
    param->simulation_time_step = tmp.dt_minutes;
  };

  Simulation sim(argc, argv, setp);
  sim.GetRandom()->SetSeed(tmp.seed);

  // Phase 2: overwrite the live SimParam with the JSON-loaded values so SP()
  // returns correct values in all behaviors.
  auto* sp = const_cast<SimParam*>(sim.GetParam()->Get<SimParam>());
  sp->LoadParams(config_file);
  if (!output_override.empty()) sp->output_dir = output_override;
  sp->PrintParams();

  // Disable mechanical forces — spatial positions are for local counting only;
  // forces would push cells out of the tumor sphere and add noise without
  // adding biological value until a proper spatial model is implemented.
  auto* scheduler = sim.GetScheduler();
  auto mech_ops = scheduler->GetOps("mechanical forces");
  if (!mech_ops.empty()) {
    scheduler->UnscheduleOp(mech_ops[0]);
  }

  // In local mode, ensure the grid box_length >= local_radius_um.
  // BioDynaMo's UniformGridEnvironment auto-sets box_length = max cell diameter
  // (12 µm), but ForEachNeighbor fatals when search_radius > box_length.
  // SetBoxLength() with is_custom_box_length_=true prevents auto-reset each step.
  if (sp->use_local_counts) {
    auto* env = dynamic_cast<UniformGridEnvironment*>(sim.GetEnvironment());
    if (env) {
      int32_t bl = static_cast<int32_t>(std::ceil(sp->local_radius_um));
      env->SetBoxLength(bl);
    }
  }

  // Open the output CSV.
  PopulationLogger::Instance().Open(sp->output_dir);

  auto* ctxt = sim.GetExecutionContext();
  auto* rng  = sim.GetRandom();

  // Helper: random position within the domain.
  const real_t lo = sp->min_bound, hi = sp->max_bound;
  auto rand_pos = [&]() -> Real3 {
    return {rng->Uniform(lo, hi), rng->Uniform(lo, hi), rng->Uniform(lo, hi)};
  };

  // Random position inside a sphere of given radius centered at domain center.
  // Uses rejection sampling; with r << domain side, average ~2 attempts per sample.
  const Real3 domain_center = {(lo + hi) / 2.0, (lo + hi) / 2.0, (lo + hi) / 2.0};
  auto rand_sphere_pos = [&](real_t radius) -> Real3 {
    while (true) {
      real_t x = rng->Uniform(-radius, radius);
      real_t y = rng->Uniform(-radius, radius);
      real_t z = rng->Uniform(-radius, radius);
      if (x * x + y * y + z * z <= radius * radius) {
        return ClampPoint({domain_center[0] + x,
                           domain_center[1] + y,
                           domain_center[2] + z}, lo, hi);
      }
    }
  };

  // Tumor (C) and PSC (P) start as a compact sphere when use_sphere_init=true.
  // This is the biologically correct initial condition for local-mode experiments:
  // the tumor mass is already formed, and immune cells infiltrate from the periphery.
  auto pos_tumor = [&]() -> Real3 {
    return sp->use_sphere_init
        ? rand_sphere_pos(sp->tumor_sphere_radius)
        : rand_pos();
  };

  // --- Seed populations ---
  for (size_t i = 0; i < sp->C0; ++i) {
    auto* c = new TumorCell(pos_tumor());
    c->AddBehavior(new TumorBehavior());
    ctxt->AddAgent(c);
  }
  for (size_t i = 0; i < sp->P0; ++i) {
    auto* p = new StellateCell(pos_tumor());
    p->AddBehavior(new PSCBehavior());
    ctxt->AddAgent(p);
  }
  for (size_t i = 0; i < sp->E0; ++i) {
    auto* e = new EffectorTCell(ClampPoint(rand_pos(), lo, hi));
    e->AddBehavior(new EffectorBehavior());
    ctxt->AddAgent(e);
  }
  for (size_t i = 0; i < sp->N0; ++i) {
    auto* n = new NKCell(ClampPoint(rand_pos(), lo, hi));
    n->AddBehavior(new NKBehavior());
    ctxt->AddAgent(n);
  }
  for (size_t i = 0; i < sp->H0; ++i) {
    auto* h = new HelperTCell(ClampPoint(rand_pos(), lo, hi));
    h->AddBehavior(new HelperBehavior());
    ctxt->AddAgent(h);
  }
  for (size_t i = 0; i < sp->R0; ++i) {
    auto* r = new TRegCell(ClampPoint(rand_pos(), lo, hi));
    r->AddBehavior(new TRegBehavior());
    ctxt->AddAgent(r);
  }

  // Reporter + constant source (invisible agent; not using AlwaysCopyToNew)
  auto* rep = new ReporterCell();
  rep->AddBehavior(new ReportPopCounts());
  rep->AddBehavior(new SourceBehavior());
  ctxt->AddAgent(rep);

  const size_t total_steps =
      static_cast<size_t>(sp->total_days * 1440.0 / sp->dt_minutes);
  sim.GetScheduler()->Simulate(total_steps);

  std::cout << "Pancreatic tumor ABM completed ("
            << sp->total_days << " days).\n";
  return 0;
}

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PANCREATIC_TUMOR_MODEL_H_
