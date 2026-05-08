// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey
// BioDynaMo collaboration. Apache-2.0 license.
//
// -----------------------------------------------------------------------------
#ifndef PANCREATIC_TUMOR_MODEL_H_
#define PANCREATIC_TUMOR_MODEL_H_

#include "biodynamo.h"
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
    csv << "step,days,C,P,E,N,H,R,total\n";
    csv.flush();
  }

  void Write(size_t step, real_t t_day,
             size_t C, size_t P, size_t E, size_t N, size_t H, size_t R) {
    csv << step << "," << t_day << ","
        << C << "," << P << "," << E << ","
        << N << "," << H << "," << R << ","
        << (C + P + E + N + H + R) << "\n";
    csv.flush();
  }
};

// ============================================================================
// Behaviors
//
// Key patterns (from CARTopiaX):
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
    const Counts cnt = GetCounts(*c);

    // --- Division (Eq. 2.1 growth terms) ---
    // (k_c + mu_c·P)·C·(1-C/K_C): mu_c·P is LINEAR in P (not Hill)
    real_t crowd  = 1.0 - Clamp(static_cast<real_t>(cnt.C) / sp->K_C, 0.0, 1.0);
    real_t boostP = sp->c_boost_from_P * static_cast<real_t>(cnt.P);  // mu_c·P
    real_t div_rate = (sp->c_base_div + boostP) * crowd;

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      c->Divide();
      c->SetCellColor(sp->color_tumor_div);
    } else {
      c->SetCellColor(sp->color_tumor);
    }

    // --- Killing (Eq. 2.1 death terms) — independent of division ---
    // b_c·N·C (bilinear in N) and d_c·E·C/(1+r1·R) (bilinear in E)
    real_t inhibit_R = 1.0 / (1.0 + sp->c_R_blocks_E * static_cast<real_t>(cnt.R));
    real_t killE = sp->c_kill_by_E * static_cast<real_t>(cnt.E) * inhibit_R;
    real_t killN = sp->c_kill_by_N * static_cast<real_t>(cnt.N);

    if (rng->Uniform(0, 1) < ProbFromRate(killE + killN, dt_day)) {
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
    const Counts cnt = GetCounts(*psc);

    // --- Death (lambda_p*P) — checked independently ---
    if (rng->Uniform(0, 1) < ProbFromRate(sp->p_base_death, dt_day)) {
      ctxt->RemoveAgent(psc->GetUid());
      return;
    }

    // --- Division (Eq. 2.2) ---
    // (k_p + f_p*C/(mu_p+C))*P*(1-a_p*P)
    real_t crowd  = 1.0 - Clamp(static_cast<real_t>(cnt.P) / sp->K_P, 0.0, 1.0);
    real_t boostC = sp->p_boost_from_C
                  * Sat(static_cast<real_t>(cnt.C), sp->p_boost_from_C_K);
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

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const Counts cnt = GetCounts(*e);

    // --- Death (Eq. 2.3 death terms): b_e·E + c_e·E·C + δ_e·R·E ---
    // c_e·C and δ_e·R are bilinear (no Hill).
    real_t die = sp->e_base_death
               + sp->e_inact_by_C * static_cast<real_t>(cnt.C)   // c_e·C
               + sp->e_suppr_by_R * static_cast<real_t>(cnt.R);  // δ_e·R

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(e->GetUid());
      return;
    }

    // --- Per-cell proliferation (Eq. 2.3, p_e·H·E/(g_e+H) term) ---
    // Constant influx a_e handled by SourceBehavior.
    real_t div_rate = sp->e_help_from_H
                    * Sat(static_cast<real_t>(cnt.H), sp->e_help_from_H_K);

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

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const Counts cnt = GetCounts(*n);

    // --- Death (Eq. 2.4 death terms): b_n·N + c_n·N·C + δ_n·R·N ---
    real_t die = sp->n_base_death
               + sp->n_inact_by_C * static_cast<real_t>(cnt.C)   // c_n·C (bilinear)
               + sp->n_suppr_by_R * static_cast<real_t>(cnt.R);  // δ_n·R (bilinear)

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(n->GetUid());
      return;
    }

    // --- Per-cell proliferation (Eq. 2.4, p_n·H·N/(g_n+H) term) ---
    // Constant influx a_n handled by SourceBehavior.
    real_t div_rate = sp->n_help_from_H
                    * Sat(static_cast<real_t>(cnt.H), sp->n_help_from_H_K);

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

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const Counts cnt = GetCounts(*h);

    // --- Death (Eq. 2.5 death terms): b_h·H + δ_h·R·H ---
    // δ_h·R is bilinear in R.
    real_t die = sp->h_base_death
               + sp->h_suppr_by_R * static_cast<real_t>(cnt.R);  // δ_h·R (bilinear)

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(h->GetUid());
      return;
    }

    // --- Per-cell self-activation (Eq. 2.5, p_h·H²/(g_h+H) = p_h·H/(g_h+H) per cell) ---
    // Constant influx a_h handled by SourceBehavior.
    real_t div_rate = sp->h_self_act
                    * Sat(static_cast<real_t>(cnt.H), sp->h_self_act_K);

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

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const Counts cnt = GetCounts(*r);

    // --- Death (Eq. 2.6 death terms): δ_r·R + r·N·R ---
    // r·N is bilinear in N.
    real_t die = sp->r_decay
               + sp->r_cleared_by_N * static_cast<real_t>(cnt.N);  // r·N (bilinear)

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(r->GetUid());
      return;
    }

    // --- Per-cell proliferation (Eq. 2.6, p_r·H·R/(g_r+H) term) ---
    // Absolute sources a, a_r·E, b_r·H handled by SourceBehavior.
    real_t div_rate = sp->r_prolif_by_H
                    * Sat(static_cast<real_t>(cnt.H), sp->r_prolif_by_H_K);

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
    const Counts cnt    = GetCounts(*agent);
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

    // Reuse the census already computed this step by the cell behaviors.
    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();
    const Counts& c = gc.Get();

    const real_t day_since_start = t_day + 7.0;  // paper starts at day 7
    std::cout << "[day " << day_since_start << "] "
              << "C=" << c.C << " P=" << c.P
              << " E=" << c.E << " N=" << c.N
              << " H=" << c.H << " R=" << c.R << "\n";

    PopulationLogger::Instance().Write(
        steps, day_since_start, c.C, c.P, c.E, c.N, c.H, c.R);
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
  SimParam tmp;
  tmp.LoadParams("params.json");

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
  *sp = tmp;
  sp->PrintParams();

  // Disable mechanical forces: this model uses global ODE replication mode where
  // spatial structure is irrelevant. Mechanics only add noise and slow the run.
  // (When local-mode spatial experiments are introduced, remove this line.)
  auto* scheduler = sim.GetScheduler();
  auto mech_ops = scheduler->GetOps("mechanical forces");
  if (!mech_ops.empty()) {
    scheduler->UnscheduleOp(mech_ops[0]);
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

  // --- Seed populations ---
  for (size_t i = 0; i < sp->C0; ++i) {
    auto* c = new TumorCell(ClampPoint(rand_pos(), lo, hi));
    c->AddBehavior(new TumorBehavior());
    ctxt->AddAgent(c);
  }
  for (size_t i = 0; i < sp->P0; ++i) {
    auto* p = new StellateCell(ClampPoint(rand_pos(), lo, hi));
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
