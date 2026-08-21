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
#include <climits>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <random>
#include <vector>

#include <unistd.h>  // readlink — locate the pvsm fixer relative to the binary

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
// Rejection-sample a random position inside a sphere (viz sphere-seed mode).
inline Real3 RandSpherePos(Random* rng, const Real3& center, real_t radius,
                           real_t lo, real_t hi) {
  while (true) {
    const real_t x = rng->Uniform(-radius, radius);
    const real_t y = rng->Uniform(-radius, radius);
    const real_t z = rng->Uniform(-radius, radius);
    if (x * x + y * y + z * z <= radius * radius) {
      return ClampPoint({center[0] + x, center[1] + y, center[2] + z}, lo, hi);
    }
  }
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

// On cell division BioDynaMo places the daughter touching the mother; with
// mechanical forces disabled and a volume split the pair can render one inside
// the other. Instead place the daughter a short fixed distance from the mother
// in a random direction: it stays a local neighbour (biologically sensible) but
// is clearly separated. Fires only for division events, not other new-agent ones.
// Visualization-friendly division for the DAUGHTER (viz_division only): pin it
// to a uniform size (BioDynaMo's Divide() otherwise splits volume and shrinks
// it) and place it a short distance from the mother so the pair doesn't overlap.
// No-op when viz_division=false → BioDynaMo's default (shrunk + adjacent).
inline void VizDivideDaughter(Cell* daughter, const NewAgentEvent& event,
                              real_t diameter) {
  if (!SP()->viz_division) return;
  if (event.GetUid() != CellDivisionEvent::kUid) return;
  daughter->SetDiameter(diameter);                    // uniform render size
  auto* rng = Simulation::GetActive()->GetRandom();
  const auto* sp = SP();
  const Real3 mpos = event.existing_agent->GetPosition();
  Real3 dir = {rng->Uniform(-1.0, 1.0), rng->Uniform(-1.0, 1.0),
               rng->Uniform(-1.0, 1.0)};
  real_t len = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  if (len < 1e-9) { dir = {1.0, 0.0, 0.0}; len = 1.0; }
  // Center-to-center = 1.1 diameters: the daughter sits right against the mother
  // (surfaces just touching, a tiny gap) in a random direction — adjacent like a
  // bud, never one inside the other. (=1 diameter would touch exactly; <1 overlaps.)
  const real_t dist = 2 * daughter->GetDiameter();
  const Real3 pos = {mpos[0] + dir[0] / len * dist,
                     mpos[1] + dir[1] / len * dist,
                     mpos[2] + dir[2] / len * dist};
  daughter->SetPosition(ClampPoint(pos, sp->min_bound, sp->max_bound));
}
// Reset the MOTHER's size after Divide() (viz_division only; no-op otherwise).
inline void VizResetMotherSize(Cell* mother, real_t diameter) {
  if (SP()->viz_division) mother->SetDiameter(diameter);
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
    SetDiameter(SP()->cell_radius_um);
    color_ = SP()->color_tumor;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    VizDivideDaughter(this, event, SP()->cell_radius_um);
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
    SetDiameter(SP()->cell_radius_um);
    color_ = SP()->color_psc;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    VizDivideDaughter(this, event, SP()->cell_radius_um);
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
    SetDiameter(SP()->cell_radius_um);
    color_ = SP()->color_eff;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    VizDivideDaughter(this, event, SP()->cell_radius_um);
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
    SetDiameter(SP()->cell_radius_um);
    color_ = SP()->color_nk;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    VizDivideDaughter(this, event, SP()->cell_radius_um);
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
    SetDiameter(SP()->cell_radius_um);
    color_ = SP()->color_helper;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    VizDivideDaughter(this, event, SP()->cell_radius_um);
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
    SetDiameter(SP()->cell_radius_um);
    color_ = SP()->color_treg;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    VizDivideDaughter(this, event, SP()->cell_radius_um);
    color_ = bdm_static_cast<TRegCell*>(event.existing_agent)->color_;
  }

  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }

 private:
  int color_ = 0;
};

// S — Cancer Stem Cell (CSC), Paper Section 6. OPTIONAL: only seeded and grown
// when csc_enable=true.  Inert otherwise (no agents created).
class CancerStemCell : public Cell {
  BDM_AGENT_HEADER(CancerStemCell, Cell, 1);

 public:
  CancerStemCell() = default;
  explicit CancerStemCell(const Real3& p) {
    SetPosition(p);
    SetDiameter(SP()->cell_radius_um);
    color_ = SP()->color_csc;
  }

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    VizDivideDaughter(this, event, SP()->cell_radius_um);
    color_ = bdm_static_cast<CancerStemCell*>(event.existing_agent)->color_;
  }

  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }

 private:
  int color_ = 0;
};

// ============================================================================
// Population census
// ============================================================================
struct Counts { size_t C=0, P=0, E=0, N=0, H=0, R=0, S=0; };  // S = Cancer Stem Cells

// Global census: computed at most once per step via double-checked locking.
// Thread-safe: multiple behaviors call RefreshIfNeeded() concurrently but only
// one executes the ForEachAgent scan, the rest read the cached result.
struct GlobalCensus {
  std::atomic<size_t> step_cached{std::numeric_limits<size_t>::max()};
  Counts              cnt;
  real_t              tumor_radius = 0.0;  // max tumor-cell dist from center (viz)
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

    // In viz sphere-seed mode also track the tumor's extent so immune cells can
    // follow the growing mass (see ImmuneRandomWalk). Skipped otherwise.
    const auto* sp = SP();
    const bool   track_r = sp->viz_sphere_seed;
    const real_t center  = (sp->min_bound + sp->max_bound) / 2.0;
    Counts c;
    real_t sum_r2 = 0.0;
    sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
      if (dynamic_cast<TumorCell*>(a)) {
        ++c.C;
        if (track_r) {
          const Real3 p = a->GetPosition();
          const real_t dx = p[0] - center, dy = p[1] - center, dz = p[2] - center;
          sum_r2 += dx * dx + dy * dy + dz * dz;
        }
      }
      else if (dynamic_cast<StellateCell*>(a))  ++c.P;
      else if (dynamic_cast<EffectorTCell*>(a)) ++c.E;
      else if (dynamic_cast<NKCell*>(a))        ++c.N;
      else if (dynamic_cast<HelperTCell*>(a))   ++c.H;
      else if (dynamic_cast<TRegCell*>(a))      ++c.R;
      else if (dynamic_cast<CancerStemCell*>(a)) ++c.S;
    });
    cnt = c;
    // Effective tumor ball radius from the RMS distance (robust vs a few
    // outliers): for a uniform ball, mean(r^2) = (3/5)R^2 -> R = rms * 1.291.
    tumor_radius = (c.C > 0) ? std::sqrt(sum_r2 / static_cast<real_t>(c.C)) * 1.291
                             : 0.0;
    step_cached.store(step, std::memory_order_release);
  }

  const Counts& Get() const { return cnt; }
  real_t TumorRadius() const { return tumor_radius; }
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
      else if (dynamic_cast<CancerStemCell*>(nb)) ++out->S;
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

// Immune random walk: one Gaussian step per step. Active in local mode (cells
// infiltrate the tumor) and in viz_sphere_seed mode. In viz mode the step is
// reflected inward at a radius that GROWS with the tumor (max of the seeding
// sphere and the live tumor extent), so immune cells expand together with the
// tumor mass rather than being stuck in the initial sphere. Always kept inside
// the domain box.
inline void ImmuneRandomWalk(Cell* cell, const SimParam* sp, Random* rng) {
  if (!sp->use_local_counts && !sp->viz_sphere_seed) return;
  const Real3 pos = cell->GetPosition();
  const real_t s  = sp->immune_step_um;
  Real3 np = {pos[0] + rng->Gaus(0.0, s),
              pos[1] + rng->Gaus(0.0, s),
              pos[2] + rng->Gaus(0.0, s)};
  if (sp->viz_sphere_seed) {
    const real_t c     = (sp->min_bound + sp->max_bound) / 2.0;
    const real_t seedR = sp->viz_seed_radius_frac * (sp->max_bound - sp->min_bound) / 2.0;
    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();
    // Follow the tumor: allow immune cells out to the current tumor extent
    // (+ a couple of steps of margin), but never less than the seed sphere.
    const real_t R = std::max(seedR, gc.TumorRadius() + 2.0 * s);
    const real_t dx = np[0] - c, dy = np[1] - c, dz = np[2] - c;
    const real_t d2 = dx * dx + dy * dy + dz * dz;
    if (d2 > R * R) {                          // stepped past the front — reflect in
      const real_t d = std::sqrt(d2);
      const real_t scale = std::max(0.0, 2.0 * R - d) / d;
      np = {c + dx * scale, c + dy * scale, c + dz * scale};
    }
  }
  cell->SetPosition(ClampPoint(np, sp->min_bound, sp->max_bound));
}

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
// DrugState  (Paper Section 5) — OPTIONAL, only active when treat_* flags set.
//
// Holds drug concentrations M_gem, M_abr and the Anti-CD47 active flag.
// Refreshed at most once per step (double-checked locking like GlobalCensus):
//   M *= exp(-γ·dt)               (exponential PK decay)
//   M += dose on injection steps  (integer-step grid, no float drift)
// When no treatment flag is set, RefreshIfNeeded() short-circuits and all
// concentrations stay 0 — so the base model path is completely unaffected.
// ============================================================================
// Diffusion-grid substance ids for the optional drug-diffusion visualization.
enum DrugSubstance : int { kGemSubstance = 0, kAbrSubstance = 1 };

struct DrugState {
  std::atomic<size_t> step_cached{std::numeric_limits<size_t>::max()};
  std::mutex mtx;

  real_t M_gem = 0.0;
  real_t M_abr = 0.0;
  bool   acd47_active = false;
  // Dose injected THIS step (0 otherwise) — drives the optional drug-diffusion
  // visualization grid (DrugDiffusionInjector). Not used by the kill dynamics.
  real_t gem_injected = 0.0;
  real_t abr_injected = 0.0;
  // True once the step is past the last active treatment's end day — triggers
  // the paper's post-treatment reduced growth rate kc_post_treat (Section 5).
  bool   post_treatment = false;
  // Step-averaged kill fraction (1-e^{-M}) over the timestep (see AvgKillFrac).
  // Behaviors use these instead of (1-e^{-M}) at the frozen peak — otherwise a
  // fast-decaying drug (gemcitabine, t½=3h) is grossly over-applied when dt is
  // not << t½ (dt=24h → ~5.6x over-kill → tumor crashes vs paper Fig. 5).
  real_t gem_kill_frac = 0.0;
  real_t abr_kill_frac = 0.0;

  static DrugState& Instance() { static DrugState ds; return ds; }

  // Time-averaged (1-e^{-M}) over one step, integrating the within-step decay
  // M(τ)=M0·e^{-γτ}. Midpoint quadrature (N=16) — matches the continuous kill
  // integral the ODE reference already uses, making the delivered dose per
  // injection dt-independent.
  static real_t AvgKillFrac(real_t M0, real_t gamma, real_t dt) {
    if (M0 <= 0.0) return 0.0;
    if (gamma <= 0.0 || dt <= 0.0) return 1.0 - std::exp(-M0);
    constexpr int N = 16;
    real_t sum = 0.0;
    for (int i = 0; i < N; ++i) {
      const real_t tau = (i + 0.5) * dt / N;
      sum += 1.0 - std::exp(-M0 * std::exp(-gamma * tau));
    }
    return sum / N;
  }

  void RefreshIfNeeded() {
    auto* sim   = Simulation::GetActive();
    size_t step = sim->GetScheduler()->GetSimulatedSteps();
    if (step_cached.load(std::memory_order_acquire) == step) return;

    std::lock_guard<std::mutex> lock(mtx);
    if (step_cached.load(std::memory_order_relaxed) == step) return;

    const auto* sp = SP();
    // Short-circuit: no treatment enabled → base model, concentrations stay 0.
    if (!sp->treat_gem && !sp->treat_abr && !sp->treat_acd47) {
      step_cached.store(step, std::memory_order_release);
      return;
    }

    const real_t dt_day  = sp->dt_minutes / 1440.0;
    const real_t day_off = 7.0;                  // paper day 7 = simulation day 0
    const real_t spd_f   = 1440.0 / sp->dt_minutes;  // steps per day

    auto s_of = [&](real_t pday) -> size_t {
      real_t sim_day = pday - day_off;
      return (sim_day <= 0.0) ? 0 : static_cast<size_t>(std::round(sim_day * spd_f));
    };
    auto is_inject = [&](size_t s0, size_t s1, size_t freq) -> bool {
      return step >= s0 && step <= s1 && ((step - s0) % freq == 0);
    };

    gem_injected = 0.0;
    if (sp->treat_gem) {
      M_gem *= std::exp(-sp->gem_gamma * dt_day);
      size_t s0   = s_of(sp->treat_start_day);
      size_t s1   = s_of(sp->gem_end_day);
      size_t freq = static_cast<size_t>(std::round(sp->gem_freq_days * spd_f));
      if (freq > 0 && is_inject(s0, s1, freq)) {
        M_gem += sp->gem_dose;
        gem_injected = sp->gem_dose;
      }
    }

    abr_injected = 0.0;
    if (sp->treat_abr) {
      M_abr *= std::exp(-sp->abr_gamma * dt_day);
      size_t s0   = s_of(sp->treat_start_day);
      size_t s1   = s_of(sp->abr_end_day);
      size_t freq = static_cast<size_t>(std::round(sp->abr_freq_days * spd_f));
      if (freq > 0 && is_inject(s0, s1, freq)) {
        M_abr += sp->abr_dose;
        abr_injected = sp->abr_dose;
      }
    }

    acd47_active = sp->treat_acd47
                && step >= s_of(sp->treat_start_day)
                && step <= s_of(sp->acd47_end_day);

    // Post-treatment begins after the LAST active treatment's end day.
    real_t last_end = 0.0;
    if (sp->treat_gem)   last_end = std::max(last_end, sp->gem_end_day);
    if (sp->treat_abr)   last_end = std::max(last_end, sp->abr_end_day);
    if (sp->treat_acd47) last_end = std::max(last_end, sp->acd47_end_day);
    post_treatment = (last_end > 0.0) && (step > s_of(last_end));

    // Step-averaged kill fractions (resolves fast PK independent of dt).
    gem_kill_frac = AvgKillFrac(M_gem, sp->gem_gamma, dt_day);
    abr_kill_frac = AvgKillFrac(M_abr, sp->abr_gamma, dt_day);

    step_cached.store(step, std::memory_order_release);
  }
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
    // Unified schema: S is 0 unless CSC enabled; M_gem/M_abr are 0 unless
    // treatment enabled. Python readers select columns by name, so the extra
    // columns are harmless for base-model runs.
    csv << "step,days,C,P,E,N,H,R,S,total,M_gem,M_abr\n";
    csv.flush();
  }

  void Write(size_t step, real_t t_day,
             size_t C, size_t P, size_t E, size_t N, size_t H, size_t R,
             size_t S = 0, real_t M_gem = 0.0, real_t M_abr = 0.0) {
    csv << step << "," << t_day << ","
        << C << "," << P << "," << E << ","
        << N << "," << H << "," << R << "," << S << ","
        << (C + P + E + N + H + R + S) << ","
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
    // Post-treatment: kc_post_treat is the FITTED EFFECTIVE growth rate ("k_c
    // after treatment", Fig. 5 titles) and is used as the TOTAL growth — the
    // mu_c*P PSC boost is dropped. Otherwise mu_c*P dominates the tiny fitted
    // rate as P relapses and C climbs instead of plateauing (abr+acd47, Fig. 5d).
    real_t div_rate;
    if (ds.post_treatment && sp->kc_post_treat > 0.0) {
      div_rate = sp->kc_post_treat * crowd;
    } else {
      real_t boostP = sp->c_boost_from_P * cnt.P;  // mu_c*P (Eq. 5.1)
      div_rate = (sp->c_base_div + boostP) * crowd;
    }

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      c->Divide();
      VizResetMotherSize(c, sp->cell_radius_um);
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

    // --- Drug kill (Eq. 5.1): c_c·(1−e^{−M}) — 0 when treatment disabled ---
    real_t drug_kill = sp->gem_c_c * ds.gem_kill_frac
                     + sp->abr_c_c * ds.abr_kill_frac;

    if (rng->Uniform(0, 1) < ProbFromRate(killE + killN + drug_kill, dt_day)) {
      ctxt->RemoveAgent(c->GetUid());
      // c->SetCellColor(0);
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

    // --- Death (lambda_p*P + drug kill Eq. 5.1) — checked independently ---
    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();
    real_t drug_kill = sp->gem_c_p * ds.gem_kill_frac
                     + sp->abr_c_p * ds.abr_kill_frac;
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
      VizResetMotherSize(psc, sp->cell_radius_um);
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
    ImmuneRandomWalk(e, sp, rng);  // move each step (bounded to sphere in viz mode)

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const EffCounts cnt(GetCounts(*e), DensityCompensation());

    // --- Death (Eq. 2.3 death terms): b_e·E + c_e·E·C + δ_e·R·E ---
    // c_e·C is local (direct tumor contact). δ_e·R uses global R (cytokine-mediated).
    const real_t R_global_e = sp->use_local_counts
        ? static_cast<real_t>(GetGlobalCounts().R) : cnt.R;
    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();
    real_t drug_kill = sp->gem_c_immune * ds.gem_kill_frac
                     + sp->abr_c_immune * ds.abr_kill_frac;
    real_t die = sp->e_base_death
               + sp->e_inact_by_C * cnt.C       // c_e·C (local)
               + sp->e_suppr_by_R * R_global_e  // δ_e·R (global)
               + drug_kill;                      // Eq. 5.1 (0 if no treatment)

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(e->GetUid());
      return;
    }

    // --- Per-cell proliferation (Eq. 2.3, p_e·H·E/(g_e+H) term) ---
    // Constant influx a_e handled by SourceBehavior.
    real_t div_rate = sp->e_help_from_H * Sat(cnt.H, sp->e_help_from_H_K);
    // Anti-CD47 (Eq. 5.x): extra CTL proliferation while active (0 otherwise).
    if (ds.acd47_active) div_rate += sp->acd47_e_boost;

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      e->Divide();
      VizResetMotherSize(e, sp->cell_radius_um);
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

    ImmuneRandomWalk(n, sp, rng);  // move each step (bounded to sphere in viz mode)

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const EffCounts cnt(GetCounts(*n), DensityCompensation());

    // --- Death (Eq. 2.4 death terms): b_n·N + c_n·N·C + δ_n·R·N ---
    const real_t R_global_n = sp->use_local_counts
        ? static_cast<real_t>(GetGlobalCounts().R) : cnt.R;
    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();
    real_t drug_kill = sp->gem_c_immune * ds.gem_kill_frac
                     + sp->abr_c_immune * ds.abr_kill_frac;
    real_t die = sp->n_base_death
               + sp->n_inact_by_C * cnt.C       // c_n·C (local, direct contact)
               + sp->n_suppr_by_R * R_global_n  // δ_n·R (global, cytokine)
               + drug_kill;                      // Eq. 5.1 (0 if no treatment)

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(n->GetUid());
      return;
    }

    // --- Per-cell proliferation (Eq. 2.4, p_n·H·N/(g_n+H) term) ---
    // Constant influx a_n handled by SourceBehavior.
    real_t div_rate = sp->n_help_from_H * Sat(cnt.H, sp->n_help_from_H_K);

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      n->Divide();
      VizResetMotherSize(n, sp->cell_radius_um);
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

    ImmuneRandomWalk(h, sp, rng);  // move each step (bounded to sphere in viz mode)

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const EffCounts cnt(GetCounts(*h), DensityCompensation());

    // --- Death (Eq. 2.5 death terms): b_h·H + δ_h·R·H ---
    const real_t R_global_h = sp->use_local_counts
        ? static_cast<real_t>(GetGlobalCounts().R) : cnt.R;
    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();
    real_t drug_kill = sp->gem_c_immune * ds.gem_kill_frac
                     + sp->abr_c_immune * ds.abr_kill_frac;
    real_t die = sp->h_base_death
               + sp->h_suppr_by_R * R_global_h  // δ_h·R (global, cytokine)
               + drug_kill;                      // Eq. 5.1 (0 if no treatment)

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(h->GetUid());
      return;
    }

    // --- Per-cell self-activation (Eq. 2.5, p_h·H²/(g_h+H) = p_h·H/(g_h+H) per cell) ---
    // Constant influx a_h handled by SourceBehavior.
    real_t div_rate = sp->h_self_act * Sat(cnt.H, sp->h_self_act_K);

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      h->Divide();
      VizResetMotherSize(h, sp->cell_radius_um);
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

    ImmuneRandomWalk(r, sp, rng);  // move each step (bounded to sphere in viz mode)

    const real_t dt_day = sp->dt_minutes / 1440.0;
    const EffCounts cnt(GetCounts(*r), DensityCompensation());

    // --- Death (Eq. 2.6 death terms): δ_r·R + r·N·R + drug kill (Eq. 5.1) ---
    // r·N is bilinear in N.
    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();
    real_t drug_kill = sp->gem_c_immune * ds.gem_kill_frac
                     + sp->abr_c_immune * ds.abr_kill_frac;
    real_t die = sp->r_decay
               + sp->r_cleared_by_N * cnt.N  // r·N (bilinear)
               + drug_kill;                   // Eq. 5.1 (0 if no treatment)

    if (rng->Uniform(0, 1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(r->GetUid());
      return;
    }

    // --- Per-cell proliferation (Eq. 2.6, p_r·H·R/(g_r+H) term) ---
    // Absolute sources a, a_r·E, b_r·H handled by SourceBehavior.
    real_t div_rate = sp->r_prolif_by_H * Sat(cnt.H, sp->r_prolif_by_H_K);

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      r->Divide();
      VizResetMotherSize(r, sp->cell_radius_um);
    }
  }
};

// ============================================================================
// CSCBehavior  (Paper Section 6, Eq. 6.7) — OPTIONAL, only runs when csc_enable.
// Self-renewal rate switches via an arctan in C:
//   C ≪ σ (tumor suppressed)  → λ ≈ λ_max  (CSCs proliferate fast → relapse)
//   C ≫ σ (tumor large)       → λ ≈ λ_min  (CSCs grow minimally)
// The (a2+2a3)·S source of new PCCs to Eq. 2.1 is handled by SourceBehavior.
// ============================================================================
class CSCBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(CSCBehavior, Behavior, 1);

 public:
  CSCBehavior() { AlwaysCopyToNew(); }

  void Run(Agent* a) override {
    if (!SP()->csc_enable) return;

    auto* s    = bdm_static_cast<CancerStemCell*>(a);
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto* rng  = Simulation::GetActive()->GetRandom();
    const auto* sp = SP();

    const real_t dt_day = sp->dt_minutes / 1440.0;

    // --- Death (δ_S·S) — paper sets δ_S = 0 for aggressive progression ---
    if (sp->csc_delta_s > 0.0 &&
        rng->Uniform(0, 1) < ProbFromRate(sp->csc_delta_s, dt_day)) {
      ctxt->RemoveAgent(s->GetUid());
      return;
    }

    // --- Self-renewal: λ(C)·(a1−a3)  [Eq. 6.7], global C (all CSCs see same C) ---
    const real_t C_global = static_cast<real_t>(GetGlobalCounts().C);
    real_t lam = -(std::atan(C_global - sp->csc_sigma) + M_PI / 2.0)
                  * (sp->csc_lambda_max - sp->csc_lambda_min) / M_PI
                  + sp->csc_lambda_max;
    real_t div_rate = lam * (sp->csc_a1 - sp->csc_a3);

    if (rng->Uniform(0, 1) < ProbFromRate(div_rate, dt_day)) {
      s->Divide();
      VizResetMotherSize(s, sp->cell_radius_um);
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

    // Viz sphere-seed: keep recruited immune cells inside the sphere too, so the
    // mass stays a spheroid rather than growing a box halo of new cells.
    const Real3 dcenter = {(lo + hi) / 2.0, (lo + hi) / 2.0, (lo + hi) / 2.0};
    const real_t vradius = sp->viz_seed_radius_frac * (hi - lo) / 2.0;
    auto rpos = [&]() -> Real3 {
      if (sp->viz_sphere_seed) return RandSpherePos(rng, dcenter, vradius, lo, hi);
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
    auto spawn_c = [&]() {  // new PCC produced by a CSC division (Eq. 6.1)
      auto* c = new TumorCell(rpos());
      c->AddBehavior(new TumorBehavior());
      ctxt->AddAgent(c);
    };

    // Poisson spawning: deterministically spawn floor(λ) cells and stochastically
    // +1 with probability frac(λ), where λ = rate * dt_day.
    // This is unbiased at any dt — eliminates Bernoulli saturation when λ ≥ 1
    // (which occurs for r_base_src at S≥1e4 even at dt=30 min).
    auto spawn_poisson = [&](auto fn, real_t rate) {
      real_t lam = rate * dt_day;
      int    n   = static_cast<int>(lam);
      if (rng->Uniform(0.0, 1.0) < (lam - static_cast<real_t>(n))) ++n;
      for (int i = 0; i < n; ++i) fn();
    };

    // a_e, a_n, a_h: constant immune influx (Eqs. 2.3–2.5)
    spawn_poisson(spawn_e, sp->e_base_birth);
    spawn_poisson(spawn_n, sp->n_base_birth);
    spawn_poisson(spawn_h, sp->h_base_birth);

    // Eq. 2.6: three Treg source terms
    spawn_poisson(spawn_r, sp->r_base_src);                                        // (1) constant a
    spawn_poisson(spawn_r, sp->r_induced_by_E * static_cast<real_t>(cnt.E));      // (2) a_r·E
    spawn_poisson(spawn_r, sp->r_induced_by_H * static_cast<real_t>(cnt.H));      // (3) b_r·H

    // Eq. 6.1: CSC→PCC source (a2+2a3)·S — 0 when CSC disabled or no CSCs.
    if (sp->csc_enable && cnt.S > 0) {
      spawn_poisson(spawn_c,
                    (sp->csc_a2 + 2.0 * sp->csc_a3) * static_cast<real_t>(cnt.S));
    }
  }
};

// ============================================================================
// DrugDiffusionInjector  (Paper Section 5) — OPTIONAL, viz_drug_diffusion only.
//
// On each dosing step, deposits the dose at a FIXED set of random points spread
// through the whole simulation volume ("vasculature" delivering drug into the
// tissue). Diffusion then smooths the blobs into overlapping clouds, so the
// field is a soft, spatially-textured haze that fills the domain (like the
// soma_clustering substance) and pulses with the Cioffi schedule, decaying
// between injections. PURELY VISUAL — the kill dynamics still use the
// well-mixed DrugState scalar, so Section-5 ODE replication is unchanged.
// The source points use their own RNG (seeded from `seed`) so enabling the
// visualization does not perturb the simulation's own random stream.
// Attached to the reporter agent (single, non-dividing).
// ============================================================================
class DrugDiffusionInjector : public Behavior {
  BDM_BEHAVIOR_HEADER(DrugDiffusionInjector, Behavior, 1);

 public:
  void Run(Agent* /*unused*/) override {
    const auto* sp = SP();
    if (!sp->viz_drug_diffusion) return;
    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();
    auto* rm = Simulation::GetActive()->GetResourceManager();

    if (!built_) BuildSources(sp);

    auto deposit = [&](int id, real_t amount) {
      auto* g = rm->GetDiffusionGrid(id);
      if (!g) return;
      for (const auto& p : sources_) g->ChangeConcentrationBy(p, amount);
    };
    if (sp->treat_gem && ds.gem_injected > 0.0) deposit(kGemSubstance, ds.gem_injected);
    if (sp->treat_abr && ds.abr_injected > 0.0) deposit(kAbrSubstance, ds.abr_injected);
  }

 private:
  bool built_ = false;
  std::vector<Real3> sources_;

  // Fixed random source points scattered through the whole domain volume
  // (systemic drug — covers the tumor as it grows/expands toward the box).
  void BuildSources(const SimParam* sp) {
    const real_t lo = sp->min_bound, hi = sp->max_bound;
    std::mt19937 gen(static_cast<unsigned>(sp->seed) + 12345u);
    std::uniform_real_distribution<real_t> U(lo, hi);
    sources_.reserve(sp->drug_n_sources);
    for (int i = 0; i < sp->drug_n_sources; ++i)
      sources_.push_back({U(gen), U(gen), U(gen)});
    built_ = true;
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

    // Drug concentrations (0 when no treatment enabled).
    auto& ds = DrugState::Instance();
    ds.RefreshIfNeeded();

    const real_t day_since_start = t_day + 7.0;  // paper starts at day 7
    std::cout << "[day " << day_since_start << "] "
              << "C=" << c.C << " P=" << c.P
              << " E=" << c.E << " N=" << c.N
              << " H=" << c.H << " R=" << c.R;
    if (sp->csc_enable)  std::cout << " S=" << c.S;
    if (sp->treat_gem || sp->treat_abr)
      std::cout << " M_gem=" << ds.M_gem << " M_abr=" << ds.M_abr;
    std::cout << "\n";

    PopulationLogger::Instance().Write(
        steps, day_since_start, c.C, c.P, c.E, c.N, c.H, c.R,
        c.S, ds.M_gem, ds.M_abr);
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
  // BDM_PARAMS env var overrides the default "params.json" path.
  const char* _params_env = std::getenv("BDM_PARAMS");
  const std::string _params_path = _params_env ? _params_env : "configs/params.json";
  SimParam tmp;
  tmp.LoadParams(_params_path);

  auto setp = [&tmp](Param* param) {
    param->bound_space          = Param::BoundSpaceMode::kClosed;
    param->min_bound            = tmp.min_bound;
    param->max_bound            = tmp.max_bound;
    param->simulation_time_step = tmp.dt_minutes;
    // Optional drug-diffusion visualization: register the concentration fields
    // for ParaView export (only the substances actually defined in Simulate()).
    // gradient=false: the extra vector array carries an L2_NORM_RANGE
    // InformationKey that some ParaView reader versions fail to parse. Color by
    // "Substance Concentration"; use ParaView's Gradient filter if a vector is
    // needed. {name, concentration, gradient}.
    if (tmp.viz_drug_diffusion) {
      if (tmp.treat_gem) param->visualize_diffusion.push_back({"Gemcitabine", true, false});
      if (tmp.treat_abr) param->visualize_diffusion.push_back({"Abraxane", true, false});
    }
  };

  Simulation sim(argc, argv, setp);
  sim.GetRandom()->SetSeed(tmp.seed);

  // Phase 2: overwrite the live SimParam with the JSON-loaded values so SP()
  // returns correct values in all behaviors.
  auto* sp = const_cast<SimParam*>(sim.GetParam()->Get<SimParam>());
  *sp = tmp;
  // Sphere-seed visualization runs in GLOBAL mode: positions become purely
  // cosmetic (dynamics depend only on total counts = the validated mean-field),
  // so seeding a spheroid changes nothing scientific.
  if (sp->viz_sphere_seed) sp->use_local_counts = false;
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

  // Optional drug-diffusion visualization grids (B-hybrid): define substances
  // with the paper PK decay so ParaView can render the drug field. The kill
  // terms still use the well-mixed DrugState scalar — this is visualization
  // only. D is clamped below BioDynaMo's FTCS stability limits (which are hard
  // Fatal aborts) so the run never crashes regardless of dt / domain / res.
  if (sp->viz_drug_diffusion && (sp->treat_gem || sp->treat_abr)) {
    const real_t dt   = sp->dt_minutes;                       // diffusion Step dt
    const real_t dt_d = sp->dt_minutes / 1440.0;
    const real_t dx   = (sp->max_bound - sp->min_bound) /
                        std::max(1, sp->drug_grid_resolution);
    // Largest D that satisfies stability (μ+12D/dx²)·dt≤2 and decay-safety
    // 1−μ−6D·dt/dx²≥0, with an 0.8 safety factor. μ is the per-step decay.
    auto safe_D = [&](real_t mu) {
      const real_t d_stab  = (2.0 / dt - mu) * dx * dx / 12.0;
      const real_t d_decay = (1.0 - mu) * dx * dx / (6.0 * dt);
      const real_t d_max   = 0.8 * std::max(0.0, std::min(d_stab, d_decay));
      return std::min(static_cast<real_t>(sp->drug_diff_coeff), d_max);
    };
    if (sp->treat_gem) {
      const real_t mu = (1.0 - std::exp(-sp->gem_gamma * dt_d)) / dt;
      const real_t D  = safe_D(mu);
      ModelInitializer::DefineSubstance(kGemSubstance, "Gemcitabine", D, mu,
                                        sp->drug_grid_resolution);
      std::cout << "[viz] Gemcitabine diffusion grid: D=" << D
                << " mu=" << mu << " res=" << sp->drug_grid_resolution << "\n";
    }
    if (sp->treat_abr) {
      const real_t mu = (1.0 - std::exp(-sp->abr_gamma * dt_d)) / dt;
      const real_t D  = safe_D(mu);
      ModelInitializer::DefineSubstance(kAbrSubstance, "Abraxane", D, mu,
                                        sp->drug_grid_resolution);
      std::cout << "[viz] Abraxane diffusion grid: D=" << D
                << " mu=" << mu << " res=" << sp->drug_grid_resolution << "\n";
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

  // Viz sphere-seed: place EVERY cell inside a sphere so the mass renders as a
  // tumor spheroid instead of the domain box. Only for global-mode viz runs
  // (positions cosmetic). Otherwise use the normal placement.
  const real_t viz_radius = sp->viz_seed_radius_frac * (hi - lo) / 2.0;
  auto seed_pos = [&](bool is_tumor) -> Real3 {
    if (sp->viz_sphere_seed)
      return RandSpherePos(rng, domain_center, viz_radius, lo, hi);
    return is_tumor ? pos_tumor() : ClampPoint(rand_pos(), lo, hi);
  };

  // --- Seed populations ---
  for (size_t i = 0; i < sp->C0; ++i) {
    auto* c = new TumorCell(seed_pos(true));
    c->AddBehavior(new TumorBehavior());
    ctxt->AddAgent(c);
  }
  for (size_t i = 0; i < sp->P0; ++i) {
    auto* p = new StellateCell(seed_pos(true));
    p->AddBehavior(new PSCBehavior());
    ctxt->AddAgent(p);
  }
  for (size_t i = 0; i < sp->E0; ++i) {
    auto* e = new EffectorTCell(seed_pos(false));
    e->AddBehavior(new EffectorBehavior());
    ctxt->AddAgent(e);
  }
  for (size_t i = 0; i < sp->N0; ++i) {
    auto* n = new NKCell(seed_pos(false));
    n->AddBehavior(new NKBehavior());
    ctxt->AddAgent(n);
  }
  for (size_t i = 0; i < sp->H0; ++i) {
    auto* h = new HelperTCell(seed_pos(false));
    h->AddBehavior(new HelperBehavior());
    ctxt->AddAgent(h);
  }
  for (size_t i = 0; i < sp->R0; ++i) {
    auto* r = new TRegCell(seed_pos(false));
    r->AddBehavior(new TRegBehavior());
    ctxt->AddAgent(r);
  }
  // Cancer Stem Cells (Section 6) — DISABLED: this build reproduces only the
  // Section-5 treatment model. Commented out (not deleted) so CSC can be restored.
  // if (sp->csc_enable) {
  //   for (size_t i = 0; i < sp->S0; ++i) {
  //     auto* s = new CancerStemCell(pos_tumor());
  //     s->AddBehavior(new CSCBehavior());
  //     ctxt->AddAgent(s);
  //   }
  // }

  // Reporter + constant source (invisible agent; not using AlwaysCopyToNew)
  auto* rep = new ReporterCell();
  rep->AddBehavior(new ReportPopCounts());
  rep->AddBehavior(new SourceBehavior());
  if (sp->viz_drug_diffusion) rep->AddBehavior(new DrugDiffusionInjector());
  ctxt->AddAgent(rep);

  size_t total_steps =
      static_cast<size_t>(sp->total_days * 1440.0 / sp->dt_minutes);
  
  
  sim.GetScheduler()->Simulate(total_steps);

  std::cout << "Pancreatic tumor ABM completed ("
            << sp->total_days << " days).\n";
  return 0;
}

// ============================================================================
// FixDrugParaviewState
// Repairs the drug colour range in BioDynaMo's generated ParaView state so
// `bdm view` renders the diffusion field correctly — automatically, however the
// binary was launched (bdm run, direct, IDE). BioDynaMo colours the substance
// from frame 0 (drug = 0 before treatment) → a degenerate [0, FLT_MIN] range →
// flat/blank field. This runs scripts/paraview/fix_drug_pvsm.py (found relative
// to the executable) which reads the real data range and rescales ONLY the
// drug's transfer functions, leaving cells and colours untouched. No-ops when
// no drug field was exported. Call AFTER Simulate() returns (state is written
// when the Simulation is destroyed at the end of Simulate()).
// ============================================================================
inline void FixDrugParaviewState(
    const std::string& viz_dir = "output/pancreatic_tumor_new") {
  namespace fs = std::filesystem;
  std::error_code ec;
  if (!fs::exists(viz_dir, ec)) return;

  bool has_drug = false;
  for (const auto& e : fs::directory_iterator(viz_dir, ec)) {
    const std::string name = e.path().filename().string();
    if (e.path().extension() == ".pvti" &&
        (name.rfind("Abraxane-", 0) == 0 || name.rfind("Gemcitabine-", 0) == 0)) {
      has_drug = true;
      break;
    }
  }
  if (!has_drug) return;  // base / non-treatment run — no drug range to fix

  const char* pv = std::getenv("ParaView_DIR");
  if (pv == nullptr) {
    std::cerr << "[viz] ParaView_DIR unset — skipping pvsm fix "
                 "(run 'pvbatch scripts/paraview/fix_drug_pvsm.py " << viz_dir
              << "' manually).\n";
    return;
  }

  // Locate the fixer script relative to this executable (<repo>/build/<bin>).
  char buf[PATH_MAX];
  const ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
  if (n <= 0) return;
  buf[n] = '\0';
  const fs::path script =
      fs::path(buf).parent_path() / ".." / "scripts" / "paraview" / "fix_drug_pvsm.py";
  if (!fs::exists(script, ec)) {
    std::cerr << "[viz] fix_drug_pvsm.py not found next to binary — skipping.\n";
    return;
  }

  const std::string cmd = "\"" + std::string(pv) + "/bin/pvbatch\" \"" +
                          script.string() + "\" \"" + viz_dir + "\" 2>/dev/null";
  std::cout << "[viz] fixing drug-diffusion ParaView state for 'bdm view'...\n";
  std::system(cmd.c_str());  // pvbatch may segfault on exit AFTER writing — fine
}

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PANCREATIC_TUMOR_MODEL_H_
