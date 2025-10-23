// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey
// BioDynaMo collaboration. Apache-2.0 license.
//
// -----------------------------------------------------------------------------
#ifndef PANCREATIC_TUMOR_MODEL_H_
#define PANCREATIC_TUMOR_MODEL_H_

#include "biodynamo.h"
#include <cmath>
#include <algorithm>
#include <fstream>

namespace bdm {
namespace pancreatic_tumor {

// ============================================================================
// Parameters
// ============================================================================
struct Params {
  // --- time & space ---
  real_t dt_minutes      = 1.0;    // scheduler step (minutes)
  real_t min_bound       = 0.0;
  real_t max_bound       = 150.0;
  real_t cell_radius_um  = 6.0;

  // --- interaction mode ---
  bool   use_local_counts = true;  // false = global interactions, true = local
  real_t local_radius_um  = 12.0;   // neighborhood radius if local mode

  // --- initial counts ---
  size_t C0 = 474;  // Tumor (C)
  size_t P0 = 2;   // PSC (P)
  size_t E0 = 431;  // CD8+ Effector (E)
  size_t N0 = 189;  // NK (N)
  size_t H0 = 839;  // Helper (H)
  size_t R0 = 65;   // Tregs (R)

  // --- carrying capacities (used in crowding) ---
  // In global mode: these are global caps. In local mode: they are recomputed.
  real_t K_C = 7500.0;
  real_t K_P = 4500.0;
  real_t K_E = 3000.0;
  real_t K_N = 3000.0;
  real_t K_H = 4000.0;
  real_t K_R = 3500.0;

  // --- local-capacity estimation knobs (only used if use_local_counts=true) ---
  real_t packing_fraction = 0.64;   // random-close packing for spheres
  real_t capacity_margin  = 0.75;   // conservative margin
  bool   lock_equal_caps  = true;   // one K for all types
  size_t min_local_K      = 10;     // never go below this

  // Recompute Ks for local mode from radius and cell size
  void RecomputeCapsForLocal() {
    if (!use_local_counts) return;
    const real_t r  = std::max<real_t>(local_radius_um, cell_radius_um);
    const real_t R  = std::max<real_t>(cell_radius_um, 1.0);
    const real_t base = packing_fraction * capacity_margin * std::pow(r / R, 3.0);
    const real_t Kest = std::max<real_t>(std::round(base), static_cast<real_t>(min_local_K));
    if (lock_equal_caps) {
      K_C = K_P = K_E = K_N = K_R = K_H = Kest;
    } else {
      K_C = Kest;
      K_P = Kest;
      K_E = Kest;
      K_N = Kest;
      K_H = std::round(1.2 * Kest);
      K_R = std::round(1.1 * Kest);
    }
  }

  // --- smooth gates / half-saturation for global signals ---
  real_t gate_C_K = 500.0; // gate strength from tumor size

  // generic half-sats if needed elsewhere
  real_t K_small = 200.0;
  real_t K_med   = 600.0;
  real_t K_big   = 1200.0;

  // ================= Tumor (C) =================
  real_t c_base_div      = 0.10;   // baseline division
  real_t c_boost_from_P  = 0.25;   // P → C proliferative boost
  real_t c_boost_from_P_K= 800.0;  // global P half-sat for boost
  real_t c_kill_by_E     = 0.012;  // E→C killing
  real_t c_kill_by_E_K   = 600.0;  // global E half-sat
  real_t c_kill_by_N     = 0.007;  // N→C killing
  real_t c_kill_by_N_K   = 600.0;  // global N half-sat
  real_t c_R_blocks_E    = 0.10;   // 1/(1 + alpha * R) reduces E-kill

  // ================= PSC (P) =================
  real_t p_base_div      = 0.12;
  real_t p_boost_from_C  = 0.10;   // strong C → P boost
  real_t p_boost_from_C_K= 500.0; // later boost onset (global)
  real_t p_base_death    = 0.05;

  // ================= Effector T (E) =================
  real_t e_base_birth    = 0.035;
  real_t e_help_from_H   = 0.22;   // H → E help
  real_t e_help_from_H_K = 500.0;  // global H half-sat
  real_t e_inact_by_C    = 0.20;  // C → E inactivation
  real_t e_inact_by_C_K  = 400.0;  // global C half-sat
  real_t e_suppr_by_R    = 0.04;  // R → E suppression (gated by C)
  real_t e_suppr_by_R_K  = 500.0;  // global R half-sat
  real_t e_base_death    = 0.1;

  // ================= NK (N) =================
  real_t n_base_birth    = 0.03;
  real_t n_help_from_H   = 0.15;
  real_t n_help_from_H_K = 500.0;
  real_t n_inact_by_C    = 0.08;
  real_t n_inact_by_C_K  = 400.0;
  real_t n_suppr_by_R    = 0.038; // gated by C
  real_t n_suppr_by_R_K  = 500.0;
  real_t n_base_death    = 0.09;

  // ================= Helper T (H) =================
  real_t h_base_birth    = 0.05;
  real_t h_self_act      = 0.09;
  real_t h_self_act_K    = 1000.0; // weak self-activation, global
  real_t h_suppr_by_R    = 0.075; // gated by C
  real_t h_suppr_by_R_K  = 700.0;
  real_t h_base_death    = 0.08;

  // ================= Tregs (R) =================
  real_t r_base_src      = 0.08;
  real_t r_induced_by_E  = 0.008;  // E → R induction
  real_t r_induced_by_E_K= 600.0;
  real_t r_induced_by_H  = 0.008;  // H → R induction
  real_t r_induced_by_H_K= 600.0;
  real_t r_cleared_by_N  = 0.003;  // N → R clearance
  real_t r_cleared_by_N_K= 500.0;
  real_t r_decay         = 0.06;

  // --- Colors (UI) ---
  struct Colors {
    int tumor          = 8;
    int psc            = 4;
    int eff            = 5;
    int nk             = 3;
    int helper         = 6;
    int treg           = 1;
    int tumor_div_tint = 7;
  } color;
};

inline Params* P() { static Params p; return &p; }

// ============================================================================
// Helpers
// ============================================================================
inline real_t Clamp(real_t v, real_t lo, real_t hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}
inline Real3 ClampPoint(const Real3& pos) {
  return Real3{Clamp(pos[0], P()->min_bound, P()->max_bound),
               Clamp(pos[1], P()->min_bound, P()->max_bound),
               Clamp(pos[2], P()->min_bound, P()->max_bound)};
}
inline real_t Sat(real_t x, real_t K) { return (x <= 0) ? 0.0 : x / (K + x); }
inline real_t ProbFromRate(real_t rate_per_day, real_t dt_day) {
  if (rate_per_day <= 0) return 0.0;
  return 1.0 - std::exp(-rate_per_day * dt_day);
}

// ============================================================================
// Agent types
// ============================================================================
class TumorCell : public Cell {
  BDM_AGENT_HEADER(TumorCell, Cell, 1);
 public:
  TumorCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.tumor; }
  explicit TumorCell(const Real3& p) : TumorCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class StellateCell : public Cell {
  BDM_AGENT_HEADER(StellateCell, Cell, 1);
 public:
  StellateCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.psc; }
  explicit StellateCell(const Real3& p) : StellateCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class EffectorTCell : public Cell {
  BDM_AGENT_HEADER(EffectorTCell, Cell, 1);
 public:
  EffectorTCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.eff; }
  explicit EffectorTCell(const Real3& p) : EffectorTCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class NKCell : public Cell {
  BDM_AGENT_HEADER(NKCell, Cell, 1);
 public:
  NKCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.nk; }
  explicit NKCell(const Real3& p) : NKCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class HelperTCell : public Cell {
  BDM_AGENT_HEADER(HelperTCell, Cell, 1);
 public:
  HelperTCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.helper; }
  explicit HelperTCell(const Real3& p) : HelperTCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

class TRegCell : public Cell {
  BDM_AGENT_HEADER(TRegCell, Cell, 1);
 public:
  TRegCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.treg; }
  explicit TRegCell(const Real3& p) : TRegCell() { SetPosition(p); }
  void SetCellColor(int c) { color_ = c; }
  int  GetCellColor() const { return color_; }
 private:
  int color_ = 0;
};

// ============================================================================
// Counting utilities: GLOBAL and LOCAL
// ============================================================================
struct Counts { size_t C=0, P=0, E=0, N=0, H=0, R=0; };

// ---- Global census (computed at most once per step) ----
struct GlobalCensus {
  size_t step_cached = std::numeric_limits<size_t>::max();
  Counts cnt;

  static GlobalCensus& Instance() { static GlobalCensus gc; return gc; }

  void RefreshIfNeeded() {
    auto* sim = Simulation::GetActive();
    auto steps = sim->GetScheduler()->GetSimulatedSteps();
    if (steps == step_cached) return;

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
    step_cached = steps;
  }
};

// ---- Local neighborhood counter (within radius) ----
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

  static Counts Around(const Agent& a, real_t radius_um) {
    Counts out;
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    Fun f(&out);
    const real_t sr = radius_um * radius_um;
    ctxt->ForEachNeighbor(f, a, sr);
    return out;
  }
};

// Unified accessor: returns either global or local counts for the calling agent
inline Counts GetRelevantCounts(const Agent& self) {
  if (!P()->use_local_counts) {
    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();
    return gc.cnt;
  } else {
    return LocalNeighborhoodCounter::Around(self, P()->local_radius_um);
  }
}

// ============================================================================
// Behaviors
// ============================================================================
class TumorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TumorBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* c   = bdm_static_cast<TumorCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    c->SetPosition(ClampPoint(c->GetPosition()));

    const real_t dt_day = P()->dt_minutes / 1440.0;
    const Counts cnt = GetRelevantCounts(*c);

    // crowding on C
    real_t crowd = 1.0 - Clamp(static_cast<real_t>(cnt.C) / std::max<real_t>(1.0, P()->K_C), 0.0, 1.0);

    // division: baseline + PSC boost (saturated by PSC count)
    real_t boostP = P()->c_boost_from_P * Sat(static_cast<real_t>(cnt.P), P()->c_boost_from_P_K);
    real_t div_rate = (P()->c_base_div + boostP) * crowd;

    if (rng->Uniform(0,1) < ProbFromRate(div_rate, dt_day)) {
      auto* d = c->Divide();
      if (auto* cd = dynamic_cast<TumorCell*>(d)) {
        cd->SetCellColor(P()->color.tumor);
        cd->AddBehavior(new TumorBehavior());
      }
      c->SetCellColor(P()->color.tumor_div_tint);
    } else {
      c->SetCellColor(P()->color.tumor);
    }

    // death: E and N killing, E reduced by R
    // In global mode we gate by tumor size; in local mode, gating naturally comes from local counts.
    real_t gateC = P()->use_local_counts ? 1.0 : Sat(static_cast<real_t>(cnt.C), P()->gate_C_K);
    real_t inhibit_E = 1.0 / (1.0 + P()->c_R_blocks_E * static_cast<real_t>(cnt.R));
    real_t killE = P()->c_kill_by_E * Sat(static_cast<real_t>(cnt.E), P()->c_kill_by_E_K) * inhibit_E * gateC;
    real_t killN = P()->c_kill_by_N * Sat(static_cast<real_t>(cnt.N), P()->c_kill_by_N_K) * gateC;
    real_t p_kill = ProbFromRate(killE + killN, dt_day);

    if (rng->Uniform(0,1) < p_kill) {
      ctxt->RemoveAgent(c->GetUid());
    }
  }
};

class PSCBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(PSCBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* psc = bdm_static_cast<StellateCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    psc->SetPosition(ClampPoint(psc->GetPosition()));

    const real_t dt_day = P()->dt_minutes / 1440.0;
    const Counts cnt = GetRelevantCounts(*psc);

    real_t crowd  = 1.0 - Clamp(static_cast<real_t>(cnt.P) / std::max<real_t>(1.0, P()->K_P), 0.0, 1.0);
    real_t boostC = P()->p_boost_from_C * Sat(static_cast<real_t>(cnt.C), P()->p_boost_from_C_K);
    real_t div_rate = (P()->p_base_div + boostC) * crowd;
    real_t die_rate = P()->p_base_death;

    if (rng->Uniform(0,1) < ProbFromRate(die_rate, dt_day)) {
      ctxt->RemoveAgent(psc->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(div_rate, dt_day)) {
      auto* d = psc->Divide();
      if (auto* pd = dynamic_cast<StellateCell*>(d)) {
        pd->SetCellColor(P()->color.psc);
        pd->AddBehavior(new PSCBehavior());
      }
    }
  }
};

class EffectorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(EffectorBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* e   = bdm_static_cast<EffectorTCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    e->SetPosition(ClampPoint(e->GetPosition()));

    const real_t dt_day = P()->dt_minutes / 1440.0;
    const Counts cnt = GetRelevantCounts(*e);

    real_t crowdE = 1.0 - Clamp(static_cast<real_t>(cnt.E) / std::max<real_t>(1.0, P()->K_E), 0.0, 1.0);
    real_t helpH  = P()->e_help_from_H * Sat(static_cast<real_t>(cnt.H), P()->e_help_from_H_K);
    real_t birth  = (P()->e_base_birth + helpH) * crowdE;

    real_t gateC = P()->use_local_counts ? 1.0 : Sat(static_cast<real_t>(cnt.C), P()->gate_C_K);
    real_t die = P()->e_base_death
               + P()->e_inact_by_C * Sat(static_cast<real_t>(cnt.C), P()->e_inact_by_C_K) * gateC
               + P()->e_suppr_by_R * Sat(static_cast<real_t>(cnt.R), P()->e_suppr_by_R_K) * gateC;

    // small nonlinearity to avoid instant wipeout at tiny counts
    die *= static_cast<real_t>(cnt.E) / (cnt.E + 1.0);

    if (rng->Uniform(0,1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(e->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth, dt_day)) {
      auto* d = e->Divide();
      if (auto* ed = dynamic_cast<EffectorTCell*>(d)) {
        ed->SetCellColor(P()->color.eff);
        ed->AddBehavior(new EffectorBehavior());
      }
    }
  }
};

class NKBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(NKBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* n   = bdm_static_cast<NKCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    n->SetPosition(ClampPoint(n->GetPosition()));

    const real_t dt_day = P()->dt_minutes / 1440.0;
    const Counts cnt = GetRelevantCounts(*n);

    real_t crowdN = 1.0 - Clamp(static_cast<real_t>(cnt.N) / std::max<real_t>(1.0, P()->K_N), 0.0, 1.0);
    real_t helpH  = P()->n_help_from_H * Sat(static_cast<real_t>(cnt.H), P()->n_help_from_H_K);
    real_t birth  = (P()->n_base_birth + helpH) * crowdN;

    real_t gateC = P()->use_local_counts ? 1.0 : Sat(static_cast<real_t>(cnt.C), P()->gate_C_K);
    real_t die = P()->n_base_death
               + P()->n_inact_by_C * Sat(static_cast<real_t>(cnt.C), P()->n_inact_by_C_K) * gateC
               + P()->n_suppr_by_R * Sat(static_cast<real_t>(cnt.R), P()->n_suppr_by_R_K) * gateC;

    die *= static_cast<real_t>(cnt.N) / (cnt.N + 1.0);

    if (rng->Uniform(0,1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(n->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth, dt_day)) {
      auto* d = n->Divide();
      if (auto* nd = dynamic_cast<NKCell*>(d)) {
        nd->SetCellColor(P()->color.nk);
        nd->AddBehavior(new NKBehavior());
      }
    }
  }
};

class HelperBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(HelperBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* h   = bdm_static_cast<HelperTCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    h->SetPosition(ClampPoint(h->GetPosition()));

    const real_t dt_day = P()->dt_minutes / 1440.0;
    const Counts cnt = GetRelevantCounts(*h);

    real_t crowdH = 1.0 - Clamp(static_cast<real_t>(cnt.H) / std::max<real_t>(1.0, P()->K_H), 0.0, 1.0);
    real_t self   = P()->h_self_act * Sat(static_cast<real_t>(cnt.H), P()->h_self_act_K);
    real_t birth  = (P()->h_base_birth + self) * crowdH;

    real_t gateC = P()->use_local_counts ? 1.0 : Sat(static_cast<real_t>(cnt.C), P()->gate_C_K);
    real_t die = P()->h_base_death
               + P()->h_suppr_by_R * Sat(static_cast<real_t>(cnt.R), P()->h_suppr_by_R_K) * gateC;

    if (rng->Uniform(0,1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(h->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth, dt_day)) {
      auto* d = h->Divide();
      if (auto* hd = dynamic_cast<HelperTCell*>(d)) {
        hd->SetCellColor(P()->color.helper);
        hd->AddBehavior(new HelperBehavior());
      }
    }
  }
};

class TRegBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TRegBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* r   = bdm_static_cast<TRegCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    r->SetPosition(ClampPoint(r->GetPosition()));

    const real_t dt_day = P()->dt_minutes / 1440.0;
    const Counts cnt = GetRelevantCounts(*r);

    real_t crowdR = 1.0 - Clamp(static_cast<real_t>(cnt.R) / std::max<real_t>(1.0, P()->K_R), 0.0, 1.0);
    real_t birth  = (P()->r_base_src
                  +  P()->r_induced_by_E * Sat(static_cast<real_t>(cnt.E), P()->r_induced_by_E_K)
                  +  P()->r_induced_by_H * Sat(static_cast<real_t>(cnt.H), P()->r_induced_by_H_K))
                  *  crowdR;

    real_t die    = P()->r_decay
                  + P()->r_cleared_by_N * Sat(static_cast<real_t>(cnt.N), P()->r_cleared_by_N_K);

    if (rng->Uniform(0,1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(r->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth, dt_day)) {
      auto* d = r->Divide();
      if (auto* rd = dynamic_cast<TRegCell*>(d)) {
        rd->SetCellColor(P()->color.treg);
        rd->AddBehavior(new TRegBehavior());
      }
    }
  }
};

// ============================================================================
// Reporter (daily CSV)
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
    auto* rm  = sim->GetResourceManager();
    const auto steps = sim->GetScheduler()->GetSimulatedSteps();
    const real_t t_day = (P()->dt_minutes * steps) / 1440.0;

    size_t C=0, Pn=0, E=0, N=0, H=0, R=0;
    rm->ForEachAgent([&](Agent* a) {
      if      (dynamic_cast<TumorCell*>(a))     ++C;
      else if (dynamic_cast<StellateCell*>(a))  ++Pn;
      else if (dynamic_cast<EffectorTCell*>(a)) ++E;
      else if (dynamic_cast<NKCell*>(a))        ++N;
      else if (dynamic_cast<HelperTCell*>(a))   ++H;
      else if (dynamic_cast<TRegCell*>(a))      ++R;
    });

    static std::ofstream csv("data-export/populations.csv");
    static bool wrote_header = false;
    if (!wrote_header) {
      csv << "step,days,C,P,E,N,H,R,total\n";
      wrote_header = true;
    }
    // write once per day
    size_t steps_per_day = static_cast<size_t>(1440.0 / std::max<real_t>(1.0, P()->dt_minutes));
    if (steps % std::max<size_t>(steps_per_day, 1) == 0) {
      std::cout << "[day " << t_day << "] C=" << C << " P=" << Pn
                << " E=" << E << " N=" << N << " H=" << H << " R=" << R << "\n";
      csv << steps << "," << t_day + 7 << ","
          << C << "," << Pn << "," << E << "," << N << "," << H << "," << R << ","
          << (C+Pn+E+N+H+R) << "\n";
      csv.flush();
    }
  }
};

// ============================================================================
// Simulate
// ============================================================================
inline int Simulate(int argc, const char** argv) {
  auto setp = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = P()->min_bound;
    param->max_bound = P()->max_bound;
    param->simulation_time_step = P()->dt_minutes; // minutes per step
    // ensure neighbor search can cover local_radius_um safely
   
  };

  Simulation sim(argc, argv, setp);

  // If using local neighborhoods, convert carrying capacities to *local* caps
  if (P()->use_local_counts) {
    P()->RecomputeCapsForLocal();
  }

  auto* ctxt = sim.GetExecutionContext();
  auto* rng  = sim.GetRandom();

  auto rand_pos = [&]() -> Real3 {
    return Real3{rng->Uniform(P()->min_bound, P()->max_bound),
                 rng->Uniform(P()->min_bound, P()->max_bound),
                 rng->Uniform(P()->min_bound, P()->max_bound)};
  };

  // Seed populations
  for (size_t i=0; i<P()->C0; ++i) {
    auto* c = new TumorCell(ClampPoint(rand_pos()));
    c->SetCellColor(P()->color.tumor);
    c->AddBehavior(new TumorBehavior());
    ctxt->AddAgent(c);
  }
  for (size_t i=0; i<P()->P0; ++i) {
    auto* p = new StellateCell(ClampPoint(rand_pos()));
    p->SetCellColor(P()->color.psc);
    p->AddBehavior(new PSCBehavior());
    ctxt->AddAgent(p);
  }
  for (size_t i=0; i<P()->E0; ++i) {
    auto* e = new EffectorTCell(ClampPoint(rand_pos()));
    e->SetCellColor(P()->color.eff);
    e->AddBehavior(new EffectorBehavior());
    ctxt->AddAgent(e);
  }
  for (size_t i=0; i<P()->N0; ++i) {
    auto* n = new NKCell(ClampPoint(rand_pos()));
    n->SetCellColor(P()->color.nk);
    n->AddBehavior(new NKBehavior());
    ctxt->AddAgent(n);
  }
  for (size_t i=0; i<P()->H0; ++i) {
    auto* h = new HelperTCell(ClampPoint(rand_pos()));
    h->SetCellColor(P()->color.helper);
    h->AddBehavior(new HelperBehavior());
    ctxt->AddAgent(h);
  }
  for (size_t i=0; i<P()->R0; ++i) {
    auto* r = new TRegCell(ClampPoint(rand_pos()));
    r->SetCellColor(P()->color.treg);
    r->AddBehavior(new TRegBehavior());
    ctxt->AddAgent(r);
  }

  // Reporter agent
  auto* rep = new ReporterCell();
  rep->AddBehavior(new ReportPopCounts());
  ctxt->AddAgent(rep);

  // Run ~100 days (1 min step → 1440 steps per day)
  sim.GetScheduler()->Simulate(1440 * 93);
  std::cout << "Pancreatic tumor (C,P,E,N,H,R) ABM (global interactions) completed.\n";
  return 0;
}

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PANCREATIC_TUMOR_MODEL_H_
