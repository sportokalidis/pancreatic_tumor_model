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

// -------------------------------------------------------
// Parameters (per day unless noted)
// -------------------------------------------------------
struct Params {
  // Time & space
  real_t dt_minutes      = 1.0;   // scheduler step in minutes
  real_t min_bound       = 0.0;
  real_t max_bound       = 100.0;
  real_t cell_radius_um  = 6.0;
  real_t local_radius_um = 8.0;   // neighborhood radius used for local counts

  // Initial AGENT counts (small ABM; we keep numbers low and scale behaviorally)
  size_t C0 = 474;  // Tumor
  size_t P0 = 30;   // PSC
  size_t E0 = 431;  // CD8+
  size_t N0 = 189;  // NK
  size_t H0 = 839;  // Helper T
  size_t R0 = 65;   // Tregs

  // Local crowding caps (agent counts) — larger caps delay early plateaus
  int cap_C = 100;
  int cap_P = 800;
  int cap_E = 80;
  int cap_N = 80;
  int cap_H = 160;
  int cap_R = 80;

  // Global regulators for when interactions become strong
  // Scale local counts before saturation: x_eff = effect_scale * x
  real_t effect_scale = 1;   // < 1 delays saturation/negative feedback activation

  // Tumor gate so immune suppressions become strong only when C is large
  real_t gate_C_K     = 60.0;   // half-saturation in agents

  // Common half-saturation helpers for Sat(x,K) (in agents) - not used now!!!
  real_t K_small = 20.0;
  real_t K_med   = 45.0;
  real_t K_big   = 80.0;

  // ================= Tumor (C) =================
  real_t c_base_div      = 0.08;             // baseline division (drives early rise)
  real_t c_boost_from_P  = 0.25;             // P → C proliferative boost
  real_t c_boost_from_P_K= 50.0;             // later-onset boost
  real_t c_kill_by_E     = 0.012;            // E→C killing
  real_t c_kill_by_E_K   = 45.0;
  real_t c_kill_by_N     = 0.007;            // N→C killing
  real_t c_kill_by_N_K   = 45.0;
  real_t c_R_blocks_E    = 0.10;             // 1/(1 + alpha * R) reduces E effectiveness

  // ================= PSC (P) =================
  real_t p_base_div      = 0.62;
  real_t p_boost_from_C  = 0.90;             // strong C → P boost
  real_t p_boost_from_C_K= 70.0;             // later boost onset (inflection ~ day 40–50)
  real_t p_base_death    = 0.05;

  // ================= Effector T (E) =================
  real_t e_base_birth    = 0.035;
  real_t e_help_from_H   = 0.22;             // H → E help
  real_t e_help_from_H_K = 40.0;
  real_t e_inact_by_C    = 0.003;            // C → E inactivation
  real_t e_inact_by_C_K  = 55.0;
  real_t e_suppr_by_R    = 0.004;            // R → E suppression (will be gated by C)
  real_t e_suppr_by_R_K  = 35.0;
  real_t e_base_death    = 0.28;

  // ================= NK (N) =================
  real_t n_base_birth    = 0.03;
  real_t n_help_from_H   = 0.20;
  real_t n_help_from_H_K = 40.0;
  real_t n_inact_by_C    = 0.0028;
  real_t n_inact_by_C_K  = 55.0;
  real_t n_suppr_by_R    = 0.0038;           // gated by C
  real_t n_suppr_by_R_K  = 35.0;
  real_t n_base_death    = 0.09;

  // ================= Helper T (H) =================
  real_t h_base_birth    = 0.05;
  real_t h_self_act      = 0.09;
  real_t h_self_act_K    = 70.0;             // weak self-activation
  real_t h_suppr_by_R    = 0.0075;           // gated by C
  real_t h_suppr_by_R_K  = 45.0;
  real_t h_base_death    = 0.28;

  // ================= Tregs (R) =================
  real_t r_base_src      = 0.18;
  real_t r_induced_by_E  = 0.008;            // E → R induction (gentle, gradual rise)
  real_t r_induced_by_E_K= 40.0;
  real_t r_induced_by_H  = 0.008;            // H → R induction
  real_t r_induced_by_H_K= 40.0;
  real_t r_cleared_by_N  = 0.003;            // N → R clearance
  real_t r_cleared_by_N_K= 35.0;
  real_t r_decay         = 0.06;

  // Colors
  struct Colors {
    int tumor   = 8;
    int psc     = 4;
    int eff     = 5;
    int nk      = 3;
    int helper  = 6;
    int treg    = 1;
    int tumor_div_tint = 7;
  } color;
};

inline Params* P() { static Params p; return &p; }

// -------------------------------------------------------
// Helpers
// -------------------------------------------------------
inline real_t Clamp(real_t v, real_t lo, real_t hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}
inline Real3 ClampPoint(const Real3& pos) {
  return Real3{Clamp(pos[0], P()->min_bound, P()->max_bound),
               Clamp(pos[1], P()->min_bound, P()->max_bound),
               Clamp(pos[2], P()->min_bound, P()->max_bound)};
}
inline real_t Sat(real_t x, real_t K) { return (x <= 0) ? 0.0 : x / (K + x); }
// Saturation with a pre-scale of local counts: x_eff = effect_scale * x
inline real_t SatEff(real_t x, real_t K) { return Sat(P()->effect_scale * x, K); }

inline real_t ProbFromRate(real_t rate_per_day, real_t dt_day) {
  if (rate_per_day <= 0) return 0.0;
  real_t v = std::exp(-rate_per_day * dt_day);
  return 1.0 - v;
}

// -------------------------------------------------------
// Agents
// -------------------------------------------------------
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

// -------------------------------------------------------
// Neighborhood counting functor (agents, local)
// -------------------------------------------------------
struct Counts6 : public Functor<void, Agent*, real_t> {
  int C=0, Pn=0, E=0, N=0, H=0, R=0;
  void operator()(Agent* a, real_t) override {
    if      (dynamic_cast<TumorCell*>(a))     { ++C; }
    else if (dynamic_cast<StellateCell*>(a))  { ++Pn; }
    else if (dynamic_cast<EffectorTCell*>(a)) { ++E; }
    else if (dynamic_cast<NKCell*>(a))        { ++N; }
    else if (dynamic_cast<HelperTCell*>(a))   { ++H; }
    else if (dynamic_cast<TRegCell*>(a))      { ++R; }
  }
};

// -------------------------------------------------------
// Behaviors
// -------------------------------------------------------
class TumorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TumorBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* c   = bdm_static_cast<TumorCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    c->SetPosition(ClampPoint(c->GetPosition()));

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *c, r2);

    const real_t dt_day = P()->dt_minutes / 1440.0;

    // Birth: baseline + P boost, moderated by local C crowding
    real_t crowd = 1.0 - Clamp(static_cast<real_t>(cnt.C) / std::max(1, P()->cap_C), 0.0, 1.0);
    real_t boostP = P()->c_boost_from_P * SatEff(cnt.Pn, P()->c_boost_from_P_K);
    real_t div_rate = (P()->c_base_div + boostP) * crowd;

    if (rng->Uniform(0,1) < ProbFromRate(div_rate, dt_day)) {
      auto* d = c->Divide();
      if (auto* cd = dynamic_cast<TumorCell*>(d)) cd->SetCellColor(P()->color.tumor);
      c->SetCellColor(P()->color.tumor_div_tint);
    } else {
      c->SetCellColor(P()->color.tumor);
    }

    // Death: E and N killing; E effectiveness reduced by R and gated by tumor size
    real_t gateC = SatEff(cnt.C, P()->gate_C_K);  // delays strong immune effects until C is large
    real_t inhibit_E = 1.0 / (1.0 + P()->c_R_blocks_E * static_cast<real_t>(cnt.R));
    real_t killE = P()->c_kill_by_E * SatEff(cnt.E, P()->c_kill_by_E_K) * inhibit_E * gateC;
    real_t killN = P()->c_kill_by_N * SatEff(cnt.N, P()->c_kill_by_N_K) * gateC;
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

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *psc, r2);

    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t crowd  = 1.0 - Clamp(static_cast<real_t>(cnt.Pn) / std::max(1, P()->cap_P), 0.0, 1.0);
    real_t boostC = P()->p_boost_from_C * SatEff(cnt.C, P()->p_boost_from_C_K);
    real_t div_rate = (P()->p_base_div + boostC) * crowd;
    real_t die_rate = P()->p_base_death;

    if (rng->Uniform(0,1) < ProbFromRate(die_rate, dt_day)) {
      ctxt->RemoveAgent(psc->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(div_rate, dt_day)) {
      auto* d = psc->Divide();
      if (auto* pd = dynamic_cast<StellateCell*>(d)) pd->SetCellColor(P()->color.psc);
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

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *e, r2);

    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t helpH = P()->e_help_from_H * SatEff(cnt.H, P()->e_help_from_H_K);
    real_t birth = P()->e_base_birth + helpH;

    // Suppressions gated by tumor size so they strengthen later in time
    real_t gateC = SatEff(cnt.C, P()->gate_C_K);
    real_t die = P()->e_base_death
               + P()->e_inact_by_C * SatEff(cnt.C, P()->e_inact_by_C_K) * gateC
               + P()->e_suppr_by_R * SatEff(cnt.R, P()->e_suppr_by_R_K) * gateC;

    if (rng->Uniform(0,1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(e->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth, dt_day)) {
      auto* d = e->Divide();
      if (auto* ed = dynamic_cast<EffectorTCell*>(d)) ed->SetCellColor(P()->color.eff);
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

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *n, r2);

    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t helpH = P()->n_help_from_H * SatEff(cnt.H, P()->n_help_from_H_K);
    real_t birth = P()->n_base_birth + helpH;

    real_t gateC = SatEff(cnt.C, P()->gate_C_K);
    real_t die = P()->n_base_death
               + P()->n_inact_by_C * SatEff(cnt.C, P()->n_inact_by_C_K) * gateC
               + P()->n_suppr_by_R * SatEff(cnt.R, P()->n_suppr_by_R_K) * gateC;

    if (rng->Uniform(0,1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(n->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth, dt_day)) {
      auto* d = n->Divide();
      if (auto* nd = dynamic_cast<NKCell*>(d)) nd->SetCellColor(P()->color.nk);
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

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *h, r2);

    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t self  = P()->h_self_act * SatEff(cnt.H, P()->h_self_act_K);
    real_t birth = P()->h_base_birth + self;

    real_t gateC = SatEff(cnt.C, P()->gate_C_K);
    real_t die = P()->h_base_death
               + P()->h_suppr_by_R * SatEff(cnt.R, P()->h_suppr_by_R_K) * gateC;

    if (rng->Uniform(0,1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(h->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth, dt_day)) {
      auto* d = h->Divide();
      if (auto* hd = dynamic_cast<HelperTCell*>(d)) hd->SetCellColor(P()->color.helper);
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

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *r, r2);

    const real_t dt_day = P()->dt_minutes / 1440.0;

    // R rise gradually: baseline source + induction by E and H (saturated)
    real_t birth = P()->r_base_src
                 + P()->r_induced_by_E * SatEff(cnt.E, P()->r_induced_by_E_K)
                 + P()->r_induced_by_H * SatEff(cnt.H, P()->r_induced_by_H_K);

    real_t die   = P()->r_decay
                 + P()->r_cleared_by_N * SatEff(cnt.N, P()->r_cleared_by_N_K);

    if (rng->Uniform(0,1) < ProbFromRate(die, dt_day)) {
      ctxt->RemoveAgent(r->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth, dt_day)) {
      auto* d = r->Divide();
      if (auto* rd = dynamic_cast<TRegCell*>(d)) rd->SetCellColor(P()->color.treg);
    }
  }
};

// -------------------------------------------------------
// Reporter (daily CSV)
// -------------------------------------------------------
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
    if (steps % static_cast<size_t>(1440.0 / std::max<real_t>(1.0, P()->dt_minutes)) == 0) {
      std::cout << "[day " << t_day << "] C=" << C << " P=" << Pn
                << " E=" << E << " N=" << N << " H=" << H << " R=" << R << "\n";
      csv << steps << "," << t_day << ","
          << C << "," << Pn << "," << E << "," << N << "," << H << "," << R << ","
          << (C+Pn+E+N+H+R) << "\n";
      csv.flush();
    }
  }
};

// -------------------------------------------------------
// Simulate
// -------------------------------------------------------
inline int Simulate(int argc, const char** argv) {
  auto setp = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = P()->min_bound;
    param->max_bound = P()->max_bound;
    param->simulation_time_step = P()->dt_minutes; // minutes
  };

  Simulation sim(argc, argv, setp);
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

  // Reporter
  auto* rep = new ReporterCell();
  rep->AddBehavior(new ReportPopCounts());
  sim.GetExecutionContext()->AddAgent(rep);

  // Run ~100 days (1 min step → 1440 steps per day)
  sim.GetScheduler()->Simulate(1440 * 100);
  std::cout << "Pancreatic tumor (C,P,E,N,H,R) simple ABM completed.\n";
  return 0;
}

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PANCREATIC_TUMOR_MODEL_H_
