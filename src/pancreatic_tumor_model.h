// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef PANCREATIC_TUMOR_MODEL_H_
#define PANCREATIC_TUMOR_MODEL_H_

#include "biodynamo.h"
#include <cmath>
#include <algorithm>   // std::min, std::max
#include <fstream>

namespace bdm {
namespace pancreatic_tumor {

// -------------------------------------------------------
// Parameters (per day unless noted)
// All half-sats & caps are in REAL CELLS (not agents).
// Each agent carries a "weight" = number of real cells.
// -------------------------------------------------------
struct Params {
  // Time & space
  real_t dt_minutes      = 1.0;     // scheduler step (minutes)
  real_t min_bound       = 0.0;
  real_t max_bound       = 100.0;   // larger box for stable densities
  real_t cell_radius_um  = 6.0;
  real_t local_radius_um = 8.0;     // neighborhood radius

  // ---------- Representative-agent weights (real cells per agent)
  // Calibrated from your CSV day-0 logs (log10 -> cells / agents).
  double wC = 9.501e6;
  double wP = 1.096e7;
  double wE = 3.349e4;
  double wN = 1.139e6;
  double wH = 2.566e5;
  double wR = 1.686e6;

  // ---------- Initial AGENT counts (not real cells)
  size_t C0 = 474;  // Tumor
  size_t P0 = 10;   // PSC
  size_t E0 = 431;  // CD8+
  size_t N0 = 189;  // NK
  size_t H0 = 839;  // Helper T
  size_t R0 = 65;   // Tregs

  // ---------- Weight housekeeping (keep agent counts bounded)
  // Agents grow by weight, not by division. Clamp weights into [w_min, w_max].
  double w_min_cells = 5e3;   // remove agent if it shrinks below this
  double w_max_cells = 2e8;   // cap per-agent weight to avoid huge representatives

  // ---------- Soft local carrying capacities (crowding; REAL CELLS)
  double cap_C_cells = 1e10;
  double cap_P_cells = 5e9;
  double cap_R_cells = 2e9;   // NEW: R logistic cap to avoid runaway induction

  // ---------- Half-saturation constants (REAL CELLS)
  double c_psc_half_sat_cells   = 5e7;   // PSC → C boost
  double p_tumor_half_sat_cells = 5e7;   // C → PSC boost
  double e_help_half_sat_cells  = 3e7;   // H → E
  double n_help_half_sat_cells  = 3e7;   // H → N
  double h_self_half_sat_cells  = 4e7;   // H → H

  // ---------------- Tumor (C) ----------------
  real_t c_base_div           = 0.06;     // ↑ slightly for faster rise
  real_t c_psc_boost          = 0.30;
  real_t c_kill_by_E          = 6.0e-3;   // ↑ mild late bending
  real_t c_kill_by_N          = 3.0e-3;
  real_t c_r_inhib_Ekill      = 0.10;

  // ---------------- PSC (P) ----------------
  real_t p_base_div           = 0.3;
  real_t p_tumor_boost        = 2.0;     // ↑ stronger C→P drive
  real_t p_base_death         = 0.01;     // ↓ slower PSC loss

  // ---------------- Effector T (E) ----------------
  real_t e_base_birth         = 0.03;
  real_t e_help_from_H        = 0.30;
  real_t e_inactivation_by_C  = 2.0e-3;   // ↑ modestly
  real_t e_suppression_by_R   = 3.0e-3;
  real_t e_base_death         = 0.05;     // ↑ smoother but faster decline

  // ---------------- NK (N) ----------------
  real_t n_base_birth         = 0.03;
  real_t n_help_from_H        = 0.25;
  real_t n_inactivation_by_C  = 2.0e-3;   // ↑
  real_t n_suppression_by_R   = 4.0e-3;
  real_t n_base_death         = 0.05;     // ↑

  // ---------------- Helper T (H) ----------------
  real_t h_base_birth         = 0.05;
  real_t h_self_activation    = 0.10;
  real_t h_suppression_by_R   = 5.0e-3;   // ↑ stronger R→H suppression
  real_t h_base_death         = 0.10;     // ↑ faster H decay (matches paper)

  // ---------------- Tregs (R) ----------------
  // Moderated vs previous to avoid agent explosion; also capped by cap_R_cells.
  real_t r_base_source        = 0.30;
  real_t r_induction_by_E     = 5.0e-3;
  real_t r_induction_by_H     = 5.0e-3;
  real_t r_decay              = 0.08;     // ↑ a bit vs 0.05
  real_t r_clear_by_N         = 2.0e-3;   // ↑ clearance

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
inline double Clamp(double v, double lo, double hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}
inline Real3 ClampPoint(const Real3& pos) {
  return Real3{Clamp(pos[0], P()->min_bound, P()->max_bound),
               Clamp(pos[1], P()->min_bound, P()->max_bound),
               Clamp(pos[2], P()->min_bound, P()->max_bound)};
}
inline double Sat(double x, double K) {   // simple Hill-like saturation x/(K+x)
  return (x <= 0.0) ? 0.0 : x / (K + x);
}
inline double ProbFromRate(double rate_per_day, double dt_day) {
  if (rate_per_day <= 0.0) return 0.0;
  double v = std::exp(-rate_per_day * dt_day);
  return 1.0 - v; // in [0,1)
}

// --- Weight update (no agent division). Grows/shrinks weights smoothly.
//     If too small => remove agent. If too big => clamp to w_max.
template <typename TCell>
inline bool UpdateWeightOrRemove(TCell* a, double birth_rate, double death_rate) {
  auto* sim  = Simulation::GetActive();
  auto* ctxt = sim->GetExecutionContext();
  const double dt_day = P()->dt_minutes / 1440.0;

  double net = birth_rate - death_rate;      // per day
  double w   = a->GetWeight();               // real cells in this agent
  w *= std::exp(net * dt_day);               // exponential Euler (stable)

  if (w < P()->w_min_cells) {
    ctxt->RemoveAgent(a->GetUid());
    return true; // removed
  }
  if (w > P()->w_max_cells) {
    w = P()->w_max_cells;                    // clamp instead of spawning
  }
  a->SetWeight(w);
  return false;
}

// -------------------------------------------------------
// Agents with weights (real cells represented per agent)
// -------------------------------------------------------
class TumorCell : public Cell {
  BDM_AGENT_HEADER(TumorCell, Cell, 1);
 public:
  TumorCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.tumor; }
  explicit TumorCell(const Real3& p) : TumorCell() { SetPosition(p); }
  void   SetCellColor(int c) { color_ = c; }
  int    GetCellColor() const { return color_; }
  void   SetWeight(double w) { weight_ = w; }
  double GetWeight() const { return weight_; }
 private:
  int    color_  = 0;
  double weight_ = 1.0;
};

class StellateCell : public Cell {
  BDM_AGENT_HEADER(StellateCell, Cell, 1);
 public:
  StellateCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.psc; }
  explicit StellateCell(const Real3& p) : StellateCell() { SetPosition(p); }
  void   SetCellColor(int c) { color_ = c; }
  int    GetCellColor() const { return color_; }
  void   SetWeight(double w) { weight_ = w; }
  double GetWeight() const { return weight_; }
 private:
  int    color_  = 0;
  double weight_ = 1.0;
};

class EffectorTCell : public Cell {
  BDM_AGENT_HEADER(EffectorTCell, Cell, 1);
 public:
  EffectorTCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.eff; }
  explicit EffectorTCell(const Real3& p) : EffectorTCell() { SetPosition(p); }
  void   SetCellColor(int c) { color_ = c; }
  int    GetCellColor() const { return color_; }
  void   SetWeight(double w) { weight_ = w; }
  double GetWeight() const { return weight_; }
 private:
  int    color_  = 0;
  double weight_ = 1.0;
};

class NKCell : public Cell {
  BDM_AGENT_HEADER(NKCell, Cell, 1);
 public:
  NKCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.nk; }
  explicit NKCell(const Real3& p) : NKCell() { SetPosition(p); }
  void   SetCellColor(int c) { color_ = c; }
  int    GetCellColor() const { return color_; }
  void   SetWeight(double w) { weight_ = w; }
  double GetWeight() const { return weight_; }
 private:
  int    color_  = 0;
  double weight_ = 1.0;
};

class HelperTCell : public Cell {
  BDM_AGENT_HEADER(HelperTCell, Cell, 1);
 public:
  HelperTCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.helper; }
  explicit HelperTCell(const Real3& p) : HelperTCell() { SetPosition(p); }
  void   SetCellColor(int c) { color_ = c; }
  int    GetCellColor() const { return color_; }
  void   SetWeight(double w) { weight_ = w; }
  double GetWeight() const { return weight_; }
 private:
  int    color_  = 0;
  double weight_ = 1.0;
};

class TRegCell : public Cell {
  BDM_AGENT_HEADER(TRegCell, Cell, 1);
 public:
  TRegCell() { SetDiameter(2.0 * P()->cell_radius_um); color_ = P()->color.treg; }
  explicit TRegCell(const Real3& p) : TRegCell() { SetPosition(p); }
  void   SetCellColor(int c) { color_ = c; }
  int    GetCellColor() const { return color_; }
  void   SetWeight(double w) { weight_ = w; }
  double GetWeight() const { return weight_; }
 private:
  int    color_  = 0;
  double weight_ = 1.0;
};

// -------------------------------------------------------
// Weighted neighborhood counting functor (sums real cells)
// -------------------------------------------------------
struct WCounts6 : public Functor<void, Agent*, real_t> {
  double C=0, Pn=0, E=0, N=0, H=0, R=0;
  void operator()(Agent* a, real_t) override {
    if      (auto* t  = dynamic_cast<TumorCell*>(a))     { C  += t->GetWeight(); return; }
    else if (auto* ps = dynamic_cast<StellateCell*>(a))  { Pn += ps->GetWeight(); return; }
    else if (auto* e  = dynamic_cast<EffectorTCell*>(a)) { E  += e->GetWeight(); return; }
    else if (auto* n  = dynamic_cast<NKCell*>(a))        { N  += n->GetWeight(); return; }
    else if (auto* h  = dynamic_cast<HelperTCell*>(a))   { H  += h->GetWeight(); return; }
    else if (auto* r  = dynamic_cast<TRegCell*>(a))      { R  += r->GetWeight(); return; }
  }
};

// -------------------------------------------------------
// Behaviors (using REAL-CELL neighbor counts, weight updates)
// -------------------------------------------------------
class TumorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TumorBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* c   = bdm_static_cast<TumorCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    c->SetPosition(ClampPoint(c->GetPosition()));

    WCounts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *c, r2);

    double Fp    = Sat(cnt.Pn, P()->c_psc_half_sat_cells);
    double crowd = 1.0 - Clamp(cnt.C / std::max(1.0, P()->cap_C_cells), 0.0, 1.0);
    double birth = (P()->c_base_div + P()->c_psc_boost * Fp) * crowd;

    double inhibit_E = 1.0 / (1.0 + P()->c_r_inhib_Ekill * cnt.R);
    double death     = P()->c_kill_by_E * cnt.E * inhibit_E
                     + P()->c_kill_by_N * cnt.N;

    UpdateWeightOrRemove(c, birth, death);
  }
};

class PSCBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(PSCBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* psc = bdm_static_cast<StellateCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    psc->SetPosition(ClampPoint(psc->GetPosition()));

    WCounts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *psc, r2);

    double Fc    = Sat(cnt.C, P()->p_tumor_half_sat_cells);
    double crowd = 1.0 - Clamp(cnt.Pn / std::max(1.0, P()->cap_P_cells), 0.0, 1.0);

    double birth = (P()->p_base_div + P()->p_tumor_boost * Fc) * crowd;
    double death = P()->p_base_death;

    UpdateWeightOrRemove(psc, birth, death);
  }
};

class EffectorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(EffectorBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* e   = bdm_static_cast<EffectorTCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    e->SetPosition(ClampPoint(e->GetPosition()));

    WCounts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *e, r2);

    double help  = P()->e_help_from_H * Sat(cnt.H, P()->e_help_half_sat_cells);
    double birth = P()->e_base_birth + help;
    double death = P()->e_base_death
                 + P()->e_inactivation_by_C * cnt.C
                 + P()->e_suppression_by_R  * cnt.R;

    UpdateWeightOrRemove(e, birth, death);
  }
};

class NKBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(NKBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* n   = bdm_static_cast<NKCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    n->SetPosition(ClampPoint(n->GetPosition()));

    WCounts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *n, r2);

    double help  = P()->n_help_from_H * Sat(cnt.H, P()->n_help_half_sat_cells);
    double birth = P()->n_base_birth + help;
    double death = P()->n_base_death
                 + P()->n_inactivation_by_C * cnt.C
                 + P()->n_suppression_by_R  * cnt.R;

    UpdateWeightOrRemove(n, birth, death);
  }
};

class HelperBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(HelperBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* h   = bdm_static_cast<HelperTCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    h->SetPosition(ClampPoint(h->GetPosition()));

    WCounts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *h, r2);

    double self  = P()->h_self_activation * Sat(cnt.H, P()->h_self_half_sat_cells);
    double birth = P()->h_base_birth + self;
    double death = P()->h_base_death
                 + P()->h_suppression_by_R * cnt.R;

    UpdateWeightOrRemove(h, birth, death);
  }
};

class TRegBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TRegBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* r   = bdm_static_cast<TRegCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    r->SetPosition(ClampPoint(r->GetPosition()));

    WCounts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *r, r2);

    // Birth from basal + induction by E and H (REAL counts), with logistic cap
    double birth = P()->r_base_source
                 + P()->r_induction_by_E * cnt.E
                 + P()->r_induction_by_H * cnt.H;

    double crowdR = 1.0 - Clamp(cnt.R / std::max(1.0, P()->cap_R_cells), 0.0, 1.0);
    birth *= crowdR;

    // Death from decay + NK clearance
    double death = P()->r_decay
                 + P()->r_clear_by_N * cnt.N;

    UpdateWeightOrRemove(r, birth, death);
  }
};


// -------------------------------------------------------
// Reporter: logs agent counts and REAL-CELL totals each day
// -------------------------------------------------------
class ReporterCell : public Cell {
  BDM_AGENT_HEADER(ReporterCell, Cell, 1);
 public:
  ReporterCell() { SetDiameter(0.1); }  // tiny, non-intrusive
};

class ReportPopCounts : public Behavior {
  BDM_BEHAVIOR_HEADER(ReportPopCounts, Behavior, 1);
 public:
  void Run(Agent* /*unused*/) override {
    auto* sim = Simulation::GetActive();
    auto* rm  = sim->GetResourceManager();
    const auto steps = sim->GetScheduler()->GetSimulatedSteps();
    const double t_day = (P()->dt_minutes * steps) / 1440.0;

    // Count by AGENT type and sum REAL-CELL totals
    size_t C=0, Pn=0, E=0, N=0, H=0, R=0;
    double RC=0, RP=0, RE=0, RN=0, RH=0, RR=0;

    rm->ForEachAgent([&](Agent* a) {
      if      (auto* x = dynamic_cast<TumorCell*>(a))     { ++C;  RC += x->GetWeight(); }
      else if (auto* x = dynamic_cast<StellateCell*>(a))  { ++Pn; RP += x->GetWeight(); }
      else if (auto* x = dynamic_cast<EffectorTCell*>(a)) { ++E;  RE += x->GetWeight(); }
      else if (auto* x = dynamic_cast<NKCell*>(a))        { ++N;  RN += x->GetWeight(); }
      else if (auto* x = dynamic_cast<HelperTCell*>(a))   { ++H;  RH += x->GetWeight(); }
      else if (auto* x = dynamic_cast<TRegCell*>(a))      { ++R;  RR += x->GetWeight(); }
    });

    static std::ofstream csv("data-export/populations.csv");
    static bool wrote_header = false;
    if (!wrote_header) {
      csv << "step,days,C,P,E,N,H,R,total\n";
      wrote_header = true;
    }

    // Optional: print daily + log to CSV daily
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
// Simulate: creates agents, attaches behaviors, runs
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

  // Seed populations (assign weights)
  for (size_t i=0; i<P()->C0; ++i) {
    auto* c = new TumorCell(ClampPoint(rand_pos()));
    c->SetCellColor(P()->color.tumor);
    c->SetWeight(P()->wC);
    c->AddBehavior(new TumorBehavior());
    ctxt->AddAgent(c);
  }
  for (size_t i=0; i<P()->P0; ++i) {
    auto* p = new StellateCell(ClampPoint(rand_pos()));
    p->SetCellColor(P()->color.psc);
    p->SetWeight(P()->wP);
    p->AddBehavior(new PSCBehavior());
    ctxt->AddAgent(p);
  }
  for (size_t i=0; i<P()->E0; ++i) {
    auto* e = new EffectorTCell(ClampPoint(rand_pos()));
    e->SetCellColor(P()->color.eff);
    e->SetWeight(P()->wE);
    e->AddBehavior(new EffectorBehavior());
    ctxt->AddAgent(e);
  }
  for (size_t i=0; i<P()->N0; ++i) {
    auto* n = new NKCell(ClampPoint(rand_pos()));
    n->SetCellColor(P()->color.nk);
    n->SetWeight(P()->wN);
    n->AddBehavior(new NKBehavior());
    ctxt->AddAgent(n);
  }
  for (size_t i=0; i<P()->H0; ++i) {
    auto* h = new HelperTCell(ClampPoint(rand_pos()));
    h->SetCellColor(P()->color.helper);
    h->SetWeight(P()->wH);
    h->AddBehavior(new HelperBehavior());
    ctxt->AddAgent(h);
  }
  for (size_t i=0; i<P()->R0; ++i) {
    auto* r = new TRegCell(ClampPoint(rand_pos()));
    r->SetCellColor(P()->color.treg);
    r->SetWeight(P()->wR);
    r->AddBehavior(new TRegBehavior());
    ctxt->AddAgent(r);
  }

  // Add one reporter agent that logs every step (daily to CSV)
  auto* rep = new ReporterCell();
  rep->AddBehavior(new ReportPopCounts());
  sim.GetExecutionContext()->AddAgent(rep);

  // Run (~100 days at dt=1 minute)
  sim.GetScheduler()->Simulate(1440 * 100);
  std::cout << "Pancreatic tumor (C,P,E,N,H,R) weighted ABM completed.\n";
  return 0;
}

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PANCREATIC_TUMOR_MODEL_H_
