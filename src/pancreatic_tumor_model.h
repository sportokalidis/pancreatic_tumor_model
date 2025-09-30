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
#include <algorithm>  // std::max, std::min

namespace bdm {
namespace pancreatic_tumor {

// -------------------------------------------------------
// Minimal, calibratable parameters (units: per day unless noted)
// -------------------------------------------------------
struct Params {
  // Time & space
  real_t dt_minutes      = 1.0;     // scheduler step = 1 minute by default
  real_t min_bound       = 0.0;
  real_t max_bound       = 100.0;
  real_t cell_radius_um  = 6.0;
  real_t local_radius_um = 6.0;     // neighborhood radius for counting

  // Initial counts (choose ratios to match figure; absolute scale is up to you)
  size_t C0 = 600;   // Tumor (PCC)
  size_t P0 = 240;   // Stellate (PSC)
  size_t E0 = 120;    // CD8+ T
  size_t N0 = 80;    // NK
  size_t H0 = 180;    // Helper T
  size_t R0 = 40;    // Tregs

  // Soft local carrying capacities (crowding caps) for division moderation
  int cap_C = 30;
  int cap_P = 20;
  int cap_E = 25;
  int cap_N = 25;
  int cap_H = 40;
  int cap_R = 30;

  // ---------------- Tumor (C) ----------------
  real_t c_base_div           = 0.08;     // basal C division
  real_t c_psc_boost          = 0.60;     // PSC -> C (strength)
  real_t c_psc_half_sat       = 8.0;      // PSC half-saturation for boost
  real_t c_kill_by_E          = 1.0e-3;   // E -> kill C (per E)
  real_t c_kill_by_N          = 8.0e-4;   // N -> kill C (per N)
  real_t c_r_inhib_Ekill      = 0.08;     // R inhibits E-kill; factor = 1/(1 + r1 * Rloc)

  // ---------------- PSC (P) ----------------
  real_t p_base_div           = 0.18;     // basal P division
  real_t p_tumor_boost        = 0.40;     // C -> P (strength)
  real_t p_tumor_half_sat     = 10.0;     // C half-saturation
  real_t p_base_death         = 0.15;     // P death

  // ---------------- Effector T (E) ----------------
  real_t e_base_birth         = 0.05;     // basal E birth / upreg
  real_t e_help_from_H        = 0.25;     // H -> E (strength)
  real_t e_help_half_sat      = 10.0;     // H half-saturation
  real_t e_inactivation_by_C  = 2.0e-3;   // C -> E death / inactivation (per C)
  real_t e_suppression_by_R   = 4.0e-3;   // R -> E death (per R)
  real_t e_base_death         = 0.10;     // baseline E death

  // ---------------- NK (N) ----------------
  real_t n_base_birth         = 0.04;     // basal N birth
  real_t n_help_from_H        = 0.20;     // H -> N (strength)
  real_t n_help_half_sat      = 10.0;     // H half-saturation
  real_t n_inactivation_by_C  = 1.5e-3;   // C -> N death (per C)
  real_t n_suppression_by_R   = 3.0e-3;   // R -> N death (per R)
  real_t n_base_death         = 0.10;     // baseline N death

  // ---------------- Helper T (H) ----------------
  real_t h_base_birth         = 0.06;     // basal H birth
  real_t h_self_activation    = 0.20;     // H -> H (auto-activation)
  real_t h_self_half_sat      = 15.0;     // half-saturation
  real_t h_suppression_by_R   = 2.0e-3;   // R -> H death (per R)
  real_t h_base_death         = 0.08;     // baseline H death

  // ---------------- Tregs (R) ----------------
  real_t r_base_source        = 0.02;     // basal R birth / influx
  real_t r_induction_by_E     = 1.5e-3;   // E -> R (per E)
  real_t r_induction_by_H     = 1.5e-3;   // H -> R (per H)
  real_t r_decay              = 0.12;     // baseline R death
  real_t r_clear_by_N         = 2.0e-3;   // N -> R death (per N)

  // Colors (ParaView integer ids, easy to retune)
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
inline real_t Sat(real_t x, real_t K) {   // simple Hill-like saturation x/(K+x)
  return (x <= 0) ? 0.0 : x / (K + x);
}
inline real_t ProbFromRate(real_t rate_per_day, real_t dt_day) {
  // hazard mapping: stable for small or large dt
  if (rate_per_day <= 0) return 0.0;
  real_t v = std::exp(-rate_per_day * dt_day);
  return 1.0 - v; // in [0,1)
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
// Neighborhood counting functor (cling/ROOT friendly)
// -------------------------------------------------------
struct Counts6 : public Functor<void, Agent*, real_t> {
  int C=0, Pn=0, E=0, N=0, H=0, R=0;
  void operator()(Agent* a, real_t) override {
    if      (dynamic_cast<TumorCell*>(a))     { ++C; return; }
    else if (dynamic_cast<StellateCell*>(a))  { ++Pn; return; }
    else if (dynamic_cast<EffectorTCell*>(a)) { ++E; return; }
    else if (dynamic_cast<NKCell*>(a))        { ++N; return; }
    else if (dynamic_cast<HelperTCell*>(a))   { ++H; return; }
    else if (dynamic_cast<TRegCell*>(a))      { ++R; return; }
  }
};

// -------------------------------------------------------
// Behaviors (simple, calibratable interaction sizes)
// -------------------------------------------------------
class TumorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TumorBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* c   = bdm_static_cast<TumorCell*>(a);
    auto* sim = Simulation::GetActive();
    auto* ctxt= sim->GetExecutionContext();
    auto* rng = sim->GetRandom();
    c->SetPosition(ClampPoint(c->GetPosition()));

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *c, r2);

    const real_t dt_day = P()->dt_minutes / 1440.0;

    // Division: (baseline + PSC boost) * C-crowding
    real_t Fp = Sat(cnt.Pn, P()->c_psc_half_sat);
    real_t crowd = 1.0 - Clamp(static_cast<real_t>(cnt.C) / std::max(1, P()->cap_C), 0.0, 1.0);
    real_t div_rate = (P()->c_base_div + P()->c_psc_boost * Fp) * crowd;
    real_t p_div = ProbFromRate(div_rate, dt_day);

    if (rng->Uniform(0,1) < p_div) {
      auto* d = c->Divide();
      if (auto* cd = dynamic_cast<TumorCell*>(d)) cd->SetCellColor(P()->color.tumor);
      c->SetCellColor(P()->color.tumor_div_tint);
    } else {
      c->SetCellColor(P()->color.tumor);
    }

    // Killing by E and N; E effectiveness reduced by Tregs R
    real_t inhibit_E = 1.0 / (1.0 + P()->c_r_inhib_Ekill * static_cast<real_t>(cnt.R));
    real_t kill_rate = P()->c_kill_by_E * static_cast<real_t>(cnt.E) * inhibit_E
                     + P()->c_kill_by_N * static_cast<real_t>(cnt.N);
    real_t p_kill = ProbFromRate(kill_rate, dt_day);

    if (rng->Uniform(0,1) < p_kill) {
      sim->GetExecutionContext()->RemoveAgent(c->GetUid());
    }
  }
};

class PSCBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(PSCBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* psc = bdm_static_cast<StellateCell*>(a);
    auto* sim = Simulation::GetActive();
    auto* ctxt= sim->GetExecutionContext();
    auto* rng = sim->GetRandom();
    psc->SetPosition(ClampPoint(psc->GetPosition()));

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *psc, r2);

    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t Fc = Sat(cnt.C, P()->p_tumor_half_sat);
    real_t crowd = 1.0 - Clamp(static_cast<real_t>(cnt.Pn) / std::max(1, P()->cap_P), 0.0, 1.0);

    real_t div_rate = (P()->p_base_div + P()->p_tumor_boost * Fc) * crowd;
    real_t die_rate = P()->p_base_death;

    // Death first (avoids immediate replace)
    if (rng->Uniform(0,1) < ProbFromRate(die_rate, dt_day)) {
      sim->GetExecutionContext()->RemoveAgent(psc->GetUid());
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
    auto* sim = Simulation::GetActive();
    auto* ctxt= sim->GetExecutionContext();
    auto* rng = sim->GetRandom();
    e->SetPosition(ClampPoint(e->GetPosition()));

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *e, r2);
    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t help = P()->e_help_from_H * Sat(cnt.H, P()->e_help_half_sat);
    real_t birth_rate = P()->e_base_birth + help;
    real_t death_rate = P()->e_base_death
                      + P()->e_inactivation_by_C * static_cast<real_t>(cnt.C)
                      + P()->e_suppression_by_R  * static_cast<real_t>(cnt.R);

    // Death first
    if (rng->Uniform(0,1) < ProbFromRate(death_rate, dt_day)) {
      sim->GetExecutionContext()->RemoveAgent(e->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth_rate, dt_day)) {
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
    auto* sim = Simulation::GetActive();
    auto* ctxt= sim->GetExecutionContext();
    auto* rng = sim->GetRandom();
    n->SetPosition(ClampPoint(n->GetPosition()));

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *n, r2);
    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t help = P()->n_help_from_H * Sat(cnt.H, P()->n_help_half_sat);
    real_t birth_rate = P()->n_base_birth + help;
    real_t death_rate = P()->n_base_death
                      + P()->n_inactivation_by_C * static_cast<real_t>(cnt.C)
                      + P()->n_suppression_by_R  * static_cast<real_t>(cnt.R);

    if (rng->Uniform(0,1) < ProbFromRate(death_rate, dt_day)) {
      sim->GetExecutionContext()->RemoveAgent(n->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth_rate, dt_day)) {
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
    auto* sim = Simulation::GetActive();
    auto* ctxt= sim->GetExecutionContext();
    auto* rng = sim->GetRandom();
    h->SetPosition(ClampPoint(h->GetPosition()));

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *h, r2);
    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t self = P()->h_self_activation * Sat(cnt.H, P()->h_self_half_sat);
    real_t birth_rate = P()->h_base_birth + self;
    real_t death_rate = P()->h_base_death
                      + P()->h_suppression_by_R * static_cast<real_t>(cnt.R);

    if (rng->Uniform(0,1) < ProbFromRate(death_rate, dt_day)) {
      sim->GetExecutionContext()->RemoveAgent(h->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth_rate, dt_day)) {
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
    auto* sim = Simulation::GetActive();
    auto* ctxt= sim->GetExecutionContext();
    auto* rng = sim->GetRandom();
    r->SetPosition(ClampPoint(r->GetPosition()));

    Counts6 cnt;
    const real_t r2 = P()->local_radius_um * P()->local_radius_um;
    ctxt->ForEachNeighbor(cnt, *r, r2);
    const real_t dt_day = P()->dt_minutes / 1440.0;

    // Birth from basal + induction by E and H (linear in local counts)
    real_t birth_rate = P()->r_base_source
                      + P()->r_induction_by_E * static_cast<real_t>(cnt.E)
                      + P()->r_induction_by_H * static_cast<real_t>(cnt.H);

    // Death from decay + NK clearance
    real_t death_rate = P()->r_decay
                      + P()->r_clear_by_N * static_cast<real_t>(cnt.N);

    if (rng->Uniform(0,1) < ProbFromRate(death_rate, dt_day)) {
      sim->GetExecutionContext()->RemoveAgent(r->GetUid());
      return;
    }
    if (rng->Uniform(0,1) < ProbFromRate(birth_rate, dt_day)) {
      auto* d = r->Divide();
      if (auto* rd = dynamic_cast<TRegCell*>(d)) rd->SetCellColor(P()->color.treg);
    }
  }
};


// A tiny, invisible agent that just logs counts each step
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
    const real_t t_day = (P()->dt_minutes * steps) / 1440.0;

    // Count by type
    size_t C=0, Pn=0, E=0, N=0, H=0, R=0;
    rm->ForEachAgent([&](Agent* a) {
      if      (dynamic_cast<TumorCell*>(a))     ++C;
      else if (dynamic_cast<StellateCell*>(a))  ++Pn;
      else if (dynamic_cast<EffectorTCell*>(a)) ++E;
      else if (dynamic_cast<NKCell*>(a))        ++N;
      else if (dynamic_cast<HelperTCell*>(a))   ++H;
      else if (dynamic_cast<TRegCell*>(a))      ++R;
    });

    // Write CSV
    static std::ofstream csv("populations.csv");
    static bool wrote_header = false;
    if (!wrote_header) {
      csv << "step,days,C,P,E,N,H,R,total\n";
      wrote_header = true;
    }
    csv << steps << "," << t_day << ","
        << C << "," << Pn << "," << E << "," << N << "," << H << "," << R << ","
        << (C+Pn+E+N+H+R) << "\n";
    csv.flush();

    // Optional: print a quick daily summary
    if (steps % static_cast<size_t>(1440.0 / std::max<real_t>(1.0, P()->dt_minutes)) == 0) {
      std::cout << "[day " << t_day << "] C=" << C << " P=" << Pn
                << " E=" << E << " N=" << N << " H=" << H << " R=" << R << "\n";
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

  // Add one reporter agent that logs every step
  auto* rep = new ReporterCell();
  rep->AddBehavior(new ReportPopCounts());
  sim.GetExecutionContext()->AddAgent(rep);


  // Run (tune steps as you like; with dt=1 min, 1440 steps ≈ 1 day)
  // sim.GetScheduler()->Simulate(1440 * 100); // ~100 days
  sim.GetScheduler()->Simulate(1500); // ~0.5 days
  std::cout << "Pancreatic tumor (C,P,E,N,H,R) simple ABM completed.\n";
  return 0;
}

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PANCREATIC_TUMOR_MODEL_H_
