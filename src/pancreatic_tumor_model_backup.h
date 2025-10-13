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
#include <mutex>

namespace {
  std::mutex cout_mtx;
}

namespace bdm {
namespace pancreatic_tumor {

// All components are now included from modular headers:
// - Params struct and P() function from ptm/params.h (with file loading)
// - Helper functions from ptm/helpers.h
// - Cell types from ptm/cell_types.h  
// - GlobalCensus from ptm/global_census.h
// - Behaviors from ptm/behaviors/*.h

// -------------------------------------------------------
  BDM_BEHAVIOR_HEADER(TumorBehavior, Behavior, 1);
 public:
  void Run(Agent* a) override {
    auto* c   = bdm_static_cast<TumorCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    c->SetPosition(ClampPoint(c->GetPosition()));

    // refresh global counts once per step
    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();

    const real_t dt_day = P()->dt_minutes / 1440.0;

    // Logistic crowding on global C
    real_t crowd = 1.0 - Clamp(static_cast<real_t>(gc.C) / std::max<real_t>(1.0, P()->K_C), 0.0, 1.0);

    // Birth: baseline + global P boost
    real_t boostP = P()->c_boost_from_P * Sat(static_cast<real_t>(gc.Pn), P()->c_boost_from_P_K);
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

    // Death: E and N global killing; E effectiveness reduced by R; gated by C
    real_t gateC = Sat(static_cast<real_t>(gc.C), P()->gate_C_K);
    real_t inhibit_E = 1.0 / (1.0 + P()->c_R_blocks_E * static_cast<real_t>(gc.R));
    real_t killE = P()->c_kill_by_E * Sat(static_cast<real_t>(gc.E), P()->c_kill_by_E_K) * inhibit_E * gateC;
    real_t killN = P()->c_kill_by_N * Sat(static_cast<real_t>(gc.N), P()->c_kill_by_N_K) * gateC;
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

    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();

    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t crowd  = 1.0 - Clamp(static_cast<real_t>(gc.Pn) / std::max<real_t>(1.0, P()->K_P), 0.0, 1.0);
    real_t boostC = P()->p_boost_from_C * Sat(static_cast<real_t>(gc.C), P()->p_boost_from_C_K);
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

    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();

    const real_t dt_day = P()->dt_minutes / 1440.0;

    // Birth: baseline + help from global H, moderated by global E crowding
    real_t crowdE = 1.0 - Clamp(static_cast<real_t>(gc.E) / std::max<real_t>(1.0, P()->K_E), 0.0, 1.0);
    real_t helpH  = P()->e_help_from_H * Sat(static_cast<real_t>(gc.H), P()->e_help_from_H_K);
    real_t birth  = (P()->e_base_birth + helpH) * crowdE;

    // Death/suppression: inactivation by C & suppression by R, gated by C
    real_t gateC = Sat(static_cast<real_t>(gc.C), P()->gate_C_K);
    real_t die = P()->e_base_death
               + P()->e_inact_by_C * Sat(static_cast<real_t>(gc.C), P()->e_inact_by_C_K) * gateC
               + P()->e_suppr_by_R * Sat(static_cast<real_t>(gc.R), P()->e_suppr_by_R_K) * gateC;

    die *= static_cast<real_t>(gc.E) / (gc.E + 1.0);
    // die *= std::pow(static_cast<real_t>(gc.E) / (gc.E + 50), 0.5);

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

    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();

    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t crowdN = 1.0 - Clamp(static_cast<real_t>(gc.N) / std::max<real_t>(1.0, P()->K_N), 0.0, 1.0);
    real_t helpH  = P()->n_help_from_H * Sat(static_cast<real_t>(gc.H), P()->n_help_from_H_K);
    real_t birth  = (P()->n_base_birth + helpH) * crowdN;

    real_t gateC = Sat(static_cast<real_t>(gc.C), P()->gate_C_K);
    real_t die = P()->n_base_death
               + P()->n_inact_by_C * Sat(static_cast<real_t>(gc.C), P()->n_inact_by_C_K) * gateC
               + P()->n_suppr_by_R * Sat(static_cast<real_t>(gc.R), P()->n_suppr_by_R_K) * gateC;

    die *= static_cast<real_t>(gc.N) / (gc.N + 1.0);
    // die *= std::pow(static_cast<real_t>(gc.N) / (gc.N + 20), 0.5);

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

    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();

    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t crowdH = 1.0 - Clamp(static_cast<real_t>(gc.H) / std::max<real_t>(1.0, P()->K_H), 0.0, 1.0);
    real_t self   = P()->h_self_act * Sat(static_cast<real_t>(gc.H), P()->h_self_act_K);
    real_t birth  = (P()->h_base_birth + self) * crowdH;

    real_t gateC = Sat(static_cast<real_t>(gc.C), P()->gate_C_K);
    real_t die = P()->h_base_death
               + P()->h_suppr_by_R * Sat(static_cast<real_t>(gc.R), P()->h_suppr_by_R_K) * gateC;


    // {
    //   std::ostringstream ss;
    //   // std::cout << "div_rate" << div_rate << ", p_kill " << p_kill << "\n";
    //   ss << "[H] step=" << Simulation::GetActive()->GetScheduler()->GetSimulatedSteps()
    //      << " birth=" << birth << " die=" << die << "\n";
    //   std::lock_guard<std::mutex> lock(cout_mtx);
    //   std::cout << ss.str();       // serialized
    //   std::cout.flush();
    // }

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

    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();

    const real_t dt_day = P()->dt_minutes / 1440.0;

    // R rises from baseline source + induction by global E and H
    real_t crowdR = 1.0 - Clamp(static_cast<real_t>(gc.R) / std::max<real_t>(1.0, P()->K_R), 0.0, 1.0);
    real_t birth  = (P()->r_base_src
                  +  P()->r_induced_by_E * Sat(static_cast<real_t>(gc.E), P()->r_induced_by_E_K)
                  +  P()->r_induced_by_H * Sat(static_cast<real_t>(gc.H), P()->r_induced_by_H_K))
                  *  crowdR;

    real_t die    = P()->r_decay
                  + P()->r_cleared_by_N * Sat(static_cast<real_t>(gc.N), P()->r_cleared_by_N_K);

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

// -------------------------------------------------------
// Reporter (daily CSV based on global counts)
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
  std::cout << "Pancreatic tumor (C,P,E,N,H,R) ABM (global interactions) completed.\n";
  return 0;
}

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PANCREATIC_TUMOR_MODEL_H_
