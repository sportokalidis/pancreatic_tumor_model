#ifndef PTM_BEHAVIORS_TUMOR_BEHAVIOR_H_
#define PTM_BEHAVIORS_TUMOR_BEHAVIOR_H_

#include "biodynamo.h"
#include "ptm/helpers.h"
#include "ptm/global_census.h"

namespace bdm {
namespace pancreatic_tumor {

class TumorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TumorBehavior, Behavior, 1);
 public:
  inline void Run(Agent* a) override {
    auto* c   = bdm_static_cast<TumorCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    c->SetPosition(ClampPoint(c->GetPosition()));

    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();

    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t crowd = 1.0 - Clamp(static_cast<real_t>(gc.C) / std::max<real_t>(1.0, P()->K_C), 0.0, 1.0);
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

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_BEHAVIORS_TUMOR_BEHAVIOR_H_
