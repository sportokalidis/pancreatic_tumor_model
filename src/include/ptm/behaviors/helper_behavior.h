#ifndef PTM_BEHAVIORS_HELPER_BEHAVIOR_H_
#define PTM_BEHAVIORS_HELPER_BEHAVIOR_H_

#include "biodynamo.h"
#include "ptm/helpers.h"
#include "ptm/global_census.h"

namespace bdm {
namespace pancreatic_tumor {

class HelperBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(HelperBehavior, Behavior, 1);
 public:
  inline void Run(Agent* a) override {
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

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_BEHAVIORS_HELPER_BEHAVIOR_H_
