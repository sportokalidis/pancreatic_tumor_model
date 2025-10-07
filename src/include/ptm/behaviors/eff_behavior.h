#ifndef PTM_BEHAVIORS_EFF_BEHAVIOR_H_
#define PTM_BEHAVIORS_EFF_BEHAVIOR_H_

#include "biodynamo.h"
#include "ptm/helpers.h"
#include "ptm/global_census.h"

namespace bdm {
namespace pancreatic_tumor {

class EffectorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(EffectorBehavior, Behavior, 1);
 public:
  inline void Run(Agent* a) override {
    auto* e   = bdm_static_cast<EffectorTCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    e->SetPosition(ClampPoint(e->GetPosition()));

    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();

    const real_t dt_day = P()->dt_minutes / 1440.0;

    real_t crowdE = 1.0 - Clamp(static_cast<real_t>(gc.E) / std::max<real_t>(1.0, P()->K_E), 0.0, 1.0);
    real_t helpH  = P()->e_help_from_H * Sat(static_cast<real_t>(gc.H), P()->e_help_from_H_K);
    real_t birth  = (P()->e_base_birth + helpH) * crowdE;

    real_t gateC = Sat(static_cast<real_t>(gc.C), P()->gate_C_K);
    real_t die = P()->e_base_death
               + P()->e_inact_by_C * Sat(static_cast<real_t>(gc.C), P()->e_inact_by_C_K) * gateC
               + P()->e_suppr_by_R * Sat(static_cast<real_t>(gc.R), P()->e_suppr_by_R_K) * gateC;

    // Smooth death as population shrinks (prevents hard crash to zero)
    die *= static_cast<real_t>(gc.E) / (gc.E + 1.0);

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

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_BEHAVIORS_EFF_BEHAVIOR_H_
