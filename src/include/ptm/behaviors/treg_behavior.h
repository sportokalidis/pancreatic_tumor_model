#ifndef PTM_BEHAVIORS_TREG_BEHAVIOR_H_
#define PTM_BEHAVIORS_TREG_BEHAVIOR_H_

#include "biodynamo.h"
#include "ptm/helpers.h"
#include "ptm/global_census.h"

namespace bdm {
namespace pancreatic_tumor {

class TRegBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TRegBehavior, Behavior, 1);
 public:
  inline void Run(Agent* a) override {
    auto* r   = bdm_static_cast<TRegCell*>(a);
    auto* ctxt= Simulation::GetActive()->GetExecutionContext();
    auto* rng = Simulation::GetActive()->GetRandom();
    r->SetPosition(ClampPoint(r->GetPosition()));

    auto& gc = GlobalCensus::Instance();
    gc.RefreshIfNeeded();

    const real_t dt_day = P()->dt_minutes / 1440.0;

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

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_BEHAVIORS_TREG_BEHAVIOR_H_
