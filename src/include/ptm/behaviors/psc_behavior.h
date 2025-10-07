#ifndef PTM_BEHAVIORS_PSC_BEHAVIOR_H_
#define PTM_BEHAVIORS_PSC_BEHAVIOR_H_

#include "biodynamo.h"
#include "ptm/helpers.h"
#include "ptm/global_census.h"

namespace bdm {
namespace pancreatic_tumor {

class PSCBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(PSCBehavior, Behavior, 1);
 public:
  inline void Run(Agent* a) override {
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

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_BEHAVIORS_PSC_BEHAVIOR_H_
