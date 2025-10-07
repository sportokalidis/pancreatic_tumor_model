#ifndef PTM_BEHAVIORS_NK_BEHAVIOR_H_
#define PTM_BEHAVIORS_NK_BEHAVIOR_H_

#include "biodynamo.h"
#include "ptm/helpers.h"
#include "ptm/global_census.h"

namespace bdm {
namespace pancreatic_tumor {

class NKBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(NKBehavior, Behavior, 1);
 public:
  inline void Run(Agent* a) override {
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

    // Smooth death at low N
    die *= static_cast<real_t>(gc.N) / (gc.N + 1.0);

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

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_BEHAVIORS_NK_BEHAVIOR_H_
