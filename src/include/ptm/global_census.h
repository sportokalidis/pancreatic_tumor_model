#ifndef PTM_GLOBAL_CENSUS_H_
#define PTM_GLOBAL_CENSUS_H_

#include "biodynamo.h"
#include "ptm/cell_types.h"
#include <limits>

namespace bdm {
namespace pancreatic_tumor {

// Global totals computed once per scheduler step.
// Header-only with inline definition.
struct GlobalCensus {
  size_t step_cached = std::numeric_limits<size_t>::max();
  size_t C=0, Pn=0, E=0, N=0, H=0, R=0;

  static GlobalCensus& Instance() {
    static GlobalCensus gc;
    return gc;
  }

  inline void RefreshIfNeeded() {
    auto* sim = Simulation::GetActive();
    auto steps = sim->GetScheduler()->GetSimulatedSteps();
    if (steps == step_cached) return;

    size_t C_=0, P_=0, E_=0, N_=0, H_=0, R_=0;
    sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
      if      (dynamic_cast<TumorCell*>(a))     ++C_;
      else if (dynamic_cast<StellateCell*>(a))  ++P_;
      else if (dynamic_cast<EffectorTCell*>(a)) ++E_;
      else if (dynamic_cast<NKCell*>(a))        ++N_;
      else if (dynamic_cast<HelperTCell*>(a))   ++H_;
      else if (dynamic_cast<TRegCell*>(a))      ++R_;
    });
    C  = C_;  Pn = P_;  E = E_;  N = N_;  H = H_;  R = R_;
    step_cached = steps;
  }
};

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_GLOBAL_CENSUS_H_
