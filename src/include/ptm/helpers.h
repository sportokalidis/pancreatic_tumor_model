#ifndef PTM_HELPERS_H_
#define PTM_HELPERS_H_

#include "biodynamo.h"
#include "ptm/params.h"
#include <cmath>
#include <algorithm>

namespace bdm {
namespace pancreatic_tumor {

// -----------------------------------------------------------------------------
// Helper functions (extracted from main header to avoid circular dependencies)
// -----------------------------------------------------------------------------
inline real_t Clamp(real_t v, real_t lo, real_t hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

inline Real3 ClampPoint(const Real3& pos) {
  return Real3{Clamp(pos[0], P()->min_bound, P()->max_bound),
               Clamp(pos[1], P()->min_bound, P()->max_bound),
               Clamp(pos[2], P()->min_bound, P()->max_bound)};
}

inline real_t Sat(real_t x, real_t K) { 
  return (x <= 0) ? 0.0 : x / (K + x); 
}

inline real_t ProbFromRate(real_t rate_per_day, real_t dt_day) {
  if (rate_per_day <= 0) return 0.0;
  real_t v = std::exp(-rate_per_day * dt_day);
  return 1.0 - v;
}

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_HELPERS_H_