#ifndef PTM_PARAMS_H_
#define PTM_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace pancreatic_tumor {

using real_t = double;

// -----------------------------------------------------------------------------
// Parameters (per day unless noted). Global interactions version.
// -----------------------------------------------------------------------------
struct Params {
  // Time & space
  real_t dt_minutes      = 1.0;   // scheduler step in minutes
  real_t min_bound       = 0.0;
  real_t max_bound       = 100.0;
  real_t cell_radius_um  = 6.0;

  // Initial agent counts (small ABM)
  size_t C0 = 474;  // Tumor
  size_t P0 = 10;   // PSC
  size_t E0 = 431;  // CD8+
  size_t N0 = 189;  // NK
  size_t H0 = 839;  // Helper T
  size_t R0 = 65;   // Tregs

  // Global carrying capacities (agents): logistic crowding per type uses these
  real_t K_C = 3000.0;
  real_t K_P = 3000.0;
  real_t K_E = 3000.0;
  real_t K_N = 3000.0;
  real_t K_H = 4000.0;
  real_t K_R = 3500.0;

  // Global gating: immune suppression becomes strong only when tumor is large
  real_t gate_C_K     = 300.0;    // half-saturation on global C

  // Generic half-sats for saturated interactions (in global agent counts)
  real_t K_small = 200.0;
  real_t K_med   = 600.0;
  real_t K_big   = 1200.0;

  // ================= Tumor (C) =================
  real_t c_base_div      = 0.05;   // baseline division
  real_t c_boost_from_P  = 0.25;   // P → C proliferative boost
  real_t c_boost_from_P_K= 800.0;  // global P half-sat for boost
  real_t c_kill_by_E     = 0.012;  // E→C killing
  real_t c_kill_by_E_K   = 600.0;  // global E half-sat
  real_t c_kill_by_N     = 0.007;  // N→C killing
  real_t c_kill_by_N_K   = 600.0;  // global N half-sat
  real_t c_R_blocks_E    = 0.10;   // 1/(1 + alpha * R) reduces E-kill

  // ================= PSC (P) =================
  real_t p_base_div      = 0.08;
  real_t p_boost_from_C  = 0.10;   // C → P boost
  real_t p_boost_from_C_K= 1000.0; // later boost onset (global)
  real_t p_base_death    = 0.05;

  // ================= Effector T (E) =================
  real_t e_base_birth    = 0.035;
  real_t e_help_from_H   = 0.22;   // H → E help
  real_t e_help_from_H_K = 600.0;  // global H half-sat
  real_t e_inact_by_C    = 0.15;   // C → E inactivation
  real_t e_inact_by_C_K  = 600.0;  // global C half-sat
  real_t e_suppr_by_R    = 0.04;   // R → E suppression (gated by C)
  real_t e_suppr_by_R_K  = 500.0;  // global R half-sat
  real_t e_base_death    = 0.10;

  // ================= NK (N) =================
  real_t n_base_birth    = 0.03;
  real_t n_help_from_H   = 0.15;
  real_t n_help_from_H_K = 600.0;
  real_t n_inact_by_C    = 0.05;
  real_t n_inact_by_C_K  = 600.0;
  real_t n_suppr_by_R    = 0.038;  // gated by C
  real_t n_suppr_by_R_K  = 500.0;
  real_t n_base_death    = 0.10;

  // ================= Helper T (H) =================
  real_t h_base_birth    = 0.05;
  real_t h_self_act      = 0.09;
  real_t h_self_act_K    = 1000.0; // weak self-activation, global
  real_t h_suppr_by_R    = 0.075;  // gated by C
  real_t h_suppr_by_R_K  = 700.0;
  real_t h_base_death    = 0.08;

  // ================= Tregs (R) =================
  real_t r_base_src      = 0.08;
  real_t r_induced_by_E  = 0.008;  // E → R induction
  real_t r_induced_by_E_K= 600.0;
  real_t r_induced_by_H  = 0.008;  // H → R induction
  real_t r_induced_by_H_K= 600.0;
  real_t r_cleared_by_N  = 0.003;  // N → R clearance
  real_t r_cleared_by_N_K= 500.0;
  real_t r_decay         = 0.06;

  // Colors
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

// Singleton accessor
inline Params* P() { static Params p; return &p; }

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_PARAMS_H_
