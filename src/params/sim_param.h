// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey
// BioDynaMo collaboration. Apache-2.0 license.
//
// -----------------------------------------------------------------------------
#ifndef PANCREATIC_SIM_PARAM_H_
#define PANCREATIC_SIM_PARAM_H_

#include "core/param/param_group.h"
#include "core/real_t.h"
#include <string>

namespace bdm {
namespace pancreatic_tumor {

// ============================================================================
// All simulation parameters for the pancreatic tumor model.
// Values are Table 1 from Akman Yıldız et al. (2021), scaled to ABM (S=1e5).
//
// Scale rules:
//   rate constant (day⁻¹)             → unchanged
//   bilinear coupling (cell day)⁻¹    → × S
//   source term (cells/day)            → ÷ S
//   half-saturation K (cells)          → ÷ S
//   carrying capacity: K_ABM = 1/(a·S)
// ============================================================================
struct SimParam : public ParamGroup {
  BDM_PARAM_GROUP_HEADER(SimParam, 1);

  // --------------------------------------------------------------------------
  // General
  // --------------------------------------------------------------------------
  int         seed        = 42;
  int         total_days  = 100;
  std::string output_dir  = "output";

  // --------------------------------------------------------------------------
  // Time & space
  // --------------------------------------------------------------------------
  real_t dt_minutes     = 1.0;
  real_t min_bound      = 0.0;
  real_t max_bound      = 150.0;
  real_t cell_radius_um = 6.0;

  // --------------------------------------------------------------------------
  // Interaction mode
  // --------------------------------------------------------------------------
  bool   use_local_counts = false;
  real_t local_radius_um  = 12.0;

  // --------------------------------------------------------------------------
  // Initial conditions (paper values ÷ S)
  // Paper: C0=4.886e7, P0=2.7362e5, E0=4.2684e7, N0=2.3531e7,
  //        H0=1.0343e8, R0=7.7570e6
  // --------------------------------------------------------------------------
  size_t C0 = 489;
  size_t P0 = 3;
  size_t E0 = 427;
  size_t N0 = 235;
  size_t H0 = 1034;
  size_t R0 = 78;

  // --------------------------------------------------------------------------
  // Carrying capacities (used in logistic crowding for C and P).
  // K_E, K_N, K_H are unused in ODE terms (paper has no cap for immune cells)
  // but kept as optional safety caps in SourceBehavior.
  // --------------------------------------------------------------------------
  real_t K_C = 8767.0;    // 1/(a_c·S)
  real_t K_P = 4384.0;    // 1/(a_p·S)
  real_t K_E = 3000.0;
  real_t K_N = 3000.0;
  real_t K_H = 4000.0;
  real_t K_R = 800000.0;

  // --------------------------------------------------------------------------
  // Tumor C  — Eq. 2.1
  // dC/dt = (k_c + mu_c·P)·C·(1-C/K_C) - b_c·N·C - d_c·E·C/(1+r1·R)
  // mu_c·P is LINEAR in P (not Hill). d_c·E and b_c·N are bilinear.
  // --------------------------------------------------------------------------
  real_t c_base_div     = 0.0886933;  // k_c  (day⁻¹)
  real_t c_boost_from_P = 2.02267e-6; // mu_c×S  (ABM cell day)⁻¹ — linear in P
  real_t c_kill_by_E    = 13.0;       // d_c×S   bilinear in E
  real_t c_kill_by_N    = 7.13e-6;    // b_c×S   bilinear in N
  real_t c_R_blocks_E   = 34500.0;    // r1×S    Treg blocks CTL 1/(1+r1·R)

  // --------------------------------------------------------------------------
  // PSC P  — Eq. 2.2
  // dP/dt = (k_p + f_p·C/(mu_p+C))·P·(1-P/K_P) - lambda_p·P
  // f_p·C/(mu_p+C) is Hill/Sat in C. lambda_p ≈ 0 → PSC grows to K_P.
  // --------------------------------------------------------------------------
  real_t p_base_div      = 0.00887;      // k_p (day⁻¹)
  real_t p_boost_from_C  = 0.154955;     // f_p (day⁻¹)
  real_t p_boost_from_C_K= 560.0;        // mu_p÷S half-sat for C
  real_t p_base_death    = 7.83296e-10;  // lambda_p (day⁻¹) ≈ 0

  // --------------------------------------------------------------------------
  // Effector T cell E  — Eq. 2.3
  // dE/dt = a_e - b_e·E - c_e·E·C + p_e·H·E/(g_e+H) - δ_e·R·E
  // c_e·C and δ_e·R are bilinear (no Hill K). p_e·H/(g_e+H) is Hill in H.
  // --------------------------------------------------------------------------
  real_t e_base_birth    = 0.13;        // a_e÷S  constant influx (cells/day)
  real_t e_help_from_H   = 0.125;       // p_e    Hill coefficient in H
  real_t e_help_from_H_K = 674490.0;    // g_e÷S  half-sat for H
  real_t e_inact_by_C    = 3.42e-5;     // c_e×S  bilinear in C
  real_t e_suppr_by_R    = 1.0e-5;      // δ_e×S  bilinear in R
  real_t e_base_death    = 0.027;       // b_e    natural death

  // --------------------------------------------------------------------------
  // NK cell N  — Eq. 2.4
  // dN/dt = a_n - b_n·N - c_n·N·C + p_n·H·N/(g_n+H) - δ_n·R·N
  // --------------------------------------------------------------------------
  real_t n_base_birth    = 0.13;
  real_t n_help_from_H   = 0.125;
  real_t n_help_from_H_K = 674490.0;
  real_t n_inact_by_C    = 1.0e-6;      // c_n×S  bilinear in C
  real_t n_suppr_by_R    = 1.0e-5;      // δ_n×S  bilinear in R
  real_t n_base_death    = 0.0412;

  // --------------------------------------------------------------------------
  // Helper T cell H  — Eq. 2.5
  // dH/dt = a_h - b_h·H + p_h·H²/(g_h+H) - δ_h·R·H
  // δ_h·R is bilinear. p_h·H/(g_h+H) is Hill per-cell proliferation.
  // --------------------------------------------------------------------------
  real_t h_base_birth  = 0.36;          // a_h÷S
  real_t h_self_act    = 0.125;         // p_h
  real_t h_self_act_K  = 674490.0;      // g_h÷S
  real_t h_suppr_by_R  = 1.0e-4;        // δ_h×S  bilinear in R
  real_t h_base_death  = 2.0833e-4;     // b_h

  // --------------------------------------------------------------------------
  // Treg R  — Eq. 2.6
  // dR/dt = a - δ_r·R + a_r·E + b_r·H + p_r·H·R/(g_r+H) - r·N·R
  // a_r·E and b_r·H are absolute source terms (handled by SourceBehavior).
  // r·N is bilinear death. p_r·H/(g_r+H) is per-cell Hill proliferation.
  // --------------------------------------------------------------------------
  real_t r_base_src      = 5.6;         // a÷S    constant source (cells/day)
  real_t r_induced_by_E  = 2.0e-4;      // a_r    rate coeff × E → abs. source
  real_t r_induced_by_H  = 4.0e-4;      // b_r    rate coeff × H → abs. source
  real_t r_prolif_by_H   = 0.125;       // p_r    per-cell Hill proliferation
  real_t r_prolif_by_H_K = 674490.0;    // g_r÷S  half-sat for H
  real_t r_cleared_by_N  = 1.0e-6;      // r×S    bilinear in N
  real_t r_decay         = 1.0e-5;      // δ_r    natural decay

  // --------------------------------------------------------------------------
  // Visualization colors
  // --------------------------------------------------------------------------
  int color_tumor     = 8;
  int color_tumor_div = 7;
  int color_psc       = 4;
  int color_eff       = 5;
  int color_nk        = 3;
  int color_helper    = 6;
  int color_treg      = 1;

  // --------------------------------------------------------------------------
  // Load / print
  // --------------------------------------------------------------------------
  void LoadParams(const std::string& filename);
  void PrintParams() const;
};

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PANCREATIC_SIM_PARAM_H_
