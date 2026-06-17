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
  bool   use_local_counts = true;
  real_t local_radius_um  = 36.0;   // 3 cell diameters — paracrine signaling range
  real_t immune_step_um   = 10.0;   // Gaussian σ per coordinate per step for immune random walk

  // --------------------------------------------------------------------------
  // Spatial initialization
  // --------------------------------------------------------------------------
  bool   use_sphere_init      = false;
  real_t tumor_sphere_radius  = 120.0;  // µm — compact tumor mass at domain center

  // --------------------------------------------------------------------------
  // Initial conditions (paper values ÷ S)
  // Paper: C0=4.886e7, P0=2.7362e5, E0=4.2684e7, N0=2.3531e7,
  //        H0=1.0343e8, R0=7.7570e6
  // --------------------------------------------------------------------------
  size_t C0 = 4890;
  size_t P0 = 30;
  size_t E0 = 4270;
  size_t N0 = 2350;
  size_t H0 = 10340;
  size_t R0 = 780;

  // --------------------------------------------------------------------------
  // Carrying capacities (used in logistic crowding for C and P).
  // K_E, K_N, K_H are unused in ODE terms (paper has no cap for immune cells)
  // but kept as optional safety caps in SourceBehavior.
  // --------------------------------------------------------------------------
  real_t K_C = 87670.0;    // 1/(a_c·S)
  real_t K_P = 43840.0;    // 1/(a_p·S)
  real_t K_E = 30000.0;
  real_t K_N = 30000.0;
  real_t K_H = 40000.0;
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
  // Drug treatment (Paper Section 5, Eqs. 5.1–5.7) — OPTIONAL, default OFF.
  //
  // Drug concentration M evolves as: dM/dt = -γ·M + v_b(t)
  // Drug kill per cell type x: rate += c_x·(1−e^{−M})
  // Anti-CD47: +α·E·v_a(t) extra CTL proliferation when active
  //
  // Schedules use paper-day coordinates (simulation day 0 = paper day 7).
  // Injection at treat_start_day, then every *_freq_days, up to *_end_day.
  // --------------------------------------------------------------------------
  bool treat_gem   = false;  // enable Gemcitabine
  bool treat_abr   = false;  // enable Abraxane
  bool treat_acd47 = false;  // enable Anti-CD47

  real_t treat_start_day = 14.0;  // paper day all treatments begin

  // Gemcitabine — Table 2 values; c_* in day⁻¹
  real_t gem_gamma    = 5.54;          // γ_Gem  drug decay rate
  real_t gem_dose     = 1.0;           // M increment per injection
  real_t gem_freq_days = 7.0;          // injection interval (days)
  real_t gem_end_day  = 56.0;          // last paper day of treatment
  real_t gem_c_c      = 10.0;          // kill rate for C
  real_t gem_c_p      = 7.5;           // kill rate for P
  real_t gem_c_immune = 1.8;           // kill rate for E, N, H, R

  // Abraxane — Table 2 values
  real_t abr_gamma    = 0.6161308;
  real_t abr_dose     = 1.0;
  real_t abr_freq_days = 7.0;
  real_t abr_end_day  = 28.0;
  real_t abr_c_c      = 10.0;
  real_t abr_c_p      = 10.0;
  real_t abr_c_immune = 6.4;

  // Anti-CD47 — α·E·v_a term: effective per-cell CTL boost rate (day⁻¹)
  real_t acd47_end_day  = 35.0;
  real_t acd47_e_boost  = 0.5;         // α·v_a_dose (tune vs Fig. 5)

  // --------------------------------------------------------------------------
  // Cancer Stem Cell (Paper Section 6, Eqs. 6.1 & 6.7) — OPTIONAL, default OFF.
  // dS/dt = λ(C)·(a1-a3)·S - δ_S·S                                (Eq. 6.7)
  //   → λ_max when C<<σ (tumor suppressed — CSCs proliferate fast)
  //   → λ_min when C>>σ (tumor large  — CSCs grow minimally)
  // dC/dt adds +(a2+2·a3)·S to Eq. 2.1.                          (Eq. 6.1)
  //   CSC asymmetric / symmetric-to-PCC divisions supply new tumor cells.
  //   a1=0.35, a2=0.65, a3=0, λ_min=0.001, λ_max=0.07, σ=1e7 cells (paper)
  //   δ_S = 0 — paper neglects natural CSC death for aggressive progression
  // --------------------------------------------------------------------------
  bool   csc_enable     = false;
  size_t S0             = 0;      // initial CSC count (usually 0 at t=0)
  real_t csc_a1         = 0.35;   // prob. symmetric division CSC→CSC+CSC   [Table 3]
  real_t csc_a2         = 0.65;   // prob. asymmetric division CSC→CSC+PCC
  real_t csc_a3         = 0.0;    // prob. symmetric division CSC→PCC+PCC (rare)
  real_t csc_lambda_max = 0.07;   // λ_max: max CSC proliferation (day⁻¹)   [Table 3]
  real_t csc_lambda_min = 0.001;  // λ_min: min CSC proliferation (day⁻¹)   [Table 3]
  real_t csc_sigma      = 100.0;  // σ÷S: tumor threshold in ABM cells (1e7÷1e5) [Table 3]
  real_t csc_delta_s    = 0.0;    // δ_S: CSC death rate (day⁻¹); 0 per paper

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
  int color_csc       = 2;

  // --------------------------------------------------------------------------
  // Load / print
  // --------------------------------------------------------------------------
  void LoadParams(const std::string& filename);
  void PrintParams() const;
};

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PANCREATIC_SIM_PARAM_H_
