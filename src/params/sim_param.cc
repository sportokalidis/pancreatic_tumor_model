// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey
// BioDynaMo collaboration. Apache-2.0 license.
//
// -----------------------------------------------------------------------------
#include "params/sim_param.h"
#include "core/param/param_group.h"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <stdexcept>

namespace bdm {
namespace pancreatic_tumor {

const ParamGroupUid SimParam::kUid = ParamGroupUidGenerator::Get()->NewUid();

void SimParam::LoadParams(const std::string& filename) {
  nlohmann::json jf;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cout << "[SimParam] params.json not found — using defaults.\n";
    return;
  }
  try {
    file >> jf;
  } catch (const std::exception& e) {
    std::cerr << "[SimParam] JSON parse error in " << filename
              << ": " << e.what() << "\n";
    return;
  }

  auto load_real = [&](const char* key, real_t& v) {
    if (jf.contains(key)) v = jf[key].get<double>();
  };
  auto load_int = [&](const char* key, int& v) {
    if (jf.contains(key)) v = jf[key].get<int>();
  };
  auto load_size = [&](const char* key, size_t& v) {
    if (jf.contains(key)) v = jf[key].get<size_t>();
  };
  auto load_bool = [&](const char* key, bool& v) {
    if (jf.contains(key)) v = jf[key].get<bool>();
  };
  auto load_str = [&](const char* key, std::string& v) {
    if (jf.contains(key)) v = jf[key].get<std::string>();
  };

  // General
  load_int("seed",        seed);
  load_int("total_days",  total_days);
  load_str("output_dir",  output_dir);

  // Time & space
  load_real("dt_minutes",     dt_minutes);
  load_real("min_bound",      min_bound);
  load_real("max_bound",      max_bound);
  load_real("cell_radius_um", cell_radius_um);

  // Interaction mode
  load_bool("use_local_counts", use_local_counts);
  load_real("local_radius_um",  local_radius_um);
  load_real("immune_step_um",   immune_step_um);

  // Spatial initialization
  load_bool("use_sphere_init",     use_sphere_init);
  load_real("tumor_sphere_radius", tumor_sphere_radius);

  // Initial counts
  load_size("C0", C0);
  load_size("P0", P0);
  load_size("E0", E0);
  load_size("N0", N0);
  load_size("H0", H0);
  load_size("R0", R0);

  // Carrying capacities
  load_real("K_C", K_C);
  load_real("K_P", K_P);
  load_real("K_E", K_E);
  load_real("K_N", K_N);
  load_real("K_H", K_H);
  load_real("K_R", K_R);

  // Tumor C (Eq. 2.1)
  load_real("c_base_div",     c_base_div);
  load_real("c_boost_from_P", c_boost_from_P);
  load_real("c_kill_by_E",    c_kill_by_E);
  load_real("c_kill_by_N",    c_kill_by_N);
  load_real("c_R_blocks_E",   c_R_blocks_E);

  // PSC P (Eq. 2.2)
  load_real("p_base_div",       p_base_div);
  load_real("p_boost_from_C",   p_boost_from_C);
  load_real("p_boost_from_C_K", p_boost_from_C_K);
  load_real("p_base_death",     p_base_death);

  // Effector E (Eq. 2.3)
  load_real("e_base_birth",    e_base_birth);
  load_real("e_help_from_H",   e_help_from_H);
  load_real("e_help_from_H_K", e_help_from_H_K);
  load_real("e_inact_by_C",    e_inact_by_C);
  load_real("e_suppr_by_R",    e_suppr_by_R);
  load_real("e_base_death",    e_base_death);

  // NK N (Eq. 2.4)
  load_real("n_base_birth",    n_base_birth);
  load_real("n_help_from_H",   n_help_from_H);
  load_real("n_help_from_H_K", n_help_from_H_K);
  load_real("n_inact_by_C",    n_inact_by_C);
  load_real("n_suppr_by_R",    n_suppr_by_R);
  load_real("n_base_death",    n_base_death);

  // Helper H (Eq. 2.5)
  load_real("h_base_birth",  h_base_birth);
  load_real("h_self_act",    h_self_act);
  load_real("h_self_act_K",  h_self_act_K);
  load_real("h_suppr_by_R",  h_suppr_by_R);
  load_real("h_base_death",  h_base_death);

  // Treg R (Eq. 2.6)
  load_real("r_base_src",      r_base_src);
  load_real("r_induced_by_E",  r_induced_by_E);
  load_real("r_induced_by_H",  r_induced_by_H);
  load_real("r_prolif_by_H",   r_prolif_by_H);
  load_real("r_prolif_by_H_K", r_prolif_by_H_K);
  load_real("r_cleared_by_N",  r_cleared_by_N);
  load_real("r_decay",         r_decay);

  // Drug treatment (Section 5) — optional, default off
  load_bool("treat_gem",   treat_gem);
  load_bool("treat_abr",   treat_abr);
  load_bool("treat_acd47", treat_acd47);
  load_real("treat_start_day", treat_start_day);
  load_real("gem_gamma",     gem_gamma);
  load_real("gem_dose",      gem_dose);
  load_real("gem_freq_days", gem_freq_days);
  load_real("gem_end_day",   gem_end_day);
  load_real("gem_c_c",       gem_c_c);
  load_real("gem_c_p",       gem_c_p);
  load_real("gem_c_immune",  gem_c_immune);
  load_real("abr_gamma",     abr_gamma);
  load_real("abr_dose",      abr_dose);
  load_real("abr_freq_days", abr_freq_days);
  load_real("abr_end_day",   abr_end_day);
  load_real("abr_c_c",       abr_c_c);
  load_real("abr_c_p",       abr_c_p);
  load_real("abr_c_immune",  abr_c_immune);
  load_real("acd47_end_day",  acd47_end_day);
  load_real("acd47_e_boost",  acd47_e_boost);

  // Cancer Stem Cell (Section 6) — optional, default off
  load_bool("csc_enable",     csc_enable);
  load_size("S0",             S0);
  load_real("csc_a1",         csc_a1);
  load_real("csc_a2",         csc_a2);
  load_real("csc_a3",         csc_a3);
  load_real("csc_lambda_max", csc_lambda_max);
  load_real("csc_lambda_min", csc_lambda_min);
  load_real("csc_sigma",      csc_sigma);
  load_real("csc_delta_s",    csc_delta_s);

  // Colors
  load_int("color_tumor",     color_tumor);
  load_int("color_tumor_div", color_tumor_div);
  load_int("color_psc",       color_psc);
  load_int("color_eff",       color_eff);
  load_int("color_nk",        color_nk);
  load_int("color_helper",    color_helper);
  load_int("color_treg",      color_treg);
  load_int("color_csc",       color_csc);
}

void SimParam::PrintParams() const {
  std::cout << "\n=== SimParam (Table 1 — Akman Yildiz 2021, S=1e5) ===\n"
            << "  seed=" << seed << "  total_days=" << total_days
            << "  output_dir=" << output_dir << "\n"
            << "  dt_minutes=" << dt_minutes
            << "  domain=[" << min_bound << "," << max_bound << "]"
            << "  cell_radius=" << cell_radius_um << " um\n"
            << "  use_local=" << use_local_counts
            << "  local_radius=" << local_radius_um << " um"
            << "  immune_step=" << immune_step_um << " um"
            << "  sphere_init=" << use_sphere_init
            << "  sphere_r=" << tumor_sphere_radius << " um\n"
            << "  C0=" << C0 << " P0=" << P0 << " E0=" << E0
            << " N0=" << N0 << " H0=" << H0 << " R0=" << R0 << "\n"
            << "  K_C=" << K_C << " K_P=" << K_P << "\n"
            << "  c_base_div=" << c_base_div
            << "  c_R_blocks_E=" << c_R_blocks_E
            << "  c_kill_by_E=" << c_kill_by_E << "\n"
            << "  p_base_death=" << p_base_death
            << "  e_base_birth=" << e_base_birth
            << "  e_base_death=" << e_base_death << "\n"
            << "  r_base_src=" << r_base_src
            << "  r_decay=" << r_decay << "\n"
            << "=== end SimParam ===\n";
}

}  // namespace pancreatic_tumor
}  // namespace bdm
