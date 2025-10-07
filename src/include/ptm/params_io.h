#ifndef PTM_PARAMS_IO_H_
#define PTM_PARAMS_IO_H_

#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include "pancreatic_tumor_model.h"

namespace bdm {
namespace pancreatic_tumor {

// Very small, dependency-free parser for "key = value" files.
// Lines starting with '#' are ignored.
// Only doubles/size_t are supported here (extend as needed).
inline bool LoadParamsFromFile(const std::string& path, Params* out) {
  if (!out) return false;
  std::ifstream in(path);
  if (!in.good()) return false;

  // Map keys to setter lambdas
  std::unordered_map<std::string, std::function<void(const std::string&)>> set;

  auto set_double = [&](const std::string& k, double* ptr) {
    set[k] = [ptr](const std::string& v) { *ptr = std::stod(v); };
  };
  auto set_size   = [&](const std::string& k, size_t* ptr) {
    set[k] = [ptr](const std::string& v) { *ptr = static_cast<size_t>(std::stoll(v)); };
  };

  // Register fields (add more if you need)
  set_double("dt_minutes", &out->dt_minutes);
  set_double("min_bound", &out->min_bound);
  set_double("max_bound", &out->max_bound);
  set_double("cell_radius_um", &out->cell_radius_um);

  set_size("C0", &out->C0); set_size("P0", &out->P0); set_size("E0", &out->E0);
  set_size("N0", &out->N0); set_size("H0", &out->H0); set_size("R0", &out->R0);

  set_double("K_C", &out->K_C); set_double("K_P", &out->K_P);
  set_double("K_E", &out->K_E); set_double("K_N", &out->K_N);
  set_double("K_H", &out->K_H); set_double("K_R", &out->K_R);
  set_double("gate_C_K", &out->gate_C_K);

  set_double("K_small", &out->K_small);
  set_double("K_med", &out->K_med);
  set_double("K_big", &out->K_big);

  set_double("c_base_div", &out->c_base_div);
  set_double("c_boost_from_P", &out->c_boost_from_P);
  set_double("c_boost_from_P_K", &out->c_boost_from_P_K);
  set_double("c_kill_by_E", &out->c_kill_by_E);
  set_double("c_kill_by_E_K", &out->c_kill_by_E_K);
  set_double("c_kill_by_N", &out->c_kill_by_N);
  set_double("c_kill_by_N_K", &out->c_kill_by_N_K);
  set_double("c_R_blocks_E", &out->c_R_blocks_E);

  set_double("p_base_div", &out->p_base_div);
  set_double("p_boost_from_C", &out->p_boost_from_C);
  set_double("p_boost_from_C_K", &out->p_boost_from_C_K);
  set_double("p_base_death", &out->p_base_death);

  set_double("e_base_birth", &out->e_base_birth);
  set_double("e_help_from_H", &out->e_help_from_H);
  set_double("e_help_from_H_K", &out->e_help_from_H_K);
  set_double("e_inact_by_C", &out->e_inact_by_C);
  set_double("e_inact_by_C_K", &out->e_inact_by_C_K);
  set_double("e_suppr_by_R", &out->e_suppr_by_R);
  set_double("e_suppr_by_R_K", &out->e_suppr_by_R_K);
  set_double("e_base_death", &out->e_base_death);

  set_double("n_base_birth", &out->n_base_birth);
  set_double("n_help_from_H", &out->n_help_from_H);
  set_double("n_help_from_H_K", &out->n_help_from_H_K);
  set_double("n_inact_by_C", &out->n_inact_by_C);
  set_double("n_inact_by_C_K", &out->n_inact_by_C_K);
  set_double("n_suppr_by_R", &out->n_suppr_by_R);
  set_double("n_suppr_by_R_K", &out->n_suppr_by_R_K);
  set_double("n_base_death", &out->n_base_death);

  set_double("h_base_birth", &out->h_base_birth);
  set_double("h_self_act", &out->h_self_act);
  set_double("h_self_act_K", &out->h_self_act_K);
  set_double("h_suppr_by_R", &out->h_suppr_by_R);
  set_double("h_suppr_by_R_K", &out->h_suppr_by_R_K);
  set_double("h_base_death", &out->h_base_death);

  set_double("r_base_src", &out->r_base_src);
  set_double("r_induced_by_E", &out->r_induced_by_E);
  set_double("r_induced_by_E_K", &out->r_induced_by_E_K);
  set_double("r_induced_by_H", &out->r_induced_by_H);
  set_double("r_induced_by_H_K", &out->r_induced_by_H_K);
  set_double("r_cleared_by_N", &out->r_cleared_by_N);
  set_double("r_cleared_by_N_K", &out->r_cleared_by_N_K);
  set_double("r_decay", &out->r_decay);

  std::string line;
  while (std::getline(in, line)) {
    // trim
    auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    auto is_space = [](unsigned char c){ return std::isspace(c); };
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), [&](char c){ return !is_space(c); }));
    line.erase(std::find_if(line.rbegin(), line.rend(), [&](char c){ return !is_space(c); }).base(), line.end());
    if (line.empty()) continue;

    auto eq = line.find('=');
    if (eq == std::string::npos) continue;

    std::string key = line.substr(0, eq);
    std::string val = line.substr(eq + 1);
    // trim key and val
    key.erase(key.begin(), std::find_if(key.begin(), key.end(), [&](char c){ return !is_space(c); }));
    key.erase(std::find_if(key.rbegin(), key.rend(), [&](char c){ return !is_space(c); }).base(), key.end());
    val.erase(val.begin(), std::find_if(val.begin(), val.end(), [&](char c){ return !is_space(c); }));
    val.erase(std::find_if(val.rbegin(), val.rend(), [&](char c){ return !is_space(c); }).base(), val.end());

    auto it = set.find(key);
    if (it != set.end()) {
      try { it->second(val); } catch (...) { /* ignore parse errors */ }
    }
  }
  return true;
}

}  // namespace pancreatic_tumor
}  // namespace bdm

#endif  // PTM_PARAMS_IO_H_
