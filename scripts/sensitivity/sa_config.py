#!/usr/bin/env python3
"""
Single source of truth for the pancreatic-tumor ABM sensitivity analysis.

Replicates the global sensitivity analysis of Akman Yildiz et al. (2021),
J Biol Syst 29(4):799-832, Section 4 (Fig. 3): Latin-Hypercube sampling of the
Table 1 parameters over [0.5x, 2x] of baseline, with PRCC of outputs C (tumor)
and E (CTL) against each parameter at several time points.

Everything that defines *what* is varied lives here so the OAT, LHS, collect and
analysis scripts all agree. All SA runs use the fixed settings baked into the
baseline config: S=1e4, dt=24h (dt_minutes=1440), seed=42.
"""
from pathlib import Path

# Repo root = two levels up from scripts/sensitivity/
REPO = Path(__file__).resolve().parent.parent.parent
CONFIGS = REPO / "configs"

# Baseline config — fixed S=1e4, dt=24h, seed=42 (user-mandated SA settings).
BASELINE_CONFIG = CONFIGS / "params_S1e4_dt24h.json"

# Perturbation range as multiplicative factors of baseline (paper uses [0.5x, 2x]).
RANGE = (0.5, 2.0)

# Time points (paper-day coordinate, matches the `days` column in populations.csv,
# which equals sim_time + 7). Paper uses 7/35/70/150; our default sim runs
# total_days=100 -> days span 7..107, so we use 100 instead of 150. To reproduce
# the paper's day-150 point, run the sweep with total_days>=143 and add 150 here.
TIMEPOINTS = [7, 35, 70, 100]

# Outputs recorded at each time point. C and E are the paper's primary targets;
# the others are recorded for completeness.
OUTPUTS = ["C", "P", "E", "N", "H", "R"]

# ---------------------------------------------------------------------------
# Table 1 biological parameters to vary (ABM name -> paper symbol).
# Structural / numerical keys are deliberately EXCLUDED:
#   seed, dt_minutes, total_days, scale_S, bounds, cell/local radii, ICs
#   (C0..R0), carrying caps K_E/K_N/K_H/K_R (safety caps, not ODE terms),
#   colors, and all treatment/CSC flags.
#
# NOTE on K_C / K_P: the paper varies a_c, a_p (inverse carrying capacities),
# and K = 1/(a*S). Varying K over [0.5x, 2x] is equivalent to varying a over
# [2x, 0.5x]; a positive PRCC on K_C therefore corresponds to a *negative*
# sensitivity to a_c. This is documented, not corrected away.
# ---------------------------------------------------------------------------
PARAMS = {
    # Tumor C — Eq. 2.1
    "c_base_div":      "k_c",
    "c_boost_from_P":  "mu_c",
    "c_kill_by_E":     "d_c",
    "c_kill_by_N":     "b_c",
    "c_R_blocks_E":    "r_1",
    "K_C":             "1/a_c",
    # PSC P — Eq. 2.2
    "p_base_div":       "k_p",
    "p_boost_from_C":   "f_p",
    "p_boost_from_C_K": "mu_p",
    "p_base_death":     "lambda_p",
    "K_P":              "1/a_p",
    # Effector E — Eq. 2.3
    "e_base_birth":    "a_e",
    "e_help_from_H":   "p_e",
    "e_help_from_H_K": "g_e",
    "e_inact_by_C":    "c_e",
    "e_suppr_by_R":    "delta_e",
    "e_base_death":    "b_e",
    # NK N — Eq. 2.4
    "n_base_birth":    "a_n",
    "n_help_from_H":   "p_n",
    "n_help_from_H_K": "g_n",
    "n_inact_by_C":    "c_n",
    "n_suppr_by_R":    "delta_n",
    "n_base_death":    "b_n",
    # Helper H — Eq. 2.5
    "h_base_birth":  "a_h",
    "h_self_act":    "p_h",
    "h_self_act_K":  "g_h",
    "h_suppr_by_R":  "delta_h",
    "h_base_death":  "b_h",
    # Treg R — Eq. 2.6
    "r_base_src":      "a",
    "r_induced_by_E":  "a_r",
    "r_induced_by_H":  "b_r",
    "r_prolif_by_H":   "p_r",
    "r_prolif_by_H_K": "g_r",
    "r_cleared_by_N":  "r",
    "r_decay":         "delta_r",
}

PARAM_NAMES = list(PARAMS.keys())


def runs_dir() -> Path:
    """Directory where SA sweep folders live: runs/sensitivity/."""
    return REPO / "runs" / "sensitivity"
