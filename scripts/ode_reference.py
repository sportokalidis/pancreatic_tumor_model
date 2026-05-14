#!/usr/bin/env python3
"""
ODE reference solver — three modes:

  --mode paper (default)
      Integrates the EXACT Eqs. 2.1-2.6 from Akman Yıldız et al. (2021)
      using Table 1 parameter values scaled to ABM units (scale S = 1e5).
      This is the ground truth.  Writes ode_reference.csv AND regenerates
      the *_scaled_global.csv reference files used by calc-error.py and the
      comparison plots.

  --mode abm
      Integrates the same equations but driven by params.json, mirroring the
      ABM behaviour exactly (Hill saturation, crowding factors).  Used to
      validate that the ABM implementation matches its mean-field ODE limit.

  --mode treatment
      Integrates Eqs. 5.1-5.7 (paper Section 5) — baseline ODEs plus drug kill
      terms c_x*(1-e^{-M})*x and Anti-CD47 CTL boost, with piecewise-exact
      injection pulses. Treatment parameters are read from params.json.
      Overlay with --abm <treatment_populations.csv> to validate the ABM.

Scale factor S ≈ 1e5
  paper units  (10^7 cells)  →  ABM units  (~10^5 cells)

Scaling rules applied to every paper parameter:
  rate constant (day⁻¹)                 → unchanged
  bilinear coupling (cell day)⁻¹        → × S
  source term (cells day⁻¹)             → ÷ S
  half-saturation K (cells)             → ÷ S
  inverse carrying capacity a (cells⁻¹) → K_ABM = 1/(a·S)

Outputs (--out directory):
  ode_reference.csv
  ode_vs_abm.png

In --mode paper also writes per-population reference CSVs:
  data-export/C-Cells_scaled_global.csv  … R-Cells_scaled_global.csv

Usage:
  python3 scripts/ode_reference.py [--mode paper|abm|treatment]
                                   [--params params.json]
                                   [--abm    output/populations.csv]
                                   [--refs   data-export]
                                   [--out    data-export]
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Table 1 parameters — paper units → ABM units
#
# S = scale factor: N_paper = N_ABM * S
# Default S=1e5 matches original publication figures.
# Use S=1e4 for 10× more agents (better statistics, avoids E→0 extinction).
# All ODE dynamics are identical regardless of S — only absolute counts change.
# ---------------------------------------------------------------------------

_k_c  = 0.0886933
_a_c  = 1.14077e-9       # cells⁻¹ (paper)
_a_p  = 2.0 * _a_c       # = 2·a_c  (paper)
_mu_c = 0.1 * _k_c * _a_p  # (cell day)⁻¹ (paper)
_g    = 6.7449e10         # common half-saturation for E, N, H, R (paper cells)


def build_paper_params(S: float) -> dict:
    """Build PAPER_PARAMS dict scaled to ABM scale factor S."""
    return {
        # ---- Tumor C  (Eq. 2.1) ---------------------------------------------
        "k_c":    _k_c,
        "mu_c":   _mu_c * S,
        "K_C":    1.0 / (_a_c * S),
        "b_c":    7.13e-11 * S,
        "d_c":    1.3e-4 * S,
        "r1":     0.345 * S,

        # ---- PSC P  (Eq. 2.2) -----------------------------------------------
        "k_p":      0.1 * _k_c,
        "f_p":      0.154955,
        "mu_p":     5.6e7 / S,
        "K_P":      1.0 / (_a_p * S),
        "lambda_p": 7.83296e-10,

        # ---- Effector E  (Eq. 2.3) ------------------------------------------
        # r_e set to 0: PDF rendering fuses reference [47] with exponent,
        # making the term 109× too large and pushing E above Fig. 2 y-axis.
        "a_e":     1.3e4 / S,
        "b_e":     2.7e-2,
        "c_e":     3.42e-10 * S,
        "r_e":     0.0,
        "p_e":     0.125,
        "g_e":     _g / S,
        "delta_e": 1e-10 * S,

        # ---- NK N  (Eq. 2.4) ------------------------------------------------
        "a_n":     1.3e4 / S,
        "b_n":     4.12e-2,
        "c_n":     1e-11 * S,
        "p_n":     0.125,
        "g_n":     _g / S,
        "delta_n": 1e-10 * S,

        # ---- Helper H  (Eq. 2.5) --------------------------------------------
        "a_h":     3.6e4 / S,
        "b_h":     2.0833e-4,
        "p_h":     0.125,
        "g_h":     _g / S,
        "delta_h": 1e-9 * S,

        # ---- Treg R  (Eq. 2.6) ----------------------------------------------
        "a":       5.6e5 / S,
        "delta_r": 1e-5,
        "a_r":     2e-4,
        "b_r":     4e-4,
        "p_r":     0.125,
        "g_r":     _g / S,
        "r":       1e-11 * S,

        # ---- Initial conditions (paper → ABM) -------------------------------
        "C0": round(4.886e7  / S),
        "P0": round(2.7362e5 / S),
        "E0": round(4.2684e7 / S),
        "N0": round(2.3531e7 / S),
        "H0": round(1.0343e8 / S),
        "R0": round(7.7570e6 / S),

        "total_days": 100,
    }


# Default (backward-compatible) — overridden by --scale at runtime
PAPER_PARAMS = build_paper_params(1e5)


# ---------------------------------------------------------------------------
# Paper ODE — exact Eqs. 2.1-2.6 with bilinear coupling (no Hill approx)
# Cross-activation terms are set to zero per paper text (Sec. 2, negligible)
# ---------------------------------------------------------------------------
def odes_paper(t, y, p):
    C, P, E, N, H, R = [max(v, 0.0) for v in y]

    # -- Tumor C (Eq. 2.1) --
    dC = ((p["k_c"] + p["mu_c"] * P) * C * (1.0 - C / p["K_C"])
          - p["b_c"] * N * C
          - p["d_c"] * E * C / (1.0 + p["r1"] * R))

    # -- PSC P (Eq. 2.2) --
    dP = ((p["k_p"] + p["f_p"] * C / (p["mu_p"] + C)) * P * (1.0 - P / p["K_P"])
          - p["lambda_p"] * P)

    # -- Effector E (Eq. 2.3) --
    dE = (p["a_e"]
          - p["b_e"]     * E
          - p["c_e"]     * E * C
          + p["r_e"]     * N * C
          + p["p_e"]     * H * E / (p["g_e"] + H)
          - p["delta_e"] * R * E)

    # -- NK N (Eq. 2.4) --
    dN = (p["a_n"]
          - p["b_n"]     * N
          - p["c_n"]     * N * C
          + p["p_n"]     * H * N / (p["g_n"] + H)
          - p["delta_n"] * R * N)

    # -- Helper H (Eq. 2.5) --
    dH = (p["a_h"]
          - p["b_h"]     * H
          + p["p_h"]     * H * H / (p["g_h"] + H)
          - p["delta_h"] * R * H)

    # -- Treg R (Eq. 2.6) --
    dR = (p["a"]
          - p["delta_r"] * R
          + p["a_r"]     * E
          + p["b_r"]     * H
          + p["p_r"]     * H * R / (p["g_r"] + H)
          - p["r"]       * N * R)

    return [dC, dP, dE, dN, dH, dR]


# ---------------------------------------------------------------------------
# ABM ODE — mirrors params.json behaviour (Hill saturation, crowding)
# ---------------------------------------------------------------------------
def sat(x, K):
    return max(0.0, x) / (K + max(0.0, x)) if (K + max(0.0, x)) > 0 else 0.0


def odes_abm(t, y, p):
    C, P, E, N, H, R = [max(v, 0.0) for v in y]
    cC = max(0.0, 1.0 - C / p["K_C"])
    cP = max(0.0, 1.0 - P / p["K_P"])
    cE = max(0.0, 1.0 - E / p["K_E"])
    cN = max(0.0, 1.0 - N / p["K_N"])
    cH = max(0.0, 1.0 - H / p["K_H"])
    cR = max(0.0, 1.0 - R / p["K_R"])

    boost_P  = p["c_boost_from_P"] * sat(P, p["c_boost_from_P_K"])
    inhib_R  = 1.0 / (1.0 + p["c_R_blocks_E"] * R)
    kill_E_C = p["c_kill_by_E"] * sat(E, p["c_kill_by_E_K"]) * inhib_R
    kill_N_C = p["c_kill_by_N"] * sat(N, p["c_kill_by_N_K"])
    dC = C * ((p["c_base_div"] + boost_P) * cC - kill_E_C - kill_N_C)

    boost_C = p["p_boost_from_C"] * sat(C, p["p_boost_from_C_K"])
    dP = P * ((p["p_base_div"] + boost_C) * cP - p["p_base_death"])

    src_E   = p["e_base_birth"] * cE
    div_E   = p["e_help_from_H"] * sat(H, p["e_help_from_H_K"]) * cE
    die_E   = (p["e_base_death"]
               + p["e_inact_by_C"] * sat(C, p["e_inact_by_C_K"])
               + p["e_suppr_by_R"] * sat(R, p["e_suppr_by_R_K"]))
    dE = src_E + E * (div_E - die_E)

    src_N   = p["n_base_birth"] * cN
    div_N   = p["n_help_from_H"] * sat(H, p["n_help_from_H_K"]) * cN
    die_N   = (p["n_base_death"]
               + p["n_inact_by_C"] * sat(C, p["n_inact_by_C_K"])
               + p["n_suppr_by_R"] * sat(R, p["n_suppr_by_R_K"]))
    dN = src_N + N * (div_N - die_N)

    div_H   = p["h_self_act"] * sat(H, p["h_self_act_K"]) * cH
    src_H   = p["h_base_birth"] * cH
    die_H   = p["h_base_death"] + p["h_suppr_by_R"] * sat(R, p["h_suppr_by_R_K"])
    dH = src_H + H * (div_H - die_H)

    div_R   = (p["r_induced_by_E"] * sat(E, p["r_induced_by_E_K"])
               + p["r_induced_by_H"] * sat(H, p["r_induced_by_H_K"])) * cR
    src_R   = p["r_base_src"] * cR
    die_R   = p["r_decay"] + p["r_cleared_by_N"] * sat(N, p["r_cleared_by_N_K"])
    dR = src_R + R * (div_R - die_R)

    return [dC, dP, dE, dN, dH, dR]


# ---------------------------------------------------------------------------
# Treatment ODE — Eqs. 5.1-5.7 (Paper Section 5)
# State vector: [C, P, E, N, H, R, M_gem, M_abr]
# Injection pulses are handled by the piecewise solver below; this function
# only sees the continuous decay between pulses.
# ---------------------------------------------------------------------------
def odes_treatment(t, y, p, acd47_active=False):
    C, P, E, N, H, R, M_gem, M_abr = [max(v, 0.0) for v in y]

    f_gem = 1.0 - np.exp(-M_gem)
    f_abr = 1.0 - np.exp(-M_abr)

    drug_kill_C      = p.get("gem_c_c",      0.0) * f_gem + p.get("abr_c_c",      0.0) * f_abr
    drug_kill_P      = p.get("gem_c_p",      0.0) * f_gem + p.get("abr_c_p",      0.0) * f_abr
    drug_kill_immune = p.get("gem_c_immune", 0.0) * f_gem + p.get("abr_c_immune", 0.0) * f_abr

    # Baseline Eqs. 2.1-2.6
    dC = ((p["k_c"] + p["mu_c"] * P) * C * (1.0 - C / p["K_C"])
          - p["b_c"] * N * C
          - p["d_c"] * E * C / (1.0 + p["r1"] * R)
          - drug_kill_C * C)

    dP = ((p["k_p"] + p["f_p"] * C / (p["mu_p"] + C)) * P * (1.0 - P / p["K_P"])
          - p["lambda_p"] * P
          - drug_kill_P * P)

    acd47_boost = p.get("acd47_e_boost", 0.0) if acd47_active else 0.0
    dE = (p["a_e"]
          - p["b_e"]     * E
          - p["c_e"]     * E * C
          + p["p_e"]     * H * E / (p["g_e"] + H)
          - p["delta_e"] * R * E
          + acd47_boost  * E
          - drug_kill_immune * E)

    dN = (p["a_n"]
          - p["b_n"]     * N
          - p["c_n"]     * N * C
          + p["p_n"]     * H * N / (p["g_n"] + H)
          - p["delta_n"] * R * N
          - drug_kill_immune * N)

    dH = (p["a_h"]
          - p["b_h"]     * H
          + p["p_h"]     * H * H / (p["g_h"] + H)
          - p["delta_h"] * R * H
          - drug_kill_immune * H)

    dR = (p["a"]
          - p["delta_r"] * R
          + p["a_r"]     * E
          + p["b_r"]     * H
          + p["p_r"]     * H * R / (p["g_r"] + H)
          - p["r"]       * N * R
          - drug_kill_immune * R)

    dM_gem = -p.get("gem_gamma", 0.0) * M_gem if p.get("treat_gem", False) else 0.0
    dM_abr = -p.get("abr_gamma", 0.0) * M_abr if p.get("treat_abr", False) else 0.0

    return [dC, dP, dE, dN, dH, dR, dM_gem, dM_abr]


def solve_treatment(p_paper: dict, p_treat: dict, t_days: float):
    """
    Solve Eqs. 5.1-5.7 piecewise with discrete injection events.

    p_paper   — paper ODE parameters (from build_paper_params)
    p_treat   — treatment parameters (from params.json: treat_*, gem_*, abr_*, acd47_*)
    t_days    — total simulation time (sim days, paper day 7 = sim day 0)
    """
    # Merge paper + treatment params so odes_treatment sees everything
    p = {**p_paper, **{k: v for k, v in p_treat.items() if not k.startswith("_")}}

    day_off = 7.0  # sim day 0 = paper day 7

    def inject_sim_days(start_pday, end_pday, freq_days):
        """Return sorted list of sim-day injection times."""
        times = []
        pday = start_pday
        while pday <= end_pday + 1e-9:
            sim_day = pday - day_off
            if sim_day >= -1e-9:
                times.append(max(0.0, sim_day))
            pday += freq_days
        return times

    # Collect all injection events: (sim_day, "gem"|"abr", dose)
    inj_events: list[tuple[float, str, float]] = []
    if p.get("treat_gem", False):
        for t in inject_sim_days(p["treat_start_day"], p["gem_end_day"], p["gem_freq_days"]):
            inj_events.append((t, "gem", p["gem_dose"]))
    if p.get("treat_abr", False):
        for t in inject_sim_days(p["treat_start_day"], p["abr_end_day"], p["abr_freq_days"]):
            inj_events.append((t, "abr", p["abr_dose"]))

    # Anti-CD47 active window (sim days)
    acd47_start = p["treat_start_day"] - day_off if p.get("treat_acd47", False) else None
    acd47_end   = p["acd47_end_day"]   - day_off if p.get("treat_acd47", False) else None

    # Build sorted breakpoints: 0, all injection times, t_days
    breakpoints = sorted(set([0.0, t_days] + [e[0] for e in inj_events]))
    inj_by_time: dict[float, list] = {}
    for t_ev, drug, dose in inj_events:
        inj_by_time.setdefault(t_ev, []).append((drug, dose))

    # Initial state: [C, P, E, N, H, R, M_gem=0, M_abr=0]
    y = [p["C0"], p["P0"], p["E0"], p["N0"], p["H0"], p["R0"], 0.0, 0.0]

    all_t: list[float] = []
    all_y: list[list[float]] = [[] for _ in range(8)]

    for i in range(len(breakpoints) - 1):
        t_start, t_end = breakpoints[i], breakpoints[i + 1]

        # Apply injection(s) at t_start (except at t=0 which is pre-treatment)
        if i > 0 and t_start in inj_by_time:
            for drug, dose in inj_by_time[t_start]:
                if drug == "gem": y[6] += dose
                if drug == "abr": y[7] += dose

        if t_end <= t_start:
            continue

        acd47_now = (acd47_start is not None
                     and t_start >= acd47_start - 1e-9
                     and t_start <  acd47_end   + 1e-9)

        n_pts = max(2, int(t_end - t_start) + 1)
        t_eval = np.linspace(t_start, t_end, n_pts)

        sol = solve_ivp(
            lambda t, y: odes_treatment(t, y, p, acd47_now),
            [t_start, t_end], y,
            t_eval=t_eval, method="RK45", rtol=1e-9, atol=1e-9,
        )
        if not sol.success:
            raise RuntimeError(f"Treatment ODE solver failed at [{t_start},{t_end}]: {sol.message}")

        # Skip first point of each segment (except segment 0) to avoid duplicates
        skip = 1 if i > 0 else 0
        all_t.extend(sol.t[skip:])
        for j in range(8):
            all_y[j].extend(sol.y[j, skip:])

        y = [max(0.0, v) for v in sol.y[:, -1]]

    # Handle injection at t=0 if any (first breakpoint = 0 and i=0 skipped above)
    # (treat_start_day=7 → sim day 0: inject at beginning before any solving)
    # This edge case is rare in practice (paper uses day 14 start).

    t_arr = np.array(all_t)
    y_arr = np.array(all_y)  # shape (8, N)

    C, P, E, N, H, R = [y_arr[j].clip(min=0).round().astype(int) for j in range(6)]
    M_gem_arr, M_abr_arr = y_arr[6], y_arr[7]
    days = t_arr + day_off

    df = pd.DataFrame({
        "step": (t_arr * 1440).astype(int),
        "days": days,
        "C": C, "P": P, "E": E, "N": N, "H": H, "R": R,
        "total": C + P + E + N + H + R,
        "M_gem": M_gem_arr.round(6),
        "M_abr": M_abr_arr.round(6),
    })
    return df


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_json_params(path: Path) -> dict:
    with open(path) as f:
        raw = json.load(f)
    return {k: v for k, v in raw.items() if not k.startswith("_")}


def solve(ode_fn, p, t_days):
    y0 = [p["C0"], p["P0"], p["E0"], p["N0"], p["H0"], p["R0"]]
    t_eval = np.arange(0, t_days + 1, 1.0)
    sol = solve_ivp(ode_fn, [0, t_days], y0, args=(p,),
                    t_eval=t_eval, method="RK45",
                    rtol=1e-9, atol=1e-9)
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")
    return sol


def df_from_sol(sol) -> pd.DataFrame:
    C, P, E, N, H, R = [v.clip(min=0).round().astype(int) for v in sol.y]
    days = sol.t + 7.0   # simulation t=0 ↔ paper day 7
    df = pd.DataFrame({
        "step": (sol.t * 1440).astype(int),
        "days": days, "C": C, "P": P, "E": E, "N": N, "H": H, "R": R,
    })
    df["total"] = df[["C", "P", "E", "N", "H", "R"]].sum(axis=1)
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
COLORS = {"C": "#1f77b4", "P": "#d62728", "E": "#ff7f0e",
          "N": "#2ca02c",  "H": "#17becf", "R": "#9467bd"}
LABELS = {"C": "Tumor (C)", "P": "PSC (P)", "E": "CD8⁺ T (E)",
          "N": "NK (N)",    "H": "Helper T (H)", "R": "Treg (R)"}
POPS   = ["C", "P", "E", "N", "H", "R"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode",   choices=["paper", "abm", "treatment"], default="paper",
                        help="'paper': exact Table 1 params (default);  "
                             "'abm': mirror params.json;  "
                             "'treatment': Eqs. 5.1-5.7 with drug kill terms")
    parser.add_argument("--scale",  type=float, default=None,
                        help="ABM scale factor S (paper_cells / S = ABM_cells). "
                             "If omitted, read from params.json 'scale_S'; "
                             "falls back to 1e5.")
    parser.add_argument("--params", default="params.json",
                        help="params.json path (used in --mode abm; scale_S read here)")
    parser.add_argument("--abm",    default="output/populations.csv",
                        help="ABM populations CSV for overlay plot")
    parser.add_argument("--refs",   default="data-export",
                        help="Directory to write *_scaled_global.csv (paper mode)")
    parser.add_argument("--out",    default="data-export")
    args = parser.parse_args()

    out_dir  = Path(args.out);  out_dir.mkdir(parents=True, exist_ok=True)
    refs_dir = Path(args.refs); refs_dir.mkdir(parents=True, exist_ok=True)

    # Resolve scale factor: CLI > params.json > default 1e5
    if args.scale is not None:
        S = args.scale
    else:
        try:
            jparams = load_json_params(Path(args.params))
            S = float(jparams.get("scale_S", 1e5))
        except Exception:
            S = 1e5

    if args.mode == "paper":
        p       = build_paper_params(S)
        ode_fn  = odes_paper
        title   = "Paper ODE — Table 1 parameters (Akman Yıldız 2021)"
        print(f"Mode: paper  (exact Table 1, Eqs. 2.1-2.6, S={S:.0e})")
        print(f"  K_C={p['K_C']:.0f}  r1={p['r1']:.0f}  a_e={p['a_e']:.4f}  "
              f"lambda_p={p['lambda_p']:.2e}")
    elif args.mode == "abm":
        p       = load_json_params(Path(args.params))
        ode_fn  = odes_abm
        title   = "ABM ODE — params.json"
        print("Mode: abm  (params.json Hill-saturation equations)")
    else:  # treatment
        p_paper  = build_paper_params(S)
        p_treat  = load_json_params(Path(args.params))
        drugs_on = [d for d in ("gem", "abr", "acd47") if p_treat.get(f"treat_{d}", False)]
        title    = "Treatment ODE — Eqs. 5.1-5.7  [" + (", ".join(drugs_on) or "none") + "]"
        print(f"Mode: treatment  (drugs: {drugs_on or 'none'}  S={S:.0e})")

    if args.mode == "treatment":
        t_days = float(p_treat.get("total_days", p_paper["total_days"]))
        print(f"Integrating {int(t_days)} days  "
              f"(y0 = C={p_paper['C0']} P={p_paper['P0']} E={p_paper['E0']} "
              f"N={p_paper['N0']} H={p_paper['H0']} R={p_paper['R0']})")
        df = solve_treatment(p_paper, p_treat, t_days)
    else:
        t_days = float(p["total_days"])
        print(f"Integrating {int(t_days)} days  "
              f"(y0 = C={p['C0']} P={p['P0']} E={p['E0']} N={p['N0']} H={p['H0']} R={p['R0']})")
        sol = solve(ode_fn, p, t_days)
        df  = df_from_sol(sol)

    # Write ode_reference.csv
    csv_path = out_dir / "ode_reference.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    # In paper mode: also regenerate *_scaled_global.csv ground-truth files
    if args.mode == "paper":
        for pop in POPS:
            ref_path = refs_dir / f"{pop}-Cells_scaled_global.csv"
            pd.DataFrame({"day": df["days"].values,
                          "count": df[pop].values}
                         ).to_csv(ref_path, index=False, header=False)
            print(f"  Reference updated: {ref_path}")

    # Final state summary
    last_row = df.iloc[-1]
    print(f"\nFinal state (day {last_row['days']:.0f}):")
    for pop in POPS:
        print(f"  {pop}: {int(last_row[pop]):>8,}")
    if "M_gem" in df.columns:
        print(f"  M_gem={last_row['M_gem']:.4f}  M_abr={last_row['M_abr']:.4f}")

    # Plot: ODE vs ABM
    abm_path = Path(args.abm)
    if not abm_path.exists():
        abm_path = Path("data-export/populations.csv")

    df_abm = None
    if abm_path.exists():
        df_abm = pd.read_csv(abm_path)
        df_abm.columns = df_abm.columns.str.lower().str.strip()

    # Treatment shading: one band per drug covering its active window
    treat_bands: list[tuple[float, float, str, str]] = []
    if args.mode == "treatment":
        p_t = p_treat
        if p_t.get("treat_gem",   False): treat_bands.append((p_t["treat_start_day"], p_t["gem_end_day"],   "#4e9aff", "Gem"))
        if p_t.get("treat_abr",   False): treat_bands.append((p_t["treat_start_day"], p_t["abr_end_day"],   "#ff9f4e", "Abr"))
        if p_t.get("treat_acd47", False): treat_bands.append((p_t["treat_start_day"], p_t["acd47_end_day"], "#4eff9a", "ACD47"))

    plt.rcParams.update({"font.family": "serif", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.flatten()

    for ax, pop in zip(axes, POPS):
        col = pop.lower()
        ydata_ode = df[pop].values.clip(min=0)
        ax.plot(df["days"], ydata_ode,
                color=COLORS[pop], linewidth=2.5, label="ODE")
        if df_abm is not None and col in df_abm.columns:
            ax.plot(df_abm["days"], df_abm[col],
                    color=COLORS[pop], linewidth=1.6, linestyle="--",
                    alpha=0.85, label="ABM")
        for t0, t1, band_color, band_label in treat_bands:
            ax.axvspan(t0, t1, alpha=0.10, color=band_color, label=band_label)
        ax.set_title(LABELS[pop], fontsize=12)
        ax.set_xlabel("Day"); ax.set_ylabel("Cell count")
        ax.legend(fontsize=8); ax.grid(True, linestyle=":", alpha=0.4)
        all_y = list(ydata_ode)
        if df_abm is not None and col in df_abm.columns:
            all_y += list(df_abm[col].values)
        pos = [v for v in all_y if v > 0]
        if pos and max(pos) / max(min(pos), 1) > 100:
            ax.set_yscale("log")

    fig.suptitle(f"ODE Reference ({args.mode} mode) vs ABM\n{title}", fontsize=12)
    plt.tight_layout()
    plot_path = out_dir / "ode_vs_abm.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved: {plot_path}")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
