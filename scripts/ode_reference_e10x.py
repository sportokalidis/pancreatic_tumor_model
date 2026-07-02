#!/usr/bin/env python3
"""
ODE reference solver with custom initial conditions (10x E cells).

Integrates the Paper ODE (Eqs. 2.1-2.6) with modified initial conditions
from params_S1e4_E10x.json for sensitivity analysis.

Usage:
  python3 scripts/ode_reference_e10x.py --scale 1e4 --out-dir OUTPUT_DIR
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# ============================================================================
# Paper ODE parameters (same as ode_reference.py)
# ============================================================================

_k_c  = 0.0886933
_a_c  = 1.14077e-9
_a_p  = 2.0 * _a_c
_mu_c = 0.1 * _k_c * _a_p
_g    = 6.7449e10


def build_paper_params(S: float) -> dict:
    """Build paper parameters scaled to scale factor S."""
    return {
        "k_c":    _k_c,
        "mu_c":   _mu_c * S,
        "K_C":    1.0 / (_a_c * S),
        "b_c":    7.13e-11 * S,
        "d_c":    1.3e-4 * S,
        "r1":     0.345 * S,

        "k_p":      0.1 * _k_c,
        "f_p":      0.154955,
        "mu_p":     5.6e7 / S,
        "K_P":      1.0 / (_a_p * S),
        "lambda_p": 7.83296e-10,

        "a_e":     1.3e4 / S,
        "b_e":     2.7e-2,
        "c_e":     3.42e-10 * S,
        "r_e":     0.0,
        "p_e":     0.125,
        "g_e":     _g / S,
        "f_e":     0.4167,
        "delta_e": 1e-10 * S,

        "a_n":     1.3e4 / S,
        "b_n":     4.12e-2,
        "c_n":     1e-11 * S,
        "p_n":     0.125,
        "g_n":     _g / S,
        "f_n":     0.4167,
        "delta_n": 1e-10 * S,

        "a_h":     3.6e4 / S,
        "b_h":     2.0833e-4,
        "p_h":     0.125,
        "g_h":     _g / S,
        "f_h":     0.4167,
        "delta_h": 1e-9 * S,

        "a":       5.6e5 / S,
        "delta_r": 1e-5,
        "a_r":     2e-4,
        "b_r":     4e-4,
        "p_r":     0.125,
        "g_r":     _g / S,
        "r":       1e-11 * S,

        "beta_E":  4.4691e-13,
        "beta_N":  4.4691e-13,
        "beta_H":  4.4691e-13,

        "total_days": 100,
    }


def odes_paper(t, y, p):
    """Paper ODE equations (Eqs. 2.1-2.6) with f_hat terms."""
    C, P, E, N, H, R = [max(v, 0.0) for v in y]

    dC = ((p["k_c"] + p["mu_c"] * P) * C * (1.0 - C / p["K_C"])
          - p["b_c"] * N * C
          - p["d_c"] * E * C / (1.0 + p["r1"] * R))

    dP = ((p["k_p"] + p["f_p"] * C / (p["mu_p"] + C)) * P * (1.0 - P / p["K_P"])
          - p["lambda_p"] * P)

    dE = (p["a_e"]
          - p["b_e"]     * E
          - p["c_e"]     * E * C
          + p["r_e"]     * N * C
          + p["p_e"]     * H * E / (p["g_e"] + H)
          + p["f_e"]     * (p["beta_E"] * E + p["beta_N"] * N + p["beta_H"] * H) * E
          - p["delta_e"] * R * E)

    dN = (p["a_n"]
          - p["b_n"]     * N
          - p["c_n"]     * N * C
          + p["p_n"]     * H * N / (p["g_n"] + H)
          + p["f_n"]     * (p["beta_E"] * E + p["beta_N"] * N + p["beta_H"] * H) * N
          - p["delta_n"] * R * N)

    dH = (p["a_h"]
          - p["b_h"]     * H
          + p["p_h"]     * H * H / (p["g_h"] + H)
          + p["f_h"]     * (p["beta_E"] * E + p["beta_N"] * N + p["beta_H"] * H) * H
          - p["delta_h"] * R * H)

    dR = (p["a"]
          - p["delta_r"] * R
          + p["a_r"]     * E
          + p["b_r"]     * H
          + p["p_r"]     * H * R / (p["g_r"] + H)
          - p["r"]       * N * R)

    return [dC, dP, dE, dN, dH, dR]


def solve(ode_fn, p, y0, t_days):
    """Solve ODE with custom initial conditions."""
    t_eval = np.arange(0, t_days + 1, 1.0)
    sol = solve_ivp(ode_fn, [0, t_days], y0, args=(p,),
                    t_eval=t_eval, method="RK45",
                    rtol=1e-9, atol=1e-9)
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")
    return sol


def df_from_sol(sol) -> pd.DataFrame:
    """Convert solution to DataFrame."""
    C, P, E, N, H, R = [v.clip(min=0).round().astype(int) for v in sol.y]
    days = sol.t + 7.0
    df = pd.DataFrame({
        "step": (sol.t * 1440).astype(int),
        "days": days, "C": C, "P": P, "E": E, "N": N, "H": H, "R": R,
    })
    df["total"] = df[["C", "P", "E", "N", "H", "R"]].sum(axis=1)
    return df


def main():
    parser = argparse.ArgumentParser(description="ODE reference with 10x E initial conditions")
    parser.add_argument("--scale", type=float, default=1e4, help="Scale factor (default 1e4)")
    parser.add_argument("--params", default="configs/params_S1e4_E10x.json",
                        help="Custom params file with modified initial conditions")
    parser.add_argument("--out-dir", default="data-export", help="Output directory")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load custom params to get initial conditions
    with open(args.params) as f:
        custom_params = json.load(f)

    # Extract initial conditions
    C0 = custom_params.get("C0", 4890)
    P0 = custom_params.get("P0", 30)
    E0 = custom_params.get("E0", 42700)
    N0 = custom_params.get("N0", 2350)
    H0 = custom_params.get("H0", 10340)
    R0 = custom_params.get("R0", 780)
    y0 = [C0, P0, E0, N0, H0, R0]

    # Build paper parameters
    p = build_paper_params(args.scale)

    # Solve ODE
    print(f"ODE Reference — 10x E cells sensitivity analysis")
    print(f"  Scale: S={args.scale:.0e}")
    print(f"  Initial conditions (from {args.params}):")
    print(f"    C0={C0}  P0={P0}  E0={E0}  N0={N0}  H0={H0}  R0={R0}")

    t_days = float(p["total_days"])
    sol = solve(odes_paper, p, y0, t_days)
    df = df_from_sol(sol)

    # Save result
    csv_path = out_dir / "ode_reference_e10x.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n✓ Saved: {csv_path}")

    # Print summary at day 100
    df_final = df[df["days"] >= 99.9].iloc[0] if len(df) > 0 else None
    if df_final is not None:
        print(f"\nSummary at day {df_final['days']:.1f}:")
        print(f"  C={df_final['C']:.0f}  P={df_final['P']:.0f}  E={df_final['E']:.0f}  "
              f"N={df_final['N']:.0f}  H={df_final['H']:.0f}  R={df_final['R']:.0f}")


if __name__ == "__main__":
    main()
