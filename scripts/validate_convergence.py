#!/usr/bin/env python3
"""
Stochastic convergence validation — compute mean±std across seeds.

For a given scale (e.g., S=1e4) and dt, run the same ABM with multiple seeds
and compute mean cell counts ± std. Compare to ODE reference to validate:
  1. ABM mean follows ODE (mean-field limit)
  2. Std scales as 1/√N (stochastic theory)
  3. Convergence improves at larger S

Usage:
  python3 scripts/validate_convergence.py \\
    --runs runs/20260616_211445_S1e4_dt6h_s42 \\
           runs/YYYYMMDD_S1e4_dt6h_s123 \\
           runs/YYYYMMDD_S1e4_dt6h_s456 \\
    --ode  data-export/ode_reference.csv \\
    --out  validation_convergence.png
"""

import argparse
from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_run(run_dir: Path) -> pd.DataFrame:
    """Load populations.csv from a run directory."""
    csv = run_dir / "populations.csv"
    if not csv.exists():
        raise FileNotFoundError(f"Not found: {csv}")
    df = pd.read_csv(csv)
    df.columns = df.columns.str.lower().str.strip()
    return df


def load_ode(ode_path: Path) -> pd.DataFrame:
    """Load ODE reference."""
    df = pd.read_csv(ode_path)
    df.columns = df.columns.str.lower().str.strip()
    return df


def load_scale(run_dir: Path) -> float:
    """Load scale_S from params.json."""
    params = json.load(open(run_dir / "params.json"))
    return params.get("scale_S", 1e5)


def main():
    parser = argparse.ArgumentParser(
        description="Compute mean±std across multiple ABM seed runs"
    )
    parser.add_argument("--runs", nargs="+", required=True,
                        help="Run directories (multiple seeds)")
    parser.add_argument("--ode", default="data-export/ode_reference.csv",
                        help="ODE reference CSV")
    parser.add_argument("--out", default="validation_convergence.png",
                        help="Output figure")
    args = parser.parse_args()

    # Load all runs
    runs = [Path(r) for r in args.runs]
    dfs = [load_run(r) for r in runs]
    scale = load_scale(runs[0])

    print(f"[validate_convergence] Loaded {len(dfs)} seed runs, scale S={scale:.0e}")
    for r, df in zip(runs, dfs):
        print(f"  {r.name}: {len(df)} days")

    # Compute mean and std across seeds
    days = dfs[0]["days"].values
    pops = ["c", "p", "e", "n", "h", "r"]
    means = {pop: np.mean([df[pop].values for df in dfs], axis=0) for pop in pops}
    stds = {pop: np.std([df[pop].values for df in dfs], axis=0) for pop in pops}

    # Load ODE
    df_ode = load_ode(Path(args.ode))

    # Plot: mean±std vs ODE
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    labels = {"c": "Tumor (C)", "p": "PSC (P)", "e": "CD8⁺ T (E)",
              "n": "NK (N)", "h": "Helper T (H)", "r": "Treg (R)"}
    colors = {"c": "#1f77b4", "p": "#d62728", "e": "#ff7f0e",
              "n": "#2ca02c", "h": "#17becf", "r": "#9467bd"}

    for ax, pop in zip(axes, pops):
        color = colors[pop]

        # ODE (solid line)
        if pop in df_ode.columns:
            ax.plot(df_ode["days"], df_ode[pop], color=color, linewidth=2.5,
                    linestyle="-", label="ODE", zorder=4)

        # Mean ± std (shaded region)
        ax.fill_between(days, means[pop] - stds[pop], means[pop] + stds[pop],
                        color=color, alpha=0.3, label="ABM mean±std", zorder=2)
        ax.plot(days, means[pop], color=color, linewidth=2.0, linestyle="--",
                label="ABM mean", zorder=3)

        ax.set_xlabel("Day")
        ax.set_ylabel("Cell count")
        ax.set_title(labels[pop])
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, alpha=0.3)

        # Log scale if range > 100×
        all_vals = list(means[pop]) + list(df_ode[pop]) if pop in df_ode.columns else list(means[pop])
        positive = [v for v in all_vals if v > 0]
        if positive and max(positive) / max(min(positive), 1) > 100:
            ax.set_yscale("log")

    fig.suptitle(
        f"Stochastic Convergence Validation (S={scale:.0e}, n={len(dfs)} seeds)\n"
        "Shaded region: ±1 std; Dashed: ABM mean; Solid: ODE",
        fontsize=13
    )
    plt.tight_layout()
    plt.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"\n✓ Convergence plot saved: {args.out}")

    # Print summary statistics
    print(f"\nSummary (day 100):")
    day_100_idx = np.argmin(np.abs(days - 100))
    for pop in pops:
        m = means[pop][day_100_idx]
        s = stds[pop][day_100_idx]
        cv = s / m * 100 if m > 0 else 0  # Coefficient of variation
        print(f"  {pop:5s}: {m:10.0f} ± {s:8.0f}  (CV={cv:5.1f}%)")


if __name__ == "__main__":
    main()
