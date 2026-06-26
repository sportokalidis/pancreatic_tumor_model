#!/usr/bin/env python3
"""
Post-process validation suite results.

Computes mean±std across seeds, generates convergence validation plots.

Usage:
  python3 scripts/post_validation_suite.py \\
    --suite-dir runs/base/20260625_153000_S1e4_dt24h \\
    --skip-metrics \\
    --skip-plots
"""

import argparse
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Post-process validation suite")
    parser.add_argument("--suite-dir", required=True, help="Suite directory path")
    parser.add_argument("--skip-metrics", action="store_true", help="Skip fit metrics")
    parser.add_argument("--skip-plots", action="store_true", help="Skip plots")
    return parser.parse_args()


def compute_fit_metrics(suite_dir: Path):
    """Compute fit metrics for each run."""
    print(f"\nComputing fit metrics...")
    run_dirs = sorted([d for d in suite_dir.iterdir() if d.is_dir() and "_s" in d.name])

    for run_dir in run_dirs:
        pop_csv = run_dir / "populations.csv"
        if not pop_csv.exists():
            continue

        ode_csv = run_dir / "ode_reference.csv"
        if not ode_csv.exists():
            continue

        print(f"  {run_dir.name}...", end="", flush=True)
        cmd = [
            "python3", "scripts/calc-error.py",
            "--sim", str(pop_csv),
            "--ode-csv", str(ode_csv),
            "--out", str(run_dir)
        ]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode == 0:
            print(" ✓")
        else:
            print(f" ✗ ({result.stderr[:50]})")


def compute_convergence_plots(suite_dir: Path):
    """Compute mean±std convergence plots."""
    print(f"\nComputing convergence plots...")

    run_dirs = sorted([d for d in suite_dir.iterdir() if d.is_dir() and "_s" in d.name])
    print(f"  Found {len(run_dirs)} runs")

    if len(run_dirs) < 2:
        print("  (Need at least 2 runs for convergence plot)")
        return

    # Load all runs
    dfs = []
    for run_dir in run_dirs:
        pop_csv = run_dir / "populations.csv"
        if pop_csv.exists():
            df = pd.read_csv(pop_csv)
            df.columns = df.columns.str.lower().str.strip()
            dfs.append(df)

    if not dfs:
        print("  No population CSV files found")
        return

    # Compute mean and std
    days = dfs[0]["days"].values
    pops = ["c", "p", "e", "n", "h", "r"]
    means = {pop: np.mean([df[pop].values for df in dfs], axis=0) for pop in pops}
    stds = {pop: np.std([df[pop].values for df in dfs], axis=0) for pop in pops}

    # Load ODE reference
    ode_csv = run_dirs[0] / "ode_reference.csv"
    if not ode_csv.exists():
        print(f"  ODE reference not found: {ode_csv}")
        return

    df_ode = pd.read_csv(ode_csv)
    df_ode.columns = df_ode.columns.str.lower().str.strip()

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
        f"Stochastic Convergence Validation ({len(dfs)} seeds)\n"
        "Shaded region: ±1 std; Dashed: ABM mean; Solid: ODE",
        fontsize=13
    )
    plt.tight_layout()

    out_path = suite_dir / "validation_convergence.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"  ✓ Plot saved: {out_path}")
    plt.close(fig)

    # Print summary
    print(f"\nSummary (day 100):")
    day_100_idx = np.argmin(np.abs(days - 100))
    for pop in pops:
        m = means[pop][day_100_idx]
        s = stds[pop][day_100_idx]
        cv = s / m * 100 if m > 0 else 0
        print(f"  {pop:5s}: {m:10.0f} ± {s:8.0f}  (CV={cv:5.1f}%)")


def main():
    args = parse_args()
    suite_dir = Path(args.suite_dir)

    if not suite_dir.exists():
        print(f"ERROR: Suite directory not found: {suite_dir}")
        return 1

    print(f"\n{'='*60}")
    print(f"Post-processing: {suite_dir.name}")
    print(f"{'='*60}")

    if not args.skip_metrics:
        compute_fit_metrics(suite_dir)

    if not args.skip_plots:
        compute_convergence_plots(suite_dir)

    print(f"\n{'='*60}")
    print(f"✓ Post-processing complete")
    print(f"{'='*60}\n")

    return 0


if __name__ == "__main__":
    exit(main())
