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

    # Use only the common time range (all runs may have different lengths if stopped early)
    min_len = min(len(df) for df in dfs)
    dfs = [df.iloc[:min_len] for df in dfs]

    # Compute mean and std
    days = dfs[0]["days"].values
    pops = ["c", "p", "e", "n", "h", "r"]
    means = {pop: np.mean([df[pop].values for df in dfs], axis=0) for pop in pops}
    stds = {pop: np.std([df[pop].values for df in dfs], axis=0) for pop in pops}

    # Try to auto-detect scale from first run's params.json
    scale_s = 1e5
    params_file = run_dirs[0] / "params.json"
    if params_file.exists():
        try:
            import json
            with open(params_file) as f:
                params = json.load(f)
                scale_s = params.get("scale_S", 1e5)
        except:
            pass

    # Load digitized paper reference (scaled to ABM's scale)
    scale_str = f"S{scale_s:.0e}".replace("+0", "").replace("e", "e")
    digitized_dir = Path("data-export/digitized")
    df_ref = {}
    has_ref = False

    for pop in pops:
        digitized_csv = digitized_dir / f"{pop.upper()}-Cells-{scale_str}.csv"
        if digitized_csv.exists():
            df = pd.read_csv(digitized_csv, header=None, names=["days", pop])
            df_ref[pop] = df
            has_ref = True

    if not has_ref:
        print(f"  Digitized references not found in {digitized_dir}")
        print(f"  (Looked for *-Cells-{scale_str}.csv)")
        return

    # Load ODE reference from first run
    df_ode = {}
    ode_csv = run_dirs[0] / "ode_reference.csv"
    if ode_csv.exists():
        ode_df = pd.read_csv(ode_csv)
        ode_df.columns = ode_df.columns.str.lower().str.strip()
        for pop in pops:
            if pop in ode_df.columns:
                df_ode[pop] = ode_df[["days", pop]].copy()

    labels = {"c": "Tumor (C)", "p": "PSC (P)", "e": "CD8⁺ T (E)",
              "n": "NK (N)", "h": "Helper T (H)", "r": "Treg (R)"}
    colors = {"c": "#1f77b4", "p": "#d62728", "e": "#ff7f0e",
              "n": "#2ca02c", "h": "#17becf", "r": "#9467bd"}

    def plot_cell_type(pop, ax, use_log=False):
        """Plot a single cell type."""
        color = colors[pop]

        # ODE reference (dotted line) — theoretical solution
        if pop in df_ode:
            ode = df_ode[pop]
            ax.plot(ode["days"], ode[pop], color=color, linewidth=3.0,
                    linestyle=":", label="ODE reference (theory)", zorder=5, alpha=0.8)

        # Paper reference (solid line) — DISABLED
        # if pop in df_ref:
        #     ref = df_ref[pop]
        #     ax.plot(ref["days"], ref[pop], color=color, linewidth=2.5,
        #             linestyle="-", label="Paper ref. (digitized)", zorder=4)

        # Mean ± std (shaded region)
        ax.fill_between(days, means[pop] - stds[pop], means[pop] + stds[pop],
                        color=color, alpha=0.3, label="ABM mean±std", zorder=2)
        ax.plot(days, means[pop], color=color, linewidth=1.0, linestyle="--",
                label="ABM mean", zorder=3)

        ax.set_xlabel("Day")
        ax.set_ylabel("Cell count")
        ax.set_title(labels[pop])
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, alpha=0.3, which="both")

        if use_log:
            ax.set_yscale("log")

    # 1. Standard 2x3 plot (linear scale with auto-log where needed)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for ax, pop in zip(axes, pops):
        plot_cell_type(pop, ax, use_log=False)
        # Auto-log if range > 100×
        all_vals = list(means[pop])
        if pop in df_ref:
            all_vals.extend(df_ref[pop][pop].tolist())
        positive = [v for v in all_vals if v > 0]
        if positive and max(positive) / max(min(positive), 1) > 100:
            ax.set_yscale("log")

    fig.suptitle(
        f"Stochastic Convergence Validation ({len(dfs)} seeds)\n"
        "Shaded region: ±1 std; Dashed: ABM mean; Dotted: ODE reference",
        fontsize=13
    )
    plt.tight_layout()

    out_path = suite_dir / "validation_convergence.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"  ✓ Plot saved: {out_path}")
    plt.close(fig)

    # 2. Individual plots for each cell type
    for pop in pops:
        fig, ax = plt.subplots(figsize=(10, 6))
        plot_cell_type(pop, ax, use_log=False)
        # Auto-log if needed
        all_vals = list(means[pop])
        if pop in df_ref:
            all_vals.extend(df_ref[pop][pop].tolist())
        positive = [v for v in all_vals if v > 0]
        if positive and max(positive) / max(min(positive), 1) > 100:
            ax.set_yscale("log")
        fig.suptitle(f"{labels[pop]} — Stochastic Convergence ({len(dfs)} seeds)", fontsize=14)
        plt.tight_layout()
        out_path = suite_dir / f"validation_{pop}_individual.png"
        plt.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)

    # 3. 2x3 plot with logarithmic y-axis for all subplots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for ax, pop in zip(axes, pops):
        plot_cell_type(pop, ax, use_log=True)

    fig.suptitle(
        f"Stochastic Convergence Validation — Logarithmic Scale ({len(dfs)} seeds)\n"
        "Shaded region: ±1 std; Dashed: ABM mean; Dotted: ODE reference",
        fontsize=13
    )
    plt.tight_layout()

    out_path = suite_dir / "validation_convergence_log.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"  ✓ Plot saved: {out_path}")
    plt.close(fig)

    # 4. Combined plot (standard linear with better aesthetics)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for ax, pop in zip(axes, pops):
        plot_cell_type(pop, ax, use_log=False)
        # Auto-log if range > 100×
        all_vals = list(means[pop])
        if pop in df_ref:
            all_vals.extend(df_ref[pop][pop].tolist())
        positive = [v for v in all_vals if v > 0]
        if positive and max(positive) / max(min(positive), 1) > 100:
            ax.set_yscale("log")

    fig.suptitle(
        f"Stochastic Convergence Validation — Combined View ({len(dfs)} seeds)\n"
        "Shaded region: ±1 std; Dashed: ABM mean; Dotted: ODE reference",
        fontsize=13
    )
    plt.tight_layout()

    out_path = suite_dir / "validation_convergence_combined.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"  ✓ Plot saved: {out_path}")
    plt.close(fig)

    print(f"  ✓ Individual plots saved for each cell type")

    # 5. All cell types on same axis (linear scale)
    fig, ax = plt.subplots(figsize=(14, 8))

    for pop in pops:
        color = colors[pop]

        # ODE reference (dotted line)
        if pop in df_ode:
            ode = df_ode[pop]
            ax.plot(ode["days"], ode[pop], color=color, linewidth=3.0,
                    linestyle=":", alpha=0.7, zorder=4)

        # Paper reference (solid line)
        if pop in df_ref:
            ref = df_ref[pop]
            ax.plot(ref["days"], ref[pop], color=color, linewidth=2.5,
                    linestyle="-", alpha=0.7, zorder=4.5)

        # Mean ± std (shaded region)
        ax.fill_between(days, means[pop] - stds[pop], means[pop] + stds[pop],
                        color=color, alpha=0.15, zorder=2)
        ax.plot(days, means[pop], color=color, linewidth=1.0, linestyle="--",
                label=labels[pop], zorder=3)

    ax.set_xlabel("Day", fontsize=12)
    ax.set_ylabel("Cell count", fontsize=12)
    ax.set_title(f"All Cell Types — Stochastic Convergence ({len(dfs)} seeds)", fontsize=14)
    ax.legend(loc="best", fontsize=11, ncol=3)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = suite_dir / "validation_all_types_linear.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"  ✓ Plot saved: {out_path}")
    plt.close(fig)

    # 6. All cell types on same axis (logarithmic scale)
    fig, ax = plt.subplots(figsize=(14, 8))

    for pop in pops:
        color = colors[pop]

        # ODE reference (dotted line)
        if pop in df_ode:
            ode = df_ode[pop]
            ax.plot(ode["days"], ode[pop], color=color, linewidth=3.0,
                    linestyle=":", alpha=0.7, zorder=4)

        # Paper reference (solid line)
        if pop in df_ref:
            ref = df_ref[pop]
            ax.plot(ref["days"], ref[pop], color=color, linewidth=2.5,
                    linestyle="-", alpha=0.7, zorder=4.5)

        # Mean ± std (shaded region)
        ax.fill_between(days, means[pop] - stds[pop], means[pop] + stds[pop],
                        color=color, alpha=0.15, zorder=2)
        ax.plot(days, means[pop], color=color, linewidth=1.0, linestyle="--",
                label=labels[pop], zorder=3)

    ax.set_xlabel("Day", fontsize=12)
    ax.set_ylabel("Cell count", fontsize=12)
    ax.set_title(f"All Cell Types — Stochastic Convergence ({len(dfs)} seeds) [Log Scale]", fontsize=14)
    ax.set_yscale("log")
    ax.legend(loc="best", fontsize=11, ncol=3)
    ax.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    out_path = suite_dir / "validation_all_types_log.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"  ✓ Plot saved: {out_path}")
    plt.close(fig)

    # 7. Min-max envelope with mean (linear scale)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for ax, pop in zip(axes, pops):
        color = colors[pop]

        # Min-Max envelope
        mins = np.min([df[pop].values for df in dfs], axis=0)
        maxs = np.max([df[pop].values for df in dfs], axis=0)
        ax.fill_between(days, mins, maxs, color=color, alpha=0.25,
                        label="Min-Max range", zorder=2)

        # ODE reference (dotted line)
        if pop in df_ode:
            ode = df_ode[pop]
            ax.plot(ode["days"], ode[pop], color=color, linewidth=2.0,
                    linestyle=":", label="ODE reference", zorder=5, alpha=0.8)

        # Paper reference (solid line) — DISABLED
        # if pop in df_ref:
        #     ref = df_ref[pop]
        #     ax.plot(ref["days"], ref[pop], color=color, linewidth=2.5,
        #             linestyle="-", label="Paper ref. (digitized)", zorder=4)

        # Mean line
        ax.plot(days, means[pop], color=color, linewidth=2.5, linestyle="--",
                label="ABM mean", zorder=3)

        ax.set_xlabel("Day")
        ax.set_ylabel("Cell count")
        ax.set_title(labels[pop])
        ax.legend(loc="best", fontsize=8)
        ax.grid(True, alpha=0.3)

        # Auto-log if needed
        all_vals = list(means[pop])
        all_vals.extend(mins)
        all_vals.extend(maxs)
        if pop in df_ref:
            all_vals.extend(df_ref[pop][pop].tolist())
        positive = [v for v in all_vals if v > 0]
        if positive and max(positive) / max(min(positive), 1) > 100:
            ax.set_yscale("log")

    fig.suptitle(
        f"Min-Max Envelope + ABM Mean ({len(dfs)} seeds)\n"
        "Shaded: min-max range across all runs; Dashed: ABM mean",
        fontsize=13
    )
    plt.tight_layout()

    out_path = suite_dir / "validation_all_runs_minmax.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"  ✓ Plot saved: {out_path}")
    plt.close(fig)

    # 8. All cell types on same axis with min-max envelope (linear scale)
    fig, ax = plt.subplots(figsize=(14, 8))

    for pop in pops:
        color = colors[pop]

        # Min-Max envelope
        mins = np.min([df[pop].values for df in dfs], axis=0)
        maxs = np.max([df[pop].values for df in dfs], axis=0)
        ax.fill_between(days, mins, maxs, color=color, alpha=0.2, zorder=2)

        # ODE reference (dotted line)
        if pop in df_ode:
            ode = df_ode[pop]
            ax.plot(ode["days"], ode[pop], color=color, linewidth=3.0,
                    linestyle=":", alpha=0.7, zorder=4)

        # Paper reference (solid line)
        if pop in df_ref:
            ref = df_ref[pop]
            ax.plot(ref["days"], ref[pop], color=color, linewidth=2.5,
                    linestyle="-", alpha=0.7, zorder=4.5)

        # Mean line
        ax.plot(days, means[pop], color=color, linewidth=1.0, linestyle="--",
                label=labels[pop], zorder=3)

    ax.set_xlabel("Day", fontsize=12)
    ax.set_ylabel("Cell count", fontsize=12)
    ax.set_title(f"All Cell Types — Min-Max Envelope + Mean ({len(dfs)} seeds)", fontsize=14)
    ax.legend(loc="best", fontsize=11, ncol=3)
    ax.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    out_path = suite_dir / "validation_all_runs_combined.png"
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
