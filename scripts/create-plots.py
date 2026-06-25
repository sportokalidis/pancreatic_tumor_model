#!/usr/bin/env python3
"""
Publication-quality comparison: ABM output vs ODE reference vs paper reference.

Generates two figures:
  1. comparison_grid.png  — 2×3 panel, one population per axis, all three sources
  2. comparison_combined.png — all 6 populations on one log-scale axis (overview)

Usage:
  python3 data-export/create-plots.py [--out data-export]
"""

import argparse
import os
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
POPS = ["C", "P", "E", "N", "H", "R"]

LABELS = {
    "C": "Tumor (C)",
    "P": "PSC (P)",
    "E": "CD8⁺ T (E)",
    "N": "NK (N)",
    "H": "Helper T (H)",
    "R": "Treg (R)",
}

COLORS = {
    "C": "#1f77b4",   # blue
    "P": "#d62728",   # red
    "E": "#ff7f0e",   # orange
    "N": "#2ca02c",   # green
    "H": "#17becf",   # cyan
    "R": "#9467bd",   # purple
}

plt.rcParams.update({
    "font.family":      "serif",
    "font.size":        11,
    "axes.titlesize":   12,
    "axes.labelsize":   11,
    "legend.fontsize":  9,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "figure.dpi":       150,
    "axes.spines.top":  False,
    "axes.spines.right":False,
})


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------
def load_abm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.lower().str.strip()
    return df


def load_ode(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.lower().str.strip()
    return df


def load_paper_ref(data_dir: Path, scale_s: float = 1e5) -> dict[str, pd.DataFrame]:
    """
    Load per-population paper reference CSVs (no header, col0=day, col1=count).
    Paper references are at S=1e5; scale them to match ABM's scale_S.
    """
    refs = {}
    scale_factor = scale_s / 1e5  # Paper refs are always S=1e5
    for pop in POPS:
        csv = data_dir / f"{pop}-Cells_scaled_global.csv"
        if csv.exists():
            df = pd.read_csv(csv, header=None, names=["day", "count"])
            df["count"] = df["count"] * scale_factor
            refs[pop] = df
    return refs


# ---------------------------------------------------------------------------
# Figure 1: 2×3 grid — one panel per population
# ---------------------------------------------------------------------------
def plot_grid(df_abm, df_ode, paper_refs, out_path: Path):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=False)
    axes = axes.flatten()

    for ax, pop in zip(axes, POPS):
        color = COLORS[pop]
        pop_lc = pop.lower()

        # Paper reference (scatter dots)
        if pop in paper_refs:
            ref = paper_refs[pop]
            ax.scatter(ref["day"], ref["count"],
                       color=color, s=18, zorder=5,
                       alpha=0.7, label="Paper ref.")

        # ODE reference (thick solid)
        if pop_lc in df_ode.columns:
            ax.plot(df_ode["days"], df_ode[pop_lc],
                    color=color, linewidth=2.5, linestyle="-",
                    label="ODE", zorder=4)

        # ABM output (dashed)
        if pop_lc in df_abm.columns:
            ax.plot(df_abm["days"], df_abm[pop_lc],
                    color=color, linewidth=1.8, linestyle="--",
                    alpha=0.9, label="ABM", zorder=3)

        ax.set_title(LABELS[pop])
        ax.set_xlabel("Day")
        ax.set_ylabel("Cell count")
        ax.legend(loc="best")
        ax.grid(True, linestyle=":", alpha=0.4)

        # Log scale if dynamic range > 100×
        all_vals = []
        if pop_lc in df_ode.columns:
            all_vals.extend(df_ode[pop_lc].tolist())
        if pop_lc in df_abm.columns:
            all_vals.extend(df_abm[pop_lc].tolist())
        if pop in paper_refs:
            all_vals.extend(paper_refs[pop]["count"].tolist())
        positive = [v for v in all_vals if v > 0]
        if positive and max(positive) / max(min(positive), 1) > 100:
            ax.set_yscale("log")

    fig.suptitle(
        "Pancreatic Tumor Model — ABM vs ODE vs Paper Reference\n"
        "(Akman Yıldız et al., 2021)",
        fontsize=13, y=1.01
    )

    # Shared legend at the bottom
    handles = [
        mlines.Line2D([], [], color="gray", linewidth=2.5, linestyle="-",  label="ODE (same params)"),
        mlines.Line2D([], [], color="gray", linewidth=1.8, linestyle="--", label="ABM"),
        mlines.Line2D([], [], color="gray", marker="o", linestyle="None",  markersize=5, label="Paper ref."),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.04), frameon=True)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Grid figure saved: {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: all-in-one overview (log scale)
# ---------------------------------------------------------------------------
def plot_combined(df_abm, df_ode, paper_refs, out_path: Path):
    fig, ax = plt.subplots(figsize=(12, 6))

    for pop in POPS:
        color = COLORS[pop]
        pop_lc = pop.lower()
        lbl = LABELS[pop]

        if pop in paper_refs:
            ref = paper_refs[pop]
            ax.scatter(ref["day"], ref["count"],
                       color=color, s=14, alpha=0.6, zorder=5)

        if pop_lc in df_ode.columns:
            ax.plot(df_ode["days"], df_ode[pop_lc],
                    color=color, linewidth=2.0, linestyle="-",
                    label=f"{lbl} ODE", zorder=4)

        if pop_lc in df_abm.columns:
            ax.plot(df_abm["days"], df_abm[pop_lc],
                    color=color, linewidth=1.4, linestyle="--",
                    alpha=0.85, label=f"{lbl} ABM", zorder=3)

    ax.set_yscale("log")
    ax.set_xlabel("Day")
    ax.set_ylabel("Cell count (log scale)")
    ax.set_title(
        "Pancreatic Tumor Model — All Populations Overview\n"
        "Solid: ODE  |  Dashed: ABM  |  Dots: Paper reference",
        fontsize=12
    )
    ax.legend(ncol=3, fontsize=8, loc="upper left")
    ax.grid(True, which="both", linestyle=":", alpha=0.35)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Combined figure saved: {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3: Individual per-cell-type plots
# ---------------------------------------------------------------------------
def plot_per_cell(df_abm, df_ode, paper_refs, out_dir: Path):
    """Generate one figure per population (ABM, ODE, paper ref on same axis)."""
    for pop in POPS:
        fig, ax = plt.subplots(figsize=(10, 6))
        color = COLORS[pop]
        pop_lc = pop.lower()
        label = LABELS[pop]

        # Paper reference (scatter)
        if pop in paper_refs:
            ref = paper_refs[pop]
            ax.scatter(ref["day"], ref["count"],
                       color=color, s=80, alpha=0.6, zorder=5,
                       marker="o", label="Paper reference", edgecolors="black", linewidth=0.5)

        # ODE (thick solid line)
        if pop_lc in df_ode.columns:
            ax.plot(df_ode["days"], df_ode[pop_lc],
                    color=color, linewidth=3.0, linestyle="-",
                    label="ODE (same params)", zorder=4)

        # ABM (dashed line)
        if pop_lc in df_abm.columns:
            ax.plot(df_abm["days"], df_abm[pop_lc],
                    color=color, linewidth=2.0, linestyle="--",
                    alpha=0.85, label="ABM", zorder=3)

        ax.set_xlabel("Day", fontsize=12)
        ax.set_ylabel("Cell count", fontsize=12)
        ax.set_title(label, fontsize=14, fontweight="bold")
        ax.legend(fontsize=11, loc="best", framealpha=0.95)
        ax.grid(True, linestyle=":", alpha=0.4)

        # Log scale if dynamic range > 100×
        all_vals = []
        if pop_lc in df_ode.columns:
            all_vals.extend(df_ode[pop_lc].tolist())
        if pop_lc in df_abm.columns:
            all_vals.extend(df_abm[pop_lc].tolist())
        if pop in paper_refs:
            all_vals.extend(paper_refs[pop]["count"].tolist())
        positive = [v for v in all_vals if v > 0]
        if positive and max(positive) / max(min(positive), 1) > 100:
            ax.set_yscale("log")

        plt.tight_layout()
        out_path = out_dir / f"plot_{pop}.png"
        plt.savefig(out_path, dpi=200, bbox_inches="tight")
        print(f"Per-cell plot saved: {out_path}")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    import json
    parser = argparse.ArgumentParser()
    parser.add_argument("--abm",    default="output/populations.csv",
                        help="ABM populations CSV (default: output/populations.csv)")
    parser.add_argument("--ode",    default="data-export/ode_reference.csv",
                        help="ODE reference CSV")
    parser.add_argument("--refs",   default="data-export",
                        help="Directory containing *-Cells_scaled_global.csv (at S=1e5, auto-scaled)")
    parser.add_argument("--out",    default="data-export",
                        help="Output directory for figures")
    parser.add_argument("--scale-s", type=float, default=None,
                        help="Scale factor (auto-detected from params.json if not provided)")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Try both the given path and the data-export copy
    abm_path = Path(args.abm)
    if not abm_path.exists():
        abm_path = Path("data-export/populations.csv")
    if not abm_path.exists():
        print(f"[ERROR] ABM output not found at {args.abm} or data-export/populations.csv")
        return 1

    ode_path = Path(args.ode)
    if not ode_path.exists():
        print(f"[ERROR] ODE reference not found: {ode_path} — run scripts/ode_reference.py first")
        return 1

    # Auto-detect scale_S from params.json if not provided
    scale_s = args.scale_s or 1e5
    params_file = abm_path.parent / "params.json"
    if params_file.exists():
        try:
            with open(params_file) as f:
                params = json.load(f)
                scale_s = params.get("scale_S", scale_s)
        except:
            pass

    df_abm = load_abm(abm_path)
    df_ode = load_ode(ode_path)
    paper_refs = load_paper_ref(Path(args.refs), scale_s=scale_s)

    print(f"ABM:   {len(df_abm)} rows, days {df_abm['days'].iloc[0]:.0f}–{df_abm['days'].iloc[-1]:.0f}")
    print(f"ODE:   {len(df_ode)} rows, days {df_ode['days'].iloc[0]:.0f}–{df_ode['days'].iloc[-1]:.0f}")
    print(f"Paper refs loaded: {sorted(paper_refs.keys())}")

    plot_grid(df_abm, df_ode, paper_refs, out_dir / "comparison_grid.png")
    plot_combined(df_abm, df_ode, paper_refs, out_dir / "comparison_combined.png")
    plot_per_cell(df_abm, df_ode, paper_refs, out_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
