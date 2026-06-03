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
import json
import os
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

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


def load_paper_ref(data_dir: Path) -> dict[str, pd.DataFrame]:
    """Load per-population paper reference CSVs (no header, col0=day, col1=count)."""
    refs = {}
    for pop in POPS:
        csv = data_dir / f"{pop}-Cells_scaled_global.csv"
        if csv.exists():
            df = pd.read_csv(csv, header=None, names=["day", "count"])
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
# Figure 3: Treatment comparison — C and E panels with window shading
# ---------------------------------------------------------------------------
def plot_treatment(df_abm, df_ode, params: dict, out_path: Path):
    # Convert sim-day columns to paper-day for display
    abm_pday = df_abm["days"] + 7
    ode_pday = df_ode["days"] + 7

    # Build list of treatment windows (paper days)
    windows = []
    start = params.get("treat_start_day", 14.0)
    if params.get("treat_gem"):
        windows.append(("GEM", start, params.get("gem_end_day", 56.0), "#e74c3c", 0.13))
    if params.get("treat_abr"):
        windows.append(("ABR", start, params.get("abr_end_day", 28.0), "#27ae60", 0.13))
    if params.get("treat_acd47"):
        windows.append(("ACD47", start, params.get("acd47_end_day", 35.0), "#2980b9", 0.10))

    proto_parts = []
    if params.get("treat_gem"):   proto_parts.append("Gemcitabine")
    if params.get("treat_abr"):   proto_parts.append("Abraxane")
    if params.get("treat_acd47"): proto_parts.append("Anti-CD47")
    title = " + ".join(proto_parts) if proto_parts else "Untreated"

    fig, (ax_c, ax_e) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f"Drug Treatment: {title}\nABM vs ODE  |  sim days shown as paper days",
                 fontsize=12)

    for ax, col, col_label, color in [
        (ax_c, "c", "Tumor cells (C)", COLORS["C"]),
        (ax_e, "e", "CD8⁺ T cells (E)", COLORS["E"]),
    ]:
        col_up = col.upper()
        handles = []

        # Treatment windows
        for w_name, w_start, w_end, w_color, w_alpha in windows:
            ax.axvspan(w_start, w_end, color=w_color, alpha=w_alpha, zorder=1)
            handles.append(mpatches.FancyArrow(0, 0, 0, 0,  # invisible placeholder
                                               width=0, head_width=0, head_length=0,
                                               color=w_color, alpha=w_alpha + 0.2,
                                               label=f"{w_name} window"))

        # ODE (solid)
        if col in df_ode.columns:
            line_ode, = ax.plot(ode_pday, df_ode[col], color=color, lw=2.2, zorder=3)
            handles.insert(0, mlines.Line2D([], [], color=color, lw=2.2, label="ODE"))

        # ABM (scatter) — columns are lowercased by load_abm
        if col in df_abm.columns:
            ax.scatter(abm_pday, df_abm[col], color=color, s=12, alpha=0.65, zorder=4)
            handles.insert(1, mlines.Line2D([], [], color=color, marker="o",
                                            linestyle="None", markersize=4, label="ABM"))

        ax.set_xlabel("Paper day")
        ax.set_ylabel(col_label)
        ax.set_title(col_label)
        ax.legend(handles=handles, fontsize=8, loc="best")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Treatment figure saved: {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--abm",       default="output/populations.csv",
                        help="ABM populations CSV")
    parser.add_argument("--ode",       default="data-export/ode_reference.csv",
                        help="ODE reference CSV")
    parser.add_argument("--refs",      default="data-export",
                        help="Directory containing *-Cells_scaled_global.csv")
    parser.add_argument("--out",       default="data-export",
                        help="Output directory for figures")
    parser.add_argument("--treatment", action="store_true",
                        help="Generate treatment comparison plot instead of baseline grid")
    parser.add_argument("--params",    default="",
                        help="params.json path — required with --treatment for window shading")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    abm_path = Path(args.abm)
    if not abm_path.exists():
        abm_path = Path("data-export/populations.csv")
    if not abm_path.exists():
        print(f"[ERROR] ABM output not found at {args.abm} or data-export/populations.csv")
        return 1

    ode_path = Path(args.ode)
    if not ode_path.exists():
        print(f"[ERROR] ODE reference not found: {ode_path}")
        return 1

    df_abm = load_abm(abm_path)
    df_ode = load_ode(ode_path)

    if args.treatment:
        params = {}
        if args.params and Path(args.params).exists():
            with open(args.params) as f:
                params = {k: v for k, v in json.load(f).items() if not k.startswith("_")}
        plot_treatment(df_abm, df_ode, params, out_dir / "treatment_comparison.png")
    else:
        paper_refs = load_paper_ref(Path(args.refs))
        print(f"ABM:   {len(df_abm)} rows, days {df_abm['days'].iloc[0]:.0f}–{df_abm['days'].iloc[-1]:.0f}")
        print(f"ODE:   {len(df_ode)} rows, days {df_ode['days'].iloc[0]:.0f}–{df_ode['days'].iloc[-1]:.0f}")
        print(f"Paper refs loaded: {sorted(paper_refs.keys())}")
        plot_grid(df_abm, df_ode, paper_refs, out_dir / "comparison_grid.png")
        plot_combined(df_abm, df_ode, paper_refs, out_dir / "comparison_combined.png")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
