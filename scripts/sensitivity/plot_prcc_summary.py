#!/usr/bin/env python3
"""
Presentation figures for the PRCC sensitivity analysis (paper Fig. 3 analog).

Reads <sweep>/prcc_results.csv (from analyze_prcc.py) and produces publication-
ready artifacts that line up directly against Akman Yildiz et al. (2021) Fig. 3:

  plots/prcc_heatmap_C.png     PRCC of C vs top params across time points (shows
  plots/prcc_heatmap_E.png       the temporal shift k_c -> a_c as tumor saturates)
  plots/prcc_combined_C_E.png  2x(time) tornado grid: C (top row), E (bottom row)
  prcc_significant.csv         tidy table of all significant (p<alpha) entries

Usage (inside the Singularity container — needs numpy/pandas/matplotlib):
    python3 scripts/sensitivity/plot_prcc_summary.py <sweep_dir> [--top 15] [--alpha 0.01]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

import sa_config as sa

PLOT_OUTPUTS = ["C", "E"]


def _label(param: str, sym) -> str:
    return f"{param} ({sym})" if isinstance(sym, str) and sym else param


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("sweep_dir")
    ap.add_argument("--top", type=int, default=15,
                    help="Max parameters shown per heatmap (by max |PRCC|)")
    ap.add_argument("--alpha", type=float, default=0.01)
    args = ap.parse_args()
    sweep = Path(args.sweep_dir).resolve()

    res = pd.read_csv(sweep / "prcc_results.csv")
    plot_dir = sweep / "plots"
    plot_dir.mkdir(exist_ok=True)

    # Time points that actually have PRCC values (day 7 = identical ICs -> none).
    tps = sorted(t for t in res["timepoint"].unique()
                 if res[res["timepoint"] == t]["prcc"].notna().any())

    # ---- tidy significant-parameters table ---------------------------------
    sig = res[res["significant"]].copy()
    sig["abs_prcc"] = sig["prcc"].abs()
    sig = sig.sort_values(["output", "timepoint", "abs_prcc"],
                          ascending=[True, True, False])
    sig[["output", "timepoint", "param", "paper_symbol", "prcc", "pvalue"]].to_csv(
        sweep / "prcc_significant.csv", index=False)

    # ---- per-output PRCC-vs-time heatmaps ----------------------------------
    for out in PLOT_OUTPUTS:
        sub = res[res["output"] == out]
        piv = sub.pivot_table(index="param", columns="timepoint", values="prcc")
        piv = piv[[t for t in tps if t in piv.columns]]
        if piv.empty:
            continue
        order = piv.abs().max(axis=1).sort_values(ascending=False).index[:args.top]
        piv = piv.loc[order]
        symmap = dict(zip(sub["param"], sub["paper_symbol"]))
        ylabels = [_label(p, symmap.get(p, "")) for p in piv.index]

        data = piv.values
        fig, ax = plt.subplots(figsize=(1.8 + 1.2 * piv.shape[1],
                                        0.45 * piv.shape[0] + 1.6))
        im = ax.imshow(data, cmap="RdBu_r",
                       norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1), aspect="auto")
        ax.set_xticks(range(piv.shape[1]))
        ax.set_xticklabels([f"day {c}" for c in piv.columns])
        ax.set_yticks(range(piv.shape[0]))
        ax.set_yticklabels(ylabels, fontsize=8)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                v = data[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                            fontsize=7, color="white" if abs(v) > 0.55 else "black")
        ax.set_title(f"PRCC of {out} vs parameters over time (top {piv.shape[0]})")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="PRCC")
        fig.tight_layout()
        fig.savefig(plot_dir / f"prcc_heatmap_{out}.png", dpi=150)
        plt.close(fig)

    # ---- combined C/E tornado grid across time -----------------------------
    fig, axes = plt.subplots(2, len(tps), figsize=(3.5 * len(tps), 7.2),
                             squeeze=False)
    for r, out in enumerate(PLOT_OUTPUTS):
        for c, t in enumerate(tps):
            ax = axes[r][c]
            s = res[(res["output"] == out) & (res["timepoint"] == t)].dropna(
                subset=["prcc"]).copy()
            s["abs_prcc"] = s["prcc"].abs()
            s = s.sort_values("abs_prcc").tail(10)
            symmap = dict(zip(s["param"], s["paper_symbol"]))
            colors = ["#c0392b" if v > 0 else "#2980b9" for v in s["prcc"]]
            bars = ax.barh(range(len(s)), s["prcc"], color=colors)
            for b, ok in zip(bars, s["significant"]):
                b.set_alpha(1.0 if ok else 0.35)
            ax.set_yticks(range(len(s)))
            ax.set_yticklabels([_label(p, symmap.get(p, "")) for p in s["param"]],
                               fontsize=7)
            ax.set_xlim(-1, 1)
            ax.axvline(0, color="k", lw=0.6)
            if r == 0:
                ax.set_title(f"day {t}")
            if c == 0:
                ax.set_ylabel(out, fontsize=13, fontweight="bold")
    fig.suptitle(f"PRCC tornado — C (top) and E (bottom) across time "
                 f"(red +, blue −; solid: p<{args.alpha:g})")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(plot_dir / "prcc_combined_C_E.png", dpi=150)
    plt.close(fig)

    print(f"Wrote: {plot_dir}/prcc_heatmap_C.png, prcc_heatmap_E.png, "
          f"prcc_combined_C_E.png")
    print(f"Wrote: {sweep}/prcc_significant.csv ({len(sig)} significant entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
