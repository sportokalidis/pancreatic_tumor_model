#!/usr/bin/env python3
"""
PRCC global sensitivity analysis (Akman Yildiz 2021, Section 4 / Fig. 3).

Reads <sweep>/lhs_samples.csv (the parameter design) and <sweep>/outputs.csv
(from collect.py), joins them on the sample id, and for each output / time point
computes the Partial Rank Correlation Coefficient (PRCC) of every parameter
against the output.

Method: rank-transform all parameters and the output, then use the inverse of
the rank correlation matrix:

    PRCC(x_i, y) = -Cinv[i, y] / sqrt(Cinv[i, i] * Cinv[y, y])

with a two-sided t-test, df = N - 2 - (k-1), k = number of parameters.

Outputs:
    <sweep>/prcc_results.csv                 param, output, timepoint, prcc, p, significant
    <sweep>/plots/prcc_tornado_<OUT>_d<T>.png   for C and E at each time point

Usage (inside the Singularity container — needs numpy/pandas/scipy/matplotlib):
    python3 scripts/sensitivity/analyze_prcc.py <sweep_dir> [--alpha 0.01]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import sa_config as sa

PLOT_OUTPUTS = ["C", "E"]  # paper's primary targets


def prcc_for(X: np.ndarray, y: np.ndarray, names):
    """PRCC of each column of X against y. Returns (names, prcc, pval, n, k)."""
    n = len(y)
    # Rank-transform (average ranks for ties).
    Xr = np.column_stack([stats.rankdata(X[:, j]) for j in range(X.shape[1])])
    yr = stats.rankdata(y)

    # Drop constant columns (zero variance -> undefined correlation).
    keep = [j for j in range(Xr.shape[1]) if np.ptp(Xr[:, j]) > 0]
    kept_names = [names[j] for j in keep]
    Xr = Xr[:, keep]
    k = Xr.shape[1]

    M = np.column_stack([Xr, yr])
    C = np.corrcoef(M, rowvar=False)
    Cinv = np.linalg.pinv(C)

    yi = k  # index of output in the matrix
    prcc, pval = [], []
    df = n - 2 - (k - 1)
    for i in range(k):
        denom = np.sqrt(Cinv[i, i] * Cinv[yi, yi])
        r = -Cinv[i, yi] / denom if denom > 0 else np.nan
        r = float(np.clip(r, -1.0, 1.0)) if np.isfinite(r) else np.nan
        prcc.append(r)
        if np.isfinite(r) and df > 0 and abs(r) < 1.0:
            t = r * np.sqrt(df / (1.0 - r * r))
            p = 2.0 * stats.t.sf(abs(t), df)
        else:
            p = np.nan
        pval.append(p)
    return kept_names, np.array(prcc), np.array(pval), n, k


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("sweep_dir")
    ap.add_argument("--alpha", type=float, default=0.01,
                    help="Significance level for highlighting (default 0.01)")
    args = ap.parse_args()
    sweep = Path(args.sweep_dir).resolve()

    samples = pd.read_csv(sweep / "lhs_samples.csv")
    outs = pd.read_csv(sweep / "outputs.csv").rename(columns={"run_id": "sample_id"})
    df = samples.merge(outs, on="sample_id", how="inner")
    if df.empty:
        print("[ERROR] no overlap between lhs_samples.csv and outputs.csv",
              file=sys.stderr)
        return 1
    print(f"Joined {len(df)} samples (of {len(samples)} designed).")

    param_names = [p for p in sa.PARAM_NAMES if p in df.columns]
    X = df[param_names].to_numpy(dtype=float)

    records = []
    for out in sa.OUTPUTS:
        for t in sa.TIMEPOINTS:
            col = f"{out}_d{t}"
            if col not in df.columns:
                continue
            y = df[col].to_numpy(dtype=float)
            if np.ptp(y) == 0:  # output constant across samples
                continue
            names, prcc, pval, n, k = prcc_for(X, y, param_names)
            for nm, r, p in zip(names, prcc, pval):
                records.append({
                    "param": nm,
                    "paper_symbol": sa.PARAMS.get(nm, ""),
                    "output": out, "timepoint": t,
                    "prcc": r, "pvalue": p,
                    "significant": bool(np.isfinite(p) and p < args.alpha),
                    "n_samples": n, "n_params": k,
                })

    res = pd.DataFrame(records)
    res_path = sweep / "prcc_results.csv"
    res.to_csv(res_path, index=False)
    print(f"PRCC results -> {res_path}")

    plot_dir = sweep / "plots"
    plot_dir.mkdir(exist_ok=True)
    for out in PLOT_OUTPUTS:
        for t in sa.TIMEPOINTS:
            sub = res[(res["output"] == out) & (res["timepoint"] == t)].copy()
            sub = sub.dropna(subset=["prcc"])
            if sub.empty:
                continue
            sub["abs"] = sub["prcc"].abs()
            sub = sub.sort_values("abs")
            labels = [f"{sym}\n({p})" if sym else p
                      for p, sym in zip(sub["param"], sub["paper_symbol"])]
            colors = ["#c0392b" if v > 0 else "#2980b9" for v in sub["prcc"]]
            # Non-significant bars drawn translucent.
            alphas = [1.0 if s else 0.35 for s in sub["significant"]]
            fig, ax = plt.subplots(figsize=(8, max(4, 0.32 * len(sub))))
            bars = ax.barh(range(len(sub)), sub["prcc"], color=colors)
            for b, a in zip(bars, alphas):
                b.set_alpha(a)
            ax.set_yticks(range(len(sub)))
            ax.set_yticklabels(labels, fontsize=7)
            ax.set_xlim(-1, 1)
            ax.axvline(0, color="k", lw=0.8)
            ax.set_xlabel("PRCC")
            ax.set_title(f"PRCC — {out} at day {t}  "
                         f"(N={int(sub['n_samples'].iloc[0])}, "
                         f"solid: p<{args.alpha})")
            fig.tight_layout()
            fig.savefig(plot_dir / f"prcc_tornado_{out}_d{t}.png", dpi=140)
            plt.close(fig)
    print(f"Plots -> {plot_dir}/prcc_tornado_*.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
