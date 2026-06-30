#!/usr/bin/env python3
"""
Analyze the One-At-a-Time (OAT) sweep: normalized local sensitivities + tornado.

Reads <sweep>/oat_plan.csv and <sweep>/outputs.csv (produced by collect.py) and
computes, for each parameter / output / time point, the normalized sensitivity
(elasticity):

    S = (Delta_output / output_baseline) / (Delta_param / param_baseline)

For each parameter both endpoints (x0.5, x2.0) are used and averaged into a
single elasticity (sign preserved); the per-endpoint values are kept in the CSV.

Outputs:
    <sweep>/oat_sensitivity.csv          param, output, timepoint, elasticity, ...
    <sweep>/plots/oat_tornado_<OUT>_d<T>.png   for C and E at each time point

Usage (inside the Singularity container — needs numpy/pandas/matplotlib):
    python3 scripts/sensitivity/analyze_oat.py <sweep_dir>
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import sa_config as sa

PLOT_OUTPUTS = ["C", "E"]  # paper's primary targets


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("sweep_dir")
    args = ap.parse_args()
    sweep = Path(args.sweep_dir).resolve()

    plan = pd.read_csv(sweep / "oat_plan.csv")
    outs = pd.read_csv(sweep / "outputs.csv").set_index("run_id")

    base_run = "oat_0000"
    if base_run not in outs.index:
        print(f"[ERROR] baseline run {base_run} missing from outputs.csv",
              file=sys.stderr)
        return 1
    base_out = outs.loc[base_run]

    out_cols = [c for c in outs.columns]
    records = []

    perturb = plan[plan["param"] != ""]
    for param, grp in perturb.groupby("param"):
        for out in sa.OUTPUTS:
            for t in sa.TIMEPOINTS:
                col = f"{out}_d{t}"
                if col not in out_cols:
                    continue
                y0 = base_out[col]
                elas_endpoints = []
                rec = {"param": param,
                       "paper_symbol": sa.PARAMS.get(param, ""),
                       "output": out, "timepoint": t,
                       "baseline_value": y0}
                for _, r in grp.iterrows():
                    rid, factor = r["run_id"], float(r["factor"])
                    if rid not in outs.index:
                        continue
                    yf = outs.loc[rid, col]
                    dp = factor - 1.0
                    if y0 != 0 and dp != 0:
                        elas = ((yf - y0) / y0) / dp
                    else:
                        elas = np.nan
                    elas_endpoints.append(elas)
                    rec[f"out_x{factor:g}"] = yf
                    rec[f"elasticity_x{factor:g}"] = elas
                finite = [e for e in elas_endpoints if np.isfinite(e)]
                rec["elasticity"] = float(np.mean(finite)) if finite else np.nan
                records.append(rec)

    res = pd.DataFrame(records)
    res_path = sweep / "oat_sensitivity.csv"
    res.to_csv(res_path, index=False)
    print(f"Sensitivities -> {res_path}")

    plot_dir = sweep / "plots"
    plot_dir.mkdir(exist_ok=True)
    for out in PLOT_OUTPUTS:
        for t in sa.TIMEPOINTS:
            sub = res[(res["output"] == out) & (res["timepoint"] == t)].copy()
            sub = sub.dropna(subset=["elasticity"])
            if sub.empty:
                continue
            sub["abs"] = sub["elasticity"].abs()
            sub = sub.sort_values("abs")
            labels = [f"{sym}\n({p})" if sym else p
                      for p, sym in zip(sub["param"], sub["paper_symbol"])]
            colors = ["#c0392b" if v > 0 else "#2980b9" for v in sub["elasticity"]]
            fig, ax = plt.subplots(figsize=(8, max(4, 0.32 * len(sub))))
            ax.barh(range(len(sub)), sub["elasticity"], color=colors)
            ax.set_yticks(range(len(sub)))
            ax.set_yticklabels(labels, fontsize=7)
            ax.axvline(0, color="k", lw=0.8)
            ax.set_xlabel("Normalized sensitivity (elasticity)")
            ax.set_title(f"OAT tornado — {out} at day {t}")
            fig.tight_layout()
            fig.savefig(plot_dir / f"oat_tornado_{out}_d{t}.png", dpi=140)
            plt.close(fig)
    print(f"Plots -> {plot_dir}/oat_tornado_*.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
