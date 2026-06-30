#!/usr/bin/env python3
"""
Collect ABM outputs from a sensitivity sweep into a single matrix.

Scans <sweep>/runs/<run_id>/populations.csv for every run, linearly interpolates
each output (C,P,E,N,H,R) at each time point (sa_config.TIMEPOINTS, in the
paper-day coordinate = the `days` column), and writes <sweep>/outputs.csv with
one row per run and columns like C_d7, C_d35, ..., R_d100.

Shared by both the OAT and PRCC pipelines.

Usage (inside the Singularity container — needs numpy/pandas):
    python3 scripts/sensitivity/collect.py <sweep_dir>
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import sa_config as sa


def interp_run(pop_csv: Path) -> dict:
    """Return {output_dXX: value} for one run, or None if unreadable/empty."""
    try:
        df = pd.read_csv(pop_csv)
    except Exception as e:  # noqa: BLE001
        print(f"  [WARN] cannot read {pop_csv}: {e}", file=sys.stderr)
        return None
    if df.empty or "days" not in df.columns:
        print(f"  [WARN] {pop_csv}: empty or no 'days' column", file=sys.stderr)
        return None

    days = df["days"].to_numpy(dtype=float)
    res = {}
    for out in sa.OUTPUTS:
        if out not in df.columns:
            continue
        y = df[out].to_numpy(dtype=float)
        for t in sa.TIMEPOINTS:
            # np.interp clamps to the endpoints outside the data range.
            res[f"{out}_d{t}"] = float(np.interp(t, days, y))
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("sweep_dir")
    ap.add_argument("--out", default=None, help="Output CSV (default <sweep>/outputs.csv)")
    args = ap.parse_args()

    sweep = Path(args.sweep_dir).resolve()
    run_root = sweep / "runs"
    if not run_root.is_dir():
        print(f"[ERROR] no runs/ under {sweep}", file=sys.stderr)
        return 1

    rows = []
    missing = []
    for run_dir in sorted(run_root.iterdir()):
        if not run_dir.is_dir():
            continue
        pop = run_dir / "populations.csv"
        if not pop.exists():
            missing.append(run_dir.name)
            continue
        vals = interp_run(pop)
        if vals is None:
            missing.append(run_dir.name)
            continue
        row = {"run_id": run_dir.name}
        row.update(vals)
        rows.append(row)

    if not rows:
        print("[ERROR] no usable runs found", file=sys.stderr)
        return 1

    out_path = Path(args.out) if args.out else sweep / "outputs.csv"
    df = pd.DataFrame(rows).sort_values("run_id")
    df.to_csv(out_path, index=False)

    print(f"Collected : {len(rows)} runs -> {out_path}")
    if missing:
        print(f"Missing   : {len(missing)} run(s) without populations.csv: "
              f"{', '.join(missing[:10])}{' ...' if len(missing) > 10 else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
