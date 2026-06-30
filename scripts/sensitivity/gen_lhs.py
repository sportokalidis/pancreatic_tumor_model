#!/usr/bin/env python3
"""
Generate Latin-Hypercube samples for the global (PRCC) sensitivity analysis.

Replicates Akman Yildiz et al. (2021) Section 4: sample all Table 1 parameters
simultaneously over [0.5x, 2x] of baseline using Latin-Hypercube sampling, run
each sample, then compute PRCC (see analyze_prcc.py).

Sampling is log-uniform over the multiplicative range by default: the factor f
is drawn uniformly in [ln(0.5), ln(2.0)], so the range is symmetric about the
baseline (f=1). Use --dist linear for uniform sampling in [0.5x, 2x] linear space.

LHS is implemented with numpy only (stratified one-point-per-bin + independent
per-dimension permutation) so it does not depend on scipy.stats.qmc.

Usage (inside the Singularity container — needs numpy):
    python3 scripts/sensitivity/gen_lhs.py [--samples 250] [--seed 12345] [--out-dir DIR]

Outputs (under runs/sensitivity/<timestamp>_lhs/ unless --out-dir given):
    configs/sample_0000.json ...   one JSON per sample
    lhs_samples.csv                sample_id + the 35 sampled parameter values
"""
import argparse
import copy
import csv
import datetime as dt
import json
import math
from pathlib import Path

import numpy as np

import sa_config as sa


def latin_hypercube(n: int, d: int, rng: np.random.RandomState) -> np.ndarray:
    """Return an (n, d) array in [0,1) — one sample per stratum per dimension."""
    out = np.empty((n, d))
    cut = np.arange(n)
    for j in range(d):
        # one uniform point inside each of the n equal strata, then shuffle
        u = (cut + rng.random_sample(n)) / n
        out[:, j] = u[rng.permutation(n)]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=250)
    ap.add_argument("--seed", type=int, default=12345,
                    help="RNG seed for the LHS design (reproducibility)")
    ap.add_argument("--dist", choices=["log", "linear"], default="log",
                    help="Sampling distribution over the [0.5x, 2x] range")
    ap.add_argument("--out-dir", default=None,
                    help="Sweep directory (default: runs/sensitivity/<ts>_lhs)")
    args = ap.parse_args()

    with open(sa.BASELINE_CONFIG) as f:
        base = json.load(f)

    for p in sa.PARAM_NAMES:
        if p not in base:
            raise KeyError(f"Param '{p}' not present in baseline config "
                           f"{sa.BASELINE_CONFIG}")

    names = sa.PARAM_NAMES
    d = len(names)
    lo, hi = sa.RANGE

    rng = np.random.RandomState(args.seed)
    unit = latin_hypercube(args.samples, d, rng)  # (n, d) in [0,1)

    # Map unit cube -> per-parameter multiplicative factor in [lo, hi].
    if args.dist == "log":
        log_lo, log_hi = math.log(lo), math.log(hi)
        factors = np.exp(log_lo + unit * (log_hi - log_lo))
    else:
        factors = lo + unit * (hi - lo)

    base0 = np.array([float(base[p]) for p in names])
    values = factors * base0  # (n, d) actual parameter values

    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    sweep = Path(args.out_dir) if args.out_dir else sa.runs_dir() / f"{stamp}_lhs"
    cfg_dir = sweep / "configs"
    run_root = sweep / "runs"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    run_root.mkdir(parents=True, exist_ok=True)

    rows = []
    for i in range(args.samples):
        sample_id = f"sample_{i:04d}"
        cfg = copy.deepcopy(base)
        for j, p in enumerate(names):
            cfg[p] = float(values[i, j])
        out = run_root / sample_id
        cfg["output_dir"] = str(out)
        with open(cfg_dir / f"{sample_id}.json", "w") as f:
            json.dump(cfg, f, indent=2)
        row = {"sample_id": sample_id}
        row.update({p: float(values[i, j]) for j, p in enumerate(names)})
        rows.append(row)

    samples_path = sweep / "lhs_samples.csv"
    with open(samples_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample_id"] + names)
        w.writeheader()
        w.writerows(rows)

    print(f"Sweep dir : {sweep}")
    print(f"Samples   : {args.samples}  params: {d}  dist: {args.dist}  "
          f"design-seed: {args.seed}")
    print(f"Samples   : {samples_path}")
    print(f"Array size: --array=0-{args.samples - 1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
