#!/usr/bin/env python3
"""
Generate One-At-a-Time (OAT) parameter variants for the sensitivity analysis.

For each Table 1 parameter (see sa_config.PARAMS), emit two config variants
(baseline x 0.5 and baseline x 2.0) plus a single unperturbed baseline run.
This is the cheap first pass (~2N+1 runs) before the full PRCC/LHS sweep — it
gives a local tornado plot and validates the whole pipeline at low cost.

Each config inherits the full baseline (S=1e4, dt=24h, seed=42) and overwrites
exactly one parameter; output_dir points at a per-run folder inside the sweep.

Usage (inside the Singularity container or any Python 3.7+):
    python3 scripts/sensitivity/gen_oat.py [--out-dir DIR]

Outputs (under runs/sensitivity/<timestamp>_oat/ unless --out-dir given):
    configs/oat_0000.json ...     one JSON per run
    oat_plan.csv                  run_id, param, paper_symbol, factor, value, config, output_dir
"""
import argparse
import copy
import csv
import datetime as dt
import json
from pathlib import Path

import sa_config as sa


def load_baseline() -> dict:
    with open(sa.BASELINE_CONFIG) as f:
        return json.load(f)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=None,
                    help="Sweep directory (default: runs/sensitivity/<ts>_oat)")
    args = ap.parse_args()

    base = load_baseline()
    lo, hi = sa.RANGE

    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    sweep = Path(args.out_dir) if args.out_dir else sa.runs_dir() / f"{stamp}_oat"
    cfg_dir = sweep / "configs"
    run_root = sweep / "runs"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    run_root.mkdir(parents=True, exist_ok=True)

    plan = []

    def emit(run_id: str, param: str, factor, value):
        cfg = copy.deepcopy(base)
        if param is not None:
            cfg[param] = value
        out = run_root / run_id
        cfg["output_dir"] = str(out)
        cfg_path = cfg_dir / f"{run_id}.json"
        with open(cfg_path, "w") as f:
            json.dump(cfg, f, indent=2)
        plan.append({
            "run_id": run_id,
            "param": param or "",
            "paper_symbol": sa.PARAMS.get(param, "") if param else "",
            "factor": "" if factor is None else factor,
            "value": "" if value is None else value,
            "config": str(cfg_path),
            "output_dir": str(out),
        })

    # run 0 = unperturbed baseline
    emit("oat_0000", None, None, None)

    idx = 1
    for param in sa.PARAM_NAMES:
        if param not in base:
            raise KeyError(f"Param '{param}' not present in baseline config "
                           f"{sa.BASELINE_CONFIG}")
        p0 = float(base[param])
        for factor in (lo, hi):
            emit(f"oat_{idx:04d}", param, factor, p0 * factor)
            idx += 1

    plan_path = sweep / "oat_plan.csv"
    with open(plan_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["run_id", "param", "paper_symbol",
                                          "factor", "value", "config",
                                          "output_dir"])
        w.writeheader()
        w.writerows(plan)

    print(f"Sweep dir : {sweep}")
    print(f"Configs   : {len(plan)} ({len(sa.PARAM_NAMES)} params x 2 + 1 baseline)")
    print(f"Plan      : {plan_path}")
    print(f"Array size: --array=0-{len(plan) - 1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
