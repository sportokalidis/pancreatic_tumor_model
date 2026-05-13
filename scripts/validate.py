#!/usr/bin/env python3
"""
Automated validation: compare simulation output to scaled ODE reference.

Exit code 0 = all populations pass thresholds.
Exit code 1 = one or more populations fail (prints details).
Exit code 2 = input files missing or unreadable.

Usage:
    python3 scripts/validate.py [--populations-csv PATH] [--data-dir PATH] [--strict]

Thresholds (current model state — update as the model improves):
    Calibrated against a full 100-day simulation run (fit_metrics_summary.csv baseline):
      Tumor: R²=0.18, MAPE=42%  | PSC: R²=-0.45, MAPE=156%
      CD8:   R²=0.77, MAPE=563% | NK:  R²=0.80,  MAPE=21%
      Helper:R²=0.28, MAPE=37%  | Treg:R²=-1.70, MAPE=49%

    Thresholds are regression guards, not aspirational targets.
    Tighten them as each known issue (KI-01 through KI-05) is fixed.
    Known issues that explain failures: KI-01 (constant sources), KI-02 (overcrowding),
    KI-03 (K scale), KI-04 (P0 too small). See docs/known_issues.md.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Thresholds — update as model improves
# ---------------------------------------------------------------------------
# Format: population_name -> {metric: (min_ok, max_ok)}
# None means "don't check this metric for now"
# Baseline from 100-day run (fit_metrics_summary.csv); buffer ≈ 20% below baseline.
THRESHOLDS = {
    "Tumor": {
        "R2":     (0.0, None),    # baseline 0.18; alert if model gets worse than flat
        "MAPE_%": (None, 100.0),  # baseline 42%; double is a regression
    },
    "PSC": {
        "R2":     (-1.0, None),   # baseline -0.45; KI-04 pending
        "MAPE_%": (None, 300.0),  # baseline 156%
    },
    "CD8": {
        "R2":     (0.5, None),    # baseline 0.77; shape match is valuable
        "MAPE_%": (None, 800.0),  # baseline 563%; scale broken until KI-01 fixed
    },
    "NK": {
        "R2":     (0.5, None),    # baseline 0.80; alert if drops below 0.5
        "MAPE_%": (None, 50.0),   # baseline 21%
    },
    "HelperT": {
        "R2":     (0.0, None),    # baseline 0.28
        "MAPE_%": (None, 80.0),   # baseline 37%
    },
    "Treg": {
        "R2":     (-3.0, None),   # baseline -1.70; KI-01 pending
        "MAPE_%": (None, 100.0),  # baseline 49%
    },
}

# Minimum number of data points the comparison must cover
MIN_POINTS = 10

# Expected minimum simulation length (rows) for thresholds to be meaningful.
# Thresholds are calibrated against a 100-day run (~100 rows in populations.csv).
# A short run will give different metrics and may produce spurious failures.
EXPECTED_MIN_ROWS = 60

# ---------------------------------------------------------------------------
# Shared metric functions (duplicated from calc-error.py to avoid import)
# ---------------------------------------------------------------------------

def _mape(y_true, y_pred):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    denom = np.where(y_true == 0, np.nan, y_true)
    return float(np.nanmean(np.abs(y_pred - y_true) / denom * 100.0))


def _r2(y_true, y_pred):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    ss_res = np.nansum((y_true - y_pred) ** 2)
    ss_tot = np.nansum((y_true - np.nanmean(y_true)) ** 2)
    if ss_tot == 0 or not np.isfinite(ss_tot):
        return np.nan
    return float(1.0 - ss_res / ss_tot)


def _rmse(y_true, y_pred):
    return float(np.sqrt(np.nanmean((np.asarray(y_true, float) - np.asarray(y_pred, float)) ** 2)))


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

POP_MAP = {
    "Tumor":   ("C-Cells_scaled_global.csv", "C"),
    "PSC":     ("P-Cells_scaled_global.csv", "P"),
    "CD8":     ("E-Cells_scaled_global.csv", "E"),
    "NK":      ("N-Cells_scaled_global.csv", "N"),
    "HelperT": ("H-Cells_scaled_global.csv", "H"),
    "Treg":    ("R-Cells_scaled_global.csv", "R"),
}


def load_reference(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, header=None, names=["day", "value"])
    df = df.apply(pd.to_numeric, errors="coerce").dropna()
    df = df.sort_values("day").drop_duplicates("day", keep="first")
    return df


def load_simulation(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip().str.lower()
    day_col = next((c for c in df.columns if c in ("day", "days", "time", "t")), df.columns[0])
    df = df.rename(columns={day_col: "day"})
    return df


def compute_metrics(ref_df, sim_df, sim_col: str) -> dict:
    ref_days = ref_df["day"].values.astype(float)
    ref_vals = ref_df["value"].values.astype(float)
    sim_days = sim_df["day"].values.astype(float)
    sim_col_lc = sim_col.lower()
    if sim_col_lc not in sim_df.columns:
        return {"error": f"column '{sim_col}' not in simulation output"}

    sim_vals = sim_df[sim_col_lc].values.astype(float)
    msk = np.isfinite(sim_days) & np.isfinite(sim_vals)
    sim_days = sim_days[msk]
    sim_vals = sim_vals[msk]

    ref_interp = np.interp(sim_days, ref_days, ref_vals, left=np.nan, right=np.nan)
    valid = np.isfinite(ref_interp) & np.isfinite(sim_vals)
    y_true = ref_interp[valid]
    y_pred = sim_vals[valid]

    if len(y_true) < MIN_POINTS:
        return {"error": f"only {len(y_true)} comparable points (need ≥{MIN_POINTS})"}

    return {
        "n_points": int(len(y_true)),
        "R2": _r2(y_true, y_pred),
        "MAPE_%": _mape(y_true, y_pred),
        "RMSE": _rmse(y_true, y_pred),
    }


# ---------------------------------------------------------------------------
# Threshold checking
# ---------------------------------------------------------------------------

def check_threshold(metric_name: str, value: float, bounds: tuple) -> tuple[bool, str]:
    lo, hi = bounds
    if lo is not None and value < lo:
        return False, f"{metric_name}={value:.4g} < min {lo}"
    if hi is not None and value > hi:
        return False, f"{metric_name}={value:.4g} > max {hi}"
    return True, f"{metric_name}={value:.4g} OK"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Validate simulation against reference ODE data")
    parser.add_argument("--populations-csv", default="output/populations.csv",
                        help="Path to simulation populations.csv (default: output/populations.csv)")
    parser.add_argument("--data-dir", default="data-export",
                        help="Directory with reference CSVs (default: data-export)")
    parser.add_argument("--strict", action="store_true",
                        help="Treat warnings as failures")
    args = parser.parse_args()

    sim_path = Path(args.populations_csv)
    data_dir = Path(args.data_dir)

    if not sim_path.exists():
        print(f"[ERROR] Simulation output not found: {sim_path}", file=sys.stderr)
        sys.exit(2)

    try:
        sim_df = load_simulation(sim_path)
    except Exception as e:
        print(f"[ERROR] Cannot read {sim_path}: {e}", file=sys.stderr)
        sys.exit(2)

    failures = []
    warnings_list = []
    passed = []

    print(f"\nValidating: {sim_path.resolve()}")
    print(f"Reference:  {data_dir.resolve()}")

    n_rows = len(sim_df)
    if n_rows < EXPECTED_MIN_ROWS:
        msg = (f"Simulation output has only {n_rows} rows (expected ≥{EXPECTED_MIN_ROWS} "
               f"for a 100-day run). Thresholds are calibrated for a full run; "
               f"short-run metrics may differ significantly.")
        print(f"\n[WARN] {msg}")
        warnings_list.append(msg)
    print()
    print(f"{'Population':<12} {'n_pts':>6} {'R²':>8} {'MAPE%':>9}  Result")
    print("-" * 60)

    for pop_name, (ref_file, sim_col) in POP_MAP.items():
        ref_path = data_dir / ref_file
        if not ref_path.exists():
            msg = f"[SKIP] {pop_name}: reference file missing ({ref_path})"
            warnings_list.append(msg)
            print(f"  {msg}")
            continue

        ref_df = load_reference(ref_path)
        metrics = compute_metrics(ref_df, sim_df, sim_col)

        if "error" in metrics:
            msg = f"[SKIP] {pop_name}: {metrics['error']}"
            warnings_list.append(msg)
            print(f"  {msg}")
            continue

        pop_thresholds = THRESHOLDS.get(pop_name, {})
        pop_failures = []
        pop_ok = []

        for metric_name, bounds in pop_thresholds.items():
            if metric_name not in metrics:
                continue
            val = metrics[metric_name]
            if not np.isfinite(val):
                pop_failures.append(f"{metric_name}=NaN")
                continue
            ok, msg = check_threshold(metric_name, val, bounds)
            if ok:
                pop_ok.append(msg)
            else:
                pop_failures.append(msg)

        r2_val = metrics.get("R2", float("nan"))
        mape_val = metrics.get("MAPE_%", float("nan"))
        n_pts = metrics.get("n_points", 0)

        if pop_failures:
            result = "FAIL: " + "; ".join(pop_failures)
            failures.append(f"{pop_name}: {result}")
        else:
            result = "pass"
            passed.append(pop_name)

        print(f"  {pop_name:<12} {n_pts:>6} {r2_val:>8.3f} {mape_val:>9.1f}  {result}")

    print("-" * 60)
    print(f"\nPassed: {len(passed)}/{len(POP_MAP)}  |  Failed: {len(failures)}  |  Skipped: {len(warnings_list)}")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)

    if args.strict and warnings_list:
        print("\nWARNINGS (--strict mode, treating as failures):")
        for w in warnings_list:
            print(f"  - {w}")
        sys.exit(1)

    print("\nAll checks passed.")
    sys.exit(0)


if __name__ == "__main__":
    main()
