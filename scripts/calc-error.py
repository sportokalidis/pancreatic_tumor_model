#!/usr/bin/env python3
"""
Compare simulation output (populations.csv) with experimental curves
from *_scaled_global.csv files that have **no headers**.

Outputs: fit_metrics_summary.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional




# --------- config ---------
BASE_DIR = Path("./data-export")  # change if needed
SIM_PATH = BASE_DIR / "populations.csv"

# Map: population -> (experimental CSV path, tokens to find the sim column)
POP_MAP = {
    "Tumor":   {"exp_csv": BASE_DIR / "C-Cells_scaled_global.csv", "tokens": ["tumor", "c", "c_cells", "c-cells"]},
    "PSC":     {"exp_csv": BASE_DIR / "P-Cells_scaled_global.csv", "tokens": ["psc", "p", "p_cells", "p-cells"]},
    "CD8":     {"exp_csv": BASE_DIR / "E-Cells_scaled_global.csv", "tokens": ["cd8", "cd8+", "e", "e_cells", "e-cells"]},
    "NK":      {"exp_csv": BASE_DIR / "N-Cells_scaled_global.csv", "tokens": ["nk", "n", "n_cells", "n-cells"]},
    "HelperT": {"exp_csv": BASE_DIR / "H-Cells_scaled_global.csv", "tokens": ["helper", "helpert", "th", "h", "h_cells", "h-cells"]},
    "Treg":    {"exp_csv": BASE_DIR / "R-Cells_scaled_global.csv", "tokens": ["treg", "tregs", "r", "r_cells", "r-cells"]},
}
# --------------------------

def read_exp_no_header(path: Path) -> pd.DataFrame:
    """
    Read experimental CSV with NO HEADER and force columns to [day, cell_population].
    Also ensures numeric types and sorts by day; drops duplicate days keeping the first.
    """
    df = pd.read_csv(path, header=None, names=["day", "cell_population"])
    df["day"] = pd.to_numeric(df["day"], errors="coerce")
    df["cell_population"] = pd.to_numeric(df["cell_population"], errors="coerce")
    df = df.dropna(subset=["day", "cell_population"])
    df = df.sort_values("day")
    # If any duplicate days exist, keep the first occurrence to avoid ambiguity
    df = df[~df["day"].duplicated(keep="first")].reset_index(drop=True)
    return df

def find_day_column(df: pd.DataFrame) -> str:
    for cand in df.columns:
        if str(cand).strip().lower() in ("day", "days", "time", "t"):
            return cand
    return df.columns[0]  # fallback

def find_sim_column(df: pd.DataFrame, tokens) -> Optional[str]:
    cols = list(df.columns)
    cols_lc = [str(c).lower() for c in cols]
    # exact match
    for i, c in enumerate(cols_lc):
        if c in tokens:
            return cols[i]
    # substring match
    for i, c in enumerate(cols_lc):
        if any(tok in c for tok in tokens):
            return cols[i]
    return None

# ---- metrics (no external deps) ----
def mae(y_true, y_pred):
    y_true = np.asarray(y_true, float); y_pred = np.asarray(y_pred, float)
    return float(np.nanmean(np.abs(y_pred - y_true)))

def rmse(y_true, y_pred):
    y_true = np.asarray(y_true, float); y_pred = np.asarray(y_pred, float)
    return float(np.sqrt(np.nanmean((y_pred - y_true)**2)))

def nrmse(y_true, y_pred, mode="range"):
    y_true = np.asarray(y_true, float); y_pred = np.asarray(y_pred, float)
    r = rmse(y_true, y_pred)
    if mode == "range":
        denom = np.nanmax(y_true) - np.nanmin(y_true)
    elif mode == "mean":
        denom = np.nanmean(y_true)
    else:
        raise ValueError("mode must be 'range' or 'mean'")
    return float(r / denom) if denom and np.isfinite(denom) else np.nan

def mape(y_true, y_pred):
    y_true = np.asarray(y_true, float); y_pred = np.asarray(y_pred, float)
    denom = np.where(y_true == 0, np.nan, y_true)
    pct = np.abs(y_pred - y_true) / denom * 100.0
    return float(np.nanmean(pct))

def r2(y_true, y_pred):
    y_true = np.asarray(y_true, float); y_pred = np.asarray(y_pred, float)
    ss_res = np.nansum((y_true - y_pred)**2)
    ss_tot = np.nansum((y_true - np.nanmean(y_true))**2)
    if ss_tot == 0 or not np.isfinite(ss_tot):
        return np.nan
    return float(1.0 - ss_res/ss_tot)
# ------------------------------------

def compute_fit_metrics(exp_df, sim_df, sim_col, sim_day_col="day"):
    """
    Interpolates experimental curve to simulation days and computes metrics.
    Discards sim days outside experimental domain.
    """
    exp_days = exp_df["day"].astype(float).values
    exp_vals = exp_df["cell_population"].astype(float).values
    sim_days = sim_df[sim_day_col].astype(float).values
    sim_vals = sim_df[sim_col].astype(float).values

    # Clean sim NaNs
    msk = np.isfinite(sim_days) & np.isfinite(sim_vals)
    sim_days = sim_days[msk]; sim_vals = sim_vals[msk]

    # Interpolate, mark outside-range as NaN then drop
    exp_interp = np.interp(sim_days, exp_days, exp_vals, left=np.nan, right=np.nan)
    valid = np.isfinite(exp_interp) & np.isfinite(sim_vals)
    y_true = exp_interp[valid]; y_pred = sim_vals[valid]

    if len(y_true) == 0:
        return {"n_points_compared": 0, "MAE": np.nan, "RMSE": np.nan,
                "NRMSE_range": np.nan, "NRMSE_mean": np.nan, "MAPE_%": np.nan, "R2": np.nan}

    return {
        "n_points_compared": int(len(y_true)),
        "MAE": mae(y_true, y_pred),
        "RMSE": rmse(y_true, y_pred),
        "NRMSE_range": nrmse(y_true, y_pred, mode="range"),
        "NRMSE_mean": nrmse(y_true, y_pred, mode="mean"),
        "MAPE_%": mape(y_true, y_pred),
        "R2": r2(y_true, y_pred),
    }

def compute_all_metrics(sim_csv: Path, refs_dir: Optional[Path] = None, ode_csv: Optional[Path] = None, scale_s: float = 1e5) -> pd.DataFrame:
    """
    Core logic: compare ABM (sim_csv) against either:
      1. ODE reference (ode_csv) — PREFERRED for validation (same scale/params)
      2. Paper reference directory (refs_dir) — scaled to ABM's scale_S

    Args:
        sim_csv: Path to ABM populations CSV
        refs_dir: Directory with paper reference CSVs (S=1e5)
        ode_csv: Path to ODE reference CSV (same scale/params as ABM)
        scale_s: Scale factor of ABM (e.g., 1e5, 1e4, 1e3). Used to scale paper refs.

    Returns a DataFrame with fit metrics.
    """
    sim_df = pd.read_csv(sim_csv)
    sim_day_col = find_day_column(sim_df)

    # Auto-detect scale_S from ABM if not provided
    if scale_s == 1e5:
        try:
            import json
            params_file = Path(sim_csv).parent / "params.json"
            if params_file.exists():
                with open(params_file) as f:
                    params = json.load(f)
                    scale_s = params.get("scale_S", 1e5)
        except:
            pass

    # If ODE reference provided, use it (proper validation)
    if ode_csv and ode_csv.exists():
        ode_df = pd.read_csv(ode_csv)
        ode_day_col = find_day_column(ode_df)

        pops = ["Tumor", "PSC", "CD8", "NK", "HelperT", "Treg"]
        tokens_map = {
            "Tumor":   ["tumor", "c", "c_cells", "c-cells"],
            "PSC":     ["psc", "p", "p_cells", "p-cells"],
            "CD8":     ["cd8", "cd8+", "e", "e_cells", "e-cells"],
            "NK":      ["nk", "n", "n_cells", "n-cells"],
            "HelperT": ["helper", "helpert", "th", "h", "h_cells", "h-cells"],
            "Treg":    ["treg", "tregs", "r", "r_cells", "r-cells"],
        }

        rows = []
        for pop in pops:
            ode_col = find_sim_column(ode_df, tokens_map[pop])
            sim_col = find_sim_column(sim_df, tokens_map[pop])

            if sim_col is None or ode_col is None:
                rows.append({
                    "population": pop, "sim_column": sim_col, "ode_column": ode_col,
                    "status": "missing", "n_points_compared": 0, "MAE": np.nan,
                    "RMSE": np.nan, "NRMSE_range": np.nan, "NRMSE_mean": np.nan,
                    "MAPE_%": np.nan, "R2": np.nan
                })
                continue

            # Create synthetic "paper ref" by treating ODE as ground truth
            exp_df = ode_df[[ode_day_col, ode_col]].copy()
            exp_df.columns = ["day", "cell_population"]

            metrics = compute_fit_metrics(exp_df, sim_df, sim_col, sim_day_col)
            rows.append({
                "population": pop, "sim_column": sim_col, "ode_column": ode_col,
                "status": "ok", **metrics
            })
        return pd.DataFrame(rows)

    # Fallback: use paper reference directory (auto-scaled to ABM's scale_S)
    if not refs_dir:
        raise ValueError("Must provide either ode_csv or refs_dir")

    # Paper references are at S=1e5; scale them to match ABM scale
    scale_factor = scale_s / 1e5

    pop_map = {
        "Tumor":   {"exp_csv": refs_dir / "C-Cells_scaled_global.csv", "tokens": ["tumor", "c", "c_cells", "c-cells"]},
        "PSC":     {"exp_csv": refs_dir / "P-Cells_scaled_global.csv", "tokens": ["psc", "p", "p_cells", "p-cells"]},
        "CD8":     {"exp_csv": refs_dir / "E-Cells_scaled_global.csv", "tokens": ["cd8", "cd8+", "e", "e_cells", "e-cells"]},
        "NK":      {"exp_csv": refs_dir / "N-Cells_scaled_global.csv", "tokens": ["nk", "n", "n_cells", "n-cells"]},
        "HelperT": {"exp_csv": refs_dir / "H-Cells_scaled_global.csv", "tokens": ["helper", "helpert", "th", "h", "h_cells", "h-cells"]},
        "Treg":    {"exp_csv": refs_dir / "R-Cells_scaled_global.csv", "tokens": ["treg", "tregs", "r", "r_cells", "r-cells"]},
    }

    rows = []
    for pop, meta in pop_map.items():
        exp_df = read_exp_no_header(meta["exp_csv"])
        # Scale the paper reference to match ABM's scale_S
        exp_df["cell_population"] = exp_df["cell_population"] * scale_factor

        sim_col = find_sim_column(sim_df, [t.lower() for t in meta["tokens"]])
        if sim_col is None:
            rows.append({
                "population": pop, "sim_column": None,
                "status": "missing in simulation",
                "n_points_compared": 0, "MAE": np.nan, "RMSE": np.nan,
                "NRMSE_range": np.nan, "NRMSE_mean": np.nan, "MAPE_%": np.nan, "R2": np.nan
            })
            continue
        metrics = compute_fit_metrics(exp_df, sim_df, sim_col, sim_day_col)
        rows.append({"population": pop, "sim_column": sim_col, "status": "ok", **metrics})

    return pd.DataFrame(rows)


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description="Compute fit metrics: ABM vs ODE reference (preferred) or paper reference (auto-scaled)")
    parser.add_argument("--sim",  default=str(BASE_DIR / "populations.csv"),
                        help="Path to simulation populations CSV")
    parser.add_argument("--ode-csv", default=None,
                        help="Path to ODE reference CSV (PREFERRED: same scale/params as ABM)")
    parser.add_argument("--refs", default=str(BASE_DIR),
                        help="Directory with *_scaled_global.csv paper references (at S=1e5)")
    parser.add_argument("--out",  default=str(BASE_DIR),
                        help="Output directory for fit_metrics_summary.csv")
    parser.add_argument("--scale-s", type=float, default=None,
                        help="Scale factor (auto-detected from params.json if not provided)")
    args = parser.parse_args()

    # Auto-detect scale_S from params.json if not provided
    scale_s = args.scale_s or 1e5
    params_file = Path(args.sim).parent / "params.json"
    if params_file.exists():
        try:
            with open(params_file) as f:
                params = json.load(f)
                scale_s = params.get("scale_S", scale_s)
        except:
            pass

    ode_path = Path(args.ode_csv) if args.ode_csv else None
    summary = compute_all_metrics(Path(args.sim), refs_dir=Path(args.refs), ode_csv=ode_path, scale_s=scale_s)

    out_path = Path(args.out) / "fit_metrics_summary.csv"
    summary.to_csv(out_path, index=False)
    print(summary.to_string(index=False))
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
