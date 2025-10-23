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

def main():
    # Load simulation with headers
    sim_df = pd.read_csv(SIM_PATH)
    sim_day_col = find_day_column(sim_df)

    rows = []
    for pop, meta in POP_MAP.items():
        exp_df = read_exp_no_header(meta["exp_csv"])
        sim_col = find_sim_column(sim_df, [t.lower() for t in meta["tokens"]])
        if sim_col is None:
            rows.append({
                "population": pop,
                "sim_column": None,
                "status": "missing in simulation",
                "n_points_compared": 0, "MAE": np.nan, "RMSE": np.nan,
                "NRMSE_range": np.nan, "NRMSE_mean": np.nan, "MAPE_%": np.nan, "R2": np.nan
            })
            continue

        metrics = compute_fit_metrics(exp_df, sim_df, sim_col, sim_day_col)
        rows.append({"population": pop, "sim_column": sim_col, "status": "ok", **metrics})

    summary = pd.DataFrame(rows)
    out = BASE_DIR / "fit_metrics_summary.csv"
    summary.to_csv(out, index=False)
    print(summary.to_string(index=False))
    print(f"\nSaved: {out}")

if __name__ == "__main__":
    main()
