#!/usr/bin/env python3
"""
Archive a completed ABM run into a structured folder and update the run index.

Baseline run creates:
  runs/YYYYMMDD_HHMMSS_s{seed}/
    params.json              ← exact parameters used
    populations.csv          ← ABM output
    ode_reference.csv        ← ODE with same params
    fit_metrics_summary.csv  ← R²/MAPE vs paper reference
    comparison_grid.png      ← 2×3 publication figure
    comparison_combined.png  ← all-in-one overview
    run_info.json            ← metadata (git, duration, note, metrics)

Treatment run (auto-detected when treat_* flags are True) creates:
  runs/YYYYMMDD_HHMMSS_{protocol}_s{seed}/
    params.json
    populations.csv
    ode_reference.csv        ← ODE in treatment mode
    treatment_comparison.png ← C + E vs ODE + window shading
    run_info.json

Updates:
  runs/index.csv  ← one row per run; all params + key metrics

Usage:
  python3 scripts/save_run.py [--note "description"] [--duration 120]
                              [--abm output/populations.csv]
                              [--params params.json]
                              [--refs data-export]
                              [--runs-dir runs]
                              [--label gem]   # optional name override
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT  = Path(__file__).resolve().parent.parent
PYTHON     = sys.executable
REFS_DIR   = REPO_ROOT / "data-export"
ODE_SCRIPT = REPO_ROOT / "scripts" / "ode_reference.py"
ERR_SCRIPT = REPO_ROOT / "scripts" / "calc-error.py"
PLT_SCRIPT = REPO_ROOT / "scripts" / "create-plots.py"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def git_info() -> dict:
    def run(cmd):
        try:
            return subprocess.check_output(cmd, cwd=REPO_ROOT,
                                           stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return "unknown"
    return {
        "commit": run(["git", "rev-parse", "--short", "HEAD"]),
        "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
    }


def params_hash(params: dict) -> str:
    numeric = {k: v for k, v in params.items() if isinstance(v, (int, float, bool))}
    blob = json.dumps(numeric, sort_keys=True).encode()
    return hashlib.md5(blob).hexdigest()[:8]


def load_params(path: Path) -> dict:
    with open(path) as f:
        raw = json.load(f)
    return {k: v for k, v in raw.items() if not k.startswith("_")}


def final_counts(populations_csv: Path) -> dict:
    df = pd.read_csv(populations_csv)
    last = df.iloc[-1]
    return {col: int(last[col]) for col in ["C", "P", "E", "N", "H", "R"]
            if col in df.columns}


def summarise_metrics(metrics_csv: Path) -> dict:
    df = pd.read_csv(metrics_csv)
    flat = {}
    for _, row in df.iterrows():
        pop = row["population"]
        flat[f"{pop}_R2"]   = round(float(row["R2"]),   4) if pd.notna(row["R2"])    else None
        flat[f"{pop}_MAPE"] = round(float(row["MAPE_%"]), 2) if pd.notna(row["MAPE_%"]) else None
    return flat


def compute_ode_metrics(abm_csv: Path, ode_csv: Path) -> dict:
    """Compute ABM vs ODE R² for each population (treatment mode: no paper reference)."""
    df_abm = pd.read_csv(abm_csv)
    df_ode = pd.read_csv(ode_csv)
    df_ode.columns = df_ode.columns.str.lower()

    metrics = {}
    for col in ["C", "P", "E", "N", "H", "R"]:
        col_lc = col.lower()
        if col not in df_abm.columns or col_lc not in df_ode.columns:
            continue
        abm_days = df_abm["days"].values.astype(float)
        abm_vals = df_abm[col].values.astype(float)
        ode_interp = np.interp(abm_days,
                               df_ode["days"].values.astype(float),
                               df_ode[col_lc].values.astype(float))
        ss_res = np.nansum((abm_vals - ode_interp) ** 2)
        ss_tot = np.nansum((abm_vals - np.nanmean(abm_vals)) ** 2)
        r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
        metrics[f"{col}_R2_vs_ODE"] = round(r2, 4)
    return metrics


def protocol_label(params: dict) -> str:
    """Build a short protocol name from active treat_* flags."""
    parts = []
    if params.get("treat_gem"):   parts.append("gem")
    if params.get("treat_abr"):   parts.append("abr")
    if params.get("treat_acd47"): parts.append("acd47")
    return "+".join(parts) if parts else "untreated"


def run_subprocess(cmd: list, label: str) -> bool:
    print(f"  [{label}] {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [{label}] FAILED:\n{result.stderr[-2000:]}")
        return False
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Archive an ABM run")
    parser.add_argument("--note",     default="",
                        help="Short human-readable description of this run")
    parser.add_argument("--duration", type=float, default=None,
                        help="Wall-clock duration of the simulation in seconds")
    parser.add_argument("--abm",      default="output/populations.csv",
                        help="Path to ABM populations CSV (default: output/populations.csv)")
    parser.add_argument("--params",   default="params.json",
                        help="Path to params JSON (default: params.json)")
    parser.add_argument("--refs",     default="data-export",
                        help="Directory with *_scaled_global.csv paper references")
    parser.add_argument("--runs-dir", default="runs",
                        help="Root directory for archived runs (default: runs/)")
    parser.add_argument("--label",    default="",
                        help="Optional label inserted into the run folder name (e.g. 'gem')")
    args = parser.parse_args()

    abm_path    = (REPO_ROOT / args.abm).resolve()
    params_path = (REPO_ROOT / args.params).resolve()
    refs_dir    = (REPO_ROOT / args.refs).resolve()
    runs_dir    = (REPO_ROOT / args.runs_dir).resolve()

    for p, name in [(abm_path, "--abm"), (params_path, "--params")]:
        if not p.exists():
            print(f"[ERROR] {name} not found: {p}")
            return 1

    # -----------------------------------------------------------------------
    # 1. Determine mode and create run folder
    # -----------------------------------------------------------------------
    params      = load_params(params_path)
    seed        = params.get("seed", 0)
    is_treatment = any(params.get(f"treat_{d}", False) for d in ["gem", "abr", "acd47"])
    proto       = protocol_label(params)
    label       = args.label or (proto if is_treatment else "")

    ts     = datetime.now()
    run_id = ts.strftime("%Y%m%d_%H%M%S")
    if label:
        run_id += f"_{label}"
    run_id += f"_s{seed}"

    run_dir = runs_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n[save_run] Archiving → {run_dir}")
    if is_treatment:
        print(f"  Mode: TREATMENT  protocol={proto}")

    # -----------------------------------------------------------------------
    # 2. Copy core files
    # -----------------------------------------------------------------------
    shutil.copy2(params_path, run_dir / "params.json")
    shutil.copy2(abm_path,    run_dir / "populations.csv")
    if not is_treatment:
        shutil.copy2(abm_path, REFS_DIR / "populations.csv")

    # -----------------------------------------------------------------------
    # 3. ODE reference
    # -----------------------------------------------------------------------
    scale_s = params.get("scale_S", 1e5)
    print("\n[1/3] Running ODE reference...")
    ode_cmd = [
        PYTHON, str(ODE_SCRIPT),
        "--params", str(run_dir / "params.json"),
        "--scale",  str(scale_s),
        "--abm",    str(run_dir / "populations.csv"),
        "--out",    str(run_dir),
    ]
    if is_treatment:
        ode_cmd += ["--mode", "treatment"]
    ok = run_subprocess(ode_cmd, "ODE")
    if not ok:
        print("  [WARN] ODE reference failed — skipping")

    # -----------------------------------------------------------------------
    # 4. Fit metrics
    # -----------------------------------------------------------------------
    ode_csv     = run_dir / "ode_reference.csv"
    metrics_flat = {}

    if is_treatment:
        print("\n[2/3] Computing ABM vs ODE metrics...")
        if ode_csv.exists():
            metrics_flat = compute_ode_metrics(run_dir / "populations.csv", ode_csv)
            print(f"  ABM vs ODE R²: " +
                  "  ".join(f"{k}={v}" for k, v in metrics_flat.items()))
        else:
            print("  [WARN] ODE csv missing — skipping metrics")
    else:
        print("\n[2/3] Computing fit metrics vs paper reference...")
        ok = run_subprocess([
            PYTHON, str(ERR_SCRIPT),
            "--sim",  str(run_dir / "populations.csv"),
            "--refs", str(refs_dir),
            "--out",  str(run_dir),
        ], "metrics")
        if not ok:
            print("  [WARN] Fit metrics failed — skipping")
        metrics_csv = run_dir / "fit_metrics_summary.csv"
        if metrics_csv.exists():
            metrics_flat = summarise_metrics(metrics_csv)

    # -----------------------------------------------------------------------
    # 5. Plots
    # -----------------------------------------------------------------------
    print("\n[3/3] Generating plots...")
    if is_treatment:
        ok = run_subprocess([
            PYTHON, str(PLT_SCRIPT),
            "--treatment",
            "--abm",    str(run_dir / "populations.csv"),
            "--ode",    str(ode_csv) if ode_csv.exists() else "",
            "--params", str(run_dir / "params.json"),
            "--out",    str(run_dir),
        ], "plots")
    else:
        ok = run_subprocess([
            PYTHON, str(PLT_SCRIPT),
            "--abm",  str(run_dir / "populations.csv"),
            "--ode",  str(ode_csv) if ode_csv.exists() else str(REFS_DIR / "ode_reference.csv"),
            "--refs", str(refs_dir),
            "--out",  str(run_dir),
        ], "plots")
    if not ok:
        print("  [WARN] Plot generation failed — skipping")

    # -----------------------------------------------------------------------
    # 6. Write run_info.json
    # -----------------------------------------------------------------------
    git    = git_info()
    counts = final_counts(run_dir / "populations.csv")

    run_info = {
        "run_id":         run_id,
        "timestamp":      ts.isoformat(timespec="seconds"),
        "note":           args.note,
        "git_commit":     git["commit"],
        "git_branch":     git["branch"],
        "duration_s":     args.duration,
        "params_hash":    params_hash(params),
        "seed":           seed,
        "total_days":     params.get("total_days"),
        "is_treatment":   is_treatment,
        "protocol":       proto,
        "final_counts":   counts,
        "fit_metrics":    metrics_flat,
    }
    with open(run_dir / "run_info.json", "w") as f:
        json.dump(run_info, f, indent=2)

    # -----------------------------------------------------------------------
    # 7. Update runs/index.csv
    # -----------------------------------------------------------------------
    index_path = runs_dir / "index.csv"

    row = {
        "run_id":        run_id,
        "timestamp":     ts.isoformat(timespec="seconds"),
        "note":          args.note,
        "git_commit":    git["commit"],
        "git_branch":    git["branch"],
        "duration_s":    args.duration,
        "params_hash":   params_hash(params),
        "is_treatment":  is_treatment,
        "protocol":      proto,
    }
    for k, v in params.items():
        if isinstance(v, (int, float, bool)):
            row[k] = v
    for pop, val in counts.items():
        row[f"{pop}_final"] = val
    row.update(metrics_flat)

    new_row_df = pd.DataFrame([row])
    if index_path.exists():
        index_df = pd.read_csv(index_path, dtype=str)
        for col in new_row_df.columns:
            if col not in index_df.columns:
                index_df[col] = np.nan
        for col in index_df.columns:
            if col not in new_row_df.columns:
                new_row_df[col] = np.nan
        index_df = pd.concat([index_df, new_row_df], ignore_index=True)
    else:
        index_df = new_row_df
    index_df.to_csv(index_path, index=False)

    # -----------------------------------------------------------------------
    # 8. Summary
    # -----------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"  Run archived: {run_id}")
    print(f"  Note: {args.note or '(none)'}")
    print(f"  Git:  {git['branch']}@{git['commit']}")
    if is_treatment:
        print(f"  Protocol: {proto}")
    print(f"  Final counts: C={counts.get('C','?')} P={counts.get('P','?')} "
          f"E={counts.get('E','?')} N={counts.get('N','?')} "
          f"H={counts.get('H','?')} R={counts.get('R','?')}")
    if metrics_flat:
        if is_treatment:
            print("  ABM vs ODE (R²):")
            for pop in ["C", "P", "E", "N", "H", "R"]:
                r2 = metrics_flat.get(f"{pop}_R2_vs_ODE", "?")
                print(f"    {pop}: R²={r2}")
        else:
            print("  Fit vs paper (R²):")
            for pop in ["Tumor", "PSC", "CD8", "NK", "HelperT", "Treg"]:
                r2  = metrics_flat.get(f"{pop}_R2",   "?")
                mpe = metrics_flat.get(f"{pop}_MAPE", "?")
                print(f"    {pop:8s}: R²={r2:>7}  MAPE={mpe:>8}%")
    print(f"  Index: {index_path}")
    print(f"{'='*60}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
