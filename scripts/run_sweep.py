#!/usr/bin/env python3
"""
Parameter sweep driver for the pancreatic tumor ABM.

Usage:
    make sweep CONFIG=scripts/sweep_config.yaml
    # or directly (with BDM env sourced):
    python3 scripts/run_sweep.py --config scripts/sweep_config.yaml

Sweep config format (YAML):
    runs:
      - name: baseline
        params: {}
      - name: high_treg_reinforce
        params:
          r_reinforce_by_H: 0.08
          r_reinforce_by_H_K: 200.0

Each run:
  1. Patches the given real_t params in src/pancreatic_tumor_model.h via regex
  2. Runs `biodynamo run` (recompiles only changed files, then executes)
  3. Runs calc-error.py to produce fit_metrics_summary.csv
  4. Archives populations.csv + fit_metrics_summary.csv → results/<name>/
  5. Restores the original header (even on failure)

After all runs, prints a side-by-side R² comparison table.
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import yaml
    import pandas as pd
except ImportError:
    print("ERROR: pyyaml and pandas required. Source the BDM environment first:")
    print("  source ~/Documents/dev/MyBDM/biodynamo/build/bin/thisbdm.sh")
    sys.exit(1)

PROJECT_ROOT = Path(__file__).parent.parent
HEADER       = PROJECT_ROOT / "src" / "pancreatic_tumor_model.h"
DATA_DIR     = PROJECT_ROOT / "data-export"
RESULTS_DIR  = PROJECT_ROOT / "results"


def patch_header(content: str, params: dict) -> str:
    """Replace real_t <name> = <old>; with real_t <name> = <new>; for each param."""
    for name, value in params.items():
        pattern     = rf'(\breal_t\s+{re.escape(name)}\s*=\s*)[\d.e+\-]+'
        replacement = rf'\g<1>{value}'
        new_content, n = re.subn(pattern, replacement, content)
        if n == 0:
            raise ValueError(f"Parameter '{name}' not found in header — check spelling")
        content = new_content
    return content


def run_biodynamo() -> bool:
    result = subprocess.run(["biodynamo", "run"], cwd=PROJECT_ROOT)
    return result.returncode == 0


def run_metrics() -> bool:
    result = subprocess.run(
        [sys.executable, str(DATA_DIR / "calc-error.py")],
        cwd=PROJECT_ROOT,
    )
    return result.returncode == 0


def archive_run(name: str):
    dest = RESULTS_DIR / name
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy(DATA_DIR / "populations.csv", dest)
    metrics_src = DATA_DIR / "fit_metrics_summary.csv"
    if metrics_src.exists():
        shutil.copy(metrics_src, dest)
    print(f"  Archived → results/{name}/")


def load_metrics(name: str) -> "pd.DataFrame | None":
    path = RESULTS_DIR / name / "fit_metrics_summary.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)[["population", "R2", "MAPE_%"]]
    df = df.rename(columns={"R2": f"{name}_R2", "MAPE_%": f"{name}_MAPE"})
    return df


def main():
    parser = argparse.ArgumentParser(description="Run ABM parameter sweep")
    parser.add_argument("--config", required=True, help="Path to YAML sweep config")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        print(f"ERROR: config not found: {config_path}")
        sys.exit(1)

    config = yaml.safe_load(config_path.read_text())
    runs   = config.get("runs", [])
    if not runs:
        print("No runs defined in config.")
        sys.exit(0)

    original = HEADER.read_text()
    completed = []

    for run in runs:
        name   = run["name"]
        params = run.get("params", {})
        print(f"\n{'='*55}")
        print(f"Run: {name}  ({len(params)} override(s))")
        for k, v in params.items():
            print(f"  {k} = {v}")
        print()

        try:
            HEADER.write_text(patch_header(original, params))

            if not run_biodynamo():
                print(f"  FAILED: biodynamo run returned non-zero")
                continue

            run_metrics()
            archive_run(name)
            completed.append(name)

        except ValueError as e:
            print(f"  SKIPPED: {e}")
        finally:
            HEADER.write_text(original)  # always restore, even on crash

    # side-by-side comparison table
    print(f"\n{'='*55}")
    print("Sweep summary (R²):\n")
    if completed:
        from functools import reduce
        dfs = [load_metrics(n) for n in completed]
        dfs = [d for d in dfs if d is not None]
        if dfs:
            merged = reduce(lambda a, b: a.merge(b, on="population"), dfs)
            print(merged.to_string(index=False))
        else:
            print("  No metrics files found.")
    else:
        print("  No runs completed.")


if __name__ == "__main__":
    main()
