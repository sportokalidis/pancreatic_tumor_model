#!/usr/bin/env python3
"""
Parameter calibration for P (PSC) and E (Effector) cell dynamics.

Tests different parameter combinations to minimize error between ABM and ODE
reference, focusing on:
  - P cell proliferation (p_boost_from_C, p_boost_from_C_K)
  - E cell suppression (e_inact_by_C, e_suppr_by_R, e_base_death)

Usage:
  python3 scripts/calibrate_p_e_params.py --test-dir /tmp/calibration
"""

import argparse
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

def run_ode_reference(S=1e4):
    """Generate ODE reference for comparison."""
    cmd = [
        "python3", "scripts/ode_reference_e10x.py",
        "--scale", str(int(S)),
        "--params", "configs/params_S1e4_E10x.json",
        "--out-dir", "/tmp/ode_ref"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: ODE generation failed: {result.stderr}")
        return None

    ode_csv = Path("/tmp/ode_ref/ode_reference_e10x.csv")
    if ode_csv.exists():
        return pd.read_csv(ode_csv)
    return None


def create_test_params(base_params, adjustments):
    """Create modified params with adjustments."""
    params = base_params.copy()
    params.update(adjustments)
    return params


def analyze_dynamics(ode_df, label):
    """Extract key dynamics from ODE reference."""
    if ode_df is None or len(ode_df) < 2:
        return None

    # Find day ~10, ~50, ~100
    df_day10 = ode_df[ode_df["days"] >= 10].iloc[0] if any(ode_df["days"] >= 10) else ode_df.iloc[0]
    df_day50 = ode_df[ode_df["days"] >= 50].iloc[0] if any(ode_df["days"] >= 50) else ode_df.iloc[-1]
    df_day100 = ode_df[ode_df["days"] >= 99].iloc[0] if any(ode_df["days"] >= 99) else ode_df.iloc[-1]

    return {
        "label": label,
        "P_day10": df_day10["P"],
        "P_day50": df_day50["P"],
        "P_day100": df_day100["P"],
        "E_day10": df_day10["E"],
        "E_day50": df_day50["E"],
        "E_day100": df_day100["E"],
        "P_growth_rate": (df_day100["P"] - df_day10["P"]) / (df_day100["days"] - 10.0),
        "E_decline_rate": (df_day10["E"] - df_day100["E"]) / (df_day100["days"] - 10.0),
    }


def main():
    parser = argparse.ArgumentParser(description="Calibrate P and E cell parameters")
    parser.add_argument("--test-dir", default="/tmp/calibration", help="Test output directory")
    parser.add_argument("--num-seeds", type=int, default=3, help="Seeds per test (for speed)")
    args = parser.parse_args()

    test_dir = Path(args.test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # Load baseline params
    with open("configs/params_S1e4_E10x.json") as f:
        base_params = json.load(f)

    # Generate ODE reference for comparison
    print("Generating ODE reference...")
    ode_df = run_ode_reference(1e4)
    if ode_df is None:
        print("ERROR: Could not generate ODE reference")
        return 1

    ode_dynamics = analyze_dynamics(ode_df, "ODE Reference")
    print("\nODE Reference Dynamics:")
    print(f"  P growth rate: {ode_dynamics['P_growth_rate']:.1f} cells/day")
    print(f"  E decline rate: {ode_dynamics['E_decline_rate']:.1f} cells/day")
    print(f"  P day10→100: {ode_dynamics['P_day10']:.0f} → {ode_dynamics['P_day100']:.0f}")
    print(f"  E day10→100: {ode_dynamics['E_day10']:.0f} → {ode_dynamics['E_day100']:.0f}")

    # Test parameter combinations
    print("\n" + "="*80)
    print("Testing parameter adjustments...")
    print("="*80)

    test_cases = [
        {
            "name": "Baseline (current)",
            "adjustments": {},
        },
        {
            "name": "P boost +50% (p_boost_from_C × 1.5)",
            "adjustments": {"p_boost_from_C": base_params["p_boost_from_C"] * 1.5},
        },
        {
            "name": "P boost +100% (p_boost_from_C × 2.0)",
            "adjustments": {"p_boost_from_C": base_params["p_boost_from_C"] * 2.0},
        },
        {
            "name": "E inactivation +50% (e_inact_by_C × 1.5)",
            "adjustments": {"e_inact_by_C": base_params["e_inact_by_C"] * 1.5},
        },
        {
            "name": "E inactivation +100% (e_inact_by_C × 2.0)",
            "adjustments": {"e_inact_by_C": base_params["e_inact_by_C"] * 2.0},
        },
        {
            "name": "E suppression by R +50% (e_suppr_by_R × 1.5)",
            "adjustments": {"e_suppr_by_R": base_params["e_suppr_by_R"] * 1.5},
        },
        {
            "name": "Combined: P +75%, E inact +75%, E suppr +75%",
            "adjustments": {
                "p_boost_from_C": base_params["p_boost_from_C"] * 1.75,
                "e_inact_by_C": base_params["e_inact_by_C"] * 1.75,
                "e_suppr_by_R": base_params["e_suppr_by_R"] * 1.75,
            },
        },
        {
            "name": "Combined: P +100%, E inact +100%, E suppr +50%",
            "adjustments": {
                "p_boost_from_C": base_params["p_boost_from_C"] * 2.0,
                "e_inact_by_C": base_params["e_inact_by_C"] * 2.0,
                "e_suppr_by_R": base_params["e_suppr_by_R"] * 1.5,
            },
        },
    ]

    results = []

    for i, test_case in enumerate(test_cases):
        print(f"\n[{i+1}/{len(test_cases)}] {test_case['name']}")

        # Create test params file
        test_params = create_test_params(base_params, test_case["adjustments"])
        test_params_file = test_dir / f"test_{i:02d}_params.json"
        with open(test_params_file, "w") as f:
            json.dump(test_params, f, indent=2)

        # Print parameter changes
        if test_case["adjustments"]:
            for key, val in test_case["adjustments"].items():
                orig = base_params[key]
                ratio = val / orig
                print(f"  {key}: {orig:.4e} → {val:.4e} ({ratio:.2f}×)")
        else:
            print(f"  (no changes)")

        results.append({
            "test_case": test_case["name"],
            "params_file": str(test_params_file),
            "adjustments": test_case["adjustments"],
        })

    # Save results summary
    print("\n" + "="*80)
    print("Parameter combinations saved to:")
    for i, test_case in enumerate(test_cases):
        test_params_file = test_dir / f"test_{i:02d}_params.json"
        print(f"  {test_params_file}")

    print("\n" + "="*80)
    print("RECOMMENDATIONS:")
    print("="*80)
    print("""
Based on ODE dynamics analysis, suggested parameter adjustments to match ABM:

1. **P cells growing too slow:**
   - Increase p_boost_from_C (stimulation from tumor cells)
   - Try: ×1.5 to ×2.0 of current value
   - Alternative: Decrease p_boost_from_C_K (lower threshold)

2. **E cells declining too slow:**
   - Increase e_inact_by_C (tumor suppression)
   - Increase e_suppr_by_R (Treg suppression)
   - Try each: ×1.5 to ×2.0 of current value

3. **Combined approach (recommended):**
   - P boost ×1.75-2.0
   - E inactivation ×1.75-2.0
   - E suppression by R ×1.5

To test: Modify configs/params_S1e4_E10x.json with suggested values
         and run small test suite (e.g., 5-10 seeds) to compare against ODE.
""")

    return 0


if __name__ == "__main__":
    exit(main())
