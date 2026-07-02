#!/usr/bin/env python3
"""
Analyze which E cell parameters have the most impact.

Test different mechanisms:
1. Carrying capacity (K_E) - prevents E from accumulating
2. Help from H reduction - less proliferation stimulus
3. More aggressive inactivation by C
4. Different combinations
"""

import json
from pathlib import Path

# Load baseline
with open("configs/params_S1e4_E10x.json") as f:
    base = json.load(f)

print("═" * 80)
print("E CELL MECHANISM ANALYSIS")
print("═" * 80)

print("\nCurrent v2 Parameters:")
print(f"  e_base_birth:      {base['e_base_birth']}")
print(f"  e_base_death:      {base['e_base_death']}")
print(f"  e_help_from_H:     {base['e_help_from_H']}")
print(f"  e_inact_by_C:      {base['e_inact_by_C']:.2e}")
print(f"  e_suppr_by_R:      {base['e_suppr_by_R']:.2e}")
print(f"  K_E:               {base['K_E']}")
print(f"  e_inact_by_C_K:    {base.get('e_inact_by_C_K', 'not set')}")
print(f"  e_suppr_by_R_K:    {base.get('e_suppr_by_R_K', 'not set')}")

print("\n" + "═" * 80)
print("RECOMMENDED PARAMETER COMBINATIONS TO TEST")
print("═" * 80)

test_cases = [
    {
        "name": "v4a: More aggressive C inactivation",
        "description": "e_inact_by_C ×4.0 (instead of ×2.5), reduces E suppression threshold",
        "adjustments": {
            "e_inact_by_C": 3.42e-06 * 4.0,
        }
    },
    {
        "name": "v4b: Reduce E help from H",
        "description": "e_help_from_H ×0.5 (less H-stimulated proliferation)",
        "adjustments": {
            "e_help_from_H": 0.125 * 0.5,
        }
    },
    {
        "name": "v4c: Lower E carrying capacity",
        "description": "K_E ÷2 (15,000 instead of 30,000) - limits E growth ceiling",
        "adjustments": {
            "K_E": 30000 / 2,
        }
    },
    {
        "name": "v4d: Combined aggressive E suppression",
        "description": "e_inact_by_C ×3.5, e_suppr_by_R ×5.0, K_E ÷2",
        "adjustments": {
            "e_inact_by_C": 3.42e-06 * 3.5,
            "e_suppr_by_R": 1.0e-06 * 5.0,
            "K_E": 30000 / 2,
        }
    },
    {
        "name": "v4e: Reduce help + increase inactivation",
        "description": "e_help_from_H ×0.6, e_inact_by_C ×3.5, e_suppr_by_R ×4.0",
        "adjustments": {
            "e_help_from_H": 0.125 * 0.6,
            "e_inact_by_C": 3.42e-06 * 3.5,
            "e_suppr_by_R": 1.0e-06 * 4.0,
        }
    },
    {
        "name": "v4f: Maximum reasonable adjustment",
        "description": "p_boost ×2.0, e_help_from_H ×0.5, e_inact_by_C ×5.0, e_suppr_by_R ×6.0",
        "adjustments": {
            "e_help_from_H": 0.125 * 0.5,
            "e_inact_by_C": 3.42e-06 * 5.0,
            "e_suppr_by_R": 1.0e-06 * 6.0,
        }
    },
]

print("\nRECOMMENDATIONS:")
print("-" * 80)

for i, test in enumerate(test_cases):
    print(f"\n{i+1}. {test['name']}")
    print(f"   {test['description']}")
    print(f"   Adjustments:")
    for key, val in test['adjustments'].items():
        base_val = base[key]
        if isinstance(val, float) and val < 1e-3:
            ratio = val / base_val
            print(f"     {key}: {base_val:.2e} → {val:.2e} (×{ratio:.1f})")
        else:
            ratio = val / base_val if base_val != 0 else 0
            print(f"     {key}: {base_val:.4f} → {val:.4f} (×{ratio:.2f})" if ratio != 0 else f"     {key}: {base_val} → {val}")

print("\n" + "═" * 80)
print("STRATEGY:")
print("═" * 80)
print("""
Test cases in order of aggressiveness:
  v4a: Surgical (only C inactivation) → least risky
  v4b: Reduces proliferation stimulus → moderate
  v4c: Carrying capacity ceiling → moderate
  v4d: Combined suppression → aggressive
  v4e: Multi-mechanism → balanced aggression
  v4f: Maximum tuning → highest risk of over-fitting

RECOMMENDATION: Test v4a and v4e first
  - v4a is minimally invasive (just boost C suppression)
  - v4e is balanced (targets proliferation + suppression)
""")
