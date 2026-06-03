#!/usr/bin/env python3
"""
Rescale params.json from one ABM scale factor S to another.

Updates all biology parameters that depend on S (bilinear couplings,
source terms, half-saturations, initial conditions, carrying capacities).
Leaves spatial, treatment, seed, and dimension-less rate params untouched.

Usage:
  python3 scripts/rescale_params.py --s 1e4 [--in params.json] [--out params_S1e4.json]
"""

import argparse
import json
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Maps: params.json key → build_paper_params key (both must exist)
# ---------------------------------------------------------------------------
RESCALE_MAP = {
    # Initial conditions (÷S)
    "C0": "C0", "P0": "P0", "E0": "E0", "N0": "N0", "H0": "H0", "R0": "R0",
    # Carrying capacities (÷S via 1/(a·S))
    "K_C": "K_C",
    "K_P": "K_P",
    # Tumor bilinear couplings (×S)
    "c_boost_from_P": "mu_c",
    "c_kill_by_E":    "d_c",
    "c_kill_by_N":    "b_c",
    "c_R_blocks_E":   "r1",
    # PSC: half-saturation for C→P boost (÷S)
    "p_boost_from_C_K": "mu_p",
    # CTL (E)
    "e_base_birth":   "a_e",       # source ÷S
    "e_help_from_H_K":"g_e",       # half-sat ÷S
    "e_inact_by_C":   "c_e",       # bilinear ×S
    "e_suppr_by_R":   "delta_e",   # bilinear ×S
    # NK (N)
    "n_base_birth":   "a_n",
    "n_help_from_H_K":"g_n",
    "n_inact_by_C":   "c_n",
    "n_suppr_by_R":   "delta_n",
    # Helper T (H)
    "h_base_birth":   "a_h",
    "h_self_act_K":   "g_h",
    "h_suppr_by_R":   "delta_h",
    # Treg (R)
    "r_base_src":     "a",
    "r_prolif_by_H_K":"g_r",
    "r_cleared_by_N": "r",
}


def main():
    parser = argparse.ArgumentParser(description="Rescale params.json to a new ABM scale S")
    parser.add_argument("--s",   type=float, required=True,
                        help="New scale factor (e.g. 1e4 = 10,000 real cells per ABM agent)")
    parser.add_argument("--in",  dest="infile",  default="params.json",
                        help="Input params JSON  (default: params.json)")
    parser.add_argument("--out", dest="outfile", default=None,
                        help="Output params JSON (default: params_S{s}.json)")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo_root / "scripts"))
    from ode_reference import build_paper_params

    in_path  = (repo_root / args.infile).resolve()
    new_S    = args.s
    out_name = args.outfile or f"params_S{new_S:.0e}.json".replace("+", "")
    out_path = (repo_root / out_name).resolve()

    with open(in_path) as f:
        p = json.load(f)

    paper = build_paper_params(new_S)

    # Apply rescaling
    updated = []
    for json_key, paper_key in RESCALE_MAP.items():
        if json_key in p and paper_key in paper:
            p[json_key] = paper[paper_key]
            updated.append(json_key)

    p["scale_S"] = new_S

    with open(out_path, "w") as f:
        json.dump(p, f, indent=2)

    print(f"Rescaled {in_path.name} → {out_path.name}  (S={new_S:.0e})")
    print(f"Updated {len(updated)} params: {', '.join(updated)}")
    print(f"Initial conditions: C0={p['C0']}  P0={p['P0']}  "
          f"E0={p['E0']}  N0={p['N0']}  H0={p['H0']}  R0={p['R0']}")


if __name__ == "__main__":
    main()
