#!/usr/bin/env python3
"""
Rescale an ABM params.json from one scale factor S to another.

Usage:
    python3 scripts/rescale_params.py params.json 1e4 params_S1e4.json
    python3 scripts/rescale_params.py params.json 1e3 params_S1e3.json

Scale rules (from sim_param.h comments):
    bilinear coupling (cell·day)^-1   :  ABM = paper × S   →  scale ×(S_new/S_old)
    source terms (cells/day)          :  ABM = paper ÷ S   →  scale ×(S_old/S_new)
    half-saturation K (cells)         :  ABM = paper ÷ S   →  scale ×(S_old/S_new)
    carrying caps K_ABM               :  = 1/(a·S)         →  scale ×(S_old/S_new)
    initial conditions                :  ABM = paper ÷ S   →  scale ×(S_old/S_new)
    rate constants (day^-1)           :  unchanged
"""

import argparse
import json
import math
from pathlib import Path

# Parameters whose ABM value = paper_value × S  (bilinear couplings)
# When S changes: new = old × (S_new / S_old)
SCALE_WITH_S = {
    "c_boost_from_P",   # mu_c × S
    "c_kill_by_E",      # d_c  × S
    "c_kill_by_N",      # b_c  × S
    "c_R_blocks_E",     # r1   × S   (Treg blocking denominator)
    "e_inact_by_C",     # c_e  × S
    "e_suppr_by_R",     # δ_e  × S
    "n_inact_by_C",     # c_n  × S
    "n_suppr_by_R",     # δ_n  × S
    "h_suppr_by_R",     # δ_h  × S
    "r_cleared_by_N",   # r    × S
}

# Parameters whose ABM value = paper_value ÷ S  (sources, half-sats, caps, ICs)
# When S changes: new = old × (S_old / S_new)
SCALE_INV_S = {
    "e_base_birth",      "n_base_birth",
    "h_base_birth",      "r_base_src",
    "p_boost_from_C_K",
    "e_help_from_H_K",   "n_help_from_H_K",
    "h_self_act_K",      "r_prolif_by_H_K",
    "K_C", "K_P", "K_E", "K_N", "K_H", "K_R",
    "C0",  "P0",  "E0",  "N0",  "H0",  "R0",
}

# Integer keys — rounded after rescaling
INTEGER_KEYS = {"C0", "P0", "E0", "N0", "H0", "R0"}


def rescale(src: Path, S_new: float, dst: Path) -> None:
    raw = json.load(open(src))
    S_old = float(raw.get("scale_S", 1e5))
    if abs(S_old - S_new) < 1:
        raise ValueError(f"Source already has scale_S={S_old:.0e}; nothing to do.")

    ratio_fwd = S_new / S_old   # for ×S params
    ratio_inv = S_old / S_new   # for ÷S params

    exp_new = int(round(math.log10(S_new)))

    out = {}
    for k, v in raw.items():
        if k.startswith("_"):
            # Update scale references in comment strings
            out[k] = v.replace(f"S=1e{int(round(math.log10(S_old)))}", f"S=1e{exp_new}") \
                       .replace(f"(S={S_old:.0e})", f"(S={S_new:.0e})")
        elif k == "scale_S":
            out[k] = int(S_new)
        elif k in SCALE_WITH_S:
            out[k] = v * ratio_fwd
        elif k in SCALE_INV_S:
            new_v = v * ratio_inv
            out[k] = round(new_v) if k in INTEGER_KEYS else new_v
        else:
            out[k] = v  # rate constants and non-physics keys: unchanged

    json.dump(out, open(dst, "w"), indent=2)

    print(f"Rescaled: {src}  S={S_old:.0e} → {dst}  S={S_new:.0e}")
    ics = {k: out[k] for k in ["C0", "P0", "E0", "N0", "H0", "R0"]}
    print(f"  ICs:  {ics}")
    print(f"  K_C={out['K_C']:.0f}  K_P={out['K_P']:.0f}")
    total_max = sum(out[k] for k in ["K_C", "K_P", "K_E", "K_N", "K_H"])
    print(f"  Max agents (sum of caps): ~{total_max:,.0f}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rescale ABM params.json to a different scale factor S.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("src",   help="source params.json (e.g. params.json)")
    parser.add_argument("S_new", type=float, help="target scale factor (e.g. 1e4)")
    parser.add_argument("dst",   help="output JSON path (e.g. params_S1e4.json)")
    args = parser.parse_args()
    rescale(Path(args.src), args.S_new, Path(args.dst))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
