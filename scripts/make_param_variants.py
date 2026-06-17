#!/usr/bin/env python3
"""
Generate treatment + CSC parameter files from the corrected base configs.

Every variant = a base config (params.json / params_S1e4.json) with a small
feature-flag overlay on top, so all variants automatically inherit the latest
corrected base body (Poisson spawner, dt, ICs, caps).

Usage:
    python3 scripts/make_param_variants.py          # regenerate all variants

Outputs (repo root):
    params_treat_acd47[_S1e4].json       Anti-CD47 alone        (Fig 5a)
    params_treat_gem[_S1e4].json         Gemcitabine alone      (Fig 5b)
    params_treat_abr[_S1e4].json         Abraxane alone         (Fig 5c)
    params_treat_abr_acd47[_S1e4].json   Abraxane + Anti-CD47   (Fig 5d)
    params_csc[_S1e4].json               CSC, untreated         (Sec 6)
    params_treat_abr_csc_S1e4.json       Abraxane + CSC relapse (Fig 6)

Feature flags default OFF in the base configs, so the base files themselves stay
the plain no-treatment / no-CSC model.
"""
import copy
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CONFIGS = REPO / "configs"   # all param JSONs live here

# Calibrated drug block (validated in the drug-treatment branch vs paper Fig 5).
# Rate constants (day^-1) — NOT scaled with S.
DRUG_BLOCK = {
    "treat_start_day": 14.0,
    "gem_gamma": 5.54, "gem_dose": 0.4, "gem_freq_days": 3.5,
    "gem_end_day": 56.0, "gem_c_c": 10.0, "gem_c_p": 7.5, "gem_c_immune": 1.8,
    "abr_gamma": 0.6161308, "abr_dose": 0.1, "abr_freq_days": 4.0,
    "abr_end_day": 28.0, "abr_c_c": 10.0, "abr_c_p": 10.0, "abr_c_immune": 6.4,
    "acd47_end_day": 35.0, "acd47_e_boost": 0.05,
}

# Protocol flag sets (Paper Section 5)
PROTOCOLS = {
    "acd47":     {"treat_gem": False, "treat_abr": False, "treat_acd47": True},
    "gem":       {"treat_gem": True,  "treat_abr": False, "treat_acd47": False},
    "abr":       {"treat_gem": False, "treat_abr": True,  "treat_acd47": False},
    "abr_acd47": {"treat_gem": False, "treat_abr": True,  "treat_acd47": True},
}

# CSC defaults (Paper Section 6, Table 3) — csc_sigma is scaled per base file.
CSC_BLOCK = {"csc_enable": True, "S0": 10, "csc_a1": 0.35, "csc_a2": 0.65,
             "csc_a3": 0.0, "csc_lambda_max": 0.07, "csc_lambda_min": 0.001,
             "csc_delta_s": 0.0}


def load(name: str) -> dict:
    with open(CONFIGS / name) as f:
        return json.load(f)


def write(name: str, data: dict) -> None:
    CONFIGS.mkdir(exist_ok=True)
    with open(CONFIGS / name, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  wrote configs/{name}")


def with_overlay(base: dict, *overlays: dict) -> dict:
    out = copy.deepcopy(base)
    for ov in overlays:
        out.update(ov)
    return out


def csc_block_for(base: dict) -> dict:
    """CSC block with csc_sigma = 1e7 / scale_S for this base."""
    S = float(base.get("scale_S", 1e5))
    blk = dict(CSC_BLOCK)
    blk["csc_sigma"] = 1e7 / S
    return blk


def main() -> int:
    bases = {"": "params.json", "_S1e4": "params_S1e4.json"}

    print("Generating treatment variants...")
    for suffix, base_name in bases.items():
        base = load(base_name)
        for proto, flags in PROTOCOLS.items():
            variant = with_overlay(base, DRUG_BLOCK, flags)
            write(f"params_treat_{proto}{suffix}.json", variant)

    print("Generating CSC variants...")
    for suffix, base_name in bases.items():
        base = load(base_name)
        write(f"params_csc{suffix}.json", with_overlay(base, csc_block_for(base)))

    # Fig 6: Abraxane + CSC relapse (S1e4 only — the working HPC scale)
    print("Generating Fig 6 (Abraxane + CSC) variant...")
    base = load("params_S1e4.json")
    fig6 = with_overlay(base, DRUG_BLOCK, PROTOCOLS["abr"], csc_block_for(base))
    write("params_treat_abr_csc_S1e4.json", fig6)

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
