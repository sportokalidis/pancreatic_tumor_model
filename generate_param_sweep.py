#!/usr/bin/env python3
"""
Parameter Sweep Generator for Pancreatic Tumor Model

This script generates multiple parameter files for running parameter sweeps
on HPC systems. It creates different scenarios by varying key parameters.

Usage:
    python3 generate_param_sweep.py [--sweep-type TYPE] [--num-jobs N] [--output-dir DIR]
"""

import os
import argparse
import itertools
import json
from pathlib import Path

# Base parameter set (default values)
BASE_PARAMS = {
    # Time & space
    "dt_minutes": 1.0,
    "min_bound": 0.0,
    "max_bound": 100.0,
    "cell_radius_um": 6.0,
    
    # Initial agent counts
    "C0": 474,  # Tumor
    "P0": 10,   # PSC
    "E0": 431,  # CD8+
    "N0": 189,  # NK
    "H0": 839,  # Helper T
    "R0": 65,   # Tregs
    
    # Global carrying capacities
    "K_C": 3000.0,
    "K_P": 3000.0,
    "K_E": 3000.0,
    "K_N": 3000.0,
    "K_H": 4000.0,
    "K_R": 3500.0,
    
    # Global gating
    "gate_C_K": 300.0,
    
    # Generic half-sats
    "K_small": 200.0,
    "K_med": 600.0,
    "K_big": 1200.0,
    
    # Tumor parameters
    "c_base_div": 0.05,
    "c_boost_from_P": 0.25,
    "c_boost_from_P_K": 800.0,
    "c_kill_by_E": 0.012,
    "c_kill_by_E_K": 600.0,
    "c_kill_by_N": 0.007,
    "c_kill_by_N_K": 600.0,
    "c_R_blocks_E": 0.10,
    
    # PSC parameters
    "p_base_div": 0.08,
    "p_boost_from_C": 0.10,
    "p_boost_from_C_K": 1000.0,
    "p_base_death": 0.05,
    
    # Effector T parameters
    "e_base_birth": 0.035,
    "e_help_from_H": 0.22,
    "e_help_from_H_K": 600.0,
    "e_inact_by_C": 0.15,
    "e_inact_by_C_K": 600.0,
    "e_suppr_by_R": 0.04,
    "e_suppr_by_R_K": 500.0,
    "e_base_death": 0.10,
    
    # NK parameters
    "n_base_birth": 0.03,
    "n_help_from_H": 0.15,
    "n_help_from_H_K": 600.0,
    "n_inact_by_C": 0.05,
    "n_inact_by_C_K": 600.0,
    "n_suppr_by_R": 0.038,
    "n_suppr_by_R_K": 500.0,
    "n_base_death": 0.10,
    
    # Helper T parameters
    "h_base_birth": 0.05,
    "h_self_act": 0.09,
    "h_self_act_K": 1000.0,
    "h_suppr_by_R": 0.075,
    "h_suppr_by_R_K": 700.0,
    "h_base_death": 0.08,
    
    # Treg parameters
    "r_base_src": 0.08,
    "r_induced_by_E": 0.008,
    "r_induced_by_E_K": 600.0,
    "r_induced_by_H": 0.008,
    "r_induced_by_H_K": 600.0,
    "r_cleared_by_N": 0.003,
    "r_cleared_by_N_K": 500.0,
    "r_decay": 0.06,
}

def write_param_file(params, filename):
    """Write parameters to a file in the expected format."""
    with open(filename, 'w') as f:
        f.write("# Pancreatic Tumor Model Parameters\n")
        f.write("# Auto-generated parameter file\n\n")
        
        # Group parameters by category
        categories = {
            "Time & Space": ["dt_minutes", "min_bound", "max_bound", "cell_radius_um"],
            "Initial Counts": ["C0", "P0", "E0", "N0", "H0", "R0"],
            "Carrying Capacities": ["K_C", "K_P", "K_E", "K_N", "K_H", "K_R", "gate_C_K"],
            "Generic Parameters": ["K_small", "K_med", "K_big"],
            "Tumor Parameters": [k for k in params.keys() if k.startswith("c_")],
            "PSC Parameters": [k for k in params.keys() if k.startswith("p_")],
            "Effector T Parameters": [k for k in params.keys() if k.startswith("e_")],
            "NK Parameters": [k for k in params.keys() if k.startswith("n_")],
            "Helper T Parameters": [k for k in params.keys() if k.startswith("h_")],
            "Treg Parameters": [k for k in params.keys() if k.startswith("r_")],
        }
        
        for category, param_list in categories.items():
            f.write(f"# {category}\n")
            for param in param_list:
                if param in params:
                    f.write(f"{param} = {params[param]}\n")
            f.write("\n")

def generate_tumor_growth_sweep(output_dir, num_jobs=10):
    """Generate parameter files varying tumor growth parameters."""
    print(f"Generating tumor growth parameter sweep with {num_jobs} jobs...")
    
    # Vary tumor division rate and PSC boost
    c_base_div_values = [0.03, 0.04, 0.05, 0.06, 0.07]
    c_boost_from_P_values = [0.15, 0.20, 0.25, 0.30, 0.35]
    
    combinations = list(itertools.product(c_base_div_values, c_boost_from_P_values))[:num_jobs]
    
    job_configs = []
    for i, (c_base_div, c_boost_from_P) in enumerate(combinations):
        params = BASE_PARAMS.copy()
        params["c_base_div"] = c_base_div
        params["c_boost_from_P"] = c_boost_from_P
        
        job_name = f"tumor_growth_{i:03d}"
        param_file = f"{output_dir}/params_{job_name}.txt"
        write_param_file(params, param_file)
        
        job_configs.append({
            "job_name": job_name,
            "param_file": param_file,
            "description": f"c_base_div={c_base_div}, c_boost_from_P={c_boost_from_P}"
        })
    
    return job_configs

def generate_immune_response_sweep(output_dir, num_jobs=10):
    """Generate parameter files varying immune response parameters."""
    print(f"Generating immune response parameter sweep with {num_jobs} jobs...")
    
    # Vary CD8+ killing efficiency and NK activity
    c_kill_by_E_values = [0.008, 0.010, 0.012, 0.014, 0.016]
    c_kill_by_N_values = [0.005, 0.006, 0.007, 0.008, 0.009]
    
    combinations = list(itertools.product(c_kill_by_E_values, c_kill_by_N_values))[:num_jobs]
    
    job_configs = []
    for i, (c_kill_by_E, c_kill_by_N) in enumerate(combinations):
        params = BASE_PARAMS.copy()
        params["c_kill_by_E"] = c_kill_by_E
        params["c_kill_by_N"] = c_kill_by_N
        
        job_name = f"immune_response_{i:03d}"
        param_file = f"{output_dir}/params_{job_name}.txt"
        write_param_file(params, param_file)
        
        job_configs.append({
            "job_name": job_name,
            "param_file": param_file,
            "description": f"c_kill_by_E={c_kill_by_E}, c_kill_by_N={c_kill_by_N}"
        })
    
    return job_configs

def generate_initial_conditions_sweep(output_dir, num_jobs=10):
    """Generate parameter files varying initial cell populations."""
    print(f"Generating initial conditions parameter sweep with {num_jobs} jobs...")
    
    # Vary initial tumor and immune cell counts
    C0_values = [300, 400, 474, 550, 650]
    E0_values = [300, 365, 431, 500, 570]
    
    combinations = list(itertools.product(C0_values, E0_values))[:num_jobs]
    
    job_configs = []
    for i, (C0, E0) in enumerate(combinations):
        params = BASE_PARAMS.copy()
        params["C0"] = C0
        params["E0"] = E0
        
        job_name = f"initial_conditions_{i:03d}"
        param_file = f"{output_dir}/params_{job_name}.txt"
        write_param_file(params, param_file)
        
        job_configs.append({
            "job_name": job_name,
            "param_file": param_file,
            "description": f"C0={C0}, E0={E0}"
        })
    
    return job_configs

def generate_treg_suppression_sweep(output_dir, num_jobs=10):
    """Generate parameter files varying Treg suppression parameters."""
    print(f"Generating Treg suppression parameter sweep with {num_jobs} jobs...")
    
    # Vary Treg suppression of effectors and helpers
    e_suppr_by_R_values = [0.02, 0.03, 0.04, 0.05, 0.06]
    h_suppr_by_R_values = [0.05, 0.06, 0.075, 0.09, 0.10]
    
    combinations = list(itertools.product(e_suppr_by_R_values, h_suppr_by_R_values))[:num_jobs]
    
    job_configs = []
    for i, (e_suppr_by_R, h_suppr_by_R) in enumerate(combinations):
        params = BASE_PARAMS.copy()
        params["e_suppr_by_R"] = e_suppr_by_R
        params["h_suppr_by_R"] = h_suppr_by_R
        
        job_name = f"treg_suppression_{i:03d}"
        param_file = f"{output_dir}/params_{job_name}.txt"
        write_param_file(params, param_file)
        
        job_configs.append({
            "job_name": job_name,
            "param_file": param_file,
            "description": f"e_suppr_by_R={e_suppr_by_R}, h_suppr_by_R={h_suppr_by_R}"
        })
    
    return job_configs

def main():
    parser = argparse.ArgumentParser(description="Generate parameter files for HPC parameter sweeps")
    parser.add_argument("--sweep-type", choices=["tumor_growth", "immune_response", "initial_conditions", "treg_suppression"], 
                       default="tumor_growth", help="Type of parameter sweep to generate")
    parser.add_argument("--num-jobs", type=int, default=10, help="Number of parameter combinations to generate")
    parser.add_argument("--output-dir", default="param_sweep", help="Directory to store parameter files")
    
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Generate parameter sweep based on type
    if args.sweep_type == "tumor_growth":
        job_configs = generate_tumor_growth_sweep(args.output_dir, args.num_jobs)
    elif args.sweep_type == "immune_response":
        job_configs = generate_immune_response_sweep(args.output_dir, args.num_jobs)
    elif args.sweep_type == "initial_conditions":
        job_configs = generate_initial_conditions_sweep(args.output_dir, args.num_jobs)
    elif args.sweep_type == "treg_suppression":
        job_configs = generate_treg_suppression_sweep(args.output_dir, args.num_jobs)
    
    # Save job configuration for use by submission scripts
    config_file = f"{args.output_dir}/job_configs.json"
    with open(config_file, 'w') as f:
        json.dump({
            "sweep_type": args.sweep_type,
            "num_jobs": len(job_configs),
            "jobs": job_configs
        }, f, indent=2)
    
    print(f"Generated {len(job_configs)} parameter files in {args.output_dir}/")
    print(f"Job configuration saved to {config_file}")
    print("\nExample job configurations:")
    for i, job in enumerate(job_configs[:3]):
        print(f"  {job['job_name']}: {job['description']}")
    if len(job_configs) > 3:
        print(f"  ... and {len(job_configs) - 3} more")

if __name__ == "__main__":
    main()