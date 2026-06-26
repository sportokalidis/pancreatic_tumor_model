#!/usr/bin/env python3
"""
Configurable validation suite runner.

Runs multiple ABM simulations with different seeds, generates ODE references,
calculates fit metrics, and produces convergence validation plots.

Usage:
  python3 scripts/run_validation_suite.py \\
    --scale 1e4 \\
    --dt-hr 24 \\
    --num-seeds 30 \\
    --seed-start 1 \\
    --output-base runs/base \\
    --skip-submit \\
    --skip-post

Configuration:
  --scale: Scale factor (1e5, 1e4, 1e3, etc.)
  --dt-hr: Time step in hours (0.25=6h, 1=24h, etc.)
  --num-seeds: Number of seeds to run
  --seed-start: First seed number (default 1)
  --output-base: Base directory for results (default: runs/base)
  --skip-submit: Don't submit jobs, just prepare config
  --skip-post: Skip post-processing (metrics, plots)
  --wait-timeout: Max seconds to wait for jobs (default 36000 = 10h)
"""

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
import shutil


def parse_args():
    parser = argparse.ArgumentParser(description="Run validation suite")
    parser.add_argument("--scale", type=float, default=1e4, help="Scale factor")
    parser.add_argument("--dt-hr", type=float, default=24, help="Time step (hours)")
    parser.add_argument("--num-seeds", type=int, default=30, help="Number of seeds")
    parser.add_argument("--seed-start", type=int, default=1, help="First seed number")
    parser.add_argument("--output-base", default="runs/base", help="Output base directory")
    parser.add_argument("--skip-submit", action="store_true", help="Don't submit jobs")
    parser.add_argument("--skip-post", action="store_true", help="Skip post-processing")
    parser.add_argument("--wait-timeout", type=int, default=36000, help="Max wait time (sec)")
    return parser.parse_args()


def format_scale_dt(scale: float, dt_hr: float) -> tuple:
    """Format scale and dt into readable names."""
    if scale == 1e5:
        s_str = "S1e5"
    elif scale == 1e4:
        s_str = "S1e4"
    elif scale == 1e3:
        s_str = "S1e3"
    else:
        s_str = f"S{scale:.0e}".replace("+0", "").replace("e", "e")

    if dt_hr == 0.25:
        dt_str = "dt6h"
    elif dt_hr == 1.0:
        dt_str = "dt24h"
    else:
        dt_str = f"dt{dt_hr:.0f}h"

    return s_str, dt_str


def create_param_file(seed: int, scale: float, dt_hr: float, output_dir: Path) -> Path:
    """Create params.json for a run."""
    # Load base params
    base_params = json.load(open("configs/params.json"))

    # Update for this run
    base_params["scale_S"] = scale
    base_params["dt_hr"] = dt_hr
    base_params["random_seed"] = seed

    # Save to run directory
    param_file = output_dir / "params.json"
    with open(param_file, "w") as f:
        json.dump(base_params, f, indent=2)

    return param_file


def submit_job(seed: int, scale: float, dt_hr: float, suite_dir: Path, repo_root: Path):
    """Submit one ABM job."""
    s_str, dt_str = format_scale_dt(scale, dt_hr)
    job_name = f"bdm-{s_str}-{dt_str}-s{seed}"

    # Create run directory
    run_dir = suite_dir / f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{s_str}_{dt_str}_s{seed}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Create params file
    create_param_file(seed, scale, dt_hr, run_dir)

    # Submit job
    cmd = [
        "sbatch",
        "--partition=cpu",
        "--mem=30G",
        "--job-name", job_name,
        "--export", f"ALL,REPO_ROOT={repo_root},SCALE={s_str}_{dt_str},SEED={seed},SKIP_BUILD=true",
        "scripts/hpc/job_base.sh"
    ]

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if result.returncode != 0:
        print(f"ERROR submitting job {seed}: {result.stderr}")
        return None

    job_id = result.stdout.strip().split()[-1]
    print(f"  Seed {seed:3d}: Job {job_id} → {run_dir.name}")
    return job_id, run_dir


def main():
    args = parse_args()
    repo_root = Path.cwd()
    s_str, dt_str = format_scale_dt(args.scale, args.dt_hr)

    # Create suite directory
    suite_name = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{s_str}_{dt_str}"
    suite_dir = Path(args.output_base) / suite_name
    suite_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Validation Suite: {s_str} {dt_str}")
    print(f"Seeds: {args.num_seeds} (starting from {args.seed_start})")
    print(f"Output: {suite_dir}")
    print(f"{'='*60}\n")

    if not args.skip_submit:
        print(f"Submitting {args.num_seeds} jobs...")
        jobs = {}
        for i in range(args.num_seeds):
            seed = args.seed_start + i
            job_id, run_dir = submit_job(seed, args.scale, args.dt_hr, suite_dir, repo_root)
            if job_id:
                jobs[seed] = (job_id, run_dir)

        print(f"\n✓ Submitted {len(jobs)} jobs")

        # Save job tracking file
        job_file = suite_dir / "jobs.json"
        with open(job_file, "w") as f:
            json.dump({str(k): [v[0], str(v[1])] for k, v in jobs.items()}, f, indent=2)
        print(f"✓ Job list saved: {job_file}")
        print(f"\nMonitor with: squeue -u $USER | grep {s_str}-{dt_str}")

    if not args.skip_post:
        print(f"\n(Post-processing will be done when jobs complete)")
        print(f"Run: python3 scripts/run_validation_suite.py --scale {args.scale} --dt-hr {args.dt_hr} --skip-submit")

    print(f"\n{'='*60}\n")


if __name__ == "__main__":
    main()
