#!/usr/bin/env python3
"""
Sensitivity analysis: 30 runs with S=1e4, dt=24h, 10x E cells.
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path


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


def submit_job(seed: int, scale: float, dt_hr: float, suite_dir: Path, repo_root: Path, custom_params_file: str):
    """Submit one ABM job with custom params."""
    s_str, dt_str = format_scale_dt(scale, dt_hr)
    job_name = f"bdm-{s_str}-{dt_str}-E10x-s{seed}"

    # Create run directory
    run_dir = suite_dir / f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{s_str}_{dt_str}_E10x_s{seed}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Create params file for this run from the custom params template
    base_params = json.load(open(custom_params_file))
    base_params["scale_S"] = scale
    base_params["dt_hr"] = dt_hr
    base_params["random_seed"] = seed

    param_file = run_dir / "params.json"
    with open(param_file, "w") as f:
        json.dump(base_params, f, indent=2)

    # Submit job
    if scale == 1e3:
        time_limit = "48:00:00"
    else:
        time_limit = "04:00:00"

    cmd = [
        "sbatch",
        "--partition=cpu",
        "--mem=30G",
        f"--time={time_limit}",
        "--job-name", job_name,
        "--export", f"ALL,REPO_ROOT={repo_root},SCALE={s_str}_{dt_str},SEED={seed},SKIP_BUILD=true,OUTPUT_DIR={run_dir}",
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
    parser = argparse.ArgumentParser(description="Sensitivity analysis: S=1e4, dt=24h, 10x E cells")
    parser.add_argument("--num-seeds", type=int, default=30, help="Number of seeds")
    parser.add_argument("--seed-start", type=int, default=1, help="First seed number")
    parser.add_argument("--output-base", default="runs/base", help="Output base directory")
    parser.add_argument("--skip-submit", action="store_true", help="Don't submit jobs")
    args = parser.parse_args()

    repo_root = Path.cwd()
    custom_params = "configs/params_S1e4_E10x.json"

    if not Path(custom_params).exists():
        print(f"ERROR: {custom_params} not found")
        return 1

    scale = 1e4
    dt_hr = 24.0
    s_str, dt_str = format_scale_dt(scale, dt_hr)

    # Create suite directory
    suite_name = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{s_str}_{dt_str}_E10x"
    suite_dir = Path(args.output_base) / suite_name
    suite_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Sensitivity Suite: {s_str} {dt_str} — 10x E cells")
    print(f"Seeds: {args.num_seeds} (starting from {args.seed_start})")
    print(f"Output: {suite_dir}")
    print(f"Params: {custom_params}")
    print(f"{'='*60}\n")

    if not args.skip_submit:
        print(f"Submitting {args.num_seeds} jobs...")
        jobs = {}
        for i in range(args.num_seeds):
            seed = args.seed_start + i
            job_id, run_dir = submit_job(seed, scale, dt_hr, suite_dir, repo_root, custom_params)
            if job_id:
                jobs[seed] = (job_id, run_dir)

        print(f"\n✓ Submitted {len(jobs)} jobs")

        # Save job tracking file
        job_file = suite_dir / "jobs.json"
        with open(job_file, "w") as f:
            json.dump({str(k): [v[0], str(v[1])] for k, v in jobs.items()}, f, indent=2)
        print(f"✓ Job list saved: {job_file}")
        print(f"\nMonitor with: squeue -u $USER | grep {s_str}-{dt_str}-E10x")

    print(f"\n{'='*60}\n")


if __name__ == "__main__":
    sys.exit(main())
