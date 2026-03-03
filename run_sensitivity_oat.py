import os
import csv
import subprocess
from pathlib import Path

EXEC = "./build/pancreatic_tumor_new"          # adjust if needed
OUT_DIR = "data-export/sens_oat"
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

DELTA = 0.05
REPLICATES = 3              # set 1 if you want very fast
BASELINE_RUN_ID_START = 0
SEED_START = 100

# The 10 “paper-style” parameters to vary
PARAMS = [
    "c_base_div",
    # "c_boost_from_P",
    # "c_kill_by_E",
    # "c_kill_by_N",
    "p_boost_from_C",
    "e_base_birth",
    # "e_inact_by_C",
    # "e_suppr_by_R",
    "r_base_src",
    # "r_decay",
]

def run_one(run_id: int, seed: int, set_arg=None):
    cmd = [EXEC, f"--run_id={run_id}", f"--seed={seed}", f"--out_dir={OUT_DIR}"]
    if set_arg is not None:
        cmd.append(f"--set={set_arg}")
    subprocess.run(cmd, check=True)

def read_baseline_value_from_params_file(run_id: int, seed: int, key: str) -> float:
    # file written by your model:
    # params_run_<run_id>_seed_<seed>.csv with rows: run_id,seed,key,value
    p = Path(OUT_DIR) / f"params_run_{run_id}_seed_{seed}.csv"
    with p.open() as f:
        r = csv.DictReader(f)
        for row in r:
            if row["key"] == key:
                return float(row["value"])
    raise RuntimeError(f"Could not find key={key} in {p}")

def main():
    run_plan_path = Path(OUT_DIR) / "run_plan.csv"
    with run_plan_path.open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["run_id", "seed", "scenario", "param", "factor", "set_arg"])

        run_id = BASELINE_RUN_ID_START
        seed = SEED_START

        # ---- Baseline replicates ----
        baseline_runs = []
        for rep in range(REPLICATES):
            print(f"[baseline] rep={rep} run_id={run_id} seed={seed}")
            run_one(run_id, seed, None)
            w.writerow([run_id, seed, "baseline", "", 1.0, ""])
            baseline_runs.append((run_id, seed))
            run_id += 1
            seed += 1

        # Use the first baseline run to get baseline parameter values (deterministic params)
        baseline_ref_run_id, baseline_ref_seed = baseline_runs[0]

        # ---- Perturb each parameter ±5% with replicates ----
        for p in PARAMS:
            p0 = read_baseline_value_from_params_file(baseline_ref_run_id, baseline_ref_seed, p)
            for factor in (1.0 - DELTA, 1.0 + DELTA):
                p_new = p0 * factor
                set_arg = f"{p}={p_new}"
                for rep in range(REPLICATES):
                    print(f"[perturb] {p} factor={factor} rep={rep} run_id={run_id} seed={seed}")
                    run_one(run_id, seed, set_arg)
                    w.writerow([run_id, seed, "perturb", p, factor, set_arg])
                    run_id += 1
                    seed += 1

    print(f"Done. Run plan written to: {run_plan_path}")

if __name__ == "__main__":
    main()