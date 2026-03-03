import csv
from pathlib import Path
from collections import defaultdict
import math

OUT_DIR = Path("data-export/sens_oat")

# Choose the endpoint exactly once here:
# Options available from your summary csv: C,P,E,N,H,R,total
ENDPOINT = "C"   # tumor population
DAYS = [7, 10, 16]

DELTA = 0.05  # matches run script

def load_summary(run_id: int, seed: int):
    path = OUT_DIR / f"summary_run_{run_id}_seed_{seed}.csv"
    data = {}
    with path.open() as f:
        r = csv.DictReader(f)
        for row in r:
            day = int(float(row["day"]))
            if day in DAYS:
                data[day] = float(row[ENDPOINT])
    # ensure all days exist
    for d in DAYS:
        if d not in data:
            raise RuntimeError(f"Missing day {d} in {path}")
    return data

def load_param_value(run_id: int, seed: int, key: str) -> float:
    path = OUT_DIR / f"params_run_{run_id}_seed_{seed}.csv"
    with path.open() as f:
        r = csv.DictReader(f)
        for row in r:
            if row["key"] == key:
                return float(row["value"])
    raise RuntimeError(f"Could not find {key} in {path}")

def mean(xs):
    return sum(xs)/len(xs) if xs else float("nan")

def main():
    run_plan = OUT_DIR / "run_plan.csv"
    rows = []
    with run_plan.open() as f:
        r = csv.DictReader(f)
        rows = list(r)

    baseline = [row for row in rows if row["scenario"] == "baseline"]
    perturb = [row for row in rows if row["scenario"] == "perturb"]

    # Group baseline replicates
    baseline_summaries = []
    for row in baseline:
        rid = int(row["run_id"]); seed = int(row["seed"])
        baseline_summaries.append(load_summary(rid, seed))

    # Baseline mean Y0(day)
    Y0 = {d: mean([s[d] for s in baseline_summaries]) for d in DAYS}

    # Use first baseline run to read p0
    base_rid = int(baseline[0]["run_id"])
    base_seed = int(baseline[0]["seed"])

    # Group perturb runs by (param, factor)
    grp = defaultdict(list)
    for row in perturb:
        grp[(row["param"], float(row["factor"]))].append(row)

    # Compute sensitivity per (param, day) using central difference from ±5%
    # S ≈ ( (Y_plus - Y_minus) / (2*DELTA*Y0) )  because dp/p0 = DELTA
    # More explicitly:
    #   dY/Y0 over dp/p0 -> (Yplus - Yminus) / (2*DELTA*Y0)
    results = []  # (param, day, S)

    params = sorted({row["param"] for row in perturb})

    for p in params:
        p0 = load_param_value(base_rid, base_seed, p)

        minus_runs = grp[(p, 1.0 - DELTA)]
        plus_runs  = grp[(p, 1.0 + DELTA)]

        # Mean Y_minus(day), Y_plus(day) across replicates
        Yminus = {d: mean([load_summary(int(r["run_id"]), int(r["seed"]))[d] for r in minus_runs]) for d in DAYS}
        Yplus  = {d: mean([load_summary(int(r["run_id"]), int(r["seed"]))[d] for r in plus_runs])  for d in DAYS}

        for d in DAYS:
            if Y0[d] == 0:
                S = float("nan")
            else:
                S = (Yplus[d] - Yminus[d]) / (2.0 * DELTA * Y0[d])
            results.append((p, d, S))

    # Write results table
    out_path = OUT_DIR / f"sensitivity_{ENDPOINT}.csv"
    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["param", "day", f"S_{ENDPOINT}"])
        for p, d, s in results:
            w.writerow([p, d, s])

    print(f"Wrote: {out_path}")

    # Print a quick ranking using |S| at the last day
    last_day = DAYS[-1]
    last = [(p, abs(s)) for (p, d, s) in results if d == last_day and not math.isnan(s)]
    last.sort(key=lambda x: x[1], reverse=True)

    print(f"\nTop sensitivities by |S| at day {last_day} for endpoint {ENDPOINT}:")
    for p, val in last[:10]:
        print(f"  {p:15s} |S|={val:.4f}")

if __name__ == "__main__":
    main()