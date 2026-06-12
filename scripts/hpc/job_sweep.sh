#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# SLURM array job — seed sweep across all 4 treatment protocols.
#
# Runs N_SEEDS × 4 combinations.  Array task ID encodes:
#   task_id = seed_index * 4 + protocol_index
#
# Results land in: runs_treatment/<TIMESTAMP>_<SCALE>_sweep/
#
# Submit:
#   N_SEEDS=5 sbatch scripts/hpc/job_sweep.sh
#
# With 5 seeds: --array=0-19  (5×4=20 tasks)
# With 10 seeds: --array=0-39
# Edit the #SBATCH --array line below to match N_SEEDS*4-1.
#
# Environment variables:
#   N_SEEDS        Number of random seeds to sweep (default: 5)
#   SEEDS          Explicit comma-separated seed list  (overrides N_SEEDS)
#                  e.g. SEEDS=42,123,456,789,1000
#   SCALE          Scale label: S1e5 (default) or S1e4
#   BDM_BUILD      Path to BioDynaMo build dir (auto-detected if unset)
#   PYTHON         Path to python3 (auto-detected if unset)
# ---------------------------------------------------------------------------

#SBATCH --job-name=tumor_sweep
#SBATCH --array=0-19        # N_SEEDS*4 - 1  →  change this when N_SEEDS changes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/slurm_sweep_%A_%a.out
#SBATCH --error=logs/slurm_sweep_%A_%a.err
##SBATCH --partition=compute

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
mkdir -p "${REPO_ROOT}/logs"

source "${REPO_ROOT}/scripts/hpc/env.sh"

SCALE="${SCALE:-S1e5}"
N_SEEDS="${N_SEEDS:-5}"
BUILD_DIR="${REPO_ROOT}/build"
BINARY="${BUILD_DIR}/pancreatic_tumor_new"

# Build seeds list
if [ -n "${SEEDS:-}" ]; then
  IFS=',' read -ra SEED_LIST <<< "${SEEDS}"
else
  # Deterministic seeds spaced by 1000 starting at 42
  SEED_LIST=()
  for i in $(seq 0 $(( N_SEEDS - 1 ))); do
    SEED_LIST+=("$(( 42 + i * 1000 ))")
  done
fi

PROTOCOLS=("acd47" "gem" "abr" "abr_acd47")
N_PROTO=4

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
SEED_IDX=$(( TASK_ID / N_PROTO ))
PROTO_IDX=$(( TASK_ID % N_PROTO ))

SEED="${SEED_LIST[$SEED_IDX]}"
PROTO="${PROTOCOLS[$PROTO_IDX]}"

# Build a per-seed config by injecting the seed into the base config
if [ "${SCALE}" = "S1e5" ]; then
  BASE_CONFIG="${REPO_ROOT}/params_treat_${PROTO}.json"
else
  BASE_CONFIG="${REPO_ROOT}/params_treat_${PROTO}_${SCALE}.json"
fi

if [ ! -f "${BASE_CONFIG}" ]; then
  echo "[ERROR] Base config not found: ${BASE_CONFIG}" >&2; exit 1
fi

# Write a tmp config with the overridden seed
TMP_CONFIG="${REPO_ROOT}/output/.tmp_seed${SEED}_${PROTO}.json"
mkdir -p "${REPO_ROOT}/output"
${PYTHON} - <<PYEOF
import json, pathlib
p = json.load(open("${BASE_CONFIG}"))
p["seed"] = ${SEED}
json.dump(p, open("${TMP_CONFIG}", "w"), indent=2)
PYEOF

if [ ! -f "${BINARY}" ]; then
  echo "[ERROR] Binary not found: ${BINARY}" >&2
  echo "        Run 'bash scripts/hpc/build_hpc.sh' first." >&2; exit 1
fi

echo "================================================"
echo "  Sweep task ${TASK_ID}: proto=${PROTO}  seed=${SEED}"
echo "  Scale: ${SCALE}   CPUs: ${SLURM_CPUS_PER_TASK:-1}"
echo "  $(date)"
echo "================================================"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

OUTPUT_DIR="${REPO_ROOT}/output/sweep_${PROTO}_s${SEED}"
mkdir -p "${OUTPUT_DIR}"

START=$(date +%s)
cd "${REPO_ROOT}"
"${BINARY}" --config "${TMP_CONFIG}" --output "${OUTPUT_DIR}"
END=$(date +%s)
echo "  Done in $(( END - START ))s"

# Archive into a sweep group folder (keyed by parent job ID)
SWEEP_GROUP="${REPO_ROOT}/runs_treatment/sweep_${SLURM_ARRAY_JOB_ID:-local}_${SCALE}"
mkdir -p "${SWEEP_GROUP}"

"${PYTHON}" "${REPO_ROOT}/scripts/save_run.py" \
  --params    "${TMP_CONFIG}" \
  --abm       "${OUTPUT_DIR}/populations.csv" \
  --refs      "${REPO_ROOT}/data-export" \
  --group-dir "${SWEEP_GROUP}" \
  --duration  "$(( END - START ))" \
  --note      "sweep seed=${SEED} proto=${PROTO}"

rm -f "${TMP_CONFIG}"
echo "  Archived → ${SWEEP_GROUP}/${PROTO}_s${SEED}/"
