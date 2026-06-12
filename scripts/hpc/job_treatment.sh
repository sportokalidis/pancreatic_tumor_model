#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# SLURM array job — run all 4 treatment protocols in parallel.
#
# Each array task (0-3) handles one protocol:
#   0 → Anti-CD47
#   1 → Gemcitabine
#   2 → Abraxane
#   3 → Abraxane + Anti-CD47
#
# Results land in: runs_treatment/<TIMESTAMP>_<SCALE>_s<SEED>/
#
# Submit:
#   sbatch scripts/hpc/job_treatment.sh
#
# Override any #SBATCH option at submit time:
#   sbatch --partition=gpu --time=2:00:00 scripts/hpc/job_treatment.sh
#
# Environment variables (set before sbatch or in #SBATCH --export):
#   SCALE          Scale label, e.g. S1e5 (default) or S1e4
#   PARAMS_PREFIX  Config prefix  (default: params_treat_)
#   NOTE           Human-readable label for this run group
#   BDM_BUILD      Path to BioDynaMo build dir (auto-detected if unset)
#   PYTHON         Path to python3 (auto-detected if unset)
# ---------------------------------------------------------------------------

#SBATCH --job-name=tumor_treatment
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/slurm_treatment_%A_%a.out
#SBATCH --error=logs/slurm_treatment_%A_%a.err
# Uncomment and set your partition:
##SBATCH --partition=compute

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
mkdir -p "${REPO_ROOT}/logs"

# Load portable environment (BioDynaMo + Python)
source "${REPO_ROOT}/scripts/hpc/env.sh"

SCALE="${SCALE:-S1e5}"
PARAMS_PREFIX="${PARAMS_PREFIX:-params_treat_}"
NOTE="${NOTE:-HPC treatment run}"
BUILD_DIR="${REPO_ROOT}/build"
BINARY="${BUILD_DIR}/pancreatic_tumor_new"

# Protocol index → name mapping
PROTOCOLS=("acd47" "gem" "abr" "abr_acd47")
PROTO="${PROTOCOLS[${SLURM_ARRAY_TASK_ID:-0}]}"

# Config file (S1e5 omits the scale suffix to match existing filenames)
if [ "${SCALE}" = "S1e5" ]; then
  CONFIG="${REPO_ROOT}/${PARAMS_PREFIX}${PROTO}.json"
else
  CONFIG="${REPO_ROOT}/${PARAMS_PREFIX}${PROTO}_${SCALE}.json"
fi

if [ ! -f "${CONFIG}" ]; then
  echo "[ERROR] Config not found: ${CONFIG}" >&2; exit 1
fi
if [ ! -f "${BINARY}" ]; then
  echo "[ERROR] Binary not found: ${BINARY}" >&2
  echo "        Run 'bash scripts/hpc/build_hpc.sh' first." >&2; exit 1
fi

echo "================================================"
echo "  Protocol: ${PROTO}   Scale: ${SCALE}"
echo "  Config:   ${CONFIG}"
echo "  CPUs:     ${SLURM_CPUS_PER_TASK:-1}"
echo "  Node:     ${SLURMD_NODENAME:-local}"
echo "  $(date)"
echo "================================================"

# BioDynaMo reads OMP_NUM_THREADS for parallelism
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

OUTPUT_DIR="${REPO_ROOT}/output/hpc_${PROTO}"
mkdir -p "${OUTPUT_DIR}"

START=$(date +%s)
cd "${REPO_ROOT}"
"${BINARY}" --config "${CONFIG}" --output "${OUTPUT_DIR}"
END=$(date +%s)
echo "  Done in $(( END - START ))s  → ${OUTPUT_DIR}/populations.csv"

# Archive into runs_treatment group (uses the shared group timestamp written
# by task 0; other tasks detect it from the group_id file).
GROUP_ID_FILE="${REPO_ROOT}/output/.hpc_group_id"
if [ "${SLURM_ARRAY_TASK_ID:-0}" = "0" ]; then
  SEED=$(${PYTHON} -c "import json; print(json.load(open('${CONFIG}')).get('seed',42))")
  GID="$(date +%Y%m%d_%H%M%S)_${SCALE}_s${SEED}"
  echo "${GID}" > "${GROUP_ID_FILE}"
  echo "  Group ID: ${GID}"
else
  # Brief back-off: wait until task 0 writes the group id
  for i in $(seq 1 30); do
    [ -f "${GROUP_ID_FILE}" ] && break
    sleep 2
  done
  GID="$(cat "${GROUP_ID_FILE}" 2>/dev/null || echo "unknown")"
fi

GROUP_DIR="${REPO_ROOT}/runs_treatment/${GID}"
mkdir -p "${GROUP_DIR}"

"${PYTHON}" "${REPO_ROOT}/scripts/save_run.py" \
  --params    "${CONFIG}" \
  --abm       "${OUTPUT_DIR}/populations.csv" \
  --refs      "${REPO_ROOT}/data-export" \
  --group-dir "${GROUP_DIR}" \
  --duration  "$(( END - START ))" \
  --note      "${NOTE}"

echo "  Archived → ${GROUP_DIR}/${PROTO}/"
