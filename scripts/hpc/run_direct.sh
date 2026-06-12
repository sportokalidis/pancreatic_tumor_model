#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_direct.sh — run all 4 treatment protocols WITHOUT SLURM.
#
# Designed for Apptainer/Singularity containers and interactive HPC nodes
# where sbatch is not available.
#
# Usage:
#   bash scripts/hpc/run_direct.sh [OPTIONS]
#
# Options:
#   --scale LABEL      S1e5 (default) or S1e4
#   --protocols LIST   Comma-separated subset (default: acd47,gem,abr,abr_acd47)
#   --threads N        OMP threads per protocol (default: auto = nproc/4)
#   --parallel         Run all protocols simultaneously in the background
#   --skip-build       Skip cmake build (binary must already exist)
#   --note TEXT        Label stored in the run archive
#
# Environment variables (override auto-detection):
#   BDM_BUILD    Path to BioDynaMo build/install dir
#   PYTHON       Path to Python >= 3.7
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO_ROOT}/scripts/hpc/env.sh"

# ---- defaults ---------------------------------------------------------------
SCALE="S1e5"
PROTOCOLS="acd47,gem,abr,abr_acd47"
THREADS=""          # empty → computed below
PARALLEL=false
SKIP_BUILD=false
NOTE="HPC direct run"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scale)      SCALE="$2";     shift 2 ;;
    --protocols)  PROTOCOLS="$2"; shift 2 ;;
    --threads)    THREADS="$2";   shift 2 ;;
    --parallel)   PARALLEL=true;  shift ;;
    --skip-build) SKIP_BUILD=true; shift ;;
    --note)       NOTE="$2";      shift 2 ;;
    *) echo "[ERROR] Unknown option: $1" >&2; exit 1 ;;
  esac
done

BUILD_DIR="${REPO_ROOT}/build"
BINARY="${BUILD_DIR}/pancreatic_tumor_new"

# ---- build ------------------------------------------------------------------
if [ "${SKIP_BUILD}" = false ]; then
  echo "[1/2] Building..."
  cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBIODYNAMO_ROOT="${BDM_BUILD}" \
    -Wno-dev 2>&1 | tail -3
  cmake --build "${BUILD_DIR}" --parallel "$(nproc)"
  echo "      Build complete."
fi

if [ ! -f "${BINARY}" ]; then
  echo "[ERROR] Binary not found: ${BINARY}" >&2
  echo "        Run without --skip-build first." >&2
  exit 1
fi

# ---- thread count -----------------------------------------------------------
NCPU=$(nproc 2>/dev/null || echo 4)
IFS=',' read -ra PROTO_LIST <<< "${PROTOCOLS}"
N_PROTO=${#PROTO_LIST[@]}

if [ -z "${THREADS}" ]; then
  # Give each protocol a fair share; at least 1
  THREADS=$(( NCPU / N_PROTO ))
  [ "${THREADS}" -lt 1 ] && THREADS=1
fi
export OMP_NUM_THREADS="${THREADS}"
echo "[2/2] Running ${N_PROTO} protocols  (OMP_NUM_THREADS=${THREADS}  parallel=${PARALLEL})"

# ---- run protocols ----------------------------------------------------------
SEED=$(${PYTHON} -c "import json; print(json.load(open('${REPO_ROOT}/params.json')).get('seed',42))" 2>/dev/null || echo 42)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
GROUP_ID="${TIMESTAMP}_${SCALE}_s${SEED}"
GROUP_DIR="${REPO_ROOT}/runs_treatment/${GROUP_ID}"
mkdir -p "${GROUP_DIR}"
echo "  Group folder: ${GROUP_DIR}"
echo ""

PIDS=()
for PROTO in "${PROTO_LIST[@]}"; do
  if [ "${SCALE}" = "S1e5" ]; then
    CONFIG="${REPO_ROOT}/params_treat_${PROTO}.json"
  else
    CONFIG="${REPO_ROOT}/params_treat_${PROTO}_${SCALE}.json"
  fi

  if [ ! -f "${CONFIG}" ]; then
    echo "[ERROR] Config not found: ${CONFIG}" >&2; exit 1
  fi

  OUTPUT_DIR="${REPO_ROOT}/output/${PROTO}"
  mkdir -p "${OUTPUT_DIR}"

  echo "  → ${PROTO}  (config: $(basename "${CONFIG}"))"

  if [ "${PARALLEL}" = true ]; then
    # Background: all protocols run at the same time
    (
      START=$(date +%s)
      "${BINARY}" --config "${CONFIG}" --output "${OUTPUT_DIR}" \
        > "${REPO_ROOT}/logs/${PROTO}.log" 2>&1
      END=$(date +%s)
      echo "  [${PROTO}] done in $(( END - START ))s"
      "${PYTHON}" "${REPO_ROOT}/scripts/save_run.py" \
        --params    "${CONFIG}" \
        --abm       "${OUTPUT_DIR}/populations.csv" \
        --refs      "${REPO_ROOT}/data-export" \
        --group-dir "${GROUP_DIR}" \
        --duration  "$(( END - START ))" \
        --note      "${NOTE}"
    ) &
    PIDS+=($!)
  else
    # Sequential
    START=$(date +%s)
    "${BINARY}" --config "${CONFIG}" --output "${OUTPUT_DIR}"
    END=$(date +%s)
    echo "  [${PROTO}] done in $(( END - START ))s"
    "${PYTHON}" "${REPO_ROOT}/scripts/save_run.py" \
      --params    "${CONFIG}" \
      --abm       "${OUTPUT_DIR}/populations.csv" \
      --refs      "${REPO_ROOT}/data-export" \
      --group-dir "${GROUP_DIR}" \
      --duration  "$(( END - START ))" \
      --note      "${NOTE}"
  fi
done

# Wait for background jobs if running in parallel
if [ "${PARALLEL}" = true ] && [ ${#PIDS[@]} -gt 0 ]; then
  echo ""
  echo "  Waiting for ${#PIDS[@]} background jobs..."
  mkdir -p "${REPO_ROOT}/logs"
  for pid in "${PIDS[@]}"; do
    wait "${pid}" || echo "[WARN] job ${pid} exited non-zero"
  done
fi

echo ""
echo "============================================================"
echo "  All protocols complete.  Results: ${GROUP_DIR}"
echo "============================================================"
