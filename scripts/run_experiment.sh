#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Full experiment pipeline: Build → Run → Archive
#
# Each invocation creates a timestamped folder under runs/ containing the
# params, output, ODE reference, fit metrics, and publication-quality plots.
# A one-row summary is appended to runs/index.csv.
#
# Usage:
#   ./scripts/run_experiment.sh [--note "description"] [--skip-build] [--skip-run]
#
# Options:
#   --note TEXT    Human-readable label for this run (recorded in index)
#   --skip-build   Reuse existing binary in ./build/
#   --skip-run     Skip simulation, archive existing output/populations.csv
#
# Exit codes: 0=success, 1=archive error, 2=build/run error
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
BINARY="${BUILD_DIR}/pancreatic_tumor_new"
OUTPUT_CSV="${REPO_ROOT}/output/populations.csv"

# Portable environment: resolves BDM_BUILD and PYTHON on any machine.
source "${REPO_ROOT}/scripts/hpc/env.sh" || exit 2

NOTE=""
SKIP_BUILD=false
SKIP_RUN=false

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --note)        NOTE="$2"; shift 2 ;;
    --skip-build)  SKIP_BUILD=true; shift ;;
    --skip-run)    SKIP_RUN=true; SKIP_BUILD=true; shift ;;
    *) echo "[ERROR] Unknown option: $1" >&2; exit 2 ;;
  esac
done

echo "============================================================"
echo "  Pancreatic Tumor ABM — Experiment Pipeline"
echo "  $(date)"
if [[ -n "${NOTE}" ]]; then
  echo "  Note: ${NOTE}"
fi
echo "============================================================"

# -- 1. Build ---------------------------------------------------------------
if [ "${SKIP_BUILD}" = false ]; then
  echo ""
  echo "[1/3] Building..."
  cd "${REPO_ROOT}"
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -Wno-dev -DBIODYNAMO_ROOT="${BIODYNAMO_ROOT:-}" 2>&1 | tail -3
  cmake --build build --parallel "$(nproc)" 2>&1 | tail -15
  echo "      Build complete."
else
  echo "[1/3] Build skipped."
fi

if [ ! -f "${BINARY}" ]; then
  echo "[ERROR] Binary not found: ${BINARY}" >&2
  echo "  Run without --skip-build first." >&2
  exit 2
fi

# -- 2. Simulate -------------------------------------------------------------
START_TIME=$(date +%s)

if [ "${SKIP_RUN}" = false ]; then
  echo ""
  echo "[2/3] Running simulation..."
  mkdir -p "${REPO_ROOT}/output"
  cd "${REPO_ROOT}"
  # Binary reads the config from BDM_PARAMS (param files now live in configs/).
  BDM_PARAMS="${REPO_ROOT}/configs/params.json" "${BINARY}"
  echo "      Simulation complete. Output: ${OUTPUT_CSV}"
else
  echo "[2/3] Simulation run skipped."
fi

END_TIME=$(date +%s)
DURATION=$(( END_TIME - START_TIME ))

if [ ! -f "${OUTPUT_CSV}" ]; then
  echo "[ERROR] No simulation output found: ${OUTPUT_CSV}" >&2
  exit 2
fi

# -- 3. Archive run ----------------------------------------------------------
echo ""
echo "[3/3] Archiving run..."
cd "${REPO_ROOT}"
"${PYTHON}" scripts/save_run.py \
  --note     "${NOTE}" \
  --duration "${DURATION}" \
  --abm      "output/populations.csv" \
  --params   "configs/params.json" \
  --refs     "data-export" \
  --runs-dir "runs"
