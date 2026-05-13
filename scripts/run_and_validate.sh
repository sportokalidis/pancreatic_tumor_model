#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Build → Run → Validate pipeline
#
# Usage:
#   ./scripts/run_and_validate.sh [--skip-build] [--skip-run] [--strict]
#
# Options:
#   --skip-build   Use existing build (must already exist in ./build/)
#   --skip-run     Skip simulation, only validate existing output/populations.csv
#   --strict       Pass --strict to validate.py (warnings become failures)
#
# Exit codes: 0=all pass, 1=validation failed, 2=build/run failed
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Source BioDynaMo environment (provides Python with numpy/pandas and the BDM libs)
BDM_ENV="${HOME}/Documents/dev/biodynamo/build/bin/thisbdm.sh"
if [ -f "${BDM_ENV}" ]; then
  source "${BDM_ENV}" 2>/dev/null || true
else
  echo "[WARN] BioDynaMo environment not found at ${BDM_ENV}" >&2
  echo "       Ensure numpy/pandas are available in the current Python environment." >&2
fi
BUILD_DIR="${REPO_ROOT}/build"
OUTPUT_CSV="${REPO_ROOT}/output/populations.csv"
VALIDATE_PY="${REPO_ROOT}/scripts/validate.py"

SKIP_BUILD=false
SKIP_RUN=false
STRICT_FLAG=""

for arg in "$@"; do
  case "$arg" in
    --skip-build) SKIP_BUILD=true ;;
    --skip-run)   SKIP_RUN=true; SKIP_BUILD=true ;;
    --strict)     STRICT_FLAG="--strict" ;;
    *) echo "[ERROR] Unknown option: $arg"; exit 2 ;;
  esac
done

echo "============================================================"
echo "  Pancreatic Tumor ABM — Build / Run / Validate"
echo "  $(date)"
echo "============================================================"

# -- 1. Build ---------------------------------------------------------------
if [ "$SKIP_BUILD" = false ]; then
  echo ""
  echo "[STEP 1/3] Building..."
  mkdir -p "${BUILD_DIR}"
  cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release -Wno-dev 2>&1 | tail -5
  cmake --build "${BUILD_DIR}" --parallel "$(nproc)" 2>&1 | tail -20
  echo "[BUILD] Done."
else
  echo "[STEP 1/3] Build skipped."
fi

# -- 2. Run simulation -------------------------------------------------------
if [ "$SKIP_RUN" = false ]; then
  echo ""
  echo "[STEP 2/3] Running simulation..."
  mkdir -p "${REPO_ROOT}/output"

  # Find the built binary (BioDynaMo puts it in build/pancreatic_tumor_model)
  BIN="${BUILD_DIR}/pancreatic_tumor_model"
  if [ ! -f "${BIN}" ]; then
    echo "[ERROR] Binary not found: ${BIN}" >&2
    echo "  Try building first or check CMakeLists.txt target name." >&2
    exit 2
  fi

  # Run from REPO_ROOT so params.json is found at ./params.json
  cd "${REPO_ROOT}"
  "${BIN}"
  echo "[RUN] Done. Output: ${OUTPUT_CSV}"
else
  echo "[STEP 2/3] Simulation run skipped."
fi

# -- 3. Validate -------------------------------------------------------------
echo ""
echo "[STEP 3/3] Validating output against reference..."

if [ ! -f "${OUTPUT_CSV}" ]; then
  echo "[ERROR] Simulation output not found: ${OUTPUT_CSV}" >&2
  echo "  Run the simulation first or use --skip-run only when output exists." >&2
  exit 2
fi

# Always compute error metrics first (updates fit_metrics_summary.csv)
cd "${REPO_ROOT}"
python3 scripts/calc-error.py 2>/dev/null || true

# Then run threshold validation
python3 "${VALIDATE_PY}" \
  --populations-csv "${OUTPUT_CSV}" \
  --data-dir "${REPO_ROOT}/data-export" \
  ${STRICT_FLAG}
