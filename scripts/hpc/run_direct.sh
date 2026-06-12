#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_direct.sh — run the ABM WITHOUT SLURM (Apptainer / interactive node).
#
# Two modes:
#   --base (default)   Single run of the base (no-treatment) model.
#                      Uses params.json or params_S<scale>.json.
#   --treatment        Run all 4 treatment protocols (treatment branch only).
#
# Usage:
#   # Base model, default scale S=1e5
#   bash scripts/hpc/run_direct.sh
#
#   # Base model, 10× more agents (S=1e4)
#   bash scripts/hpc/run_direct.sh --scale S1e4
#
#   # Base model, 100× more agents (S=1e3) — needs more RAM/time
#   bash scripts/hpc/run_direct.sh --scale S1e3
#
#   # Treatment protocols (4 runs, sequential)
#   bash scripts/hpc/run_direct.sh --treatment
#
#   # Treatment protocols, parallel background jobs
#   bash scripts/hpc/run_direct.sh --treatment --parallel
#
# Options:
#   --scale LABEL      S1e5 (default) | S1e4 | S1e3
#   --treatment        Run all 4 treatment protocols instead of base
#   --protocols LIST   With --treatment: comma-separated subset
#                      (default: acd47,gem,abr,abr_acd47)
#   --threads N        OMP_NUM_THREADS (default: auto)
#   --parallel         Run protocols in background simultaneously
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
MODE="base"
PROTOCOLS="acd47,gem,abr,abr_acd47"
THREADS=""
PARALLEL=false
SKIP_BUILD=false
NOTE="HPC direct run"
SEED_OVERRIDE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scale)      SCALE="$2";          shift 2 ;;
    --treatment)  MODE="treatment";    shift ;;
    --protocols)  PROTOCOLS="$2";      shift 2 ;;
    --threads)    THREADS="$2";        shift 2 ;;
    --parallel)   PARALLEL=true;       shift ;;
    --skip-build) SKIP_BUILD=true;     shift ;;
    --note)       NOTE="$2";           shift 2 ;;
    --seed)       SEED_OVERRIDE="$2";  shift 2 ;;
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
  echo "        Remove --skip-build to rebuild." >&2; exit 1
fi

# ---- thread count -----------------------------------------------------------
NCPU=$(nproc 2>/dev/null || echo 4)
if [ -z "${THREADS}" ]; then
  if [ "${MODE}" = "base" ]; then
    THREADS="${NCPU}"
  else
    IFS=',' read -ra _TMP <<< "${PROTOCOLS}"
    N_PROTO=${#_TMP[@]}
    THREADS=$(( NCPU / N_PROTO ))
    [ "${THREADS}" -lt 1 ] && THREADS=1
  fi
fi
export OMP_NUM_THREADS="${THREADS}"

SEED=$(${PYTHON} -c "import json; print(json.load(open('${REPO_ROOT}/params.json')).get('seed',42))" 2>/dev/null || echo 42)
[ -n "${SEED_OVERRIDE}" ] && SEED="${SEED_OVERRIDE}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# ============================================================
# BASE MODEL
# ============================================================
if [ "${MODE}" = "base" ]; then
  # Select config file
  if [ "${SCALE}" = "S1e5" ]; then
    CONFIG="${REPO_ROOT}/params.json"
  else
    CONFIG="${REPO_ROOT}/params_${SCALE}.json"
  fi

  # Auto-generate rescaled params if missing
  if [ ! -f "${CONFIG}" ]; then
    echo "[INFO] ${CONFIG} not found — generating via rescale_params.py..."
    S_NUM="${SCALE#S}"   # e.g. S1e4 → 1e4
    ${PYTHON} "${REPO_ROOT}/scripts/rescale_params.py" \
      "${REPO_ROOT}/params.json" "${S_NUM}" "${CONFIG}"
  fi

  OUTPUT_DIR="${REPO_ROOT}/output/base_${SCALE}"
  mkdir -p "${OUTPUT_DIR}"

  echo "[2/2] Running base model  (scale=${SCALE}  threads=${THREADS})"
  echo "  Config: ${CONFIG}"

  # The binary reads "params.json" from CWD via LoadParams(); BDM_PARAMS overrides
  # that path, and output_dir is patched to an absolute path so it is CWD-independent.
  TMP_CFG=$(mktemp /tmp/bdm_cfg_XXXXXX.json)
  ${PYTHON} -c "
import json
cfg = json.load(open('${CONFIG}'))
cfg['output_dir'] = '${OUTPUT_DIR}'
if '${SEED_OVERRIDE}':
    cfg['seed'] = int('${SEED_OVERRIDE}')
json.dump(cfg, open('${TMP_CFG}', 'w'), indent=2)
"
  START=$(date +%s)
  BDM_PARAMS="${TMP_CFG}" "${BINARY}"
  END=$(date +%s)
  rm -f "${TMP_CFG}"
  echo "  Done in $(( END - START ))s  → ${OUTPUT_DIR}/populations.csv"

  # Archive
  GROUP_DIR="${REPO_ROOT}/runs/${TIMESTAMP}_${SCALE}_s${SEED}"
  mkdir -p "${GROUP_DIR}"
  ${PYTHON} "${REPO_ROOT}/scripts/save_run.py" \
    --params   "${CONFIG}" \
    --abm      "${OUTPUT_DIR}/populations.csv" \
    --refs     "${REPO_ROOT}/data-export" \
    --runs-dir "${GROUP_DIR}" \
    --duration "$(( END - START ))" \
    --note     "${NOTE}"

  echo ""
  echo "================================================================"
  echo "  Base run complete.  Scale: ${SCALE}  Seed: ${SEED}"
  echo "  Output: ${OUTPUT_DIR}/populations.csv"
  echo "  Archive: ${GROUP_DIR}/"
  echo "================================================================"
  exit 0
fi

# ============================================================
# TREATMENT MODE  (requires treatment branch)
# ============================================================
IFS=',' read -ra PROTO_LIST <<< "${PROTOCOLS}"
N_PROTO=${#PROTO_LIST[@]}
echo "[2/2] Running ${N_PROTO} treatment protocols  (scale=${SCALE}  threads=${THREADS}  parallel=${PARALLEL})"

GROUP_DIR="${REPO_ROOT}/runs_treatment/${TIMESTAMP}_${SCALE}_s${SEED}"
mkdir -p "${GROUP_DIR}"
mkdir -p "${REPO_ROOT}/logs"
echo "  Group folder: ${GROUP_DIR}"

PIDS=()
for PROTO in "${PROTO_LIST[@]}"; do
  if [ "${SCALE}" = "S1e5" ]; then
    CONFIG="${REPO_ROOT}/params_treat_${PROTO}.json"
  else
    CONFIG="${REPO_ROOT}/params_treat_${PROTO}_${SCALE}.json"
  fi
  [ ! -f "${CONFIG}" ] && { echo "[ERROR] Config not found: ${CONFIG}" >&2; exit 1; }

  OUTPUT_DIR="${REPO_ROOT}/output/${PROTO}"
  mkdir -p "${OUTPUT_DIR}"
  echo "  → ${PROTO}  config: $(basename "${CONFIG}")"

  _run_proto() {
    local proto="$1" cfg="$2" outdir="$3"
    local start end tmp_cfg
    tmp_cfg=$(mktemp /tmp/bdm_cfg_XXXXXX.json)
    ${PYTHON} -c "
import json
c = json.load(open('$cfg'))
c['output_dir'] = '$outdir'
json.dump(c, open('$tmp_cfg', 'w'), indent=2)
"
    start=$(date +%s)
    BDM_PARAMS="${tmp_cfg}" "${BINARY}" \
      > "${REPO_ROOT}/logs/${proto}.log" 2>&1
    end=$(date +%s)
    rm -f "${tmp_cfg}"
    echo "  [${proto}] done in $(( end - start ))s"
    "${PYTHON}" "${REPO_ROOT}/scripts/save_run.py" \
      --params    "${cfg}" \
      --abm       "${outdir}/populations.csv" \
      --refs      "${REPO_ROOT}/data-export" \
      --group-dir "${GROUP_DIR}" \
      --duration  "$(( end - start ))" \
      --note      "${NOTE}"
  }

  if [ "${PARALLEL}" = true ]; then
    ( _run_proto "${PROTO}" "${CONFIG}" "${OUTPUT_DIR}" ) &
    PIDS+=($!)
  else
    _TMP=$(mktemp /tmp/bdm_cfg_XXXXXX.json)
    ${PYTHON} -c "
import json
c = json.load(open('${CONFIG}'))
c['output_dir'] = '${OUTPUT_DIR}'
json.dump(c, open('${_TMP}', 'w'), indent=2)
"
    START=$(date +%s)
    BDM_PARAMS="${_TMP}" "${BINARY}"
    END=$(date +%s)
    rm -f "${_TMP}"
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

if [ "${PARALLEL}" = true ] && [ ${#PIDS[@]} -gt 0 ]; then
  echo "  Waiting for ${#PIDS[@]} background jobs..."
  for pid in "${PIDS[@]}"; do wait "${pid}" || echo "[WARN] job ${pid} exited non-zero"; done
fi

echo ""
echo "================================================================"
echo "  All protocols complete.  Results: ${GROUP_DIR}"
echo "================================================================"
