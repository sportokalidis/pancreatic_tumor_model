#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Run all 4 treatment protocols and archive them into ONE group folder.
#
# Creates:
#   runs_treatment/<TIMESTAMP>_<SCALE>_s<SEED>/
#     acd47/          ← Anti-CD47 run
#     gem/            ← Gemcitabine run
#     abr/            ← Abraxane run
#     abr+acd47/      ← Abraxane + Anti-CD47 run
#     group_summary.json
#
# Usage:
#   ./scripts/run_treatment_group.sh [OPTIONS]
#
# Options:
#   --scale LABEL      Scale label in the folder name (default: S1e5)
#   --params-prefix P  Prefix for config files (default: params_treat_)
#                      e.g. --params-prefix params_treat_ --scale S1e4
#                      will use params_treat_{proto}_S1e4.json
#   --note TEXT        Human-readable note for all runs
#   --skip-build       Reuse existing binary
#   --protocols LIST   Comma-separated subset (default: acd47,gem,abr,abr+acd47)
#
# Exit codes: 0=success, non-zero on first failure
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
BINARY="${BUILD_DIR}/pancreatic_tumor_new"
RUNS_TREATMENT="${REPO_ROOT}/runs_treatment"

# Portable environment: resolves BDM_BUILD and PYTHON on any machine.
source "${REPO_ROOT}/scripts/hpc/env.sh" || exit 1

SCALE="S1e5"
PARAMS_PREFIX="params_treat_"
NOTE=""
SKIP_BUILD=false
PROTOCOLS="acd47,gem,abr,abr+acd47"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scale)           SCALE="$2"; shift 2 ;;
    --params-prefix)   PARAMS_PREFIX="$2"; shift 2 ;;
    --note)            NOTE="$2"; shift 2 ;;
    --skip-build)      SKIP_BUILD=true; shift ;;
    --protocols)       PROTOCOLS="$2"; shift 2 ;;
    *) echo "[ERROR] Unknown option: $1" >&2; exit 1 ;;
  esac
done

# Derive seed from the base params file (first protocol config or params.json)
SEED=$(${PYTHON} -c "
import json, sys, pathlib
prefix = '${PARAMS_PREFIX}'
protos = '${PROTOCOLS}'.split(',')
# Try first protocol's config file, fall back to params.json
for p in protos:
    key = p.replace('+', '_')
    suffix = '_${SCALE}' if '${SCALE}' != 'S1e5' else ''
    fname = pathlib.Path('${REPO_ROOT}') / f'{prefix}{key}{suffix}.json'
    if not fname.exists():
        fname = pathlib.Path('${REPO_ROOT}') / 'params.json'
    if fname.exists():
        d = json.load(open(fname))
        print(d.get('seed', 42))
        sys.exit(0)
print(42)
" 2>/dev/null || echo "42")

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
GROUP_ID="${TIMESTAMP}_${SCALE}_s${SEED}"
GROUP_DIR="${RUNS_TREATMENT}/${GROUP_ID}"

echo "============================================================"
echo "  Treatment Group Experiment"
echo "  Group folder: ${GROUP_DIR}"
echo "  Scale: ${SCALE}   Seed: ${SEED}"
echo "  Protocols: ${PROTOCOLS}"
echo "  $(date)"
echo "============================================================"

# -- 1. Build ----------------------------------------------------------------
if [ "${SKIP_BUILD}" = false ]; then
  echo ""
  echo "[1/N] Building..."
  cd "${REPO_ROOT}"
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -Wno-dev 2>&1 | tail -3
  cmake --build build --parallel "$(nproc)" 2>&1 | tail -10
  echo "      Build complete."
else
  echo "[1/N] Build skipped."
fi

if [ ! -f "${BINARY}" ]; then
  echo "[ERROR] Binary not found: ${BINARY}" >&2; exit 1
fi

mkdir -p "${GROUP_DIR}"
echo "  Created group folder: ${GROUP_DIR}"

# -- 2. Run each protocol ----------------------------------------------------
IFS=',' read -ra PROTO_LIST <<< "${PROTOCOLS}"
N_PROTO=${#PROTO_LIST[@]}
IDX=1

for PROTO in "${PROTO_LIST[@]}"; do
  # Derive config filename: replace + with _ for filename
  PROTO_KEY="${PROTO//+/_}"
  if [ "${SCALE}" = "S1e5" ]; then
    CONFIG="${REPO_ROOT}/${PARAMS_PREFIX}${PROTO_KEY}.json"
  else
    CONFIG="${REPO_ROOT}/${PARAMS_PREFIX}${PROTO_KEY}_${SCALE}.json"
  fi

  if [ ! -f "${CONFIG}" ]; then
    echo "[ERROR] Config not found: ${CONFIG}" >&2; exit 1
  fi

  # Fail fast if post-treatment tumor-growth switch is missing.
  # This parameter is critical for reproducing Section 5 trajectories.
    if ! ${PYTHON} -c 'import json,sys; cfg=sys.argv[1]; p=json.load(open(cfg));\
  sys.exit((print(f"[ERROR] Missing kc_post_treat in {cfg}"),2)[1]) if "kc_post_treat" not in p else None;\
  _x=p.get("kc_post_treat");\
  sys.exit((print(f"[ERROR] Invalid kc_post_treat in {cfg}: {_x}"),3)[1]) if (not isinstance(_x,(int,float))) else None' "${CONFIG}"; then
    exit 1
  fi

  # BioDynaMo treats --output as a directory; the CSV lands at <dir>/populations.csv
  OUTPUT_DIR="${REPO_ROOT}/output/populations_${PROTO_KEY}.csv"
  OUTPUT_CSV="${OUTPUT_DIR}/populations.csv"

  echo ""
  echo "[$(( IDX + 1 ))/$((N_PROTO + 1))] Protocol: ${PROTO}"
  echo "  Config: ${CONFIG}"

  mkdir -p "${REPO_ROOT}/output"
  START=$(date +%s)
  cd "${REPO_ROOT}"
  "${BINARY}" --config "${CONFIG}" --output "${OUTPUT_DIR}"
  END=$(date +%s)
  DURATION=$(( END - START ))
  echo "  Simulation done in ${DURATION}s → ${OUTPUT_CSV}"

  echo "  Archiving → ${GROUP_DIR}/${PROTO}/"
  "${PYTHON}" scripts/save_run.py \
    --params    "${CONFIG}" \
    --abm       "${OUTPUT_CSV}" \
    --refs      "data-export" \
    --group-dir "${GROUP_DIR}" \
    --duration  "${DURATION}" \
    --note      "${NOTE}"

  IDX=$(( IDX + 1 ))
done

# -- 3. Summary --------------------------------------------------------------
echo ""
echo "============================================================"
echo "  Group complete: ${GROUP_ID}"
echo "  Folder: ${GROUP_DIR}"
echo ""
"${PYTHON}" -c "
import json, pathlib
p = pathlib.Path('${GROUP_DIR}') / 'group_summary.json'
if p.exists():
    s = json.load(open(p))
    print('  Protocol results:')
    for proto, info in s.get('protocols', {}).items():
        c = info.get('final_counts', {}).get('C', '?')
        print(f'    {proto:15s}: C={c}')
"
echo "============================================================"
