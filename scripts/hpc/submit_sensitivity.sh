#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit_sensitivity.sh — generate SA configs and submit the SLURM array.
#
# Host-side convenience wrapper around gen_{oat,lhs}.py + job_sensitivity.sh.
# Generators run INSIDE Singularity (LHS needs numpy). Prints the post-run
# collect + analyze command to run once the array completes.
#
# Usage (from repo root, on the login node):
#   bash scripts/hpc/submit_sensitivity.sh oat
#   bash scripts/hpc/submit_sensitivity.sh lhs [N_SAMPLES]   # default 250
#
# Build the binary first if needed:
#   singularity exec --cleanenv --bind /cephhome $BDM_SIF bash scripts/hpc/build_hpc.sh
#
# Env vars:
#   BDM_SIF        Singularity image (default: /ceph/hpc/home/eustavrosp/biodynamo/Singularity.sif)
#   PYTHON         container python (default: /opt/.pyenv/versions/3.9.1/bin/python3)
#   PARTITION      SLURM partition (default: cpu-st — longcpu was down for HW maint.)
#   MAX_CONCURRENT cap simultaneous array tasks, e.g. 30 (default: unset = no cap)
#   NO_SUBMIT      set to 1 to only generate configs (skip sbatch)
# ---------------------------------------------------------------------------
set -euo pipefail

MODE="${1:-}"
N_SAMPLES="${2:-250}"
[ "${MODE}" = "oat" ] || [ "${MODE}" = "lhs" ] || {
  echo "Usage: $0 {oat|lhs} [n_samples]" >&2; exit 1; }

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
SIF="${BDM_SIF:-/ceph/hpc/home/eustavrosp/biodynamo/Singularity.sif}"
PYTHON="${PYTHON:-/opt/.pyenv/versions/3.9.1/bin/python3}"
PARTITION="${PARTITION:-cpu-st}"
[ -e "${SIF}" ] || { echo "[ERROR] SIF not found: ${SIF}" >&2; exit 1; }

STAMP="$(date +%Y%m%d_%H%M%S)"
SWEEP_DIR="${REPO_ROOT}/runs/sensitivity/${STAMP}_${MODE}"
mkdir -p "${SWEEP_DIR}" "${REPO_ROOT}/logs"

_sing() {
  singularity exec --cleanenv --bind /cephhome \
    --env "PYTHONUSERBASE=/cephhome/eustavrosp/.local" \
    "${SIF}" "$@"
}

echo "[1/2] Generating ${MODE} configs -> ${SWEEP_DIR}"
if [ "${MODE}" = "oat" ]; then
  _sing "${PYTHON}" "${REPO_ROOT}/scripts/sensitivity/gen_oat.py" --out-dir "${SWEEP_DIR}"
else
  _sing "${PYTHON}" "${REPO_ROOT}/scripts/sensitivity/gen_lhs.py" \
    --samples "${N_SAMPLES}" --out-dir "${SWEEP_DIR}"
fi

N_CFG=$(ls -1 "${SWEEP_DIR}/configs"/*.json | wc -l)
LAST=$(( N_CFG - 1 ))
echo "      ${N_CFG} configs."

ANALYZE="analyze_oat.py"
[ "${MODE}" = "lhs" ] && ANALYZE="analyze_prcc.py"
POST_CMD="singularity exec --cleanenv --bind /cephhome --env PYTHONUSERBASE=/cephhome/eustavrosp/.local ${SIF} \\
  ${PYTHON} ${REPO_ROOT}/scripts/sensitivity/collect.py ${SWEEP_DIR} && \\
singularity exec --cleanenv --bind /cephhome --env PYTHONUSERBASE=/cephhome/eustavrosp/.local ${SIF} \\
  ${PYTHON} ${REPO_ROOT}/scripts/sensitivity/${ANALYZE} ${SWEEP_DIR}"

ARRAY_SPEC="0-${LAST}"
[ -n "${MAX_CONCURRENT:-}" ] && ARRAY_SPEC="${ARRAY_SPEC}%${MAX_CONCURRENT}"

if [ "${NO_SUBMIT:-0}" = "1" ]; then
  echo "[2/2] NO_SUBMIT=1 — skipping sbatch."
else
  echo "[2/2] Submitting array ${ARRAY_SPEC} on partition '${PARTITION}'"
  sbatch --array="${ARRAY_SPEC}" --partition="${PARTITION}" \
    --export=ALL,REPO_ROOT="${REPO_ROOT}",SWEEP_DIR="${SWEEP_DIR}",BDM_SIF="${SIF}" \
    "${REPO_ROOT}/scripts/hpc/job_sensitivity.sh"
fi

echo ""
echo "================================================================"
echo "  Sweep dir: ${SWEEP_DIR}"
echo "  When the array finishes, collect + analyze with:"
echo ""
echo "${POST_CMD}"
echo "================================================================"
