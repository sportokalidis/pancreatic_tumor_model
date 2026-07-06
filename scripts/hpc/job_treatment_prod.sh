#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# job_treatment_prod.sh — SLURM job for ONE treatment protocol × ONE seed,
# run in the same per-seed "production" style as job_base.sh.
#
# Normally submitted by scripts/run_treatment_production.py (one job per
# protocol per seed), but can be launched directly:
#
#   CONFIG=configs/params_treat_gem_S1e4_E10x_dt6h.json \
#   OUTPUT_DIR=runs/treatment_prod/mysuite/gem_s1 \
#   SEED=1 sbatch scripts/hpc/job_treatment_prod.sh
#
# Env vars:
#   CONFIG       (required) self-contained protocol config (see make_treatment_configs.py)
#   OUTPUT_DIR   (required) per-seed run directory to archive into
#   SEED         Random seed (default: from CONFIG)
#   NOTE         Label stored in the run archive
#   BDM_SIF      Path to Singularity image (default: Vega BioDynaMo SIF)
#   THREADS      OMP_NUM_THREADS (default: cpus-per-task)
#   SKIP_BUILD   true (default) | false
# ---------------------------------------------------------------------------
#SBATCH --job-name=bdm-treat-prod
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm-treatprod-%j.out
#SBATCH --error=logs/slurm-treatprod-%j.err
#SBATCH --export=ALL

for _m in /etc/profile.d/modules.sh /usr/share/Modules/init/bash \
          /opt/modules/init/bash /usr/local/Modules/init/bash; do
  [ -f "${_m}" ] && { source "${_m}" 2>/dev/null; break; }
done

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}}"
REPO_ROOT="$(cd "${REPO_ROOT}" && pwd -P)"   # resolve symlinks -> /cephhome/...
SIF="${BDM_SIF:-/ceph/hpc/home/eustavrosp/biodynamo/Singularity.sif}"

if [ -z "${CONFIG:-}" ]; then
  echo "[ERROR] CONFIG env var is required (path to protocol config)." >&2; exit 1
fi
if [ -z "${OUTPUT_DIR:-}" ]; then
  echo "[ERROR] OUTPUT_DIR env var is required." >&2; exit 1
fi
# Resolve CONFIG relative to the repo if not absolute
case "${CONFIG}" in
  /*) : ;;
  *)  CONFIG="${REPO_ROOT}/${CONFIG}" ;;
esac
[ -f "${CONFIG}" ] || { echo "[ERROR] Config not found: ${CONFIG}" >&2; exit 1; }

NOTE="${NOTE:-SLURM treatment-prod $(basename "${CONFIG}") seed=${SEED:-cfg}}"

if [ ! -e "${SIF}" ]; then
  echo "[ERROR] Singularity image not found: ${SIF}" >&2
  echo "        Set BDM_SIF=/path/to/image.sif before submitting." >&2
  exit 1
fi

_find_sing() {
  command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || echo ""
}
SING="$(_find_sing)"
if [ -z "${SING}" ] && command -v module &>/dev/null; then
  module load singularity 2>/dev/null || module load apptainer 2>/dev/null || true
  SING="$(_find_sing)"
fi
if [ -z "${SING}" ]; then
  echo "[ERROR] singularity/apptainer not found. Load the module and resubmit." >&2
  exit 1
fi

mkdir -p "${REPO_ROOT}/logs" "${OUTPUT_DIR}"

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
RUN_ARGS=(--config-file "${CONFIG}" --output-dir "${OUTPUT_DIR}"
          --threads "${THREADS}" --note "${NOTE}")
[ -n "${SEED:-}" ] && RUN_ARGS+=(--seed "${SEED}")
[ "${SKIP_BUILD:-true}" = true ] && RUN_ARGS+=(--skip-build)

echo "[job_treatment_prod] cfg=$(basename "${CONFIG}")  seed=${SEED:-cfg}  node=$(hostname)"
echo "                     out=${OUTPUT_DIR}  SLURM_JOB_ID=${SLURM_JOB_ID:-local}"

"${SING}" exec --cleanenv --bind /cephhome \
  --env "LD_PRELOAD=${REPO_ROOT}/scripts/hpc/fake_numa.so" \
  --env "PYTHONUSERBASE=/cephhome/eustavrosp/.local" \
  "${SIF}" \
  /bin/bash "${REPO_ROOT}/scripts/hpc/run_direct.sh" \
  "${RUN_ARGS[@]}"
