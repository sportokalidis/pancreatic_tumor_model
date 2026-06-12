#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# job_base.sh — SLURM job for a single base-model run via Singularity.
#
# Submit from the HOST (outside the container):
#
#   sbatch scripts/hpc/job_base.sh
#   SCALE=S1e3 sbatch scripts/hpc/job_base.sh
#   SCALE=S1e4 SEED=123 NOTE="rerun" sbatch scripts/hpc/job_base.sh
#
# Build first (once):
#   singularity exec --cleanenv Singularity.sif bash scripts/hpc/build_hpc.sh
#
# Env vars (all optional):
#   SCALE     S1e5 | S1e4 (default) | S1e3
#   SEED      Override random seed (default: from params.json)
#   NOTE      Label stored in run archive
#   BDM_SIF   Path to Singularity image (default: <repo>/Singularity.sif)
#   THREADS   OMP_NUM_THREADS (default: all cpus-per-task)
# ---------------------------------------------------------------------------
#SBATCH --job-name=bdm-base
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm-base-%j.out
#SBATCH --error=logs/slurm-base-%j.err

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SIF="${BDM_SIF:-/ceph/hpc/home/eustavrosp/biodynamo}"
SCALE="${SCALE:-S1e4}"
NOTE="${NOTE:-SLURM base scale=${SCALE}}"

if [ ! -f "${SIF}" ]; then
  echo "[ERROR] Singularity image not found: ${SIF}" >&2
  echo "        Set BDM_SIF=/path/to/image.sif before submitting." >&2
  exit 1
fi

mkdir -p "${REPO_ROOT}/logs"

EXTRA_ARGS=()
[ -n "${SEED:-}"    ] && EXTRA_ARGS+=(--seed    "${SEED}")
[ -n "${THREADS:-}" ] && EXTRA_ARGS+=(--threads "${THREADS}")

echo "[job_base] scale=${SCALE}  sif=$(basename "${SIF}")  node=$(hostname)"
echo "           SLURM_JOB_ID=${SLURM_JOB_ID:-local}"

singularity exec --cleanenv "${SIF}" \
  bash "${REPO_ROOT}/scripts/hpc/run_direct.sh" \
  --scale      "${SCALE}" \
  --skip-build \
  --note       "${NOTE}" \
  "${EXTRA_ARGS[@]}"
