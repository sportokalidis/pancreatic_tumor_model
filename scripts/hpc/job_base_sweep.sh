#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# job_base_sweep.sh — SLURM array job: run the base model with N seeds
#                     in parallel via Singularity.
#
# Submit from the HOST (outside the container):
#
#   # 5 seeds (array tasks 0-4) at S1e4
#   sbatch scripts/hpc/job_base_sweep.sh
#
#   # 5 seeds at S1e3 (more agents)
#   SCALE=S1e3 sbatch scripts/hpc/job_base_sweep.sh
#
#   # Custom number of seeds (override --array on the command line)
#   sbatch --array=0-9 scripts/hpc/job_base_sweep.sh
#
# Build first (once):
#   singularity exec --cleanenv Singularity.sif bash scripts/hpc/build_hpc.sh
#
# Env vars (all optional):
#   SCALE     S1e5 | S1e4 (default) | S1e3
#   NOTE      Label prefix (seed appended automatically)
#   BDM_SIF   Path to Singularity image (default: <repo>/Singularity.sif)
# ---------------------------------------------------------------------------
#SBATCH --job-name=bdm-sweep
#SBATCH --array=0-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm-sweep-%A-%a.out
#SBATCH --error=logs/slurm-sweep-%A-%a.err

# SLURM copies the script to a spool dir, so BASH_SOURCE is unreliable.
# REPO_ROOT is injected by submit_base_scales.sh; fall back to SLURM_SUBMIT_DIR.
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}}"
# Resolve any symlinks so the path matches /cephhome/... inside the container
REPO_ROOT="$(cd "${REPO_ROOT}" && pwd -P)"
SIF="${BDM_SIF:-/ceph/hpc/home/eustavrosp/biodynamo/Singularity.sif}"
SCALE="${SCALE:-S1e4}"

# One seed per array task — extend the list for more tasks
SEEDS=(42 123 456 789 1001 2024 3141 5926 9999 7777)
SEED="${SEEDS[$SLURM_ARRAY_TASK_ID]}"
NOTE="${NOTE:-sweep ${SCALE} seed=${SEED} task=${SLURM_ARRAY_TASK_ID}}"

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

mkdir -p "${REPO_ROOT}/logs"
echo "[job_base_sweep] array=${SLURM_ARRAY_JOB_ID:-local}[${SLURM_ARRAY_TASK_ID:-0}]  scale=${SCALE}  seed=${SEED}  sing=${SING}"

"${SING}" exec --cleanenv --bind /cephhome \
  --env "LD_PRELOAD=${REPO_ROOT}/scripts/hpc/fake_numa.so" \
  --env "PYTHONUSERBASE=/cephhome/eustavrosp/.local" \
  "${SIF}" \
  /bin/bash "${REPO_ROOT}/scripts/hpc/run_direct.sh" \
  --scale      "${SCALE}" \
  --seed       "${SEED}" \
  --threads    "${SLURM_CPUS_PER_TASK:-16}" \
  --skip-build \
  --note       "${NOTE}"
