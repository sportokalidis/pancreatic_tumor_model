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
#   SCALE        S1e5 | S1e4 (default) | S1e3 | S1e4_dt4h | S1e4_dt6h | ...
#   SEED         Override random seed (default: from params.json)
#   NOTE         Label stored in run archive
#   BDM_SIF      Path to Singularity image (default: <repo>/Singularity.sif)
#   THREADS      OMP_NUM_THREADS (default: all cpus-per-task)
#   SKIP_BUILD   true (default) | false — set false to rebuild binary before run
# ---------------------------------------------------------------------------
#SBATCH --job-name=bdm-base
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm-base-%j.out
#SBATCH --error=logs/slurm-base-%j.err
#SBATCH --export=ALL

# Ensure module system is initialised on the compute node
for _m in /etc/profile.d/modules.sh /usr/share/Modules/init/bash \
          /opt/modules/init/bash /usr/local/Modules/init/bash; do
  [ -f "${_m}" ] && { source "${_m}" 2>/dev/null; break; }
done

# SLURM copies the script to a spool dir, so BASH_SOURCE is unreliable.
# REPO_ROOT is injected by submit_base_scales.sh; fall back to SLURM_SUBMIT_DIR.
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}}"
# Resolve any symlinks so the path matches /cephhome/... inside the container
REPO_ROOT="$(cd "${REPO_ROOT}" && pwd -P)"
SIF="${BDM_SIF:-/ceph/hpc/home/eustavrosp/biodynamo/Singularity.sif}"
SCALE="${SCALE:-S1e4}"
NOTE="${NOTE:-SLURM base scale=${SCALE}}"

if [ ! -e "${SIF}" ]; then
  echo "[ERROR] Singularity image not found: ${SIF}" >&2
  echo "        Set BDM_SIF=/path/to/image.sif before submitting." >&2
  exit 1
fi

# Locate singularity/apptainer — compute nodes may need a module load
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

EXTRA_ARGS=()
[ -n "${SEED:-}" ] && EXTRA_ARGS+=(--seed "${SEED}")
# Default to SLURM allocation; override with THREADS env var
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
EXTRA_ARGS+=(--threads "${THREADS}")

echo "[job_base] scale=${SCALE}  sif=$(basename "${SIF}")  node=$(hostname)"
echo "           sing=${SING}  SLURM_JOB_ID=${SLURM_JOB_ID:-local}"

RUN_ARGS=(--scale "${SCALE}" --note "${NOTE}" "${EXTRA_ARGS[@]}")
# Add explicit output directory if provided (from validation suite runner)
[ -n "${OUTPUT_DIR:-}" ] && RUN_ARGS+=(--output-dir "${OUTPUT_DIR}")
# Default to skipping build (binary pre-built); set SKIP_BUILD=false to rebuild.
[ "${SKIP_BUILD:-true}" = true ] && RUN_ARGS+=(--skip-build)

"${SING}" exec --cleanenv --bind /cephhome \
  --env "LD_PRELOAD=${REPO_ROOT}/scripts/hpc/fake_numa.so" \
  --env "PYTHONUSERBASE=/cephhome/eustavrosp/.local" \
  "${SIF}" \
  /bin/bash "${REPO_ROOT}/scripts/hpc/run_direct.sh" \
  "${RUN_ARGS[@]}"
