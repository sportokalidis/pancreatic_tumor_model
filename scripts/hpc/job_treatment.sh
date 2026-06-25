#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# job_treatment.sh — SLURM job: run drug-treatment protocols via Singularity.
#
# Runs the unified flag-gated binary through run_direct.sh --treatment, which
# selects each protocol's config in configs/, patches output_dir, runs the
# binary (BDM_PARAMS), and archives every protocol under ONE group folder in
# runs/treatment/.  All protocols run in parallel inside the single job.
#
# Submit from the HOST (outside the container):
#
#   sbatch scripts/hpc/job_treatment.sh                       # Fig 5: 4 protocols, S1e4
#   PROTOCOLS=abr_csc NOTE="Fig 6 relapse" sbatch scripts/hpc/job_treatment.sh
#   SCALE=S1e5 sbatch scripts/hpc/job_treatment.sh
#
# Env vars (all optional):
#   SCALE        S1e4 (default) | S1e5 | S1e3 — selects configs/params_treat_<p>_<SCALE>.json
#   PROTOCOLS    Comma list (default: acd47,gem,abr,abr_acd47).
#                Use "abr_csc" for the Fig 6 (Abraxane + CSC relapse) scenario.
#   NOTE         Label stored in each run archive
#   BDM_SIF      Path to Singularity image (default: Vega BioDynaMo SIF)
#   SKIP_BUILD   true (default) | false — set false to rebuild binary first
# ---------------------------------------------------------------------------
#SBATCH --job-name=bdm-treat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --time=16:00:00
#SBATCH --output=logs/slurm-treat-%j.out
#SBATCH --error=logs/slurm-treat-%j.err

# SLURM copies the script to a spool dir, so BASH_SOURCE is unreliable.
# REPO_ROOT is injected at submit time; fall back to SLURM_SUBMIT_DIR.
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}}"
REPO_ROOT="$(cd "${REPO_ROOT}" && pwd -P)"   # resolve symlinks → /cephhome/...
SIF="${BDM_SIF:-/ceph/hpc/home/eustavrosp/biodynamo/Singularity.sif}"
SCALE="${SCALE:-S1e4}"
PROTOCOLS="${PROTOCOLS:-acd47,gem,abr,abr_acd47}"
NOTE="${NOTE:-SLURM treatment scale=${SCALE} protocols=${PROTOCOLS}}"

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

echo "[job_treatment] scale=${SCALE}  protocols=${PROTOCOLS}  node=$(hostname)"
echo "                sif=$(basename "${SIF}")  SLURM_JOB_ID=${SLURM_JOB_ID:-local}"

# run_direct.sh --treatment runs every protocol in parallel and archives them
# into one runs/treatment/<group> folder.
RUN_ARGS=(--treatment --parallel
          --scale "${SCALE}"
          --protocols "${PROTOCOLS}"
          --note "${NOTE}")
# Default to skipping build (binary pre-built); set SKIP_BUILD=false to rebuild.
[ "${SKIP_BUILD:-true}" = true ] && RUN_ARGS+=(--skip-build)

"${SING}" exec --cleanenv --bind /cephhome \
  --env "LD_PRELOAD=${REPO_ROOT}/scripts/hpc/fake_numa.so" \
  --env "PYTHONUSERBASE=/cephhome/eustavrosp/.local" \
  "${SIF}" \
  /bin/bash "${REPO_ROOT}/scripts/hpc/run_direct.sh" \
  "${RUN_ARGS[@]}"
