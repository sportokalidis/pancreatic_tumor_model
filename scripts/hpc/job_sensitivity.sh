#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# job_sensitivity.sh — SLURM ARRAY job for the sensitivity-analysis sweep.
#
# Each array task runs ONE pre-generated config from a sweep directory:
#   <SWEEP_DIR>/configs/{oat,sample}_NNNN.json
# (configs are produced by scripts/sensitivity/gen_oat.py / gen_lhs.py).
#
# The config already carries the fixed SA settings (S=1e4, dt=24h, seed=42) and
# its own absolute output_dir, so the task just runs the binary under
# Singularity with BDM_PARAMS pointing at its config. No per-run post-processing
# — collect.py + analyze_*.py run ONCE after the whole array finishes.
#
# Submit from the HOST (outside the container) — usually via submit_sensitivity.sh:
#   SWEEP_DIR=/abs/path/to/sweep sbatch --array=0-70 scripts/hpc/job_sensitivity.sh
#
# Env vars:
#   SWEEP_DIR   (required) absolute path to the sweep dir holding configs/
#   BDM_SIF     Singularity image (default: /ceph/hpc/home/eustavrosp/biodynamo/Singularity.sif)
#   REPO_ROOT   repo root (default: SLURM_SUBMIT_DIR)
#   THREADS     OMP_NUM_THREADS (default: SLURM cpus-per-task)
# ---------------------------------------------------------------------------
#SBATCH --job-name=bdm-sens
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=46G
#SBATCH --time=01:00:00
#SBATCH --partition=cpu-st
# Partition: cpu-st (longcpu was down for HW maintenance 2026-06 — cn0012/cn0014
# draining; switch back with `sbatch --partition=longcpu` or PARTITION= once it
# returns). Override via the PARTITION env var in submit_sensitivity.sh.
# Memory: baseline S=1e4 dt=24h uses ~27 GB; SA perturbs growth/capacity params
# up to 2x, so high-growth samples need more headroom (~40 GB). 46G/24 cores =
# 1962 MB/core stays under the cluster's 2048 MB/core policy cap (32G/16=2048
# tripped it). cpu/cpu-st nodes are 256 cores / ~251 GB, so 24-core tasks pack fine.
#SBATCH --output=logs/slurm-sens-%A_%a.out
#SBATCH --error=logs/slurm-sens-%A_%a.err
#SBATCH --export=ALL
# Partition: longcpu by default (cpu partition is often busy with base runs).
# Each task is short (~6-10 min); override with `sbatch --partition=...` or the
# PARTITION env var in submit_sensitivity.sh.

set -euo pipefail

# Ensure module system is available on the compute node
for _m in /etc/profile.d/modules.sh /usr/share/Modules/init/bash \
          /opt/modules/init/bash /usr/local/Modules/init/bash; do
  [ -f "${_m}" ] && { source "${_m}" 2>/dev/null; break; }
done

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}}"
REPO_ROOT="$(cd "${REPO_ROOT}" && pwd -P)"
SIF="${BDM_SIF:-/ceph/hpc/home/eustavrosp/biodynamo/Singularity.sif}"
BINARY="${REPO_ROOT}/build/pancreatic_tumor_new"

if [ -z "${SWEEP_DIR:-}" ]; then
  echo "[ERROR] SWEEP_DIR not set." >&2; exit 1
fi
SWEEP_DIR="$(cd "${SWEEP_DIR}" && pwd -P)"
CFG_DIR="${SWEEP_DIR}/configs"
[ -d "${CFG_DIR}" ] || { echo "[ERROR] no configs/ in ${SWEEP_DIR}" >&2; exit 1; }
[ -f "${BINARY}" ]  || { echo "[ERROR] binary not found: ${BINARY} (build first)" >&2; exit 1; }
[ -e "${SIF}" ]     || { echo "[ERROR] SIF not found: ${SIF}" >&2; exit 1; }

# Map array task id -> Nth config (sorted, stable).
mapfile -t CONFIGS < <(ls -1 "${CFG_DIR}"/*.json | sort)
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
if [ "${TASK_ID}" -ge "${#CONFIGS[@]}" ]; then
  echo "[ERROR] task ${TASK_ID} >= ${#CONFIGS[@]} configs" >&2; exit 1
fi
CONFIG="${CONFIGS[$TASK_ID]}"

# Locate singularity/apptainer
_find_sing() { command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || echo ""; }
SING="$(_find_sing)"
if [ -z "${SING}" ] && command -v module &>/dev/null; then
  module load singularity 2>/dev/null || module load apptainer 2>/dev/null || true
  SING="$(_find_sing)"
fi
[ -n "${SING}" ] || { echo "[ERROR] singularity/apptainer not found." >&2; exit 1; }

mkdir -p "${REPO_ROOT}/logs"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"

echo "================================================"
echo "  sensitivity task ${TASK_ID}/${#CONFIGS[@]}"
echo "  config : $(basename "${CONFIG}")"
echo "  threads: ${THREADS}   node: $(hostname)   $(date)"
echo "================================================"

START=$(date +%s)
"${SING}" exec --cleanenv --bind /cephhome \
  --env "LD_PRELOAD=${REPO_ROOT}/scripts/hpc/fake_numa.so" \
  --env "PYTHONUSERBASE=/cephhome/eustavrosp/.local" \
  --env "BDM_PARAMS=${CONFIG}" \
  --env "OMP_NUM_THREADS=${THREADS}" \
  --env "OMP_PROC_BIND=true" \
  "${SIF}" \
  /bin/bash -c "cd '${REPO_ROOT}' && '${BINARY}'"
END=$(date +%s)
echo "  done in $(( END - START ))s"
