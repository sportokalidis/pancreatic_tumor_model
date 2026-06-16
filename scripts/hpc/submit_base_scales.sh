#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit_base_scales.sh — submit base-model seed sweeps for S=1e4 and S=1e3.
#
# Usage (from repo root):
#   bash scripts/hpc/submit_base_scales.sh              # 5 seeds × 2 scales
#   bash scripts/hpc/submit_base_scales.sh --seeds 10   # 10 seeds × 2 scales
#   bash scripts/hpc/submit_base_scales.sh --s1e4-only
#   bash scripts/hpc/submit_base_scales.sh --s1e3-only
#
# Results land in:  runs/<timestamp>_<SCALE>_s<seed>/
# Log files:        logs/slurm-sweep-<JOBID>-<TASKID>.{out,err}
#
# Resource notes:
#   S=1e4  ~23 k initial cells — 32 G / 8 h (ran ~60 days in 4 h; 8 h is safe)
#   S=1e3  ~230 k initial cells, dt=30 min (3× speedup) — 128 CPUs / 24 h / 240 G
#           Estimated ~20 h for 100 days (was 130 h at dt=10 with 64 CPUs).
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
JOB_SWEEP="${REPO_ROOT}/scripts/hpc/job_base_sweep.sh"

N_SEEDS=5
RUN_S1E4=true
RUN_S1E3=true

while [[ $# -gt 0 ]]; do
  case "$1" in
    --seeds)      N_SEEDS="$2"; shift 2 ;;
    --s1e4-only)  RUN_S1E3=false; shift ;;
    --s1e3-only)  RUN_S1E4=false; shift ;;
    *) echo "[ERROR] Unknown option: $1"; exit 1 ;;
  esac
done

ARRAY_END=$(( N_SEEDS - 1 ))

mkdir -p "${REPO_ROOT}/logs"
echo "Submitting base-model sweeps: N_SEEDS=${N_SEEDS}  (array 0-${ARRAY_END})"
echo ""

# ---- S=1e4 -----------------------------------------------------------------
if [ "${RUN_S1E4}" = true ]; then
  JID4=$(sbatch \
    --job-name=bdm-S1e4 \
    --array="0-${ARRAY_END}" \
    --cpus-per-task=16 \
    --mem=32G \
    --time=08:00:00 \
    --output="${REPO_ROOT}/logs/slurm-S1e4-%A-%a.out" \
    --error="${REPO_ROOT}/logs/slurm-S1e4-%A-%a.err" \
    --export=ALL,SCALE=S1e4,REPO_ROOT="${REPO_ROOT}" \
    "${JOB_SWEEP}" | awk '{print $NF}')
  echo "  S=1e4  → job ${JID4}  (${N_SEEDS} tasks, 32G, 8h, cpu)"
fi

# ---- S=1e3 -----------------------------------------------------------------
# dt=30 min (3× fewer steps) + 128 CPUs = ~6× total speedup over previous runs.
# Estimated ~20 h for 100 days; cpu partition (48 h max) is sufficient.
if [ "${RUN_S1E3}" = true ]; then
  JID3=$(sbatch \
    --job-name=bdm-S1e3 \
    --array="0-${ARRAY_END}" \
    --cpus-per-task=128 \
    --mem=240G \
    --time=24:00:00 \
    --output="${REPO_ROOT}/logs/slurm-S1e3-%A-%a.out" \
    --error="${REPO_ROOT}/logs/slurm-S1e3-%A-%a.err" \
    --export=ALL,SCALE=S1e3,REPO_ROOT="${REPO_ROOT}" \
    "${JOB_SWEEP}" | awk '{print $NF}')
  echo "  S=1e3  → job ${JID3}  (${N_SEEDS} tasks, 240G, 24h, 128 CPUs)"
fi

echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Logs in:       ${REPO_ROOT}/logs/"
echo "Results in:    ${REPO_ROOT}/runs/"
