#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Build the ABM binary on HPC.  Run this once after cloning/pulling.
#
# Usage:
#   bash scripts/hpc/build_hpc.sh [--clean]
#
# Options:
#   --clean    Remove the build/ directory before building
#
# Environment variables:
#   BDM_BUILD  Path to BioDynaMo build dir (auto-detected if unset)
#   BUILD_TYPE Release (default) | Debug
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO_ROOT}/scripts/hpc/env.sh"

BUILD_TYPE="${BUILD_TYPE:-Release}"
CLEAN=false

for arg in "$@"; do
  case "$arg" in
    --clean) CLEAN=true ;;
    *) echo "[ERROR] Unknown option: $arg"; exit 1 ;;
  esac
done

if [ "${CLEAN}" = true ]; then
  echo "[build] Removing old build directory..."
  rm -rf "${REPO_ROOT}/build"
fi

echo "[build] Configuring (${BUILD_TYPE})..."
cd "${REPO_ROOT}"
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DBIODYNAMO_ROOT="${BDM_BUILD}" \
  -Wno-dev \
  2>&1 | tail -5

echo "[build] Compiling ($(nproc) cores)..."
cmake --build build --parallel "$(nproc)"

BINARY="${REPO_ROOT}/build/pancreatic_tumor_new"
if [ -f "${BINARY}" ]; then
  echo "[build] Success: ${BINARY}"
  "${BINARY}" --version 2>/dev/null || true
else
  echo "[build] ERROR: binary not found after build" >&2; exit 1
fi
