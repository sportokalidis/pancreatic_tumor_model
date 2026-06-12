#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# env.sh — portable environment setup for local dev and HPC clusters.
#
# Source this file at the top of every run script:
#   source "$(dirname "${BASH_SOURCE[0]}")/hpc/env.sh"
#
# Exports: BDM_BUILD  PYTHON  BDM_THISBDM
#
# Override any variable before sourcing to skip auto-detection:
#   export BDM_BUILD=/opt/biodynamo/build
#   source scripts/hpc/env.sh
# ---------------------------------------------------------------------------

# ---------- BioDynaMo installation -----------------------------------------
if [ -z "${BDM_BUILD:-}" ]; then
  # Priority order: common paths, then module system
  for candidate in \
      "${HOME}/biodynamo/build" \
      "${HOME}/biodynamo-build" \
      "/opt/biodynamo/build" \
      "/usr/local/biodynamo/build" \
      "${HOME}/Documents/dev/biodynamo/build"; do
    if [ -f "${candidate}/bin/thisbdm.sh" ]; then
      BDM_BUILD="${candidate}"
      break
    fi
  done

  # Last resort: ask the module system
  if [ -z "${BDM_BUILD:-}" ] && command -v module &>/dev/null; then
    module load biodynamo 2>/dev/null || true
    if command -v biodynamo &>/dev/null; then
      BDM_BUILD="$(dirname "$(dirname "$(command -v biodynamo)")")"
    fi
  fi
fi

if [ -z "${BDM_BUILD:-}" ]; then
  echo "[env.sh] ERROR: Cannot find BioDynaMo installation." >&2
  echo "         Set BDM_BUILD=/path/to/biodynamo/build before sourcing env.sh" >&2
  return 1 2>/dev/null || exit 1
fi
export BDM_BUILD

BDM_THISBDM="${BDM_BUILD}/bin/thisbdm.sh"
if [ ! -f "${BDM_THISBDM}" ]; then
  echo "[env.sh] ERROR: thisbdm.sh not found at ${BDM_THISBDM}" >&2
  return 1 2>/dev/null || exit 1
fi
export BDM_THISBDM

{ set +eu; source "${BDM_THISBDM}" >/dev/null 2>&1; set -eu; } || true
echo "[env.sh] BioDynaMo: ${BDM_BUILD}"

# ---------- Python ≥ 3.7 ---------------------------------------------------
if [ -z "${PYTHON:-}" ]; then
  for candidate in \
      "$(which python3 2>/dev/null)" \
      "$(which python 2>/dev/null)" \
      "${HOME}/.pyenv/versions/3.9.1/bin/python3" \
      "/usr/bin/python3" \
      "/usr/local/bin/python3"; do
    [ -z "${candidate}" ] && continue
    if "${candidate}" -c "import sys; sys.exit(0 if sys.version_info>=(3,7) else 1)" 2>/dev/null; then
      PYTHON="${candidate}"
      break
    fi
  done

  if [ -z "${PYTHON:-}" ] && command -v module &>/dev/null; then
    module load python 2>/dev/null || module load python3 2>/dev/null || true
    PYTHON="$(which python3 2>/dev/null || which python 2>/dev/null || echo '')"
  fi
fi

if [ -z "${PYTHON:-}" ]; then
  echo "[env.sh] ERROR: Cannot find Python >= 3.7." >&2
  echo "         Set PYTHON=/path/to/python3 before sourcing env.sh" >&2
  return 1 2>/dev/null || exit 1
fi
export PYTHON
echo "[env.sh] Python: ${PYTHON}  ($(${PYTHON} --version 2>&1))"
