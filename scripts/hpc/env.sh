#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# env.sh — portable environment setup for local dev, HPC, and Apptainer.
#
# Source this file at the top of every run script:
#   source "$(dirname "${BASH_SOURCE[0]}")/hpc/env.sh"
#
# Exports: BDM_BUILD  PYTHON
#
# Override before sourcing to skip auto-detection:
#   export BDM_BUILD=/opt/biodynamo/build
#   source scripts/hpc/env.sh
# ---------------------------------------------------------------------------

# ---------- Apptainer / Singularity detection ------------------------------
# Inside a container, BioDynaMo shared libs are already in LD_LIBRARY_PATH
# and the binary is on PATH.  We only need to locate BDM_BUILD for cmake
# (build_hpc.sh) and for thisbdm.sh (optional, only needed if libs aren't
# already loaded).
_IN_CONTAINER=false
if [ -n "${APPTAINER_CONTAINER:-}" ] || [ -n "${SINGULARITY_CONTAINER:-}" ]; then
  _IN_CONTAINER=true
fi

# ---------- BioDynaMo installation -----------------------------------------
if [ -z "${BDM_BUILD:-}" ]; then
  # Priority: Apptainer standard prefix, then user builds, then modules
  for candidate in \
      "/biodynamo/build" \
      "/biodynamo" \
      "/opt/biodynamo" \
      "/opt/biodynamo/build" \
      "/usr/local/biodynamo" \
      "${HOME}/biodynamo/build" \
      "${HOME}/biodynamo-build" \
      "${HOME}/Documents/dev/biodynamo/build"; do
    # Accept the path if it has either thisbdm.sh (user build) or the BDM
    # cmake config (system install with no thisbdm.sh needed)
    if [ -f "${candidate}/bin/thisbdm.sh" ] || \
       [ -f "${candidate}/lib/cmake/BioDynaMo/BioDynaMoConfig.cmake" ] || \
       [ -f "${candidate}/cmake/BioDynaMoConfig.cmake" ]; then
      BDM_BUILD="${candidate}"
      break
    fi
  done

  # Module system fallback
  if [ -z "${BDM_BUILD:-}" ] && command -v module &>/dev/null; then
    module load biodynamo 2>/dev/null || true
    if command -v biodynamo &>/dev/null; then
      BDM_BUILD="$(dirname "$(dirname "$(command -v biodynamo)")")"
    fi
  fi
fi

if [ -z "${BDM_BUILD:-}" ]; then
  echo "[env.sh] ERROR: Cannot find BioDynaMo installation." >&2
  echo "         Set BDM_BUILD=/path/to/biodynamo before sourcing env.sh" >&2
  return 1 2>/dev/null || exit 1
fi
export BDM_BUILD

# Always source thisbdm.sh when it exists — needed even inside a container
# when launched with --cleanenv, which strips LD_LIBRARY_PATH and ROOTSYS.
_THISBDM="${BDM_BUILD}/bin/thisbdm.sh"
if [ -f "${_THISBDM}" ]; then
  { set +eu; source "${_THISBDM}" >/dev/null 2>&1; set -eu; } || true
  echo "[env.sh] BioDynaMo (sourced thisbdm.sh): ${BDM_BUILD}"
else
  echo "[env.sh] BioDynaMo (no thisbdm.sh, assuming pre-loaded): ${BDM_BUILD}"
fi

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
