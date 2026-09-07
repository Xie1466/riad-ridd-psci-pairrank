#!/usr/bin/env bash
set -euo pipefail
HERE="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  echo "Usage: bash 04_scripts/run_full_reproduction.sh TRAJECTORY_ROOT"
  echo "Fail-closed trajectory-to-PSCI entry point; exact trajectory deposit is not bundled."
  exit 0
fi
if [[ $# -lt 1 ]]; then
  echo "FULL_TRAJECTORY_REPRODUCTION_NOT_AVAILABLE: provide a future standardized trajectory DOI/mount path" >&2
  exit 2
fi
"${PYTHON:-python3}" "$HERE/01_build_psci_features.py" --trajectory-root "$1"
