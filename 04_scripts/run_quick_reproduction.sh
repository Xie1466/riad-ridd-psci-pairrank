#!/usr/bin/env bash
set -euo pipefail
HERE="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
ROOT="$HERE/.."
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  echo "Usage: bash 04_scripts/run_quick_reproduction.sh [OUTPUT_DIR]"
  echo "Exact mode uses the pinned Python/NumPy/SciPy/Matplotlib/PyYAML environment."
  echo "Compatible environments use a 1e-7 probability tolerance while requiring identical pair directions and manuscript metrics."
  exit 0
fi
OUT="${1:-$ROOT/reproduction_output}"
PYTHON_BIN="${PYTHON:-python3}"
EXACT="$($PYTHON_BIN "$ROOT/01_environment/check_runtime.py" --exact-flag)"
if [[ "$EXACT" == "1" ]]; then
  TOL="1e-12"
  MODE="EXACT_PINNED_ENV"
else
  TOL="1e-7"
  MODE="COMPATIBILITY_ENV"
  echo "WARNING: runtime differs from the pinned reference environment; using compatibility tolerance $TOL." >&2
  "$PYTHON_BIN" "$ROOT/01_environment/check_runtime.py" >&2
fi
"$PYTHON_BIN" "$HERE/00_validate_inputs.py"
"$PYTHON_BIN" "$HERE/01_build_psci_features.py" --output-dir "$OUT"
"$PYTHON_BIN" "$HERE/02_fit_literature_prior.py" --output-dir "$OUT"
"$PYTHON_BIN" "$HERE/03_run_pairrank_loso.py" --output-dir "$OUT"
"$PYTHON_BIN" "$HERE/04_reproduce_metrics.py" --output-dir "$OUT"
"$PYTHON_BIN" "$HERE/05_reproduce_figures.py" --output-dir "$OUT"
"$PYTHON_BIN" "$HERE/06_verify_reference_outputs.py" --output-dir "$OUT" --tolerance "$TOL"
echo "QUICK_REPRODUCTION_PASS mode=$MODE tolerance=$TOL"
