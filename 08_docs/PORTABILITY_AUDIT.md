# Numerical portability audit

The frozen scientific reference was generated and package-tested under the pinned runtime recorded in `01_environment/` (Python 3.12.3, NumPy 1.26.4, SciPy 1.11.4, Matplotlib 3.6.3, PyYAML 6.0.1). The release-generation record reports quick reproduction PASS and exact reconstruction of all frozen M2 probabilities.

An independent public-package audit was additionally run on a deliberately different modern stack (Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8, PyYAML 6.0.3). All manuscript metrics and all 96 M0–M3 pair directions were unchanged. The largest probability-level floating-point deviation from the frozen reference was 1.0754e-08.

Accordingly, the public wrapper implements two verification modes without changing any model, data, features, penalties or folds:

- **EXACT_PINNED_ENV:** tolerance 1e-12.
- **COMPATIBILITY_ENV:** tolerance 1e-7, while requiring identical pair directions and identical reported manuscript metrics.

This distinction prevents harmless BLAS/NumPy floating-point variation from being misreported as a scientific reproduction failure while preserving the strict exact-environment gate.
