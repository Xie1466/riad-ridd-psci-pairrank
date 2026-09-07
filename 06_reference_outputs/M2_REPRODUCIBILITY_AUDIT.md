# M2 reproducibility audit

Status: **EXACT_REPRODUCIBLE**.

The frozen M2 was reconstructed without model changes: features `L_lit,D,C,F1,F4,F5,L`, fold-local scale-only population-SD scaling, no intercept, L2 penalty 1.0, and the original strict four-fold LOSO split. The public-history architecture prior was rebuilt per fold with the original held-out-system overlap exclusions. No tuning, feature substitution, imputation, label change, or candidate reselection occurred.

- Compared predictions: 24/24.
- Direction-exact predictions: 24/24.
- Reconstructed correct comparisons: 19/24 (expected 19/24).
- Maximum absolute probability difference from `07_PAIRWISE_DATASET.csv`: `0.000e+00`.
- Frozen metrics cross-check: pair accuracy `0.791667` versus `19/24 = 0.791667` in `08_PAIRWISE_V2_FINAL_METRICS.md`.
- Fold artifacts: coefficients and scalers are in the two requested CSVs; training/held-out systems, pair probabilities, and construct rankings are in the two additional fold CSVs.

The recorded difference is within deterministic decimal serialization/floating-point tolerance and does not change any probability, direction call, correctness flag, or ranking.
