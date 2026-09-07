# Final public-release audit

## Status

**CODE/DATA RELEASE READY FOR ACS/Zenodo DEPOSIT: YES, subject only to DOI/URL insertion.**

## Scientific integrity

- Frozen PSCI values: unchanged.
- Literature prior inputs: unchanged.
- M0–M3 model definitions/hyperparameters/splits: unchanged.
- Frozen predictions/reference probabilities: unchanged.
- Reported M2 metrics: unchanged (19/24; 79.2%; Spearman 0.654846; Kendall tau-b 0.521749; Top1 3/4; true-best-in-predicted-Top2 3/4; exact Top2 set 2/4).

## Corrections made in v1.0.1

1. Current ACS Figure 1–6 source-data mapping replaced stale earlier-draft figure numbering.
2. Historical absolute server paths in one stored reference output were sanitized.
3. The README now makes the pinned Python 3.12.3 environment explicit.
4. Exact-environment and cross-version numerical portability are separated rather than treating ~1e-8 floating-point drift as a scientific failure.
5. A stale reference to a missing clean-room report was removed.
6. Software/data/literature-metadata license boundaries were made explicit.
7. Author-verified `n=3` and `(G4S)3` experimental reporting annotations were added without synthesizing raw replicate values, SDs or missing phenotypes.

## Independent portability check

On Python 3.13.5 / NumPy 2.3.5 / SciPy 1.17.0, the unmodified scientific pipeline reproduced all manuscript metrics and all 96 M0–M3 pair directions; maximum probability deviation was 1.0754e-08. This is documented as compatibility-level floating-point variation.

## Remaining placeholders

Only public-deposit metadata remain: GitHub URL, software Zenodo DOI, article DOI and (if deposited separately) trajectory DOI. These do not alter the frozen science.
