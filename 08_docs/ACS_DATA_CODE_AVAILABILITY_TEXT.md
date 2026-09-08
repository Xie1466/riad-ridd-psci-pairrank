# ACS manuscript availability text

## Code Availability

Custom code used to assemble PSCI descriptors from deposited descriptor inputs, fit the literature-derived architecture prior, perform whole-enzyme-system leave-one-system-out PSCI-PairRank evaluation, reproduce the reported computational metrics, and generate computational figures is available at https://github.com/Xie1466/riad-ridd-psci-pairrank, release `v1.0.2-acs-reproducibility`, commit `70fb08f246f88302c6797b5ebccc9774a409637b`, and is permanently archived at Zenodo: https://doi.org/10.5281/zenodo.22656378. The archive includes the pinned environment, exact run instructions, frozen model configuration, deterministic seed policy, tests, reference outputs, and SHA-256 manifests.

In the pinned Python 3.12.3, NumPy 1.26.4, and SciPy 1.11.4 environment, frozen probabilities reproduce at the strict `1e-12` gate. Compatible modern numerical stacks use a `1e-7` probability tolerance and must retain identical pair directions and manuscript-level metrics.

## Data Availability

The archived release contains Literature V3.2 architecture tables (19 studies, 20 systems, 91 constructs, 102 phenotype records, and 381 validated pairs), the Internal16 construct/rank freeze (four systems, 16 topologies, and 24 pair relations), experimental endpoint summaries, MD-QC summaries, the exact PSCI feature dictionary and matrix, all 24 frozen M2 out-of-fold predictions, M0-M3 comparisons, fold coefficients/scalers, Figure 1-6 source tables, provenance, and manifests. Individual experimental replicate measurements and SDs are included only where recoverable from authoritative primary records and are not reconstructed from means.

`ARTICLE_DOI_STATUS = NOT_YET_AVAILABLE`. `TRAJECTORY_DOI_STATUS = NOT_YET_AVAILABLE`. Full production MD trajectories are not included in the current compact software archive; full trajectory-to-feature reproduction is `NOT_AVAILABLE`.
