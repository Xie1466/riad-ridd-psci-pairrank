# Code Availability

Custom code used to assemble PSCI descriptors from deposited descriptor inputs, fit the literature-derived architecture prior, perform whole-enzyme-system leave-one-system-out PSCI-PairRank evaluation, reproduce the reported computational metrics, and generate computational figures is available at https://github.com/Xie1466/riad-ridd-psci-pairrank (release `v1.0.2-acs-reproducibility`). The versioned Zenodo archive is linked from the GitHub release; the stable software concept DOI is https://doi.org/10.5281/zenodo.22638943.

The release includes the frozen computational environment, exact run instructions, model configuration, deterministic seed policy, reference outputs, automated tests, and SHA-256 manifests. Under Python 3.12.3, NumPy 1.26.4 and SciPy 1.11.4, frozen probabilities reproduce at the strict `1e-12` gate. Compatible modern environments use `1e-7` while still requiring identical pair directions and manuscript-level metrics.
