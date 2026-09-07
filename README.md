# RIAD/RIDD PSCI-PairRank ACS reproducibility release v1.0.1

This is the compact code-and-data release for the computational component of the RIAD/RIDD protein-scaffold study **“Productive spatial coupling shapes topology-dependent performance of modular enzyme scaffolds.”** It contains Literature V3.2 (19 studies, 20 systems, 91 constructs, 102 phenotype records and 381 validated pairs), the Internal16 four-system/16-topology freeze, the structured PSCI representation, the frozen literature architecture prior and the final M2 PSCI-PairRank evaluation.

Repository: https://github.com/Xie1466/riad-ridd-psci-pairrank

## What this release reproduces

The quick workflow starts from frozen trajectory-derived PSCI values and public literature-architecture tables and re-executes:

1. PSCI feature-table assembly/validation;
2. the literature architecture prior;
3. M0–M3 whole-enzyme-system LOSO evaluation;
4. construct-level mean-win-probability ranking;
5. the reported computational metrics; and
6. compact computational check figures.

Expected M2 results are: **19/24 pair directions (79.2%)**, Spearman **0.654846**, Kendall tau-b **0.521749**, Top1 **3/4**, true-best-in-predicted-Top2 **3/4**, and exact Top2 set **2/4**.

## Recommended exact reproduction

Use the pinned environment whenever possible. With Conda/Mamba:

```bash
conda env create -f 01_environment/environment.yml
conda activate riad-ridd-psci-reproduction
bash 04_scripts/run_quick_reproduction.sh
PYTHONPATH=07_tests python -m unittest discover -s 07_tests -p 'test_*.py'
```

A `venv` is also acceptable **only when the interpreter is Python 3.12.3**:

```bash
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r 01_environment/requirements-lock.txt
bash 04_scripts/run_quick_reproduction.sh
PYTHONPATH=07_tests python -m unittest discover -s 07_tests -p 'test_*.py'
```

`01_environment/check_runtime.py` reports whether the exact pinned runtime is active. The pinned environment uses the strict 1e-12 probability gate. Compatible modern numerical stacks use a 1e-7 probability gate while still requiring identical pair directions and manuscript metrics; see `08_docs/PORTABILITY_AUDIT.md`.

## Current ACS Figure 1–6 source-data mapping

`05_data/figures/` is synchronized to the current ACS manuscript:

- Fig. 1 — Literature V3.2 scale + controlled reciprocal topology grammar.
- Fig. 2 — experimental reciprocal-topology atlas.
- Fig. 3 — proximity/contact paradox and geometric-extreme audit.
- Fig. 4 — PSCI feature landscape and MD QC.
- Fig. 5 — PSCI-PairRank probabilities, held-out system metrics and M0–M3 challenge.
- Fig. 6 — scaffoldability regimes and experimental-prioritization workflow.

## Experimental reporting annotations

The compact release does **not** synthesize missing replicate values or SDs. `05_data/internal16/author_verified_reporting_annotations.csv` records the manuscript-level author verification that the experimental series used three independent replicates (`n = 3`) and the common flexible `(G4S)3` experimental linker. Individual replicate measurements remain absent from this compact computational archive unless separately deposited from authoritative laboratory records.

## Full trajectory reproduction

Production trajectories are intentionally not included in this small archive. `04_scripts/run_full_reproduction.sh` fails closed until a separately archived standardized-trajectory DOI/mount is supplied. The required 16-trajectory/reference/topology/index inventory is described under `05_data/md_reduced_or_manifest/`. This compact release therefore supports exact/portable **model-level** reproduction from the frozen PSCI matrix, not a claim of trajectory-to-feature recreation.

## Licensing and provenance

Software is MIT licensed. Project-authored derived data are CC BY 4.0; literature metadata retain source-specific terms (`DATA_LICENSE.md`). Frozen repository commits and per-file provenance are under `09_provenance/`. The release contains no article PDFs or copied third-party figures.

No article, software, or trajectory DOI is claimed in this pre-publication staging release. The software DOI will be added only after Zenodo issues a validated identifier; an article DOI can be linked after publication.
