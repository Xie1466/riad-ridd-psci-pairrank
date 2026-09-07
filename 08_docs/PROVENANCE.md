# Provenance

The code lineage uses five frozen repository states without rewriting their history:

| Repository role | Frozen commit | Current HEAD at release extraction | Relationship |
|---|---|---|---|
| Historical MD/source lineage | `6fef6ba77965e504837e39cf2f5ef0205d86f567` | same | immutable source lineage; worktree contains unrelated untracked later material |
| TrueFit PSCI resolution | `c3fb6030860fa7d19bc0bb3782aa371487f24a93` | same | final structured PSCI and Internal16 tables |
| Pairwise V2/release source | `b895eb013f87f87da3f1e9a3e175bca77f2b30dd` | same | current frozen PairRank plus provenance patch |
| Literature V3.2 C-line | `0e818da437056b0488a4033eac696beb91ada0ea` | same | final literature curation |
| B-line scientific acceptance | `f2f215ed4fa8b5ee86c36b2ec0c8fdbc8c965279` | same | immutable import and final acceptance audit; worktree has unrelated untracked handoff material |

The earlier manuscript bundle was generated at Pairwise commit `b9ad43dfd7930824ecd0e2d452804d41e3bc4b23`; the later `b895eb0` commit adds the provenance/robustness patch without changing frozen M2 results. The current release therefore uses `b895eb0` while retaining bundle tables generated from unchanged scientific inputs.

Exact per-file lineage and hashes are in `09_provenance/SOURCE_FILE_PROVENANCE.csv` and `MANIFEST.sha256`.
