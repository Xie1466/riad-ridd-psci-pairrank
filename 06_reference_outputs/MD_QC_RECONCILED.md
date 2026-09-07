# MD QC reconciled

This is a factual reconciliation only. No trajectory classification was changed. The final bundle table and both upstream records (`metadata/analysis_trajectory_manifest.csv` and `metadata/basic_md_summary.csv`) agree for all 16 topologies.

| Classification | Count | Topology IDs |
|---|---:|---|
| `not_converged_within_100ns` | 3 | DPA_T01, GPTRC_T01, RHLARHLB_T04 |
| `slow_drift` | 12 | DPA_T02, DPA_T03, DPA_T04, GIDAE_T01, GIDAE_T02, GIDAE_T04, GPTRC_T02, GPTRC_T03, GPTRC_T04, RHLARHLB_T01, RHLARHLB_T02, RHLARHLB_T03 |
| `late_transition` | 1 | GIDAE_T03 |
| `equilibrium_confirmed=true` | 0 | none |

## Why the old narrative emphasized DPA_T01

DPA_T01 had a dedicated nojump/whole/center/compact PBC-processing recovery and therefore received a topology-specific resolution note. That note described **technical trajectory recovery**, not convergence. GPTRC_T01 and RHLARHLB_T04 did not require that special processing narrative, but both are independently classified `not_converged_within_100ns`. The earlier emphasis was therefore incomplete editorial emphasis, not a discrepancy in the underlying CSV.

All mechanistic use remains descriptive: technical QC passed for 16/16, none has a confirmed equilibrium window, and three are explicitly not converged within 100 ns.
