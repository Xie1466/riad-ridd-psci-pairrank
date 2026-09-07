#!/usr/bin/env python3
import argparse,csv,math
from pathlib import Path
from _bootstrap import default_data,default_work

def read(p): return list(csv.DictReader(Path(p).open(encoding="utf-8")))
def main():
 p=argparse.ArgumentParser(description="Fail closed unless probabilities and manuscript metrics match frozen references.");p.add_argument("--data-root",type=Path,default=default_data());p.add_argument("--output-dir",type=Path,default=default_work());p.add_argument("--tolerance",type=float,default=1e-12);a=p.parse_args()
 got=read(a.output_dir/"pairwise_probabilities.csv");ref=read(a.data_root/"model/frozen_all_model_predictions.csv"); key=lambda r:(r["model"],r["pair_id"]);g={key(r):r for r in got if key(r)[0] in {"M0_GLOBAL_PSCI_PAIRWISE","M1_LITERATURE_PRIOR","M2_PRIOR_PLUS_PSCI","M3_SYSTEM_CONDITIONED_PSCI"}};f={key(r):r for r in ref if key(r)[0] in {"M0_GLOBAL_PSCI_PAIRWISE","M1_LITERATURE_PRIOR","M2_PRIOR_PLUS_PSCI","M3_SYSTEM_CONDITIONED_PSCI"}}
 if set(g)!=set(f):raise SystemExit("PAIR_KEY_MISMATCH")
 errors=[abs(float(g[k]["probability_i_gt_j"])-float(f[k]["probability_i_gt_j"])) for k in g]
 direction_mismatch=sum((float(g[k]["probability_i_gt_j"])>0.5)!=(float(f[k]["probability_i_gt_j"])>0.5) for k in g)
 if direction_mismatch:raise SystemExit(f"PAIR_DIRECTION_MISMATCH N={direction_mismatch}")
 if max(errors)>a.tolerance:raise SystemExit(f"PROBABILITY_MISMATCH max={max(errors)} tolerance={a.tolerance}")
 metrics={r["model"]:r for r in read(a.output_dir/"metrics_summary.csv")};m2=metrics["M2_PRIOR_PLUS_PSCI"]
 expected={"pair_accuracy":19/24,"rho":0.6548461875980991,"tau":0.521749194749951,"Top1_exact":.75,"true_best_in_predicted_top2":.75,"Top2_exact_set":.5}
 for name,value in expected.items():
  if not math.isclose(float(m2[name]),value,abs_tol=a.tolerance):raise SystemExit(f"M2_METRIC_MISMATCH {name}")
 reported={r["model"]:float(r["pair_accuracy"]) for r in read(a.data_root/"model/reported_model_comparison.csv")};
 for model in reported:
  if model in metrics and not math.isclose(float(metrics[model]["pair_accuracy"]),reported[model],abs_tol=a.tolerance):raise SystemExit(f"MODEL_ACCURACY_MISMATCH {model}")
 print(f"PASS probabilities=96/96 directions_identical=96/96 M2=24/24 max_abs_error={max(errors):.3e} tolerance={a.tolerance:.1e}")
if __name__=="__main__":main()
