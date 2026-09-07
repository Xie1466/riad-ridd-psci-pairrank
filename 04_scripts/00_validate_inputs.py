#!/usr/bin/env python3
import argparse, hashlib, json
from pathlib import Path
from _bootstrap import RELEASE, default_data

def main():
    p=argparse.ArgumentParser(description="Fail-closed validation of bundled frozen inputs.")
    p.add_argument("--release-root",type=Path,default=RELEASE); p.add_argument("--data-root",type=Path,default=default_data()); a=p.parse_args()
    required={"internal16/INTERNAL16_PSCI_FINAL_VALUES.csv":17,"internal16/pair_definitions_24.csv":25,"literature_prior/construct_architecture_91.csv":92,"literature_prior/validated_pairs_381.csv":382,"literature_prior/eligible_non_tie_pairs_377.csv":378,"model/frozen_all_model_predictions.csv":145}
    checks=[]
    for rel, lines in required.items():
        path=a.data_root/rel
        if not path.is_file(): raise SystemExit(f"MISSING_REQUIRED_INPUT: {path}")
        observed=sum(1 for _ in path.open(encoding="utf-8"))
        if observed!=lines: raise SystemExit(f"ROW_COUNT_MISMATCH: {path}: {observed} != {lines}")
        checks.append({"path":str(path.relative_to(a.release_root)),"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"lines":observed})
    print(json.dumps({"status":"PASS","inputs":checks},indent=2))
if __name__=="__main__": main()
