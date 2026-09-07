#!/usr/bin/env python3
import argparse, csv, shutil
from pathlib import Path
from _bootstrap import default_data, default_work

def main():
    p=argparse.ArgumentParser(description="Validate/copy frozen PSCI matrix; trajectory mode is fail-closed until a DOI deposit is supplied.")
    p.add_argument("--data-root",type=Path,default=default_data()); p.add_argument("--output-dir",type=Path,default=default_work()); p.add_argument("--trajectory-root",type=Path); a=p.parse_args()
    if a.trajectory_root is not None: raise SystemExit("FULL_TRAJECTORY_REPRODUCTION_NOT_AVAILABLE: no exact trajectory deposit is bundled")
    src=a.data_root/"internal16/INTERNAL16_PSCI_FINAL_VALUES.csv"; rows=list(csv.DictReader(src.open(encoding="utf-8")))
    if len(rows)!=16 or {r["primary_window"] for r in rows}!={"0-75 ns"}: raise SystemExit("PSCI_FREEZE_VALIDATION_FAILED")
    a.output_dir.mkdir(parents=True,exist_ok=True); shutil.copyfile(src,a.output_dir/"psci_features.csv"); print("PASS frozen PSCI matrix: 16 topologies")
if __name__=="__main__": main()
