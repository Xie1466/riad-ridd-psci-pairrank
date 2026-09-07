#!/usr/bin/env python3
import argparse
from pathlib import Path
from _bootstrap import default_data, default_work
from pairrank.pipeline import run_loso
def main():
 p=argparse.ArgumentParser(description="Reproduce frozen M0-M3 whole-system LOSO probabilities.");p.add_argument("--data-root",type=Path,default=default_data());p.add_argument("--output-dir",type=Path,default=default_work());p.add_argument("--prior-artifact",type=Path);a=p.parse_args(); prior=a.prior_artifact or a.output_dir/"literature_prior_folds.json";rows=run_loso(a.data_root,a.output_dir,prior);print(f"PASS PairRank LOSO: {len(rows)} model-pair rows")
if __name__=="__main__":main()
