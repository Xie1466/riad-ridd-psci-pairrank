#!/usr/bin/env python3
import argparse
from pathlib import Path
from _bootstrap import default_data, default_work
from pairrank.pipeline import rankings_and_metrics
def main():
 p=argparse.ArgumentParser(description="Reproduce M0-M3 rankings and manuscript metrics.");p.add_argument("--data-root",type=Path,default=default_data());p.add_argument("--output-dir",type=Path,default=default_work());p.add_argument("--predictions",type=Path);a=p.parse_args(); pred=a.predictions or a.output_dir/"pairwise_probabilities.csv";rows=rankings_and_metrics(a.data_root,pred,a.output_dir);print(f"PASS metrics: {len(rows)} models")
if __name__=="__main__":main()
