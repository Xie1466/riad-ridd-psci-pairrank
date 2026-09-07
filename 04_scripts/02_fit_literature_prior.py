#!/usr/bin/env python3
import argparse
from pathlib import Path
from _bootstrap import default_data, default_work
from pairrank.pipeline import fit_priors
def main():
 p=argparse.ArgumentParser(description="Reproduce frozen fold-specific literature architecture prior.");p.add_argument("--data-root",type=Path,default=default_data());p.add_argument("--output-dir",type=Path,default=default_work());a=p.parse_args(); rows=fit_priors(a.data_root,a.output_dir);print(f"PASS literature prior: {len(rows)} Internal16 pair logits")
if __name__=="__main__":main()
