#!/usr/bin/env python3
import argparse
from pathlib import Path
from _bootstrap import RELEASE, default_work
import sys
sys.path.insert(0, str(RELEASE / "03_code"))
from figures.render import render_check_figures

def main():
 p=argparse.ArgumentParser(description="Reproduce compact computational metric and M2 probability figures.");p.add_argument("--output-dir",type=Path,default=default_work());a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
 paths=render_check_figures(a.output_dir);print(f"PASS figures: {len(paths)}")
if __name__=="__main__":main()
