import csv, hashlib, json, os, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"reproduction_output"
def read(path):
 with Path(path).open(encoding="utf-8", newline="") as handle:
  return list(csv.DictReader(handle))
def ensure_run():
 if not (OUT/"metrics_summary.csv").is_file(): subprocess.run(["bash",str(ROOT/"04_scripts/run_quick_reproduction.sh")],check=True,cwd=ROOT)
