from pathlib import Path
import sys

RELEASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RELEASE / "03_code"))

def default_data(): return RELEASE / "05_data"
def default_work(): return RELEASE / "reproduction_output"
