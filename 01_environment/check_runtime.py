#!/usr/bin/env python3
import argparse, json, platform, sys
try:
 import numpy as np
except Exception:
 np=None
try:
 import scipy
except Exception:
 scipy=None
try:
 import matplotlib
except Exception:
 matplotlib=None
try:
 import yaml
except Exception:
 yaml=None
expected={"python":"3.12.3","numpy":"1.26.4","scipy":"1.11.4","matplotlib":"3.6.3","yaml":"6.0.1"}
actual={"python":platform.python_version(),"numpy":getattr(np,"__version__",None),"scipy":getattr(scipy,"__version__",None),"matplotlib":getattr(matplotlib,"__version__",None),"yaml":getattr(yaml,"__version__",None)}
exact=all(actual[k]==expected[k] for k in expected)
parser=argparse.ArgumentParser();parser.add_argument("--exact-flag",action="store_true");a=parser.parse_args()
if a.exact_flag: print("1" if exact else "0")
else: print(json.dumps({"exact_pinned_environment":exact,"expected":expected,"actual":actual,"recommended_probability_tolerance":1e-12 if exact else 1e-7},indent=2))
