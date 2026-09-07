import platform, unittest
import numpy as np
from _util import ROOT,OUT,read,ensure_run
class TestReconstruction(unittest.TestCase):
 def _error(self):
  ensure_run();models={"M0_GLOBAL_PSCI_PAIRWISE","M1_LITERATURE_PRIOR","M2_PRIOR_PLUS_PSCI","M3_SYSTEM_CONDITIONED_PSCI"};key=lambda r:(r["model"],r["pair_id"]);a={key(r):r for r in read(OUT/"pairwise_probabilities.csv")};b={key(r):r for r in read(ROOT/"05_data/model/frozen_all_model_predictions.csv") if r["model"] in models};self.assertEqual(set(a),set(b));errs=[abs(float(a[k]["probability_i_gt_j"])-float(b[k]["probability_i_gt_j"])) for k in a];directions=sum((float(a[k]["probability_i_gt_j"])>0.5)!=(float(b[k]["probability_i_gt_j"])>0.5) for k in a);return max(errs),directions
 def test_portable_probability_equivalence(self):
  err,directions=self._error();self.assertEqual(directions,0);self.assertLessEqual(err,1e-7)
 def test_exact_pinned_environment(self):
  exact=(platform.python_version()=="3.12.3" and np.__version__=="1.26.4")
  if not exact:self.skipTest("strict 1e-12 test requires pinned Python 3.12.3 / NumPy 1.26.4 runtime")
  err,directions=self._error();self.assertEqual(directions,0);self.assertLessEqual(err,1e-12)
if __name__=="__main__":unittest.main()
