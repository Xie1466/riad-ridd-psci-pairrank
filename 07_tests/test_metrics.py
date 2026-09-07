import unittest
from _util import OUT,read,ensure_run
class TestMetrics(unittest.TestCase):
 def test_m2(self):
  ensure_run();r=next(x for x in read(OUT/"metrics_summary.csv") if x["model"]=="M2_PRIOR_PLUS_PSCI");expected={"pair_accuracy":19/24,"rho":0.6548461875980991,"tau":0.521749194749951,"Top1_exact":.75,"true_best_in_predicted_top2":.75,"Top2_exact_set":.5};
  for k,v in expected.items():self.assertAlmostEqual(float(r[k]),v,places=12)
if __name__=="__main__":unittest.main()
