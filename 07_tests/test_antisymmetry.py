import unittest
from _util import OUT,read,ensure_run
class TestAntisymmetry(unittest.TestCase):
 def test_probability_complements(self):
  ensure_run();self.assertTrue(all(abs(float(r["probability_i_gt_j"])+float(r["probability_j_gt_i"])-1)<1e-15 for r in read(OUT/"pairwise_probabilities.csv")))
if __name__=="__main__":unittest.main()
