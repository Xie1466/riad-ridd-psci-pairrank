import unittest
from _util import ROOT
class TestOrder(unittest.TestCase):
 def test_order(self):
  text=(ROOT/"03_code/pairrank/pipeline.py").read_text();self.assertIn('FEATURES = ("D", "C", "F1", "F4", "F5", "L")',text);self.assertIn('np.column_stack([train_prior, train_dx])',text)
if __name__=="__main__":unittest.main()
