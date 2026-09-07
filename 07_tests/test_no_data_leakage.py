import unittest
from _util import ROOT,read
class TestLeakage(unittest.TestCase):
 def test_whole_system_fold(self):
  rows=read(ROOT/"05_data/internal16/pair_definitions_24.csv")
  for heldout in {r["system"] for r in rows}:
   self.assertEqual(len([r for r in rows if r["system"]==heldout]),6);self.assertEqual(len([r for r in rows if r["system"]!=heldout]),18)
 def test_overlap_audit_present(self):self.assertTrue((ROOT/"05_data/literature_prior/overlap_exclusion_audit.csv").is_file())
if __name__=="__main__":unittest.main()
