import hashlib,json,unittest
from pathlib import Path
from _util import ROOT
class TestManifest(unittest.TestCase):
 def test_manifest(self):
  for line in (ROOT/"MANIFEST.sha256").read_text().splitlines():
   digest,rel=line.split("  ",1);self.assertEqual(hashlib.sha256((ROOT/rel).read_bytes()).hexdigest(),digest)
if __name__=="__main__":unittest.main()
