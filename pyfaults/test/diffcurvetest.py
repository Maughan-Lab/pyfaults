import unittest
from pyfaults.diffCalc import diffCurve
from pyfaults import importFile

class TestDiffCurve(unittest.TestCase):
    exptQ, exptInts = importFile('./testdata/', 'expt_XRD')
    simQ, simInts = importFile('./testdata/', 'sim_XRD')
    
    diffQ, diffInts = diffCurve(exptQ, simQ, exptInts, simInts)
    
    def test_diffCurve(self):
        self.assertIsNotNone(self.diffQ)
        self.assertIsNotNone(self.diffInts)
        
if __name__=='__main__':
    unittest.main()
        