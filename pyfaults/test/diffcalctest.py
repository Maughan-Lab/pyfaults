import unittest
import pyfaults.diffCalc as dc
from pyfaults import importSim

class TestDiffCalc(unittest.TestCase):
    exptQ, exptInts = importSim('./testdata/', 'expt_XRD')
    simQ, simInts = importSim('./testdata/', 'sim_XRD')
    
    diffQ, diffInts = dc.diffCurve(exptQ, simQ, exptInts, simInts)
    
    def test_diffCurve(self):
        self.assertIsNotNone(self.diffQ)
        self.assertIsNotNone(self.diffInts)
        
if __name__=='__main__':
    unittest.main()
        