import unittest

from pyfaults.simXRD import fullSim

class TestSimXRD(unittest.TestCase):
    
    q, ints = fullSim('./testdata/', 'testunitcell', 1.54, 20)
    
    def test_fullSim(self):
        self.assertIsNotNone(self.q)
        self.assertIsNotNone(self.ints)
        
if __name__=='__main__':
    unittest.main()