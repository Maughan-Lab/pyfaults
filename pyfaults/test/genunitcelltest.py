import unittest
from pyfaults.genUnitCell import genUnitCell

class TestGenUnitCell(unittest.TestCase):

    lattParams = [1.0, 1.0, 1.0, 90, 90, 90]
    unitcell = genUnitCell('./testdata/', 'testGetLayers', lattParams, ['A', 'B'])
    
    def test_genUnitCell(self):
        self.assertIsNotNone(self.unitcell.layers)
        self.assertIsNotNone(self.unitcell.lattice)
        
if __name__=='__main__':
    unittest.main()