import unittest
from pyfaults.unitcell import Unitcell, genUnitCell
from pyfaults.layer import Layer
from pyfaults.layerAtom import LayerAtom
from pyfaults.lattice import Lattice

class TestUnitcell(unittest.TestCase):
    
    latt = Lattice(1.0, 1.0, 1.0, 90, 90, 90)
    a = LayerAtom('A', 'H1', 'H+', [0,0,0], 1.0, latt)
    lyr1 = Layer([a], latt, 'A')
    lyr2 = lyr1.genChildLayer('A2', [0.5, 0.5, 0.5])
    
    unitcell = Unitcell('testcell', [lyr1, lyr2], latt)
    
    genCell = genUnitCell('testgencell', './testdata/', 'testGetLayers', latt, 
                          ['A', 'B'], 'c', childLayers=[('A2', 'A', [0.1, 0.1, 0.1])])
    
    def test__init__(self):
        self.assertIsNotNone(self.unitcell.name)
        self.assertIsNotNone(self.unitcell.layers)
        self.assertIsNotNone(self.unitcell.lattice)
        
    def test_setParam(self):
        self.assertEqual('testcell', self.unitcell.name)
        self.assertEqual([self.lyr1, self.lyr2], self.unitcell.layers)
        self.assertEqual(self.latt, self.unitcell.lattice)
        
    def test_genUnitCell(self):
        self.assertEqual('testgencell', self.genCell.name)
        self.assertEqual(3, len(self.genCell.layers))
        self.assertEqual(self.latt, self.genCell.lattice)
        
if __name__=='__main__':
    unittest.main()