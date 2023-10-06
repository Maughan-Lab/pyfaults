import unittest

from pyfaults.supercell import Supercell
from pyfaults.unitcell import genUnitCell
from pyfaults.lattice import Lattice

class TestSupercell(unittest.TestCase):
    
    latt = Lattice(1.0, 1.0, 1.0, 90, 90, 90)
    unitcell = genUnitCell('testgencell', './testdata/', 'testGetLayers', latt, 
                          ['A', 'B'], 'c')
    
    supercell1 = Supercell(unitcell, 100)
    supercell2 = Supercell(unitcell, 100, fltLayer='B', stackVec=[0.1, 0.1, 0.1],
                           stackProb=0.2)
    
    def test__init__(self):
        self.assertIsNotNone(self.supercell1.unitcell)
        self.assertIsNotNone(self.supercell1.nStacks)
        self.assertIsNotNone(self.supercell1.layers)
        
        self.assertIsNotNone(self.supercell2.unitcell)
        self.assertIsNotNone(self.supercell2.nStacks)
        self.assertIsNotNone(self.supercell2.layers)
        self.assertIsNotNone(self.supercell2.fltLayer)
        self.assertIsNotNone(self.supercell2.stackVec)
        self.assertIsNotNone(self.supercell2.stackProb)
        
    def test_setParam(self):
        self.assertEqual(100, self.supercell1.nStacks)
        self.assertEqual(100, self.supercell2.nStacks)
        
    def test_setLayers(self):
        self.assertEqual(self.unitcell, self.supercell1.unitcell)
        self.assertIsNotNone(self.supercell1.layers)
        
        self.assertEqual(self.unitcell, self.supercell2.unitcell)
        self.assertIsNotNone(self.supercell2.layers)
        self.assertEqual('B', self.supercell2.fltLayer)
        self.assertEqual([0.1, 0.1, 0.1], self.supercell2.stackVec)
        self.assertEqual(0.2, self.supercell2.stackProb)
        
if __name__=='__main__':
    unittest.main()