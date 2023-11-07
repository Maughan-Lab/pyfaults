import unittest

from pyfaults.supercell import Supercell, genSupercells
from pyfaults.lattice import Lattice
from pyfaults.layerAtom import LayerAtom
from pyfaults.layer import Layer
from pyfaults.unitcell import Unitcell

class TestSupercell(unittest.TestCase):
    
    latt = Lattice(1.0, 1.0, 1.0, 90, 90, 90)
    Sc1 = LayerAtom('A', 'Sc1', 'Sc3+', [0,0,0], 1.0, latt)
    Sc2 = LayerAtom('A', 'Sc2', 'Sc3+', [0.5,0.5,0], 1.0, latt)
    Li1 = LayerAtom('B', 'Li1', 'Li+', [0.5,0,0.5], 1.0, latt)
    Li2 = LayerAtom('B', 'Li2', 'Li+', [0,0.5,0.5], 1.0, latt)
    
    A = Layer([Sc1, Sc2], latt, 'A')
    B = Layer([Li1, Li2], latt, 'B')
    
    unitcell = Unitcell('testcell', [A, B], latt)
    
    supercell1 = Supercell(unitcell, 100)
    supercell2 = Supercell(unitcell, 100, fltLayer='B', stackVec=[0.1, 0.1, 0.1],
                           stackProb=0.2)
    
    df = genSupercells(unitcell, 100, 'B', [0.1, 0.2], 
                       [[0.25, 0, 0], [0, 0.25, 0]], './output')
    
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
        
    def test_show_faults(self):
        self.assertIsNone(self.supercell1.show_faults())
        
if __name__=='__main__':
    unittest.main()