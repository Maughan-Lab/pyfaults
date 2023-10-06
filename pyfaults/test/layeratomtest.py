import unittest
from pyfaults.layerAtom import LayerAtom
from pyfaults.lattice import Lattice


class TestLayerAtom(unittest.TestCase):
    
    latt = Lattice(1.0, 1.0, 1.0, 90, 90, 90)
    a = LayerAtom('A', 'H1', 'H+', [0,0,0], 1.0, latt)
    
    def test__init__(self):
        self.assertIsNotNone(self.a.layerName)
        self.assertIsNotNone(self.a.atomLabel)
        self.assertIsNotNone(self.a.element)
        self.assertIsNotNone(self.a.xyz)
        self.assertIsNotNone(self.a.occupancy)
        self.assertIsNotNone(self.a.lattice)
        
    def test_setParam(self):
        self.assertEqual('A', self.a.layerName)
        self.assertEqual('H1_A', self.a.atomLabel)
        self.assertEqual('H+', self.a.element)
        self.assertEqual([0,0,0], self.a.xyz)
        self.assertEqual(1.0, self.a.occupancy)
        self.assertEqual(self.latt, self.a.lattice)
        
if __name__=='__main__':
    unittest.main()