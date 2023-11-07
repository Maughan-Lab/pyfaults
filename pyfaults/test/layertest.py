import unittest
from pyfaults.layer import Layer
from pyfaults.layerAtom import LayerAtom
from pyfaults.lattice import Lattice

class TestLayer(unittest.TestCase):
    
    latt = Lattice(1.0, 2.0, 3.0, 80, 90, 100)
    a = LayerAtom('A', 'a', 'H+', [0,0,0], 1.0, latt)
    lyr = Layer([a], latt, 'A')
    lyr2 = lyr.genChildLayer('A2', [0.5, 0.5, 0.5])
    
    def test__init__(self):
        self.assertIsNotNone(self.lyr.atoms)
        self.assertIsNotNone(self.lyr.lattice)
        self.assertIsNotNone(self.lyr.layerName)
        
    def test_setParam(self):
        self.assertEqual([self.a], self.lyr.atoms)
        self.assertEqual(self.latt, self.lyr.lattice)
        self.assertEqual('A', self.lyr.layerName)
    
    def test_display(self):
        self.assertIsNone(self.lyr.display())
        
    def test_genChildLayer(self):
        self.assertIsNot(self.lyr.atoms, self.lyr2.atoms)
        self.assertEqual('A2', self.lyr2.layerName)
        self.assertEqual(self.latt, self.lyr2.lattice)
    
    def test_getLayers(self):
        from pyfaults.layer import getLayers
        from pyfaults import importCSV

        data = importCSV('./testdata/', 'testGetLayers')
        self.assertIsNotNone(getLayers(data, self.latt, ['A', 'B'], 'c'))

if __name__=='__main__':
    unittest.main()