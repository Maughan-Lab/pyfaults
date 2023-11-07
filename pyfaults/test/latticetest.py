import unittest
from pyfaults.lattice import Lattice

class TestLattice(unittest.TestCase):
    
    latt = Lattice(1.0, 2.0, 3.0, 80, 90, 100)
    
    def test__init__(self):
        self.assertIsNotNone(self.latt.a)
        self.assertIsNotNone(self.latt.b)
        self.assertIsNotNone(self.latt.c)
        self.assertIsNotNone(self.latt.alpha)
        self.assertIsNotNone(self.latt.beta)
        self.assertIsNotNone(self.latt.gamma)
        
    def test_setParam(self):
        self.assertEqual(1.0, self.latt.a)
        self.assertEqual(2.0, self.latt.b)
        self.assertEqual(3.0, self.latt.c)
        self.assertEqual(80, self.latt.alpha)
        self.assertEqual(90, self.latt.beta)
        self.assertEqual(100, self.latt.gamma)
        
    def test_display(self):
        expectedStr = '\n'.join(('a : 1.0', 
                         'b : 2.0', 
                         'c : 3.0',
                         'alpha : 80', 
                         'beta : 90',
                         'gamma : 100'))
        self.assertEqual(expectedStr, self.latt.display())
    
if __name__=='__main__':
    unittest.main()