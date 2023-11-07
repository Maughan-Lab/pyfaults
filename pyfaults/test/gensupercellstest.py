import unittest
from pyfaults.genUnitCell import genUnitCell
from pyfaults.genSupercells import genSupercells

class TestGenSupercells(unittest.TestCase):
    
    lattParams = [1.0, 1.0, 1.0, 90, 90, 90]
    unitcell = genUnitCell('./testdata/', 'testGetLayers', lattParams, ['A', 'B'])
    
    pList = [0.1, 0.2]
    sList = [[0.25, 0, 0], [0, 0.25, 0]]
    supercells = genSupercells(unitcell, 10, 'A', pList, sList, './testoutput/')
    
    def test_genSupercells(self):
        expModel = ['Unfaulted', 'S1_P10', 'S2_P10', 'S1_P20', 'S2_P20']
        for i in range(0, 5):
            self.assertEqual(expModel[i], self.supercells['Model'][i])

        expSx = [0, 0.25, 0, 0.25, 0]
        for i in range(0, 5):
            self.assertEqual(expSx[i], self.supercells['S_x'][i])
        
        expSy = [0, 0, 0.25, 0, 0.25]
        for i in range(0, 5):
            self.assertEqual(expSy[i], self.supercells['S_y'][i])
        
        expSz = [0, 0, 0, 0, 0]
        for i in range(0, 5):
            self.assertEqual(expSz[i], self.supercells['S_z'][i])
        
        expP = [0, 0.1, 0.1, 0.2, 0.2]
        for i in range(0, 5):
            self.assertEqual(expP[i], self.supercells['P'][i])
        
if __name__=='__main__':
    unittest.main()