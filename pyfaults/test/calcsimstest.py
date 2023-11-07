import unittest
from pyfaults.calcSims import calcSims
from pyfaults import importFile

class TestCalcSims(unittest.TestCase):
    
    calcSims('./testoutput/supercell_CIFs/', 'model_info', 0.459744, 30, 0.01, './testoutput/')
    
    dataPath = './testdata/'
    UF_q, UF_ints = importFile(dataPath, 'Unfaulted_sim')
    S1_P10_q, S1_P10_ints = importFile(dataPath, 'S1_P10_sim')
    S2_P10_q, S2_P10_ints = importFile(dataPath, 'S2_P10_sim')
    S1_P20_q, S1_P20_ints = importFile(dataPath, 'S1_P20_sim')
    S2_P20_q, S2_P20_ints = importFile(dataPath, 'S2_P20_sim')
    
    def test_calcSims(self):
        UF_q_sim, UF_ints_sim = importFile('./testoutput/sims/', 'Unfaulted_sim')
        S1_P10_q_sim, S1_P10_ints_sim = importFile('./testoutput/sims/', 'S1_P10_sim')
        S2_P10_q_sim, S2_P10_ints_sim = importFile('./testoutput/sims/', 'S2_P10_sim')
        S1_P20_q_sim, S1_P20_ints_sim = importFile('./testoutput/sims/', 'S1_P20_sim')
        S2_P20_q_sim, S2_P20_ints_sim = importFile('./testoutput/sims/', 'S2_P20_sim')
        
        for i in range(len(self.UF_q)):
            self.assertEqual(self.UF_q[i], UF_q_sim[i])
            self.assertEqual(self.S1_P10_q[i], S1_P10_q_sim[i])
            self.assertEqual(self.S2_P10_q[i], S2_P10_q_sim[i])
            self.assertEqual(self.S1_P20_q[i], S1_P20_q_sim[i])
            self.assertEqual(self.S2_P20_q[i], S2_P20_q_sim[i])
            
        for i in range(len(self.UF_ints)):
            self.assertEqual(self.UF_ints[i], UF_ints_sim[i])
            self.assertEqual(self.S1_P10_ints[i], S1_P10_ints_sim[i])
            self.assertEqual(self.S2_P10_ints[i], S2_P10_ints_sim[i])
            self.assertEqual(self.S1_P20_ints[i], S1_P20_ints_sim[i])
            self.assertEqual(self.S2_P20_ints[i], S2_P20_ints_sim[i])

if __name__=='__main__':
    unittest.main()