import unittest

import sys
import os.path
sys.path.append(
            os.path.abspath(os.path.join(os.path.dirname(__file__),
                os.path.pardir)))

from recombination import recombination_locus, recombination_fx
import random
import numpy as np
import pandas as pd

random.seed(100)

class Test_Recombination_Locus(unittest.TestCase):
    """Testing recombination_locus output
    """
    def setUp(self):
        self.basepairs = 10
        self.h1 = [1, 5, 6]
        self.h2 =  [2, 8]
        self.crossover_pos = 4

    def test_recombination_norm(self):
        h1, h2 = recombination_locus(self.h1,
                self.h2, self.crossover_pos)
        np.testing.assert_equal(h1, [1, 8])
        np.testing.assert_equal(h2, [2, 5, 6])

    def test_end(self):
        """Test when recombination happens at the end
        """
        h1, h2 = recombination_locus(self.h1,
                self.h2, 9)
        np.testing.assert_equal(h1, self.h1)
        np.testing.assert_equal(h2, self.h2)

    def test_start(self):
        h1, h2 = recombination_locus(self.h1,
                self.h2, 0)
        np.testing.assert_equal(h1, self.h2)
        np.testing.assert_equal(h2, self.h1)


    def test_multiple_crossovers(self):
        h1, h2 = recombination_locus(self.h1,
                self.h2, 4)
        h1, h2 = recombination_locus(h1, h2, 7)
        np.testing.assert_equal(self.h1, [1,5,6])





class Test_Recombination_Fx(unittest.TestCase):
    def setUp(self):
        villages = [0, 0, 0]
        sex = ['M', 'F', 'F']
        hostidx = ['v0h1', 'v0h1', 'v0h1']
        R0net = [0.5299222, 0.658231, 0.444]
        fec = [0, 10, 2]

        M_hap1= [1, 3, 9]
        M_hap2= [2, 4, 6]

        F_hap1 = [2,7]
        F_hap2 = [5,9]

        df_adult = pd.DataFrame({
            'village': villages,
            'sex' : sex,
            'hostidx' : hostidx,
            'R0net' : R0net,
            'fec' : fec,
            'locus_0_h1' : [M_hap1, F_hap1, M_hap1],
            'locus_0_h2' : [M_hap2, F_hap2, F_hap1],
            })

        self.adult = df_adult

        

    def test_recombination_fx(self):
        df_adult_mf = recombination_fx(1, self.adult, [0.5], [10])




if __name__ == '__main__':
    unittest.main()
