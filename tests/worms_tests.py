import unittest
import sys
import cProfile
import copy

import numpy as np
import pandas as pd
from pstats import Stats
from IPython import embed

from figs.worm import Worms

np.random.seed(20)
class Test_Worms(unittest.TestCase):
    def setUp(self):
        #:TODO load this from a pickle?
        villages = [0, 0, 0, 0, 1]
        sex = ['M', 'M', 'F', 'F', 'F']
        hostidx = ['v0h1', 'v0h1', 'v0h1', 'v0h1', 'v0h2']
        R0net = [0.66, 0.5299222, 0.658231, 0.444, 0.222]
        fec = [0, 0, 10, 2, 20]
        positions = {
                '0' : np.array([20, 30], dtype=np.uint64),
                '1' : np.array([1, 10, 50, 100], dtype=np.int64)
                }

        loc0 = np.array([
            [0, 1],
            [0, 0],
            [1, 1],
            [1, 0], 
            [2, 2]], dtype=np.uint8)

        hap1 = np.array([
            [0, 3, 0, 0],
            [4, 3, 4, 4], 
            [6, 6, 6, 6],
            [1, 0, 0, 0],
            [2, 2, 2, 2]], dtype = np.uint8)
        hap2 = np.array([
            [1, 1, 0, 1],
            [0, 0, 0, 0], 
            [5, 5, 5, 5], 
            [1, 1, 0, 0], 
            [2, 2, 2, 2]], dtype = np.uint8)

        meta = pd.DataFrame({
            'village': villages, 
            'sex': sex, 
            'hostidx': hostidx, 
            'fec': fec,
            'R0net': R0net})
        worms = Worms(meta, haplotype1={'0' : loc0, '1' : hap1},
            haplotype2={'1': hap2}, positions=positions)

        positions_2 = {
                '0' : np.array([10, 30], dtype=np.uint64),
                '1' : np.array([1, 5, 29, 100], dtype=np.int64)
                }
        worms2 = Worms(meta, 
                haplotype1={'0' : np.copy(loc0), '1' : np.copy(hap1)},
            haplotype2={'1': np.copy(hap2)}, 
            positions=positions_2)
        self.worms = worms
        self.worms2 = worms2 


    def test_add_worms(self):
        #### Testing regular add
        self.worms.add_worms(self.worms2, index = [0,1])
        self.assertEqual(self.worms.pos['1'].shape[0] , self.worms.h1['1'].shape[1]) 
        np.testing.assert_equal(self.worms.h1['1'][:,1],
                np.array([0, 0, 0, 0 , 0, 3, 3], dtype=np.uint8))

    def test_add_worms_empty(self):
        pass

    def test_drop_worms(self):
        self.worms2.drop_worms([0, 1])
        np.testing.assert_equal(self.worms2.h1['1'].shape[0], 3)
        np.testing.assert_equal(self.worms2.meta.shape[0], 3)



if __name__ == '__main__':
    unittest.main()

