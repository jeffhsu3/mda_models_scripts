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
        villages = [0, 0, 0, 0, 1, 0, 0]
        sex = ['M', 'M', 'F', 'F', 'F', 'M', 'F']
        stage = ['A', 'A', 'A', 'A', 'A', 'J', 'J']
        age = [9, 10, 11, 12, 13, 12, 0]
        hostidx = ['v0h1', 'v0h1', 'v0h1', 'v0h1', 'v0h2',
                'v0h3', 'v0h4']
        R0net = [0.66, 0.5299222, 0.658231, 0.444, 0.222, 0.111, 0.23]
        fec = [0, 0, 10, 2, 20, 23, 40]
        positions = {
                '0' : np.array([20, 30], dtype=np.uint64),
                '1' : np.array([1, 10, 50, 100], dtype=np.int64)
                }

        new_positions = []
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
        ng_h1 = np.array([
            [7, 7, 8, 7, 0],
            [2, 3, 4, 5, 6]])
        ng_h2 = np.array([
            [9, 9, 9, 9, 9],
            [2, 3, 4, 5, 6]])
        meta = pd.DataFrame({
            'village': villages, 
            'stage': stage,
            'sex': sex, 
            'hostidx': hostidx, 
            'fec': fec,
            'age': age,
            'R0net': R0net})
        worms = Worms(meta, haplotype1={'0' : loc0, '1' : hap1},
            haplotype2={'1': hap2}, positions=positions)
        new_index = [np.nan, np.nan, np.nan, np.nan, np.nan, 0, 0]
        worms.meta['ng_index'] = new_index
        worms.ng_h1.append({'0': np.array([[0, 1], [3,3]]), '1': ng_h1})
        worms.ng_h2.append({'1': ng_h2})
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


    def test_add_worms_reg(self):
        #### Testing regular add
        self.worms.add_worms(self.worms2, index = [0,1])
        self.assertEqual(self.worms.pos['1'].shape[0] , self.worms.h1['1'].shape[1]) 

        #### WORK IN PROGRESS
        '''
        np.testing.assert_equal(self.worms.h1['1'][:,1],
                np.array([0, 0, 0, 0 , 0, 3, 3], dtype=np.uint8))
        ######
        '''

    def test_drop_worms_regular(self):
        self.worms2.drop_worms([0, 1])
        np.testing.assert_equal(self.worms2.h1['1'].shape[0], 3)
        #np.testing.assert_equal(self.worms2.meta.shape[0], 3)

    def test_drop_worms_mfs(self):
        self.worms.drop_worms([0, 6])
        np.testing.assert_equal(self.worms.ng_h1[0]['1'].shape[0], 1)

    def test_age_worms(self):
        # Kills 1 juvenile
        print(self.worms.meta)
        self.worms.age_worms(0.70, 1, 1)
        surviving_juvenile = self.worms.meta.query('stage == "J"')
        np.testing.assert_equal(surviving_juvenile.shape[0], 0)




if __name__ == '__main__':
    unittest.main()

