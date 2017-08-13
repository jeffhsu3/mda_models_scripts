"""
"""
import unittest
import numpy as np
import pandas as pd

from figs.worm import Worms
from figs.survival import kill_adults

np.random.seed(10)

from IPython import embed

class Test_Kill_Adult_Worms(unittest.TestCase):
    """
    """
    def setUp(self):
        villages = [0, 0, 0, 0, 1]
        sex = ['M', 'M', 'F', 'F', 'F']
        hostidx = ['v0h1', 'v0h1', 'v0h1', 'v0h1', 'v0h2']
        stage = ['A', 'A', 'A', 'A', 'A', 'A']
        R0net = [0.66, 0.5299222, 0.658231, 0.444, 0.222]
        fec = [0, 0, 10, 2, 20]
        age = [7, 8, 6, 11, 9, 12]
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
            [0, 0, 0, 0],
            [4, 4, 4, 4], 
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
            'age': age,
            'stage': stage,
            'R0net': R0net})
        worms = Worms(meta, haplotype1={'0' : loc0, '1' : hap1},
            haplotype2={'1': hap2}, positions=positions)
        self.worms = worms

        hostidx = ['v0h1', 'v0h2']
        village = [0, 1]
        agedeath = [50, 67]
        age = [50, 65]

        dfHost = pd.DataFrame({'hostidx': hostidx, 'village' : village, 
            'age': age, 'agedeath': agedeath})
        self.hosts = dfHost

    def test_kill_adult_worms(self):
        newhost = kill_adults(self.worms, self.hosts, 12, 2, 2, 2)
        embed()

    def test_mda_survival(self):
        pass



if __name__ == '__main__':
    unittest.main()

