import unittest

import numpy as np
import pandas as pd

from figs.worm import Worms
from figs.mutation import mutation_fx, mutation_at_segsite

from IPython import embed

np.random.seed(20)

class Test_Mutation(unittest.TestCase):
    def setUp(self):
        villages = [0, 0, 0]
        sex = ['M', 'F', 'F']
        hostidx = ['v0h1', 'v0h1', 'v0h1']
        R0net = [0.5299222, 0.658231, 0.444]
        fec = [0, 10, 2]
        positions = np.array([1, 10, 50, 100],
                dtype=np.uint64)
        hap1 = np.array([[0, 1, 0, 1], 
                         [0, 0, 1, 0], 
                         [1, 0, 0, 0]], 
                         dtype = np.uint8)
        hap2 = np.array([[0, 0, 0, 0], 
                         [0, 1, 0, 1], 
                         [1, 1, 0, 0]], 
                         dtype = np.uint8)
        meta = pd.DataFrame({
            'village': villages, 
            'sex': sex, 
            'hostidx': hostidx, 
            'fec': fec})
        hname = str(0)
        worms = Worms(meta, haplotype1 = {hname : hap1},
            haplotype2 = {hname : hap2}, positions={hname : positions})
        self.worms = worms


    def test_one_mutation(self):
        # New position should be 138 with seed
        # Mutation is in first worm and on hap1
        orig_shape = self.worms.h1['0'].shape[1]
        worms, newpositions = mutation_fx(1, self.worms, [0.001],
            [1e-9], [200])
        self.assertEqual(len(newpositions), 
                worms.h1['0'].shape[1] - orig_shape) 
        np.testing.assert_equal(worms.h1['0'][:, -1],
                np.array([1, 0, 0], dtype=np.uint8))


    def test_multiple_mutations(self):
        # Mainly test speed and insertion orders
        orig_nsegsites = self.worms.h1['0'].shape[1]
        worms, newpositions = mutation_fx(1, self.worms, [0.005],
            [1e-9], [200])
        print(orig_nsegsites)


    def test_mutation_at_segsites(self):
        mutation_at_segsite(10, '0', self.worms, 2)
        # Whap is 0
        np.testing.assert_equal(self.worms.h1['0'][2, 1], 1)


if __name__ == '__main__':
    unittest.main()
