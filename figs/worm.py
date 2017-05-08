"""
"""
from numpy import delete as ndelete
from numpy import vstack, hstack
import numpy as np
import pandas as pd

class Worms(object):
    def __init__(self, meta, haplotype1=None, haplotype2=None,
            positions = None, selection=None, cds_coords=None):
        self.meta = meta
        if haplotype1:
            self.h1 = haplotype1
        else:
            self.h1= {}
        if haplotype2:
            self.h2 = haplotype2
        else:
            self.h2 = {}
        if positions:
            self.pos = positions
        else:
            self.pos = {}
        if selection:
            self.sel = selection
        else:
            self.sel = {}
        if cds_coords:
            self.coord = cds_coords
        else:
            self.coord = {}

    def _merge_positions(self, loc, oworm):
        assert self.h1[loc].shape[1] == len(self.pos[loc])
        pos1 = np.copy(self.pos[loc])
        pos2 = np.copy(oworm.pos[loc])
        m1 = np.setdiff1d(pos1, pos2)
        m2 = np.setdiff1d(pos2, pos1)
        n1 = self.h1[loc].shape[0]


        self.h1[loc] = hstack((self.h1[loc],
            np.zeroes((n1, len(m2)), dtype=np.uint8)))
        iix = np.searchsorted(self.pos[loc], m2)

        np.arange(
        






    def add_worms(self, oworms, new_pos, index, 
            update=False):
        """
        Parameters
        ----------
        oworms : figs.worm.Worms object
            other Worms object to add worms from
        index : int list
            numerical index from the other Worms object to add
        new_pos : dict of lists
        """
        if len(new_pos) != 0 and self.meta.shape[0] !=0:
            self.meta = pd.concat([self.meta, oworms.meta.ix[index,:]],
                    ignore_index=True)
            self.meta.reset_index(drop=True, inplace=True)
            for i in oworms.h1.keys():
                if np.array_equal(self.pos[i],  oworms.pos[i]):
                    self.h1[i] = vstack((self.h1[i], oworms.h1[i]))
                    if i in oworms.h2.keys():
                        self.h2[i] = vstack((self.h2[i],
                            oworms.h2[i]))
                    else:
                        pass
                else:
                    self.h1[i] = hstack((self.h1[i],
                        np.zeros((self.h1[i].shape[0], len(new_pos[i])), 
                            dtype = np.uint8)))
                    self.pos[i] = np.append(self.pos[i], new_pos[i])
                    iix = np.argsort(self.pos[i])
                    self.h1[i] = self.h1[i][:, iix]
                    self.h1[i] = vstack((self.h1[i], oworms.h1[i]))
                    self.pos[i] = self.pos[i][iix]
                    try:
                        self.h2[i] = hstack((self.h2[i],
                            np.zeros((self.h2[i].shape[0], len(new_pos[i])),
                                dtype=np.uint8)))
                        self.h2[i] = self.h2[i][:, iix]
                        self.h2[i] = vstack((self.h2[i],
                            oworms.h2[i]))
                    except KeyError:
                        pass
                if update:
                    oworms.h1[i] = oworms
                    oworms.pos[i] =  self.pos[i]
        elif self.meta.shape[0] == 0 and len(new_pos) != 0:
            self.meta = oworms.meta
            self.meta.reset_index(drop=True, inplace=True)
            for i in oworms.h1.keys():
                self.h1[i] = oworms.h1[i]
                self.pos[i] = oworms.pos[i]
            for i in oworms.h2.keys():
                self.h2[i] = oworms.h2[i]
        else:
            self.meta = pd.concat([self.meta, oworms.meta],
                    ignore_index=True)
            self.meta.reset_index(drop=True, inplace=True)
            print("Nothing to add")




    def drop_worms(self, index):
        if len(index) != 0 and self.meta.shape[0] != 0:
            self.meta.drop(index, inplace=True)
            self.meta.reset_index(drop=True, inplace=True)
            for i in self.h1.keys():
                self.h1[i] = ndelete(self.h1[i], index, axis=0)
            for i in self.h2.keys():
                self.h2[i] = ndelete(self.h2[i], index, axis=0)
        else:
            print('No worms to drop')
            pass


    def calc_allele_frequencies(self, host=None, village=None, loci=None):
        """ Calculate allele frequencies
        """
        if loci == None:
            loci = self.h1.keys()
        else:
            loci = loci
        all_loci_shape = [self.h1[i].shape[1] for i in loci]
        allele_freqs = np.zeros(np.sum(all_loci_shape), dtype=np.float64)
        c=0
        nworms = self.meta.shape[0]
        for loc in loci:
            nsites = self.h1[loc].shape[1]
            if loc in self.h2.keys():
                allele_freqs[c:c+nsites] = np.sum(self.h1[loc] +\
                        self.h2[loc], dtype=np.float64, axis=0)/2*nworms
            else:
                allele_freqs[c: c+nsites] = np.sum(self.h1[loc],
                        dtype=np.float64, axis=0)/nworms
            c += self.h1[loc].shape[1]
        return(allele_freqs)
