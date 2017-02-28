"""
""" 
from numpy import delete as ndelete 
from numpy import vstack
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


    def _merge_positions(self, loc, oworm, newpos = None):
        # Not the fastest
        assert self.h1[loc].shape[1] == len(self.pos[loc])
        pos1 = np.copy(self.pos[loc])
        pos2 = np.copy(oworm.pos[loc])
        common = np.intersect1d(self.pos[loc], oworm.pos[loc])
        m1 = [i for i in pos1 if i not in common]
        m2 = [i for i in pos2 if i not in common]
        n1 = self.h1[loc].shape[0]
        for i in m2:
            iix = np.argmax(pos1 > i)
            self.h1[loc] = np.insert(self.h1[loc], iix, 
                    np.zeros(n1, dtype=np.uint8), axis=1)
            try:
                self.h2[loc] = np.insert(self.h2[loc], iix, 
                        np.zeros(n1, dtype=np.uint8), axis=1)
            except KeyError:
                pass
            pos1 = np.insert(pos1, iix, i)
        self.pos[loc] = pos1
        n2 = oworm.h1[loc].shape[0]
        for i in m1:
            iix = np.argmax(pos2 > i)
            oworm.h1[loc] = np.insert(oworm.h1[loc], iix, 
                    np.zeros(n2, dtype=np.uint8), axis=1)
            try:
                oworm.h2[loc] = np.insert(oworm.h2[loc], iix, 
                        np.zeros(n2, dtype=np.uint8), axis=1)
            except KeyError:
                pass
            pos2 = np.insert(pos2, iix, i)
        return(oworm)

            

    def add_worms(self, oworms, index):
        """
        Parameters
        ----------
        oworms : figs.worm.Worms object
            other Worms object to add worms from
        index : int list
            numerical index from the other Worms object to add
        """
        if len(index) != 0 and self.meta.shape[0] !=0:
            self.meta = pd.concat([self.meta, oworms.meta.ix[index, :]], 
                    ignore_index=True)
            self.meta.reset_index(drop=True, inplace=True) 
            for i in oworms.h1.keys():
                if np.array_equal(self.pos[i],  oworms.pos[i]):
                    self.h1[i] = vstack((self.h1[i], oworms.h1[i][index,:]))
                    if i in oworms.h2.keys():
                        self.h2[i] = vstack((self.h2[i],
                            oworms.h2[i][index,:]))
                    else:
                        pass
                else:
                    _oworm = self._merge_positions(i, oworms)
                    self.h1[i] = vstack((self.h1[i], _oworm.h1[i][index,:]))
                    try:
                        self.h2[i] = vstack((self.h2[i],
                            _oworm.h2[i][index,:]))
                    except KeyError:
                        pass
        elif self.meta.shape[0] == 0 and len(index) != 0:
            self.meta = oworms.meta.ix[index, :]
            self.meta.reset_index(drop=True, inplace=True)
            for i in oworms.h1.keys():
                self.h1[i] = oworms.h1[i][index, :]
                self.pos[i] = oworms.pos[i]
            for i in oworms.h2.keys():
                self.h2[i] = oworms.h2[i][index, :]

        else:
            self.meta = pd.concat([self.meta, oworms.meta.ix[index, :]], 
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


    def calc_allele_frequencies(self, host=None, village=None):
        """ Calculate allele frequencies
        """
        all_loci_shape = [i.shape[0] for _, i in self.h1.iteritems()]
        allele_freqs = np.empty(np.sum(all_loci_shape), dtype=np.float64)
        c=0
        nworms = self.meta.shape[0]
        for loc in self.h1.keys():
            if loc in self.h2.keys():
                nsites = self.h1[loc].shape[1]
                allele_freqs[c:c+nsites] = np.sum(self.h1[loc] +\
                        self.h2[loc])/float(2*nworms)
            else:
                allele_freqs[c: c+nsites] = np.sum(self.h1[loc])/float(nworms)
            c += self.h1[loc].shape[1]
        return(allele_freqs)
