"""
"""

from numpy import delete as ndelete
from numpy import vstack
import pandas as pd


class Worms(object):
    def __init__(self, meta, haplotype1=None, haplotype2=None,
            positions = None):
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

    def add_worms(self, df, index):
        """
        Parameters
        ----------
        df : figs.worm.Worms object
            other Worms object to add worms from
        index : int list
            numerical index from the other Worms object to add
        """
        self.meta = pd.concat([self.meta, df.meta.ix[index, :]], ignore_index=True)
        self.meta.reset_index(inplace=True)
        for i in df.h1.keys():
            self.h1[i] = vstack(self.h1[i], df.h1[i][index,:])
        for i in df.h2.keys():
            self.h2[i] = vstack(self.h2[i], df.h2[i][index,:])


    def drop_worms(self, index):
        self.meta.drop(index)
        self.meta = self.meta.reset_index(drop=True)
        for i in self.h1.keys():
            self.h1[i] = ndelete(self.h1[i], index, axis=0)
        for i in self.h2.keys():
            self.h2[i] = ndelete(self.h2[i], index, axis=0)
        