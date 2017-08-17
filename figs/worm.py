"""
"""
from numpy import delete as ndelete
from numpy import vstack, hstack
import numpy as np
import pandas as pd

from scipy.stats import weibull_min

class Worms(object):
    def __init__(self, meta, haplotype1=None, haplotype2=None,
            positions=None, selection=None, cds_coords=None,
            ng_index=None):
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
        if ng_index:
            self.meta['ng_index'] = ng_index
        else:
            self.meta['ng_index'] = np.repeat(np.nan, self.meta.shape[0])
        self.ng_h1 = []
        self.ng_h2 = []
        self.ng_pos = []
        self.new_pos= []
        self.new_pos_iix = []

    def _merge_positions_2(self, loc, oh1, oh2, pos):
        """
        """
        #assert(len(self.pos[loc]) > oworm.shape[1])
        # Equalize both
        pos1 = np.copy(self.pos[loc])
        pos2 = pos
        m1 = np.setdiff1d(pos1, pos2)
        m2 = np.setdiff1d(pos2, pos1)
        n1 = self.h1[loc].shape[0]
        n2 = oh1.shape[0]
        iix = np.searchsorted(self.pos[loc], m2)
        iixf = [i for i in range(len(pos1) + len(m2)) if i not in iix]
        iixf.extend(iix)
        self.h1[loc] = hstack((self.h1[loc],
            np.zeros((n1, len(m2)), dtype=np.uint8)))
        self.h1[loc] = self.h1[loc][:, iixf]
        try:
            self.h2[loc] = hstack((self.h2[loc],
                np.zeros((n1, len(m2)), dtype=np.uint8)))
            self.h2[loc] = self.h2[loc][:, iixf]
        except KeyError:
            pass
        
        self.pos[loc] = np.concatenate((pos1, m2))
        self.pos[loc] = self.pos[loc][iixf]

        oh1 = hstack((oh1, 
            np.zeros((n2, len(m1)), dtype=np.uint8)))
        iix = np.searchsorted(pos, m1)
        iixf = [i for i in range(len(pos2) + len(m1)) if i not in iix]
        iixf.extend(iix)
        oh1 = oh1[:, iixf]
        import ipdb
        ipdb.set_trace()
        try:
            assert oh1.shape[1] == self.h1[loc].shape[1]
        except AssertionError:
            import ipdb
            ipdb.set_trace()
        self.h1[loc] = vstack((self.h1[loc], oh1))
        try:
            oh2 = hstack((oh2,
                np.zeros((n2, len(m1)), dtype=np.uint8)))
            oh2 = oh2[:, iixf]
            self.h2[loc] = vstack((self.h2[loc], oh2))
            self.h2[loc] = self.h2[loc].copy(order='C')
        except KeyError:
            pass
        except ValueError:
            print(loc)
        self.h1[loc] = self.h1[loc].copy(order='C')

    def _merge_positions(self, loc, oh1, oh2, index):
        assert self.h1[loc].shape[1] == len(self.pos[loc])
        pos1 = np.copy(self.pos[loc])
        pos2 = np.copy(oworm.pos[loc])
        m1 = np.setdiff1d(pos1, pos2)
        m2 = np.setdiff1d(pos2, pos1)
        n1 = self.h1[loc].shape[0]
        n2 = oworm.h1[loc].shape[0]
        self.h1[loc] = hstack((self.h1[loc],
            np.zeros((n1, len(m2)), dtype=np.uint8)))
        iix = np.searchsorted(self.pos[loc], m2)
        iixf = [i for i in range(len(pos1) + len(m2)) if i not in iix]
        iixf.extend(iix)
        self.h1[loc] = self.h1[loc][:, iixf]
        try:
            self.h2[loc] = hstack((self.h2[loc],
                np.zeros((n1, len(m2)), dtype=np.uint8)))
            self.h2[loc] = self.h2[loc][:, iixf]
        except KeyError:
            pass
        
        self.pos[loc] = np.concatenate((pos1, m2))
        self.pos[loc] = self.pos[loc][iixf]

        oh1 = oworm.h1[loc].copy()
        oh1 = hstack((oh1, 
            np.zeros((n2, len(m1)), dtype=np.uint8)))
        iix = np.searchsorted(oworm.pos[loc], m1)
        iixf = [i for i in range(len(pos2) + len(m1)) if i not in iix]
        iixf.extend(iix)
        oh1 = oh1[index][:, iixf]
        
        assert oh1.shape[1] == self.h1[loc].shape[1]
        self.h1[loc] = vstack((self.h1[loc], oh1))
        try:
            oh2 = oworm.h2[loc].copy()
            oh2 = hstack((oh2, 
                np.zeros((n2, len(m1)), dtype=np.uint8)))
            oh2 = oh2[index][: , iixf]
            self.h2[loc] = vstack((self.h2[loc], oh2))
            self.h2[loc] = self.h2[loc].copy(order='C')
        except KeyError:
            pass
        self.h1[loc] = self.h1[loc].copy(order='C')


    def _prune_genotypes(self, ng_index):
        ipos = self.new_pos[int(ng_index)]
        i = self.new_pos_iix[int(ng_index)]
        for loc in i.keys():
            iix = i[loc]
            to_remove = np.all(self.ng_h1[ng_index][loc][:, iix]==0, axis=0)
            to_delete = iix[to_remove]
            ndelete(self.ng_h1[ng_index][loc], to_delete, axis=0)
            i[loc] = [x for (x, v) in zip(iix, to_remove) if v]
            ipos[loc] = [x for (x, v) in zip(ipos[loc], to_remove) if v]
        self.new_pos_iix[int(ng_index)] = i
        self.new_pos[int(ng_index)] = ipos


    def _add_worms(self, ng_index, wormsindex):
        """
        """
        for loc in self.h1.keys():
            opos = self.ng_pos[ng_index][loc]
            #new_pos = self.new_pos[ng_index][loc]
            h1 = self.ng_h1[ng_index][loc][wormsindex, :]
            try:
                h2 = self.ng_h2[ng_index][loc][wormsindex, :]
            except KeyError:
                h2 = None
            self._merge_positions_2(loc, h1, h2, opos)
            self.ng_h1[ng_index][loc] = ndelete(self.ng_h1[ng_index][loc], wormsindex, 0)
            try:
                self.ng_h2[ng_index][loc] = ndelete(self.ng_h2[ng_index][loc], wormsindex, 0)
            except KeyError:
                pass




    def add_worms(self, oworms=None, index=None, update=False):
        """
        Add worms from one worms object to another object
        
        Parameters
        ----------
        oworms : figs.worm.Worms object
            other Worms object to add worms from
        index : int list
            numerical index from the other Worms object to add
        new_pos : dict of lists
        """
        if oworms.meta.shape[0] != 0 and self.meta.shape[0] !=0:
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
                    self._merge_positions(i, oworms, index)
        elif self.meta.shape[0] == 0 and oworms.meta.shape[0] != 0:
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
        """
        if len(index) != 0 and self.meta.shape[0] != 0:
            self.meta.drop(index, inplace=True)
            self.meta.reset_index(drop=True, inplace=True)
            for i in self.h1.keys():
                self.h1[i] = ndelete(self.h1[i], index, axis=0)
            for i in self.h2.keys():
                self.h2[i] = ndelete(self.h2[i], index, axis=0)
            for i in self.ng_h1:
                for j in i.keys():
                    i[j] = ndelete(i[j], index, axis=0)
            for i in self.ng_h2:
                for j in i.ng_h2.keys():
                    i[j] = ndelete(i[j], index, axis=0)
        else:
            print('No worms to drop')
            pass
        """
        # Need to rework this. Or simply move adults up
        index = np.array(index)
        if len(index) != 0 and self.meta.shape[0] !=0:
            adults = self.meta.ng_index.isnull().sum()
            shift = self.meta.ng_index.value_counts().values
            shift = np.insert(shift, 0, adults)
            shift = np.cumsum(shift)
            ad_ix = index[index <= adults]
            for i in self.h1.keys():
                self.h1[i] = ndelete(self.h1[i], ad_ix, axis=0)
            for i in self.h2.keys():
                self.h2[i] = ndelete(self.h2[i], ad_ix, axis=0)
            for j, k in enumerate(self.ng_h1):
                mf_ix = index[np.logical_and(index > shift[j], 
                    index < shift[j+1])] - shift[j] 
                for i in k.keys():
                    k[i] = ndelete(k[i], mf_ix, axis=0)
            for j, k in enumerate(self.ng_h2):
                mf_ix = index[np.logical_and(index > shift[j], 
                    index <= shift[j+1])] 
                for i in k.keys():
                    k[i] = ndelete(k[i], mf_ix, axis=0)
            self.meta.drop(index, inplace=True)
            self.meta.reset_index(drop=True, inplace=True)


    def age_worms(self):
        # Move juveniles to age l3 to adult age 1
        juviix12 = self.meta.query('(age > 12) & (stage == "J")').index.values
        print('juviix12 stuff')
        print(self.meta.loc[juviix12, 'ng_index'].unique())
        # Run specialized addworms here
        # ADD in juveniles to main haplotypes
        ng_indexes = self.meta.loc[juviix12, 'ng_index'].unique()
        for i in ng_indexes:
            print('Aging worms from ng_index: {0!s}'.format(i))
            print(i)
            i = int(i)
            ng_set = self.meta.ix[np.isclose(self.meta.ng_index , i),:]
            ng_set_index = np.arange(ng_set.shape[0])
            ng_hap_index = ng_set_index[(ng_set.age > 12) &\
                    (ng_set.stage == "J")]
            # Not sure if this part 
            # Prune only the about to be added set
            # self._prune_genotypes(i)
            self._add_worms(i, ng_hap_index)
            # ng_index get's reset to nothing if added to adults
            self.meta.loc[self.meta.ng_index == i, 'ng_index'] = np.nan

        self.meta.loc[juviix12, 'stage'] = "A"
        self.meta.loc[juviix12, 'R0net'] += 1
        self.meta.loc[juviix12, 'age'] = 1
        # Resort so that only adults are at the beginning of meta
        if len(juviix12) > 0:
            self.meta = pd.concat([self.meta.query("stage == 'A'"),
                self.meta.query("stage != 'A'")])
            import ipdb
            ipdb.set_trace()
        else: pass

        


    def _partial_add(self, other_worms, new_positions, new_pos_iix):
        """
        """
        other_worms.meta['ng_index'] = np.repeat(len(self.ng_h1), other_worms.meta.shape[0])
        
        self.ng_h1.append(other_worms.h1)
        self.ng_h2.append(other_worms.h2)
        self.new_pos.append(new_positions)
        self.ng_pos.append(other_worms.pos)
        self.new_pos_iix.append(new_pos_iix)
        self.meta = pd.concat([self.meta, other_worms.meta], ignore_index=True)



