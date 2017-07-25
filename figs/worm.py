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
        '''
        if ng_index:
            self.meta['ng_index'] = ng_index
        else:
            self.meta['ng_index'] = np.repeat(np.nan, self.meta.shape[0])
        '''
        self.ng_h1 = []
        self.ng_h2 = []
        self.new_positions = []
        self.new_pos_iix = []

    def _merge_positions(self, loc, oworm, index):
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
        


    def add_worms(self, oworms, index, update=False):
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
        # Reimplantations
        # Alpha sum
        index = np.array(index)
        if len(index) != 0 and self.meta.shape[0] !=0:
            adults = self.meta.ng_index.isnull().sum()
            shift = self.meta.ng_index.value_counts().cumsum() + adults
            ad_ix = index[index <= adults]
            for i in self.h1.keys():
                self.h1[i] = ndelete(self.h1[i], ad_ix, axis=0)
            for i in self.h2.keys():
                self.h2[i] = ndelete(self.h2[i], ad_ix, axis=0)
            for j, k in enumerate(len(self.ng_h1)):

                pass
            for j, k in enumerate(len(self.ng_h2)):
                pass
            self.meta.drop(index, inplace=True)
            self.meta.reset_index(drop=True, inplace=True)
        




    def _kill_mf(self, shapeMF, scaleMF, increment_age=True):
        mfiix = self.meta[self.meta.stage == "M"].index.values
        kill_mfrand = np.random.random(mfiix.shape[0])
        try:
            kill_mffxage = weibull_min.cdf(self.meta.ix[mfiix].age,
                    shapeMF, loc=0, scale=scaleMF)
        except TypeError:
            kill_mffxage = weibull_min.cdf(0, shapeMF, loc=0, scale=scaleMF)
        dieMF = mfiix[np.where(kill_mfrand < kill_mffxage)]
        if increment_age:
            self.meta.ix[mfiix, 'age'] += 1
        else: pass
        return(dieMF)


    def _kill_juvenile(self, surv_Juv, increment_age=False):
        ''' Returns bool of worms to kill. Also increments age
        '''
        juviix = self.meta[self.meta.stage == "J"].index.values
        kill_juvrand = np.random.random(juviix.shape[0])
        dieJuv = juviix[np.where(kill_juvrand > surv_Juv)]
        if increment_age:
            self.meta.loc[juviix, 'age'] += 1
        else: pass
        # Dataframe mapppning meta to each list
        return(dieJuv)

    def age_worms(self, surv_Juv, shapeMF, scaleMF):
        dieJuv = self._kill_juvenile(surv_Juv, increment_age=True)
        dieMF = self._kill_mf(shapeMF, scaleMF, increment_age=True)
        self.drop_worms(np.append(dieJuv, dieMF))
        #age_juvenile(i)
        #add_only_variants(i, self)
