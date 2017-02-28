#!/u fecsr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:59:27 2017

@author: scott
"""
import numpy as np
import random
from collections import defaultdict

def mutation_at_segsite(newsite, loc, worms, randmf):
    iix = np.where(worms.pos[str(loc)] == newsite)
    whap = str(np.random.randint(1, 3))
    try:
        hap = getattr(worms, "h" + whap)[str(loc)][:, iix]
    except KeyError:
        hap = getattr(worms, "h1")[str(loc)][:, iix]
        whap = str(1)
    if hap[randmf] == 0:
        hap[randmf] = 1
    else:
        hap[randmf] = 0
    getattr(worms,"h"+ whap)[str(loc)][:, iix]  = hap 


def mutation_fx(locus,
                worms,
                mutation_rate,
                recombination_rate,
                basepairs):
    '''Calculates number of mutations, changes haplotypes

    Parameters
    ---------
    locus: int
         number of loci
    dfAdult_mf : figs.Worms object
          New larval parasites
    mutation_rate : float, list
          mutation rates for each locus
    recombination_rate : float, list
        recombination rates for each locus
    basepairs : int, list
          length of each locus in basepairs
    cds_coordinates : list
        list of coding positions

    Returns
    ------
    dfAdult_mf : pandas df
         df of larval parasites with mutated positions
    mutations : int, list
         list of positions of new mutations
    '''
    new_positions = defaultdict(list)
    nworms = worms.meta.shape[0]
    for loc in range(locus):
        if recombination_rate[loc] == 0:
            mut_coef = 1
        else:
            mut_coef = 2
        num_muts = np.random.binomial(mut_coef * nworms,
                basepairs[loc] * mutation_rate[loc])
        positions = np.copy(worms.pos[str(loc)])
        max_seg = positions[-1]
        for mut in range(num_muts):
            iix = 0
            randmf = np.random.randint(0, nworms)
            newsite = np.random.randint(1, basepairs[loc])
            whap = np.random.randint(1, 3)
            if newsite in worms.pos[str(loc)]: 
                mutation_at_segsite(newsite, loc, worms, randmf)
            else:
                narray = np.zeros(nworms, np.uint8)
                narray[randmf] = 1
                if newsite > max_seg:
                    iix = len(positions)
                else:
                    iix = np.argmax(positions > newsite)
                if recombination_rate[loc] == 0:
                    worms.h1[str(loc)] =\
                            np.insert(worms.h1[str(loc)],
                                    iix, narray, axis=1)
                else:
                    oarray = np.zeros(nworms, np.uint8)
                    whap = np.random.randint(1, 3)
                    if whap == 1: whap2 =str(2)
                    else: whap2 = str(1)
                    whap = str(whap)
                    hap = getattr(worms,"h"+ whap)[str(loc)]
                    hap = np.insert(hap, iix, narray, axis=1)
                    getattr(worms,"h"+ whap)[str(loc)] = hap
                    ohap = getattr(worms, "h"+ whap2)[str(loc)]
                    ohap = np.insert(ohap, iix, oarray, axis=1)
                    getattr(worms, "h"+ whap2)[str(loc)] = ohap
                positions = np.insert(positions, iix, newsite)
                new_positions[str(loc)].append(newsite)
        worms.pos[str(loc)] = positions
        #assert worms.pos[str(loc)].shape[0] == dfAdult_mf.h1[str(loc)].shape[1]
    return(worms, new_positions)
