#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import copy


def fitness_fx(dfAdult_mf,
               dfworm,
               cdslist):
    ''' calculates mean fitness for each individual by summing fitness effects
    from dfSel for each position across all loci.
    selection coefficient, s, is found by fit - 1
    Parameters
    ---------
    locus : int
        number of loci
    dfAdult_mf : df
         data frame of adult worms containing genotype information
    dfSel : df
         data fram of fitness benefit for each allele
    Returns
    ------
    dfAdult_mf : df
         updated df for mf
    '''

    ninds = len(dfAdult_mf.meta)
    dom = cdslist[4]
    fitF_ind = np.zeros(ninds)
    fitS_ind = np.zeros(ninds)
    Floc = 0
    Sloc = 0
    for locus in dfAdult_mf.h2.keys():
        count_sites = dfAdult_mf.h1[locus] + dfAdult_mf.h2[locus]
        if dom == "codom":
        ##CODOMINANT: AA = 1 + 2s, Aa = 1 + hs, aa = 1, where h is 1
            h = 1
        elif dom == "dom":
        ##DOMINANT: AA = 1 + 2s, Aa = 1 + hs, aa = 1, where h is 2
            h = 2
            count_sites[count_sites > 0] = 2
        else:
        ##RECESSIVE: AA = 1 + 2s, Aa = 1 + hs, aa = 1, where h is 0
            h = 0
            count_sites[count_sites < 2] = 0

        sum_selsites_S = np.dot(count_sites, dfworm.sel[locus + "S"])
        sum_selsites_F = np.dot(count_sites, dfworm.sel[locus + "F"])

        intsites_S = copy.copy(dfworm.sel[locus + "S"])
        intsites_S[intsites_S > 0] = 1
        intsites_F = copy.copy(dfworm.sel[locus + "F"])
        intsites_F[intsites_F > 0] = 1

        cds_sites_S = np.dot(dfAdult_mf.h1[locus], intsites_S) \
            + np.dot(dfAdult_mf.h2[locus], intsites_S)
        cds_sites_F = np.dot(dfAdult_mf.h1[locus], intsites_F) \
            + np.dot(dfAdult_mf.h2[locus], intsites_F)
        fitS_ind += (( (dfworm.sel[locus + "St"] * 2) - cds_sites_S) + sum_selsites_S) / (dfworm.sel[locus + "St"] * 2)
        fitF_ind += (( (dfworm.sel[locus + "Ft"] * 2) - cds_sites_F) + sum_selsites_F) / (dfworm.sel[locus + "Ft"] * 2)
        Floc += dfworm.coord[locus + "F"].shape[0]
        Sloc += dfworm.coord[locus + "S"].shape[0]
    #import ipdb; ipdb.set_trace()
    dfAdult_mf.meta["fitS"] = fitS_ind / Sloc
    dfAdult_mf.meta["fitF"] = fitF_ind / Floc
    return(dfAdult_mf)

def selection_fx(dfworm,
                 dfAdult_mf,
                 new_positions,
                 cdslist):
    '''recalculates DFE for new mutations and phenotype for new mf
    Parameters
    ---------
    dfSel : df
         updates this dataframe
    dfAdult_mf : df
         adds phenotype info
    dfMuts : df
         gives positions of new mutations for addition to dfSel
    locus : int
         number of loci
    Returns
    ------
    dfSel : df
         now updates with new mutations
    dfAdult_mf : df
         updated with phenotype
    '''
    for loc in dfworm.h2.keys(): #since this wont include 0
        new_positions[loc].sort()
        for pos in new_positions[loc]: #this is the dict of positions
            iix = np.argmax(dfworm.pos[loc] == pos) #this is the position it is inserted
            if any([i <= pos <= j for i,j in dfworm.coord[loc + "F"]]):
                 #shape = 4, mean = 1, scale = mean/shape
                 #here mean is mean_fitness, wildtype is assumed to be 1
                 selF = np.random.gamma(4, scale=0.25)
                 selS = 0
            elif any([i <= pos <= j for i,j in dfworm.coord[loc + "S"]]):
                     selS = np.random.gamma(4, scale=0.25)
                     selF = 0
            else: #not in a cds
                selS = 0
                selF = 0
            dfworm.sel[loc + "S"] = np.insert(dfworm.sel[loc + "S"], iix, selS)
            dfworm.sel[loc + "F"] = np.insert(dfworm.sel[loc + "F"], iix, selF)

    dfAdult_mf = fitness_fx(dfAdult_mf, dfworm, cdslist)
    return(dfAdult_mf, dfworm)
