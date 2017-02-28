#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import pandas as pd
from collections import Counter

def R0net_fx(dfAdult, dfMF, dfJuv):
    '''Calculates the reproductive number, R0, by counting the uniqueness
    of R0net per village and taking the mean counts

    Parameters
    ----------
    dfJuv : df
        dataframe of juvenilles age 13

    Returns
    --------
    R0 : float, list
        reproductive rate of each village
    '''
    #cant remember whether to calculate this as Juv age 13 or adult age 0
    dfJuv[dfJuv.age == 13].groupby(["village", "R0net"]).size()[1].mean()
    R0 = 1
    return(R0)


def allelefreq_fx(dfAdult, dfSel):
    ''' Calculates allele freq for mutations in dfSel by village

    Parameters
    ---------
    dfAdult : df
        data frame of adult worms
    dfSel : df
        data frame of mutation positions at each locus

    Returns
    -------
    dfFreq : list, float
        freq of mutations in dfSel for each village
    '''
    locus = []
    positions = []
    freqv1 = []
    freqv2 = []
    freqv3 = []
    freqv4 = []
    #list of things to count from dfSel per locus, position
    pos = dfSel.groupby("locus").position.apply(lambda x: x.values.tolist())
    #counts of each mutation in dfAdult by village for each locus
    #locus = len(positions)
    for loc in range(1,len(pos) + 1):
        hap1 = dfAdult.groupby("village")["locus_" + str(loc) + "_h1"].apply(lambda x: [i for l in
                       x.values.tolist() for i in l])
        hap2 = dfAdult.groupby("village")["locus_" + str(loc) + "_h2"].apply(lambda x: [i for l in
                       x.values.tolist() for i in l])        #filter hap1,hap2 by positions in dfSel
        #combine counts for same positions in hap1 and hap2 at each locus
        hap3 = hap1 + hap2
        villcounts = dfAdult.groupby("village").size() * 2
        ##To Do: Generalize to any number of villages, although current limit is 4 and loci
        locus.extend(np.repeat(loc, len(pos[loc])))

        freqv1.extend([round(hap3[0].count(i) / float(villcounts[0]),2) for i in pos[loc]])
        freqv2.extend([round(hap3[1].count(i) / float(villcounts[1]),2) for i in pos[loc]])
        positions.extend(pos[loc])

    dfFreq = pd.DataFrame({"locus" : np.repeat(range(1,len(pos)+1), [len(i) for i in pos]),
                       "position" : positions,
                       "freqv1" : freqv1,
                       "freqv2" : freqv2
                       })
    dfFreq = dfFreq.loc[:, ['locus', 'position', 'freqv1','freqv2']]
    return(dfFreq)
