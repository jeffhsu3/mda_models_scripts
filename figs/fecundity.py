#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
from .recombination import recombination_fx
from .mutation import mutation_fx
from .selection import selection_fx

def fecunditybase_fx(fecund,
                     dfworm,
                     locus,
                     mutation_rate,
                     recombination_rate,
                     basepairs,
                     selection,
                     densitydep_fec,
                     cdslist):
    '''Base fecundity function, simpliest scenario

    Parameters
    ---------
    fecund: int
        rate of the poisson distribution for offspring number
    dfAdult: figs.Worms
        dataframe containing adult parasites
    locus : int
    mutation_rate : list, float recombination_rate : list, float basepairs : list, int
    selection : boolean
    dfSel : df
    cds_coordinates : [list], int

    Returns
    ------
    dfAdult_mf : df
        figs.worms object
    dfSel : df

    '''
    adiix = dfworm.meta[dfworm.meta.stage == "A"].index.values
    young = adiix[np.where(dfworm.meta.ix[adiix].age < 6)]
    old = adiix[np.where(dfworm.meta.ix[adiix].age >= 6)]
    dfworm.meta.ix[young, 'fec'] = np.random.poisson(fecund, young.shape[0])
    #linear function defining decline in fecundity with age
    m = float(0 - fecund) / (21 - 6)
    b = 0 - m * 21
    #assign fecundity value based on age function
    positive_lambda = (dfworm.meta.ix[old].age.values * m) + b
    positive_lambda[positive_lambda < 0] = 0
    dfworm.meta.ix[old, 'fec'] = np.random.poisson(positive_lambda).astype(np.int64)
    #sex, recombination, mutation
    dfAdult_mf = recombination_fx(locus, dfworm, adiix, recombination_rate, basepairs)
    dfAdult_mf.meta.sex = np.random.choice(['M', 'F'] , size=len(dfAdult_mf.meta))
    dfAdult_mf.meta.age = 1
    dfAdult_mf, new_positions = mutation_fx(locus, dfAdult_mf,
         mutation_rate, recombination_rate, basepairs)
    dfworm.add_worms(dfAdult_mf, new_positions)
    if selection: #dfAdult.sel will be updated here to same length as dfAdult_mf.pos
        dfAdult_mf, dfworm = selection_fx(dfworm, dfAdult_mf, new_positions, cdslist)
    return(dfAdult_mf, dfworm)
