#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
#import ipdb

from .recombination import recombination_fx
from .mutation import mutation_fx
from .selection import selection_fx

def fecunditybase_fx(fecund,
                     dfAdult,
                     locus,
                     mutation_rate,
                     recombination_rate,
                     basepairs,
                     selection,
                     densitydep_fec):
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
    dfAdult.meta.loc[dfAdult.meta.age < 6, "fec"] = np.random.poisson(fecund,
           len(dfAdult.meta[dfAdult.meta.age < 6]))
    #linear function defining decline in fecundity with age
    m = float(0 - fecund) / (21 - 6)
    b = 0 - m * 21
    #assign fecundity value based on age function
    # :TODO this errors out
    dfAdult.meta.loc[dfAdult.meta.age >= 6, "fec"] = np.random.poisson(m
           * dfAdult.meta.loc[dfAdult.meta.age >= 6,"age"] + b)
    #sex, recombination, mutation
    dfAdult_mf = recombination_fx(locus, dfAdult, recombination_rate, basepairs)
    # Positions is just the new positions
    dfAdult_mf, new_positions = mutation_fx(locus, dfAdult_mf,
         mutation_rate, recombination_rate, basepairs)
    if selection: #dfAdult.sel will be updated here to same length as dfAdult_mf.pos
        dfAdult_mf, dfAdult = selection_fx(dfAdult, dfAdult_mf, new_positions)
    return(dfAdult_mf, dfAdult)
