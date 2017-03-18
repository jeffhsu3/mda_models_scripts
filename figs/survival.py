#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""

import numpy as np
from scipy.stats import weibull_min

from figs.fecundity import fecunditybase_fx
from figs.host_migration import hostmigration_fx


def kill_adults(dfworm, dfHost, month, shapeAdult, scaleAdult,
        village):
    """
    """
    adiix = dfworm.meta.index[dfworm.meta.stage == "A"].values
    #Adult survival is based on weibull cdf
    kill_adult_rand = np.random.random(adiix.shape[0])
    try:
        kill_adult_fx_age = weibull_min.cdf(dfworm.meta.age[adiix], shapeAdult,
                loc=0, scale=scaleAdult)
    except TypeError:
        kill_adult_fx_age = weibull_min.cdf(0, shapeAdult, loc=0, scale=scaleAdult)
    dieAdult = adiix[np.where(kill_adult_rand < kill_adult_fx_age)]
    dfworm.meta.ix[adiix, "age"] += 1 #2 - 21
    ##host survival is from act table
    dfHost = dfHost.query("age < agedeath")
    diehost = dfHost.hostidx.values
    dead_worms = np.append(dieAdult,
            dfworm.meta[~dfworm.meta.hostidx.isin(diehost)].index.values)
    dfworm.drop_worms(dead_worms)
    return(dfHost, dfworm)


def survivalbase_fx(month,
                    village,
                    surv_Juv,
                    shapeMF,
                    scaleMF,
                    shapeAdult,
                    scaleAdult,
                    fecund,
                    locus,
                    mutation_rate,
                    recombination_rate,
                    basepairs,
                    selection,
                    hostmigrate,
                    mdalist,
                    densitydep_surv,
                    densitydep_fec,
                    dfHost,
                    dfworm,
                    R0netlist,
                    cdslist):
    '''Base survival function
    Parameters
    ---------
    month: int
         time in months of the simulation
    surv_Juv: float
         survival prob of juvenille life stage
    shapeMF: float
         shape parameter for weibull distribution of MF
    scaleMF: int
         scale parameter for weibull distribution of MF
    shapeAdult: float
         shape parameter for weibull distribution of Adults
    scaleAdult: int
         scale parameter for weibull distribution of MF
    dfMF : df
        dataframe of MF
    dfAdult : df
        dataframe of adult worms
    dfJuv : df
        dataframe of juvenille worms
    dfHost : df
        dataframe of hosts
    basepairs : int, list
        length of loci
    hostmigrate : float
        rate of migration per year between villages
    selection : boolean
        T/F for selection
    dfSel : df
        dataframe of cds positions and fitness
    cds_coordinates : list
        list of coding seq coordinates
    mdalist : list
        list of mda parameters
    densitydep_surv : boolean
    densitydep_fec : boolean

    Returns
    -------
    dfworm
    dfHost

    '''
    if month%12 == 0:
        ##stats
        x = dfworm.meta.groupby(["village","stage"]).apply(lambda y: y[(y.R0net < (len(R0netlist['R0']) + 1))
                            & (y.R0net > len(R0netlist['R0']))]).R0net[:,'A']
        R0netlist['R0'].append([len(x[i]) for i in range(len(x.index.levels[0]))])
        R0netlist['repoavg'].append([np.mean((np.unique(x[i],return_counts=True)[1])) for i in range(len(x.index.levels[0]))])
        R0netlist['repovar'].append([np.var((np.unique(x[i],return_counts=True)[1])) for i in range(len(x.index.levels[0]))])
        ##
        dfHost, dfworm = kill_adults(dfworm, dfHost, month, shapeAdult, scaleAdult,
                village)

        dfHost.age = dfHost.age + 1
        if hostmigrate != 0:
            dfHost = hostmigration_fx(village, dfHost, hostmigrate)
    else: pass

    ##Juv is exponential 0.866; surv_Juv
    juviix = dfworm.meta[dfworm.meta.stage == "J"].index.values
    kill_juvrand = np.random.random(juviix.shape[0])
    dieJuv = juviix[np.where(kill_juvrand > surv_Juv)]
    dfworm.meta.ix[juviix,'age'] += 1

    ##MF is weibull cdf
    mfiix = dfworm.meta[dfworm.meta.stage == "M"].index.values
    kill_mfrand = np.random.random(mfiix.shape[0])
    try:
        kill_mffxage = weibull_min.cdf(dfworm.meta.ix[mfiix].age,shapeMF,loc=0,scale=scaleMF)
    except TypeError:
        kill_mffxage = weibull_min.cdf(0,shapeMF,loc=0,scale=scaleMF)
    dieMF = mfiix[np.where(kill_mfrand < kill_mffxage)]
    dfworm.meta.ix[mfiix, 'age'] += 1

    ##move Juv age 13 to adult age 1
    juviix12 = dfworm.meta.ix[juviix].query('age > 12').index.values
    if any(juviix12):
        #reset age to adult
        dfworm.meta.ix[juviix12,'age'] = 1
        #increase R0net for next gen
        dfworm.meta.ix[juviix12,'R0net'] += 1
        dfworm.meta.ix[juviix12,'stage'] = "A"

    dfworm.drop_worms(np.append(dieJuv, dieMF))
    #fecundity calls mutation/recombination
    dfAdult_mf, dfworm = fecunditybase_fx(fecund, dfworm, locus, mutation_rate,
                                         recombination_rate, basepairs, selection,
                                         densitydep_fec, cdslist)

    return(dfHost, dfworm, R0netlist)
