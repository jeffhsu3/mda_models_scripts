#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
#import ipdb
import numpy as np
import pandas as pd
from scipy.stats import weibull_min
import random

from figs.fecundity import fecunditybase_fx
from figs.host_migration import hostmigration_fx

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
                    dfAdult,
                    dfJuv,
                    dfMF):


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
    dfMF
    dfAdult
    dfJuv
    dfHost
    dfSel

    '''

    #adult worms and hosts are only evaluated per year
    if month%12 == 0:
        #Adult survival is based on weibull cdf
        kill_adult_rand = np.random.random(dfAdult.meta.shape[0])
        try:
            kill_adultfxage = weibull_min.cdf(dfAdult.meta.age, shapeAdult,
                    loc=0, scale=scaleAdult)
        except TypeError:
            kill_adultfxage = weibull_min.cdf(0, shapeAdult, loc=0, scale=scaleAdult)
        dieAdult = np.where(kill_adult_rand < kill_adultfxage)
        dfAdult.drop_worms(dieAdult[0])

        dfAdult.meta.age = dfAdult.meta.age + 1 #2 - 21
        ##host survival is from act table
        dfHost = dfHost[dfHost.age < dfHost.agedeath]
        #remove all worms with dead host.hostidx from all dataframes
        dfAdult.drop_worms(dfAdult.meta[~dfAdult.meta.hostidx.isin(dfHost.hostidx)].index.values)
        dfJuv.drop_worms(dfJuv.meta[~dfJuv.meta.hostidx.isin(dfHost.hostidx)].index.values)
        dfMF.drop_worms(dfMF.meta[~dfMF.meta.hostidx.isin(dfHost.hostidx)].index.values)
        #add 1 year to all ages of hosts
        dfHost.age = dfHost.age + 1
        if hostmigrate != 0:
            dfHost = hostmigration_fx(village, dfHost, hostmigrate)

    print('Post killing')
    print(dfAdult.pos['1'].shape[0])
    print(dfAdult.h1['1'].shape[1])
    print('Nworms :{0!s}'.format(dfAdult.h1['1'].shape[0]))

    ##Juv is exponential 0.866; surv_Juv
    #dont include age 0 which just moved from transmission fx
    dfJuv.meta.age += 1
    kill_juvrand = np.random.random(dfJuv.meta.shape[0])
    dieJuv = np.where(kill_juvrand > surv_Juv)
    dfJuv.drop_worms(dieJuv[0])

    ##MF is weibull cdf
    kill_mfrand = np.random.random(dfMF.meta.shape[0])
    try:
        kill_mffxage = weibull_min.cdf(dfMF.meta.age,shapeMF,loc=0,scale=scaleMF)
    except TypeError:
        kill_mffxage = weibull_min.cdf(0,shapeMF,loc=0,scale=scaleMF)
    dieMF = np.where(kill_mfrand < kill_mffxage)
    dfMF.drop_worms(dieMF[0])
    dfMF.meta.age = dfMF.meta.age + 1 #2 - 12
    dfMF.drop_worms(dfMF.meta.ix[dfMF.meta.age > 12].index.values) #hard cutoff at 12 months
    try:
        assert dfMF.pos['0'].shape[0] == dfMF.h1['0'].shape[1]
    except KeyError:
        pass

    ##move Juv age 13 to adult age 1
    juv_rows = dfJuv.meta[dfJuv.meta.age > 12].index.values
    try:
        #reset age to adult
        dfJuv.meta.ix[juv_rows,"age"] = 1
        #increase R0net for next gen
        dfJuv.meta.ix[juv_rows, "R0net"] += 1
    except TypeError:
        print("dfJuv empty")
    dfJuv.drop_worms(juv_rows)
    dfAdult.add_worms(dfJuv, juv_rows)

    #fecundity calls mutation/recombination
    dfAdult_mf, dfAdult = fecunditybase_fx(fecund, dfAdult, locus, mutation_rate,
                                         recombination_rate, basepairs, selection,
                                         densitydep_fec)
    dfAdult_mf.meta.sex = [random.choice("MF") for i in range(len(dfAdult_mf.meta))]
    dfAdult_mf.meta.age = 1
    try:
        dfMF.add_worms(dfAdult_mf, dfAdult_mf.meta.index.values)
    except ValueError:
        from IPython import embed
        embed()


    return(dfHost, dfAdult, dfJuv, dfMF)
