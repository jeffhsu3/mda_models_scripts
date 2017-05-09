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
    ##host survival is from act table dfHost = dfHost.query("age < agedeath")
    diehost = dfHost.hostidx.values
    dead_worms = np.append(dieAdult,
            dfworm.meta[~dfworm.meta.hostidx.isin(diehost)].index.values)
    dfworm.drop_worms(dead_worms)
    return(dfHost, dfworm)


def kill_juvenile(dfworm, surv_Juv, increment_age=False):
    ''' Returns bool of worms to kill. Also increments age
    '''
    juviix = dfworm.meta[dfworm.meta.stage == "J"].index.values
    kill_juvrand = np.random.random(juviix.shape[0])
    dieJuv = juviix[np.where(kill_juvrand > surv_Juv)]
    if increment_age:
        dfworm.meta.loc[juviix, 'age'] += 1
    else: pass
    return(dieJuv)


def kill_mf(dfworm, shapeMF, scaleMF, increment_age=False):
    ''' Kills MF
    '''
    mfiix = dfworm.meta[dfworm.meta.stage == "M"].index.values
    kill_mfrand = np.random.random(mfiix.shape[0])
    try:
        kill_mffxage = weibull_min.cdf(dfworm.meta.ix[mfiix].age,
                shapeMF, loc=0, scale=scaleMF)
    except TypeError:
        kill_mffxage = weibull_min.cdf(0, shapeMF, loc=0, scale=scaleMF)
    dieMF = mfiix[np.where(kill_mfrand < kill_mffxage)]
    if increment_age:
        dfworm.meta.ix[mfiix, 'age'] += 1
    else: pass
    return(dieMF)


def age_juvenile(dfworm):
    ''' Ages juveniles
    '''
    juviix = dfworm.meta[dfworm.meta.stage == "J"].index.values
    juviix12 = dfworm.meta.ix[juviix].query('age > 12').index.values
    if any(juviix12):
        #reset age to adult
        dfworm.meta.ix[juviix12,'age'] = 1
        #increase R0net for next gen
        dfworm.meta.ix[juviix12,'R0net'] += 1
        dfworm.meta.ix[juviix12,'stage'] = "A"
    else:pass


def add_only_variants(dfworm, dfworm_to_add):
    to_add = dfworm.meta.stage == 'A'
    if sum(to_add) >= 1:
        pass
    else: pass
    # :TODO check if nothing
    if any(to_add):
        dfworm_to_add.add_worms(dfworm, to_add) 
        dfworm.drop_worms(to_add)
    else:
        pass

    
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
                    cdslist, 
                    juvs):
    '''
    Base survival function
    
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
        length of locI
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
    dfworm :
    dfHost :
    new_pos_iix :
    R0netlist :


    '''
    ### Need to move this to function above
    new_worms = juvs['worms']
    new_indexes = juvs['indexes']
    new_positions = juvs['pos']
    if month%12 == 0:
        ##stats
        x = dfworm.meta.groupby(["village","stage"]).apply(
                lambda y: y[(y.R0net < (len(R0netlist['R0']) + 1))\
                        & (y.R0net > len(R0netlist['R0']))]).R0net[:,'A']
        R0netlist['R0'].append([len(x[i]) for i in range(len(x.index.levels[0]))])
        R0netlist['repoavg'].append([np.mean((np.unique(x[i],return_counts=True)[1]))\
                 for i in range(len(x.index.levels[0]))])
        R0netlist['repovar'].append([np.var((np.unique(x[i],return_counts=True)[1]))\
                 for i in range(len(x.index.levels[0]))])
        ##
        dfHost, dfworm = kill_adults(dfworm, dfHost, month, 
                shapeAdult, scaleAdult, village)
        dfHost.age = dfHost.age + 1
        hostmignumb = np.random.poisson(hostmigrate)
        dfHost = hostmigration_fx(village, dfHost, hostmignumb)
    else: pass


    #### df worm #########################
    dieJuv = kill_juvenile(dfworm, surv_Juv, increment_age=True)
    dieMF = kill_mf(dfworm, shapeMF, scaleMF, increment_age=True)
    dfworm.drop_worms(np.append(dieJuv, dieMF))
    age_juvenile(dfworm)

    for i, j in zip(new_worms, new_positions):
        dieJuv = kill_juvenile(i, surv_Juv, increment_age=True)
        dieMF = kill_mf(i, shapeMF, scaleMF, increment_age=True)
        i.drop_worms(np.append(dieJuv, dieMF))
        age_juvenile(i)
        add_only_variants(i, dfworm)


    #fecundity calls mutation/recombination
    df_new_worms, dfworm, new_pos = fecunditybase_fx(fecund, dfworm, locus, mutation_rate,
                                         recombination_rate, basepairs, selection,
                                         densitydep_fec, cdslist)
    ##### NEW WORMS ######################
    juvs['worms'].append(df_new_worms)
    juvs['pos'].append(new_pos)
    return(dfHost, dfworm, juvs, R0netlist)
