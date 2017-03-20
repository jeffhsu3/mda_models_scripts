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
from figs.hostmda import hostmda_fx
from figs.fecundity_mda import fecunditymda_fx
from figs.fecundity_mda import fecunditymda_sel1_fx
from figs.fecundity_mda import fecunditymda_sel2_fx
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

def survivalmda_fx(month,
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

    '''base survival function
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
    mda : boolean
         mda in village or not
    mda_start : list, int
         time in months that mda began in each village
    mda_num : list, int
         how many mdas to simulate
    mda_freq : int
         how often
    mda_coverage : float
         what percentage of the population get mda
    mda_macro: float
         percent of adults killed by drug
    mda_micro: float
         percent of MF killed by drug
    mda_juvcide: float
         percent of juveniile killed by drug; could be 0 or avg micro/macro
    mda_sterile : float
         percent of adult worms perm sterilized
    mda_clear : int
         clearance time

    Returns
    -------
    dfMF
    dfAdult
    dfJuv
    dfHost
    dfSel

   '''
    mda_start = mdalist[0]
    mda_num = mdalist[1]
    mda_freq = mdalist[2]
    mda_coverage = mdalist[3]
    mda_macro = mdalist[4]
    mda_juvicide = mdalist[5]
    mda_micro = mdalist[6]
    mda_sterile = mdalist[7]
    mda_clear = mdalist[8]

    ##calculates time since mda as clear_count
    if month < mda_start or month > (mda_start + mda_freq * mda_num):
        clear_count = 0
    else:
        if month%mda_freq == 0:
             clear_count = 1
             print("\n\nStart of MDA\n\n")
        elif (month - mda_clear)%mda_freq == 0:
             clear_count = mda_clear + 1
        else:
             if ((month-mda_freq)%mda_freq + 1) <= mda_clear:
                  clear_count = ((month-mda_freq)%mda_freq + 1)
             else:
                  clear_count = 0
    print(clear_count)
#######################
    ##assign MDA to hosts
    if clear_count == 1:
        ##add mda to dfHost per coverage
        dfHost = hostmda_fx(village, dfHost, mda_coverage)
        ##apply MDA effect to host populations of worm
        hostmda = dfHost[dfHost.MDAstate == 1].hostidx.values
        wormmda = dfworm.meta[dfworm.meta.hostidx.isin(hostmda)].index.values
        MFmda = dfworm.meta.ix[wormmda].query("stage == 'M'").index.values
        Juvmda = dfworm.meta.ix[wormmda].query("stage == 'J'").index.values
        Adultmda = dfworm.meta.ix[wormmda].query("stage == 'A'").index.values
        mfdiemda = dfworm.meta.ix[MFmda].groupby("hostidx").apply(lambda y: y.sample(frac = mda_micro).index.values)
        juvdiemda = dfworm.meta.ix[Juvmda].groupby("hostidx").apply(lambda y: y.sample(frac = mda_juvicide).index.values)
        adultdiemda = dfworm.meta.ix[Adultmda].groupby("hostidx").apply(lambda y: y.sample(frac = mda_macro).index.values)
        diemda = [i for i in mfdiemda.values] + [i for i in juvdiemda.values] + [i for i in adultdiemda.values]
        dfworm.drop_worms(np.concatenate(diemda).ravel())
    elif clear_count == 0:
        hostmda = ''
    else:
        hostmda = dfHost[dfHost.MDAstate == 1].hostidx.values
#######################
    if month%12 == 0:
        #stats
        x = dfworm.meta.groupby(["village","stage"]).apply(lambda y: y[(y.R0net < (len(R0netlist['R0']) + 1))
                            & (y.R0net > len(R0netlist['R0']))]).R0net[:,'A']
        R0netlist['R0'].append([len(x[i]) for i in range(len(x.index.levels[0]))])
        R0netlist['repoavg'].append([np.mean((np.unique(x[i],return_counts=True)[1])) for i in range(len(x.index.levels[0]))])
        R0netlist['repovar'].append([np.var((np.unique(x[i],return_counts=True)[1])) for i in range(len(x.index.levels[0]))])

        dfHost, dfworm = kill_adults(dfworm, dfHost, month, shapeAdult, scaleAdult,
                village)

        dfHost.age = dfHost.age + 1
        hostmignumb = np.random.poisson(hostmigrate)
        if hostmignumb != 0:
            dfHost = hostmigration_fx(village, dfHost, hostmignumb)
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
    else:pass
    dfworm.drop_worms(np.append(dieJuv, dieMF))
    #fecundity calls mutation/recombination
    dfAdult_mf, dfworm = fecunditymda_fx(fecund, dfworm, locus, mutation_rate,
                                         recombination_rate, basepairs, selection,
                                         densitydep_fec, cdslist, mda_sterile, clear_count, mda_clear,
                                         hostmda)
    return(dfHost, dfworm, R0netlist)

def survivalmda_sel1_fx(month,
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

    '''

   '''
    mda_start = mdalist[0]
    mda_num = mdalist[1]
    mda_freq = mdalist[2]
    mda_coverage = mdalist[3]
    mda_macro = mdalist[4]
    mda_juvicide = mdalist[5]
    mda_micro = mdalist[6]
    mda_sterile = mdalist[7]
    mda_clear = mdalist[8]

    ##calculates time since mda as clear_count
    if month < mda_start or month > (mda_start + mda_freq * mda_num):
            clear_count = 0
    else:
        if month%mda_freq == 0:
             clear_count = 1
             print("\n\nStart of MDA\n\n")
        elif (month - mda_clear)%mda_freq == 0:
             clear_count = mda_clear + 1
        else:
             if ((month-mda_freq)%mda_freq + 1) <= mda_clear:
                  clear_count = ((month-mda_freq)%mda_freq + 1)
             else:
                  clear_count = 0
##################################
    ##assign MDA to hosts
    if clear_count == 1:
        ##add mda to dfHost per coverage
        dfHost = hostmda_fx(village, dfHost, mda_coverage)
        ##apply MDA effect to host populations of worms
        hostmda = dfHost[dfHost.MDAstate == 1].hostidx.values
        wormmda = dfworm.meta[dfworm.meta.hostidx.isin(hostmda)].index.values
        MFmda = dfworm.meta.ix[wormmda].query("stage == 'M'").index.values
        Juvmda = dfworm.meta.ix[wormmda].query("stage == 'J'").index.values
        Adultmda = dfworm.meta.ix[wormmda].query("stage == 'A'").index.values
        #import ipdb; ipdb.set_trace()
        mfdiemda = dfworm.meta.ix[MFmda].groupby("hostidx").apply(lambda y: y[np.random.random(y.shape[0]) < (mda_micro ** y.fitS)].index.values)
        juvdiemda = dfworm.meta.ix[Juvmda].groupby("hostidx").apply(lambda y: y[np.random.random(y.shape[0]) < (mda_juvicide ** y.fitS)].index.values)
        adultdiemda = dfworm.meta.ix[Adultmda].groupby("hostidx").apply(lambda y: y[np.random.random(y.shape[0]) < (mda_macro ** y.fitS)].index.values)
        diemda = [i for i in mfdiemda.values] + [i for i in juvdiemda.values] + [i for i in adultdiemda.values]
        dfworm.drop_worms(np.concatenate(diemda).ravel())
    elif clear_count == 0:
        hostmda = ''
    else:
        hostmda = dfHost[dfHost.MDAstate == 1].hostidx.values
#########################################
    if month%12 == 0:
        x = dfworm.meta.groupby(["village","stage"]).apply(lambda y: y[(y.R0net < (len(R0netlist['R0']) + 1))
                            & (y.R0net > len(R0netlist['R0']))]).R0net[:,'A']
        R0netlist['R0'].append([len(x[i]) for i in range(len(x.index.levels[0]))])
        R0netlist['repoavg'].append([np.mean((np.unique(x[i],return_counts=True)[1])) for i in range(len(x.index.levels[0]))])
        R0netlist['repovar'].append([np.var((np.unique(x[i],return_counts=True)[1])) for i in range(len(x.index.levels[0]))])

        dfHost, dfworm = kill_adults(dfworm, dfHost, month, shapeAdult, scaleAdult,
                village)

        dfHost.age = dfHost.age + 1
        hostmignumb = np.random.poisson(hostmigrate)
        if hostmignumb != 0:
            dfHost = hostmigration_fx(village, dfHost, hostmignumb)
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
    else:pass
    dfworm.drop_worms(np.append(dieJuv, dieMF))
    #fecundity calls mutation/recombination
    dfAdult_mf, dfworm = fecunditymda_sel1_fx(fecund, dfworm, locus, mutation_rate,
                                         recombination_rate, basepairs, selection,
                                         densitydep_fec, cdslist, mda_sterile, clear_count, mda_clear,
                                         hostmda)
    return(dfHost, dfworm, R0netlist)

def survivalmda_sel2_fx(month,
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

    '''base survival function

   '''
    mda_start = mdalist[0]
    mda_num = mdalist[1]
    mda_freq = mdalist[2]
    mda_coverage = mdalist[3]
    mda_macro = mdalist[4]
    mda_juvicide = mdalist[5]
    mda_micro = mdalist[6]
    mda_sterile = mdalist[7]
    mda_clear = mdalist[8]

    ##calculates time since mda as clear_count
    if month < mda_start or month > (mda_start + mda_freq * mda_num):
            clear_count = 0
    else:
        if month%mda_freq == 0:
             clear_count = 1
             print("\n\nStart of MDA\n\n")
        elif (month - mda_clear)%mda_freq == 0:
             clear_count = mda_clear + 1
        else:
             if ((month-mda_freq)%mda_freq + 1) <= mda_clear:
                  clear_count = ((month-mda_freq)%mda_freq + 1)
             else:
                  clear_count = 0
##################################
    ##assign MDA to hosts
    if clear_count == 1:
        ##add mda to dfHost per coverage
        dfHost = hostmda_fx(village, dfHost, mda_coverage)
        ##apply MDA effect to host populations of worms
        hostmda = dfHost[dfHost.MDAstate == 1].hostidx.values
        wormmda = dfworm.meta[dfworm.meta.hostidx.isin(hostmda)].index.values
        MFmda = dfworm.meta.ix[wormmda].query("stage == 'M'").index.values
        Juvmda = dfworm.meta.ix[wormmda].query("stage == 'J'").index.values
        Adultmda = dfworm.meta.ix[wormmda].query("stage == 'A'").index.values
        #import ipdb; ipdb.set_trace()
        mfdiemda = dfworm.meta.ix[MFmda].groupby("hostidx").apply(lambda y: y[np.random.random(y.shape[0]) < (mda_micro ** y.fitS)].index.values)
        juvdiemda = dfworm.meta.ix[Juvmda].groupby("hostidx").apply(lambda y: y[np.random.random(y.shape[0]) < (mda_juvicide ** y.fitS)].index.values)
        adultdiemda = dfworm.meta.ix[Adultmda].groupby("hostidx").apply(lambda y: y[np.random.random(y.shape[0]) < (mda_macro ** y.fitS)].index.values)
        diemda = [i for i in mfdiemda.values] + [i for i in juvdiemda.values] + [i for i in adultdiemda.values]
        dfworm.drop_worms(np.concatenate(diemda).ravel())
    elif clear_count == 0:
        hostmda = ''
    else:
        hostmda = dfHost[dfHost.MDAstate == 1].hostidx.values
########################################
    if month%12 == 0:
        x = dfworm.meta.groupby(["village","stage"]).apply(lambda y: y[(y.R0net < (len(R0netlist['R0']) + 1))
                            & (y.R0net > len(R0netlist['R0']))]).R0net[:,'A']
        R0netlist['R0'].append([len(x[i]) for i in range(len(x.index.levels[0]))])
        R0netlist['repoavg'].append([np.mean((np.unique(x[i],return_counts=True)[1])) for i in range(len(x.index.levels[0]))])
        R0netlist['repovar'].append([np.var((np.unique(x[i],return_counts=True)[1])) for i in range(len(x.index.levels[0]))])

        adiix = dfworm.meta.index[dfworm.meta.stage == "A"].values
        #Adult survival is based on weibull cdf
        kill_adult_rand = np.random.random(adiix.shape[0])
        try:
            kill_adult_fx_age = weibull_min.cdf(dfworm.meta.age[adiix], shapeAdult,
                    loc=0, scale=scaleAdult)
        except TypeError:
            kill_adult_fx_age = weibull_min.cdf(0, shapeAdult, loc=0, scale=scaleAdult)
        dieAdult = adiix[np.where(kill_adult_rand < (kill_adult_fx_age ** (1-abs(1-dfworm.meta.ix[adiix, "fitS"]))))]
        dfworm.meta.ix[adiix, "age"] += 1 #2 - 21
        ##host survival is from act table
        dfHost = dfHost.query("age < agedeath")
        diehost = dfHost.hostidx.values

        dead_worms = np.append(dieAdult,
                dfworm.meta[~dfworm.meta.hostidx.isin(diehost)].index.values)
        dfworm.drop_worms(dead_worms)

        dfHost.age = dfHost.age + 1
        hostmignumb = np.random.poisson(hostmigrate)
        if hostmignumb != 0:
            dfHost = hostmigration_fx(village, dfHost, hostmignumb)
    else: pass

    ##Juv is exponential 0.866; surv_Juv
    juviix = dfworm.meta[dfworm.meta.stage == "J"].index.values
    kill_juvrand = np.random.random(juviix.shape[0])
    dieJuv = juviix[np.where(kill_juvrand > (surv_Juv ** (1-abs(1-dfworm.meta.ix[juviix, "fitS"]))))]
    dfworm.meta.ix[juviix,'age'] += 1

    ##MF is weibull cdf
    mfiix = dfworm.meta[dfworm.meta.stage == "M"].index.values
    kill_mfrand = np.random.random(mfiix.shape[0])
    try:
        kill_mffxage = weibull_min.cdf(dfworm.meta.ix[mfiix].age,shapeMF,loc=0,scale=scaleMF)
    except TypeError:
        kill_mffxage = weibull_min.cdf(0,shapeMF,loc=0,scale=scaleMF)
    dieMF = mfiix[np.where(kill_mfrand < (kill_mffxage ** (1-abs(1-dfworm.meta.ix[mfiix, "fitS"]))))]
    dfworm.meta.ix[mfiix, 'age'] += 1
############################################
    ##move Juv age 13 to adult age 1
    juviix12 = dfworm.meta.ix[juviix].query('age > 12').index.values
    if any(juviix12):
        #reset age to adult
        dfworm.meta.ix[juviix12,'age'] = 1
        #increase R0net for next gen
        dfworm.meta.ix[juviix12,'R0net'] += 1
        dfworm.meta.ix[juviix12,'stage'] = "A"
    else:pass
    dfworm.drop_worms(np.append(dieJuv, dieMF))
    #fecundity calls mutation/recombination
    dfAdult_mf, dfworm = fecunditymda_sel2_fx(fecund, dfworm, locus, mutation_rate,
                                         recombination_rate, basepairs, selection,
                                         densitydep_fec, cdslist, mda_sterile, clear_count, mda_clear,
                                         hostmda)
    return(dfHost, dfworm, R0netlist)
