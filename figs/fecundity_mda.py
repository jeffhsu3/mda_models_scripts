#!/usr/bin/env python2
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

def fecunditymda_fx(fecund,
                    dfworm,
                    locus,
                    mutation_rate,
                    recombination_rate,
                    basepairs,
                    selection,
                    densitydep_fec,
                    cdslist,
                    mda_sterile,
                    clear_count,
                    mda_clear,
                    hostmda):

    '''function for reduced fecundity under mda
    conditions: mda=True, selection=False

    Parameters
    ----------
    fecund: int
         rate of the poisson distribution for offspring number
    sterile_p: float
         proportion of the adult parasite that are permanently sterilized from drugs
    clear_time: int
         number of months the drugs is effective or present in the host
    clear_count: int
         count of months since drug administered
    dfAdult: dataframe
         pandas dataframe containing adult parasites
    dfMF: dataframe
          pandas dataframe containing larval stage parasites

    Returns
    ------
    dfAdult
    dfMF
    '''
####################################################
    if clear_count == 1: #permanent sterility
        adultmda = dfworm.meta[(dfworm.meta["stage"] == "A") &
                               (dfworm.meta.hostidx.isin(hostmda))].index.values
        adultsterile = dfworm.meta.ix[adultmda].groupby("hostidx").apply(lambda y: y.sample(frac = mda_sterile).index.values)
        steriix = np.concatenate([i for i in adultsterile.values]).ravel()
        dfworm.meta.ix[steriix,"sex"] = "S"
    else: pass

##reduced fertility only for hosts with MDA
    if clear_count > 0 and clear_count <= mda_clear: #Drugs cause temporary sterility over clear_time
        adiix = dfworm.meta[dfworm.meta.stage == "A"].index.values
        mda_true = dfworm.meta.ix[adiix][dfworm.meta.ix[adiix].hostidx.isin(hostmda)].index.values
        mda_false = dfworm.meta.ix[adiix][~dfworm.meta.ix[adiix].hostidx.isin(hostmda)].index.values
        #mda affected fecundity
        mdayoung = mda_true[np.where(dfworm.meta.ix[mda_true].age < 6)]
        mdaold = mda_true[np.where(dfworm.meta.ix[mda_true].age >= 6)]
        #linear function defining fecundity during drug clearance
        mmda = float(fecund - 1) / (mda_clear - 1 )
        bmda = 1 - mmda * 1
        #new base fecundity under drugs
        sterile_t = (mmda * clear_count + bmda)
        dfworm.meta.ix[mdayoung, 'fec'] = np.random.poisson(sterile_t, mdayoung.shape[0])
       #linear function defining decline in fecundity with age
        m = float(0 - sterile_t) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        positive_lambda = (dfworm.meta.ix[mdaold].age.values * m) + b
        positive_lambda[positive_lambda < 0] = 0
        dfworm.meta.ix[mdaold, 'fec'] = np.random.poisson(positive_lambda).astype(np.int64)

        #normal fecunidty in individual with no mda
        young = mda_false[np.where(dfworm.meta.ix[mda_false].age < 6)]
        old = mda_false[np.where(dfworm.meta.ix[mda_false].age >= 6)]
        dfworm.meta.ix[young, 'fec'] = np.random.poisson(fecund, young.shape[0])
        #linear function defining decline in fecundity with age
        m = float(0 - fecund) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        positive_lambda = (dfworm.meta.ix[old].age.values * m) + b
        positive_lambda[positive_lambda < 0] = 0
        dfworm.meta.ix[old, 'fec'] = np.random.poisson(positive_lambda).astype(np.int64)
########################################################
    else: #base fecundity when no drugs
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
    if selection: #dfAdult.sel will be updated here to same length as dfAdult_mf.pos
        dfAdult_mf, dfworm = selection_fx(dfworm, dfAdult_mf, new_positions, cdslist)
    dfworm.add_worms(dfAdult_mf, new_positions)
    return(dfAdult_mf, dfworm)

def fecunditymda_sel1_fx(fecund,
                    dfworm,
                    locus,
                    mutation_rate,
                    recombination_rate,
                    basepairs,
                    selection,
                    densitydep_fec,
                    cdslist,
                    mda_sterile,
                    clear_count,
                    mda_clear,
                    hostmda):
    '''function for reduced fecundity under mda option 1
    option 1 simplifies that when no MDA or selective event all phenotypes
    are essetially wildtype, so fitness is not evaluated
    conditions: mda=True, selection=True, 1

    '''
############################################################
    if clear_count == 1: #permanent sterility
        adultmda = dfworm.meta[(dfworm.meta["stage"] == "A") &
                               (dfworm.meta.hostidx.isin(hostmda))].index.values
        adultsterile = dfworm.meta.ix[adultmda].groupby("hostidx").apply(lambda y:
                                y[np.random.random(y.shape[0]) < (mda_sterile ** y.fitF)].index.values)
        steriix = np.concatenate([i for i in adultsterile.values]).ravel()
        dfworm.meta.ix[steriix,"sex"] = "S"
    else: pass

    if clear_count > 0 and clear_count <= mda_clear: #Drugs cause temporary sterility over clear_time
        adiix = dfworm.meta[dfworm.meta.stage == "A"].index.values
        mda_true = dfworm.meta.ix[adiix][dfworm.meta.ix[adiix].hostidx.isin(hostmda)].index.values
        mda_false = dfworm.meta.ix[adiix][~dfworm.meta.ix[adiix].hostidx.isin(hostmda)].index.values
        #mda affected fecundity
        mdayoung = mda_true[np.where(dfworm.meta.ix[mda_true].age < 6)]
        mdaold = mda_true[np.where(dfworm.meta.ix[mda_true].age >= 6)]

        #linear function defining fecundity during drug clearance
        mmda = float(fecund - 1) / (mda_clear - 1 )
        bmda = 1 - mmda * 1
        #new base fecundity under drugs
        sterile_t = (mmda * clear_count + bmda)
        sterile_tarray = np.repeat(sterile_t, mdayoung.shape[0])
        dfworm.meta.ix[mdayoung, 'fec'] = np.random.poisson((sterile_tarray ** dfworm.meta.ix[mdayoung, "fitF"])).astype(np.int64)

       #linear function defining decline in fecundity with age and mda
        m = float(0 - sterile_t) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        positive_lambda = (dfworm.meta.ix[mdaold].age.values * m) + b
        positive_lambda[positive_lambda < 0] = 0
        dfworm.meta.ix[mdaold, 'fec'] = np.random.poisson((positive_lambda ** dfworm.meta.ix[mdaold,"fitF"])).astype(np.int64)

        #normal fecunidty in individual with no mda
        young = mda_false[np.where(dfworm.meta.ix[mda_false].age < 6)]
        old = mda_false[np.where(dfworm.meta.ix[mda_false].age >= 6)]
        dfworm.meta.ix[young, 'fec'] = np.random.poisson(fecund, young.shape[0])
        #linear function defining decline in fecundity with age
        m = float(0 - fecund) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        positive_lambda = (dfworm.meta.ix[old].age.values * m) + b
        positive_lambda[positive_lambda < 0] = 0
        dfworm.meta.ix[old, 'fec'] = np.random.poisson(positive_lambda).astype(np.int64)
###############################################################
    else: #base fecundity when no drugs
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
    if selection: #dfAdult.sel will be updated here to same length as dfAdult_mf.pos
        dfAdult_mf, dfworm = selection_fx(dfworm, dfAdult_mf, new_positions, cdslist)
    dfworm.add_worms(dfAdult_mf, new_positions)
    return(dfAdult_mf, dfworm)

def fecunditymda_sel2_fx(fecund,
                    dfworm,
                    locus,
                    mutation_rate,
                    recombination_rate,
                    basepairs,
                    selection,
                    densitydep_fec,
                    cdslist,
                    mda_sterile,
                    clear_count,
                    mda_clear,
                    hostmda):
    '''function for reduced fecundity under mda option 2
    option 2 is when the mutant are less fit then the wildtype when no mda
    is being applied.
    conditions: mda=True, selection=True, 2
    '''
    if clear_count == 1: #permanent sterility
        adultmda = dfworm.meta[(dfworm.meta["stage"] == "A") &
                               (dfworm.meta.hostidx.isin(hostmda))].index.values
        adultsterile = dfworm.meta.ix[adultmda].groupby("hostidx").apply(lambda y:
                                y[np.random.random(y.shape[0]) < (mda_sterile ** y.fitF)].index.values)
        steriix = np.concatenate([i for i in adultsterile.values]).ravel()
        dfworm.meta.ix[steriix,"sex"] = "S"
    else: pass

    if clear_count > 0 and clear_count <= mda_clear: #Drugs cause temporary sterility over clear_time
        adiix = dfworm.meta[dfworm.meta.stage == "A"].index.values
        mda_true = dfworm.meta.ix[adiix][dfworm.meta.ix[adiix].hostidx.isin(hostmda)].index.values
        mda_false = dfworm.meta.ix[adiix][~dfworm.meta.ix[adiix].hostidx.isin(hostmda)].index.values
        #mda affected fecundity
        mdayoung = mda_true[np.where(dfworm.meta.ix[mda_true].age < 6)]
        mdaold = mda_true[np.where(dfworm.meta.ix[mda_true].age >= 6)]

        #linear function defining fecundity during drug clearance
        mmda = float(fecund - 1) / (mda_clear - 1 )
        bmda = 1 - mmda * 1
        #new base fecundity under drugs
        sterile_t = (mmda * clear_count + bmda)
        sterile_tarray = np.repeat(sterile_t, mdayoung.shape[0])
        dfworm.meta.ix[mdayoung, 'fec'] = np.random.poisson((sterile_tarray ** dfworm.meta.ix[mdayoung, "fitF"])).astype(np.int64)

       #linear function defining decline in fecundity with age and mda
        m = float(0 - sterile_t) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        positive_lambda = (dfworm.meta.ix[mdaold].age.values * m) + b
        positive_lambda[positive_lambda < 0] = 0
        dfworm.meta.ix[mdaold, 'fec'] = np.random.poisson((positive_lambda ** dfworm.meta.ix[mdaold,"fitF"])).astype(np.int64)

        #normal fecunidty in individual with no mda
        young = mda_false[np.where(dfworm.meta.ix[mda_false].age < 6)]
        old = mda_false[np.where(dfworm.meta.ix[mda_false].age >= 6)]
        dfworm.meta.ix[young, 'fec'] = np.random.poisson(fecund, young.shape[0])
        #linear function defining decline in fecundity with age
        m = float(0 - fecund) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        positive_lambda = (dfworm.meta.ix[old].age.values * m) + b
        positive_lambda[positive_lambda < 0] = 0
        dfworm.meta.ix[old, 'fec'] = np.random.poisson(positive_lambda).astype(np.int64)
############################################
    else: #base fecundity when no drugs
        adiix = dfworm.meta[dfworm.meta.stage == "A"].index.values
        young = adiix[np.where(dfworm.meta.ix[adiix].age < 6)]
        old = adiix[np.where(dfworm.meta.ix[adiix].age >= 6)]
        fecund_tarray = np.repeat(fecund, young.shape[0])
        dfworm.meta.ix[young, 'fec'] = np.random.poisson((fecund_tarray ** (1 - abs(1 - dfworm.meta.ix[young, "fitF"])))).astype(np.int64)
        #linear function defining decline in fecundity with age
        m = float(0 - fecund) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        positive_lambda = (dfworm.meta.ix[old].age.values * m) + b
        positive_lambda[positive_lambda < 0] = 0
        dfworm.meta.ix[old, 'fec'] = np.random.poisson((positive_lambda ** 1 - abs(1 - dfworm.meta.ix[old,"fitF"]))).astype(np.int64)
############################################
    #sex, recombination, mutation
    dfAdult_mf = recombination_fx(locus, dfworm, adiix, recombination_rate, basepairs)
    dfAdult_mf.meta.sex = np.random.choice(['M', 'F'] , size=len(dfAdult_mf.meta))
    dfAdult_mf.meta.age = 1
    dfAdult_mf, new_positions = mutation_fx(locus, dfAdult_mf,
         mutation_rate, recombination_rate, basepairs)
    if selection: #dfAdult.sel will be updated here to same length as dfAdult_mf.pos
        dfAdult_mf, dfworm = selection_fx(dfworm, dfAdult_mf, new_positions, cdslist)
    dfworm.add_worms(dfAdult_mf, new_positions)
    return(dfAdult_mf, dfworm)
