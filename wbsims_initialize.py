#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""

import math
import subprocess
import numpy as np
import pandas as pd
import random
from agehost import agehost_fx
from filtercoords import filtercoords_fx
import pickle

def host_fx(villages, infhost, muTrans, sizeTrans):
    '''Creates a transmission matrix for locations of infected hosts
       default_test : host_fx(2, [100, 300], 100, 1)

    Parameters
    ----------
    villages : int
        number of villages in the simulation
    infhost : list of pop. sizes
         number of infected hosts at intial time. prev*hostpopsize
    muTrans : float
        average distance in meters between hosts
    sizeTrans : float
        aggregation parameter

    Returns
    -------
    dfHost : dataframe
    '''
    deathdict = pickle.load( open( "./acttable.p", "rb" ) )
    assert villages == len(infhost)
    coordinates = []
    host_idx = []
    for vill in range(villages):
         #list of host positions
         coordinates.extend(np.random.negative_binomial(sizeTrans, sizeTrans
                             / float((sizeTrans+muTrans)), (infhost[vill],
                                                           2)))
         for host in range(infhost[vill]):
             host_idx.append("v" + str(vill) + "h" + str(host + 1))
    sex = [random.choice("01") for i in range(sum(infhost))]
    age_death = [agehost_fx(i, deathdict) for i in sex]

    dfHost = pd.DataFrame({
        'village': np.repeat(range(villages), infhost),
        'hostidx': host_idx,
        'sex': sex,
        'age': [i[0] for i in age_death],
        'agedeath': [i[1] for i in age_death],
        'coordinates': coordinates,
        'MDA': np.zeros(sum(infhost)),
        'MDA_cum': np.zeros(sum(infhost))
    })
    dfHost = dfHost.loc[:, ['village', 'hostidx', 'sex',
            'age', 'agedeath', 'coordinates', 'MDA', 'MDA_cum']]
    return(dfHost)


def coalsims_migmat_fx(villages, initial_migration, initial_distance_m, thetaN0,
                     basepairs, mutation_rate):
    '''Creates a string that represents a migration matrix between
    metapopulations. Migration here is a stepping stone not island model

    Format specified by ms (hudson 2000). Uses euclidian distances. The value
    of 2Nm is weighted by the distance from the next population as an
    exponential random variable.  The highest this can be is 2Nm/1

    Parameters
    ----------
    thetaN0 : float
         theta from the first village of the locus

    Returns
    -------
    Migration matrix : string
        migration matrix for scrm
    '''
    ne = thetaN0 / (4 * mutation_rate * basepairs)

    if villages > 4:
        raise ValueError("only handles 4 villages ATM")
    elif villages < 4:
        if len(initial_distance_m) != ((villages) * (villages - 1) / 2):
            raise ValueError(("there are not adequate pairwise comparisons in"
                              "distance_m to match villages"))
        mig = []  # initiate blank migration list
        for meters in initial_distance_m:
            mig.append((initial_migration) / (np.random.exponential(meters)))
        if villages == 2:
            m1 = 4 * ne * mig[0]  # 4Nm
            return("{}".format(m1))  # mig_matrix is symmetrical and island
        elif villages == 3:
            m1 = 4 * ne * mig[0]
            m2 = 4 * ne * mig[1]
            m3 = 4 * ne * mig[2]
            # mig_matrix is symmetrical
            return("{} {} {} {} {} {} {} {} {}".format(
                0, m1, m2, m1, 0, m3, m2, m3, 0))
        elif villages == 4:
            m1 = 4 * ne * mig[0]
            m2 = 4 * ne * mig[1]
            m3 = 4 * ne * mig[2]
            m4 = 4 * ne * mig[3]
            m5 = 4 * ne * mig[4]
            m6 = 4 * ne * mig[5]
            return("{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(
                0, m1, m2, m3, m1, 0, m4, m5, m2, m4, 0, m6, m3, m5, m6, 0))


def parse_coalsims_fx(msout, ploidy):
    """Parse output from ms or scrm.
    Parameters
    ---------
    msout : iterator,file
         stdout from scrm
    ploidy : int
         ploidy of the current locus

    Returns
    ------
    gt_array
         columns to be added to dfAdult
    pos
         list of positions for the mutaitons
    """
    for line in iter(msout.stdout.readline, ''):
        line = line.decode('utf-8')
        if line.startswith("positions"):
            # collisions can result here when theta is high
            pos = np.round(
                np.array(
                    line.strip().split()[1:],
                    dtype=np.float64))
            prev = 0
            # :TODO double check if this works correctly
            for idx, item in enumerate(pos, start=0):
                while prev >= item:
                    item += 1
                pos[idx] = item
                prev = pos[idx]
            break
        else:
            pass
    pos = pos.astype(np.int64)
    if ploidy == 1:
        gt_array = []
        #:TODO Seriously think about returning the whole genotype array
        for line in iter(msout.stdout.readline, ''):
            line = line.decode('utf-8')
            hap = np.array(list(line.strip()), dtype=int)
            try:
                gt = hap * pos
                gt_array.append(gt[gt != 0])
            except ValueError:
                # Last line
                break
        #print(gt_array)
        return(gt_array, pos)
    elif ploidy == 2:
        gt_array = []
        gt_array2 = []
        for line in iter(msout.stdout.readline, ''):
            line = line.decode('utf-8')
            hap = np.array(list(line.strip()), dtype=int)
            try:
                gt = hap * pos
                gt_array.append(gt[gt != 0])
                hap2_temp = next(iter(msout.stdout.readline, ''))
                hap2 = np.array(list(hap2_temp.strip()), dtype=int)
                gt2 = hap2 * pos
                gt_array2.append(gt2[gt2 != 0])
            except ValueError:
                # End of output
                break
        return(gt_array, gt_array2, pos)


def coalsims_fx(worm_popsize, villages, initial_migration, initial_distance_m,
        theta, basepairs, mutation_rate, recombination_rate, time2Ancestral, thetaRegional,
        time_join):
    '''External call to ms (Hudson 2000) or scrm.

    Calls the function migration_matrix() if using more than 1 village.
    This function will then read in the stdout from ms/scrm. The
    location of segsites are assigned a random number then scaled
    by the length of the desired sequence in basepairs.

    Parameters
    ----------
    worm_popsize: (list, int)
        how many worms to simulate
    villages: int
        number of villages/metapopulations
    theta: (list,float)
        list of theta values for each locus. With more than 1 village
        it is list of lists [[locus1_meta1,locus1_meta2],
        [locus2_meta1,locus2_meta2]]
    basepairs:(list,int)
        list of basepairs (lengths) for each locus
    mutation: (list,float)
        mutation rate for each locus as probability per base per generation
    recombination:(list,float)
        recombination rate for each locus as probability per base
        per generation
    thetaRegional: float
        theta for regional pops to the ratio of N0 e.g.,23 times larger
    time2Ancestral: int
        time in generations to the ancestral population
    time_join: int
        time in generations for joining/splitting villages

    Returns
    -------
    gt : array
         genotype data for first haplotype
    gt2 : array
         genotype data for second haplotype
    mutations : array
         array of mutation positions
    '''
    # theta for first village
    thetaN0 = theta[0]
    #recombination rate for locus
    rho = thetaN0 * (recombination_rate / mutation_rate)
    rho is 0
    # time to the ancestral population
    tA = basepairs * (time2Ancestral / (thetaN0 / mutation_rate))
    # time to join villages
    tjoin = basepairs * (time_join / (thetaN0 / mutation_rate))
    #check ploidy for total haplotypes
    if rho == 0:
        ploidy = 1
    else:
        worm_popsize = [x * 2 for x in worm_popsize]
        ploidy = 2
    # parameters for ms or scrm
    ms_params = {
        'nhaps': sum(worm_popsize),
        'theta': thetaN0,
        'rho': rho,
        'basepair': basepairs - 1,
        'exp_growth': (-1 / tA) * math.log(thetaRegional),
        'time_growth': tA,
        'sig_digits': 12,
    }

    # order matters for the command call
    if villages == 1:
        scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
                     "-G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                     " {sig_digits} ")
        mscmd = scrm_base.format(**ms_params)
    else:  # ms setup for >1 villages
        num_subpops = len(worm_popsize)  # -I num_pops
        sub_pop = " ".join(map(str, worm_popsize))  # -I X i j ...
        mm = coalsims_migmat_fx(
           villages,
           initial_migration,
           initial_distance_m,
           thetaN0,
           basepairs,
           mutation_rate)
        if villages == 2:
            scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
                         "{sub_pop} {present_pop} {join} -G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                         " {sig_digits} ")
            ms_params[
                'present_pop'] = '-n 1 1 -n 2 {}'.format(float(theta[1]) / thetaN0)
            ms_params[
                'sub_pop'] = '-I {} {} {}'.format(num_subpops, sub_pop, mm)
            ms_params['join'] = '-ej {} 1 2'.format(tjoin)
            mscmd = scrm_base.format(**ms_params)
        else:
            scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
                         "{sub_pop} {present_pop} {ma} {join} -G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                         " {sig_digits} ")
            ms_params['sub_pop'] = '-I {} {}'.format(num_subpops, sub_pop)
            ms_params['ma'] = '-ma {} '.format(mm)
            joinstr = ''
            subpopstr = ''

            for village_ix in range(villages):
                present_pop = float(theta[village_ix]) / thetaN0
                subpopstr += '-n {} {} '.format(village_ix + 1, present_pop)
                if village_ix != villages - 1:
                    joinstr += '-ej {0} {1} {2} '.format(
                        tjoin,
                        village_ix + 1,
                        village_ix + 2)
            ms_params['present_pop'] = subpopstr
            ms_params['join'] = joinstr
            mscmd = scrm_base.format(**ms_params)
    print(mscmd)
    msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
    #msout = subprocess.check_output(mscmd.split(" ")).decode('utf-8')
    #:TODO refactor this, less variation on returns
    if ploidy == 1:
        gt, mutations = parse_coalsims_fx(msout, ploidy)
        return(gt, mutations)
    elif ploidy == 2:
        gt, gt2, mutations = parse_coalsims_fx(msout, ploidy)
        return(gt, gt2, mutations)


def sel_fx(locus, positions, basepairs, perc_locus, cds_length, intgen_length):
    '''Initializes the distribution of fitness effects for each mutation

    Parameters
    ---------
    positions : list
      list from coalsims_fx, the the line "positions" in the scrm/ms output

    Returns
    ------
    dfSel : df
        df containing mutations and fitness effect
    cds_coordinates : nested array
        coordinates (start, end) of coding sequences

    '''
    cds_positions, cds_coordinates = filtercoords_fx(positions, basepairs, perc_locus,
                                        cds_length, intgen_length)
    selF =[]
    selS =[]
    #cds_positions = [item for sublist in cds_positions for item in sublist]
    #below functs assume the list of cds_positions contains mutliple loci
    if locus > 2:
        numpos = [len(i) for i in cds_positions]
        pos = 0
        while pos < sum(numpos):
            if random.choice("SF") is "F":
                 #shape = 4, mean = 1, scale = mean/shape
                 #here mean is mean_fitness, wildtype is assumed to be 1
                 selF.append(np.random.gamma(4, scale=0.25))
                 selS.append(1)
            else:
                 selS.append(np.random.gamma(4, scale=0.25))
                 selF.append(1)
            pos += 1
        dfSel = pd.DataFrame({
                            'locus' : np.repeat(range(1, locus), numpos),
                            'position' : [item for sub in cds_positions for item in sub],
                            'selF' : selF,
                            'selS' : selS
                            })
        dfSel = dfSel.loc[:, ['locus', 'position', 'selF',
             'selS']]
    else: #list only contains a single locus
        numpos = len(cds_positions[0])
        pos = 0
        while pos < numpos:
            if random.choice("SF") is "F":
                #shape = 4, mean = 1, scale = mean/shape
                #here mean is mean_fitness, wildtype is assumed to be 1
                selF.append(np.random.gamma(4, scale=0.25))
                selS.append(1)
            else:
                selS.append(np.random.gamma(4, scale=0.25))
                selF.append(1)
            pos += 1
        dfSel = pd.DataFrame({
                            'locus' : np.repeat(range(1, locus), numpos),
                            'position' : cds_positions[0],
                            'selF' : selF,
                            'selS' : selS,
                             })
        dfSel = dfSel.loc[:, ['locus', 'position', 'selF',
             'selS']]
    return(dfSel, cds_coordinates)


def fit_fx(locus, dfAdult, dfSel):
     ''' Calculates mean fitness for each individual by summing fitness effects
     from dfSel for each position across all loci

     Parameters
     ----------
     dfAdult : df
          data frame of adult worms containing genotype information
     dfSel : df
          data fram of fitness benefit for each allele

     Returns
     -------
     fitF : array
          array filling selF column for fecundity fitness
     fitS : array
          array filling selS column for survival fitness

     '''
     ##fitness of individual in dfAdult from values in dfSel
     fitS_ind = []
     fitF_ind = []
     fitS = []
     fitF = []
     for index, row in dfAdult.iterrows():
          for loc in range(1, locus):
               fitS_ind.extend(dfSel.loc[dfSel["position"].isin(row.ix
                                         ["locus_" + str(loc) + "_h1"])]
                                         ['selS'][dfSel["locus"] == loc])

               fitF_ind.extend(dfSel.loc[dfSel["position"].isin(row.ix
                                         ["locus_" + str(loc) + "_h1"])]
                                         ['selF'][dfSel["locus"] == loc])
          fitS.append(round(np.mean(fitS_ind), 5))
          fitF.append(round(np.mean(fitF_ind), 5))

     return(fitS, fitF)

def wormdf_fx(villages, infhost, muWormBurden, sizeWormBurden, locus,
              initial_migration, initial_distance_m, theta, basepairs, mutation,
              recombination, time2Ancestral, thetaRegional, time_join, selection,
              perc_locus, cds_length, intgen_length):
     '''Function to generate df and coalescent histories
     Parameters
     ---------
     infhost : int, list
        number of infectious hosts per village, calculate by hostpopsize * prevalence
     Returns
     -------
     dfAdult : df
         dataframe of adult wb worms, if selection is True contains the columns
         selF and selS, denoting average fitness for fecundity and selection
     dfSel : df
         if selection is True, returns list of mutation and fitness effects
     cds_coordinates : int, array
         if selection is True, returns list of coding seq locations
     '''
     #create parasite burden per host
     popinit = []
     for mu, size, numworms in zip(muWormBurden, sizeWormBurden, infhost):
         #the number of worms per infected host
         wb_burden = np.random.negative_binomial(size, size/float(mu+size), numworms)
         #[[10,12,1,15 ... infhostvill1],[19,1,5,4, ... infhostvill2]]
         popinit.append(np.array(wb_burden).tolist())

     #total worms per villages
     wormpopsize = [sum(i) for i in popinit]
     host_idx = []
     for vill in range(villages):
          for host in range(len(popinit[vill])):
               host_idx.extend(["v" + str(vill) + "h" + str(host + 1)] * popinit[vill][host])

     dfAdult = pd.DataFrame({
                      'village' : np.repeat(range(villages), wormpopsize),
                      'hostidx' : host_idx,
                      'age' : 1,
                      'sex' : [random.choice("MF") for i in range(sum(wormpopsize))],
                      'R0net' : np.random.random(sum(wormpopsize)),
                      'fec' : 0
                      })
     dfAdult = dfAdult.loc[:, ['village', 'hostidx', 'age',
            'sex', 'R0net', 'fec']]
     # Add genetic data
     posSel = []
     for loc in range(locus):
         if recombination[loc] == 0:
             gt_array, mutations = coalsims_fx(wormpopsize, villages, 
                     initial_migration, initial_distance_m, 
                     theta[loc], basepairs[loc], mutation[loc],
                     recombination[loc], time2Ancestral, thetaRegional, 
                     time_join)
             dfAdult["locus_" + str(loc)] = gt_array
         elif recombination[loc] > 0:
             gt_array, gt_array2, mutations = coalsims_fx(wormpopsize, 
                     villages, initial_migration, initial_distance_m, 
                     theta[loc], basepairs[loc], mutation[loc],
                     recombination[loc], time2Ancestral, thetaRegional, time_join)
             dfAdult["locus_" + str(loc) + "_h1"] = gt_array
             dfAdult["locus_" + str(loc) + "_h2"] = gt_array2
             posSel.append(mutations)
     # Create dfSel
     if selection:
          dfSel, cds_coordinates = sel_fx(locus, posSel, basepairs, 
                  perc_locus, cds_length, intgen_length)
          fitS, fitF = fit_fx(locus, dfAdult, dfSel)
          dfAdult["fitF"] = fitF
          dfAdult["fitS"] = fitS
          return(dfAdult, dfSel, cds_coordinates)
     else:
          return(dfAdult)


def wbsims_init(villages, hostpopsize, prevalence, muTrans, sizeTrans, muWormBurden,
                sizeWormBurden, locus, initial_migration, initial_distance_m, theta,
                basepairs, mutation_rate, recombination_rate, time2Ancestral, thetaRegional,
                time_join, selection, cdslist):
    '''This function initializes all dataframes for wbsims.py by using coalescent
    simulations to create a population genetic history.

    Parameters
    ---------
    villages : int
      number of villages to simulate
    hostpopsize : int, list
      maximum population size of each village
    prevalence : float, list
      prevalence of disease in each village
    muTrans : int
     average distance between hosts in village
    sizeTrans : int
     aggregation of hosts in village
    muWormBurden : int, list
     average number of worms per host
    sizeWormBurden : int, list
     aggregation of worms in hosts
    locus : int
     number of loci to simulate, if 1 is always assumed to be haploid
    initial_migration : float
     migration for ms/scrm as the fraction of subpopulation i
     which is made up of migrants from subpopulation j each generation
    initial_distance_m : int, list
     distances between pairs of villages in meters
    theta : float, list
     population mutation rate of each locus in each village
    basepairs : int, list
     length of each locus
    mutation_rate : float, list
     mutation rate per basepair for each locus
    recombination_rate : float, list
     recombination rate between basepairs for each locus, a value of 0 is
     no recombination
    time2Ancestral : int
     time in generation to the ancestral population
    thetaRegional : int
     size of the ancestral population compared to current population in theta
    time_join : int
     time of split between villages
    selection : boolean
     use selection in the model, runs sel_fx, and fit_fx to return dfSel
    perc_locus  : float, list
     percent of locus that is coding
    cds_length : int
     average length of coding regions
    intgen_length : int
     distance between coding regions

    Returns
    ------
    dfAdult : df
      initialize df for adult worms
    dfHost : df
      initialize df for host
    dfSel :df
      initialized df for selection
    dfJuv : df
      emtpy and initialized for main functions
    dfMF : df
      empty and intialized for main functions

    Example of parameters 2 villages with selection:
    villages = 2
    hostpopsize = [100, 200]
    prevalence = [0.1, 0.3]
    muTrans = 100
    sizeTrans = 1
    muWormBurden = [5, 5]
    sizeWormBurden = [50, 50]
    locus = 2
    initial_migration = 0.0001
    initial_distance_between = [1000]
    theta = [[5, 5], [1, 1]]
    basepairs = [13000, 200000]
    mutation_rate = [7.6E-8, 2.9E-9]
    recombination_rate = [0, 2.9E-9]
    time2Ancestral = 1800
    thetaRegional = 23
    time_join = 240
    selection = True
    perc_locus = [0, 0.18]
    cds_length = 1100
    intgen_length = 2500
    '''
    hostpopsize = np.array(hostpopsize)
    prevalence = np.array(prevalence)
    infhost = np.round(hostpopsize * prevalence).astype(np.int64)
    dfHost = host_fx(villages, infhost, muTrans, sizeTrans)
    perc_locus = cdslist[0]
    cds_length = cdslist[1]
    intgen_length = cdslist[2]
    if selection:
         dfAdult, dfSel, cds_coordinates= wormdf_fx(villages, infhost, muWormBurden, sizeWormBurden,
                              locus, initial_migration, initial_distance_m, theta,
                              basepairs, mutation_rate, recombination_rate, time2Ancestral, thetaRegional,
                              time_join, selection, perc_locus, cds_length, intgen_length)
    else: #no selection
        dfAdult = wormdf_fx(villages, infhost, muWormBurden, sizeWormBurden,
                    locus, initial_migration, initial_distance_m, theta,
                    basepairs, mutation_rate, recombination_rate, time2Ancestral, thetaRegional,
                    time_join, selection, perc_locus, cds_length, intgen_length)

        dfSel = pd.DataFrame({
            'locus' : [],
            'position' : [],
            'selF' : [],
            'selS' : [],
            })

        cds_coordinates = []

    dfJuv = pd.DataFrame({}, columns = dfAdult.columns)
    dfMF = pd.DataFrame({}, columns = dfAdult.columns)
    return(dfHost, dfAdult, dfJuv, dfMF, dfSel, cds_coordinates)

##3 loci, 2 villages, with selection
#if __name__ == '__main__':
#    dfHost, dfAdult, dfJuv, dfMF, dfSel, cds_coordinates=wbsims_init(2, [100, 200],
#    [0.1, 0.3], 100, 1, [5, 5], [50, 50], 3, 0.0001, [1000], [[5, 5], [10, 10],[10, 10]],
#    [13000, 200000, 100000],[7.6E-8, 2.9E-9, 2.9E-9], [0, 2.9E-9, 2.9E-9], 1800, 23, 240,
#    True, [[0, 0.18, 0.20], 1100, 2500])
#    from IPython import embed
#    embed()
