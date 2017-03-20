#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
#import ipdb
import math
import subprocess
import numpy as np
import pandas as pd
import random
import pickle
import copy

from .agehost import agehost_fx
from .worm import Worms

def host_fx(village, infhost):
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
    deathdict = pickle.load( open( "../figs/data/acttable.p", "rb" ) )
    assert len(village) == len(infhost)
    coordinates = []
    host_idx = []

    for vill in range(len(village)):
         #list of host positions
         coordinates.extend( np.random.negative_binomial(village[vill].size, village[vill].size
                             / float((village[vill].size+village[vill].mu)), (infhost[vill],
                                                           2)) + village[vill].dist)

         for host in range(infhost[vill]):
             host_idx.append("v" + str(vill) + "h" + str(host + 1))
    sex = [random.choice("01") for i in range(sum(infhost))]
    age_death = [agehost_fx(i, deathdict) for i in sex]

    dfHost = pd.DataFrame({
        'village': np.repeat(range(len(village)), infhost),
        'hostidx': host_idx,
        'sex': sex,
        'age': [i[0] for i in age_death],
        'agedeath': [i[1] for i in age_death],
        'coordinates': coordinates,
        'MDAstate': np.zeros(sum(infhost)),
        'MDAcum': [[] for i in range(sum(infhost))]
    })
    dfHost = dfHost.loc[:, ['village', 'hostidx', 'sex',
            'age', 'agedeath', 'coordinates', 'MDAstate', 'MDAcum']]
    return(dfHost)


def coalsims_migmat_fx(numvillages, initial_migration, initial_distance_m, thetaN0,
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

    if numvillages > 4:
        raise ValueError("only handles 4 villages ATM")
    elif numvillages <= 4:
        assert len(initial_distance_m) == ((numvillages * (numvillages - 1)) / 2)
        mig = []  # initiate blank migration list
        for meters in initial_distance_m:
            mig.append((initial_migration) / (np.random.exponential(meters)))
        if numvillages == 2:
            m1 = 4 * ne * mig[0]  # 4Nm
            return("{}".format(m1))  # mig_matrix is symmetrical and island
        elif numvillages == 3:
            m1 = 4 * ne * mig[0]
            m2 = 4 * ne * mig[1]
            m3 = 4 * ne * mig[2]
            # mig_matrix is symmetrical
            return("{} {} {} {} {} {} {} {} {}".format(
                0, m1, m2, m1, 0, m3, m2, m3, 0))
        elif numvillages == 4:
            m1 = 4 * ne * mig[0]
            m2 = 4 * ne * mig[1]
            m3 = 4 * ne * mig[2]
            m4 = 4 * ne * mig[3]
            m5 = 4 * ne * mig[4]
            m6 = 4 * ne * mig[5]
            return("{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(
                0, m1, m2, m3, m1, 0, m4, m5, m2, m4, 0, m6, m3, m5, m6, 0))


def parse_coalsims_fx(msout, ploidy, nind):
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
         columns to be added to dfworm
    pos
         list of positions for the mutaitons
    """
    print('NINDVS ' + str(nind))
    nind = int(nind)
    # Parsing positions
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
        gt_array = np.zeros((nind, pos.shape[0]) , dtype=np.uint8)
        #:TODO Seriously think about returning the whole genotype array
        cix = 0
        for line in iter(msout.stdout.readline, ''):
            line = line.decode('utf-8')
            line = list(line.strip())
            try:
                gt_array[cix, :] = np.array(line, dtype=np.uint8)
            except IndexError:
                break
            cix += 1
        return(gt_array, pos)
    elif ploidy == 2:
        cix = 0
        gt_array = np.zeros((nind, pos.shape[0]), dtype=np.uint8)
        gt_array2 = np.zeros((nind, pos.shape[0]) , dtype=np.uint8)
        for line in iter(msout.stdout.readline, ''):
            line = line.decode('utf-8')
            hap = np.array(list(line.strip()), dtype=np.uint8)
            gt_array[cix, :] = hap
            hap2 = next(iter(msout.stdout.readline, ''))
            gt_array2[cix, :] =  np.array(list(hap2.strip()),
                    dtype=np.uint8)
            cix += 1
            if cix == nind:
                break
        return(gt_array, gt_array2, pos)


def coalsims_fx(worm_popsize, numvillages, initial_migration, initial_distance_m,
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
        worm_popsize = [x*2 for x in worm_popsize]
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
    if numvillages == 1:
        scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
                     "-G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
                     " {sig_digits} ")
        mscmd = scrm_base.format(**ms_params)
    else:  # ms setup for >1 villages
        num_subpops = len(worm_popsize)  # -I num_pops
        sub_pop = " ".join(map(str, worm_popsize))  # -I X i j ...
        mm = coalsims_migmat_fx(
           numvillages,
           initial_migration,
           initial_distance_m,
           thetaN0,
           basepairs,
           mutation_rate)
        if numvillages == 2:
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

            for village_ix in range(numvillages):
                present_pop = float(theta[village_ix]) / thetaN0
                subpopstr += '-n {} {} '.format(village_ix + 1, present_pop)
                if village_ix != numvillages - 1:
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
        gt, positions = parse_coalsims_fx(msout, ploidy, sum(worm_popsize))
        return(gt, positions)
    elif ploidy == 2:
        gt, gt2, positions = parse_coalsims_fx(msout, ploidy,
                sum(worm_popsize)/2)
        return(gt, gt2, positions)


def sel_fx(dfworm, basepairs, cdslist):
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
#    #size parameter
#    size = 3
#    #average distance between in bp
#    mu = intgen_length
#    #last coordinate is approx num_cds * mu; so if num_cds is too large or mu is too long
#    #genes will run over the locus length
#    #dfworm.pos
#    for locus in dfworm.h2.keys():
#        num_cds = int(round((perc_locus[int(locus)]*basepairs[int(locus)]) /
#            cds_length))
#        size_cds = np.round(np.random.gamma(4, 0.25, num_cds) *
#                cds_length)
#
#         #clustered like operons in c elegans
##        cds_F = range(0, int(round(num_cds / 2.0)) + 1)
##        cds_S = range(cds_F.pop(), num_cds)
#
#        #random
#        rand_cds = np.arange(num_cds)
#        np.random.shuffle(rand_cds)
#        cds_F = rand_cds[0:int(round(num_cds / 2.0))]
#        cds_S = rand_cds[len(cds_F): num_cds]
#
#        dfworm.sel[locus + "Ft"] = sum(size_cds[cds_F])
#        dfworm.sel[locus + "St"] = sum(size_cds[cds_S])
#        #r = size
#        #m = mean
#        #p = r / (  r + m )
#        cds_between = np.random.negative_binomial(size, size/float(mu+size), num_cds)
#        cds_stop = 0
#        cds_coords = []
#
#        for i, j in zip(map(np.int, cds_between), map(np.int,size_cds)):
#            #[i + cds_stop, i + cds_stop + j]
#            if (i + cds_stop > basepairs[int(locus)]) or (i + j + cds_stop > basepairs[int(locus)]):
#                break
#            else:
#                cds_coords.append([i + cds_stop, i + j + cds_stop])
#                cds_stop += (i + j)
        #[[1000,600],[]]

    surv_length = cdslist[0]
    surv_position = cdslist[1]
    fec_length = cdslist[2]
    fec_position = cdslist[3]

    for loc in range(len(dfworm.h2.keys())):
        bp = basepairs[loc + 1]
        dfworm.coord[str(loc + 1) + "F"] = [(bp * fec_position[loc][i], bp
                    * fec_position[loc][i] + pos) for i, pos in enumerate(fec_length[loc])]
        dfworm.coord[str(loc + 1) + "S"] = [(bp * surv_position[loc][i], bp
                    * surv_position[loc][i] + pos) for i, pos in enumerate(surv_length[loc])]
        selS = []
        selF = []
        for position in dfworm.pos[str(loc + 1)]:
            if any([i <= position <= j for i, j in dfworm.coord[str(loc + 1) + "F"]]):
                 #shape = 4, mean = 1, scale = mean/shape
                 #here mean is mean_fitness, wildtype is assumed to be 1
                 selF.append(np.random.gamma(4, scale=0.25))
                 selS.append(0)
            elif any([i <= position <= j for i,j in dfworm.coord[str(loc + 1) + "S"]]):
                     selS.append(np.random.gamma(4, scale=0.25))
                     selF.append(0)
            else: #not in a cds
                selS.append(0)
                selF.append(0)
        dfworm.sel[str(loc + 1) + "F"] = np.array(selS)
        dfworm.sel[str(loc + 1) + "S"] = np.array(selF)
        dfworm.sel[str(loc + 1) + "Ft"] = sum(fec_length[loc])
        dfworm.sel[str(loc + 1) + "St"] = sum(surv_length[loc])
    return(dfworm)


def fit_fx(dfworm, cdslist):
    ''' Calculates mean fitness for each individual by summing fitness effects
    from dfSel for each position across all loci. Fitness here is 1 + s, where s
    is found by fit - 1.

    Parameters
    ----------
    dfworm : df
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
    ninds = len(dfworm.meta)
    dom = cdslist[4]
    fitF_ind = np.zeros(ninds)
    fitS_ind = np.zeros(ninds)
    Floc = 0
    Sloc = 0
    for locus in dfworm.h2.keys():
        count_sites = dfworm.h1[locus] + dfworm.h2[locus]
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

        cds_sites_S = np.dot(dfworm.h1[locus], intsites_S) \
            + np.dot(dfworm.h2[locus], intsites_S)
        cds_sites_F = np.dot(dfworm.h1[locus], intsites_F) \
            + np.dot(dfworm.h2[locus], intsites_F)
        fitS_ind += (( (dfworm.sel[locus + "St"] * 2) - cds_sites_S) + sum_selsites_S) / (dfworm.sel[locus + "St"] * 2)
        fitF_ind += (( (dfworm.sel[locus + "Ft"] * 2) - cds_sites_F) + sum_selsites_F) / (dfworm.sel[locus + "Ft"] * 2)
        Floc += 1
        Sloc += 1
    return(fitF_ind / Floc, fitS_ind / Sloc)


def wormdf_fx(village, infhost, muWormBurden, sizeWormBurden, locus,
              initial_migration, initial_distance_m, theta, basepairs, mutation,
              recombination, time2Ancestral, thetaRegional, time_join, selection,
              cdslist):
     '''Function to generate df and coalescent histories
     Parameters
     ---------
     infhost : int, list
        number of infectious hosts per village, calculate by hostpopsize * prevalence
     Returns
     -------
     dfworm : df
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
     worm_popsize = [sum(i) for i in popinit]
     host_idx = []
     for vill in range(len(village)):
          for host in range(len(popinit[vill])):
               host_idx.extend(["v" + str(vill) + "h" + str(host + 1)] * popinit[vill][host])

     dfworm = pd.DataFrame({
                      'village' : np.repeat(range(len(village)), worm_popsize),
                      'hostidx' : host_idx,
                      'age' : 1,
                      'sex' : [random.choice("MF") for i in range(sum(worm_popsize))],
                      'stage' : "A",
                      'R0net' : np.random.random(sum(worm_popsize)),
                      'fec' : 0
                      })
     dfworm = dfworm.loc[:, ['village', 'hostidx', 'age',
            'sex', 'stage', 'R0net', 'fec']]
     # Add genetic data
     dfworm = Worms(dfworm)
     posSel = []
     for loc in range(locus):
         if recombination[loc] == 0:
             gt_array, mutations = coalsims_fx(worm_popsize, len(village),
                     initial_migration, initial_distance_m,
                     theta[loc], basepairs[loc], mutation[loc],
                     recombination[loc], time2Ancestral, thetaRegional,
                     time_join)
             dfworm.h1[str(loc)] = gt_array
             dfworm.pos[str(loc)] = mutations
         elif recombination[loc] > 0:
             gt_array, gt_array2, mutations = coalsims_fx(worm_popsize,
                     len(village), initial_migration, initial_distance_m,
                     theta[loc], basepairs[loc], mutation[loc],
                     recombination[loc], time2Ancestral, thetaRegional, time_join)
             dfworm.h1[str(loc)] = gt_array
             dfworm.h2[str(loc)] = gt_array2
             dfworm.pos[str(loc)] = mutations
             posSel.append(mutations)
     # Create dfSel
     if selection:
         dfworm = sel_fx(dfworm, basepairs, cdslist)
         fitS, fitF = fit_fx(dfworm, cdslist)
         dfworm.meta["fitF"] = fitF
         dfworm.meta["fitS"] = fitS
     return(dfworm)


def wbsims_init(village, hostpopsize, prevalence, muWormBurden,
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
    dfworm : df
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
    dfHost = host_fx(village, infhost)

    dfworm = wormdf_fx(village, infhost, muWormBurden, sizeWormBurden,
                              locus, initial_migration, initial_distance_m, theta,
                              basepairs, mutation_rate, recombination_rate, time2Ancestral, thetaRegional,
                              time_join, selection, cdslist)
#    ipdb.set_trace()
    return(dfHost, dfworm)

####3 loci, 2 villages, with selection
#if __name__ == '__main__':
#    village, dfHost, dfworm, dfJuv, dfMF=wbsims_init(2, [100, 200],
#    [0.1, 0.3], 100, 1, [5, 5], [50, 50], 3, 0.0001, [1000], [[5, 5], [10, 10],[10, 10]],
#    [13000, 200000, 100000],[7.6E-8, 2.9E-9, 2.9E-9], [0, 2.9E-9, 2.9E-9], 1800, 23, 240,
#    True, [[0, 0.18, 0.20], 1100, 2500])
