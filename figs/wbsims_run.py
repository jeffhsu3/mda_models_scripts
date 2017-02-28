#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
#import python modules
import random
import argparse
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
import pandas as pd
#import program functions
import wbsims_initialize
import transmission
import survival
import record

def get_args():
    '''setup configparser()'''
    parser = argparse.ArgumentParser()
    ## initialize
    #migration_matrix
    parser.add_argument('-im', '--initial_migration', type=float, default=.0001, 
            help=("migration rate between villages/metapopulations in the" 
            "model.This is strictly for initial conditions and can be" 
            "changed for the forward-in-time portion"))
    parser.add_argument('-idm', '--initial_distance_m', type=list, 
            help="initial_distance_m is list [1000] such that distance_m[0] is between 1 & 2")
    parser.add_argument('-v', '--villages', type=int, default=1, 
            help="sets the intial number of villages.")
    parser.add_argument('-t', '--theta', type=list, required=True, 
            help=("observed theta value of worm populations for" 
            "each locus in format [[meta1_locus1,meta2_locus1],"
            "[meta1_locus2, meta2_locus2]]"))
    parser.add_argument('-bp', '--basepairs', type=list, default=13000, 
            help="length in basepairs of each locus")
    parser.add_argument('-u', '--mutation', type=list, default=7.6E-8, 
            help="expected as list, mutation rate per bp per generation for each locus")
    #ms_outcall
    parser.add_argument('-t12', '--time_join12', type=int, default=240, 
            help="generations until time of joining for pop1 and pop2")
    parser.add_argument('-t23', '--time_join23', type=int, default=240, 
            help="generations until time of joining for pop2 and pop3")
    parser.add_argument('-t34', '--time_join34', type=int, default=240, 
            help="generations until time of joining for pop3 and pop4")
    parser.add_argument('-t2a', '--time2Ancestral', type=int, default=1800, 
            help="generations until ancestral population for PNG is 1800 generations for Africa/Haiti 500 generations")
    parser.add_argument('-at', '--thetaAncestral', type=int, default=344, 
            help="ancestral theta before")
    parser.add_argument('-ar', '--thetaRegional', type=int, default=23, 
            help="regional theta")
    parser.add_argument('-r', '--recombination', type=list, default=0, 
            help="recombination for each locus, if 0 assumes haploid is 1")
    #trans_init
    parser.add_argument('-mt', '--muTrans', type=int, default=100, help="mu for neg bino in transmission, distances between hosts")
    parser.add_argument('-st', '--sizeTrans', type=int, default=1, help="size for neg bino in transmission, distance between hosts")
    parser.add_argument('-dp', '--sigma', type=int, default=2, help="sigma, for dispersal")
    #worm_burden
    parser.add_argument('-prev', '--prev', type=list, required=True, help="prevalance of Wb in host populations; should be a list")
    parser.add_argument('-host', '--hostpopsize', type=list, required=True, help="size of host population")
    parser.add_argument('-mw', '--muWormBurden', type=list, help="mu for negative binomial in worm burden per host, add value or defined by uniform")
    parser.add_argument('-sw', '--sizeWormBurden', type=list, help="size for negative binomial in worm burden per host")
    ## wb_sims
    parser.add_argument('-ng', '--numberGens', type=int, required=True, help="total number of generation to run the simulation")
    parser.add_argument('-bpp', '--bites_person', type=int, default=20, help="number of bites recieved per person per hour")
    parser.add_argument('-h2b', '--hours2bite', type=int, default=6, help="number of hours exposed to mosquitoes")
    #wblifecycles
    parser.add_argument('-ma', '--mortalityAdult', type=float, default=0.8863848717161292, help="adult worms dies at rate 0.01 per month or survive at .99")
    parser.add_argument('-mj', '--mortalityJuv', type=float, default=0.865880173963825, help="prob that juv makes it to adult is 0.2 or 0.8177651 per month")
    parser.add_argument('-mf', '--mortalityMF', type=float, default=0.90, help="MF die at rate 0.10 per month or survive at 0.90")
    parser.add_argument('-mh', '--mortalityHost', type=float, default=0.014286, help="Host death per year")
    parser.add_argument('-f', '--fecundity', type=int, default=20, help="mean number of MF born to each Adult per month")
    parser.add_argument('-Dd', '--density_dependence', action="store_true", help="use density dependence")
    #not used
    parser.add_argument('-gtime', '--generations', type=int, default=0.125, help="generation time in years")
    parser.add_argument('-hostmig', '--host_migration_rates', help="list of host migration rates between villages per month")
    args = parser.parse_args()
    return args

def wb_sims(numberGens, burnin, config_file):
    '''main function for simulations
    Parameters
    ---------
    numberGens : int
         how many months to run the simulation
    burning : int
         burning before recording data
    config_file : file
         file with options to be read by configparser()
    Returns
    -------
    
    '''
    #load needed actuarial table from file
    deathdict = {}
    with open("act.tbl",'r') as tbl:
         for line in tbl:
              line = line.strip()
              deathdict["{}".format(line.split()[0])] = list(map(float,
                  line.split()[1:]))
    #initialize
    dfAdult, dfHost, dfSel = wbsims_initialize()
    #set counters
    month = 0
    sim_time = numberGens + burnin
    while month <= sim_time:
        dfJuv, dfHost = transmission() 
        dfAdult, dfJuv, dfMF, dfHost, dfSel = survival(month) 
        if month > burnin:
             record()
        month += 1

#def main():
#    wb_sims(args.*)
#if __name__ == '__main__':
#    main()

