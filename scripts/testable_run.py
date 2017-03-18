# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
import subprocess
import math
import random
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from scipy.stats import weibull_min
import matplotlib.pyplot as plt
from collections import defaultdict
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import os.path
import pickle
from figs.worm import Worms
import figs.wbsims_initialize as wbinit
import figs.transmissionKDtree as trans
from figs.village import Village
from figs.figs2stats import output_tables_fx
#import figs.figs2stats as figstats

def wb_sims(config_file):
    '''main function for simulations
    Parameters
    ---------
    numberGens : int
         how many months to run the simulation
    burn-in : int
         burn-in before recording data
    config_file : file
         file with options to be read by configparser()
    Returns
    -------
    '''

    config = configparser.ConfigParser()
    config.read(config_file)

    # simulation
    sh = "simulation"
    numberGens = config.getint(sh, "numberGens")
    burn_in = config.getint(sh, "burn_in")
    seed = config.getint(sh, "seed")

    # host_demography
    sh = 'host_demography'
    villages = config.getint(sh, 'villages')
    hostpopsize = list(map(int, config.get(sh, 'hostpopsize').split(",")))
    prevalence = list(map(float, config.get(sh, 'prevalence').split(",")))
    muWormBurden = list(map(int, config.get(sh, 'muWormBurden').split(",")))
    sizeWormBurden = list(
        map(int, config.get(sh, 'sizeWormBurden').split(",")))
    assert villages == len(hostpopsize)
    muTrans = config.getint(sh, 'muTrans')
    sizeTrans = config.getint(sh, 'sizeTrans')

    # host_demography - between village parameters
    hostmigrate = config.getfloat(sh, 'hostmigrate')
    initial_distance_m = list(
        map(int, config.get(sh, 'initial_distance_m').split(",")))
    assert len(initial_distance_m) == (villages * (villages - 1)) / 2

    # vector
    sh = 'vector'
    sigma = config.getfloat(sh, 'sigma')
    bitesPperson = list(map(int, config.get(sh, 'bitesPperson').split(",")))
    hours2bite = list(map(int, config.get(sh, 'hours2bite').split(",")))
    densitydep_uptake = config.getboolean(sh, "densitydep_uptake")

    # parasite
    sh = 'parasite'
    fecund = config.getint(sh, 'fecund')
    surv_Juv = config.getfloat(sh, 'surv_Juv')
    shapeMF = config.getfloat(sh, 'shapeMF')
    scaleMF = config.getfloat(sh, 'scaleMF')
    shapeAdult = config.getfloat(sh, 'shapeAdult')
    scaleAdult = config.getfloat(sh, 'scaleAdult')
    densitydep_surv = config.getboolean(sh, 'densitydep_surv')
    densitydep_fec = config.getboolean(sh, 'densitydep_fec')

    # genetic
    sh = 'genetic'
    locus = config.getint(sh, 'locus')
    initial_migration = config.getfloat(sh, 'initial_migration')
    theta = [list(map(int,i.split(","))) for i in config.get(sh, 'theta').split()]
    basepairs = list(map(int, config.get(sh, 'basepairs').split(",")))
    mutation_rate = list(map(float, config.get(sh, 'mutation_rate').split(",")))
    recombination_rate = list(map(float, config.get(sh, 'recombination_rate').split(",")))
    time2Ancestral = config.getint(sh, 'time2Ancestral')
    thetaRegional = config.getfloat(sh, 'thetaRegional')
    time_join = config.getint(sh, 'time_join')

    #selection
    sh = 'selection'
    selection = config.getboolean(sh, 'selection')
    fitness = config.getint(sh, 'fitness')
    surv_length = [list(map(int,i.split(","))) for i in config.get(sh, 'surv_length').split()]
    surv_position = [list(map(float,i.split(","))) for i in config.get(sh, 'surv_position').split()]
    fec_length = [list(map(int,i.split(","))) for i in config.get(sh, 'fec_length').split()]
    fec_position = [list(map(float,i.split(","))) for i in config.get(sh, 'fec_position').split()]
    dominance = config.get(sh, 'dominance')

    # treatment
    sh = 'treatment'
    bednets = config.getboolean(sh, 'bednets')
    bnstart = list(map(int, config.get(sh, 'bnstart').split(",")))
    bnstop = list(map(int, config.get(sh, 'bnstop').split(",")))
    bncoverage = list(map(float, config.get(sh, 'bncoverage').split(",")))
    mda = config.getboolean(sh, 'mda')
    mda_start = config.getint(sh, 'mda_start')
    mda_num = config.getint(sh, 'mda_num')
    mda_freq = config.getint(sh, 'mda_freq')
    mda_coverage = list(map(float, config.get(sh, 'mda_coverage').split(",")))
    mda_macro = config.getfloat(sh, 'mda_macro')
    mda_juvicide= config.getfloat(sh, 'mda_juvicide')
    mda_micro = config.getfloat(sh, 'mda_micro')
    mda_sterile = config.getfloat(sh, 'mda_sterile')
    mda_clear = config.getint(sh, 'mda_clear')

    # outfiles
    sh = 'outfiles'
    logTime=config.getint(sh, 'logdf2file')
    sample_size=config.getfloat(sh, 'sample_size')
    window_length=config.getint(sh, 'window_length')
    num_windows=config.getint(sh, 'num_windows')
    wb2vcf=config.getboolean(sh, 'wb2vcf')
    figs2scikit=config.getboolean(sh, 'figs2scikit')

    outstats = [basepairs, sample_size, window_length, num_windows, wb2vcf, figs2scikit, burn_in]

    #start intialize
    if mda:
        if selection and fitness == 1:
             from figs.survival_mda import survivalmda_sel1_fx as survfx
             print("\nUsing MDA and Selection, fitness is 1")
        elif selection and fitness == 2:
             from figs.survival_mda import survivalmda_sel2_fx as survfx
             print("\nUsing MDA and Selection, fitness is 2")
        else:
             from figs.survival_mda import survivalmda_fx as survfx
             print("\nUsing MDA and NO Selection")
    else:
        from figs.survival import survivalbase_fx as survfx
        print("\nUsing Base Model and Selection is {}".format(selection))

    #print("\n\nSelection is {}\nMDA is {}\nFitness is {}\n\n".format(selection, mda, fitness))
    dist = [0]
    dist.extend(initial_distance_m)
    distvill = [sum(dist[:i+1]) for i in range(len(dist))]
    village=[]

    if os.path.isfile("dfworm_burnin.pkl"):
        with open('dfworm_burnin.pkl', 'rb') as input:
            dfworm = pickle.load(input)
        with open("dfHost_burnin.pkl",'rb') as input:
            dfHost = pickle.load(input)
        sim_time = numberGens
        burn_in = 0
        for i in range(villages):
            village.append(Village(i,hostpopsize[i],prevalence[i],distvill[i], hours2bite[i],
                               bitesPperson[i], bednets, bnstart[i] + burn_in, bnstop[i] + burn_in,
                               bncoverage[i], muTrans, sizeTrans))

    else:
        sim_time = numberGens + burn_in
        cdslist = [surv_length, surv_position, fec_length, fec_position, dominance]
        for i in range(villages):
            village.append(Village(i,hostpopsize[i],prevalence[i],distvill[i], hours2bite[i],
                                   bitesPperson[i], bednets, bnstart[i] + burn_in, bnstop[i] + burn_in,
                                   bncoverage[i], muTrans, sizeTrans))
        dfHost, dfworm =\
                 wbinit.wbsims_init(village,
                                   hostpopsize,
                                   prevalence,
                                   muWormBurden,
                                   sizeWormBurden,
                                   locus,
                                   initial_migration,
                                   initial_distance_m,
                                   theta,
                                   basepairs,
                                   mutation_rate,
                                   recombination_rate,
                                   time2Ancestral,
                                   thetaRegional,
                                   time_join,
                                   selection,
                                   cdslist)

    mdalist = [mda_start + burn_in, mda_num, mda_freq, mda_coverage, mda_macro, mda_juvicide,
            mda_micro, mda_sterile, mda_clear]
    L3transdict = defaultdict(list)
    R0netlist = defaultdict(list)

####start sims
    for month in range(1,sim_time):
        print("\nmonth is {}\n".format(month))
        village, dfHost, dfworm, L3transdict = trans.transmission_fx(month,
                                                            village,
                                                            sigma,
                                                            densitydep_uptake,
                                                            dfHost,
                                                            dfworm,
                                                            L3transdict)
        dfHost, dfworm, R0netlist = survfx(month,
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
                                cdslist)
        #print(dfworm.meta.shape[0])
        if dfworm.meta.shape[0] == 0:
            break
            print("\n\n****POPULATION IS EXTINCT***\n\n")
        #store intialized after burnin
        if month == burn_in:
            with open('dfworm_burnin.pkl', 'wb') as output:
                pickler = pickle.Pickler(output, -1)
                dfworm_burnin = Worms(dfworm.meta, dfworm.h1, dfworm.h2, dfworm.pos,
                                       dfworm.sel, dfworm.coord)
                pickler.dump(dfworm_burnin)
            del dfworm_burnin
            with open('dfHost_burnin.pkl', 'wb') as output:
                pickle.dump(dfHost, output, -1)
        #log data
        elif month > burn_in and month%logTime == 0:
            with open('dfworm_{}.pkl'.format(month - burn_in), 'wb') as output:
                pickler = pickle.Pickler(output, -1)
                dfworm_x = Worms(dfworm.meta, dfworm.h1, dfworm.h2, dfworm.pos,
                                       dfworm.sel, dfworm.coord)
                pickler.dump(dfworm_x)
            del dfworm_x
            with open('dfHost_{}.pkl'.format(month - burn_in), 'wb') as output:
                pickle.dump(dfHost, output, -1)
        else: pass

    #print out last rep
    with open('dfworm_final.pkl', 'wb') as output:
        pickler = pickle.Pickler(output, -1)
        dfworm_x = Worms(dfworm.meta, dfworm.h1, dfworm.h2, dfworm.pos,
                               dfworm.sel, dfworm.coord)
        pickler.dump(dfworm_x)
    del dfworm_x
    with open('dfHost_final.pkl', 'wb') as output:
        pickle.dump(dfHost, output, -1)
    with open('L3transdict.pkl'.format(month), 'wb') as output:
        pickle.dump(L3transdict, output, -1)
    with open('R0netlist.pkl'.format(month), 'wb') as output:
        pickle.dump(R0netlist, output, -1)

    #start stats
    import ipdb; ipbd.set_trace()
    print("\n\ncalculating stats\n\n")
    output_tables_fx(logTime, numberGens, outstats)

    return(seed, "Selection is {}\nMDA is {}\nFitness is {}".format(selection, mda, fitness))

if __name__ == '__main__':
     # this probably needs to be run for at least 240 - 360 months to get away from starting conditions
     seed, model = wb_sims('../tests/wbsims.cfg')
     print("\nSEED\t{}\n".format(seed))
     print(model)