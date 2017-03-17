#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 11:33:55 2017

@author: scott
"""
import numpy as np
import pandas as pd
#from .plotting import plot_coordinates_host

def prevTrans_fx(L3transdict, logTime):
    '''calculates the mean prev and transmission for summary table

    Parameters
    ---------

    Returns
    -------
    '''
    print("prevTrans")
    #trans_prev
    vill = len(L3transdict['prev'][0])
    time = len(L3transdict['prev'])
    trans = [e for l in L3transdict['trans'] for e in l]
    prev = np.round([e for l in L3transdict['prev'] for e in l],2)
    month = np.repeat(range(1, time + 1), vill)
    village = range(vill) * time

    transTable = pd.DataFrame({"month" : month,
                                "village" : village,
                                "trans" : trans,
                                "prev" : prev})
    transTable = transTable.loc[:, ['month','village', 'trans', 'prev']]
    meantran = [transTable.groupby("village").rolling(logTime).mean().dropna()['trans'][v::logTime][v] for v in range(vill)]
    vartran = [transTable.groupby("village").rolling(logTime).var().dropna()['trans'][v::logTime][v] for v in range(vill)]
    meanprev = [transTable.groupby("village").rolling(logTime).mean().dropna()['prev'][v::logTime][v] for v in range(vill)]
    varprev = [transTable.groupby("village").rolling(logTime).var().dropna()['prev'][v::logTime][v] for v in range(vill)]

#    avgTrans = [val for pair in zip(*meantran)[] for val in pair]
#    varTrans = [val for pair in zip(*vartran)[] for val in pair]
#    avgPrev = [val for pair in zip(*meanprev)[] for val in pair]
#    varPrev = [val for pair in zip(*varprev)[] for val in pair]
    transTable = transTable.round(4)
    #write to out
    transTable.to_csv("transTable.csv")

    return(meantran, meanprev, vartran, varprev, vill, month)

def R0net_fx(R0netlist, month):
    '''calculate the reproductive rate, R0 for summary table

    Parameters
    ---------

    Returns
    -------
    '''
    print("R0")
    #R0 by village
    R0 = []
    R0_t = zip(*R0netlist['R0'])
    for v in range(len(R0_t)):
        R0.append([float(j) / R0_t[v][i-1] for i,j in enumerate(R0_t[v])])

    #reproductive mean and var
    ravg = R0netlist['repoavg']
    rvar = R0netlist['repovar']

    R0Table = pd.DataFrame({"R0" : [j for i in zip(*R0) for j in i],
                            "ravg" : [j for i in ravg for j in i],
                            "rvar" : [j for i in rvar for j in i]
                            })
    R0Table = R0Table.loc[:, ['R0', 'ravg', 'rvar']]
    R0Table = R0Table.round(4)
    R0Table.to_csv("R0Table.csv")
    return(R0, ravg, rvar)

def demo_stats_fx(dfworm, vill):
    '''calculates various demographic statistics for summary table

    Parameters
    ---------

    Returns
    -------
    '''
    print("demostats")
    adult = dfworm.meta.groupby(["stage","village","hostidx"]).size()['A'][vill].mean()
    vadult = dfworm.meta.groupby(["stage","village","hostidx"]).size()['A'][vill].var()
    juv = dfworm.meta.groupby(["stage","village","hostidx"]).size()['J'][vill].mean()
    vjuv = dfworm.meta.groupby(["stage","village","hostidx"]).size()['J'][vill].var()
    mf = dfworm.meta.groupby(["stage","village","hostidx"]).size()['M'][vill].mean()
    vmf = dfworm.meta.groupby(["stage","village","hostidx"]).size()['M'][vill].var()

    return(adult, vadult, juv, vjuv, mf, vmf)

#################################
def demo_hoststats_fx(dfworm, dfHost, mon):
    '''calculates various demographic statistics for each hostidx

    Parameters
    ---------

    Returns
    -------
    '''
    print("demohoststats")
    adult = dfworm.meta.groupby(["stage","hostidx"])
    juv = dfworm.meta.groupby(["stage","hostidx"])
    mf = dfworm.meta.groupby(["stage","hostidx"])
    hostidx = dfworm.meta.hostidx.unique()
    nadult = []
    njuv = []
    nmf = []
    for host in hostidx:
        try:
            nadult.append(adult.size()['A'][host])
        except KeyError:
            nadult.append(0)
        try:
            njuv.append(juv.size()['J'][host])
        except KeyError:
            njuv.append(0)
        try:
            nmf.append(mf.size()['M'][host])
        except KeyError:
            nmf.append(0)

    demohostTable = pd.DataFrame({"month" : mon,
                                  "hostidx" : hostidx,
                                  "adult" : nadult,
                                  "juv" : njuv,
                                  "mf" : nmf
                                  })
    demohostTable = demohostTable.loc[:, ['month', 'hostidx', 'adult', 'juv', 'mf']]

    demohostTable.to_csv('demohostTable.csv')
    return(None)
###################################
def host_stats_fx(dfHost,villages):
    '''calculates stats for the specific host

    Parameters
    ---------

    Returns
    -------
    '''
    print("infhost")
    infhost = [dfHost.groupby("village").size()[i] for i in range(villages)]
#    plot_coordinates_host(dfHost, thetaHost)
    return(infhost)