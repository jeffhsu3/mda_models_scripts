#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import pandas as pd
import math
import random
import pickle
import numpy as np
from scipy.spatial import cKDTree
import math
#import warnings
#warnings.simplefilter("error", pd.core.common.SettingWithCopyWarning)
from .agehost import agehost_fx

deathdict = pickle.load(open('../figs/data/acttable.p', "rb"))

def vectorbite_fx(vill,
                  prev_t,
                  month,
                  village,
                  densitydep_uptake,
                  avgMF):

    '''Counts the number of successful infectious mosquito bites

    Parameters
    --------
    bitesPperson : int
        rate of bites per person per unit time, here hours
    hours2bite : int
        number of hours mosq bite per night/day
    hostpopsize : int
        total population of village
    prev_t : float
        percent of people infected
    densitydep_uptake : Boolean
        use the density dependece function for developing L3
    avgMF : float
         average number of MF per host per village
    bednets : Boolean, list
        use values from bednets
    bnstart : int, list
        month to start bednets
    bnstop : int, list
        mont to stop bednets
    bncoverage : int, list
        percent of population using bednets
    month : int
        current time, month
    Returns
    ------
    L3trans : int
         all new juv are age class 0
    '''

    bitesPperson = village[vill].bpp
    hours2bite = village[vill].h2b
    hostpopsize = village[vill].hostpopsize
    bednets = village[vill].bn
    bnstart = village[vill].bnstr
    bnstop = village[vill].bnstp
    bncoverage =  village[vill].bncov
    prev_t = prev_t
    #print("vectorbite")
    if bednets and (month > bnstart and month < bnstop):
        print("Start BedNets, vilage {} has coverage {}".format(vill,bncoverage))
        totalbites = ((1 - bncoverage) * bitesPperson * hours2bite * 30
                               * hostpopsize)
    elif bednets and (month > bnstop):
        print("BedNets stop")
        totalbites = (bitesPperson * hours2bite * 30
                              * hostpopsize)
    else:
        totalbites = (bitesPperson * hours2bite * 30
                      * hostpopsize)
    # 0.37 is prob of bite on infected host picking up MF
    infbites = np.random.binomial(totalbites, (prev_t * 0.37))
    if densitydep_uptake: #values for anopheles from doi:10.1371/journal.pone.0002874.s001
       #number of MF in 20ul of blood
       #235ml is 5% of total host blood
       #there are 50 units of 20ul in 1ml
       mfBlood = avgMF / 50.0
       # 0.414 is proportion of  L3 that leave mosquito per bite
       # 0.32 proportion of L3 that make it into the host
       L3trans = round(infbites * (4.395 * ((1 -\
               math.exp( -1* (0.055 * mfBlood / 4.395))) ** 2)) *\
               (0.414 * 0.32))
    else:
       # 0.414 is proportion of  L3 that leave mosquito per bite
       # 0.32 proportion of L3 that make it into the host
       L3trans = np.random.binomial(infbites, (0.414 * 0.32))
    return(int(L3trans))


def new_infection_fx(dispersal,
                     mfhostidx,
                     dfHost):
    '''A transmission event infecting a naive host

    Parameters
    ----------
    dispersal: float
         dispersal distance as 2*sigma
    transMF
        row of transmitted MF
    dfHost: df
        host dataframe

    Returns
    ------
    dfHost: df
         updated with new host information
    newhostidx: str
          new host index
    '''
    #how close
    mindist = 5
    #how far
    maxdist = dispersal
    #assign new position
    maxc = int( math.sqrt(( maxdist ** 2 / 2 )))
    #makes the below statement true
    newpts = dfHost[dfHost.hostidx == mfhostidx].coordinates.values[0]

    while next((True for elem in dfHost['coordinates'].values if np.array_equal(elem, newpts)),False) is True:
         #shift of new points
         x_new = random.randint( -maxc, maxc )
         y_new = random.randint( -maxc, maxc )
         #test that new shifts are min distance
         if math.sqrt(x_new ** 2 + y_new ** 2) > mindist:
              newhostptX = dfHost.coordinates[dfHost.hostidx == mfhostidx].values[0][0] + x_new
              newhostptY = dfHost.coordinates[dfHost.hostidx == mfhostidx].values[0][1] + y_new
              newpts = np.array([newhostptX, newhostptY])
    #copy village
    vill = int(mfhostidx[:mfhostidx.rfind('h')][-1])
    #new host index
    old_hostidx = dfHost[dfHost.village == vill].hostidx.values
    max_hostidx = max(map(int,[x.split("h")[1] for x in old_hostidx]))
    new_hostidx = 'v' + str(vill) + 'h' + str(max_hostidx + 1)
    #radnom sex
    sex = random.choice("01")
    #age and agedeath function from wbsims_initialize
    age, agedeath = agehost_fx(sex, deathdict)
    #mda list
    mdalist = len(dfHost[dfHost.hostidx == mfhostidx].MDAcum.values[0])
    tmdalist = [0] * mdalist
    #add to dfHost at bottom
    newhostlist = [[vill, new_hostidx, sex, age, agedeath, newpts, 0, tmdalist ]]
    #setting with enlargement
    dfHost.ix[dfHost.index[-1] + 1] = newhostlist[0]
    return(dfHost, new_hostidx)


def temp_function(dfworm, transMF, new_hostidx, hi,  dispersal=0):
    """
    hi : dictionary
       host info
    """
    tcount = ''
    transMFidx = dfworm.meta.ix[transMF].hostidx.values
    transMFhostidx = [hi['dfHost'].query('hostidx in @x').index.values \
            for x in transMFidx]

    for mfhostidx, transhostidx in zip(transMFidx, transMFhostidx):
        if mfhostidx != tcount:
            transhost = hi['tree'].query_ball_point(
                    hi['dfHost'].ix[transhostidx[0]].coordinates, 
                    dispersal, n_jobs=-1)
            tcount = mfhostidx
        if hi['infhost'] < hi['hostpopsize']:
             prob_newinfection = 1.0 / (len(transhost) + 1)
        else: #everyone is already infected
             prob_newinfection = 0
        trand = np.random.random()
        if trand < prob_newinfection:
            hi['dfHost'], rehostidx = new_infection_fx(dispersal, 
                    mfhostidx, hi['dfHost'])
            new_hostidx.append(rehostidx)
            #new host so have to resort and rebuild KDTree
            hi['hostcoords'] = np.concatenate([hi['hostcoords'], 
                [hi['dfHost'].ix[hi['dfHost'].index[-1]].coordinates]])
            hi['tree'] = cKDTree(hi['hostcoords'])
            tcount = ''
            hi['infhost'] += 1
        else:
            try:
                rehostidx = np.random.choice(transhost)
                new_hostidx.append(hi['dfHost'].ix[rehostidx,'hostidx'])
            except:
                import ipdb; ipdb.set_trace()
    return(tcount, new_hostidx)

#@profile
def transmission_fx(month,
                    village,
                    sigma,
                    densitydep_uptake,
                    dfHost,
                    dfworm,
                    L3transdict, 
                    juvs):
    '''Transmission events resolved as either reinfection or new infection

    Parameters
    --------
    villages : int
        number of villages
    hostpopsize : int, list
        total population size of village
    sigma : float
        dispersal mean
    bitesPperson : int, list
        rate of bites per person per unit time, here hours
    hours2bite : int, list
        number of hours mosq bite per night/day
    densitydep_uptake : Boolean
        use the density dependece function for developing L3
    bednets : Boolean, list
        use values from bednets
    bnstart : int, list
        month to start bednets
    bnstop : int, list
        mont to stop bednets
    bncoverage : int, list
        percent of population using bednets
    month : int
        current time, month
    dfHost: df
         chooses the donating and accepting hosts for transmission
    juvs: dict

    Returns
    ------
    dfworm : figs.worm.Worms object
        updated object 
    L3trans : int
        number of transmitted MF
    '''
    #print("transmission_fx")
    hostlist =  dfworm.meta["hostidx"].unique()
    if dfHost[dfHost.hostidx.isin(hostlist)].shape[0] != dfHost.shape[0]:
        dfHost = dfHost.loc[dfHost.hostidx.isin(hostlist)]
    dfHost.reset_index(inplace=True, drop=True)
    dispersal = 2 * sigma
    new_hostidx = []
    trans_t = []
    prev_t = []
    hostcoords = np.vstack(dfHost.coordinates)
    tree = cKDTree(hostcoords)
    new_juv = []

    for vill in range(len(village)):
        infhost = (dfHost.village == vill).sum()
        prev = infhost / float(village[vill].hostpopsize)
        prev_t.append(prev)
        
        mfmany = []
        L3_weight_factors = []
        chosen_mfs = []

        ####### Original Worm #######################################
        mfiix_vill = dfworm.meta[(dfworm.meta["stage"] == "M")\
                & (dfworm.meta["village"] == vill)].index.values
        L3_weight_factors.append(mfiix_vill.shape[0])

        for mfs in juvs['worms']:
            mfmany.append(mfs.meta.query(
            "stage=='M' & village==@vill").index.values)
            L3_weight_factors.append(mfiix_vill.shape[0])
        
        avgMF = sum(L3_weight_factors)/float(infhost)
        L3trans = vectorbite_fx(vill, prev, month, village, densitydep_uptake, avgMF)
        trans_t.append(L3trans)
        print("village is %i transmitted is %i" %(vill, L3trans))
        if L3trans != 0:
            if L3trans > sum(L3_weight_factors):  #more transmision events than possible MF
                transMF = mfiix_vill
                chosen_mfs = mfmany
            else:
                transMF = np.random.choice(mfiix_vill,
                        math.ceil(L3_weight_factors[0]/sum(L3_weight_factors)*L3trans), 
                        replace=False)
                for iix, mfs in enumerate(mfmany):
                    c_mf = np.random.choice(mfs, 
                            math.ceil(L3_weight_factors[iix
                                +1]/sum(L3_weight_factors)*L3trans),
                            replace=False)
                    c_mf.sort()
                    chosen_mfs.append(c_mf)
            transMF.sort()
            new_juv.extend(transMF)
            host_info = {'dfHost': dfHost, 
                         'infhost': infhost,
                         'hostcoords': hostcoords,
                         'tree': tree,
                         'hostpopsize': village[vill].hostpopsize 
                         }

            tcount, new_hostidx = temp_function(dfworm, transMF,
                    new_hostidx, host_info, dispersal=dispersal)

            '''
            if len(chosen_mfs) >= 1:
                import ipdb
                ipdb.set_trace()
            '''

            for mfs, tmf in zip(juvs['worms'], chosen_mfs):
                new_hostidx_t = []
                tcount, new_hostidx_t =temp_function(mfs, tmf, 
                        new_hostidx_t,  
                        host_info,
                        dispersal=dispersal)
                mfs.meta.ix[tmf, 'stage'] = "J"
                mfs.meta.ix[tmf, 'hostidx'] = new_hostidx_t

        else:
            print("dfMF is empty")
    L3transdict['prev'].append(prev_t)
    L3transdict['trans'].append(trans_t)
    dfworm.meta.ix[new_juv, 'stage'] = "J"
    dfworm.meta.ix[new_juv, 'hostidx'] = new_hostidx
    return(village, dfHost, dfworm, L3transdict)
