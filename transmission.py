#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import math
from sklearn.metrics import pairwise_distances
import random
from agehost import agehost_fx

def vectorbite_fx(bitespperson, 
                  hours2bite, 
                  hostpopsize,
                  prev_t,
                  densityDep, 
                  avgMF, 
                  bednets, 
                  bnstart, 
                  bnstop, 
                  bncoverage, 
                  month):
    '''Counts the number of successful infectious mosquito bites

    Parameters
    --------
    bitespperson : int
        rate of bites per person per unit time, here hours
    hours2bite : int
        number of hours mosq bite per night/day
    hostpopsize : int
        total population of village 
    prev_t : float
        percent of people infected 
    densityDep : Boolean
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
    print("vectorbite")
    if bednets:
        if month > bnstart and month < bnstop:
             totalbites = ((1 - bncoverage) * bitespperson * hours2bite * 30 
                                   * hostpopsize)
        else:
             totalbites = (bitespperson * hours2bite * 30 
                                   * hostpopsize)
    else:
        totalbites = (bitespperson * hours2bite * 30 
                              * hostpopsize)
    # 0.37 is prob of bite on infected host picking up MF    
    infbites = np.random.binomial(totalbites, (prev_t * 0.37))
    if densityDep: #values for anopheles from CITE
       #number of MF in 20ul of blood
       #235ml is 5% of total host blood 
       #there are 50 units of 20ul in 1ml
       mfBlood = avgMF / 50.0 
       # 0.414 is proportion of  L3 that leave mosquito per bite 
       # 0.32 proportion of L3 that make it into the host
       L3trans = round(infbites * (4.395 * (1 - math.exp( -(0.055 * (mfBlood)) 
           / 4.395)) ** 2) * (0.414 * 0.32)) 
    else:
       # 0.414 is proportion of  L3 that leave mosquito per bite 
       # 0.32 proportion of L3 that make it into the host
       L3trans = np.random.binomial(infbites, (0.414 * 0.32))
    print(int(L3trans))
    return(int(L3trans))


def new_infection_fx(dispersal, transMF, disthost, dfHost):
    '''A transmission event infecting a naive host

    Parameters
    ----------
    mindist : int
         minimum distance
    dispersal: float
         dispersal distance as 2*sigma
    dfHost: df
        host dataframe 
    transMF
        row of transmitted MF
    Returns
    ------
    dfHost: df
         updated with new host information
    newhostidx: str
          new host index     
    '''
    print("newfx")
    #how close
    mindist = 5
    #how far
    maxdist = dispersal
    #assign new position
    maxc = int( math.sqrt(( maxdist ** 2 / 2 )))
    #makes the below statement true
    newpts = dfHost[dfHost.hostidx == transMF.hostidx].coordinates.values[0]

    while next((True for elem in dfHost['coordinates'].values if np.array_equal(elem, newpts)),False) is True:
         #shift of new points
         print("loop")
         x_new = random.randint( -maxc, maxc )
         y_new = random.randint( -maxc, maxc )
         #test that new shifts are min distance
         if math.sqrt(x_new ** 2 + y_new ** 2) > mindist:
              newhostptX = dfHost.coordinates[dfHost.hostidx == transMF.hostidx].values[0][0] + x_new
              newhostptY = dfHost.coordinates[dfHost.hostidx == transMF.hostidx].values[0][1] + y_new
              newpts = np.array([newhostptX, newhostptY])
              print(newpts)
                      #[12, 12], [36, 36], 1, dfMF, dfJuv, dfHost, deathdict) 
    print("outloop")          
    #copy village
    vill = transMF.village
    #new host index
    old_hostidx = dfHost[dfHost.village == transMF.village].hostidx.iloc[-1]
    new_hostidx = old_hostidx[:old_hostidx.rfind('h')] + 'h' + str(int(old_hostidx.split('h')[1]) + 1)
    #radnom sex
    sex = random.choice("01")
    #age and agedeath function from wbsims_initialize
    age, agedeath = agehost_fx(sex)
    #add to dfHost at bottom
    dfHost.loc[len(dfHost) + 1] = [vill, new_hostidx, sex, age, agedeath, newpts, 0]
    
    return dfHost, new_hostidx

def transmission_fx(villages, hostpopsize, sigma, bitesPperson, hours2bite, 
        densitydep_uptake, bednets, bnstart, bnstop, bncoverage, month, 
        dfMF, dfJuv, dfHost):
    '''Transmission events resolved as either reinfection or new infection

    Parameters
    --------
    villages : int
        number of villages
    hostpopsize : int, list
        total population size of village      
    sigma : float
        dispersal mean
    bitespperson : int, list
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
    dfMF : df
         donating MF genotype to transmit
    dfJuv : df
         donated MF move to Juvenille age class
    deathdict : dictionary
         actuarial table, loaded as dictionary
    Returns
    ------
    dfJuv : df
         all new juv are age class 0
    dfHost : df
         in the case of newinfection, new host
    '''
    print("transmission")
    dispersal = 2 * sigma
    from IPython import embed
    embed()
    for vill in range(villages):
        infhost = (dfHost.village == vill).sum()
        prev_t = infhost / float(hostpopsize[vill])

        avgMF = len(dfMF[dfMF.village == vill])/float(infhost)
        avgMF = 200
        L3trans = vectorbite_fx(bitesPperson[vill], hours2bite[vill], hostpopsize[vill], 
                                 prev_t, densitydep_uptake, avgMF, bednets[vill], 
                                 bnstart[vill], bnstop[vill], bncoverage[vill], month)    
        if L3trans > len(dfMF[dfMF.village == vill]):
              transMF = dfMF[dfMF.village == vill]     
        else:
            transMF = dfMF[dfMF.village == vill].sample(L3trans)
        distMat = pairwise_distances(np.vstack(dfHost[dfHost.village == vill].coordinates))
###alternative way only returns points whose distance is <= dispersal
         #meaning these are the only hosts that can transmit to each other and 
         #disMat is a dict with keys as tuples
         #from scipy.spatial import KDTree
         #tree = KDTree(np.vstack(dfHost[dfHost.village == vill].coordinates))
         #distMat = tree.sparse_distance_matrix(np.vstack(dfHost[dfHost.village == vill].coordinates), dispersal)
         #rehost = [i for i, x in enumerate(distMat.keys()) if x[0] == index]                                    
               
        for index, row in transMF.iterrows():
             dfdistHost = dfHost[dfHost.village == vill]
             #new infection
             disthost = np.where((dfdistHost["hostidx"] == row.hostidx))[0]
             if len(dfHost[dfHost.village == vill]) < hostpopsize[vill]:
                  prob_newinfection = 1.0/(len(distMat[disthost] <= dispersal) + 1)
             else: #everyone is already infected
                  prob_newinfection = 0
              
             print(prob_newinfection, row, index) 
             if np.random.random() < prob_newinfection:     
                  print("new loop")
                  #new host
                  dfHost, newidx = new_infection_fx(dispersal, row, disthost, dfHost)
                  row.hostidx = newidx
                  row.age = 0
                  dfJuv = dfJuv.append(row)
                  #need to update distMat to include new host
                  distMat = pairwise_distances(np.vstack(dfHost[dfHost.village == vill].coordinates)) 
             else: #reinfection
                  print("reinfect loop")
                  rehost = dfHost.iloc[random.choice(np.where((distMat[disthost] <= dispersal)[0])[0])]
                  row.hostidx = rehost.hostidx
                  row.age = 0
                  dfJuv = dfJuv.append(row)                          
    return dfHost, dfJuv
