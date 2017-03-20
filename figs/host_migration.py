#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np


def hostmigration_fx(village,
                     dfHost,
                     hostmignumb):
     ''' allows host to move between villages
     Parameters
     ---------
     dfHost : df
          df of hosts
     hostmigrate : float
          host migration rate between villages per year
     Returns
     --------
     dfHost : df
     '''
     move_host = hostmignumb
     i = 0
     while i < move_host:
         print("hostmigration")
         migrant = np.random.choice(dfHost.index, move_host)
         #stepping stone
         for mv in migrant:
             if dfHost.ix[mv, "village"] == min(dfHost.village):
                  #less than max can go up or down
                  dfHost.ix[mv, "village"] += 1
             elif dfHost.ix[mv, "village"] == max(dfHost.village):
                  dfHost.ix[mv, "village"] -= 1
             else: #at max can only go down
                 dfHost.ix[mv, "village"] += np.random.choice([1,-1])
             vill = dfHost.ix[mv].village
             #new coordinates
             newcoords = (np.random.negative_binomial(village[vill].size,
                      village[vill].size / float((village[vill].size+village[vill].mu)), 2)
                      + village[vill].dist)
             dfHost.loc[mv, "coordinates"] = np.vstack(newcoords)
         i += 1
     return(dfHost)
