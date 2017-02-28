#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np


def hostmda_fx(village,
               dfHost,
               mda_coverage):
     ''' Assigns MDA to hosts per village based on mda_coverage
     Parameters
     ----------
     villages : int
          number of villages
     dfHost : df
          dataframe containing host information
     mda_coverage : list, float
          coverage of mda for each village

     Returns
     -------
     dfHost : df
          now will mda column as identity 0 or 1, and updates mda_sum column
     '''
     ##add to MDA_cum
     dfHost.ix[dfHost.MDA==1, "MDA_cum"] += 1
     ##rezero MDA identifier
     dfHost.MDA = 0
     ##assign MDA randomly
     for vill in range(len(village)):
          sub = round(mda_coverage[vill] * len(dfHost[dfHost.village == vill]))
          dfHost.ix[np.random.choice(dfHost.index, sub, replace = False), "MDA"] = 1
     return(dfHost)
