#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import bisect

def agehost_fx(sex, deathdict):
    '''Calculate age of host and age of death from acturial table
    Parameters
    ---------
    sex : str
         male or female

    Returns
    ------
    age : int
         age of host
    death : int
         age of death
    '''
    ageclass = [0.24, 0.46, 0.64, 0.78, 0.88, 0.94, 0.98]
    agerange = [[0,10],[11,20],[21,30],[31,40],[41,50],[51,60],[61,70],[71,80]]
    ageint = bisect.bisect_left(ageclass, np.random.random())
    age = np.random.randint(agerange[ageint][0], agerange[ageint][1])
    
    #when do they die
    death = deathdict[str(age)][int(sex)] + np.random.normal(0,6) + age
         
    return(age, round(death))
