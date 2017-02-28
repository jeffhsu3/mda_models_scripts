#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:11:26 2017

@author: scott
"""

import cPickle as pickle

def pickdeathdict():
#dictionary from actuarial tables
    deathdict = {}
    with open("../act.tbl",'r') as tbl:
         for line in tbl:
              line = line.strip()
              deathdict["{}".format(line.split()[0])] = list(map(float,
                  line.split()[1:]))
    pickle.dump( deathdict, open( "../acttable.p", "wb" ) )
