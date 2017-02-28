#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 16:33:14 2017

@author: scott
"""

def hapfreq(hap):
    ''' a very sloppy way to get haplotype freq from column in pandas df
    Parameters
    ---------
    hap : list, array
        dfAdult.locus_0
    Returns
    -------
    counts : dict values
        counts of each haplotype
    '''
    #TO DO: when/if haps are full gt array, many of these step will be uneeded
    a=[]
    [a.append(i) for i in hap]
    b=[]
    [b.append(str(list(i))) for i in a]
    counts = dict()
    for i in b:
        counts[i] = counts.get(i,0) + 1
    return(counts.values())

def hapcount(hap):
    ''' counts the number of unique haplotypes in pandas column
    Parameters
    ----------
    hap : list, array
        dfAdult.locus_0
    Returns
    ------
    uniqhap : int
        int of unique haps
    '''
    uniqhap = len(set([tuple(v) for v in hap]))
    return(uniqhap)
