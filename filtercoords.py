#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np

def filtercoords_fx(positions,
                    basepairs,
                    perc_locus,
                    cds_length,
                    intgen_length):
    ''' Filters coordinates to keep only those in defined cds

    Parameters
    ---------
    locus : int
    positions : list
    basepairs : list, int
    perc_locus : list, float
    cds_length : int
    intgen_length : int

    Returns
    -------
    cds_positions : list
        filterd positions to only include cds
    cds_coordinates : lisit
        coordinates (start,end) of cds

    '''
    cds_positions = []
    cds_coordinates = []
    #size parameter
    size = 3
    #average distance between in bp
    mu = intgen_length
    #last coordinate is approx num_cds * mu; so if num_cds is too large or mu is too long
    #genes will run over the locus length
    for loc in range(len(positions)):
        num_cds = int(round((perc_locus[loc + 1]*basepairs[loc + 1]) /
            cds_length[loc +1]))
        size_cds = np.round(np.random.gamma(4, 0.25, num_cds) *
                cds_length[loc+1])
        #r = size
        #m = mean
        #p = r / (  r + m )
        cds_between = np.random.negative_binomial(size, size/float(mu+size), num_cds)
        cds_stop = 0
        cds_coords = []
        for i, j in zip(cds_between, size_cds):
            #[i + cds_stop, i + cds_stop + j]
            if (i + cds_stop > basepairs[loc + 1]) or (i + j + cds_stop > basepairs[loc + 1]):
                break
            else:
                cds_coords.append([i + cds_stop, i + j + cds_stop])
                cds_stop += (i + j)

        keep_pos = []
        for start, end in cds_coords:
            keep_pos.extend(positions[loc][np.where(np.logical_and(positions[loc] >= start, positions[loc] <= end))])
        cds_positions.append(keep_pos) #this is a nested list for each locus
        cds_coordinates.append(cds_coords) #this is a nested list for each locus
    return(cds_positions, cds_coordinates)   
