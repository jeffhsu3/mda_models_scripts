#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
FiGS Copyright (C) 2017 Scott T. Small
This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
This is free software, and you are welcome to redistribute it
under certain conditions; type `show c' for details.
'''
import numpy as np
cimport numpy as np
import random
import pandas as pd
import cython
from cpython cimport array
from libc.stdlib cimport rand, RAND_MAX
#from cython.parallel import parallel, prange
#from libc.stdlib cimport abort, malloc, free

DTYPE = np.uint8
ctypedef np.uint8_t DTYPE_t

cdef long[:] sorted_random_ints(long[:] pos, int size, float[:] weight_array):
    cdef long[:] random_ints = np.random.choice(pos, size=size, p=weight_array)
    return(np.sort(random_ints))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef float[:] weighted_random_index(int basepairs, unsigned long[:] pos):
    cdef np.intp_t i
    cdef float[:] weight_array = np.empty(pos.shape[0] + 1, dtype=np.float32)
    cdef int prev_value
    prev_value = 0
    for i in range(pos.shape[0]):
        print(pos[i] - prev_value)
        weight_array[i] = (pos[i] - prev_value)/float(basepairs)  
        prev_value = pos[i]
    # The last interval
    weight_array[pos.shape[0] + 1] = (basepairs - pos[pos.shape[0]])
    return(np.sort(weight_array))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef np.ndarray[dtype=np.uint8_t, ndim=2] mate_worms(
        long[:] mate_array, 
        long[:] fec, 
        unsigned long[:] pos,
        int basepairs,
        float recomb_rate,
        np.ndarray[DTYPE_t, ndim=2, mode='c'] fem,
        np.ndarray[DTYPE_t, ndim=2, mode='c'] males):
    """ Mates and recombines at a given loci

    Parameters
    ----------
    mate_array : array of longs
        matches females with males
    fec : number of children each female has
    pos : array
    """
    # :TODO need to check max integer
    cdef np.intp_t i, j, l, prev_break, c_break
    cdef np.int64_t outsize
    cdef int k
    outsize = np.sum(fec)
    cdef int mnworms, fnworms
    cdef int hapc, recomb_pos, ohapc 
    cdef long[:] iix_ma = np.repeat(mate_array, fec)
    cdef long[:] femindex = np.arange(fem.shape[0]/2, dtype=np.int64)
    cdef float[:] weight_array 
    cdef long[:] posarray = np.arange(fem.shape[1] + 1, dtype=np.int64)
    cdef np.ndarray iix_fem = np.repeat(femindex, fec)
    cdef np.ndarray mnum_recomb = np.random.poisson(
            recomb_rate * basepairs, outsize)
    cdef np.ndarray fnum_recomb = np.random.poisson(
            recomb_rate * basepairs, outsize)
    cdef DTYPE_t[:, ::1] hout = np.empty((2*outsize, fem.shape[1]), dtype=np.uint8)
    mnworms = males.shape[0]/2
    fnworms = fem.shape[0]/2
    # Pos must be sorted
    weight_array = weighted_random_index(basepairs, pos)
    print(np.sum(weight_array))
    for i in range(outsize):
        print('Number of recombinations')
        print(mnum_recomb[i])
        print(fnum_recomb[i])
        hapc = np.int(rand()/RAND_MAX)
        if hapc == 0: 
            ohapc = 1
        if mnum_recomb[i] == 0:
            hout[i, :] = males[iix_ma[i] + mnworms * hapc, :]
        else:
            k = 0
            cpos = sorted_random_ints(posarray, mnum_recomb[i], weight_array)
            mhapc = np.int(rand()/RAND_MAX)
            prev_break = 0
            c_break = 0
            male_index = i
            while k < mnum_recomb[i]:
                c_break = cpos[k]
                if c_break == posarray[pos.shape[0] + 1]:
                    continue
                else:
                    hout[i, prev_break:c_break] = males[iix_ma[i] + mnworms *
                            hapc, prev_break:c_break]
                    hout[i, c_break: ] = males[iix_ma[i] + mnworms * ohapc, :]
                    prev_break = c_break
                    hapc = ohapc
                    if hapc == 1:
                        ohapc = 0
                    else:
                        ohapc = 1
                    k += 1
        '''        
        hapc = np.int(rand()/RAND_MAX)
        if hapc == 0: 
            ohapc = 1
        if fnum_recomb[i] == 0:
            hout[i + outsize, :] = fem[iix_fem[i] + fnworms * hapc, :]
        else:
            k = 0
            while k <= mnum_recomb[i]:
                k += 1
                recomb_pos = np.int(rand()/RAND_MAX*basepairs)
                for l in range(len(pos)):
                    if recomb_pos > pos[l]:
                        break
        '''
    return(hout)


def recombination_fx(locus,
                     dfAdult,
                     list recombination_rate,
                     list basepairs):
    """Calculate number of recombination events and rearranges haplotypes
    :TODO add for recombination map

    Parameters
    ---------
    locus: int
        number of loci
    dfAdult_mf : figs.worm.Worms object
        Worms containing new larval parasites
    dfAdult : figs.worm.Worms object
        dataframe containing reproducing adults
    recombination_rate : float, list
        recombination rate for each locus
    basepairs : int, list
        length of each locus in basepairs

    Returns
    -------
    dfAdult_mf : pd df

    """
    hosts = dfAdult.meta.hostidx.unique()
    cdef str host
    # How to type this?
    #cdef bool[:] ahost, females, males
    cdef np.ndarray out_array
    cdef Py_ssize_t loc
    #cdef np.ndarray fec
    cdef float rr
    #cdef int[:] mate_array = np.empty(np.sum(females))
    for host in hosts:
        #chost = dfAdult.meta.index[dfAdult.meta.hostidx == host].values
        ahost = dfAdult.meta.hostidx == host
        females = np.logical_and(ahost, dfAdult.meta.sex == 'F').values
        males = np.logical_and(ahost, dfAdult.meta.sex == 'M').values
        if np.sum(males) == 0 or np.sum(females) == 0:
            print('Either there are 0 males in host or zero females in host')
            continue
        else:
            fec = dfAdult.meta.fec[females].values
            mate_array = np.random.randint(0, np.sum(males), np.sum(females),
                    dtype=np.int64)
            # Parallelize this
            for loc in range(locus):
                rr = recombination_rate[loc]
                if rr == 0:
                    pass
                else:
                    cfemales = np.vstack((dfAdult.h1[str(loc)][females, :],
                        dfAdult.h2[str(loc)][females, :]))
                    cmales = np.vstack((dfAdult.h1[str(loc)][males, :],
                        dfAdult.h2[str(loc)][males, :]))
                    out_array = mate_worms(mate_array,
                            fec,
                            dfAdult.pos[str(loc)],
                            basepairs[loc],
                            rr,
                            cfemales,
                            cmales)
    #new_meta = pd.DataFrame({})
    #dfAdult_mf = Worms(meta = pd
    return(out_array)
