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
from collections import defaultdict
from figs.worm import Worms
#from cython.parallel import parallel, prange
#from libc.stdlib cimport abort, malloc, free

DTYPE = np.uint8
ctypedef np.uint8_t DTYPE_t

cdef long[:] sorted_random_ints(long[:] pos, int size, double[:] weight_array):
    cdef long[:] random_ints = np.random.choice(pos, size=size, p=weight_array)
    return(np.sort(random_ints))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double[:] weighted_random_index(int basepairs, long[:] pos):
    cdef np.intp_t i
    cdef size_t weight_shape = pos.shape[0] + 1
    cdef size_t pos_len = pos.shape[0]
    cdef double[:] weight_array = np.empty(weight_shape, dtype=np.float64)
    cdef int prev_value
    prev_value = 0
    for i in range(pos.shape[0]):
        weight_array[i] = (pos[i] - prev_value)/float(basepairs)  
        prev_value = pos[i]
    # The last interval
    weight_array[i + 1] = (basepairs - pos[i])/float(basepairs)
    return(np.sort(weight_array))


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef np.ndarray[dtype=np.uint8_t, ndim=2] mate_worms(
        long[:] mate_array, 
        long[:] fec, 
        long[:] pos,
        int basepairs,
        float recomb_rate,
        np.ndarray[DTYPE_t, ndim=2, mode='c'] fem,
        np.ndarray[DTYPE_t, ndim=2, mode='c'] males):
    """ Mates and recombines at a given loci

    Parameters
    ----------
    mate_array : array long
        matches females with males
    fec : array long 
        number of children each female has pos : array long
        genomic positions array same size as axis-1 of fem/males
    recomb_rate: float
    fem : np.uint8 ndim 2
    males : np.uint8 ndim 2
        this array is actually all genotypes male and females
     """
    # :TODO need to check max integer
    cdef np.intp_t i, prev_break, c_break
    cdef np.intp_t hout_index
    cdef np.int64_t outsize
    cdef int k
    outsize = np.sum(fec)
    cdef int mnworms, fnworms
    cdef size_t pos_len = pos.shape[0]
    # Haplotype chooser 
    cdef int hapc, ohapc 
    # Father mate array
    cdef long[:] iix_ma = np.repeat(mate_array, fec)
    cdef long[:] femindex = np.arange(fem.shape[0]/2, dtype=np.int64)
    cdef double[:] weight_array 
    cdef long[:] posarray = np.arange(fem.shape[1] + 1, dtype=np.int64)
    cdef np.ndarray iix_fem = np.repeat(femindex, fec)
    cdef np.ndarray mnum_recomb = np.random.poisson(
            recomb_rate * basepairs, outsize)
    cdef np.ndarray fnum_recomb = np.random.poisson(
            recomb_rate * basepairs, outsize)
    cdef np.ndarray[DTYPE_t, ndim=2] hout = np.empty((2*outsize, fem.shape[1]), 
            dtype=np.uint8)
    mnworms = males.shape[0]/2
    fnworms = fem.shape[0]/2
    # Pos must be sorted
    weight_array = weighted_random_index(basepairs, pos)
    for i in range(outsize):
        hapc = np.int(rand()/RAND_MAX)
        if hapc == 0: 
            ohapc = 1
        if mnum_recomb[i] == 0:
            hout[i, :] = males[iix_ma[i] + mnworms * hapc, :]
        else:
            cpos = sorted_random_ints(posarray, mnum_recomb[i], weight_array)
            prev_break = 0
            c_break = 0
            k = 0
            while k < mnum_recomb[i]:
                c_break = cpos[k]
                if c_break == pos_len + 1:
                    break
                else:
                    hout[i, prev_break:c_break] = males[iix_ma[i] + mnworms *
                            hapc, prev_break:c_break]
                    hout[i, c_break: ] = males[iix_ma[i] + mnworms * ohapc,
                            c_break:]
                    prev_break = c_break
                    hapc = ohapc
                    if hapc == 1: ohapc = 0
                    else: ohapc = 1
                    k += 1
        hapc = np.int(rand()/RAND_MAX)
        
        hout_index = i + outsize
        if hapc == 0: 
            ohapc = 1
        if fnum_recomb[i] == 0:
            hout[hout_index, :] = fem[iix_fem[i] + fnworms * hapc, :]
        else:
            cpos = sorted_random_ints(posarray, fnum_recomb[i], weight_array)
            prev_break = 0
            c_break = 0
            k = 0
            while k < fnum_recomb[i]:
                c_break = cpos[k]
                if c_break == pos_len + 1:
                    break
                else:
                    hout[hout_index, prev_break:c_break] = fem[iix_fem[i] + fnworms *
                            hapc, prev_break:c_break]
                    hout[hout_index, c_break: ] = fem[iix_fem[i] + fnworms * ohapc,
                            c_break:]
                    prev_break = c_break
                    hapc = ohapc
                    if hapc == 1: ohapc = 0
                    else: ohapc = 1
                    k += 1
    return(hout)


@cython.boundscheck(False)
@cython.wraparound(False)
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
    dfAdult_mf : figs.worm.Worms object
        new worms 
    """
    h1t = {}
    h2t = {}
    cdef str host
    cdef np.ndarray out_array
    cdef Py_ssize_t loc
    cdef np.ndarray[long, ndim=1] fec
    cdef double rr
    # Need to make copy or otherwise awkward downstream
    # behavior
    cdef np.ndarray[DTYPE_t, ndim=2] h1_orig
    cdef np.ndarray[DTYPE_t, ndim=2] h2_orig  

    ###############################################
    # Size of full meta along axis-0
    ###############################################
    cdef np.ndarray[np.uint8_t, cast=True] fem_bool =\
            (dfAdult.meta.sex == 'F').values
    cdef np.ndarray[np.uint8_t, cast=True] mal_bool =\
            (dfAdult.meta.sex == 'M').values

    # Sex boolean arrays
    cdef np.ndarray[np.uint8_t, cast=True] m_host_bool 
    cdef np.ndarray[np.uint8_t, cast=True] f_host_bool 
    cdef np.ndarray[np.uint8_t, cast=True] bad_host_bool =\
            np.ones(dfAdult.meta.shape[0], dtype=np.uint8)  
    

    ###############################################
    # Size of Females only
    ###############################################
    cdef int mother_n = np.sum(fem_bool)
    cdef np.ndarray[np.int64_t] mate_array = np.empty(np.sum(fem_bool), 
            dtype=np.int64)
    cdef np.ndarray[np.uint8_t, cast=True] g_fem_ix
    cdef long total_offspring = dfAdult.meta['fec'][fem_bool].sum()
    
    ###############################################
    # Input arrays axis-0 and number
    # of segregating sties
    ###############################################
    cdef np.ndarray[DTYPE_t, ndim=2] cmales 
    cdef np.ndarray[DTYPE_t, ndim=2] cfemales 
    # Generate mate array
    hosts = dfAdult.meta.hostidx.unique()
    for host in hosts:
        m_host_bool = np.logical_and(mal_bool,
                (dfAdult.meta.hostidx == host).values)
        f_host_bool = np.logical_and(fem_bool,
                (dfAdult.meta.hostidx == host).values)
        if (np.sum(m_host_bool) != 0) and (np.sum(m_host_bool) != 0):
            mate_array[f_host_bool[fem_bool]] =\
                    np.random.choice(dfAdult.meta.index[m_host_bool].values, 
                    size = np.sum(f_host_bool))
        else: 
            bad_host_bool[f_host_bool] = False 
    # :TODO This should be made more simple
    g_fem_ix = np.logical_and((dfAdult.meta.fec > 0).values, 
            bad_host_bool)
    g_fem_ix = np.logical_and(g_fem_ix, fem_bool)
    fec = dfAdult.meta.fec[g_fem_ix].values
    total_offspring= np.sum(fec)
    # Removing bad hosts and zero_fecundity females
    mate_array = mate_array[g_fem_ix[fem_bool]]
    assert np.sum(g_fem_ix) == mate_array.shape[0]

    # :TODO parallelize this
    for loc in range(locus):
        rr = recombination_rate[loc]
        if rr == 0:
            h1t[str(loc)] = np.repeat(dfAdult.h1[str(loc)][g_fem_ix,:],
                        fec, axis=0)
        else:
            cfemales = np.vstack((dfAdult.h1[str(loc)][g_fem_ix, :],
                dfAdult.h2[str(loc)][g_fem_ix, :]))
            cmales = np.vstack((dfAdult.h1[str(loc)][:, :],
                dfAdult.h2[str(loc)][:, :]))
            # :TODO passing in male genotypes instead of all
            out_array = mate_worms(
                    mate_array,
                    fec,
                    dfAdult.pos[str(loc)],
                    basepairs[loc],
                    rr,
                    cfemales,
                    cmales)
            h1t[str(loc)] = out_array[0:total_offspring, :]
            h2t[str(loc)] = out_array[total_offspring:, :]
    new_meta = pd.DataFrame({
        'village' : np.repeat(dfAdult.meta.village[g_fem_ix], fec),
        'sex' : np.random.choice(['M', 'F'], size = total_offspring),
        'hostidx' : np.repeat(dfAdult.meta.hostidx[g_fem_ix], fec),
        'fec' : np.repeat(0, total_offspring),
        'R0net' : np.repeat(dfAdult.meta['R0net'][g_fem_ix], fec),
        'age' : np.repeat(0, total_offspring),
        })
    new_meta.reset_index(drop=True, inplace=True)
    dfAdult_mf = Worms(meta = new_meta.ix[:, dfAdult.meta.columns], 
            haplotype1 = h1t, haplotype2=h2t,
            positions=dfAdult.pos.copy())
    return(dfAdult_mf)
