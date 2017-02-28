#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
FiGS Copyright (C) 2017 Scott T. Small
This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
This is free software, and you are welcome to redistribute it
under certain conditions; type `show c' for details.
'''
import numpy as np
import random
import pandas as pd
from functools import partial
from IPython import embed

from locus import Locus


def recombination_locus(h1, h2, crossover_pos):
    """Calculates the recombination at a given locus
    """
    h1_ix = len(h1)
    h2_ix = len(h2)
    for i, j in enumerate(h1):
        if j > crossover_pos:
            h1_ix = i
            break
        else:pass
    for i, j in enumerate(h2):
        if j > crossover_pos:
            h2_ix = i
            break
        else:pass
    h1_new = np.append(h1[0:h1_ix], h2[h2_ix:])
    h2_new = np.append(h2[0:h2_ix], h1[h1_ix:])
    return(h1_new, h2_new)


def _temp(df, loci):
    females = df.query('sex == "F" and fec > 0')
    males = df.query('sex == "M"')
    if males.shape[0] == 0:
        return(males)
    elif females.shape[0] == 0:
        return(females)
    outs = []
    for _, f_row in females.iterrows():
        nr = f_row.copy()
        male = males.sample(1).iloc[0, :]
        nr = f_row.copy()
        for mf in range(f_row.fec):
            for loc in loci:
                lid = loc.idx + '_h{0!s}'
                num_recomb = np.random.poisson(
                        loc.recombination_rate * loc.basepairs * 2)
                if num_recomb == 0:
                    nr[lid.format(1)] = f_row[lid.format(random.choice("12"))]
                    nr[lid.format(2)] = male[lid.format(random.choice("12"))]
                else:
                    sex_xing = random.choice("MF")
                    h1m = male[lid.format(1)].copy()
                    h2m = male[lid.format(2)].copy()
                    h1f = f_row[lid.format(1)].copy()
                    h2f = f_row[lid.format(2)].copy()
                    for _ in range(num_recomb):
                        crossover_pos = random.randint(0,
                                loc.basepairs)
                        if sex_xing == "M":
                            h1m, h2m = recombination_locus(h1m, h2m,
                                    crossover_pos)
                        elif sex_xing == "F":
                            h1f, h2f = recombination_locus(h1f, h2f,
                                    crossover_pos)
                    nr[lid.format(1)] = random.choice([h1f, h2f])
                    nr[lid.format(2)] = random.choice([h1m, h2m])
            outs.append(nr)
    outdf = pd.concat(outs, axis=1, ignore_index=True).T
    return(outdf)


def recombination_fx(locus,
                     dfAdult,
                     recombination_rate,
                     basepairs):
    """calculate number of recombination events and rearranges haplotypes
    :TODO add for recombination map

    Parameters
    ---------
    locus: int
        number of loci
    dfAdult_mf : pandas dataframe
        dataframe containing new larval parasites
    dfAdult : pd df
        dataframe containing reproducing adults
    recombination_rate : float, list
        recombination rate for each locus
    basepairs : int, list
        length of each locus in basepairs

    Returns
    -------
    dfAdult_mf : pd df

    """
    lid = "locus_{0!s}"
    hosts = dfAdult.groupby('hostidx')
    loci = [Locus(lid.format(i), recombination_rate = j, basepairs = k)
            for i, j, k in zip(range(locus), recombination_rate, basepairs) if 
            j != 0]
    _temp2 = partial(_temp, loci=loci)
    dfAdult_mf = hosts.apply(_temp2)
    dfAdult_mf.reset_index(drop=True, inplace=True)
    return dfAdult_mf
