#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import pandas as pd
from itertools import combinations
from libsequence.polytable import simData
from libsequence.summstats import polySIM
from libsequence.summstats import garudStats
from libsequence.summstats import ld
from libsequence.summstats import lhaf
from libsequence.windows import Windows
from libsequence.parallel import scheduler_init
init = scheduler_init(10)
from libsequence.fst import fst
#plotting
#from .plotting import plot_allele_trace #for sel_trace_allele
#from .plotting import plot_allele_frequency #for SFS and jSFS
#from .plotting import plot_hapnetwork #for mitochondrial networks on MF
#from .plotting import plot_pairwise

def site_freqspec_fx(dfworm, mon, mf, locus):
    '''calculates the site frequency spectrum. If >1 pop also does jSFS

    Parameters
    ---------
    dfworm : df
        padas df with meta, h1/h2 etc ...
    mon : int
        month of logTime
    mf : array
        list of indexes corresponding to stage == M
    locus : str
        locus to query
    Returns
    -------
    plotSFS : fig
        plot of site freq spectrum

    '''
    print("sfs")
    #calculate SFS
    mfsfs = np.vstack([dfworm.h1[locus][mf], dfworm.h2[locus][mf]])
    freqsum = np.sum(mfsfs, axis=0)
    sfs = np.unique(np.sort(freqsum), return_counts = True)
#    sfsx = sfs[0]
#    sfsy = sfs[1]
    #np.histogram(mfsfs, bins = 10)
#    plot_allele_frequency(sfs)
    print("{}".format(sfs))
    #pass to sel_trace
    sel_trace_fx(dfworm, freqsum, locus, mf)

    return(None)

def sel_trace_fx(dfworm, freqsum, locus, mf):
    '''calculates the change in allele frequency through time for mutations with
    fitness affects

    Parameters
    ---------
    dfworm : df
        padas df with meta, h1/h2 etc ...
    freqsum : array
        counts of all allele as sfs
    mf : array
        list of indexes corresponding to stage == M
    locus : str
        locus to query
    Returns
    -------
    plotSelTrace : fig
        plot of allele freq trace for selected sites
    '''
    print("seltrace")
    if bool(dfworm.sel):
        #records trace of selected alleles
        nind = mf.shape[0]
        ftrace = np.where(dfworm.sel[locus + 'F'] > 0)[0]
        f_freq = freqsum[ftrace] / (2.0 * nind)
        strace = np.where(dfworm.sel[locus + 'S'] > 0)[0]
        s_freq = freqsum[strace] / (2.0 * nind)
#        plot_allele_trace(f_freq, s_freq)
        #avg fitness
        avgF = dfworm.meta.loc[mf,"fitF"].mean()
        avgS = dfworm.meta.loc[mf, "fitS"].mean()
        #dfworm.meta.loc[mf,"fitF"].hist()
        #dfworm.meta.loc[mf,"fitS"].hist()
    else: pass
    return(None)

def haplotype_net_fx(dfworm):
    '''constructs haplotype network table for plotting, vcf2hap

    Parameters
    ---------
    dfworm : df
        padas df with meta, h1/h2 etc ...

    Returns
    -------
    hapnetwork : fig
        haplotype network
    HapConfig : list
        haplotype configuration
    '''
    print("hapnet")
    sample_size = 100
#    K = np.vstack({tuple(row) for row in dfworm.h1['0']})
    mf = dfworm.meta[(dfworm.meta["stage"] == "M")].index.values
    mfhaps = dfworm.meta.ix[mf].groupby("village").apply(lambda y: y.sample(sample_size).index.values)
    #haplotype configuration
    H = []
    M = []
    K = []
    C = []
    for v in range(len(mfhaps)):
        hap_sample = dfworm.h1['0'][[(mfhaps[v])]]
        uniqhaps = np.array([np.array(x) for x in set(tuple(x) for x in hap_sample)])
        hapfreq = np.array([len(hap_sample[np.all(hap_sample==x, axis=1)]) for x in uniqhaps],dtype=int)
        n = sum(hapfreq)
        C_freq, C_count = np.unique(hapfreq,return_counts=True)
        C_v = np.zeros(n)
        C_v[C_freq-1] = C_count #full hap configuration
        C.append(C_v)
        #haplotype diversity
        H.append(1 - sum([(((i+1)/float(n))**2) * c for i,c in enumerate(C_v)]))
        M.append(max(C_v))
        K.append(sum(C_v))

    #haplotype networks
    #TO DO: calculate haplotype networks
    #vcf2networks, fitchi, pegas/adgenet in R ???
    #plot_hapnetwork(hapnetobj)
    return(None)

def figs2vcf_fx(dfworm):
    '''Takes worm obj and outputs vcf file

    Parameters
    ---------
    dfworm : df
        padas df with meta, h1/h2 etc ...

    Returns
    -------
    figs2vcf : file
        vcf formatted file
    '''
    print("figs2vcf")


    return(None)

def figs2scikit_fx(dfworm, sample_size, vill, locus):
    '''takes worm object and output format to load into scikit-allel

    Parameters
    ---------
    dfworm : df
        padas df with meta, h1/h2 etc ...
    sample_size : int
        number of worms to use for stats calculation
    vill : int
    locus : str

    Returns
    -------
    scikit object
        calculate stats using modules in scikit-allel
    '''
    print("scikit")
    import allel
    sample_size = 30
    vill = 1
    #groupby: village == vill, stage == "M", hostidx
    mf = dfworm.meta[(dfworm.meta["stage"] == "M") & (dfworm.meta["village"] == vill)].index.values
    mf_subsample = mf[np.random.choice(mf.shape[0], sample_size, replace = False)]
    haps1 = dfworm.h1[locus][mf_subsample]
    haps2 = dfworm.h2[locus][mf_subsample]
    #import into allel
    positions = dfworm.pos[locus]
    haparray_ind = np.dstack((haps1, haps2))
    haparray_pop = np.concatenate(haparray_ind, axis = 1)
    h = allel.HaplotypeArray(haparray_pop)
    #g = h.to_genotypes(ploidy=2)
    ac = h.count_alleles()
#    mf_mask = dfworm.meta.ix[mf].groupby("hostidx").size() < sample_size
#    if any(mf_mask):
#        mf = dfworm.meta.ix[mf][~dfworm.meta["hostidx"].isin(mf_mask[mf_mask].index.values)].index.values
#    mf_pairs = dfworm.meta.ix[mf].groupby("hostidx").apply(lambda y: y.sample(sample_size).index.values)
#    hostidx = [h for h in mf_pairs.keys()]

    return(None)

def pairwise_div_fx(dfworm, mon, vill, basepairs, sample_size):
    '''calculates the pairwise FST between hosts

    Parameters
    ---------
    dfworm : df
        padas df with meta, h1/h2 etc ...
    mon : int
        month of logTime
    vill : int
        which village
    basepairs : list
        list of locus sizes in basepairs
    sample_size : int
        number of worms to use for stats calculation

    Returns
    -------
    fst_t : list
        slatkin 1991 calculation of Fst
    dxy : list
        pairwise diversity between 2 populations
    da : list
        pairwise diversity normalize by pi
    popgenTable : df
        summary stats for all host in each village at time mon
    '''
    print("pairwise")
    fst_t = []
    dxy = []
    da = []
    theta = []
    tajimasd = []
    thetapi = []
    fulid = []
    fulidstar = []
    fulif = []
    fulifstar = []
    hprime = []
    numexternalmutations = []
    nummutations = []
    numpoly = []
    numsingletons = []
    #pi = []
#####
    #groupby: village == vill, stage == "M", hostidx
    mf = dfworm.meta[(dfworm.meta["stage"] == "M") & (dfworm.meta["village"] == vill)].index.values
    mf_mask = dfworm.meta.ix[mf].groupby("hostidx").size() < sample_size
    if any(mf_mask):
        mf = dfworm.meta.ix[mf][~dfworm.meta["hostidx"].isin(mf_mask[mf_mask].index.values)].index.values
        import ipdb; ipdb.set_trace()
    mf_pairs = dfworm.meta.ix[mf].groupby("hostidx").apply(lambda y: y.sample(sample_size).index.values)
    hostidx = [h for h in mf_pairs.keys()]
#####
    for loc in range(1, len(basepairs)):
        locus = str(loc)
        #print sfs and allele traces
        site_freqspec_fx(dfworm, mon, mf, locus)
        print("back2hostgen")
        seq = float(basepairs[loc])
        pos = dfworm.pos[locus] / seq
        for host in hostidx:
            #pylibseq
            pop1 = np.vstack([dfworm.h1[locus][mf_pairs[host]], dfworm.h2[locus][mf_pairs[host]]])
            gtpop1 = [''.join(str(n) for n in y) for y in pop1]
            #host stats
            sdpop1 = simData()
            sdpop1.assign_sep(pos, gtpop1)
            pspop1 = polySIM(sdpop1)
            print("hoststats {}".format(host))
            theta.append(pspop1.thetaw() / seq)
            tajimasd.append(pspop1.tajimasd())
            thetapi.append(pspop1.thetapi() / seq)
            #pi.append(sum([(2.0*i*(size-i)) / (size*(size-1)) for i in np.sum(pop1)]))
            fulid.append(pspop1.fulid())
            fulidstar.append(pspop1.fulidstar())
            fulif.append(pspop1.fulif())
            fulifstar.append(pspop1.fulifstar())
            hprime.append(pspop1.hprime())
            numexternalmutations.append(pspop1.numexternalmutations())
            nummutations.append(pspop1.nummutations())
            numpoly.append(pspop1.numpoly())
            numsingletons.append(pspop1.numsingletons())
            #this is super slow
            size = pop1.shape[0]
#            wallsb = pspop1.wallsb()
#            wallsbprime = pspop1.wallsbprime()
#            wallsq = pspop1.wallsq()
#            garudStats_t = garudStats(sdpop1) #garud 2015
#            lhaf_t = lhaf(sdpop1,1) #1-HAF is most common; Ronen 2016
#            ld_t = ld(sdpop1, haveOutgroup = False, mincount = .05, maxDist = 5000)

#############this part is just slow, maybe parallel ??
        for hostX, hostY in combinations(hostidx, 2):
            #print("fst")
            popX = np.vstack([dfworm.h1[locus][mf_pairs[hostX]], dfworm.h2[locus][mf_pairs[hostX]]])
            popY = np.vstack([dfworm.h1[locus][mf_pairs[hostY]], dfworm.h2[locus][mf_pairs[hostY]]])
            sdfst = simData()
            geno_fst = np.vstack([popX, popY])
            gtpop_fst = [''.join(str(n) for n in y) for y in geno_fst]
            sdfst.assign_sep(pos, gtpop_fst)
            size = [popX.shape[0], popY.shape[0]]
            f1 = fst(sdfst, size)
            fst_t.append(f1.slatkin())
            #good summary stats for ABC analysis
#            fst_trad = f1.hsm()
#            gst = f1.hbk()
#            shared_sites = f1.shared()
#            private_sites = f1.priv()
#            fixed_sites = f1.fixed()
            #dxy, sample sizes must be equal
            sizedxy = min(size)
            #is this total pairwise or product of 2 pops? Diploid?
            #print("dxy")
            dxy.append(sum([sum((x + y) == 1) for x, y in zip(popX, popY)]) / float((sizedxy)) / seq) #approximate
            #dxy_t = sum([sum((x + y) == 1) for x in popX for y in popY]) / (sizedxy**2.0) #more precise, but really slow
#############

        print("da")
        da_t = [dx - np.mean(pi_xy) for dx, pi_xy in zip(dxy, combinations(thetapi,2))]
        da.append(da_t)

    popgenTable = pd.DataFrame({"month" : [mon] * len(hostidx),
                                "hostidx" : hostidx,
                                "theta" : theta,
                                "tajD" : tajimasd,
                                "thetapi" : thetapi,
                                "fulid" : fulid,
                                "fulidstar" : fulidstar,
                                "fulif" : fulif,
                                "fulifstar" : fulifstar,
                                "hprime" : hprime,
                                "numexternalmutations" : numexternalmutations,
                                "nummutations" : nummutations,
                                "numpoly" : numpoly,
                                "numsingletons" : numsingletons
                                })
    popgenTable = popgenTable.loc[:, ['month', 'hostidx', 'theta', 'tajD', 'thetapi', 'fulid',
                              'fulidstar', 'fulif', 'fulifstar', 'hprime', 'numexternalmutations', 'nummutations', 'numpoly',
                              'numsingletons']]
    popgenTable = popgenTable.round(5)
    popgenTable.to_csv("popgenHostTable_v{}_m{}.csv".format(vill, mon))
    #maybe these are nested lists then to numpy then to ndarray.tofile('fstmatrix_)
#    #writes full matrix for analysis
#    write2file.csv("fstmatrix_v{}_m{}".format(vill, mon))
#    write2file.csv("dxymatrix_v{}_m{}".format(vill, mon))
#    write2file.csv("damatrix_v{}_m{}".format(vill, mon))
#    #surface plot of pairwise
#    plot_pairwise(fst)
#    plot_pairwise(dxy)
#    plot_pairwise(da)
    return(fst_t, dxy, da)

def villpopgen_fx(dfworm, outstats, vill, mon):
    '''calculates popgen statistic for each village
    [basepairs, sample_size, window_length, num_windows, wb2vcf, figs2scikit]
    Parameters
    ---------
    dfworm : df
        padas df with meta, h1/h2 etc ...
    mon : int
        month of logTime
    vill : int
        which village
    outstats : list
        contains [basepairs, sample_size, window_length, num_windows, wb2vcf, figs2scikit]

    Returns
    -------

    '''
    print("villagePopgen")
    #vary these for sample size and locus size info
    sample_size = outstats[1] #number of sequences to sample
    num_windows = outstats[3]
    basepairs = outstats[0]

    ##example of 1 window of 5kb from the start of a 13kb locus
    #window_length = outstats[2] #5000
    #win_size = window_length / basepairs[loc] #5000 / 13000.0
    #num_windows = 1
    #win_step = 1.0 / num_windows #1.0
    #w = Windows(sd, window_size = 0.38, step_len = 1, starting_pos = 0., ending_pos = 1.0)
    #w = Windows(sd, window_size = 0.38, step_len = 0.38, starting_pos = 0.38, ending_pos = 0.76)

    thetaA_t = []
    thetaJ_t = []
    thetaM_t = []
    tajD_t = []

    for loc in range(1, len(basepairs)):
        locus = str(loc)
        seq = float(basepairs[loc])
        pos = dfworm.pos[locus] / seq
#        win_size = float(outstats[2]) / float(basepairs[loc])
#        win_step = 1.0/num_windows
        adpop = dfworm.meta[(dfworm.meta.village == vill) & (dfworm.meta.stage == "A")].sample(sample_size).index.values
        jvpop = dfworm.meta[(dfworm.meta.village == vill) & (dfworm.meta.stage == "J")].sample(sample_size).index.values
        mfpop = dfworm.meta[(dfworm.meta.village == vill) & (dfworm.meta.stage == "M")].sample(sample_size).index.values

        mfgeno = np.vstack([dfworm.h1[locus][mfpop], dfworm.h2[locus][mfpop]])
        jvgeno = np.vstack([dfworm.h1[locus][jvpop], dfworm.h2[locus][jvpop]])
        adgeno = np.vstack([dfworm.h1[locus][adpop], dfworm.h2[locus][adpop]])
        mfgeno1 = [''.join(str(n) for n in y) for y in mfgeno]
        jvgeno1 = [''.join(str(n) for n in y) for y in jvgeno]
        adgeno1 = [''.join(str(n) for n in y) for y in adgeno]

        #haplotype 1
        sdad1 = simData()
        sdad1.assign_sep(pos, adgeno1)
        sdjv1 = simData()
        sdjv1.assign_sep(pos, jvgeno1)
        sdmf1 = simData()
        sdmf1.assign_sep(pos, mfgeno1)
        psad1 = polySIM(sdad1)
        psjv1 = polySIM(sdjv1)
        psmf1 = polySIM(sdmf1)

        #stats
        thetaA_t.append(psad1.thetaw() / seq)
        thetaJ_t.append(psjv1.thetaw() / seq)
        thetaM_t.append(psmf1.thetaw() / seq)
        tajD_t.append(psmf1.tajimasd() / seq)

    #call paired functions
    fst, dxy, da = pairwise_div_fx(dfworm, mon, vill, basepairs, sample_size) #returns mean fst_slatkin
    fst_slat_t = np.mean(fst)
    dxy_t = np.mean(dxy)
    da_t = np.mean(da)

    return(np.mean(thetaA_t), np.mean(thetaJ_t), np.mean(thetaM_t), fst_slat_t,
           da_t, dxy_t, np.mean(tajD_t))