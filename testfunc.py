#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 16:17:00 2016

@author: scott
"""
import numpy as np
from math import pi, exp, sqrt
import math
import subprocess, re, random
from sklearn.metrics.pairwise import euclidean_distances
import argparse
import copy
prev=[.1,.1,.05]
hostpopsize=[100,100,100]


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
def migration_matrix(villages):
    '''calculated a migration matrix for use with ms (hudson 2000) from euclidian
    distances. The value of 2Nm is weighted by the distance from the next population
    as an exponential random variable. The highest this can be is 2Nm/1'''  
    #distance_m is list [1000,1500] such that distance_m[0] is between 1 & 2
    #distance_m[1] is between 1 & 3 etc ...
    #villages is int 3; such that len(distance_m) == ((villages)*(villages-1)/2)
    m = 0.0001
    Ne = 10000
    distance_m = [1000,1500,500]
    if villages > 4: 
        raise ValueError("only handles 4 villages ATM")
    elif villages < 4:
        if len(distance_m) != ((villages)*(villages-1)/2): raise ValueError("there are not adequate pairwise comparisons")
        mig = []
        for meters in distance_m:
            mig.append((m)/(np.random.exponential(meters)))
        if villages == 2:
            M1=4*Ne*mig[0]
            return "{} {} {} {}".format(0,M1,0,M1) #mig_matrix is symmetrical 
        elif villages == 3: 
            M1=4*Ne*mig[0]
            M2=4*Ne*mig[1]
            M3=4*Ne*mig[2]      
            return "{} {} {} {} {} {} {} {} {}".format(0,M1,M2,M1,0,M3,M2,M3,0) #mig_matrix is symmetrical
        elif villages == 4:
            M1=4*Ne*mig[0]
            M2=4*Ne*mig[1]
            M3=4*Ne*mig[2]
            M4=4*Ne*mig[3]       
            M5=4*Ne*mig[4]
            M6=4*Ne*mig[5]
            return "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(0,M1,M2,M3,M1,0,M4,M5,M2,M4,0,M6,M3,M5,M6,0)

def ms_outcall(worm_popsize):
    '''external call to ms (Hudson 2000) using the migration matrix calculated
above from the euclidian distances. This function will then read in the stdout from
ms. The continious location of Seg given by ms is transformed to discrete coordinates. 
It then builds a list of lists from the discrete coordinates which is then pointer
referenced in the dictionary recording the metapopulation information
    metapop_init = [100,100] #2 metapops
    theta = [metapop1,metapop2]    
    ''' 
    theta = [5,6,2]
    theta_anc = 12
    ta = 0.05
    t12 = 0.05
    t23 = 0.05
    t34 = 0.05
    if len(worm_popsize) == 1: #1 village
        mscmd = "ms {} 1 -t {} -eN {} {} > temp_file".format(worm_popsize[0],theta[0],ta,float(theta_anc/theta[0]))        
    else: #set up for 2 villages
        num_subpops = len(worm_popsize) #-I num_pops
        total_inds = sum(worm_popsize) # how many
        sub_pop = " ".join(map(str,worm_popsize))#-I X i j ...
        #t = 0.05 #2000gens/4*N0  #merge time   
        #-ej the villages split at 2000 gens in the past which is 2000 years
        if len(worm_popsize) == 2:        
            mscmd = "ms {} 1 -t {} -I {} {} -n 1 {} -n 2 {} -ma {} -ej {} 1 2 -eN {} {} > temp_file".format(total_inds,theta[0],num_subpops,sub_pop,1,float(theta[1])/theta[0],migration_matrix(len(worm_popsize)),t12,ta,theta_anc)       
        #after the split the pops are half the size of ancestral and share migrants at rate -ma
        #pop continues to go back until ancestral time    
        elif len(worm_popsize)==3:
            mscmd = "ms {} 1 -t {} -I {} {} -n 1 {} -n 2 {} -n 3 {} -ma {} -ej {} 1 2 -ej {} 2 3 -eN {} 1 > temp_file".format(total_inds,theta[0],num_subpops,sub_pop,1,float(theta[1])/theta[0],float(theta[2])/theta[0],migration_matrix(len(worm_popsize)),t12,t23,ta,theta_anc)       
        elif len(worm_popsize)==4:
            mscmd = "ms {} 1 -t {} -I {} {} -n 1 {} -n 2 {} -n 3 {} -n4 {} -ma {} -ej {} 1 2 -ej {} 2 3 -ej {} 3 4 -eN {} 1> temp_file".format(total_inds,theta[0],num_subpops,sub_pop,1,float(theta[1])/theta[0],float(theta[2])/theta[0],float(theta[3])/theta[0],migration_matrix(len(worm_popsize)),t12,t23,t34,ta,theta_anc)       
    print mscmd
    proc = subprocess.Popen(mscmd, shell=True)
    proc.wait()

    #parses ms output    
    hap_pop = []  
    with open("temp_file",'r') as ms:
        for line in ms:
            if line.startswith("positions"):             
                for line in ms:            
                    hap = line.rstrip("\n")                    
                    hap_pop.append([m.start() for m in re.finditer("1", hap)])    
    return hap_pop
    
def trans_init(prev,hostpopsize):
    '''initializes locations for above infections within
    prevalence = [.8,.2]
    hostpop_size = [100,100]    
    '''
    #from sklearn.metrics.pairwise import euclidean_distances
    sigma = 2
    dispersal = 4*pi*sigma**2

    size = 1

    mu = 100
    metan = 1
    transmission_mat = AutoVivification()     
    for i,j in zip(prev,hostpopsize):    
        x1 = np.random.negative_binomial(size,size/float((size+mu)),round(i*j)) # number of successes, prob of success (size/size+mu),number of values to return
        x2 = np.random.negative_binomial(size,size/float((size+mu)),round(i*j)) # number of successes, prob of success (size/size+mu),number of values to return
        X = np.vstack((x1,x2)).T    
        dist = euclidean_distances(X,X) 
        for pop in range(0,len(x1)):
            transmission_mat["meta_{}".format(metan)]["pop_{}".format(pop+1)] = ["{},{}".format(x1[pop],x2[pop]), [1 if item < dispersal else 0 for item in dist[pop,:]]]
        metan += 1
    return transmission_mat, dispersal    
  
def worm_burden(prev,hostpopsize):
    '''
    num_pops = [15,25] #number of infections in meta1 and meta2; prevalance * population
    mu = [5,1] #avg_burden, average number of adult female worms in infections;
    size = [1,50] #size = dispersion, size parameter for negative binomial distribution
    hap_pop is from ms_outcall    
    '''  
    # number of adult worms per infection    
    num_inf = []
    for i,j in zip(prev,hostpopsize):
        num_inf.append(round(i*j))
    pops = len(hostpopsize)

    mu = np.random.uniform(1,5,pops) #return num equal to

    size = np.random.uniform(1,50,pops) #return num equal to
    #worm burden    
    pop_init=[]
    for i,j,k in zip(mu,size,num_inf):
        wb_burden = np.random.negative_binomial(i,i/float(i+j),k) # number of successes, prob of success (size/size+mu),number of values to return
        pop_init.append(np.array(wb_burden).tolist())
    
    #make populations with haplotypes
    worm_popsize = []
    for meta in pop_init:    
        worm_popsize.append(sum(meta))
    hap_pop = ms_outcall(worm_popsize)  
        
   #meta_popdict[meta1][pop][age][list]
    meta_popdict = AutoVivification()            
    pop = 1 #initial value
    meta = 1 #initial value
    k = 0 #initial start of hap_pop
    for metapop in pop_init:
        for wb_a1 in metapop:
            j = 0 #initializes         
            meta_popdict["meta_" + str(meta)]["pop_" + str(pop)]["A_1"]=[]
            while j < wb_a1:
                meta_popdict["meta_" + str(meta)]["pop_" + str(pop)]["A_1"].append([np.random.uniform(), k])
                j += 1 #counts the A1 in population
                k += 1 #counts the haps in hap_pop[]
            pop += 1 #advances the pop counter
        meta += 1 #advances the meta counter
        pop = 1 #resets pop counter for new meta
    transmission_mat,dispersal=trans_init(prev,hostpopsize)
    return meta_popdict, hap_pop, transmission_mat #dictionary, max(max(hap_pop))
    
#-------    
    
def maturation(meta_popdict, hap_pop, month):
    '''step all individuals forward, after deaths. Kill them then move, then births'''    

    # density dependent mortality    

#    b = 1.0/K #K is carrying capacity
#    S = 0.99
#    a = (S*b)/exp(-1) #where S is maximum survival
#    mort_A = a * sum(len(A_1 ... A_8)) * exp(-b * sum(len(A_1 ... A_8))) #Ricker fx 

# fixed mortality actually prob of surviving  

    mort_A = .88 #year
    mort_J = .86 #month
    mort_M = .90 #month
    
# density dep fecundity   

#    b = 1.0/K #K is carrying capacity
#    S = 0.99
#    a = (S*b)/exp(-1) #where S is maximum survival
#    fecund = 20 * a * sum(len(A_1 ... A_8)) * exp(-b * sum(len(A_1 ... A_8))) #Ricker fx 

# fixed fecundity     

    fecund = 20   

# mutation

    bp = 13000
    pmut = 7.6E-8
    frac_gen = 0.125 #gen is 8 months so 1 month is .125 of a generation
    
    rzero_freq = []
    hap_freq = []
   #since this month to month 
    if month%12 is 0: #this denotes 1 year has passed so adults mature to next age class 
        for mpop in meta_popdict.keys(): #villages
            for npop in meta_popdict[mpop].keys(): #inf
                meta_popdict[mpop][npop]["A_8"] = random.sample(meta_popdict[mpop][npop]["A_7"],int(round(len(meta_popdict[mpop][npop]["A_7"])*mort_A)))
                meta_popdict[mpop][npop]["A_7"] = random.sample(meta_popdict[mpop][npop]["A_6"],int(round(len(meta_popdict[mpop][npop]["A_6"])*mort_A)))
                meta_popdict[mpop][npop]["A_6"] = random.sample(meta_popdict[mpop][npop]["A_5"],int(round(len(meta_popdict[mpop][npop]["A_5"])*mort_A)))
                meta_popdict[mpop][npop]["A_5"] = random.sample(meta_popdict[mpop][npop]["A_4"],int(round(len(meta_popdict[mpop][npop]["A_4"])*mort_A)))
                meta_popdict[mpop][npop]["A_4"] = random.sample(meta_popdict[mpop][npop]["A_3"],int(round(len(meta_popdict[mpop][npop]["A_3"])*mort_A)))
                meta_popdict[mpop][npop]["A_3"] = random.sample(meta_popdict[mpop][npop]["A_2"],int(round(len(meta_popdict[mpop][npop]["A_2"])*mort_A)))
                meta_popdict[mpop][npop]["A_2"] = random.sample(meta_popdict[mpop][npop]["A_1"],int(round(len(meta_popdict[mpop][npop]["A_1"])*mort_A)))
                # first year adults, later add from different age classes with maturity an evolvable trait            
                meta_popdict[mpop][npop]["A_1"] = random.sample(meta_popdict[mpop][npop]["J_12"],int(round(len(meta_popdict[mpop][npop]["J_12"])*mort_J)))               
                # juvenille
                meta_popdict[mpop][npop]["J_12"] = random.sample(meta_popdict[mpop][npop]["J_11"],int(round(len(meta_popdict[mpop][npop]["J_11"])*mort_J)))
                meta_popdict[mpop][npop]["J_11"] = random.sample(meta_popdict[mpop][npop]["J_10"],int(round(len(meta_popdict[mpop][npop]["J_10"])*mort_J)))
                meta_popdict[mpop][npop]["J_10"] = random.sample(meta_popdict[mpop][npop]["J_9"],int(round(len(meta_popdict[mpop][npop]["J_9"])*mort_J)))
                meta_popdict[mpop][npop]["J_9"] = random.sample(meta_popdict[mpop][npop]["J_8"],int(round(len(meta_popdict[mpop][npop]["J_8"])*mort_J)))
                meta_popdict[mpop][npop]["J_8"] = random.sample(meta_popdict[mpop][npop]["J_7"],int(round(len(meta_popdict[mpop][npop]["J_7"])*mort_J)))
                meta_popdict[mpop][npop]["J_7"] = random.sample(meta_popdict[mpop][npop]["J_6"],int(round(len(meta_popdict[mpop][npop]["J_6"])*mort_J)))
                meta_popdict[mpop][npop]["J_6"] = random.sample(meta_popdict[mpop][npop]["J_5"],int(round(len(meta_popdict[mpop][npop]["J_5"])*mort_J)))
                meta_popdict[mpop][npop]["J_5"] = random.sample(meta_popdict[mpop][npop]["J_4"],int(round(len(meta_popdict[mpop][npop]["J_4"])*mort_J)))
                meta_popdict[mpop][npop]["J_4"] = random.sample(meta_popdict[mpop][npop]["J_3"],int(round(len(meta_popdict[mpop][npop]["J_3"])*mort_J)))
                meta_popdict[mpop][npop]["J_3"] = random.sample(meta_popdict[mpop][npop]["J_2"],int(round(len(meta_popdict[mpop][npop]["J_2"])*mort_J)))
                meta_popdict[mpop][npop]["J_2"] = random.sample(meta_popdict[mpop][npop]["J_1"],int(round(len(meta_popdict[mpop][npop]["J_1"])*mort_J)))
                # microfilaria                
                meta_popdict[mpop][npop]["MF_12"] = random.sample(meta_popdict[mpop][npop]["MF_11"],int(round(len(meta_popdict[mpop][npop]["MF_11"])*mort_M)))
                meta_popdict[mpop][npop]["MF_11"] = random.sample(meta_popdict[mpop][npop]["MF_10"],int(round(len(meta_popdict[mpop][npop]["MF_10"])*mort_M)))
                meta_popdict[mpop][npop]["MF_10"] = random.sample(meta_popdict[mpop][npop]["MF_9"],int(round(len(meta_popdict[mpop][npop]["MF_9"])*mort_M)))
                meta_popdict[mpop][npop]["MF_9"] = random.sample(meta_popdict[mpop][npop]["MF_8"],int(round(len(meta_popdict[mpop][npop]["MF_8"])*mort_M)))
                meta_popdict[mpop][npop]["MF_8"] = random.sample(meta_popdict[mpop][npop]["MF_7"],int(round(len(meta_popdict[mpop][npop]["MF_7"])*mort_M)))
                meta_popdict[mpop][npop]["MF_7"] = random.sample(meta_popdict[mpop][npop]["MF_6"],int(round(len(meta_popdict[mpop][npop]["MF_6"])*mort_M)))
                meta_popdict[mpop][npop]["MF_6"] = random.sample(meta_popdict[mpop][npop]["MF_5"],int(round(len(meta_popdict[mpop][npop]["MF_5"])*mort_M)))
                meta_popdict[mpop][npop]["MF_5"] = random.sample(meta_popdict[mpop][npop]["MF_4"],int(round(len(meta_popdict[mpop][npop]["MF_4"])*mort_M)))
                meta_popdict[mpop][npop]["MF_4"] = random.sample(meta_popdict[mpop][npop]["MF_3"],int(round(len(meta_popdict[mpop][npop]["MF_3"])*mort_M)))
                meta_popdict[mpop][npop]["MF_3"] = random.sample(meta_popdict[mpop][npop]["MF_2"],int(round(len(meta_popdict[mpop][npop]["MF_2"])*mort_M)))
                meta_popdict[mpop][npop]["MF_2"] = random.sample(meta_popdict[mpop][npop]["MF_1"],int(round(len(meta_popdict[mpop][npop]["MF_1"])*mort_M)))
                #count A_1                
                rzero_freq.append([i[0]for i in meta_popdict[mpop][npop]["A_1"]])
                #reset A_1 since they are new from Juv12
                for subl in meta_popdict[mpop][npop]["A_1"]:
                    subl[0] = random.random()
                #MF_1
                mf1 = []                
                for i in range(1,8):
                    for j in meta_popdict[mpop][npop]["A_{}".format(i)]:
                        mf1.append([j]*np.random.poisson(fecund)) 
                mf = sum(mf1,[])
                #set up mutation        
                num_muts = np.random.binomial(len(mf), bp * pmut * frac_gen)
                if num_muts != 0:
                    mf,hap_pop = mutation(mf,hap_pop,num_muts)
                meta_popdict[mpop][npop]["MF_1"] = mf
                hap_freq.append([i[1]for i in mf])
                
    else: #a year has not passed on months, juveniles and MF move to next age class
        for mpop in meta_popdict.keys():        
            for npop in meta_popdict[mpop].keys(): #inf
                # first year adults, later add from different age classes with maturity an evolvable trait            
                meta_popdict[mpop][npop]["A_1"].append(random.sample(meta_popdict[mpop][npop]["J_12"],int(round(len(meta_popdict[mpop][npop]["J_12"])*mort_J))))
                # juvenilles                
                meta_popdict[mpop][npop]["J_12"] = random.sample(meta_popdict[mpop][npop]["J_11"],int(round(len(meta_popdict[mpop][npop]["J_11"])*mort_J)))
                meta_popdict[mpop][npop]["J_11"] = random.sample(meta_popdict[mpop][npop]["J_10"],int(round(len(meta_popdict[mpop][npop]["J_10"])*mort_J)))
                meta_popdict[mpop][npop]["J_10"] = random.sample(meta_popdict[mpop][npop]["J_9"],int(round(len(meta_popdict[mpop][npop]["J_9"])*mort_J)))
                meta_popdict[mpop][npop]["J_9"] = random.sample(meta_popdict[mpop][npop]["J_8"],int(round(len(meta_popdict[mpop][npop]["J_8"])*mort_J)))
                meta_popdict[mpop][npop]["J_8"] = random.sample(meta_popdict[mpop][npop]["J_7"],int(round(len(meta_popdict[mpop][npop]["J_7"])*mort_J)))
                meta_popdict[mpop][npop]["J_7"] = random.sample(meta_popdict[mpop][npop]["J_6"],int(round(len(meta_popdict[mpop][npop]["J_6"])*mort_J)))
                meta_popdict[mpop][npop]["J_6"] = random.sample(meta_popdict[mpop][npop]["J_5"],int(round(len(meta_popdict[mpop][npop]["J_5"])*mort_J)))
                meta_popdict[mpop][npop]["J_5"] = random.sample(meta_popdict[mpop][npop]["J_4"],int(round(len(meta_popdict[mpop][npop]["J_4"])*mort_J)))
                meta_popdict[mpop][npop]["J_4"] = random.sample(meta_popdict[mpop][npop]["J_3"],int(round(len(meta_popdict[mpop][npop]["J_3"])*mort_J)))
                meta_popdict[mpop][npop]["J_3"] = random.sample(meta_popdict[mpop][npop]["J_2"],int(round(len(meta_popdict[mpop][npop]["J_2"])*mort_J)))
                meta_popdict[mpop][npop]["J_2"] = random.sample(meta_popdict[mpop][npop]["J_1"],int(round(len(meta_popdict[mpop][npop]["J_1"])*mort_J)))
                # microfilaria                
                meta_popdict[mpop][npop]["MF_12"] = random.sample(meta_popdict[mpop][npop]["MF_11"],int(round(len(meta_popdict[mpop][npop]["MF_11"])*mort_M)))
                meta_popdict[mpop][npop]["MF_11"] = random.sample(meta_popdict[mpop][npop]["MF_10"],int(round(len(meta_popdict[mpop][npop]["MF_10"])*mort_M)))
                meta_popdict[mpop][npop]["MF_10"] = random.sample(meta_popdict[mpop][npop]["MF_9"],int(round(len(meta_popdict[mpop][npop]["MF_9"])*mort_M)))
                meta_popdict[mpop][npop]["MF_9"] = random.sample(meta_popdict[mpop][npop]["MF_8"],int(round(len(meta_popdict[mpop][npop]["MF_8"])*mort_M)))
                meta_popdict[mpop][npop]["MF_8"] = random.sample(meta_popdict[mpop][npop]["MF_7"],int(round(len(meta_popdict[mpop][npop]["MF_7"])*mort_M)))
                meta_popdict[mpop][npop]["MF_7"] = random.sample(meta_popdict[mpop][npop]["MF_6"],int(round(len(meta_popdict[mpop][npop]["MF_6"])*mort_M)))
                meta_popdict[mpop][npop]["MF_6"] = random.sample(meta_popdict[mpop][npop]["MF_5"],int(round(len(meta_popdict[mpop][npop]["MF_5"])*mort_M)))
                meta_popdict[mpop][npop]["MF_5"] = random.sample(meta_popdict[mpop][npop]["MF_4"],int(round(len(meta_popdict[mpop][npop]["MF_4"])*mort_M)))
                meta_popdict[mpop][npop]["MF_4"] = random.sample(meta_popdict[mpop][npop]["MF_3"],int(round(len(meta_popdict[mpop][npop]["MF_3"])*mort_M)))
                meta_popdict[mpop][npop]["MF_3"] = random.sample(meta_popdict[mpop][npop]["MF_2"],int(round(len(meta_popdict[mpop][npop]["MF_2"])*mort_M)))
                meta_popdict[mpop][npop]["MF_2"] = random.sample(meta_popdict[mpop][npop]["MF_1"],int(round(len(meta_popdict[mpop][npop]["MF_1"])*mort_M)))
                #count A_1                
                rzero_freq.append([i[0]for i in meta_popdict[mpop][npop]["A_1"]])
                #reset A_1
                for subl in meta_popdict[mpop][npop]["A_1"]:
                    subl[0] = random.random()
                # MF_1
                mf1 = []                
                for i in range(1,8):
                    for j in meta_popdict[mpop][npop]["A_{}".format(i)]:
                        mf1.append([j]*np.random.poisson(fecund)) 
                mf = sum(mf1,[])
                #set up mutation        
                num_muts = np.random.binomial(len(mf), bp * pmut * frac_gen)
                if num_muts != 0:
                    mf,hap_pop = mutation(mf,hap_pop,num_muts)
                meta_popdict[mpop][npop]["MF_1"] = mf
                hap_freq.append([i[1]for i in mf])

    #this will calculate for the enitre meta and each inf
    return meta_popdict, hap_pop, {x:rzero_freq.count(x) for x in rzero_freq},{y:hap_freq.count(y) for y in hap_freq}  
    
def mutation(mf,hap_pop,num_muts):
   '''this is run every time the prob of a mutation is true, updates seq_base'''   
   #calculate the number of mutations expected
   mut_mf = [random.randrange(len(mf)) for i in range(num_muts)] #choose random index in mf list
       #add new sequence to hap_pop
   for m in mut_mf:                 
       new_hap = copy.copy(hap_pop[mf[m][1]])
       new_allele = (max([max(a) for a in hap_pop])) + 1
       new_hap.append(new_allele)
       hap_pop.append(new_hap)
       mf[m][1] = len(hap_pop)-1
   return mf, hap_pop

   
   
   
   
meta_popdict,hap_pop,transmission_mat=worm_burden(prev,hostpopsize)   
   
def transmission(transmission_mat,meta_p,meta_popdict,dispersal):    
    '''this has a continious prob using a gillepsie algorithm to get the waiting time between transmission events'''    
    #this executes if tranmission is true     
    #pick a random donating pop:
    dpop = random.choice(meta_popdict[meta_p].keys())    
    while [len(meta_popdict[meta_p][dpop][i]) is 0 for i in meta_popdict[meta_p][dpop] if "MF" in i].count(True) is 12:
        dpop = random.choice(meta_popdict[meta_p].keys())
    r = random.uniform(0,1)    
    #new infection
    if r > 1 - (float(1)/len(transmission_mat[meta_p][dpop][1])):
        newpop = "pop_{}".format(len(meta_popdict[meta_p].keys())+1)
        meta_popdict[meta_p][newpop]["J_1"] = []
        rmf = int(round(random.uniform(.51,12)))
        while len(meta_popdict[meta_p][dpop]["MF_{}".format(rmf)]) is 0:
            rmf = int(round(random.uniform(.51,12)))
        dmf = random.choice(meta_popdict[meta_p][dpop]["MF_{}".format(rmf)])      
        meta_popdict[meta_p][newpop]["J_1"].append(dmf)
        print "new"
        transmission_mat = new_infection(transmission_mat, meta_p,dpop,newpop,dispersal) #find new transmission position
    #reinfection
    else:     
        rpop = random.choice([i for i, x in enumerate(transmission_mat[meta_p][dpop][1]) if x == 1]) + 1 #choose a random value that is 1 as the receiving pop        
        if meta_popdict[meta_p]["pop_{}".format(rpop)]["J_1"] is not list:
            meta_popdict[meta_p]["pop_{}".format(rpop)]["J_1"] = []
        rmf = int(round(random.uniform(.51,12)))
        while len(meta_popdict[meta_p][dpop]["MF_{}".format(rmf)]) is 0:
            rmf = int(round(random.uniform(.51,12)))
        dmf = random.choice(meta_popdict[meta_p][dpop]["MF_{}".format(rmf)])      
        meta_popdict[meta_p]["pop_{}".format(rpop)]["J_1"].append(dmf)
    return meta_popdict, transmission_mat

def new_infection(transmission_mat, meta_p,dpop,newpop,dispersal):
    '''this is run every time a new individual is infected to rebuild the transmission prob matrix, updates'''        
    x1n = np.random.uniform(-1*sqrt(dispersal), sqrt(dispersal)) 
    x2n = np.random.uniform(-1*sqrt(dispersal), sqrt(dispersal))    
    transmission_mat[meta_p][newpop] = ["{},{}".format(int(transmission_mat[meta_p][dpop][0].split(",")[0]) + x1n, int(transmission_mat[meta_p][dpop][0].split(",")[1]) + x2n),[]]
    for tpop in range(1,len(transmission_mat[meta_p].keys())):
        dist = sqrt((float(transmission_mat[meta_p]["pop_{}".format(tpop)][0].split(",")[0]) - x1n)**2 + (float(transmission_mat[meta_p]["pop_{}".format(tpop)][0].split(",")[1]) - x2n)**2)        
        if dist < dispersal:
            dist = 1
        else: 
            dist = 0
        transmission_mat[meta_p][newpop][1].append(dist)
        transmission_mat[meta_p]["pop_{}".format(tpop)][1].append(dist)
    transmission_mat[meta_p][newpop][1].append(1)   
    #if pop dies it is N in the transmission_mat, it is removed from meta_popdict
     
    return transmission_mat

