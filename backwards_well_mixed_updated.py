# -*- coding: utf-8 -*-
"""
Well-mixed sexual population simulation

Not using np.repeat to update the leaf counts
"""

import numpy as np
import random
from os import path
#from multiprocessing import Pool
#import sys

#N = int(sys.argv[1]) # population size
#s = float(sys.argv[2]) # selection coefficient
#N_forward_sim = int(sys.argv[3])
#N_back_sim_per_forward = int(sys.argv[4]) # number of simulation (forward + backward)
#nsample = int(sys.argv[5]) # number of individuals sampled for the backward part
#T_after_fix = int(sys.argv[6]) # time between fixation and sampling
#r = float(sys.argv[7]) # recombination rate


N = 10 ** 7
s = 0.05
N_forward_sim = 1
N_back_sim_per_forward = 5
nsample = 100
T_after_fix = 0
mutationSite = 50000000
startIndx = 1
endIndx = 100000000
step = 0

r_per_base = 1e-8

N_sim = N_forward_sim * N_back_sim_per_forward

def moving_average(a, n = 3, start_smooth = 100):
    a_start = a[:start_smooth]
    a_end = a[start_smooth:]
    ret = np.cumsum(a_end, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a_start, ret[n - 1:] / n))

def writeOutput(position, mutationCount, outputFileName):
    if path.exists(outputFileName) == False:
        output = open(outputFileName, "w")
        output.write("position\t" + "x\t" + "n\t" + "folded\n")
    output = open(outputFileName, "a")
    output.write(str(position) + "\t" + str(mutationCount) + "\t" + str(nsample) + "\t" + "0\n")
    output.close()

def recalcRecombination(position):
    relDistance = abs(mutationSite - position)
    genDistance = float(relDistance * r_per_base)
    newRecomb = 1 - (np.exp(-2 * genDistance)) 
    newRecomb *= 0.5
    return newRecomb

def updateStep(position):
    relDistance = abs(mutationSite - position)

    if relDistance < 10000:
        return 1000
    elif relDistance <= 100000:
        return 10000
    elif relDistance <= 1000000:
        return 100000
    else:
        return 1000000

def coalescent(parents, leaf_counts_offspring):
    '''
    Sum the leaf counts if two offsprings have the same parent (coalescence)
    return the set of parents (no duplicate), and leaf counts for the parents
    '''
    num_parents = len(parents)
    parents_set = np.array([])
    leaf_counts_parents = np.array([])
    for i in range(num_parents):
        parent = parents[i]
        leaf_count_offspring = leaf_counts_offspring[i]
        if len(parents_set) == 0:
            parents_set = np.array([parent])
            leaf_counts_parents = np.array([leaf_count_offspring])
        elif (parent ==  parents_set).all(1).any():
            idx = np.where((parents_set == parent).all(1))[0][0]
            leaf_counts_parents[idx] += leaf_count_offspring
        else:
            parents_set = np.append(parents_set, [parent], axis = 0)
            leaf_counts_parents = np.append(leaf_counts_parents, leaf_count_offspring)
    
    return parents_set, leaf_counts_parents

def backward_sim(forward_idx):
    SFS = np.zeros(nsample)
    n_selected_series = np.loadtxt(
            'well_mixed_forward_N={}_s={:.2f}_{}.txt'.format(
                N, s, forward_idx))
    
    T_sweep = len(n_selected_series)
    t_after_fix = 0
    leaf_counts = [1 for _ in range(nsample)]
    leaf_counts = np.array(leaf_counts).astype(np.int64)
    hist, bin_edges = np.histogram(leaf_counts, 
                                       bins = np.arange(1, nsample + 2))
    SFS += hist
    n_mut_tracked = nsample
    while t_after_fix < T_after_fix and n_mut_tracked > 1:
        T2 = np.random.exponential(
                N / ((n_mut_tracked - 1) * n_mut_tracked / 2))
        t_after_fix += T2
        
        hist, bin_edges = np.histogram(leaf_counts,
                                   bins = np.arange(1, nsample + 2))
        if t_after_fix < T_after_fix:
            SFS += hist * T2
    
            coal_inds = np.random.choice(range(n_mut_tracked), 
                                      size = 2, replace = False)
            leaf_counts_copy = leaf_counts
            leaf_counts_copy[coal_inds[0]] += leaf_counts[coal_inds[1]]
            leaf_counts_copy[coal_inds[1]] = 0
            leaf_counts = leaf_counts_copy[leaf_counts_copy > 0]
            n_mut_tracked -= 1


    individuals = np.array([[1, i] for i in range(n_mut_tracked)])
    individuals = individuals.astype(np.int64)
    n_tracked = len(individuals)
    leaf_counts = np.array(leaf_counts).astype(np.int64)

    t_backward_during_sweep = 1
    
    while t_backward_during_sweep < T_sweep and n_tracked > 1:
        t_backward_during_sweep += 1
        mut_types, idxs = individuals.T
        #starts at the second to the last number of the series
        Ne_parent = n_selected_series[-t_backward_during_sweep] 
        Ne = n_selected_series[-t_backward_during_sweep + 1]
        mut_types_parents = mut_types
        p_vals = np.random.random(n_tracked)
        rec_with_mut_idxs = np.where(p_vals < r * Ne / N)[0]
        mut_types_parents[rec_with_mut_idxs] = 1
        rec_with_wt_idxs = np.where(np.logical_and(
                p_vals > r * Ne / N, 
                p_vals < r))[0]
        mut_types_parents[rec_with_wt_idxs] = 0
        mut_types_parents = np.array(mut_types_parents).astype(np.int64)
        
        indices_parents = np.floor(np.random.random(n_tracked) 
        * (Ne_parent * mut_types_parents 
           + (N - Ne_parent) * (np.ones(n_tracked) - mut_types_parents)))
        
        individuals_parents = np.vstack((mut_types_parents, indices_parents)).T
        individuals_parents = (individuals_parents).astype(np.int64)
        individuals, leaf_counts = coalescent(individuals_parents, leaf_counts)
        n_tracked = len(individuals)
        hist, bin_edges = np.histogram(leaf_counts, 
                                           bins = np.arange(1, nsample + 2))
        SFS += hist

    while n_tracked > 1:
        T2 = np.random.exponential(
                N / ((n_tracked - 1) * n_tracked / 2))

        hist, bin_edges = np.histogram(leaf_counts,
                                   bins = np.arange(1, nsample + 2))
        SFS += hist * T2

        coal_inds = np.random.choice(range(n_tracked), 
                                  size = 2, replace = False)
        leaf_counts_copy = leaf_counts
        leaf_counts_copy[coal_inds[0]] += leaf_counts[coal_inds[1]]
        leaf_counts_copy[coal_inds[1]] = 0
        leaf_counts = leaf_counts_copy[leaf_counts_copy > 0]
        n_tracked -= 1
    
    SFS *= nsample
    return SFS

def calculateSFS():
    SFS_avg = np.zeros(nsample)
    for n_forward_sim in range(N_forward_sim):
        for n_back_sim_per_forward in range(N_back_sim_per_forward):
            SFS_avg += backward_sim(n_forward_sim) / N_sim
    return SFS_avg


mutationRate = r_per_base * 10
left = mutationSite
right = mutationSite
hi = endIndx
r = recalcRecombination(right)

while right < hi:
    step = updateStep(right)
    right += step
    left -= step

    SFS = calculateSFS()
    SFS *= mutationRate

    SFS_List = SFS.tolist()

    zeroProbability = 1 - sum(SFS_List)
    print(sum(SFS_List))
    print(zeroProbability)
    SFS_List.append(zeroProbability)

    mutations = []
    for x in range(1, nsample + 1):
        mutations.append(x)
    mutations.append(0)
    
    resultListLeft = random.choices(mutations, SFS_List, k = step)
    resultListRight = random.choices(mutations, SFS_List, k = step)

    targetLociLeft = random.sample(range(left, left + step), step)
    targetLociRight = random.sample(range(right - step, right), step)

    for i in range(step):
        writeOutput(targetLociLeft[i], resultListLeft[i], "output_well_mixed.txt")
        writeOutput(targetLociRight[i], resultListRight[i], "output_well_mixed.txt")

    r = recalcRecombination(right)
