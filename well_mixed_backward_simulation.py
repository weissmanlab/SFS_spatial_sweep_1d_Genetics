# -*- coding: utf-8 -*-
"""
Well-mixed sexual population simulation

Not using np.repeat to update the leaf counts
"""

import numpy as np
#from multiprocessing import Pool
#import sys

#N = int(sys.argv[1]) # population size
#s = float(sys.argv[2]) # selection coefficient
#N_forward_sim = int(sys.argv[3])
#N_back_sim_per_forward = int(sys.argv[4]) # number of simulation (forward + backward)
#nsample = int(sys.argv[5]) # number of individuals sampled for the backward part
#T_after_fix = int(sys.argv[6]) # time between fixation and sampling
#r = float(sys.argv[7]) # recombination rate


N = 10 ** 5
s = 0.05
N_forward_sim = 100
N_back_sim_per_forward = 5
nsample = 10 ** 4
T_after_fix = 0
r = 0.0005

N_sim = N_forward_sim * N_back_sim_per_forward

def moving_average(a, n = 3, start_smooth = 100):
    a_start = a[:start_smooth]
    a_end = a[start_smooth:]
    ret = np.cumsum(a_end, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a_start, ret[n - 1:] / n))


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
            'forward simulation/well_mixed_forward_N={}_s={:.2f}_{}.txt'.format(
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

SFS_avg = np.zeros(nsample)

for n_forward_sim in range(N_forward_sim):
    for n_back_sim_per_forward in range(N_back_sim_per_forward):
        SFS_avg += backward_sim(n_forward_sim) / N_sim
    
        print(n_forward_sim * N_back_sim_per_forward + n_back_sim_per_forward)

np.savetxt('SFS_well_mixed_N={}_r={:.2e}_s={:.2e}_Navg={}_tfix={}.txt'.format(
        N, r, s, N_sim, T_after_fix), SFS_avg)
import matplotlib.pyplot as plt
freq = np.arange(1, nsample + 1) / nsample
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 60, weight = 'bold')
plt.figure(figsize = (24, 18))
plt.xlabel(r'$logit(f)$')
plt.ylabel(r'$P(f)$')

plt.semilogy(np.log(moving_average(freq, n = 3) / 
                  (1 - moving_average(freq, n = 3))), 
moving_average(SFS_avg, n = 3), linewidth = 3)
plt.semilogy(np.log(freq / (1 - freq)), 2 * N / freq, linestyle = '--', linewidth = 3)
plt.semilogy(np.log(freq / (1 - freq)), 1 / (s * freq ** 2), linewidth = 3)
plt.savefig('SFS_well_mixed_N={}_r={:.2e}_s={:.2e}_Navg={}_tfix={}.png'.format(N, r, s, N_sim, T_after_fix))


#if __name__ == '__main__':
#    p = Pool(5)
#    
#    SFS = np.sum(p.map(backward_sim, range(N_sim)), axis = 0) / N_sim
#    np.savetxt('expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}_r={:.2e}_nsim={}.txt'.format(N, T_after_fix, s, r, N_sim), SFS)
#
