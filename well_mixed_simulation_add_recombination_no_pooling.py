# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:46:16 2021

@author: jim903
"""

import numpy as np
import sys

N = int(sys.argv[1]) # population size
s = float(sys.argv[2]) # selection coefficient
N_sim = int(sys.argv[3]) # number of simulation (forward + backward)
nsample = int(sys.argv[4]) # number of individuals sampled for the backward part
T_after_fix = int(sys.argv[5]) # time between fixation and sampling
r = float(sys.argv[6]) # recombination rate

#N = 10 ** 4
#s = 0.05
#N_sim = 10000
#nsample = 10 ** 3
#T_after_fix = 0
#r = 5 * 10 ** (-5)


def forward_sim():
    n_selected = 1
    n_selected_series = [n_selected]
    while n_selected < N and n_selected > 0:
        n_selected = np.random.binomial(N, n_selected / N 
                                        + s * n_selected / N
                                         * (1 - n_selected / N))
        n_selected_series.append(n_selected)
    return n_selected_series

def backward_sim(idx):
    SFS = np.zeros(nsample)
    n_selected_final = 0
    while n_selected_final < N:
        n_selected_series = forward_sim()
        n_selected_final = n_selected_series[-1]
    
    T_sweep = len(n_selected_series)
    t_after_fix = 0
    leaf_counts = [1 for _ in range(nsample)]
    leaf_counts = np.array(leaf_counts).astype(np.int64)
    hist, bin_edges = np.histogram(leaf_counts, 
                                       bins = np.arange(1, nsample + 2))
    SFS += hist
    n_mut_tracked = nsample
    while t_after_fix < T_after_fix:
        t_after_fix += 1
        idx_mut_new = np.random.randint(N, size = n_mut_tracked)
        idx_mut_new_repeat = np.repeat(idx_mut_new, leaf_counts, axis = 0)
        unique, leaf_counts = np.unique(idx_mut_new_repeat, 
                                        return_counts = True)
        hist, bin_edges = np.histogram(leaf_counts, 
                                       bins = np.arange(1, nsample + 2))
        SFS += hist
        n_mut_tracked = len(unique)


    individuals = np.array([[1, i] for i in range(n_mut_tracked)])
    individuals = individuals.astype(np.int64)
    n_tracked = len(individuals)
    leaf_counts = np.array(leaf_counts).astype(np.int64)

    t_backward_during_sweep = 1
    
    while t_backward_during_sweep < T_sweep:
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
        individuals_parents = []
        for i in range(n_tracked):
            if mut_types_parents[i] == 0:
                individuals_parents.append([0, np.random.randint(N - Ne_parent)])
            else:
                individuals_parents.append([1, np.random.randint(Ne_parent)])
        individuals_parents = np.array(individuals_parents).astype(np.int64)
        individuals_parents_repeat = np.repeat(individuals_parents, 
                                               leaf_counts, axis = 0)
        unique, leaf_counts = np.unique(individuals_parents_repeat, 
                                        return_counts = True, axis = 0)
        individuals = unique
        n_tracked = len(individuals)
        hist, bin_edges = np.histogram(leaf_counts, 
                                           bins = np.arange(1, nsample + 2))
        SFS += hist

    while n_tracked > 1:
        individuals = np.random.randint(N, size = n_tracked)
        T2 = np.random.exponential(
                N / ((n_tracked - 1) * n_tracked / 2))

        hist, bin_edges = np.histogram(leaf_counts,
                                   bins = np.arange(1, nsample + 2))
        SFS += hist * T2
    
        coal_idxs = np.random.choice(range(n_tracked), 
                              size = 2, replace = False)

        individuals[coal_idxs[0]] = individuals[coal_idxs[1]] 

        individuals_parents = np.repeat(individuals, leaf_counts, axis = 0)

        unique, leaf_counts = np.unique(individuals_parents, 
                                    axis = 0, return_counts = True)
        individuals = unique
        n_tracked = len(individuals)
    SFS *= nsample
    if np.mod(idx, 5) == 0:
        np.savetxt('progress_for_N={}_Tfix={}_s={:.2f}_r={:.2e}_{}.txt'.format(
                N, T_after_fix, s, r, idx), SFS[:5])    
    return SFS

if __name__ == '__main__':
    SFS = np.zeros(nsample)
    for idx in range(N_sim):
        SFS += backward_sim(idx)
    SFS /= N_sim
    np.savetxt('expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}_r={:.2e}_nsim={}.txt'.format(N, T_after_fix, s, r, N_sim), SFS)

#f = np.arange(1, nsample + 1) / nsample
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif', size = 60, weight = 'bold')
#plt.figure(figsize = (24, 18))
#plt.xlabel(r'$f$', fontsize = 75)
#plt.ylabel(r'$P(f)$', fontsize = 75)
#
#plt.loglog(moving_average(f, 40, 30), moving_average(SFS, 40, 30), linewidth = 2)
#plt.loglog(f, (1 + 2 * N * r) / f ** 2 / s, label = r'$P(f) = U_n(1 + 2Nr) / sf^2$')
#plt.loglog(f, 2 * N / f, label = r'$P(f) = 2 N U_n /f$')
#plt.legend(fontsize = 'medium', loc = 'upper right')
#plt.savefig('expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}_r={:.1e}_uptick.png'.format(N, T_after_fix, s, r))


