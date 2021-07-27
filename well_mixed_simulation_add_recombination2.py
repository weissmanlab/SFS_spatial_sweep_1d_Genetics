# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:46:16 2021

@author: jim903
"""

import numpy as np
import matplotlib.pyplot as plt
N = 10 ** 5
s = 0.05
N_sim = 500
nsample = 10 ** 3
T_after_fix = 0
r = 5 * 10 ** (-5)

n_sim = 0
SFS = np.zeros(nsample)

while n_sim < N_sim:
    n_selected = 1
    n_selected_series = [n_selected]
    while n_selected < N and n_selected > 0:
        n_selected = np.random.binomial(N, n_selected / N 
                                    + s * n_selected / N 
                                    * (1 - n_selected / N))
        n_selected_series.append(n_selected)
#    print(n_selected)
    if n_selected > 0:
        n_sim += 1
        print(n_sim)
        
        # coalescent simulation on a successful sweep data
        T_sweep = len(n_selected_series)
        t_after_fix = 0
        t = 0
        leaf_counts = [1 for _ in range(nsample)]
        leaf_counts_mut = leaf_counts
        leaf_counts_wt = []
        n_mut_inds = nsample
        n_wt_inds = 0
        while t_after_fix < T_after_fix:
            t_after_fix += 1
            inds_mut_new = np.random.randint(N, size = n_mut_inds)
            inds_mut_new_repeat = np.repeat(inds_mut_new, leaf_counts_mut, axis = 0)
            unique, leaf_counts_mut = np.unique(inds_mut_new_repeat, 
                                            return_counts = True)
            leaf_counts = np.append(leaf_counts_mut, leaf_counts_wt)
            hist, bin_edges = np.histogram(leaf_counts,
                                           bins = np.arange(1, nsample + 2))
            SFS += hist
            n_mut_inds = len(unique)
            n_wt_inds = 0
        individuals = np.array([[1, i] for i in range(n_mut_inds)])
        while t < T_sweep - 1:
            t += 1
            mut_types, idxs = individuals.T
            Ne = n_selected_series[-t - 1]
            Ne_current = n_selected_series[-t]
            
            mut_types_next = mut_types
            
            n_recom = np.random.poisson(N * r)
            idx_recom = np.random.randint(N, size = n_recom)
            
            for idx in idx_recom:
                idx_of_individuals = np.where(idx in idxs)[0]
                if len(idx_of_individuals) == 1:
                    p = np.random.random()
                    if mut_types_next[idx_of_individuals] < 1 and p < Ne_current / N:
                        mut_types_next[idx_of_individuals] = 1
                    elif mut_types_next[idx_of_individuals] > 0 and p < (N - Ne_current) / N:
                        mut_types_next[idx_of_individuals] = 0
                    
            
#            p_vals = np.random.random(len(individuals))
#            for i in range(len(individuals)):
#                if mut_types_next[i] < 1: # WT
#                    if p_vals[i] < r * Ne_current / N:
#                        mut_types_next[i] = 1
#                else:
#                    if p_vals[i] < r * (N - Ne_current) / N:
#                        mut_types_next[i] = 0
            
#            mut_types_next = mut_types
#            n_recom = np.random.poisson(len(individuals) * r)
#            idx_recom = np.random.randint(len(individuals), size = n_recom)
#            for i in range(n_recom):
#                p = np.random.random()
#                if p < Ne_current / N:
#                    mut_types_next[idx_recom[i]] = 1
#                else:
#                    mut_types_next[idx_recom[i]] = 0

            
#            mut_types_next = np.ones_like(mut_types)
##            print(t)
#            p_vals = np.random.random(len(individuals))
#            
#            zero_idxs = np.where(np.logical_and(mut_types == 1, 
#                                                p_vals < r * (N - Ne_current) / (N - 1)))
#            zero_idxs2 = np.where(np.logical_and(mut_types == 0, 
#                                                 p_vals > r * Ne_current / (N - 1)))
#            mut_types_next[zero_idxs] = 0
#            mut_types_next[zero_idxs2] = 0
            
#            n_mut_next = sum(mut_types_next > 0)
#            n_wt_next = len(individuals) - n_mut_next
            
            individuals2 = []
            for i in range(len(individuals)):
                if mut_types_next[i] > 0:
                    individuals2.append([1, np.random.randint(Ne)])
                else:
                    individuals2.append([0, np.random.randint(N - Ne)])
            
            individuals2_repeat = np.repeat(np.array(individuals2)
                                            , leaf_counts, axis = 0)
            unique, leaf_counts = np.unique(individuals2_repeat
                                            , return_counts = True, axis = 0)
            
            individuals = unique
            
            hist, bin_edges = np.histogram(leaf_counts, 
                                           bins = np.arange(1, nsample + 2))
            SFS += hist
            
        n_remains = len(individuals)
        while n_remains > 1:
            inds_new = np.random.randint(N, size = n_remains)
            inds_new_repeat = np.repeat(inds_new, leaf_counts, axis = 0)
            unique, leaf_counts = np.unique(inds_new_repeat, 
                                            return_counts = True)
            hist, bin_edges = np.histogram(leaf_counts,
                                           bins = np.arange(1, nsample + 2))
            SFS += hist
            n_remains = len(unique)
            
            
SFS /= N_sim
SFS *= nsample

def moving_average(a, n = 3, start_smooth = 100):
    a_start = a[:start_smooth]
    a_end = a[start_smooth:]
    ret = np.cumsum(a_end, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a_start, ret[n - 1:] / n))

f = np.arange(1, nsample + 1) / nsample
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 60, weight = 'bold')
plt.figure(figsize = (24, 18))
plt.xlabel(r'$f$', fontsize = 75)
plt.ylabel(r'$P(f)$', fontsize = 75)

plt.loglog(moving_average(f, 20, 100), moving_average(SFS, 20, 100), linewidth = 2)
plt.loglog(f, (1 + 2 * N * r) / f ** 2 / s, label = r'$P(f) = U_n(1 + 2Nr) / sf^2$')
plt.loglog(f, 2 * N / f, label = r'$P(f) = 2 N U_n /f$')
plt.legend(fontsize = 'medium', loc = 'upper right')
plt.savefig('expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}_r={:.1e}_uptick.png'.format(N, T_after_fix, s, r))

#np.savetxt('expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}_r={:.2e}_uptick.txt'.format(N, T_after_fix, s, r), SFS)

