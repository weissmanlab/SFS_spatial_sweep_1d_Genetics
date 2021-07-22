# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:46:16 2021

@author: jim903
"""

import numpy as np
import matplotlib.pyplot as plt
N = 2 * 10 ** 5
s = 0.05
N_sim = 500
nsample = 10 ** 4
T_after_fix = 0

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
        n_inds = nsample
        while t_after_fix < T_after_fix:
            t_after_fix += 1
            inds_new = np.random.randint(N, size = n_inds)
            inds_new_repeat = np.repeat(inds_new, leaf_counts, axis = 0)
            unique, leaf_counts = np.unique(inds_new_repeat, 
                                            return_counts = True)
            hist, bin_edges = np.histogram(leaf_counts,
                                           bins = np.arange(1, nsample + 2))
            SFS += hist
            n_inds = len(unique)
        while t < T_sweep - 1:
#            print(t)
            t += 1
            Ne = n_selected_series[-t - 1]
            inds_new = np.random.randint(Ne, size = n_inds)
            inds_new_repeat = np.repeat(inds_new, leaf_counts, axis = 0)
            unique, leaf_counts = np.unique(inds_new_repeat, 
                                            return_counts = True)
            hist, bin_edges = np.histogram(leaf_counts, 
                                           bins = np.arange(1, nsample + 2))
            SFS += hist
            n_inds = len(unique)
            
            
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
plt.loglog(f, 1 / f ** 2 / s, label = r'$P(f) = U_n / sf^2$')
plt.loglog(f, 2 * N / f, label = r'$P(f) = 2 N U_n / f$')
plt.legend(fontsize = 'medium', loc = 'upper right')
plt.savefig('expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}.png'.format(N, T_after_fix, s))

np.savetxt('expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}.txt'.format(N, T_after_fix, s), SFS)

