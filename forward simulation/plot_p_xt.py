# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:58:46 2021

@author: jim903
"""

import numpy as np
import matplotlib.pyplot as plt

L = 500 # number of demes
N = 20000 # deme capacity
s = 0.05 # selection coef
m = 0.25 # migration rate
tfinal = 1000000 # sweep time
Un = 1 # rate of neutral mutation
nbase = 1000 # sample size
N_SFS = 1 # number of coalescent simulation we run.
L_start = 450 # sampling range 
L_end = 500

cmap = plt.get_cmap('viridis')
ancestry = np.loadtxt('ancestry_tree_L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_nsample={}_Lstart={}_Lend={}_sample_uniform_0.txt'.format(L,
           N, s, m, tfinal, nbase, L_start, L_end))

plt.figure(figsize = (10, 10))
plt.xlim((0, L))
plt.title('L = {}, N = {}, sample range = [{}, {}]'.format(L, N, L_start, L_end))
plt.xlabel('x')
plt.ylabel('counts')
for i in np.arange(1, len(ancestry), 100):
    plt.hist(ancestry[-i], alpha = 0.5, color = cmap(i / len(ancestry)))

plt.savefig('ancentor_loc_L = {}, N = {}, sample range = [{}, {}].png'.format(L, N, L_start, L_end))
    