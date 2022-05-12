# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 18:14:22 2021

@author: jim903
"""
rho = 5000
L = 2000
N = rho * L
m = 0.25

n_forward = 100
T_after_fix = 0
tfinal = 100000

s = 0.05
rlist = [5 * 10 ** -4]
# rlist = [10 ** (-2), 5 * 10 ** (-4), 10 ** -4, 10 ** (-5)]
Nforwardlist = [1000, 1000, 1000, 3000]
N_forward_1d = 1000
n = 10000
n_sim_well_mixed = 1000
import numpy as np

for rind in range(len(rlist)):
    r = rlist[rind]
    SFS_1d = np.zeros(n)

    for i in np.arange(0, n_forward):
        SFS_1d += n * np.loadtxt(
        'backward_simulation_data/expected_SFS_L=' 
        + '{}_N={}_s={:.3f}_m={:.2f}_r={:.2e}_nsample={}_t_after_fix={}_Nback={}_Nforw={}.txt'.format(L, 
                 rho, s, m, r, n, T_after_fix, N_forward_1d, i))    
        
    SFS_1d /= n_forward 
    np.savetxt('backward_simulation_data/expected_SFS_L=' 
    + '{}_N={}_s={:.3f}_m={:.2f}_r={:.2e}_nsample={}_t_after_fix={}_Nback={}_avged.txt'.format(L, 
             rho, s, m, r, n, T_after_fix, N_forward_1d), SFS_1d)
