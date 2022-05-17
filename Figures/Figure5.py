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
rlist = [10 ** -2, 5 * 10 ** (-4), 10 ** -4, 10 ** (-5)]
n = 10000
n_sim_well_mixed = 10000
nback_sim_1d_list = [1000, 1000, 1000, 1000, 3000]


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 60, weight = 'bold')

viridis_cmap = cm.get_cmap('viridis')


def moving_average(a, n = 3, start_smooth = 100):
    a_start = a[:start_smooth]
    a_end = a[start_smooth:]
    ret = np.cumsum(a_end, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a_start, ret[n - 1:] / n))
def logit(x):
    return np.log(x / (1 - x))



freq_file = open(
       'forward_simulation_data/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, rho, s, m, tfinal))
lines = np.loadtxt(freq_file, dtype=np.int64)
tf_real = len(lines)
v_list = []
for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
    line1 = lines[t]
    line2 = lines[t + 1]
    psum1 = sum(line1) / rho
    psum2 = sum(line2) / rho
    v_list.append(psum2 - psum1)
    
v = np.average(v_list)

f_BSC_left = np.linspace(10 ** -2.8, 10 ** -1.2)
f_BSC_right = np.linspace(0.9, 0.995)
f_BSC_left_wellmixed = np.linspace(10 ** -2.8, 0.5)
f_BSC_right_wellmixed = np.linspace(0.5, 0.995)
f_uniform = np.linspace(10 ** -1, 0.9)

for rind in range(len(rlist)):
    plt.figure(figsize = (24, 18))
    plt.xlabel(r'$logit(f)$')
    plt.ylabel(r'$P(f)$')
    r = rlist[rind]

    plt.title(r'$r = ${:.2e}'.format(r))
    nback_sim_1d = nback_sim_1d_list[rind]
    smooth_start_1d = 70
    window_size_1d = 100
    f_1d = np.arange(1 / n, 1 + 1 / n, 1 / n)
    f_short_1d = moving_average(f_1d, window_size_1d, smooth_start_1d)
    f_short_wellmixed = moving_average(f_1d, window_size_1d, smooth_start_1d)
    SFS_1d = np.loadtxt('backward_simulation_data/expected_SFS_L=' 
    + '{}_N={}_s={:.3f}_m={:.2f}_r={:.2e}_nsample={}_t_after_fix={}_Nback={}_avged.txt'.format(L, 
             rho, s, m, r, n, T_after_fix, nback_sim_1d))
    plt.semilogy(logit(f_short_1d), moving_average(SFS_1d, 
                  window_size_1d, smooth_start_1d), linewidth = 3, 
        label = 'r = {:.2e}'.format(r), color = 'k')
    if r < 10 ** -3:
        plt.semilogy(logit(f_uniform), 
                      (1 + 2 * N * r) * L / v * np.ones(len(f_uniform)), 
                  linewidth = 5, linestyle = 'dotted', 
                  label = r'$U_n(1 + 2Nr) L / v$', 
                  color = '#ff7f00')
        plt.semilogy(logit(f_BSC_left), 
                  (1 + 2 * N * r) * np.log(N * s * f_BSC_left) / s / f_BSC_left ** 2 / 2, 
                  linewidth = 5, linestyle = '-.', 
                  label = r'$U_n(1 + 2Nr) \ln(Nsf) / 2 s f^2$', 
                  color = '#377eb8')

        plt.semilogy(logit(f_BSC_right), 
                  (2 * N * r) * np.log(N * s * (1 - f_BSC_right)) / s / (1 - f_BSC_right) ** 2 / 2, 
                  linewidth = 5, linestyle = '-.', 
                  color ='#4daf4a', 
                  label = r'$U_n N r \ln(Ns(1 - f)) / s (1 - f)^2$')

    if r < 10 ** -2:
        plt.semilogy(logit(f_BSC_left_wellmixed), 
                  (1 + 2 * N * r) / s / f_BSC_left_wellmixed ** 2, 
                  linewidth = 5, alpha = 0.5, linestyle = '-.', 
                  label = r'$U_n(1 + 2Nr) / s f^2$', 
                  color = '#377eb8')

        plt.semilogy(logit(f_BSC_right_wellmixed), 
                  (2 * N * r) / s / (1 - f_BSC_right_wellmixed) ** 2, 
                  linewidth = 5, linestyle = '-.', alpha = 0.5, 
                  color ='#4daf4a', 
                  label = r'$2 U_n N r / s (1 - f)^2$')
    SFS_well_mixed = np.loadtxt(
            'backward_simulation_data/'
            + 'expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}_r={:.2e}_nsample={}_nsim={}.txt'.format(
                    N, T_after_fix, s, r, n, n_sim_well_mixed))
    plt.semilogy(logit(f_short_wellmixed), moving_average(SFS_well_mixed,
                 window_size_1d, smooth_start_1d), linewidth = 5, alpha = 0.4, 
    color = 'k')
    plt.semilogy(logit(f_short_1d), 2 * N / f_short_1d, 
                 linewidth = 5, linestyle = '--', label = r'$2 N U_n / f$', 
                 color = '#f781bf')
    plt.vlines(logit(np.log(N * s) * r / s), 10 ** 6, 10 ** 11, 
               label = r'$f = \ln (Ns) r / s$')
    plt.vlines(-logit(np.log(N * s) * r / s), 10 ** 6, 10 ** 11, 
               label = r'$f = 1 - \ln (Ns) r / s')

    plt.savefig('Figure5_{}.pdf'.format(rind), format = 'pdf', bbox_inches = 'tight')
