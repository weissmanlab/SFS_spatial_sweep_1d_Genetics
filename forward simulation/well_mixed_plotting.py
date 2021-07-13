# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 18:14:22 2021

@author: jim903
"""
N = 10 ** 7
s = 0.05

import numpy as np
import matplotlib.pyplot as plt

well_mixed_SFS = np.loadtxt('expected_SFS_well_mixed.txt')


f = np.arange(1/ len(well_mixed_SFS), 
              1 + 1 / len(well_mixed_SFS), 
              1 / len(well_mixed_SFS))
plt.rc('text', usetex=True)

def moving_average(a, n = 3, start_smooth = 100):
    a_start = a[:start_smooth]
    a_end = a[start_smooth:]
    ret = np.cumsum(a_end, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a_start, ret[n - 1:] / n))

plt.figure(figsize = (24, 18))
plt.loglog(moving_average(f, 200, 1000), moving_average(well_mixed_SFS, 200, 1000), linewidth = 2)

plt.loglog(f, 2 * np.log(N) / f, label = r'$P(f) = 2 \ln(N) U_n / f$')
plt.loglog(f, 1 / f ** 2 / N / s, label = r'$P(f) = U_n / N sf^2$')

plt.legend(fontsize = 'medium', loc = 'upper right')