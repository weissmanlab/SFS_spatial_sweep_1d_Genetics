# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 18:14:22 2021

@author: jim903
"""
Nlist = [2 * 10 ** 5, 4 * 10 ** 5, 6 * 10 ** 5, 8 * 10 ** 5, 
         10 ** 6, 2 * 10 ** 6, 4 * 10 ** 6]
s = 0.05
Tfix = 0

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

viridis_cmap = cm.get_cmap('viridis')

def moving_average(a, n = 3, start_smooth = 100):
    a_start = a[:start_smooth]
    a_end = a[start_smooth:]
    ret = np.cumsum(a_end, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a_start, ret[n - 1:] / n))

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 60, weight = 'bold')
plt.figure(figsize = (24, 18))

f = np.arange(10 ** (-6), 1 + 10 ** (-6), 10 ** (-6))

plt.loglog(f, 1 / f / s, label = r'$1 / sf$', 
               linewidth = 2, linestyle = '-.', color = 'k')

for Nind in range(len(Nlist)):
    N = Nlist[Nind]
    well_mixed_SFS = np.loadtxt(
        'expected_SFS_well_mixed_N={}.txt'.format(N))

    f = np.arange(1/ len(well_mixed_SFS), 
              1 + 1 / len(well_mixed_SFS), 
              1 / len(well_mixed_SFS))



    plt.loglog(moving_average(f, 200, 300), moving_average(well_mixed_SFS, 200
               , 300), linewidth = 2, color = viridis_cmap(Nind * 0.1)
                , label = r'$N =${:.0e}'.format(N))

    plt.loglog(f, 10 / f ** 2 / N / s, color = viridis_cmap(Nind * 0.1),
               linewidth = 2, linestyle = '--')

plt.legend(fontsize = 'medium', loc = 'lower left')