# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 18:14:22 2021

@author: jim903
"""
N = 10 ** 6
Tfix = 0
s = 0.05
rlist = [10 ** (-7), 10 ** (-6), 10 ** (-5), 10 ** (-4)]

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
plt.xlabel(r'$f$')
plt.ylabel(r'$P(f)$')



for rind in range(len(rlist)):
    r = rlist[rind]
    well_mixed_SFS = np.loadtxt(
        'expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}_r={:.1e}.txt'.format(N, Tfix, s, r))

    f = np.arange(1/ len(well_mixed_SFS), 
              1 + 1 / len(well_mixed_SFS), 
              1 / len(well_mixed_SFS))



    plt.loglog(moving_average(f, 50, 70), moving_average(well_mixed_SFS, 50
               , 70), linewidth = 4, color = viridis_cmap(rind * 0.2 + 0.2)
                , label = r'$r =${:.0e}'.format(r))


plt.loglog(f, 1 / f ** 2 / s, label = r'$U_n / sf^2$', 
               linewidth = 3, linestyle = '-.', color = 'k')
plt.loglog(f, 2 * N / f, label = r'$2 N U_n / f$', 
               linewidth = 3, linestyle = '--', color = 'k')

plt.legend(fontsize = 'small', loc = 'lower left')



