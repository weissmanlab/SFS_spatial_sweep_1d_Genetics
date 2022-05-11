'''
Todo : Figure3a - close up of Figure 2 around the intermediate f

Figure3b, plot P(f) s / ln(Nsf) for various m and s
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from labellines import labelLines
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple


def moving_average(a, n = 3, start_smooth = 100):
    a_start = a[:start_smooth]
    a_end = a[start_smooth:]
    ret = np.cumsum(a_end, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a_start, ret[n - 1:] / n))


from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 60, weight = 'bold')
Un = 1
r = 0
n_forward = 100
tfix = 0


# Change s - sample everywhere
n = 100000
nSFS = 1000

blue_cmap = cm.get_cmap('Blues')
red_cmap = cm.get_cmap('Reds')
grey_cmap = cm.get_cmap('Greys')

slist = np.arange(0.02, 0.07, 0.01)
mlist = np.arange(0.2, 0.55, 0.05)
Lrholist = [[500, 2000], [1000, 10000], [1000, 1000], [2000, 5000]]
tfinallist = [10000, 100000, 100000, 100000]

xpositions = [10**(-3), 10**(-3), 10**(-3), 6 * 10 ** (-5), 6 * 10 ** (-5)]
ypositions = [10 ** 5, 10 ** 4, 10 ** 3, 10 ** 4, 10 ** 3]

fitlist = []

f = np.arange(1, n + 1) / n
navg = 30
start_smooth = 100
f_short = moving_average(f, navg, start_smooth)

fit_range_ind = np.arange(200,1000)
fig, ax = plt.subplots(figsize = (24, 18))
vary_s_list = []

for sind in range(len(slist)):
    s = slist[sind]
    m = 0.25
    L = 500
    rho = 20000
    N = rho * L
    tfinal = 1000000

    fstart_ind = int(10 / (rho * np.sqrt(m * s)) * n)
    SFS = np.loadtxt(
        'backward_simulation_data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_avged.txt'.format(L, 
             rho, s, m, r, tfinal, n, tfix, nSFS))
    y_smooth = moving_average(SFS * s / np.log(N * s * f), 
                                navg, start_smooth)

    s_plot = ax.scatter(f_short[fstart_ind::30], y_smooth[fstart_ind::30], 
                 s = 250, marker = 'v', 
                 color = blue_cmap(100 - sind * 10))
    vary_s_list.append(s_plot)
    
vary_s_tuple = tuple(vary_s_list)

vary_m_list = []
for mind in range(len(mlist)):
    s = 0.05
    m = mlist[mind]
    L = 500
    rho = 20000
    N = rho * L
    tfinal = 1000000

    fstart_ind = int(10 / (rho * np.sqrt(m * s)) * n)    
    SFS = np.loadtxt(
        'backward_simulation_data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_avged.txt'.format(L, 
             rho, s, m, r, tfinal, n, tfix, nSFS))
    y_smooth = moving_average(SFS * s / np.log(N * s * f), 
                                navg, start_smooth)

    m_plot = ax.scatter(f_short[fstart_ind::30], y_smooth[fstart_ind::30], 
                 s = 250, marker = 'o', 
                 color = red_cmap(100 - mind * 10))
    vary_m_list.append(m_plot)
vary_m_tuple = tuple(vary_m_list)

vary_Lrho_list = []    
for Lrhoind in range(len(Lrholist)):
    s = 0.05
    m = 0.25
    L, rho = Lrholist[Lrhoind]
    N = L * rho
    tfinal = tfinallist[Lrhoind]
    fstart_ind = int(10 / (rho * np.sqrt(m * s)) * n)

    SFS = np.loadtxt(
        'backward_simulation_data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_avged.txt'.format(L, 
             rho, s, m, r, tfinal, n, tfix, nSFS))
    y_smooth = moving_average(SFS * s / np.log(N * s * f), 
                                navg, start_smooth)

    Lrho_plot = ax.scatter(f_short[fstart_ind::15], y_smooth[fstart_ind::15], 
                 s = 250, marker = '*', 
                 color = grey_cmap(100 - Lrhoind * 10))
    vary_Lrho_list.append(Lrho_plot)
vary_Lrho_tuple = tuple(vary_Lrho_list)

ax.loglog()
ax.set_xlim((3 * 10 ** -3, 10 ** -1))
ax.set_ylim((5 * 10, 5 * 10 ** 4))
ax.set_xlabel('Frequency, ' + r'$\boldmath{f}$', fontsize = 100)
ax.set_ylabel(r'$\boldmath{P(f) \cdot s / \ln(Nsf)}$', fontsize = 100)
ax.legend([vary_s_tuple, vary_m_tuple, vary_Lrho_tuple], [r'vary $s$, '
          + r'$m = 0.25, L = 500, \rho = 2\times 10^4$', 
          r'vary $m$, '
          + r' $s = 0.05, L = 500, \rho = 2\times 10^4$', 
          r'vary $\rho$ and $L$, ' 
          + r'$s = 0.05, m = 0.25$'], fontsize = 30, 
    handler_map={tuple: HandlerTuple(ndivide=None)})
plt.savefig('Figure3b_new.pdf', format = 'pdf', bbox_inches = 'tight')
