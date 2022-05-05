'''
Todo : Figure3a - close up of Figure 2 around the intermediate f

Figure3b, plot P(f) s / ln(Nsf) for various m and s
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from labellines import labelLines
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
plt.figure(figsize = (24, 18))


for sind in range(len(slist)):
    s = slist[sind]
    m = 0.25
    L = 500
    rho = 20000
    N = rho * L
    tfinal = 1000000

    fstart_ind = int(10 / (rho * np.sqrt(m * s)) * n)
    SFS = np.zeros(n)

    for i in range(n_forward):
        SFS += n * np.loadtxt(
        'backward_simulation_data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
                 rho, s, m, r, tfinal, n, tfix, nSFS, i))

    
    SFS /= n_forward
    y_smooth = moving_average(SFS * s / np.log(N * s * f), 
                                navg, start_smooth)

    plt.scatter(f_short[fstart_ind::15], y_smooth[fstart_ind::15], 
                 s = 250, marker = 'v', alpha = 0.5, 
                 color = blue_cmap(100 - sind * 5), 
                 label = '$s = ${:.2f}'.format(s))
for mind in range(len(mlist)):
    s = 0.05
    m = mlist[mind]
    L = 500
    rho = 20000
    N = rho * L
    tfinal = 1000000

    fstart_ind = int(10 / (rho * np.sqrt(m * s)) * n)    
    SFS = np.zeros(n)

    for i in range(n_forward):
        SFS += n * np.loadtxt(
        'backward_simulation_data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
                 rho, s, m, r, tfinal, n, tfix, nSFS, i))

    
    SFS /= n_forward
    y_smooth = moving_average(SFS * s / np.log(N * s * f), 
                                navg, start_smooth)

    plt.scatter(f_short[fstart_ind::15], y_smooth[fstart_ind::15], 
                 s = 250, marker = 'o', alpha = 0.5, 
                 color = red_cmap(100 - mind * 5), 
                 label = '$m = ${:.2f}'.format(m))

    
for Lrhoind in range(len(Lrholist)):
    s = 0.05
    m = 0.25
    L, rho = Lrholist[Lrhoind]
    N = L * rho
    tfinal = tfinallist[Lrhoind]
    fstart_ind = int(10 / (rho * np.sqrt(m * s)) * n)

    SFS = np.zeros(n)

    for i in range(n_forward):
        SFS += n * np.loadtxt(
        'backward_simulation_data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
                 rho, s, m, r, tfinal, n, tfix, nSFS, i))

    
    SFS /= n_forward
    y_smooth = moving_average(SFS * s / np.log(N * s * f), 
                                navg, start_smooth)

    plt.scatter(f_short[fstart_ind::15], y_smooth[fstart_ind::15], 
                 s = 250, marker = '*', alpha = 0.5, 
                 color = grey_cmap(100 - Lrhoind * 5), 
                 label = '$L = ${}'.format(L) + '$\rho = ${}'.format(rho))
   
plt.loglog()
plt.xlim((3 * 10 ** -3, 10 ** -1))
plt.ylim((5 * 10, 5 * 10 ** 4))
plt.xlabel('Frequency, ' + r'$\boldmath{f}$', fontsize = 100)
plt.ylabel(r'$\boldmath{P(f) \cdot s / \ln(Nsf)}$', fontsize = 100)
#plt.legend(fontsize = 'medium', loc = 'upper right')
plt.savefig('Figure3b_new.pdf', format = 'pdf', bbox_inches = 'tight')


########################################3333
#plt.figure(figsize = (24, 18))
#plt.plot(slist, fitlist, 'o',  ms = 60, label = 'simulation')
#slist2 = np.linspace(slist[0], slist[-1], 100)
#plt.plot(slist2, Un / 2 / np.array(slist2), label = '$k = U_n / 2s$', 
#           linewidth = 10)
#plt.xlabel('Selection coefficient, ' + r'$\boldmath{s}$', fontsize = 100)
#plt.ylabel('Fitting coefficient, ' + r'$\boldmath{k}$', fontsize = 100)
#plt.legend()
#plt.title('fitting AFS to $P(f) = k \ln( N s f) / f^2$')
#plt.savefig('Figure3b.pdf', format = 'pdf')