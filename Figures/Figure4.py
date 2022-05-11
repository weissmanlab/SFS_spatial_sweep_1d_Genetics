import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax
def moving_average(a, n = 3, start_smooth = 100):
    a_start = a[:start_smooth]
    a_end = a[start_smooth:]
    ret = np.cumsum(a_end, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a_start, ret[n - 1:] / n))

def power_law_fit(x, a):
    ''' Use this to fit f^-2 (intermediate f) '''
    return a * x ** (-2)

from scipy.optimize import curve_fit
from matplotlib import cm

from labellines import labelLines
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 60, weight = 'bold')
plt.figure(figsize = (24, 18))
plt.xlabel(r'Frequency, $f$', fontsize = 75)
plt.ylabel(r'Number of alleles, $P(f)$', fontsize = 75)

L = 500
rho = 20000
N = L * rho
s = 0.05
m = 0.25
tfinal = 1000000
Un = 1
r = 0
n_forward = 100
tfix = 0


# Change m - sample everywhere
n = 100000
nSFS = 1000

tfixlist = [1000, 10000, 100000]
viridis_cmap = cm.get_cmap('viridis')

fitlist = []

f = np.arange(1, n + 1) / n
navg = 30
start_smooth = 100
f_short = moving_average(f, navg, start_smooth)

fit_range_ind = np.arange(200,1000)


freq_file = open(
   'forward_simulation_data/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(
           L, rho, s, m, tfinal))
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
#xpositions = [6 * 10 ** (-5), 10 ** (-3), 10 ** (-2)]
xpositions = [5 * 10 ** (-5), 8 * 10 ** (-4), 1.5 * 10 ** (-2)]

plt.text(2 * 10 ** -4, 200, r'$\boldmath{f = t / N + 1 / (\rho v)}$', fontsize = 80)
flat, = plt.loglog(f_short, np.ones(len(f_short)) * L / v, linewidth = 5, 
           linestyle = 'dotted'
              , label = r'$ U_n L / v$', color = '#d55e00')
neutral, = plt.loglog(f_short, 
           2 * Un * N * np.ones(len(f_short)) / f_short,
           label = r'$2 N U_n / f$', linestyle = '-.',
           linewidth = 5, color = '#cc79a7')
bsc, = plt.loglog(f_short, 
           Un / (2 * s * f_short ** 2) * np.log(N * s * f_short), 
           label = r'$U_n \ln(Nsf) / (2 s f^2)$', 
           linestyle = '-.', linewidth = 5, color = 'b')

xvals = [5 * 10 ** (-5), 0.1, 3 * 10 ** -5]
labelLines(plt.gca().get_lines(), xvals = xvals, fontsize = 75)

for tind in range(len(tfixlist)):
    tfix = tfixlist[tind]
    # Find v from lines

    SFS = np.zeros(n)

    for i in range(n_forward):
        SFS += n * np.loadtxt(
        'backward_simulation_data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
                 rho, s, m, r, tfinal, n, tfix, nSFS, i))

    
    SFS /= n_forward


    plt.loglog(f_short, 
                 moving_average(SFS, navg, start_smooth), 
                 linewidth = 5, 
                 color = viridis_cmap(tind * 0.4), 
                 label = '$t =${:.2e}'.format(tfix))
    
    plt.vlines(tfix / N + 1 / (v * rho), 10 ** 3, 10 ** 11, linestyle = 'dotted',
               linewidth = 5, color = viridis_cmap(tind * 0.4))
    plt.text(xpositions[tind], 5 * 10 ** 11, 
             r'$\boldmath{t = }$' + '{:.1e}'.format(tfix), 
             color = viridis_cmap(tind * 0.4), fontsize = 60)    
    popt, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS[fit_range_ind])
    fitlist.append(popt[0])


plt.xlim((2 * 10 ** -6, 0.99))
plt.savefig('Figure4.pdf', format = 'pdf', bbox_inches = 'tight')
