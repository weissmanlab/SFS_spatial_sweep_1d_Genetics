import numpy as np
import matplotlib.pyplot as plt

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
N = 20000
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

colorlist = ['y', 'r', 'c', 'g', 'm', 'b', 'k']
fitlist = []

f = np.arange(1, n + 1) / n
navg = 30
start_smooth = 100
f_short = moving_average(f, navg, start_smooth)

fit_range_ind = np.arange(200,1000)


freq_file = open(
   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s, m, tfinal))
lines = np.loadtxt(freq_file, dtype=np.int64)
tf_real = len(lines)
v_list = []
for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
    line1 = lines[t]
    line2 = lines[t + 1]
    psum1 = sum(line1) / N
    psum2 = sum(line2) / N
    v_list.append(psum2 - psum1)

v = np.average(v_list)
xpositions = [6 * 10 ** (-5), 10 ** (-3), 10 ** (-2)]

for tind in range(len(tfixlist)):
    tfix = tfixlist[tind]
    # Find v from lines

    SFS = np.zeros(n)

    for i in range(n_forward):
        SFS += n * np.loadtxt(
        'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
                 N, s, m, r, tfinal, n, tfix, nSFS, i))

    
    SFS /= n_forward


    plt.loglog(f_short, 
                 moving_average(SFS, navg, start_smooth), 
                 linewidth = 5, 
                 color = viridis_cmap(tind * 0.4))

    plt.vlines((tfix + L / v) / (N * L), 10 ** 3, 10 ** 11, linestyle = 'dotted',
               linewidth = 5, color = viridis_cmap(tind * 0.4))
    plt.text(xpositions[tind], 200, 
             r'$t = $' + '{:.0e}'.format(tfix), 
             color = viridis_cmap(tind * 0.4), fontsize = 50)    
    popt, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS[fit_range_ind])
    fitlist.append(popt[0])


plt.text(10 ** (-4), 5 * 10 ** 11, 
         r'$f = t / N + 1 / (\rho v), t = 10^3, 10^4, 10^5$', color = 'k')
plt.loglog(f_short, 
           Un / s / f_short ** 2, 
           label = r'$U / (s f^2)$', 
           linestyle = '-.', linewidth = 5, color = '#0072b2')

plt.loglog(f_short, np.ones(len(f_short)) * L / v, linewidth = 5, 
           linestyle = 'dotted'
              , label = r'$ U L / v$', color = '#d55e00')


plt.loglog(f_short, 
           2 * Un * N * L * np.ones(len(f_short)) / f_short,
           label = r'$2 N U / f$', linestyle = '-.',
           linewidth = 5, color = '#cc79a7')
xvals = [5 * 10 ** (-5), 5 * 10 ** (-5), 0.2]
labelLines(plt.gca().get_lines() ,xvals = xvals, fontsize = 75)



#plt.title('L = {}, '.format(L)  +
#          r'$\rho =$' + '{}, s = {:.2f}, m = {:.2f}, r = {:.2f}, sample uniformly everywhere, 100 forward sims'.format(N, s, m, r))

#plt.legend(loc = 'upper right')


#plt.figure(figsize = (8, 6))
#plt.loglog(slist, fitlist, 'o', label = 'simulation')
#plt.loglog(slist, 2 / np.array(slist), label = '$k = 2 U_n / s$')
#plt.xlabel('s')
#plt.ylabel('k')
#plt.legend()
#plt.title('fitting SFS to $f = k / f^2$')
