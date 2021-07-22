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

from matplotlib import cm

from labellines import labelLines
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 60, weight = 'bold')
plt.figure(figsize = (24, 18))
plt.xlabel(r'Frequency, $f$', fontsize = 75)
plt.ylabel(r'Number of alleles, $P(f)$', fontsize = 75)

N = 10000000
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

tfinallist = [1000000, 100000, 100000, 100000]
rholist = [20000, 10000, 5000, 2000]
viridis_cmap = cm.get_cmap('viridis')

colorlist = ['y', 'r', 'c', 'g', 'm', 'b', 'k']
fitlist = []

f = np.arange(1, n + 1) / n
navg = 30
start_smooth = 100
f_short = moving_average(f, navg, start_smooth)

fit_range_ind = np.arange(200,1000)


xpositions = [6 * 10 ** (-5), 10 ** (-4), 10 ** (-3), 10 ** (-2)]

for rhoind in range(len(rholist)):
    rho = rholist[rhoind]
    L = int(N / rho)
    tfinal = tfinallist[rhoind]
    # Find v from lines

    SFS = np.zeros(n)

    for i in range(n_forward):
        SFS += n * np.loadtxt(
        'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
                 rho, s, m, r, tfinal, n, tfix, nSFS, i))

    
    SFS /= n_forward

    freq_file = open(
       'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, rho, s, m, tfinal))
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

    plt.loglog(f_short, 
                 moving_average(SFS, navg, start_smooth), 
                 linewidth = 5, label = r'$\rho = $ {:.0e}, $L = $ {:.0e}'.format(rho, L), 
                 color = viridis_cmap(rhoind * 0.4))

    plt.vlines((tfix + L / v) / N, 10 ** 3, 10 ** 11, linestyle = 'dotted',
               linewidth = 5, color = viridis_cmap(rhoind * 0.4))
#    plt.text(xpositions[rhoind], 200, 
#             r'$\rho = $ {:.0e}, $L = $ {:.0e}'.format(rho, L), 
#             color = viridis_cmap(rhoind * 0.4), fontsize = 50)    

    plt.loglog(f_short, np.ones(len(f_short)) * L / v, linewidth = 5, 
               linestyle = 'dashed', 
                  color = viridis_cmap(rhoind * 0.4))


plt.legend(loc = 'upper right', fontsize = 40)

plt.text(5 * 10 ** (-4), 5 * 10 ** 11, 
         r'$f = 1 / (\rho v)$', color = 'k')
plt.text(0.15, 1.5 * 10 ** 5, 
         r'$ U L / v$', color = 'k')

plt.loglog(f_short, 
           Un / s / f_short ** 2, 
           label = r'$U / (s f^2)$', 
           linestyle = '-.', linewidth = 5, color = '#0072b2')



plt.loglog(f_short, 
           2 * Un * N * np.ones(len(f_short)) / f_short,
           label = r'$2 N U / f$', linestyle = '-.',
           linewidth = 5, color = '#cc79a7')
xvals = [5 * 10 ** (-5), 0.2]

lines = plt.gca().get_lines()

labelLines(lines[-2:] ,xvals = xvals, fontsize = 75)






