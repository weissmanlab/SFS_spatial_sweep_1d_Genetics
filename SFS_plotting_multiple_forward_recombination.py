import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

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

from labellines import labelLines

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 30, weight = 'bold')

#L = 500
#N = 20000
#s = 0.05
#m = 0.25
#tfinal = 1000000
#Un = 1
#r = 0
#n_forward = 100

L = 500
rho = 20000
N = rho * L
s = 0.05
m = 0.25
tfinal = 1000000
Un = 1
r = 0
n_forward = 100
tfix = 0


# Change r - sample everywhere
n = 100000


rlist = [0.005, 0.0005, 0.0001, 0.00003, 0.00001, 0.000001, 0]

f = np.arange(1, n + 1) / n
navg = 200
start_smooth = 400
f_short = moving_average(f, navg, start_smooth)
# Find v from lines

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

v0 = np.average(v_list)

plasma_cmap = cm.get_cmap('plasma')


plt.figure(figsize = (24, 18))
plt.xlabel(r'Frequency, $f$', fontsize = 75)
plt.ylabel(r'Number of alleles, $P(f)$', fontsize = 75)

plt.loglog(f_short, 
           2 * Un * N * np.ones(len(f_short)) / f_short,
           label = r'$p(f) = 2 N U / f$', linestyle = '--', linewidth = 6, 
           color = '#cc79a7')
plt.loglog(f_short, np.ones(len(f_short)) * L / v0, linewidth = 6, linestyle = '--'
              , label = r'$p(f) = U L / v$', color = '#d55e00')

f_short2 = np.linspace(1 / (rho * v0), 1, 100)
plt.vlines(1 / (rho * v0), 10 ** 3, 10 ** 11, linestyle = 'dotted',
           linewidth = 6, color = '#009e73', label = r'$f = 1 / \rho v$')


for rind in range(len(rlist)):
    r = rlist[rind]
    Uneff = Un * (1 + 2 * N * r) 
    if rind == 6:
        nSFS = 1000
    elif rind == 5:
        nSFS = 10000
    else:
        nSFS = 2000

    SFS = np.zeros(n)

    for i in range(n_forward):
            SFS += n * np.loadtxt(
    'backward simulation data/expected_SFS_L=' 
    + '{}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
             rho, s, m, r, tfinal, n, tfix, nSFS, i))
    
    SFS /= n_forward
    plt.loglog(f_short, 
             moving_average(SFS, navg, start_smooth), 
             label = '$r =$ {:.6f}'.format(r), linewidth = 2, color = 
             plasma_cmap(rind / len(rlist)), alpha = 0.8)
    if r < 5 * 10 ** (-4):
        plt.loglog(f_short2, 
           2.5 * Uneff / s / f_short2 ** 2, 
           linestyle = '-.', color = plasma_cmap(rind / len(rlist)), 
           linewidth = 6, alpha = 0.8)



plt.legend(fontsize = 'medium', loc = 'lower left')

