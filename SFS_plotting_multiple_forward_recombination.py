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
N = 20000
s = 0.05
m = 0.25
tfinal = 1000000
Un = 1
r = 0
n_forward = 100
tfix = 0


# Change r - sample everywhere
n = 100000
nSFS = 1000
nSFS2 = 2000
nSFS3 = 10000

r0 = 0
r1 = 0.000001
r2 = 0.00001
r3 = 0.00003
r4 = 0.0001
r5 = 0.0005
r6 = 0.005

f = np.arange(1, n + 1) / n
navg = 30
navg2 = 200
start_smooth = 50
f_short = moving_average(f, navg, start_smooth)
f_short2 = moving_average(f, navg2, start_smooth)
# Find v from lines

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

v0 = np.average(v_list)





SFS0 = np.zeros(n)
SFS1 = np.zeros(n)
SFS2 = np.zeros(n)
SFS3 = np.zeros(n)
SFS4 = np.zeros(n)
SFS5 = np.zeros(n)
SFS6 = np.zeros(n)

for i in range(n_forward):
    SFS0 += n * np.loadtxt(
    'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
             N, s, m, r0, tfinal, n, tfix, nSFS, i))
    SFS1 += n * np.loadtxt(
    'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
             N, s, m, r1, tfinal, n, tfix, nSFS3, i))
    SFS2 += n * np.loadtxt(
    'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
             N, s, m, r2, tfinal, n, tfix, nSFS2, i))
    SFS3 += n * np.loadtxt(
    'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
             N, s, m, r3, tfinal, n, tfix, nSFS2, i))
    SFS4 += n * np.loadtxt(
    'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
             N, s, m, r4, tfinal, n, tfix, nSFS2, i))
    SFS5 += n * np.loadtxt(
    'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
             N, s, m, r5, tfinal, n, tfix, nSFS2, i))
    SFS6 += n * np.loadtxt(
    'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
             N, s, m, r6, tfinal, n, tfix, nSFS2, i))

    
SFS0 /= n_forward
SFS1 /= n_forward
SFS2 /= n_forward
SFS3 /= n_forward
SFS4 /= n_forward
SFS5 /= n_forward
SFS6 /= n_forward




plasma_cmap = cm.get_cmap('plasma')


plt.figure(figsize = (24, 18))
plt.xlabel(r'Frequency, $f$', fontsize = 75)
plt.ylabel(r'Number of alleles, $P(f)$', fontsize = 75)

plt.loglog(f_short, 
             moving_average(SFS0, navg, start_smooth), 
             label = '$r_0 =$ {:.6f}'.format(r0), linewidth = 2, color = 
             plasma_cmap(0))
plt.loglog(f_short2, 
             moving_average(SFS1, navg2, start_smooth), 
             label = '$r_1 =$ {:.6f}'.format(r1), linewidth = 2, 
             color = plasma_cmap(0.15))
plt.loglog(f_short2, 
             moving_average(SFS2, navg2, start_smooth), 
             label = '$r_2 =$ {:.6f}'.format(r2), linewidth = 2,
             color = plasma_cmap(0.3))
plt.loglog(f_short2, 
             moving_average(SFS3, navg2, start_smooth), 
             label = '$r_3 =$ {:.6f}'.format(r3), linewidth = 2, 
             color = plasma_cmap(0.45))
plt.loglog(f_short2, 
             moving_average(SFS4, navg2, start_smooth), 
             label = '$r_4 =$ {:.6f}'.format(r4), linewidth = 2, 
             color = plasma_cmap(0.6))
plt.loglog(f_short2, 
             moving_average(SFS5, navg2, start_smooth), 
             label = '$r_5 =$ {:.6f}'.format(r5), linewidth = 2, 
             color = plasma_cmap(0.75))

plt.loglog(f_short2, 
             moving_average(SFS6, navg2, start_smooth), 
             label = '$r_6 =$ {:.6f}'.format(r6), linewidth = 2, 
             color = plasma_cmap(0.9))

plt.loglog(f_short, 
           2 * Un * N * L * np.ones(len(f_short)) / f_short,
           label = r'$p(f) = 2 N U / f$', linestyle = '--', linewidth = 6, 
           color = '#cc79a7')
plt.loglog(f_short, 
           Un / s / f_short ** 2, label = r'$p(f) = U / (s f^2)$', 
           linestyle = '-.', color = '#0072b2', linewidth = 6)
plt.loglog(f_short, np.ones(len(f_short)) * L / v0, linewidth = 6, linestyle = '--'
              , label = r'$p(f) = U L / v$', color = '#d55e00')

plt.vlines(1 / (N * v0), 10 ** 3, 10 ** 11, linestyle = 'dotted',
           linewidth = 6, color = '#009e73', label = r'$f = 1 / \rho v$')


#plt.title('L = {}, '.format(L) + 
#          r'$\rho =$ {}, m = {:.2f}, s = {:.2f}, sample uniformly everywhere, 100 forward sims'.format(N, m, s))

plt.legend(fontsize = 'medium', loc = 'lower left')
