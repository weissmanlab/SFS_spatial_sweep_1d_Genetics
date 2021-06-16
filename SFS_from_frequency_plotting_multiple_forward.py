'''
Plot AFS with one set of parameter values.
Used this to make a graphical summary in Adobe Illustrator.
'''



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


from labellines import labelLines
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 170, weight = 'bold')
plt.figure(figsize = (48, 36))
plt.xlabel(r'Frequency, $f$', fontsize = 180)
plt.ylabel(r'Number of alleles, $P(f)$', fontsize = 180)

L = 500
N = 20000
s = 0.05
m = 0.25
tfinal = 1000000
Un = 1
r = 0
n_forward = 100
n = 100000
nSFS = 1000
tfix = 0

plt.xlim((1/n, 1))
plt.ylim((2, 5 * 10 ** 14))
f = np.arange(1, n + 1) / n
navg = 10
start_smooth = 100
f_short = moving_average(f, navg, start_smooth)

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

v = np.average(v_list)

SFS = np.zeros(n)
for i in range(n_forward):
    SFS += n * np.loadtxt(
   'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
             N, s, m, r, tfinal, n, tfix, nSFS, i))

SFS /= n_forward

#####################################
#plt.loglog(f_short, 
#             moving_average(SFS, navg, start_smooth), 
#             linewidth = 12, color = 'k')
#
##plt.vlines(1 / (N * v), 10 ** (-2), 10 ** 15, linestyle = 'dotted',
##           linewidth = 10, color = '#009e73')
#
##vline = plt.axvline(1 / (N * v), 10 ** 3, 10 ** 11, linestyle = 'dotted',
##           linewidth = 10)
#
##plt.text(0.0004, 50, r'$\boldmath{f = 1 / \rho v}$', 
##         color = '#009e73', fontsize = 150)
#
#
#neutral, = plt.loglog(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           label = r'$\boldmath{2 N U / f}$', linestyle = '--', linewidth = 12,
#           color = '#cc79a7')
#
#
#bsc, = plt.loglog(f_short, 
#           Un / s / f_short ** 2, label = r'$\boldmath{U / (s f^2)}$', 
#           linestyle = '-.', linewidth = 12, color = '#0072b2')
#
#
#flat, = plt.loglog(f_short, np.ones(len(f_short)) * L / v, linewidth = 12, 
#           linestyle = 'dotted'
#              , label = r'$\boldmath{U L / v}$', color = '#d55e00')
#
#xvals = [0.002, 0.0002, 0.01]
#labelLines(plt.gca().get_lines() ,xvals = xvals, fontsize = 150)
#
#
#
#
#
#
##vline_legend = plt.legend(handles = [vline], loc = (0.02, 0.03), 
##                          fontsize = 'small')
##plt.gca().add_artist(vline_legend)
#
#plt.show()

##########################3
plt.rc('font', family='serif', size = 150, weight = 'bold')
plt.figure(figsize = (24, 12))
#plt.xlabel(r'Frequency, $f$')
#plt.ylabel(r'Number of alleles, $P(f)$')

plt.semilogy(f_short, 
             moving_average(SFS, navg, start_smooth), 
             linewidth = 10, color = 'k')
flat, = plt.semilogy(f_short, np.ones(len(f_short)) * L / v, linewidth = 30, 
           linestyle = 'dotted' 
              , label = r'$U_n L / v$', color = '#d55e00')
plt.locator_params(axis='x', nbins = 3)
plt.ylim((10, 5 * 10 ** 6))

plt.show()


