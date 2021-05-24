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
plt.rc('font', family='serif', size = 70)
plt.figure(figsize = (24, 18))
plt.xlabel(r'Frequency, $f$')
plt.ylabel(r'Number of alleles, $P(f)$')

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


#plt.loglog(f_short, 
#             moving_average(SFS, navg, start_smooth), 
#             linewidth = 6)
#
#neutral, = plt.loglog(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           label = r'$2 N U_n / f$', linestyle = '--', linewidth = 6,
#           alpha = 0.7)
#
#
#bsc, = plt.loglog(f_short, 
#           Un / s / f_short ** 2, label = r'$U_n / (s f^2)$', 
#           linestyle = '-.', linewidth = 6, alpha = 0.7)
#
#
#flat, = plt.loglog(f_short, np.ones(len(f_short)) * L / v, linewidth = 6, 
#           linestyle = 'dotted' , alpha = 0.7
#              , label = r'$U_n L / v$')
#
#xvals = [0.002, 0.0001, 0.018]
#labelLines(plt.gca().get_lines(),xvals=xvals)
#
#
#plt.axvline(1 / (N * v), linestyle = 'dotted', linewidth = 5)
#vline = plt.axvline(1 / (N * v), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 5, label = r'$f = 1 / \rho v$')
#
#plt.text(0.0004, 10 ** 13, r'$f = 1 / \rho v$', color = vline.get_color())
##vline_legend = plt.legend(handles = [vline], loc = (0.02, 0.03), 
##                          fontsize = 'small')
##plt.gca().add_artist(vline_legend)
#
#plt.show()
#plt.close()

##########################3
plt.figure(figsize = (24, 18))
#plt.xlabel(r'Frequency, $f$')
#plt.ylabel(r'Number of alleles, $P(f)$')

plt.semilogy(f_short, 
             moving_average(SFS, navg, start_smooth), 
             linewidth = 6)
flat, = plt.semilogy(f_short, np.ones(len(f_short)) * L / v, linewidth = 10, 
           linestyle = 'dotted' , alpha = 0.7
              , label = r'$U_n L / v$', color = 'r')
plt.show()


