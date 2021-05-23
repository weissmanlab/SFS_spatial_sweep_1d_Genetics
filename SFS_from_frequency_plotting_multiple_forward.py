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

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=50)
#plt.rcParams.update({'font.size': 50})
#axfont = {'fontname' : 'Helvetica'}
plt.figure(figsize = (24, 18))
plt.xlabel(r'Frequency, $f$')
plt.ylabel(r'Number of alleles, $P(f)$')

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
## 0) 
n = 100000
nSFS = 1000
tfix = 0

plt.xlim((1/n, 1))

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


plt.loglog(f_short, 
             moving_average(SFS, navg, start_smooth), 
             linewidth = 6)

neutral, = plt.loglog(f_short, 
           2 * Un * N * L * np.ones(len(f_short)) / f_short,
           label = r'$P(f) = 2 N U_n / f$', linestyle = '--', linewidth = 6,
           alpha = 0.7)

neutral_legend = plt.legend(handles = [neutral], loc = (0.61, 0.72), 
                            fontsize = 'small')
plt.gca().add_artist(neutral_legend)

bsc, = plt.loglog(f_short, 
           Un / s / f_short ** 2, label = r'$P(f) = U_n / (s f^2)$', 
           linestyle = '-.', linewidth = 6, alpha = 0.7)

bsc_legend = plt.legend(handles = [bsc], loc = (0.54, 0.48), 
                        fontsize = 'small')
plt.gca().add_artist(bsc_legend)

flat, = plt.loglog(f_short, np.ones(len(f_short)) * L / v, linewidth = 6, 
           linestyle = 'dotted' , alpha = 0.7
              , label = r'$P(f) = U_n L / v$')
flat_legend = plt.legend(handles = [flat], loc = (0.31, 0.25), 
                         fontsize = 'small')
plt.gca().add_artist(flat_legend)

plt.axvline(1 / (N * v), linestyle = 'dotted', linewidth = 5)
vline = plt.axvline(1 / (N * v), 10 ** 3, 10 ** 11, linestyle = 'dotted',
           linewidth = 5, label = r'$f = 1 / \rho v$')
vline_legend = plt.legend(handles = [vline], loc = (0.02, 0.03), 
                          fontsize = 'small')
plt.gca().add_artist(vline_legend)



#plt.semilogy(f_short, 
#             moving_average(SFS, navg, start_smooth), 
#             linewidth = 0.9)
#plt.semilogy(f_short, np.ones(len(f_short)) * L / v, linewidth = 2, 
#           linestyle = '--', alpha = 0.5
#              , label = '$p(f) = U_n L / v$', color = 'r')


plt.show()
#plt.title(
#  'L = {}, '.format(L) + r'$\rho$' 
#  + ' = {}, r = {:.2f}, m = {:.2f}, s = {:.2f}, # forward sims = {}, # backward sims = {}'.format(N, r, m, s, n_forward, nSFS))
##########################################################################
# 1)  Change r
#n = 1000
#nSFS = 5000
#r2 = 0.1
#r3 = 0.001
#r4 = 0.0001
#r5 = 0.00001
#f = np.arange(1, n + 1) / n
#navg = 5
#start_smooth = 900
#f_short = moving_average(f, navg, start_smooth)
#
## Find v from lines
#
#freq_file = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s, m, tfinal))
#lines = np.loadtxt(freq_file, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v = np.average(v_list)
#
#SFS = n * np.loadtxt(
#  'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#            N, s, m, r, tfinal, n, Un, nSFS))
#SFS2 = n * np.loadtxt(
#  'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#            N, s, m, r2, tfinal, n, Un, nSFS))
#SFS3 = n * np.loadtxt(
#  'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#            N, s, m, r3, tfinal, n, Un, nSFS))
#SFS4 = n * np.loadtxt(
#  'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#            N, s, m, r4, tfinal, n, Un, nSFS))
#SFS5 = n * np.loadtxt(
#  'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#            N, s, m, r5, tfinal, n, Un, nSFS))
#
#plt.loglog(f_short, 
#             moving_average(SFS, navg, start_smooth), 
#             label = 'r = {:.5f}'.format(r), linewidth = 0.9)
#plt.loglog(f_short, 
#             moving_average(SFS2, navg, start_smooth), 
#             label = 'r = {:.5f}'.format(r2), linewidth = 0.9)
#plt.loglog(f_short, 
#             moving_average(SFS3, navg, start_smooth), 
#             label = 'r = {:.5f}'.format(r3), linewidth = 0.9)
#plt.loglog(f_short, 
#             moving_average(SFS4, navg, start_smooth), 
#             label = 'r = {:.5f}'.format(r4), linewidth = 0.9)
#plt.loglog(f_short, 
#             moving_average(SFS5, navg, start_smooth), 
#             label = 'r = {:.5f}'.format(r5), linewidth = 0.9)
#
#plt.loglog(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2)
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = 0.5 U_n L / v$')
#plt.legend()
#plt.title(
#        'L = {}, N = {}, s = {:.2f}, m = {:.2f}, sample middle ([0.25 L, 0.75 L])'.format(L, N, s, m))
##############################################################################
### 2) Change s
#n = 100000
#nSFS = 5000
#s2 = 0.03
#s3 = 0.02
#s4 = 0.01
#
#f = np.arange(1, n + 1) / n
#navg = 200
#start_smooth = 300
#f_short = moving_average(f, navg, start_smooth)
#
## Find v from lines
#
#freq_file = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s, m, tfinal))
#lines = np.loadtxt(freq_file, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v = np.average(v_list)
#
#freq_file2 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s2, m, tfinal))
#lines = np.loadtxt(freq_file2, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v2 = np.average(v_list)
#
#
#freq_file3 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s3, m, tfinal))
#lines = np.loadtxt(freq_file3, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v3 = np.average(v_list)
#
#
#freq_file4 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s4, m, tfinal))
#lines = np.loadtxt(freq_file4, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v4 = np.average(v_list)
#
#
#SFS = n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m, r, tfinal, n, Un, nSFS))
#
#SFS2 = n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N, s2, m, r, tfinal, n, Un, nSFS))
#SFS3 = n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N, s3, m, r, tfinal, n, Un, nSFS))
#SFS4 = n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N, s4, m, r, tfinal, n, Un, nSFS))
#
#plt.loglog(f_short, 
#             moving_average(SFS, navg, start_smooth), 
#             label = '$s_1 =$ {:.2f}'.format(s), linewidth = 0.9,
#             color = 'b')
#plt.loglog(f_short, 
#             moving_average(SFS2, navg, start_smooth), 
#             label = '$s_2 =$ {:.2f}'.format(s2), linewidth = 0.9,
#             color = 'k')
#plt.loglog(f_short, 
#             moving_average(SFS3, navg, start_smooth), 
#             label = '$s_3 =$ {:.2f}'.format(s3), linewidth = 0.9,
#             color = 'g')
#plt.loglog(f_short, 
#             moving_average(SFS4, navg, start_smooth), 
#             label = '$s_4 =$ {:.2f}'.format(s4), linewidth = 0.9,
#             color = 'r')
#plt.loglog(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2,
#           alpha = 0.5)
#plt.loglog(f_short, 
#           Un / s / f_short ** 2, label = '$p(f) = U_n / (s_1 f^2)$', 
#           linestyle = '--', linewidth = 2, color = 'b', alpha = 0.5)
#plt.loglog(f_short, 
#           Un / s2 / f_short ** 2, label = '$p(f) = U_n / (s_2 f^2)$', 
#           linestyle = '--', linewidth = 2, color = 'k', alpha = 0.5)
#plt.loglog(f_short, 
#           Un / s3 / f_short ** 2, label = '$p(f) = U_n / (s_3 f^2)$', 
#           linestyle = '--', linewidth = 2, color = 'g', alpha = 0.5)
#plt.loglog(f_short, 
#           Un / s4 / f_short ** 2, label = '$p(f) = U_n / (s_4 f^2)$', 
#           linestyle = '--', linewidth = 2, color = 'r', alpha = 0.5)
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v, linewidth = 2, linestyle = '-.'
#              , label = '$p(f) = 0.5 U_n L / v$', color = 'b', alpha = 0.5)
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v2, linewidth = 2, linestyle = '-.'
#              , label = '$p(f) = 0.5 U_n L / v_2$', color = 'k', alpha = 0.5)
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v3, linewidth = 2, linestyle = '-.'
#              , label = '$p(f) = 0.5 U_n L / v_3$', color = 'g', alpha = 0.5)
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v4, linewidth = 2, linestyle = '-.'
#              , label = '$p(f) = 0.5 U_n L / v_4$', color = 'r', alpha = 0.5)
#
##plt.semilogy(f_short, 
##             moving_average(SFS, navg, start_smooth), 
##             label = '$s_1 =$ {:.2f}'.format(s), linewidth = 0.9,
##             color = 'b')
##plt.semilogy(f_short, 
##             moving_average(SFS2, navg, start_smooth), 
##             label = '$s_2 =$ {:.2f}'.format(s2), linewidth = 0.9,
##             color = 'k')
##plt.semilogy(f_short, 
##             moving_average(SFS3, navg, start_smooth), 
##             label = '$s_3 =$ {:.2f}'.format(s3), linewidth = 0.9,
##             color = 'g')
##plt.semilogy(f_short, 
##             moving_average(SFS4, navg, start_smooth), 
##             label = '$s_4 =$ {:.2f}'.format(s4), linewidth = 0.9,
##             color = 'r')
##plt.semilogy(f_short, 
##           2 * Un * N * L * np.ones(len(f_short)) / f_short,
##           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2,
##           alpha = 0.5)
##plt.semilogy(f_short, 
##           Un / s / f_short ** 2, label = '$p(f) = U / (s_1 f^2)$', 
##           linestyle = '--', linewidth = 2, color = 'b')
##plt.semilogy(f_short, 
##           Un / s2 / f_short ** 2, label = '$p(f) = U / (s_2 f^2)$', 
##           linestyle = '--', linewidth = 2, color = 'k')
##plt.semilogy(f_short, 
##           Un / s3 / f_short ** 2, label = '$p(f) = U / (s_3 f^2)$', 
##           linestyle = '--', linewidth = 2, color = 'g')
##plt.semilogy(f_short, 
##           Un / s4 / f_short ** 2, label = '$p(f) = U / (s_4 f^2)$', 
##           linestyle = '--', linewidth = 2, color = 'r')
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v, linewidth = 2, linestyle = '-.'
##              , label = '$p(f) = 0.5 U_n L / v$', color = 'b')
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v2, linewidth = 2, linestyle = '-.'
##              , label = '$p(f) = 0.5 U_n L / v_2$', color = 'k')
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v3, linewidth = 2, linestyle = '-.'
##              , label = '$p(f) = 0.5 U_n L / v_3$', color = 'g')
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v4, linewidth = 2, linestyle = '-.'
##              , label = '$p(f) = 0.5 U_n L / v_4$', color = 'r')
#
#plt.vlines(1 / (N * v), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'b', label = r'$f = 1 / \rho v_1$')
#plt.vlines(1 / (N * v2), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'k', label = r'$f = 1 / \rho v_2$')
#plt.vlines(1 / (N * v3), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'g', label = r'$f = 1 / \rho v_3$')
#plt.vlines(1 / (N * v4), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'r', label = r'$f = 1 / \rho v_4$')
#plt.title(
#  'L = {}, '.format(L) 
#  + r'$\rho =$' 
#  + ' {}, r = {:.2f}, m = {:.2f}, sample middle ([0.25 L, 0.75 L])'.format(N, r, m))
#

##############################################################################

### 2-2) Change s, sample uniformly from all demes
#n = 100000
#nSFS = 1000
#s2 = 0.04
#s3 = 0.03
#s4 = 0.02
#s5 = 0.01
#
#f = np.arange(1, n + 1) / n
#navg = 10
#start_smooth = 100
#f_short = moving_average(f, navg, start_smooth)
#
## Find v from lines
#
#freq_file = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s, m, tfinal))
#lines = np.loadtxt(freq_file, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v = np.average(v_list)
#
#freq_file2 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s2, m, tfinal))
#lines = np.loadtxt(freq_file2, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v2 = np.average(v_list)
#
#
#freq_file3 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s3, m, tfinal))
#lines = np.loadtxt(freq_file3, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v3 = np.average(v_list)
#
#
#freq_file4 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s4, m, tfinal))
#lines = np.loadtxt(freq_file4, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v4 = np.average(v_list)
#
#
#freq_file5 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s5, m, tfinal))
#lines = np.loadtxt(freq_file5, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v5 = np.average(v_list)
#
#SFS = np.zeros(n)
#SFS2 = np.zeros(n)
#SFS3 = np.zeros(n)
#SFS4 = np.zeros(n)
#SFS5 = np.zeros(n)
#
#for i in range(n_forward):
#    SFS += n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s, m, r, tfinal, n, Un, nSFS, i))
#
#    SFS2 += n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s2, m, r, tfinal, n, Un, nSFS, i))
#    SFS3 += n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s3, m, r, tfinal, n, Un, nSFS, i))
#    SFS4 += n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s4, m, r, tfinal, n, Un, nSFS, i))
#    SFS5 += n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s5, m, r, tfinal, n, Un, nSFS, i))
#SFS /= n_forward
#SFS2 /= n_forward
#SFS3 /= n_forward
#SFS4 /= n_forward
#SFS5 /= n_forward
#
#
##plt.loglog(f_short, 
##             moving_average(SFS, navg, start_smooth), 
##             label = '$s_1 =$ {:.2f}'.format(s), linewidth = 0.9,
##             color = 'b')
##plt.loglog(f_short, 
##             moving_average(SFS2, navg, start_smooth), 
##             label = '$s_2 =$ {:.2f}'.format(s2), linewidth = 0.9,
##             color = 'k')
##plt.loglog(f_short, 
##             moving_average(SFS3, navg, start_smooth), 
##             label = '$s_3 =$ {:.2f}'.format(s3), linewidth = 0.9,
##             color = 'g')
##plt.loglog(f_short, 
##             moving_average(SFS4, navg, start_smooth), 
##             label = '$s_4 =$ {:.2f}'.format(s4), linewidth = 0.9,
##             color = 'r')
##plt.loglog(f_short, 
##             moving_average(SFS5, navg, start_smooth), 
##             label = '$s_5 =$ {:.2f}'.format(s5), linewidth = 0.9,
##             color = 'y')
##
##
##plt.loglog(f_short, 
##           2 * Un * N * L * np.ones(len(f_short)) / f_short,
##           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2,
##           alpha = 0.5)
##plt.loglog(f_short, 
##           Un / s / f_short ** 2, label = '$p(f) = U_n / (s_1 f^2)$', 
##           linestyle = '-.', linewidth = 2, color = 'b', alpha = 0.5)
##plt.loglog(f_short, 
##           Un / s2 / f_short ** 2, label = '$p(f) = U_n / (s_2 f^2)$', 
##           linestyle = '-.', linewidth = 2, color = 'k', alpha = 0.5)
##plt.loglog(f_short, 
##           Un / s3 / f_short ** 2, label = '$p(f) = U_n / (s_3 f^2)$', 
##           linestyle = '-.', linewidth = 2, color = 'g', alpha = 0.5)
##plt.loglog(f_short, 
##           Un / s4 / f_short ** 2, label = '$p(f) = U_n / (s_4 f^2)$', 
##           linestyle = '-.', linewidth = 2, color = 'r', alpha = 0.5)
##plt.loglog(f_short, 
##           Un / s5 / f_short ** 2, label = '$p(f) = U_n / (s_5 f^2)$', 
##           linestyle = '-.', linewidth = 2, color = 'y', alpha = 0.5)
##
##plt.loglog(f_short, np.ones(len(f_short)) * L / v, linewidth = 2, 
##           linestyle = '--', alpha = 0.5
##              , label = '$p(f) = U_n L / v$', color = 'b')
##plt.loglog(f_short, np.ones(len(f_short)) * L / v2, linewidth = 2, 
##           linestyle = '--', alpha = 0.5
##              , label = '$p(f) = U_n L / v_2$', color = 'k')
##plt.loglog(f_short, np.ones(len(f_short)) * L / v3, linewidth = 2, 
##           linestyle = '--', alpha = 0.5
##              , label = '$p(f) = U_n L / v_3$', color = 'g')
##plt.loglog(f_short, np.ones(len(f_short)) * L / v4, linewidth = 2,
##           linestyle = '--', alpha = 0.5
##              , label = '$p(f) = U_n L / v_4$', color = 'r')
##plt.loglog(f_short, np.ones(len(f_short)) * L / v5, linewidth = 2,
##           linestyle = '--', alpha = 0.5
##              , label = '$p(f) = U_n L / v_5$', color = 'y')
#
#plt.semilogy(f_short, 
#             moving_average(SFS, navg, start_smooth), 
#             label = '$s_1 =$ {:.2f}'.format(s), linewidth = 0.9,
#             color = 'b')
#plt.semilogy(f_short, 
#             moving_average(SFS2, navg, start_smooth), 
#             label = '$s_2 =$ {:.2f}'.format(s2), linewidth = 0.9,
#             color = 'k')
#plt.semilogy(f_short, 
#             moving_average(SFS3, navg, start_smooth), 
#             label = '$s_3 =$ {:.2f}'.format(s3), linewidth = 0.9,
#             color = 'g')
#plt.semilogy(f_short, 
#             moving_average(SFS4, navg, start_smooth), 
#             label = '$s_4 =$ {:.2f}'.format(s4), linewidth = 0.9,
#             color = 'r')
#plt.semilogy(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2,
#           alpha = 0.5)
#plt.semilogy(f_short, 
#           Un / s / f_short ** 2, label = '$p(f) = U / (s_1 f^2)$', 
#           linestyle = '--', linewidth = 2, color = 'b', alpha = 0.5)
#plt.semilogy(f_short, 
#           Un / s2 / f_short ** 2, label = '$p(f) = U / (s_2 f^2)$', 
#           linestyle = '--', linewidth = 2, color = 'k', alpha = 0.5)
#plt.semilogy(f_short, 
#           Un / s3 / f_short ** 2, label = '$p(f) = U / (s_3 f^2)$', 
#           linestyle = '--', linewidth = 2, color = 'g', alpha = 0.5)
#plt.semilogy(f_short, 
#           Un / s4 / f_short ** 2, label = '$p(f) = U / (s_4 f^2)$', 
#           linestyle = '--', linewidth = 2, color = 'r', alpha = 0.5)
#plt.semilogy(f_short, np.ones(len(f_short)) * L / v, linewidth = 2, linestyle = '-.'
#              , label = '$p(f) = U_n L / v$', color = 'b', alpha = 0.5)
#plt.semilogy(f_short, np.ones(len(f_short)) * L / v2, linewidth = 2, linestyle = '-.'
#              , label = '$p(f) = U_n L / v_2$', color = 'k', alpha = 0.5)
#plt.semilogy(f_short, np.ones(len(f_short)) * L / v3, linewidth = 2, linestyle = '-.'
#              , label = '$p(f) = U_n L / v_3$', color = 'g', alpha = 0.5)
#plt.semilogy(f_short, np.ones(len(f_short)) * L / v4, linewidth = 2, linestyle = '-.'
#              , label = '$p(f) = U_n L / v_4$', color = 'r', alpha = 0.5)
#
#plt.vlines(1 / (N * v), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'b', label = r'$f = 1 / \rho v_1$')
#plt.vlines(1 / (N * v2), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'k', label = r'$f = 1 / \rho v_2$')
#plt.vlines(1 / (N * v3), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'g', label = r'$f = 1 / \rho v_3$')
#plt.vlines(1 / (N * v4), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'r', label = r'$f = 1 / \rho v_4$')
#plt.vlines(1 / (N * v5), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'y', label = r'$f = 1 / \rho v_5$')
#
#plt.legend(fontsize = 'xx-small', loc = 'lower left')
#plt.title(
#  'L = {}, '.format(L) + r'$\rho$' 
#  + ' = {}, r = {:.2f}, m = {:.2f}, sample everywhere uniformly, 100 forward sims'.format(N, r, m))
#
#
#fit_range_ind = np.arange(200,1000)
#popt, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS[fit_range_ind])
#popt2, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS2[fit_range_ind])
#popt3, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS3[fit_range_ind])
#popt4, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS4[fit_range_ind])
#popt5, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS5[fit_range_ind])
#
#slst = np.array([s, s2, s3, s4, s5])
#plt.figure()
#plt.loglog(slst, [popt[0], popt2[0], popt3[0], popt4[0], popt5[0]], 'o', 
#           label = 'simulation')
#plt.loglog(slst, 1 / slst, label = '$k = U_n / s$')
#plt.xlabel('s')
#plt.ylabel('k')
#plt.title('fitting SFS to $f = k / f^2$')

#############################################################################
##3) Change m
#n = 10000
#nSFS = 1000
#
#m2 = 0.15
#m3 = 0.35
#m4 = 0.45
#
#f = np.arange(1, n + 1) / n
#navg = 30
#start_smooth = 40
#f_short = moving_average(f, navg, start_smooth)
#
## Find v from lines
#
#freq_file = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s, m, tfinal))
#lines = np.loadtxt(freq_file, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v = np.average(v_list)
#
#freq_file2 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s, m2, tfinal))
#lines = np.loadtxt(freq_file2, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v2 = np.average(v_list)
#
#
#freq_file3 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s, m3, tfinal))
#lines = np.loadtxt(freq_file3, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v3 = np.average(v_list)
#
#
#freq_file4 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s, m4, tfinal))
#lines = np.loadtxt(freq_file4, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v4 = np.average(v_list)
#
#
#SFS = n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m, r, tfinal, n, Un, nSFS))
#
#SFS2 = n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m2, r, tfinal, n, Un, nSFS))
#SFS3 = n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m3, r, tfinal, n, Un, nSFS))
#SFS4 = n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m4, r, tfinal, n, Un, nSFS))
#plt.loglog(f_short, 
#             moving_average(SFS, navg, start_smooth), 
#             label = '$m_1 =$ {:.2f}'.format(m), linewidth = 0.9, color = 'b')
#plt.loglog(f_short, 
#             moving_average(SFS2, navg, start_smooth), 
#             label = '$m_2 =$ {:.2f}'.format(m2), linewidth = 0.9, color = 'k')
#plt.loglog(f_short, 
#             moving_average(SFS3, navg, start_smooth), 
#             label = '$m_3 =$ {:.2f}'.format(m3), linewidth = 0.9, color = 'g')
#
#plt.loglog(f_short, 
#             moving_average(SFS4, navg, start_smooth), 
#             label = '$m_4 =$ {:.2f}'.format(m4), linewidth = 0.9, color = 'r')
#plt.loglog(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2)
#plt.loglog(f_short, 
#           Un / s / f_short ** 2, label = '$p(f) = U_n / (s f^2)$', 
#           linestyle = '--', linewidth = 2)
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = 0.5 U_n L / v$', color = 'b')
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v2, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = 0.5 U_n L / v_2$', color = 'k')
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v3, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = 0.5 U_n L / v_3$', color = 'g')
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v4, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = 0.5 U_n L / v_4$', color = 'r')
#
#
##plt.semilogy(f_short, 
##             moving_average(SFS, navg, start_smooth), 
##             label = '$m_1 =$ {:.2f}'.format(m), linewidth = 0.9, color = 'b')
##plt.semilogy(f_short, 
##             moving_average(SFS2, navg, start_smooth), 
##             label = '$m_2 =$ {:.2f}'.format(m2), linewidth = 0.9, color = 'k')
##plt.semilogy(f_short, 
##             moving_average(SFS3, navg, start_smooth), 
##             label = '$m_3 =$ {:.2f}'.format(m3), linewidth = 0.9, color = 'g')
##
##plt.semilogy(f_short, 
##             moving_average(SFS4, navg, start_smooth), 
##             label = '$m_4 =$ {:.2f}'.format(m4), linewidth = 0.9, color = 'r')
##plt.semilogy(f_short, 
##           2 * Un * N * L * np.ones(len(f_short)) / f_short,
##           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2)
##plt.semilogy(f_short, 
##           Un / s / f_short ** 2, label = '$p(f) = U_n / (s f^2)$', 
##           linestyle = '--', linewidth = 2)
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v, linewidth = 2, linestyle = '--'
##              , label = '$p(f) = 0.5 U_n L / v$', color = 'b')
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v2, linewidth = 2, linestyle = '--'
##              , label = '$p(f) = 0.5 U_n L / v_2$', color = 'k')
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v3, linewidth = 2, linestyle = '--'
##              , label = '$p(f) = 0.5 U_n L / v_3$', color = 'g')
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v4, linewidth = 2, linestyle = '--'
##              , label = '$p(f) = 0.5 U_n L / v_4$', color = 'r')
##
#
#
#plt.vlines(1 / (N * v), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'b', label = '$f = 1 / N v_1$')
#plt.vlines(1 / (N * v2), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'k', label = '$f = 1 / N v_2$')
#plt.vlines(1 / (N * v3), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'g', label = '$f = 1 / N v_3$')
#plt.vlines(1 / (N * v4), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'r', label = '$f = 1 / N v_4$')
#
#plt.title('L = {}, N = {}, r = {:.2f}, s = {:.2f}, sample middle ([0.25 L, 0.75 L])'.format(L, N, r, s))
############################################################################

#
##3-2) Change m - sample everywhere
#n = 100000
#nSFS = 1000
#
#m0 = 0.15
#m2 = 0.35
#m3 = 0.45
#
#f = np.arange(1, n + 1) / n
#navg = 10
#start_smooth = 100
#f_short = moving_average(f, navg, start_smooth)
#
## Find v from lines
#
#freq_file = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s, m0, tfinal))
#lines = np.loadtxt(freq_file, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v0 = np.average(v_list)
#
#freq_file2 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s, m, tfinal))
#lines = np.loadtxt(freq_file2, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v1 = np.average(v_list)
#
#
#freq_file3 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s, m2, tfinal))
#lines = np.loadtxt(freq_file3, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v2 = np.average(v_list)
#
#
#freq_file4 = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s, m3, tfinal))
#lines = np.loadtxt(freq_file4, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v3 = np.average(v_list)
#
#SFS0 = np.zeros(n)
#SFS1 = np.zeros(n)
#SFS2 = np.zeros(n)
#SFS3 = np.zeros(n)
#for i in range(n_forward):
#    SFS0 += n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s, m0, r, tfinal, n, Un, nSFS, i))
#
#    SFS1 += n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s, m, r, tfinal, n, Un, nSFS, i))
#    SFS2 += n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s, m2, r, tfinal, n, Un, nSFS, i))
#    SFS3 += n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s, m3, r, tfinal, n, Un, nSFS, i))
#    
#SFS0 /= n_forward
#SFS1 /= n_forward
#SFS2 /= n_forward
#SFS3 /= n_forward
#
#plt.loglog(f_short, 
#             moving_average(SFS0, navg, start_smooth), 
#             label = '$m_0 =$ {:.2f}'.format(m0), linewidth = 0.9, color = 'b')
#plt.loglog(f_short, 
#             moving_average(SFS1, navg, start_smooth), 
#             label = '$m_1 =$ {:.2f}'.format(m), linewidth = 0.9, color = 'k')
#plt.loglog(f_short, 
#             moving_average(SFS2, navg, start_smooth), 
#             label = '$m_2 =$ {:.2f}'.format(m2), linewidth = 0.9, color = 'g')
#
#plt.loglog(f_short, 
#             moving_average(SFS3, navg, start_smooth), 
#             label = '$m_3 =$ {:.2f}'.format(m3), linewidth = 0.9, color = 'r')
#plt.loglog(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2)
#plt.loglog(f_short, 
#           Un / s / f_short ** 2, label = '$p(f) = U_n / (s f^2)$', 
#           linestyle = '--', linewidth = 2)
#plt.loglog(f_short, np.ones(len(f_short)) * L / v0, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = U_n L / v_0$', color = 'b')
#plt.loglog(f_short, np.ones(len(f_short)) * L / v1, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = U_n L / v_1$', color = 'k')
#plt.loglog(f_short, np.ones(len(f_short)) * L / v2, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = U_n L / v_2$', color = 'g')
#plt.loglog(f_short, np.ones(len(f_short)) * L / v3, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = U_n L / v_3$', color = 'r')
#
#
##plt.semilogy(f_short, 
##             moving_average(SFS0, navg, start_smooth), 
##             label = '$m_0 =$ {:.2f}'.format(m0), linewidth = 0.9, color = 'b')
##plt.semilogy(f_short, 
##             moving_average(SFS1, navg, start_smooth), 
##             label = '$m_1 =$ {:.2f}'.format(m), linewidth = 0.9, color = 'k')
##plt.semilogy(f_short, 
##             moving_average(SFS2, navg, start_smooth), 
##             label = '$m_2 =$ {:.2f}'.format(m2), linewidth = 0.9, color = 'g')
##
##plt.semilogy(f_short, 
##             moving_average(SFS3, navg, start_smooth), 
##             label = '$m_3 =$ {:.2f}'.format(m3), linewidth = 0.9, color = 'r')
##plt.semilogy(f_short, 
##           2 * Un * N * L * np.ones(len(f_short)) / f_short,
##           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2)
##plt.semilogy(f_short, 
##           Un / s / f_short ** 2, label = '$p(f) = U_n / (s f^2)$', 
##           linestyle = '--', linewidth = 2)
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v0, linewidth = 2, linestyle = '--'
##              , label = '$p(f) = 0.5 U_n L / v_0$', color = 'b')
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v1, linewidth = 2, linestyle = '--'
##              , label = '$p(f) = 0.5 U_n L / v_1$', color = 'k')
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v2, linewidth = 2, linestyle = '--'
##              , label = '$p(f) = 0.5 U_n L / v_2$', color = 'g')
##plt.semilogy(f_short, np.ones(len(f_short)) * L * 0.50 / v3, linewidth = 2, linestyle = '--'
##              , label = '$p(f) = 0.5 U_n L / v_3$', color = 'r')
#
#
#plt.vlines(1 / (N * v0), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'b', label = '$f = 1 / N v_0$')
#plt.vlines(1 / (N * v1), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'k', label = '$f = 1 / N v_1$')
#plt.vlines(1 / (N * v2), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'g', label = '$f = 1 / N v_2$')
#plt.vlines(1 / (N * v3), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'r', label = '$f = 1 / N v_3$')
#
#
#plt.title('L = {}, N = {}, r = {:.2f}, s = {:.2f}, sample uniformly everywhere, 100 forward sims'.format(L, N, r, s))
#
#plt.legend(fontsize = 'small', loc = 'lower left')
#
#fit_range_ind = np.arange(200,1000)
#popt0, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS0[fit_range_ind])
#popt1, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS1[fit_range_ind])
#popt2, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS2[fit_range_ind])
#popt3, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS3[fit_range_ind])
#
#
#mlst = np.array([m0, m, m2, m3])
#plt.figure()
#plt.loglog(mlst, [popt0[0], popt1[0], popt2[0], popt3[0]], 'o', 
#           label = 'simulation')
#plt.xlabel('m')
#plt.ylabel('k')
#plt.title('fitting SFS to $f = k / f^2$')
#########################################################################
## 4) Change N, L
#n = 10000
#nSFS = 10000
#L2 = 5000
#N2 = 2000
#N3 = 200
#tfinal2 = 10000
#tfinal3 = 100000
#f = np.arange(1, n + 1) / n
#navg = 30
#start_smooth = 100
#f_short = moving_average(f, navg, start_smooth)
#
#
## Find v from lines
#
#freq_file = open('forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s, m, tfinal))
#lines = np.loadtxt(freq_file, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v = np.average(v_list)
#
#freq_file2 = open('forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N2, s, m, tfinal2))
#lines = np.loadtxt(freq_file2, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N2
#    psum2 = sum(line2) / N2
#    v_list.append(psum2 - psum1)
#
#v2 = np.average(v_list)
#
#
#freq_file3 = open('forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N3, s, m, tfinal2))
#lines = np.loadtxt(freq_file3, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N3
#    psum2 = sum(line2) / N3
#    v_list.append(psum2 - psum1)
#
#v3 = np.average(v_list)
#
#SFS = n * np.loadtxt(
#        'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m, r, tfinal, n, Un, nSFS))
#
#SFS2 = n * np.loadtxt(
#        'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N2, s, m, r, tfinal2, n, Un, nSFS))
#SFS3 = n * np.loadtxt(
#        'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N3, s, m, r, tfinal2, n, Un, nSFS))
#SFS4 = n * np.loadtxt(
#        'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L2, 
#             N2, s, m, r, tfinal3, n, Un, nSFS))
#SFS5 = n * np.loadtxt(
#        'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L2, 
#             N3, s, m, r, tfinal3, n, Un, nSFS))
#plt.loglog(f_short, 
#             moving_average(SFS, navg, start_smooth), 
#             label = 'L = {}, '.format(L) + r'$\rho$' 
#             + ' = {}'.format(N), linewidth = 2, color = 'r')
#plt.loglog(f_short, 
#             moving_average(SFS2, navg, start_smooth), 
#             label = 'L = {}, '.format(L) + r'$\rho$' 
#             + ' = {}'.format(N2), linewidth = 2,
#             linestyle = 'dashdot', color = 'r')
#plt.loglog(f_short, 
#             moving_average(SFS3, navg, start_smooth), 
#             label = 'L = {}, '.format(L) + r'$\rho$' 
#             + ' = {}'.format(N3), linewidth = 2,
#             linestyle = 'dotted', color = 'r')
#
#plt.loglog(f_short, 
#             moving_average(SFS4, navg, start_smooth), 
#             label = 'L = {}, '.format(L2) + r'$\rho$' 
#             + ' = {}'.format(N2), linewidth = 2,
#             linestyle = 'dashdot', color = 'b')
#plt.loglog(f_short, 
#             moving_average(SFS5, navg, start_smooth), 
#             label = 'L = {}, '.format(L2) + r'$\rho$' 
#             + ' = {}'.format(N3), linewidth = 2,
#             linestyle = 'dotted', color = 'b')
#
##plt.loglog(f_short[:2000], 
##           2 * Un * N * L * np.ones(len(f_short[:2000])) / f_short[:2000],
##           label = '$p(f) = 2 N U_n / f, N = 10^7$', 
##           linewidth = 1.5, color = 'r', alpha = 0.5)
##plt.loglog(f_short[:2000], 
##           2 * Un * N2 * L * np.ones(len(f_short[:2000])) / f_short[:2000],
##           label = '$p(f) = 2 N U_n / f, N = 10^6$', 
##           linewidth = 1.5, color = 'g', alpha = 0.5)
##
##plt.loglog(f_short[:2000], 
##           2 * Un * N3 * L * np.ones(len(f_short[:2000])) / f_short[:2000],
##           label = '$p(f) = 2 N U_n / f, N = 10^5$', 
##           linewidth = 1.5, color = 'b', alpha = 0.5)
#plt.loglog(f_short, 
#           Un / s / f_short ** 2 / (1 - np.pi ** 2 / 2 / np.log(N * np.sqrt(m * s / 2)) ** 2 ), label = r'$p(f) = U_n / (s f^2), \rho = ${}'.format(N), 
#           linestyle = '--', linewidth = 1.5, alpha = 0.5, color = 'k')
##plt.loglog(f_short, 
##           Un / s / f_short ** 2 * np.log(N2 * s), label = r'$p(f) = \ln(\rho s) U_n / (s f^2), \rho = ${}'.format(N2), 
##           linestyle = '--', linewidth = 1.5, alpha = 0.5, color = 'g')
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v, linewidth = 2, 
#           linestyle = '--', color = 'r', alpha = 0.5
#              , label = '$p(f) = 0.5 U_n L / v, L = 500$')
#plt.loglog(f_short, np.ones(len(f_short)) * L2 * 0.50 / v, linewidth = 2, 
#           linestyle = '--', color = 'b', alpha = 0.5
#              , label = '$p(f) = 0.5 U_n L / v, L = 5000$')
##plt.vlines(1 / (N * v), 10 ** 3, 10 ** 11, 
##           linestyle = 'solid', linewidth = 2, color = 'k', alpha = 0.5, 
##           label = 'cut off for' + r'$\rho = 200$')
##plt.vlines(1 / (N2 * v2), 10 ** 3, 10 ** 11, 
##           linestyle = 'dashdot', linewidth = 2, color = 'k', alpha = 0.5,
##           label = 'cut off for' + r'$\rho = 2000$')
##plt.vlines(1 / (N3 * v3), 10 ** 3, 10 ** 11, 
##           linestyle = 'dotted', linewidth = 2, color = 'k', alpha = 0.5,
##           label = 'cut off for' + r'$\rho = 20000$')
#plt.title(
#  's = {:.2f}, m = {:.2f}, r = {:.2f}, sample middle ([0.25 L, 0.75 L])'.format(s, m, r))


#############################################################################
##
##
## 5) Change tfix
#n = 10000
#nSFS = 1000
#tfix = 1000
#tfix2 = 10000
#tfix3 = 100000
#tfix4 = 1000000
#
#f = np.arange(1, n + 1) / n
#navg = 50
#start_smooth = 50
#f_short = moving_average(f, navg, start_smooth)
#
## Find v from lines
#
#freq_file = open(
#  'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s, m, tfinal))
#lines = np.loadtxt(freq_file, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v = np.average(v_list)
#
#
#SFS = n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m, r, tfinal, n, Un, nSFS))
#SFS2 = n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m, r, tfinal, n, tfix, nSFS))
#SFS3 = n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m, r, tfinal, n, tfix2, nSFS))
#SFS4 = n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m, r, tfinal, n, tfix3, nSFS))
#SFS5 = n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_middle_navg={}.txt'.format(L, 
#             N, s, m, r, tfinal, n, tfix4, nSFS))
#
#plt.loglog(f_short, 
#             moving_average(SFS, navg, start_smooth), 
#             label = 'sampled immediately after fixation',
#             linewidth = 0.9, color = 'b')
#plt.loglog(f_short, 
#             moving_average(SFS2, navg, start_smooth), 
#             label = 'sampled {} gens later'.format(tfix),
#             linewidth = 0.9, color = 'k')
#plt.loglog(f_short, 
#             moving_average(SFS3, navg, start_smooth), 
#             label = 'sampled {} gens later'.format(tfix2),
#             linewidth = 0.9, color = 'r')
#plt.loglog(f_short, 
#             moving_average(SFS4, navg, start_smooth), 
#             label = 'sampled {} gens later'.format(tfix3),
#             linewidth = 0.9, color = 'g')
#plt.loglog(f_short, 
#             moving_average(SFS5, navg, start_smooth), 
#             label = 'sampled {} gens later'.format(tfix4),
#             linewidth = 0.9, color = 'y')
#plt.loglog(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2)
#plt.loglog(f_short, 
#           Un / s / f_short ** 2, label = '$p(f) = U_n / (s f^2)$', 
#           linestyle = '--', linewidth = 2)
#plt.loglog(f_short, np.ones(len(f_short)) * L * 0.50 / v, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = 0.5 U_n L / v$')
#
#plt.vlines(tf_real / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'b')
#plt.vlines((tf_real + tfix) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'k')
#plt.vlines((tf_real + tfix2) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'r')
#plt.vlines((tf_real + tfix3) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'g')
#plt.vlines((tf_real + tfix4) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'y')
#plt.title('L = {}, '.format(L)
#        + r'$\rho$' + 
#        ' = {}, r = {:.2f}, m = {:.2f}, s = {:.2f}, sample middle ([0.25 L, 0.75 L])'.format(N, r, m, s))
#

##############################################################################
##
##
## 5-2) Change tfix sample uniform
#n = 100000
#nSFS = 1000
#tfix = 1000
#tfix2 = 10000
#tfix3 = 100000
##tfix4 = 1000000
#
#f = np.arange(1, n + 1) / n
#navg = 20
#start_smooth = 500
#f_short = moving_average(f, navg, start_smooth)
#
## Find v from lines
#
#freq_file = open(
#  'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s, m, tfinal))
#lines = np.loadtxt(freq_file, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v = np.average(v_list)
#
#SFS0 = np.zeros(n)
#SFS1 = np.zeros(n)
#SFS2 = np.zeros(n)
#SFS3 = np.zeros(n)
#for i in range(n_forward):
#
#    SFS0 += n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s, m, r, tfinal, n, Un, nSFS, i))
#    SFS1 += n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s, m, r, tfinal, n, tfix, nSFS, i))
#    SFS2 += n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s, m, r, tfinal, n, tfix2, nSFS, i))
#    SFS3 += n * np.loadtxt(
#   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s, m, r, tfinal, n, tfix3, nSFS, i))
##    SFS4 = n * np.loadtxt(
##   'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
##             N, s, m, r, tfinal, n, tfix4, nSFS, i))
#
#SFS0 /= n_forward
#SFS1 /= n_forward
#SFS2 /= n_forward
#SFS3 /= n_forward
##SFS4 /= n_forward
#
#
##plt.loglog(f_short, 
##             moving_average(SFS0, navg, start_smooth), 
##             label = 'sampled immediately after fixation',
##             linewidth = 0.9, color = 'b')
##plt.loglog(f_short, 
##             moving_average(SFS1, navg, start_smooth), 
##             label = 'sampled {} gens later'.format(tfix),
##             linewidth = 0.9, color = 'k')
##plt.loglog(f_short, 
##             moving_average(SFS2, navg, start_smooth), 
##             label = 'sampled {} gens later'.format(tfix2),
##             linewidth = 0.9, color = 'r')
##plt.loglog(f_short, 
##             moving_average(SFS3, navg, start_smooth), 
##             label = 'sampled {} gens later'.format(tfix3),
##             linewidth = 0.9, color = 'g')
###plt.loglog(f_short, 
###             moving_average(SFS4, navg, start_smooth), 
###             label = 'sampled {} gens later'.format(tfix4),
###             linewidth = 0.9, color = 'y')
##plt.loglog(f_short, 
##           2 * Un * N * L * np.ones(len(f_short)) / f_short,
##           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2)
##plt.loglog(f_short, 
##           Un / s / f_short ** 2, label = '$p(f) = U_n / (s f^2)$', 
##           linestyle = '--', linewidth = 2)
##plt.loglog(f_short, np.ones(len(f_short)) * L / v, linewidth = 2, linestyle = '--'
##              , label = '$p(f) = U_n L / v$')
##
##plt.vlines(tf_real / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
##           linewidth = 1, color = 'b')
##plt.vlines((tf_real + tfix) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
##           linewidth = 1, color = 'k')
##plt.vlines((tf_real + tfix2) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
##           linewidth = 1, color = 'r')
##plt.vlines((tf_real + tfix3) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
##           linewidth = 1, color = 'g')
###plt.vlines((tf_real + tfix4) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
###           linewidth = 1, color = 'y')
##plt.title('L = {}, '.format(L)
##        + r'$\rho$' + 
##        ' = {}, r = {:.2f}, m = {:.2f}, s = {:.2f}, sample everywhere'.format(N, r, m, s))
#
#plt.semilogy(f_short, 
#             moving_average(SFS0, navg, start_smooth), 
#             label = 'sampled immediately after fixation',
#             linewidth = 0.9, color = 'b')
#plt.semilogy(f_short, 
#             moving_average(SFS1, navg, start_smooth), 
#             label = 'sampled {} gens later'.format(tfix),
#             linewidth = 0.9, color = 'k')
#plt.semilogy(f_short, 
#             moving_average(SFS2, navg, start_smooth), 
#             label = 'sampled {} gens later'.format(tfix2),
#             linewidth = 0.9, color = 'r')
#plt.semilogy(f_short, 
#             moving_average(SFS3, navg, start_smooth), 
#             label = 'sampled {} gens later'.format(tfix3),
#             linewidth = 0.9, color = 'g')
##plt.semilogy(f_short, 
##             moving_average(SFS4, navg, start_smooth), 
##             label = 'sampled {} gens later'.format(tfix4),
##             linewidth = 0.9, color = 'y')
#plt.semilogy(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2)
#plt.semilogy(f_short, 
#           Un / s / f_short ** 2, label = '$p(f) = U_n / (s f^2)$', 
#           linestyle = '--', linewidth = 2)
#plt.semilogy(f_short, np.ones(len(f_short)) * L / v, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = U_n L / v$')
#
#plt.vlines(tf_real / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'b')
#plt.vlines((tf_real + tfix) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'k')
#plt.vlines((tf_real + tfix2) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'r')
#plt.vlines((tf_real + tfix3) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'g')
##plt.vlines((tf_real + tfix4) / N / L, 10**3, 10 ** 11, linestyle = 'dotted',
##           linewidth = 1, color = 'y')
#plt.title('L = {}, '.format(L)
#        + r'$\rho$' + 
#        ' = {}, r = {:.2f}, m = {:.2f}, s = {:.2f}, sample everywhere'.format(N, r, m, s))
############################################################################
#
#
##6) Change r - sample everywhere
#L = 200
#N = 8000
#tfinal = 10000
#n = 10000
#nSFS = 100000
#
#r2 = 0.00001
#r3 = 0.00002
#r4 = 0.00004
#r5 = 0.0001
#
#f = np.arange(1, n + 1) / n
#navg = 100
#start_smooth = 35
#f_short = moving_average(f, navg, start_smooth)
#
## Find v from lines
#
#freq_file = open(
#   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s, m, tfinal))
#lines = np.loadtxt(freq_file, dtype=np.int64)
#tf_real = len(lines)
#v_list = []
#for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#    line1 = lines[t]
#    line2 = lines[t + 1]
#    psum1 = sum(line1) / N
#    psum2 = sum(line2) / N
#    v_list.append(psum2 - psum1)
#
#v = np.average(v_list)
#
#
#SFS = n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_no_mixing.txt'.format(L, 
#             N, s, m, r, tfinal, n, Un, nSFS))
#
#SFS2 = n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_no_mixing.txt'.format(L, 
#             N, s, m, r2, tfinal, n, Un, nSFS))
#SFS3 = n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_no_mixing.txt'.format(L, 
#             N, s, m, r3, tfinal, n, Un, nSFS))
#SFS4 = n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_no_mixing.txt'.format(L, 
#             N, s, m, r4, tfinal, n, Un, nSFS))
#SFS5 = n * np.loadtxt(
#    'backward simulation data/SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_uniform_navg={}_no_mixing.txt'.format(L, 
#             N, s, m, r5, tfinal, n, Un, nSFS))
#
#plt.loglog(f_short, 
#             moving_average(SFS, navg, start_smooth), 
#             label = '$r_1 =$ {:.5f}'.format(r), linewidth = 0.9, color = 'b')
#plt.loglog(f_short, 
#             moving_average(SFS2, navg, start_smooth), 
#             label = '$r_2 =$ {:.5f}'.format(r2), linewidth = 0.9, color = 'k')
#plt.loglog(f_short, 
#             moving_average(SFS3, navg, start_smooth), 
#             label = '$r_3 =$ {:.5f}'.format(r3), linewidth = 0.9, color = 'g')
#
#plt.loglog(f_short, 
#             moving_average(SFS4, navg, start_smooth), 
#             label = '$r_4 =$ {:.5f}'.format(r4), linewidth = 0.9, color = 'r')
#plt.loglog(f_short, 
#             moving_average(SFS5, navg, start_smooth), 
#             label = '$r_5 =$ {:.5f}'.format(r5), linewidth = 0.9, color = 'y')
#
#plt.loglog(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2)
#plt.loglog(f_short, 
#           Un / s / f_short ** 2, label = '$p(f) = U_n / (s f^2)$', 
#           linestyle = '--', linewidth = 2)
#plt.loglog(f_short, np.ones(len(f_short)) * L/ v, linewidth = 2, linestyle = '--'
#              , label = '$p(f) = U_n L / v$', color = 'b')
#
#
#
##plt.semilogy(f_short, 
##             moving_average(SFS, navg, start_smooth), 
##             label = '$r_1 =$ {:.5f}'.format(r), linewidth = 0.9, color = 'b')
##plt.semilogy(f_short, 
##             moving_average(SFS2, navg, start_smooth), 
##             label = '$r_2 =$ {:.5f}'.format(r2), linewidth = 0.9, color = 'k')
##plt.semilogy(f_short, 
##             moving_average(SFS3, navg, start_smooth), 
##             label = '$r_3 =$ {:.5f}'.format(r3), linewidth = 0.9, color = 'g')
##
##plt.semilogy(f_short, 
##             moving_average(SFS4, navg, start_smooth), 
##             label = '$r_4 =$ {:.5f}'.format(r4), linewidth = 0.9, color = 'r')
##plt.semilogy(f_short, 
##             moving_average(SFS4, navg, start_smooth), 
##             label = '$r_5 =$ {:.5f}'.format(r5), linewidth = 0.9, color = 'y')
##
##plt.semilogy(f_short, 
##           2 * Un * N * L * np.ones(len(f_short)) / f_short,
##           label = '$p(f) = 2 N U_n / f$', linestyle = '--', linewidth = 2)
##plt.semilogy(f_short, 
##           Un / s / f_short ** 2, label = '$p(f) = U_n / (s f^2)$', 
##           linestyle = '--', linewidth = 2)
##plt.semilogy(f_short, np.ones(len(f_short)) * L/ v, linewidth = 2, linestyle = '--'
##              , label = '$p(f) = U_n L / v$', color = 'b')
##
#
#
#plt.vlines(1 / (N * v), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, color = 'b', label = '$f = 1 / N v$')
#
#plt.title('L = {}, '.format(L) + 
#          r'$\rho = $' + 
#          ' = {}, m = {:.2f}, s = {:.2f}, sample uniformly everywhere'.format(N, m, s))
#






