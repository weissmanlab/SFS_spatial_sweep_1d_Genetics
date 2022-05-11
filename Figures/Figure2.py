'''
Generate plots for Figure 2:
    SFS of asexual 1D and well-mixed populations.

And Figure 3a which is a closer-up version of Figure 2 
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
rho = 20000
N = L * rho
s = 0.05
m = 0.25
tfinal = 1000000
Un = 1
r = 0
n_forward = 100
n = 100000
n_wellmixed = 10000
nSFS = 1000
nsim_wellmixed = 1000
tfix = 0

plt.xlim((1/n, 1))
plt.ylim((2, 5 * 10 ** 14))
f = np.arange(1, n + 1) / n
f_wellmixed = np.arange(1, n_wellmixed + 1) / n_wellmixed
navg = 30
start_smooth = 100
f_short = moving_average(f, navg, start_smooth)

# Find v from lines

freq_file = open(
   'forward_simulation_data/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, rho, s, m, tfinal))
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

SFS = np.zeros(n)
#### For tfix = 0 
#for i in range(n_forward):
#    SFS += n * np.loadtxt(
#   'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
#             N, s, m, r, tfinal, n, tfix, nSFS, i))

for i in range(n_forward):
    SFS += n * np.loadtxt(
    'backward_simulation_data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
                 rho, s, m, r, tfinal, n, tfix, nSFS, i))

SFS /= n_forward

SFS_well_mixed = np.loadtxt('backward_simulation_data/expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}_r={:.2e}_nsample={}_nsim={}.txt'.format(
        N, tfix, s, r, n_wellmixed, nsim_wellmixed))
SFS_well_mixed2 = np.loadtxt('backward_simulation_data/expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}_r={:.2e}_nsample={}_nsim={}.txt'.format(
        N, tfix, s, r, n, 20))

f_wellmixed_combined = np.append(f[:78], f_wellmixed[8:])
SFS_wellmixed_combined = np.append(SFS_well_mixed2[:78], SFS_well_mixed[8:])
#####################################

plt.vlines(1 / (rho * v), 10 ** (-2), 10 ** 15, linestyle = 'dotted',
           linewidth = 10, color = '#009e73')


plt.text(5 * 10 ** (-4), 10 ** 6, r'$\boldmath{f = 1 / \rho v}$', 
         color = '#009e73', fontsize = 150, rotation = 90, 
         backgroundcolor = 'w')

plt.text(10 ** (-1), 10 ** 5, 'spatial', color = 'k', 
         fontsize = 150, backgroundcolor = 'w')
plt.text(10 ** (-2), 60, 'well-mixed', color = 'k', alpha = 0.5, 
         fontsize = 150, backgroundcolor = 'w')

#plt.text(0.002, 3 * 10 ** 5, r'$\boldmath{U_n \ln( N s f) / (s f^2)}$'
#        , color = '#0072b2', fontsize = 150, rotation = -25)
#        


plt.loglog(moving_average(f_wellmixed_combined, navg, start_smooth) , 
             moving_average(SFS_wellmixed_combined, navg, start_smooth), 
             linewidth = 12, color = 'k', alpha = 0.4)

plt.loglog(f_short, 
             moving_average(SFS, navg, start_smooth), 
             linewidth = 12, color = 'k')

neutral, = plt.loglog(f_short, 
           2 * Un * rho * L * np.ones(len(f_short)) / f_short,
           label = r'$\boldmath{2 N U_n / f}$', linestyle = '--', 
           linewidth = 12,
           color = '#cc79a7')


bsc, = plt.loglog(f_short, 
           Un / s / f_short ** 2, 
           linestyle = '-.', linewidth = 12, color = '#0072b2', 
           label = r'$\boldmath{U_n / (sf^2)}$')


flat, = plt.semilogy(f_short, np.ones(len(f_short)) * L / v, linewidth = 12, 
           linestyle = 'dotted'
              , label = r'$\boldmath{U_n L / v}$', color = '#d55e00')

xvals = [0.02, 0.00005, 0.007]
#xvals = [0.6]
#plt.text(0.5, 10 ** 4, r'$\boldmath{U_n L / v}$', color = '#d55e00', 
#         fontsize = 180)
#
#plt.text(0.5, 150, r'$\boldmath{U_n / (sf^2)}$', color = '#0072b2',
#                                fontsize = 180)
plt.ylim((10, 2 * 10 ** 12))
labelLines(plt.gca().get_lines() ,xvals = xvals, fontsize = 180)

plt.savefig('Figure2_loglog.pdf', format = 'pdf', bbox_inches = 'tight')
##########################
#
#
plt.rc('font', family='serif', size = 150, weight = 'bold')
plt.figure(figsize = (48, 36))
plt.xlabel(r'Frequency, $f$')
plt.ylabel(r'Number of alleles, $P(f)$')
#
plt.loglog(moving_average(f_wellmixed, navg, start_smooth) , 
             moving_average(SFS_well_mixed, navg, start_smooth), 
             linewidth = 12, color = 'k', alpha = 0.4)

plt.loglog(f_short, 
             moving_average(SFS, navg, start_smooth), 
             linewidth = 12, color = 'k')
bsc, = plt.loglog(f_short, 
           Un / s / f_short ** 2, 
           linestyle = '-.', linewidth = 12, color = '#0072b2', 
           label = r'$\boldmath{U_n / (sf^2)}$')
bsc2, = plt.loglog(f_short, 
           Un * np.log(N * s * f_short) / s / f_short ** 2 / 2, 
           linestyle = '-.', linewidth = 12, color = 'b', 
           label = r'$\boldmath{U_n \ln(Nsf) / (2 sf^2)}$')

plt.xlim((5 * 10 ** -4, 2 * 10 ** -1))
plt.ylim((200, 2 * 10 ** 8))
plt.text(4 * 10 ** (-2), 2 * 10 ** 5, 'spatial', color = 'k', 
         fontsize = 150, backgroundcolor = 'w')
plt.text(2 * 10 ** (-2), 5 * 10 ** 2, 'well-mixed', color = 'k', alpha = 0.5, 
         fontsize = 150, backgroundcolor = 'w')
plt.text(5 * 10 ** (-3), 10 ** 7, r'$\boldmath{U_n \ln(Nsf) / (2 sf^2)}$',  
         color = 'b', fontsize = 150, backgroundcolor = 'w')
plt.text(1.2 * 10 ** (-3), 4 * 10 ** 5, r'$\boldmath{U_n / (sf^2)}$', 
         color = '#0072b2', 
         fontsize = 150, backgroundcolor = 'w')

plt.savefig('Figure3a.pdf', format = 'pdf', bbox_inches = 'tight')

#
