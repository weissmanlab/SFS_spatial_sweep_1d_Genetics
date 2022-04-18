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


# Change s - sample everywhere
n = 100000
nSFS = 1000

plasma_cmap = cm.get_cmap('plasma')
viridis_cmap = cm.get_cmap('viridis')

slist = [0.02, 0.03, 0.04, 0.05, 0.06]
mlist = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

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
    # Find v from lines

#    freq_file = open(
#       'forward_simulation_data/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, rho, s, m, tfinal))
#    lines = np.loadtxt(freq_file, dtype=np.int64)
#    tf_real = len(lines)
#    v_list = []
#    for t in np.arange(int(tf_real / 4), int(tf_real * 3 / 4)):
#        line1 = lines[t]
#        line2 = lines[t + 1]
#        psum1 = sum(line1) / N
#        psum2 = sum(line2) / N
#        v_list.append(psum2 - psum1)
#    
#    v = np.average(v_list)
    SFS = np.zeros(n)

    for i in range(n_forward):
        SFS += n * np.loadtxt(
        'backward_simulation_data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
                 rho, s, m, r, tfinal, n, tfix, nSFS, i))

    
    SFS /= n_forward


    plt.loglog(f_short, 
                 moving_average(SFS * s / np.log(N * s * f), 
                                navg, start_smooth), 
                 linewidth = 5, 
                 color = plasma_cmap(sind * 0.2), 
                 label = '$s = ${:.2f}'.format(s))

#    plt.loglog(f_short, np.ones(len(f_short)) * L / v, linewidth = 2, 
#               linestyle = '-.'
#                  , label = '$p(f) = U_n L / v, s = ${:.2f}'.format(s), 
#                  color = colorlist[sind])
#    plt.vlines((tfix + L / v) / (N * L), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#               linewidth = 5, color = plasma_cmap(sind * 0.2))
#    plt.text(xpositions[sind], ypositions[sind], 
#             r'$s = $' + '{:.2f}'.format(s), 
#             color = plasma_cmap(sind * 0.2), fontsize = 50)    
#    def power_law_fit(x, a):
#        ''' Use this to fit log(rho sqrt(ms) f) f^-2 (intermediate f) '''
#        return a * np.log(rho * L * s * x) * x ** (-2)
#    
#    popt, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS[fit_range_ind])
#    height = popt[0]
#    fitlist.append(height)
#    plt.loglog(f_short, height * np.log(rho * L * s * f_short) * f_short ** (-2), 
#               linewidth = 5, color = plasma_cmap(sind * 0.2), linestyle = '--')

for mind in range(len(mlist)):
    s = 0.05
    m = mlist[mind]
    SFS = np.zeros(n)

    for i in range(n_forward):
        SFS += n * np.loadtxt(
        'backward_simulation_data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
                 rho, s, m, r, tfinal, n, tfix, nSFS, i))

    
    SFS /= n_forward


    plt.loglog(f_short, 
                 moving_average(SFS * s / np.log(N * s * f), 
                                navg, start_smooth), 
                 linewidth = 5, 
                 color = viridis_cmap(mind * 0.12), 
                 label = '$m = ${:.2f}'.format(m))

    


plt.xlim((10 ** -3, 10 ** -1))
plt.ylim((5 * 10, 2 * 10 ** 6))
plt.xlabel('Frequency, ' + r'$\boldmath{f}$', fontsize = 100)
plt.ylabel(r'$\boldmath{P(f) \cdot s / \ln(Nsf)}$', fontsize = 100)
#plt.legend(fontsize = 'medium', loc = 'upper right')
plt.savefig('Figure3b_new.pdf', format = 'pdf')


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