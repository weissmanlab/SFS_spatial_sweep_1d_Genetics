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
plt.rc('font', family='serif', size = 60, weight = 'bold')
plt.figure(figsize = (24, 18))
plt.xlabel(r'$f$', fontsize = 75)
plt.ylabel(r'$P(f)$', fontsize = 75)
L = 500
N = 20000
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

s0 = 0.02
s1 = 0.03
s2 = 0.04
s3 = 0.05
s4 = 0.06

slist = [0.02, 0.03, 0.04, 0.05, 0.06]
colorlist = ['y', 'r', 'c', 'g', 'm', 'b', 'k']
fitlist = []

f = np.arange(1, n + 1) / n
navg = 30
start_smooth = 100
f_short = moving_average(f, navg, start_smooth)

fit_range_ind = np.arange(200,1000)


for sind in range(len(slist)):
    s = slist[sind]
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
                 label = '$s =$ {:.2f}'.format(s), linewidth = 0.9, 
                 color = colorlist[sind])

#    plt.loglog(f_short, np.ones(len(f_short)) * L / v, linewidth = 2, 
#               linestyle = '-.'
#                  , label = '$p(f) = U_n L / v, s = ${:.2f}'.format(s), 
#                  color = colorlist[sind])
    plt.vlines(1 / (N * v), 10 ** 3, 10 ** 11, linestyle = 'dotted',
               linewidth = 1, color = colorlist[sind], 
               label = r'$f = 1 / \rho v, s = ${:.2f}'.format(s))
    
    popt, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS[fit_range_ind])
    fitlist.append(popt[0])

plt.loglog(f_short, 
           Un / s / f_short ** 2, 
           label = '$p(f) = U_n / (s f^2), s = ${:.2f}'.format(s), 
           linestyle = '--', linewidth = 2, color = colorlist[sind])
plt.loglog(f_short, 
           2 * Un * N * L * np.ones(len(f_short)) / f_short,
           label = '$p(f) = 2 N U_n / f$', linestyle = 'dotted', linewidth = 2)



plt.title('L = {}, '.format(L) + r'$\rho =$' 
          + ' {}, m = {:.2f}, r = {:.2f}, sample uniformly everywhere, 100 forward sims'.format(N, m, r))

plt.legend(fontsize = 'medium', loc = 'upper right')


plt.figure(figsize = (24, 18))
plt.loglog(slist, fitlist, 'o',  ms = 60, label = 'simulation')
plt.loglog(slist, 2.5 / np.array(slist), label = '$k = 2.5 U / s$', 
           linewidth = 10)
plt.xlabel(r'$s$')
plt.ylabel(r'$k$')
plt.legend()
plt.title('fitting AFS to $P(f) = k / f^2$')
