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


plt.figure(figsize = (24, 18))
plt.xlabel(r'$f$', fontsize = 75)
plt.ylabel(r'$P(f)$', fontsize = 75)
plt.xlim((10**(-3), 10 ** (-1)))
plt.ylim((10 ** 3, 10 ** 9))
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

plasma_cmap = cm.get_cmap('plasma')

mlist = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
fitlist = []

f = np.arange(1, n + 1) / n
navg = 30
start_smooth = 100
f_short = moving_average(f, navg, start_smooth)

fit_range_ind = np.arange(500,2000)


for mind in range(len(mlist)):
    m = mlist[mind]
    # Find v from lines

    freq_file = open(
       'forward_simulation_data/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s, m, tfinal))
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
        'backward_simulation_data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
                 N, s, m, r, tfinal, n, tfix, nSFS, i))

    
    SFS /= n_forward


    plt.loglog(f_short, 
                 moving_average(SFS, navg, start_smooth), 
                 label = '$m =$ {:.2f}'.format(m), linewidth = 5, 
                 color = plasma_cmap(mind * 0.12))
#    plt.semilogy(f_short, np.ones(len(f_short)) * L / v, linewidth = 2, 
#               linestyle = '-.'
#                  , label = '$p(f) = U_n L / v, m = ${:.2f}'.format(m), 
#                  color = colorlist[mind])
#    plt.vlines(1 / (N * v), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#               linewidth = 5, color = plasma_cmap(mind * 0.1), 
#               label = r'$f = 1 / \rho v, m = ${:.2f}'.format(m))
    
    def power_law_fit(x, a):
        ''' Use this to fit log(N s f) f^-2 (intermediate f) '''
        return a * np.log(N * L * s * x) * x ** (-2)
    
    popt, pcov = curve_fit(power_law_fit, f[fit_range_ind], SFS[fit_range_ind])
    height = popt[0]
    fitlist.append(height)
    plt.loglog(f_short, 
               height * np.log(N * L * s * f_short) / f_short ** 2,  
               linestyle = '--', linewidth = 5, color = plasma_cmap(mind * 0.12))


#plt.semilogy(f_short, 
#           2 * Un * N * L * np.ones(len(f_short)) / f_short,
#           linestyle = 'dotted', linewidth = 5)
#plt.text(10 ** (-4), 10 ** 4, r'$p(f) = U_n \ln(Nsf) / (2sf^2)$', color = 'k')
#plt.text(10 ** (-4) * 3, 10 ** 11, r'$p(f) = 2 N U_n / f$', color = 'k')


#plt.title('L = {}, '.format(L) + r'$\rho =$ {}, s = {:.2f}, r = {:.2f}, sample uniformly everywhere, 100 forward sims'.format(N, s, r))

plt.legend(fontsize = 'small', loc = 'upper right')
plt.savefig('Figure3c.pdf', format = 'pdf')

plt.figure(figsize = (24, 18))
plt.plot(mlist, fitlist, 'o', ms = 60, label = 'simulation')
plt.plot(mlist, 1 / 2 / s * np.ones(len(mlist)), label = '$k = U_n / 2 s$', 
         linewidth = 10)
plt.xlabel('$m$')
plt.ylabel('$k$')
plt.ylim((5, 15))
plt.legend()
plt.title('fitting SFS to ' + '$P(f) = k \ln(Nsf) / f^2$')
plt.savefig('Figure3d.pdf', format = 'pdf')
