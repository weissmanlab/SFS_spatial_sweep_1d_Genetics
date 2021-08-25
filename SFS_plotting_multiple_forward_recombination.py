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

L = 2000
rho = 5000
N = rho * L
s = 0.05
m = 0.25
tfinal = 100000
Un = 1
r = 0
n_forward = 100
tfix = 0
nSFS = 10000

# Change r - sample everywhere
n = 100000
start_smooth = 100
navg = 100

rlist = [0.000001, 0]

f = np.arange(1, n + 1) / n

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

cividis_cmap = cm.get_cmap('cividis')


plt.figure(figsize = (24, 18))
plt.xlabel(r'Frequency, $f$', fontsize = 75)
plt.ylabel(r'Number of alleles, $P(f)$', fontsize = 75)

plt.loglog(f, 
           2 * Un * N * np.ones(len(f)) / f,
           label = r'$P(f) = 2 N U / f$', linestyle = '--', linewidth = 6, 
           color = '#cc79a7')


f_short2 = np.linspace(3 / (rho * v0), 1, 100)
plt.vlines(1 / (rho * v0), 10 ** 3, 10 ** 11, linestyle = 'dotted',
           linewidth = 6, color = '#009e73', label = r'$f = 1 / \rho v$')
plt.text(10 ** (-3) / 2, 8 * 10 ** 11, r'$P(f) = 2.5 U_{eff} / s f^2, U_{eff} = U (1 + 2 N r) $')

for rind in range(len(rlist)):
    r = rlist[rind]
    Uneff = Un * (1 + 2 * N * r) 
#    if rind == 7:
#        start_smooth = 300
#        nSFS = 1000
#        navg = 30
#    elif rind > 2 :
#        start_smooth = 100
#        nSFS = 10000
#        navg = 200
#    else:
#        start_smooth = 300
#        nSFS = 2000
#        navg = 100

    SFS = np.zeros(n)
    f_short = moving_average(f, navg, start_smooth)
    for i in range(n_forward):
            SFS += n * np.loadtxt(
    'backward simulation data/expected_SFS_L=' 
    + '{}_N={}_s={:.6f}_m={:.6f}_r={:.7f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L, 
             rho, s, m, r, tfinal, n, tfix, nSFS, i))
    
    SFS /= n_forward
    plt.loglog(f_short, 
             moving_average(SFS, navg, start_smooth), 
             label = '$r =$ {:.1e}'.format(r), linewidth = 2, color = 
             cividis_cmap(rind / len(rlist)), alpha = 0.8)
    print('{:.1e}'.format(r))
    print(SFS[-2:])
    if r < 5 * 10 ** (-3):
        plt.loglog(f_short2, 
           2.5 * Uneff / s / f_short2 ** 2, 
           linestyle = '-.', color = cividis_cmap(rind / len(rlist)), 
           linewidth = 6, alpha = 0.8)
        plt.loglog(f_short, Uneff * np.ones(len(f_short)) * L / v0,
                   linewidth = 6, linestyle = '--'
                   , color = cividis_cmap(rind / len(rlist)), alpha = 0.8)



plt.legend(fontsize = 'small', loc = 'lower left')

