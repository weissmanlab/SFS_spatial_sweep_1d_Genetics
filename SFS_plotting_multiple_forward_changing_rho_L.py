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

plt.rcParams.update({'font.size': 15})
plt.figure(figsize = (12, 9))
plt.xlabel('f')
plt.ylabel('p(f)')

L = 5001
rho = 200
s = 0.05
m = 0.25
tfinal = 100000
Un = 1
r = 0
n_forward = 5
tfix = 0
dx = 250

l0list = np.arange(0, L, dx)
# Change l0
n = 100000
nSFS = 1000

N = L * rho

f = np.arange(1, n + 1) / n
navg = 30
start_smooth = 100
f_short = moving_average(f, navg, start_smooth)



# Find v from lines
l0 = 0
freq_file = open(
   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_l0={}_0.txt'.format(L, rho, s, m, tfinal, l0))
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
for l0 in l0list:
    for i in range(n_forward):
        SFS += n * np.loadtxt(
                'backward simulation data/expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_l0={}_navg={}_{}.txt'.format(L, 
                    rho, s, m, r, tfinal, n, tfix, l0, nSFS, i))


SFS /= n_forward * len(l0list)

SFS_smooth = moving_average(SFS, navg, start_smooth)
plt.plot(f_short[-95000:], 
             SFS_smooth[-95000:], 
             linewidth = 0.9)
plt.plot(f_short[-95000:], 2 * Un * (1 - f_short[-95000:]) * L / v, linewidth = 2, 
           linestyle = '-.'
              , label = r'$p(f) = 2 U_n L (1 - f) / v$')

#plt.vlines(1 / (rho * v), 10 ** 3, 10 ** 11, linestyle = 'dotted',
#           linewidth = 1, 
#           label = r'$f = 1 / \rho v$')
#
#plt.plot(f_short[-10000:], 
#           Un / s / f_short[-10000:] ** 2, 
#           label = '$p(f) = U_n / (s f^2)$', 
#           linestyle = '--', linewidth = 2)

#plt.plot(f_short[-1000:], 
#           2 * Un * N / f_short[-1000:],
#           label = '$p(f) = 2 N U_n / f$', linestyle = 'dotted', linewidth = 2)



plt.title(r'$L = {:.0e}, \rho= {:.0e}, s = {:.2f}, m = {:.2f}, r = {:.2f}$, sample uniformly everywhere, 100 forward sims'.format(L, rho, s, m, r))

plt.legend(fontsize = 'small', loc = 'upper right')

