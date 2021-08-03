import numpy as np
import matplotlib.pyplot as plt

def moving_average(a, n = 100, start_smooth = 1000):
    a_start = a[:start_smooth]
    a_end = a[start_smooth:]
    ret = np.cumsum(a_end, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a_start, ret[n - 1:] / n))

def power_law_fit(x, a):
    ''' Use this to fit f^-2 (intermediate f) '''
    return a * x ** (-2)

L = 500*500
rho = 500
N = rho * L
s = 0.05
m = 0.25
Un = 1
r = 0
n_forward = 1
n = 100000
f = np.arange(1, n + 1) / n

v0 = 0.0049
plt.figure(figsize = (24, 18))
plt.xlabel('Frequency, f', fontsize = 25)
plt.ylabel('Number of alleles, P(f)', fontsize = 25)
plt.ylim((10e4,10e13))

plt.loglog(f, 
           2 * Un * N * np.ones(len(f)) / f,
           label = r'P(f) = 2 N U / f', linestyle = '--', linewidth = 3, 
           color = '#cc79a7')

#f_short2 = np.linspace(3 / (rho * v0), 1, 100)


Uneff = Un * (1 + 2 * N * r) 
start_smooth = 100
navg = 1

SFS = np.zeros(n)
f_short = moving_average(f, navg, start_smooth)
power = power_law_fit(SFS,1)
print(len(power))


for i in np.arange(0, n_forward):
    SFS += n*np.loadtxt(
    'expected_SFS_L=500_N=500_s=0.050_m=0.25_nsample=100000_tfix=30_sample_uniform_navg=4.txt') 
    SFS /= n_forward
    plt.loglog(f_short, moving_average(SFS, navg, start_smooth), label = 'Asexual Data', linewidth = 2)
    #plt.loglog(f_short2, Uneff / s / f_short2 ** 2, linestyle = '-.', linewidth = 1, alpha = 0.8)
    plt.loglog(f_short, Uneff * np.ones(len(f_short)) * L / v0, linewidth = 3, label = 'Constant UL/v', linestyle = '--')
    plt.loglog(f_short, (Uneff/s)* power_law_fit(f_short,1), linewidth = 3, linestyle = '--', label = '(1/f^2)')

plt.legend(fontsize = 'small', loc = 'lower left')
