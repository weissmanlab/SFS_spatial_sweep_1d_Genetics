# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 19:01:27 2021

@author: jim903
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=100)


L = 500
N = 20000
s = 0.05
m = 0.25
tfinal = 1000000
freq_file = open(
   'forward simulation/L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s, m, tfinal))
lines = np.loadtxt(freq_file, dtype=np.int64)
tf_real = len(lines)
v_list = []

plt.figure(figsize = (20, 8))
plt.plot(lines[1100] / N, linewidth = 10, color = 'k')
plt.xlim((90, 160))
plt.ylim((-0.01, 1.02))
#plt.gca().axes.get_xaxis().set_visible(False)
plt.xticks([])
plt.xlabel(r'$x$')
plt.ylabel(r'$p(x)$')
plt.show()
#for t in np.arange(0, tf_real, 100):
#    plt.plot(lines[t])
