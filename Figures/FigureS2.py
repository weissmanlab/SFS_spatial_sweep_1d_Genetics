# -*- coding: utf-8 -*-
"""
Sample n (1000) individuals from [l-dl/2, l + dl/2] 
at the end of the sweep and track the location
of the lineages backward in time until they all
coalesce.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

viridis_cmap = cm.get_cmap('viridis')
n = 1000
L = 500
rho = 20000
s = 0.05
m = 0.25
dl = 50
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 70, weight = 'bold')

plt.figure(figsize = (24, 18))
plt.xlabel(r'Generations before fixation, $t$', fontsize = 100)
plt.ylabel(r'Average loacation of lineages, $x_{avg}$', fontsize = 60)

llist = np.arange(25, 500, 50)
final_T_list = []
final_loc_list = []
for l in llist:
    loc_lineages = np.loadtxt('backward_simulation_data/ancestry_tree_' + 
            'L={}_N={}_s={:.6f}_m={:.6f}'.format(L, rho, s, m) + 
            '_tfinal=1000000_' + 
            'nsample={}_Lstart={}_Lend={}'.format(n, int(l - dl/2), int(l + dl/2)) + 
            '_sample_uniform_0.txt')
    T = len(loc_lineages)
    avg_locs = []


    for t in range(T):
        avg_locs.append(np.average(loc_lineages[t]))
#        plt.scatter(np.ones(int(n / 10)) * t, loc_lineages[t, ::10], 
#                    alpha = 0.002, color = viridis_cmap(l / 500), 
#                    s = 3)

    plt.plot(range(T), avg_locs, alpha = 0.7, 
             linewidth = 2, color = viridis_cmap(l / 500))
    final_T_list.append(T - 1)
    final_loc_list.append(avg_locs[-1])

plt.scatter(final_T_list, final_loc_list, s = 350, marker = '*', 
            color = viridis_cmap(llist / 500), edgecolors = 'r')

for i in np.arange(5,len(llist)):
    label = '$x_0 \in $[{}, {}]'.format(int(llist[i] - dl / 2), int(llist[i] + dl / 2))
    plt.annotate(label, (final_T_list[i], final_loc_list[i]), 
                 textcoords="offset points", 
                 xytext = (-20, 0), fontsize = 30 )

plt.savefig('FigureS2.pdf', format = 'pdf', bbox_inches = 'tight')