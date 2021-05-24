# -*- coding: utf-8 -*-
"""
Created on Thu May 20 15:36:20 2021

@author: jim903
"""

import numpy as np
import numpy.random as random

N = 10000000
nsample = 100000
left_individuals = nsample
individuals = np.arange(left_individuals)
leaf_counts = [1 for _ in range(left_individuals)]


SFS = np.zeros(nsample)

while left_individuals > 1:
    print(left_individuals)
    if left_individuals > 100:
        T2 = 1
        parents = random.choice(range(N), size = left_individuals)
        individuals2 = np.repeat(parents, leaf_counts, axis = 0)
    else:
        
        T2 = 1 + random.exponential(N / 
                (left_individuals * (left_individuals - 1) / 2))
        print(T2)
        coal_inds = random.choice(range(left_individuals), 
                              size = 2)
        individuals[coal_inds[0]] = individuals[coal_inds[1]] 
    
        individuals2 = np.repeat(individuals, leaf_counts, axis = 0)

    unique, leaf_counts = np.unique(individuals2, 
                                    axis = 0, return_counts = True)
    hist, bin_edges = np.histogram(leaf_counts * T2,
                                   bins = np.arange(1, nsample + 2))
    SFS += hist
    individuals = unique
    left_individuals = len(individuals)


import matplotlib.pyplot as plt
f = np.arange(1, nsample + 1) / nsample
plt.loglog(f, SFS)
plt.loglog(f, 2 * nsample/f)

