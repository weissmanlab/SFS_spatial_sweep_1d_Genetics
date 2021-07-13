# -*- coding: utf-8 -*-
"""
Created on Thu May 20 15:36:20 2021

@author: jim903
"""

import numpy as np
import numpy.random as random


def moving_average(a, n = 3, start_smooth = 30):
    a_start = a[:start_smooth]
    a_end = a[start_smooth:]
    ret = np.cumsum(a_end, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a_start, ret[n - 1:] / n))

N = 10000000
nsample = 10000
nsim = 10

def main():
    left_individuals = nsample
    individuals = np.arange(left_individuals)
    leaf_counts = np.array([1 for _ in range(left_individuals)])
    
    
    SFS = np.zeros(nsample)
    while left_individuals > 1:
        
        T2 = random.exponential(N / ((left_individuals - 1) * left_individuals / 2))

#        log_prob_no_coal = np.sum([np.log((N - j) / N) 
#                            for j in np.arange(1, left_individuals)])
#        p = 1 - np.exp(log_prob_no_coal)
#        T2 = random.geometric(p)        
        num_coal = 1
    
        
    
#        if T2 < 2:
#            T2 = 1
#            x = random.random()
#            p_num_coal = 1
#            while x < p_num_coal:
#                p_num_coal *= ((nsample - 2 * num_coal) * 
#                               (nsample - 2 * num_coal - 1) / 
#                               (2 * N) * (N - num_coal) / N)
#                num_coal += 1
#            
        if np.mod(left_individuals, 10000) == 0:
            print(left_individuals)
            print(T2)
#        if num_coal == 0:
#            num_coal = 1
        
        hist, bin_edges = np.histogram(leaf_counts,
                                       bins = np.arange(1, nsample + 2))
    
        SFS += hist * T2


        coal_inds = random.choice(range(left_individuals), 
                                  size = int(2 * num_coal), replace = False)
        for i in range(num_coal):
            individuals[coal_inds[2 * i]] = individuals[coal_inds[2 * i + 1]] 
        individuals2 = np.repeat(individuals, leaf_counts, axis = 0)
        unique, leaf_counts = np.unique(individuals2, 
                                        axis = 0, return_counts = True)
        individuals = unique
        left_individuals = len(individuals)
    
    SFS *= nsample
    return SFS


SFS_avg = np.zeros(nsample)
for _ in range(nsim):
    print(_)
    SFS_avg += main() / nsim
import matplotlib.pyplot as plt
f = np.arange(1, nsample + 1) / nsample
plt.figure(figsize = (12, 8))
plt.loglog(moving_average(f, 20, 20), moving_average(SFS_avg, 20, 20))
plt.loglog(f, 4 * N / f)

'''
I think the issue is that there can only be 2 branches merging in a single generation.
But if lef_individuals is large, there should be more than 1 pair merging in a generation.
This problem can be by-passed if we take a continuum limit to the population and 
draw T2 from an exponential distribution. That way 0 < T2 < 1 is possible.
'''
