# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 14:00:10 2021

@author: jim903
"""

'''
Well mixed population
forward simulation
'''
import numpy as np

N = 10 ** 7
s = 0.05
n_forward = 100

def forward_sim():
    n_selected = 1
    n_selected_series = [n_selected]
    while n_selected < N and n_selected > 0:
        n_selected = np.random.binomial(N, n_selected / N 
                                        + s * n_selected / N
                                         * (1 - n_selected / N))
        n_selected_series.append(n_selected)
    return n_selected_series

for idx in range(n_forward):
    data = forward_sim()
    np.savetxt('well_mixed_forward_N={}_s={:.2f}_{}.txt'.format(
                N, s, idx), data) 