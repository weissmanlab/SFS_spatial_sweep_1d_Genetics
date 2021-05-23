#!/usr/bin/env python3
"""
Coalescent simulation with recombination.
open the frequency file read backward in time until n individuals all coalesce.
As we go, we will record the number of leaves for each node. Then we sample
some of the leaf counts with rate Un, since a neutral mutation arises as a
Poisson process. We append the counts until we reach a single common ancestor.
The histogram of the counts is SFS




Todo: fix after time zero. (get parent pre sweep)


"""



import numpy as np
from multiprocessing import Pool
import sys
from numpy import random


L = 5000 # number of demes
N = 2000 # deme capacity
s = 0.05 # selection coef
m = 0.25 # migration rate
tfinal = 100000 # sweep time
Un = 1 # rate of neutral mutation
nbase = 1000 # sample size
N_SFS = 1 # number of coalescent simulation we run.
L_start = 0 # sampling range 
L_end = 500

#L = int(sys.argv[1]) # number of demes
#N = int(sys.argv[2]) # deme capacity
#s = float(sys.argv[3]) # selection coef
#m = float(sys.argv[4]) # migration rate
#tfinal = int(sys.argv[5]) # sweep time
#Un = float(sys.argv[6]) # rate of neutral mutation
#nbase = int(sys.argv[7]) # sample size
#N_SFS = int(sys.argv[8]) # number of coalescent simulation we run.
#L_start = int(sys.argv[9]) # sampling range 
#L_end = int(sys.argv[10])

fname = 'L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_0.txt'.format(L, N, s, m, tfinal)
freq_file = open(fname)
lines = np.loadtxt(fname, dtype=np.int64)

# a/b, except = 0 when b = 0
# https://stackoverflow.com/questions/26248654/how-to-return-0-with-divide-by-zero
def safe_divide(a, b, val=0):
    return np.divide(a, b, out=np.full_like(a, val), where=b != 0)

# Coalescence before the sweep (neutral L * N pop)

def get_parent_presweep_arr_new(inds):


    mut_types, deme_arr, ind_in_deme_arr = inds.T
    mut_types_new = (np.zeros_like(mut_types)).astype(np.int64)
    deme_arr_new = np.zeros_like(deme_arr)
    ind_in_deme_arr_new = np.zeros_like(ind_in_deme_arr)

    len_inds = len(inds)
    which_deme_rand = random.random(len_inds)
    choice_rand = random.random(len_inds)
    left_range_prob = m / 2
    mid_range_prob = 1 - m / 2

    left_idxs = np.where(which_deme_rand < left_range_prob)[0]
    deme_arr_new[left_idxs] = (deme_arr[left_idxs] - 1).astype(np.int64)

    mid_idxs = np.where(np.logical_and(
        which_deme_rand > left_range_prob,
        which_deme_rand < mid_range_prob))[0]
    deme_arr_new[mid_idxs] = (deme_arr[mid_idxs]).astype(np.int64)

    right_idxs = np.where(which_deme_rand > mid_range_prob)[0]
    deme_arr_new[right_idxs] = (deme_arr[right_idxs] + 1).astype(np.int64)

    left_edge_idxs = np.where(deme_arr_new < 0)[0]
    deme_arr_new[left_edge_idxs] = 0

    right_edge_idxs = np.where(deme_arr_new > L - 1)[0]
    deme_arr_new[right_edge_idxs] = L - 1

    ind_in_deme_arr_new = (np.floor(choice_rand * N)).astype(np.int64)

    inds2 = np.vstack((mut_types_new,
                              deme_arr_new,
                              ind_in_deme_arr_new)).T
    return inds2

def get_individuals2_new(Ne, Ne_parent, individuals):
    r = 0
    '''
    for each bucket, choose between neighboring/own buckets w/ relative
    probabilities [m/2, 1-m] respectively.
    '''
    # mut_types, i_arr = individuals.T
    mut_types, deme_arr, ind_in_deme_arr = individuals.T
    deme_arr = (deme_arr).astype(np.int64)

    Nwt = (N - Ne).astype(np.int64)
    Nwt_parent = (N - Ne_parent).astype(np.int64)
    mut_types_next = np.ones_like(mut_types)
    p_vals = random.random(len(mut_types))

    # There are two ways in which an individual be a wt:
    # First, it could be already a mutant and recombine with a WT
    zero_idx = np.where(np.logical_and(
        mut_types == 1,
        p_vals < r * np.take(Nwt, deme_arr) / (N - 1)))[0]
    mut_types_next[zero_idx] = 0

    # Second, it can be a WT originally and NOT recombine with a mutant in the same deme
    zero_idx2 = np.where(np.logical_and(
        mut_types == 0,
        p_vals > r * np.take(Ne, deme_arr) / (N - 1)))[0]
    mut_types_next[zero_idx2] = 0


    deme_arr_next = np.zeros_like(deme_arr)
    ind_in_deme_arr_next = np.zeros_like(ind_in_deme_arr)
    Ne_parent_extended = np.concatenate(([Ne_parent[0]], Ne_parent, [Ne_parent[-1]]))
    Nwt_parent_extended = np.concatenate(([Nwt_parent[0]], Nwt_parent, [Nwt_parent[-1]]))

    len_inds = len(deme_arr_next)

    # For mutant first, find the parent's deme and then its index inside the deme
    left_parent_prob = m / 2 * np.take(Ne_parent_extended, deme_arr)
    mid_parent_prob = (1 - m) * np.take(Ne_parent_extended, deme_arr + 1)
    right_parent_prob = m / 2 * np.take(Ne_parent_extended, deme_arr + 2)
    total_prob = (left_parent_prob + mid_parent_prob + right_parent_prob)

    # Set the cumulative probability
    mid_parent_prob = safe_divide(
        left_parent_prob + mid_parent_prob,
        total_prob,
        val=1,
    )
    left_parent_prob = safe_divide(left_parent_prob, total_prob)


    which_parent_rand = random.random(len_inds) # choose btw left/mid/right
    choice_rand = random.random(len_inds) # used for index within the deme


    left_parent_idxs = np.where(np.logical_and(
        which_parent_rand < left_parent_prob,
        mut_types_next == 1))[0]
    deme_arr_next[left_parent_idxs] = (deme_arr[left_parent_idxs] - 1).astype(np.int64)



    mid_parent_idxs = np.where(np.logical_and.reduce((
        which_parent_rand > left_parent_prob,
        which_parent_rand < mid_parent_prob,
        mut_types_next == 1)))[0]
    deme_arr_next[mid_parent_idxs] = (deme_arr[mid_parent_idxs]).astype(np.int64)


    right_parent_idxs = np.where(np.logical_and(
        which_parent_rand > mid_parent_prob,
        mut_types_next ==1))[0]
    deme_arr_next[right_parent_idxs] = (deme_arr[right_parent_idxs] + 1).astype(np.int64)

    left_edge_idxs = np.where(deme_arr_next < 0)[0]
    deme_arr_next[left_edge_idxs] = (np.zeros_like(left_edge_idxs)).astype(np.int64)

    right_edge_idxs = np.where(deme_arr_next > L - 1)[0]
    deme_arr_next[right_edge_idxs] = (np.ones_like(right_edge_idxs) *
                                      (L - 1)).astype(np.int64)

    mut_idxs = np.concatenate((left_parent_idxs, mid_parent_idxs, right_parent_idxs))
    ind_in_deme_arr_next[mut_idxs] = (np.floor(
        choice_rand[mut_idxs] * np.take(Ne_parent,
                              deme_arr_next[mut_idxs]))).astype(np.int64)

    # then same for wt
    left_parent_prob = m / 2 * np.take(Nwt_parent_extended, deme_arr)
    mid_parent_prob = (1 - m) * np.take(Nwt_parent_extended, deme_arr + 1)
    right_parent_prob = m / 2 * np.take(Nwt_parent_extended, deme_arr + 2)
    total_prob = (left_parent_prob + mid_parent_prob + right_parent_prob)
    mid_parent_prob = safe_divide(
        left_parent_prob + mid_parent_prob,
        total_prob,
        val=1,
    )
    left_parent_prob = safe_divide(left_parent_prob, total_prob)

    which_parent_rand = random.random(len_inds) # choose btw left/mid/right
    choice_rand = random.random(len_inds) # used for index within the deme


    left_parent_idxs = np.where(np.logical_and(
        which_parent_rand < left_parent_prob,
        mut_types_next == 0))[0]
    deme_arr_next[left_parent_idxs] = (deme_arr[left_parent_idxs] - 1).astype(np.int64)

    mid_parent_idxs = np.where(np.logical_and.reduce((
        which_parent_rand > left_parent_prob,
        which_parent_rand < mid_parent_prob,
        mut_types_next == 0)))[0]
    deme_arr_next[mid_parent_idxs] = (deme_arr[mid_parent_idxs]).astype(np.int64)


    right_parent_idxs = np.where(np.logical_and(
        which_parent_rand > mid_parent_prob,
        mut_types_next == 0))[0]
    deme_arr_next[right_parent_idxs] = (deme_arr[right_parent_idxs] + 1).astype(np.int64)

    left_edge_idxs = np.where(deme_arr_next < 0)[0]
    deme_arr_next[left_edge_idxs] = (np.zeros_like(left_edge_idxs)).astype(np.int64)

    right_edge_idxs = np.where(deme_arr_next > L - 1)[0]
    deme_arr_next[right_edge_idxs] = (np.ones_like(right_edge_idxs) * (L - 1)).astype(np.int64)

    wt_idxs = np.concatenate((left_parent_idxs, mid_parent_idxs, right_parent_idxs))
    ind_in_deme_arr_next[wt_idxs] = (np.floor(
        choice_rand[wt_idxs] * np.take(Nwt_parent,
                              deme_arr_next[wt_idxs]))).astype(np.int64)
    individuals2 = np.vstack((mut_types_next,
                              deme_arr_next,
                              ind_in_deme_arr_next)).T
    return individuals2

def runner(idx):
    Ne_0 = lines[-1]
    n = nbase
    ancestry = []
#    SFS = np.zeros(n)

    Ne = Ne_0
    Ne = (Ne).astype(np.int64)
    # Ne = [range(Ne_0[i]) for i in range(len(Ne_0))]
    individuals = [] # format will be [mut_type, deme index, individual index (inside the deme)]
    if n < (L_end - L_start) * N:
        individuals_location = random.choice(np.arange(L_start * N
        , L_end * N), size = n, replace = False)

    else:
        individuals_location = np.arange(L_start * N
        , L_end * N)
        n = (L_end - L_start) * N
    ind_inside_deme = np.mod(individuals_location, N)
    deme_ind = (individuals_location - ind_inside_deme) // N

    for k in range(n):
        # 1 is for having beneficial mutation
        individuals.append([1, deme_ind[k], ind_inside_deme[k]])
    individuals = np.array(individuals)
    unique, leaf_counts = np.unique(individuals, axis = 0, return_counts = True)
    
    hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))

    line_num = -1
    while (len(individuals) > 1) and (line_num > -len(lines)):
        line_num -= 1
        Ne_parent = (lines[line_num]).astype(np.int64)
        individuals2 = get_individuals2_new(Ne, Ne_parent, individuals)
        individuals2 = np.repeat(individuals2, leaf_counts, axis = 0)
        unique, leaf_counts = np.unique(individuals2, axis = 0,
                                                        return_counts = True)
        mut_types, deme_arr, ind_in_deme_arr = unique.T
        ancestry.append(list(np.repeat(deme_arr, leaf_counts)))
        
#        neutal_mut_counts = random.poisson(Un, len(unique))
#        descendents_counts = np.repeat(leaf_counts, neutal_mut_counts)
#        hist, bin_edges = np.histogram(descendents_counts,
#                                       bins = np.arange(1, n + 2))
#        SFS += hist
        individuals = unique
        Ne = Ne_parent
    print(len(ancestry))    
    np.savetxt('ancestry_tree_L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_nsample={}_Lstart={}_Lend={}_sample_uniform_{}.txt'.format(L,
           N, s, m, tfinal, n, L_start, L_end, idx), ancestry)


    
    return ancestry

if __name__ == '__main__':

        # this is the true sample number in case Ne < nbase.
    # print(individuals)
#    p = Pool(8)
#    
#    ret = p.map(runner, range(N_SFS))
    
    
    ancestry = runner(0)
    print(ancestry[-1][0])
