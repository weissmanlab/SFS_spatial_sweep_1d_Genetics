#!/usr/bin/env python3
"""
Coalescent simulation with recombination.
open the frequency file read backward in time until n individuals all coalesce.
As we go, we will record the number of leaves for each node. 
If the neutral mutation arises as a Poisson process with rate Un = 1, 
the expected AFS is the histogram of the leaf counts 

Update : pre-sweep coalescent process changed to account for more than one 
mergers in a single generation when there are many individuals are left. 
(compared to sqrt(2 * N))


"""
import numpy as np
from multiprocessing import Pool
import sys
from numpy import random

"""
L = int(sys.argv[1]) # number of demes
N = int(sys.argv[2]) # deme capacity
s = float(sys.argv[3]) # selection coef
m = float(sys.argv[4]) # migration rate
r = float(sys.argv[5]) # recombination rate
tfinal = int(sys.argv[6]) # sweep time
nbase = int(sys.argv[7]) # sample size
N_SFS = int(sys.argv[8]) # number of coalescent simulation we run.
T_after_fix = int(sys.argv[9]) # number of generations between fixation and sampling
n_forward = int(sys.argv[10]) # forward time simulation number
"""
L = 50 # number of demes in 1 D
N = 100 # deme capacity
s = 0.5 # selection coef
m = 0.1 # migration rate
r = 0 # recombination rate
tfinal = 528 # sweep time
nbase = 10 # sample size
N_SFS = 1 # number of coalescent simulation we run.
T_after_fix = 10 # number of generations between fixation and sampling
n_forward = 1 # forward time simulation number

fname = 'Data_Reshape_528.txt'
freq_file = open(fname)
lines = np.loadtxt(fname, dtype=np.int64)
lines = lines.reshape(tfinal, L,L)  ###lines = Data from forward simulation
"""Flatten lines[i] into 1-D array here"""

# a/b, except = 0 when b = 0
# https://stackoverflow.com/questions/26248654/how-to-return-0-with-divide-by-zero
def safe_divide(a, b, val=0):
    return np.divide(a, b, out=np.full_like(a, val), where=b != 0)


# Coalescence before the sweep (neutral L * N pop)

def get_parent_presweep_arr_new(inds):
    """
    Function picks up a particular individiual, allows it to randomly migrate and tracks it to the parent till before the sweep happened
    
    Input Arguement: inds (created in main function runner from the results of the forward simulation)
    inds = [type of mutation (0 for wt), which deme number, individual within the deme]
    mut_type = 1 for neutral mutation, 0 for wt
    
    """
    mut_types, deme_arr, ind_in_deme_arr = inds.T
    mut_types_new = (np.zeros_like(mut_types)).astype(np.int64)
    deme_arr_new = np.zeros_like(deme_arr)  ##Where the mutation goes after migration
    ind_in_deme_arr_new = np.zeros_like(ind_in_deme_arr)


    len_inds = len(inds)
    # array of random numbers used to determine which deme the parent is in.
    which_deme_rand = random.random(len_inds)
    # random numbers used to determine which individual to pick in each deme.
    choice_rand = random.random(len_inds)  
    """Should be rho?"""
    
    
    """Creating Cummulative Probablities of each direction"""
    left_range_prob = m / 4. # prob of going to the left deme
    right_range_prob = 2*m / 4. # prob of going to the right deme
    top_range_prob = 3*m / 4. # prob of going to the upper deme
    bottom_range_prob = 4*m / 4. # prob of going to the bottom deme
    mid_range_prob = 1. - m # Prob of remaining in the deme
    # Sum of these 5 probablities is 1


    #Monte Carlo like method of choosing which deme the migration happens to.
    #Only one of these cases will evaluate to true for a given which_deme_rand and hence only one case will be executed. 

    left_idxs = np.where(which_deme_rand < left_range_prob)[0]
    deme_arr_new[left_idxs] = (deme_arr[left_idxs] + [(L-1) if deme_arr[left_idx]%L == 0 else -1][0] ).astype(np.int64)
    
    right_idxs = np.where(np.logical_and(
        which_deme_rand > left_range_prob,
        which_deme_rand < right_range_prob))[0]
    deme_arr_new[right_idxs] = (deme_arr[right_idxs] + [-(L-1) if (deme_arr[right_idx]+1)%L == 0 else 1][0] ).astype(np.int64)

    top_idxs = np.where(np.logical_and(
        which_deme_rand > right_range_prob,
        which_deme_rand < top_range_prob))[0]
    deme_arr_new[top_idxs] = ((deme_arr[top_idxs]-L)%(L*L)).astype(np.int64)

    bottom_idxs = np.where(np.logical_and(
        which_deme_rand > top_range_prob,
        which_deme_rand < bottom_range_prob))[0]
    deme_arr_new[bottom_idxs] = ((deme_arr[bottom_idxs]+L)%(L*L)).astype(np.int64)

    mid_idxs = np.where(which_deme_rand > bottom_range_prob)[0]
    deme_arr_new[mid_idxs] = (deme_arr[mid_idxs]).astype(np.int64)

    """
    Note: The addition and subtraction of indexes is to do migration in 2-D space in a 1-D array. 
    For index i:
    Top = (i-L)%(L*L)
    Bottom = (i+L)%(L*)
    Right = i + [-(L-1) if (i+1)%L ==0 else 1][0]
    Left = i + [(L-1) if (i)%L ==0 else -1][0]
    
    This takes care of edge cases also
    """

    ind_in_deme_arr_new = (np.floor(choice_rand * N)).astype(np.int64)

    inds2 = np.vstack((mut_types_new,
                              deme_arr_new,
                              ind_in_deme_arr_new)).T
    return inds2



def get_individuals2_new(Ne, Ne_parent, individuals):
    '''
    for each bucket, choose between neighboring/own buckets w/ relative
    probabilities [m/2, 1-m] respectively.
    
    Input Arguements: Ne, Ne_Parents (Arrays of populations), individuals which is same as inds previously and has three parameters (mutation, which deme, number within deme)
    '''
    
    mut_types, deme_arr, ind_in_deme_arr = individuals.T
    deme_arr = (deme_arr).astype(np.int64)

    mut_types_next = np.ones_like(mut_types)
    p_vals = random.random(len(mut_types))

    deme_arr_next = np.zeros_like(deme_arr)
    ind_in_deme_arr_next = np.zeros_like(ind_in_deme_arr)
    Ne_parent_extended = np.concatenate(([Ne_parent[0]], Ne_parent, [Ne_parent[-1]]))
    """Why this step??"""


    len_inds = len(deme_arr_next)

    # For mutant first, find the parent's deme and then its index inside the deme
    left_parent_prob = m / 2 * np.take(Ne_parent_extended, deme_arr)
    mid_parent_prob = (1 - m) * np.take(Ne_parent_extended, deme_arr + 1)
    right_parent_prob = m / 2 * np.take(Ne_parent_extended, deme_arr + 2)
    total_prob = (left_parent_prob + mid_parent_prob + right_parent_prob)

    # Set the cumulative probability
    mid_parent_prob = safe_divide(left_parent_prob + mid_parent_prob,total_prob,val=1,)
    left_parent_prob = safe_divide(left_parent_prob, total_prob)


    which_parent_rand = random.random(len_inds) # to choose btw left/mid/right/top/bottom
    choice_rand = random.random(len_inds) # used for index within the deme
    """Should be rho?"""

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
    
    individuals2 = np.vstack((mut_types_next,
                              deme_arr_next,
                              ind_in_deme_arr_next)).T
    return individuals2


"""
def runner(idx):
    Ne_0 = lines[-1]
    n = nbase
    SFS = np.zeros(n)

    Ne = Ne_0
    Ne = (Ne).astype(np.int64)
    # Ne = [range(Ne_0[i]) for i in range(len(Ne_0))]
    individuals = [] # format will be [mut_type, deme index, individual index (inside the deme)]
    if n < round(sum(Ne)):
        individuals_location = random.choice(np.arange(0
        , round(sum(Ne))), size = n, replace = False)

    else:
        individuals_location = np.arange(0
        , round(sum(Ne)))
        n = sum(Ne)
    ind_inside_deme = np.mod(individuals_location, N)
    deme_ind = (individuals_location - ind_inside_deme) // N

    for k in range(n):
        # 1 is for having beneficial mutation
        individuals.append([1, deme_ind[k], ind_inside_deme[k]])
        
    individuals = np.array(individuals)
    unique, leaf_counts = np.unique(individuals, axis = 0, return_counts = True)
    hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))


    t_after_fix = 0
    while t_after_fix < T_after_fix:
        t_after_fix += 1
        individuals2 = get_individuals2_new(Ne, Ne, individuals)
        individuals2 = np.repeat(individuals2, leaf_counts, axis = 0)
        unique, leaf_counts = np.unique(individuals2, axis = 0,
                                         return_counts = True)
        
        hist, bin_edges = np.histogram(leaf_counts,
                                       bins = np.arange(1, n + 2))
        
        SFS += hist
        individuals = unique



    line_num = -1
    while (len(individuals) > 1) and (line_num > -len(lines)):
        line_num -= 1
        Ne_parent = (lines[line_num]).astype(np.int64)
        individuals2 = get_individuals2_new(Ne, Ne_parent, individuals)
        individuals2 = np.repeat(individuals2, leaf_counts, axis = 0)
        unique, leaf_counts = np.unique(individuals2, axis = 0,
                                                        return_counts = True)
        hist, bin_edges = np.histogram(leaf_counts,
                                       bins = np.arange(1, n + 2))

        SFS += hist
        individuals = unique
        Ne = Ne_parent

    # The first round of coalescence simulation ends when we get to the time
    # when the beneficial mmutation arose, and all the left over individuals
    # are WT. From this point on, the coalescence will be extremely slow.
    # Therefore, we will run the same kind of simulation until the individuals
    # disperse for L^2 / m / N. We speed up the simulation by recording the number
    # of generations between the coalescence events.
    left_individuals = len(individuals)
    branch_len = 0 # number of generations until first merging event
    extra_gen = 0 # extra run time before stopping the coalescent simulation.

    while left_individuals > 1 and extra_gen < int(L ** 2 / m / N):

        branch_len += 1
        extra_gen += 1

        individuals2 = get_parent_presweep_arr_new(individuals)
        individuals2 = np.repeat(individuals2, leaf_counts, axis = 0)

        unique, leaf_counts = np.unique(individuals2,
                                        axis = 0, return_counts = True)
        current_individuals_counts = len(unique)
        individuals2 = unique

        if left_individuals == current_individuals_counts:
            individuals = individuals2
        else:
            hist, bin_edges = np.histogram(leaf_counts * branch_len,
                                       bins = np.arange(1, n + 2))

            SFS += hist
            individuals = unique
            branch_len = 0
        left_individuals = len(individuals)
        
        
        
    # After the individuals disperse, we may ignore the spatial structure.
    # If there are more than sqrt(2 * N * L) individuals left, it is likely
    # to have more than 1 merging event in a single generation. Thus, we manually
    # draw parent for every individual and find if some of them are the same.
    # Once there are much smaller number of individuals left, it will take a long
    # time to coalesce. Thus, we draw T2 (time until the first coalescent) from
    # geometric prob. distribution. (This is found from prob. of choosing different
    # parent for every individual left.)
    
    while left_individuals > 1:
        T2 = 1 + random.exponential(
                2 * N * L / ((left_individuals - 1) * left_individuals / 2))
        
        hist, bin_edges = np.histogram(leaf_counts,
                                       bins = np.arange(1, n + 2))
        SFS += hist * T2
        
        coal_inds = random.choice(range(left_individuals), 
                                  size = 2, replace = False)

        individuals[coal_inds[0]] = individuals[coal_inds[1]] 

        individuals2 = np.repeat(individuals, leaf_counts, axis = 0)

        unique, leaf_counts = np.unique(individuals2, 
                                        axis = 0, return_counts = True)
        individuals = unique
        left_individuals = len(individuals)
        # print(left_individuals)


    f = np.arange(1, n + 1) / n
    H = np.sum(2 * f * (1 - f) * SFS) / np.sum(SFS)


    if np.mod(idx, 500) == 0:
        np.savetxt('expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_n={}_{}.txt'.format(L,
           N, s, m, r, tfinal, n, T_after_fix, n_forward, idx), SFS)
    return SFS, H

if __name__ == '__main__':

        # this is the true sample number in case Ne < nbase.
    # print(individuals)
    p = Pool(20)
    
    ret = p.map(runner, range(N_SFS))
    SFS_items = [r[0] for r in ret]
    H_items = [r[1] for r in ret]
    SFS = np.sum(SFS_items, axis=0)
    SFS /= N_SFS
    np.savetxt('expected_SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L,
                N, s, m, r, tfinal, nbase, T_after_fix, N_SFS, n_forward), SFS)
    np.savetxt('expected_H_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L,
                N, s, m, r, tfinal, nbase, T_after_fix, N_SFS, n_forward), H_items)
    
    
    # SFS_sum = runner(0)
    # for i in range(N_SFS - 1):
    #     SFS = runner(0)
    #     SFS_sum += SFS
    # SFS_avg = SFS_sum / N_SFS
    # np.savetxt('SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L,
    #             N, s, m, r, tfinal, nbase, Un, N_SFS), SFS)
    
    
"""