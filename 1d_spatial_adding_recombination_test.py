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


# L = 500
# rho = 5000
# s = 0.05
# m = 0.25
# r = 0
# tfinal = 10000
# nbase = 1000
# N_SFS = 10


L = int(sys.argv[1]) # number of demes
rho = int(sys.argv[2]) # deme capacity
s = float(sys.argv[3]) # selection coef
m = float(sys.argv[4]) # migration rate
r = float(sys.argv[5]) # recombination rate
tfinal = int(sys.argv[6]) # sweep time
nbase = int(sys.argv[7]) # sample size
N_SFS = int(sys.argv[8]) # number of coalescent simulation we run.
T_after_fix = int(sys.argv[9]) # number of generations between fixation and sampling
n_forward = int(sys.argv[10]) # forward time simulation number
fname = 'L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_{}.txt'.format(L, 
           rho, s, m, tfinal, n_forward)
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
    # array of random numbers used to determine which deme the parent is in.
    which_deme_rand = random.random(len_inds)
    # random numbers used to determine which individual to pick in each deme.
    choice_rand = random.random(len_inds)
    left_range_prob = m / 2 # prob of going to the left
    mid_range_prob = 1 - m / 2 # P(left) + P(middle)
    # P(left) + P(middle) + P(right) = 1

    left_idxs = np.where(which_deme_rand < left_range_prob)[0]
    deme_arr_new[left_idxs] = (deme_arr[left_idxs] - 1).astype(np.int64)

    mid_idxs = np.where(np.logical_and(
        which_deme_rand > left_range_prob,
        which_deme_rand < mid_range_prob))[0]
    deme_arr_new[mid_idxs] = (deme_arr[mid_idxs]).astype(np.int64)

    right_idxs = np.where(which_deme_rand > mid_range_prob)[0]
    deme_arr_new[right_idxs] = (deme_arr[right_idxs] + 1).astype(np.int64)


    # Taking care of the edge cases
    
    left_edge_idxs = np.where(deme_arr_new < 0)[0]
    deme_arr_new[left_edge_idxs] = 0

    right_edge_idxs = np.where(deme_arr_new > L - 1)[0]
    deme_arr_new[right_edge_idxs] = L - 1

    ind_in_deme_arr_new = (np.floor(choice_rand * rho)).astype(np.int64)

    inds2 = np.vstack((mut_types_new,
                              deme_arr_new,
                              ind_in_deme_arr_new)).T
    return inds2

def get_individuals2_new(rho_e, rho_e_parent, individuals):
    '''
    1. Recombination - count the number of mutants & WT tracked in each bucket,
    Choose a random number for each -> if < r, recombine
    They can recombine with anyone in the same bucket
    --> count the number of mutants & WT that are being tracked after the recombination
    2. Find the parents - accounts for migration. The parent has to be of 
    the same genotype(WT or MUT)
    '''
    
    # 1. recombination

    mut_types, deme_arr, ind_in_deme_arr = individuals.T
    
    mut_types = (mut_types).astype(np.int64)
    deme_arr = (deme_arr).astype(np.int64)

    rho_wt = (rho - rho_e).astype(np.int64)
    rho_wt_parent = (rho - rho_e_parent).astype(np.int64)

    mut_tracking_counts = (np.zeros(L)).astype(np.int64)
    wt_tracking_counts = (np.zeros(L)).astype(np.int64)
    
    for i in range(len(individuals)):
        if mut_types[i] == 1:
            mut_tracking_counts[deme_arr[i]] += 1
        else:
            wt_tracking_counts[deme_arr[i]] += 1
    

    n_recom_mut_out = [random.poisson(r * mut_tracking_counts[i]) for i in range(L)]
    n_recom_wt_out = [random.poisson(r * wt_tracking_counts[i]) for i in range(L)]
    
    n_recom_mut_in = [sum(random.choice(np.append(np.zeros(rho_wt[i]), 
                                                   np.ones(rho_e[i])), 
    size = n_recom_mut_out[i] + n_recom_wt_out[i])) for i in range(L)]

    n_recom_mut_out = (np.array(n_recom_mut_out)).astype(np.int64)
    n_recom_wt_out = (np.array(n_recom_wt_out)).astype(np.int64)
    n_recom_mut_in = (np.array(n_recom_mut_in)).astype(np.int64)
    
    mut_tracking_counts += -n_recom_mut_out + n_recom_mut_in
    wt_tracking_counts -= -n_recom_mut_out + n_recom_mut_in
    
    mut_types_next = []
    deme_arr_after_recom = []
    for i in range(L):
        mut_types_next.extend([1 for _ in range(mut_tracking_counts[i])])
        deme_arr_after_recom.extend([i for _ in range(mut_tracking_counts[i])])
        mut_types_next.extend([0 for _ in range(wt_tracking_counts[i])])
        deme_arr_after_recom.extend([i for _ in range(wt_tracking_counts[i])])
    
    
    mut_types_next = (np.array(mut_types_next)).astype(np.int64)
    deme_arr_after_recom = (np.array(deme_arr_after_recom)).astype(np.int64)        


# 2. migration and assigning parents
    deme_arr_next = np.zeros_like(deme_arr_after_recom)
    ind_in_deme_arr_next = np.zeros_like(ind_in_deme_arr)
    rho_e_parent_extended = np.concatenate(([rho_e_parent[0]], rho_e_parent, [rho_e_parent[-1]]))
    rho_wt_parent_extended = np.concatenate(([rho_wt_parent[0]], rho_wt_parent, [rho_wt_parent[-1]]))

    len_inds = len(deme_arr_next)

    # For mutant first, find the parent's deme and then its index inside the deme
    left_parent_prob = m / 2 * np.take(rho_e_parent_extended, deme_arr_after_recom)
    mid_parent_prob = (1 - m) * np.take(rho_e_parent_extended, deme_arr_after_recom + 1)
    right_parent_prob = m / 2 * np.take(rho_e_parent_extended, deme_arr_after_recom + 2)
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
    deme_arr_next[left_parent_idxs] = (deme_arr_after_recom[left_parent_idxs] - 1).astype(np.int64)



    mid_parent_idxs = np.where(np.logical_and.reduce((
        which_parent_rand > left_parent_prob,
        which_parent_rand < mid_parent_prob,
        mut_types_next == 1)))[0]
    deme_arr_next[mid_parent_idxs] = (deme_arr_after_recom[mid_parent_idxs]).astype(np.int64)


    right_parent_idxs = np.where(np.logical_and(
        which_parent_rand > mid_parent_prob,
        mut_types_next ==1))[0]
    deme_arr_next[right_parent_idxs] = (deme_arr_after_recom[right_parent_idxs] + 1).astype(np.int64)

    left_edge_idxs = np.where(deme_arr_next < 0)[0]
    deme_arr_next[left_edge_idxs] = (np.zeros_like(left_edge_idxs)).astype(np.int64)

    right_edge_idxs = np.where(deme_arr_next > L - 1)[0]
    deme_arr_next[right_edge_idxs] = (np.ones_like(right_edge_idxs) *
                                      (L - 1)).astype(np.int64)

    mut_idxs = np.concatenate((left_parent_idxs, mid_parent_idxs, right_parent_idxs))
    ind_in_deme_arr_next[mut_idxs] = (np.floor(
        choice_rand[mut_idxs] * np.take(rho_e_parent,
                              deme_arr_next[mut_idxs]))).astype(np.int64)

    # then same for wt
    left_parent_prob = m / 2 * np.take(rho_wt_parent_extended, deme_arr)
    mid_parent_prob = (1 - m) * np.take(rho_wt_parent_extended, deme_arr + 1)
    right_parent_prob = m / 2 * np.take(rho_wt_parent_extended, deme_arr + 2)
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
    deme_arr_next[left_parent_idxs] = (deme_arr_after_recom[left_parent_idxs] - 1).astype(np.int64)

    mid_parent_idxs = np.where(np.logical_and.reduce((
        which_parent_rand > left_parent_prob,
        which_parent_rand < mid_parent_prob,
        mut_types_next == 0)))[0]
    deme_arr_next[mid_parent_idxs] = (deme_arr_after_recom[mid_parent_idxs]).astype(np.int64)


    right_parent_idxs = np.where(np.logical_and(
        which_parent_rand > mid_parent_prob,
        mut_types_next == 0))[0]
    deme_arr_next[right_parent_idxs] = (deme_arr_after_recom[right_parent_idxs] + 1).astype(np.int64)

    left_edge_idxs = np.where(deme_arr_next < 0)[0]
    deme_arr_next[left_edge_idxs] = (np.zeros_like(left_edge_idxs)).astype(np.int64)

    right_edge_idxs = np.where(deme_arr_next > L - 1)[0]
    deme_arr_next[right_edge_idxs] = (np.ones_like(right_edge_idxs) * (L - 1)).astype(np.int64)

    wt_idxs = np.concatenate((left_parent_idxs, mid_parent_idxs, right_parent_idxs))
    ind_in_deme_arr_next[wt_idxs] = (np.floor(
        choice_rand[wt_idxs] * np.take(rho_wt_parent,
                              deme_arr_next[wt_idxs]))).astype(np.int64)
    individuals2 = np.vstack((mut_types_next,
                              deme_arr_next,
                              ind_in_deme_arr_next)).T
    return individuals2

def runner(idx):
    rho_e_0 = lines[-1]
    n = nbase
    SFS = np.zeros(n)

    rho_e = rho_e_0
    rho_e = (rho_e).astype(np.int64)
    # Ne = [range(Ne_0[i]) for i in range(len(Ne_0))]
    individuals = [] # format will be [mut_type, deme index, individual index (inside the deme)]
    if n < round(sum(rho_e)):
        individuals_location = random.choice(np.arange(0
        , round(sum(rho_e))), size = n, replace = False)

    else:
        individuals_location = np.arange(0
        , round(sum(rho_e)))
        n = sum(rho_e)
    ind_inside_deme = np.mod(individuals_location, rho)
    deme_ind = (individuals_location - ind_inside_deme) // rho

    for k in range(n):
        # 1 is for having beneficial mutation
        individuals.append([1, deme_ind[k], ind_inside_deme[k]])
        
    individuals = np.array(individuals)
    unique, leaf_counts = np.unique(individuals, axis = 0, return_counts = True)
    hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))


    t_after_fix = 0
    while t_after_fix < T_after_fix:
        t_after_fix += 1
        individuals2 = get_individuals2_new(rho_e, rho_e, individuals)
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
        rho_e_parent = (lines[line_num]).astype(np.int64)
        individuals2 = get_individuals2_new(rho_e, rho_e_parent, individuals)
        individuals2 = np.repeat(individuals2, leaf_counts, axis = 0)
        unique, leaf_counts = np.unique(individuals2, axis = 0,
                                                        return_counts = True)
        hist, bin_edges = np.histogram(leaf_counts,
                                       bins = np.arange(1, n + 2))

        SFS += hist
        individuals = unique
        rho_e = rho_e_parent

    # The first round of coalescence simulation ends when we get to the time
    # when the beneficial mmutation arose, and all the left over individuals
    # are WT. From this point on, the coalescence will be extremely slow.
    # Therefore, we will run the same kind of simulation until the individuals
    # disperse for L^2 / m / N. We speed up the simulation by recording the number
    # of generations between the coalescence events.
    left_individuals = len(individuals)
    branch_len = 0 # number of generations until first merging event
    extra_gen = 0 # extra run time before stopping the coalescent simulation.

    while left_individuals > 1 and extra_gen < int(L ** 2 / m / rho):

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
        T2 = random.exponential(
                rho * L / ((left_individuals - 1) * left_individuals / 2))
        
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




    if np.mod(idx, 1000) == 0:
        np.savetxt('expected_SFS_L={}_rho={}_s={:.2e}_m={:.2e}_r={:.2e}_tfinal={}_nsample={}_tfix={}_sample_uniform_n={}_{}.txt'.format(L,
           rho, s, m, r, tfinal, n, T_after_fix, n_forward, idx), SFS)
    return SFS

if __name__ == '__main__':

        # this is the true sample number in case Ne < nbase.
    # print(individuals)
    p = Pool(25)
    
    SFS = np.sum(p.map(runner, range(N_SFS)), axis = 0) / N_SFS
    np.savetxt('expected_SFS_L={}_rho={}_s={:.2e}_m={:.2e}_r={:.2e}_tfinal={}_nsample={}_tfix={}_sample_uniform_navg={}_{}.txt'.format(L,
                rho, s, m, r, tfinal, nbase, T_after_fix, N_SFS, n_forward), SFS)
    
    
    # SFS_sum = runner(0)
    # for i in range(N_SFS - 1):
    #     SFS = runner(0)
    #     SFS_sum += SFS
    # SFS_avg = SFS_sum / N_SFS
    # np.savetxt('SFS_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_Un={:.6f}_sample_middle_navg={}.txt'.format(L,
    #             N, s, m, r, tfinal, nbase, Un, N_SFS), SFS)
