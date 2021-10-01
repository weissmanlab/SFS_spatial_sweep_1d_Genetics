#!/usr/bin/env python3
'''
Coalescent simulation with recombination.
open the frequency file read backward in time until n individuals all coalesce.
As we go, we will record the number of leaves for each node. 
If the neutral mutation arises as a Poisson process with rate Un = 1, 
the expected AFS is the histogram of the leaf counts 

Update : pre-sweep coalescent process changed to account for more than one 
mergers in a single generation when there are many individuals are left. 
(compared to sqrt(2 * N))

Update : Convert the 2-D array into a 1-D array with migration 'patterns' calculated using a single index. Details in the function. 

'''
import numpy as np
from multiprocessing import Pool
from multiprocessing import Process
import sys
from numpy import random
from functions_combined import *
import matplotlib.pyplot as plt
import time

start = time.time()

L = 100
rho = 500
s = 0.050
m = 0.25
l0 = 50
Nforw = 1
nbase = 30000
N_SFS = 5
T_after_fix = 0
dimension = 2
extra_gen_cutoff = int(L ** 2 / m / rho / 100)
r = 0



# L = int(sys.argv[1]) # number of demes
# rho = int(sys.argv[2]) # deme capacity
# s = format(float(sys.argv[3]),'0.3f') # selection coef
# m = float(sys.argv[4]) # migration rate
# m_file = format(m,'0.2f')
# l0 = int(sys.argv[5])##l0 for the file name string
# nbase = int(sys.argv[6]) # sample size
# N_SFS = int(sys.argv[7]) # number of coalescent simulation we run for ensemble average.
# T_after_fix = int(sys.argv[8]) # number of generations between fixation and sampling
# Nforw = int(sys.argv[9]) ##Forwrad simulation Number
# dimension = int(sys.argv[10]) ##To check if to run in 1-D or 2-D 
# extra_gen_cutoff = sys.argv[11] if len(sys.argv) >= 12 else int(L ** 2 / m / rho)
# r = float(sys.argv[12]) # recombination rate

## python backwards_combined.py L N s m l0 n_base N_SFS T_after_fix Nforw dimensions
## python backwards_combined.py 500 20000 0.05 0.250 1 10000 4 3000 1 1


if (dimension == 1):
    print ('1_D')
    fname = 'L={}_N={}_s={:.6f}_m={:.6f}_tfinal=1000000_{}.txt'.format(L, rho, s, m, Nforw)
    print(fname)
    lines = np.loadtxt(fname, dtype = np.int64)

elif(dimension == 2):
    print('2-D')
    ##Read input file and flatten lines[i] into 1-D array here for use within existing framework
    fname = 'D:\SFS_spatial_sweep\Combined_codes\L={}_N={}_s={:.3f}_m={:.2f}_l0={}_Nforw={}.txt'.format(L, rho, s, m, l0, Nforw)
    # fname = 'L='+str(L)+'_N='+str(rho)+'_s='+str(s)+'_m='+str(m_file)+'_l0='+str(l0)+'_Nforw='+str(Nforw)+'.txt'
    #print(fname)
    lines = np.loadtxt(fname, dtype = np.int64)
    #print(len(lines))
    tfinal = int(len(lines) / L)
    lines = lines.reshape(tfinal, L, L)  ###lines = Data from forward simulation
    lines = [lines[i].flatten() for i in range(len(lines))]
    print(tfinal)    
else:
    print("Dimension not recognised")
    quit()

def runner(idx):
    '''
    idx: argumemt to spawn multiple jobs. Not relevant for the body of the simulation. 
    
    Picks a generation, starting from the last one. Tracks the sampled indviduals backwards in time 
    till all mergers have happened and then calculates the branch length to get the SFS.
    '''
    rho_e = lines[-1]    ##Mutants in a deme. We pick the absolute last generation we simulated forward in time till to start simulating backwards. 
    n = nbase
    SFS = np.zeros(n)
    rho_e = (rho_e).astype(np.int64)

    '''Creating the primary data structure by sampling'''
    #individuals format will be [mut_type, deme index, individual index (inside the deme)] in a row and as many rows as individuals we sample
    individuals = sample_data(rho_e, n, rho)
    unique, leaf_counts = np.unique(individuals, axis = 0, return_counts = True) #To avoid overcounting of repeated locations and ensure we actually sample what we want
    ###The AFS is the histogram of these leaf counts counted over all generations
    hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))


    '''Starting simulations backward in time now'''    
 
    t_after_fix = 0    
    '''Going backwards in time for user given generations after sweep has fixed'''
    #Since the population doesnt change once sweep has fixed, we dont need to change the input data for this
    
    while t_after_fix < T_after_fix:
        #print ('sampling after fixation')
        t_after_fix += 1
        individuals_post_recombination = recombination(rho_e, rho_e, individuals, rho, r)
        individuals_post_migration = migration(rho_e, individuals_post_recombination, rho, m, L, dimension)

        individuals, leaf_counts = coalescent(individuals_post_migration, leaf_counts)
        hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))
        SFS += hist

    '''Going back in time during the sweep, all the way back to the time when mutation is first seeded'''
    line_num = -1
    
    while (len(individuals) > 1) and (line_num > -len(lines)):  ##The maximum we can go backward in time is till T2
        #print('sweep')

        line_num -= 1
        rho_e_parent = (lines[line_num]).astype(np.int64) ##Getting the parent generation from each time step of our forward simulation
        individuals_post_recombination = recombination(rho_e, rho_e_parent, individuals, rho, r)
        individuals_post_migration = migration(rho_e_parent, individuals_post_recombination, rho, m, L, dimension)
        individuals, leaf_counts = coalescent(individuals_post_migration, leaf_counts)
        hist, bin_edges = np.histogram(leaf_counts,bins = np.arange(1, n + 2))

        SFS += hist
        rho_e = rho_e_parent
    

    '''
    The first round of coalescence simulation ends when we get to the time
    when the beneficial mmutation arose, and all the left over individuals
    are WT. 
    If there is no recombination, the simulation should end here since the MRCA 
    will be the benficial mutation that gave rise to the sweep. 
    
    From this point on, the coalescence will be extremely slow.
    Therefore, we will run the same kind of simulation until the individuals
    disperse for L^2 / m / N. We speed up the simulation by recording the number
    of generations between the coalescence events.
    '''
    

    left_individuals = len(individuals)  ##Individuas left to still coalesce in time before the mutation arose
    branch_len = 0 # number of generations until first merging event
    extra_gen = 0 # extra run time before stopping the coalescent simulation.


    while left_individuals > 1 and extra_gen < extra_gen_cutoff:  ###Will it be same for 2-D
        #print('extra gens')
        branch_len += 1
        extra_gen += 1
        rho_e_parent_pre_sweep = np.zeros_like(rho_e)
        rho_e_parent_pre_sweep = (rho_e_parent_pre_sweep).astype(np.int64)

        # Before the beneficial mutation appeared, everyone is WT. Thus, rho_e_parent is just an array of zeros.
        # Also, we can skip the recombination step, since the genotype will always be zero for everyone.
        individuals_post_migration = migration(rho_e_parent_pre_sweep, individuals, rho, m, L, dimension)
        individuals, leaf_counts = coalescent(individuals_post_migration, leaf_counts)
        current_individuals_counts = len(individuals)
        if current_individuals_counts < left_individuals:
            hist, bin_edges = np.histogram(leaf_counts * branch_len, bins = np.arange(1, n + 2))
            SFS += hist
            branch_len = 0
        left_individuals = current_individuals_counts
        
        
    '''    
    After the individuals disperse, we may ignore the spatial structure.
    If there are more than sqrt(2 * N * L) individuals left, it is likely
    to have more than 1 merging event in a single generation. Thus, we manually
    draw parent for every individual and find if some of them are the same.
    Once there are much smaller number of individuals left, it will take a long
    time to coalesce. Thus, we draw T2 (time until the first coalescent) from
    geometric prob. distribution. (This is found from prob. of choosing different
    parent for every individual left.)
    '''
    
    while left_individuals > 1:
        #print('approximation')
        T2 = 1 + random.exponential(2 * rho * L / ((left_individuals - 1) * left_individuals / 2)) ###Drawing the T2 from geometric distribution
        hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))
        SFS += hist * T2 ##That many branches will exist
        # Now choose two random branches that will coalesce in T2 generations. 
        coal_inds = random.choice(range(left_individuals), size = 2, replace = False)
        # Like in the coalescent function, as two branches merge, two leaf counts are summed up and the number of remaining individuals decreases by 1.
        leaf_counts_copy = leaf_counts
        leaf_counts_copy[coal_inds[0]] += leaf_counts[coal_inds[1]]
        leaf_counts_copy[coal_inds[1]] = 0
        leaf_counts = leaf_counts_copy[leaf_counts_copy > 0]
        left_individuals -= 1


    print ('caclulations')
    f = np.arange(1, n + 1) / n  ##Calculating frquency of mutant
    H = np.sum(2 * f * (1 - f) * SFS) / np.sum(SFS)  ##Calculating Heterezygosity
    return SFS, H




if __name__ == '__main__':
    p = Pool(1)    ##Nunmber of cores you want to spawn jobs on
    ret = p.map(runner, range(N_SFS))  ###Number of simulations you want to average over    
    SFS_items = [r[0] for r in ret]
    H_items = [r[1] for r in ret]
    SFS = np.sum(SFS_items, axis=0)
    SFS /= N_SFS
    np.savetxt('expected_SFS_L={}_N={}_s={:.3f}_m={:.2f}_nsample={}_tfix={}_sample_uniform_navg={}_Nforw={}.txt'.format(int(L),
                int(rho), float(s), float(m), int(nbase), int(T_after_fix), int(N_SFS), int(Nforw)), SFS)


    test = np.arange(1, nbase)
    n = len(SFS)
    f = np.arange(1, n + 1) / n
    '''Plotting SFS to see in log scale'''
    plt.xlabel('frequency')
    plt.ylabel('Number of alleles')
    plt.plot(f, SFS)
    plt.yscale('log')
    plt.xscale('log')
    #plt.savefig('Test.jpeg')
    end = time.time()
    print(start-end) 
    #print(psutil.virtual_memory())
    plt.show()
    
