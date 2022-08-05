'''
List of functions being used to create the backward time simulation for sweeps
Tracks migration and parent lineage as they coalesce back in time 
'''

import numpy as np
from numpy import random


def safe_divide(a, b, val = 0):
    return np.divide(a, b, out = np.full_like(a, val), where = b != 0)
    

def recombination(rho_e, individuals, rho, r):
    '''
    Let individuals recombine with others in the same deme. This will change the mutation types of the individuals.
    Output - individuals with the updated genotypes after recombination events.
    '''
    mut_types, deme_arr, ind_in_deme_arr = individuals.T
    mut_types = (mut_types).astype(np.int64)
    deme_arr = (deme_arr).astype(np.int64)


    post_rec_mut_types = mut_types
    # Draw a random number between 0 and 1 for every individual. 
    # Those numbers will be compared to the CDF of recombination with WT or a mutant.
    rand_for_recom = random.random(len(mut_types))


    corresponding_rho_e = rho_e[deme_arr]
    # Three possible events - recombine with a mutant, recombine with WT, no recombination
    # CDF of first event = r * rho_e / rho
    # CDF of second event = r * rho_e / rho + r * rho_wt / rho = r
    # CDF of no recombination = 1

    # a list of indices of the individuals that recombine with a mutant 
    idxs_recom_mut = np.where(rand_for_recom < r * corresponding_rho_e / rho)[0]

    # a list of indices of the individuals that recombine with a WT 
    idxs_recom_wt = np.where(np.logical_and(
        rand_for_recom > r * corresponding_rho_e / rho, 
        rand_for_recom < r
    ))[0]

    post_rec_mut_types[idxs_recom_mut] = 1
    post_rec_mut_types[idxs_recom_wt] = 0
    post_rec_mut_types = (post_rec_mut_types).astype(np.int64)
    individuals_post_rec = np.vstack((post_rec_mut_types, deme_arr, ind_in_deme_arr)).T

    del post_rec_mut_types
    del deme_arr
    del mut_types
    del ind_in_deme_arr
    return individuals_post_rec



def migration(rho_e_parent, individuals, rho , m, L, dimension):
    '''
    Find parent for each individual after the recombination step. 
    We consider the weighted probability for finding a parent in any of the neighboring deme, from the
    number of individuals of the same genotype in each deme and migration rate m. 

    Input :
    rho_e = array of number of mutant individuals in each deme
    individuals: array of individuals. Each individual has three numbers as their identifier
    1) mut_type = 0 if WT, 1 if a mutant
    2) deme_arr = index for the deme. Ranges from 0 to L - 1 or (L * L) - 1 depending on the dimension
    3) ind_in_deme_arr = index within the deme. Ranges from 0 to rho_e - 1 for a mutant, 0 to rho_WT - 1 for a WT
    rho = number of individuals in a deme ( = rho_e + rho_WT)
    m = migration rate
    L = number of demes in one direction
    dimension = 1 or 2

    Output : individuals after the migration step. Same format as the input individuals. 
    '''   
    mut_types, deme_arr, ind_in_deme_arr = individuals.T
    deme_arr = (deme_arr).astype(np.int64)
    mut_types_next = np.copy(mut_types)
    deme_arr_next = np.zeros_like(deme_arr)
    ind_in_deme_arr_next = np.zeros_like(ind_in_deme_arr)
    rho_wt_parent = (rho - rho_e_parent).astype(np.int64)
    len_inds = len(deme_arr)
    which_parent_rand = random.random(len_inds) # to choose btw left/mid/right/top/bottom

    
    if(dimension == 1):  

        rho_e_parent_extended = np.concatenate(([rho_e_parent[0]], rho_e_parent, [rho_e_parent[-1]]))
        rho_wt_parent_extended = np.concatenate(([rho_wt_parent[0]], rho_wt_parent, [rho_wt_parent[-1]]))
        
        #### Start with the mutant type #######        
        
        # Find the parent's deme and then its index inside the deme
        left_parent_prob_mut = m / 2 * np.take(rho_e_parent_extended, deme_arr)
        mid_parent_prob_mut = (1 - m) * np.take(rho_e_parent_extended, deme_arr + 1)
        right_parent_prob_mut = m / 2 * np.take(rho_e_parent_extended, deme_arr + 2)
        total_prob_mut = (left_parent_prob_mut + mid_parent_prob_mut + right_parent_prob_mut)

        # Set the cumulative probability
        mid_parent_prob_mut_cumulative = np.divide(left_parent_prob_mut + mid_parent_prob_mut,
        total_prob_mut)
        left_parent_prob_mut_cumulative = np.divide(left_parent_prob_mut, 
                                                      total_prob_mut)
         #Right parent probablity is determined by the two numbers above. 
       
        # The location of individuals where probablty to move left is found is in left_parent_idxs_mut/wt depending on type.
        # For those indivduals, move the location. Repeat for others


        left_parent_idxs_mut = np.where(np.logical_and(
            which_parent_rand < left_parent_prob_mut_cumulative,
            mut_types_next == 1))[0]
        deme_arr_next[left_parent_idxs_mut] = (deme_arr[left_parent_idxs_mut] 
                                               - 1).astype(np.int64)

        mid_parent_idxs_mut = np.where(np.logical_and.reduce((
            which_parent_rand > left_parent_prob_mut_cumulative,
            which_parent_rand < mid_parent_prob_mut_cumulative,
            mut_types_next == 1)))[0]
        deme_arr_next[mid_parent_idxs_mut] = (deme_arr[mid_parent_idxs_mut]).astype(np.int64)


        right_parent_idxs_mut = np.where(np.logical_and(
            which_parent_rand > mid_parent_prob_mut_cumulative,
            mut_types_next == 1))[0]
        deme_arr_next[right_parent_idxs_mut] = (deme_arr[right_parent_idxs_mut]
                                                + 1).astype(np.int64)

        left_edge_idxs_mut = np.where(np.logical_and(deme_arr_next < 0, 
                                                 mut_types_next == 1))[0]
        deme_arr_next[left_edge_idxs_mut] = (np.zeros_like(left_edge_idxs_mut)
                                             ).astype(np.int64)

        right_edge_idxs_mut = np.where(np.logical_and(deme_arr_next > L - 1, 
                                                  mut_types_next == 1))[0]
        deme_arr_next[right_edge_idxs_mut] = (np.ones_like(right_edge_idxs_mut) *
                                      (L - 1)).astype(np.int64)

        mut_idxs = np.concatenate((left_parent_idxs_mut, 
                               mid_parent_idxs_mut, 
                               right_parent_idxs_mut))
        ind_in_deme_arr_next[mut_idxs] = random.randint(0, np.take(rho_e_parent, 
                        deme_arr_next[mut_idxs]))
        
        #### Now repeat the same thing for the wild-type (WT) population
        left_parent_prob_wt = m / 2 * np.take(rho_wt_parent_extended, deme_arr)
        mid_parent_prob_wt = (1 - m) * np.take(rho_wt_parent_extended, deme_arr + 1)
        right_parent_prob_wt = m / 2 * np.take(rho_wt_parent_extended, deme_arr + 2)
        total_prob_wt = (left_parent_prob_wt + mid_parent_prob_wt + right_parent_prob_wt)

        # Set the cumulative probability
        mid_parent_prob_wt_cumulative = safe_divide(left_parent_prob_wt + mid_parent_prob_wt,
        total_prob_wt, val = 1)
        left_parent_prob_wt_cumulative = safe_divide(left_parent_prob_wt, 
                                                     total_prob_wt, val = 1)


        left_parent_idxs_wt = np.where(np.logical_and(
            which_parent_rand < left_parent_prob_wt_cumulative,
            mut_types_next == 0))[0]
        deme_arr_next[left_parent_idxs_wt] = (deme_arr[left_parent_idxs_wt] 
                                               - 1).astype(np.int64)

        mid_parent_idxs_wt = np.where(np.logical_and.reduce((
            which_parent_rand > left_parent_prob_wt_cumulative,
            which_parent_rand < mid_parent_prob_wt_cumulative,
            mut_types_next == 0)))[0]
        deme_arr_next[mid_parent_idxs_wt] = (deme_arr[mid_parent_idxs_wt]).astype(np.int64)


        right_parent_idxs_wt = np.where(np.logical_and(
            which_parent_rand > mid_parent_prob_wt_cumulative,
            mut_types_next == 0))[0]
        deme_arr_next[right_parent_idxs_wt] = (deme_arr[right_parent_idxs_wt]
                                                + 1).astype(np.int64)

        left_edge_idxs_wt = np.where(np.logical_and(deme_arr_next < 0, 
                                                 mut_types_next == 0))[0]
        deme_arr_next[left_edge_idxs_wt] = (np.zeros_like(left_edge_idxs_wt)
                                             ).astype(np.int64)

        right_edge_idxs_wt = np.where(np.logical_and(deme_arr_next > L - 1, 
                                                  mut_types_next == 0))[0]
        deme_arr_next[right_edge_idxs_wt] = (np.ones_like(right_edge_idxs_wt) *
                                      (L - 1)).astype(np.int64)

        wt_idxs = np.concatenate((left_parent_idxs_wt, 
                               mid_parent_idxs_wt, 
                               right_parent_idxs_wt))
        ind_in_deme_arr_next[wt_idxs] = random.randint(0, np.take(rho_wt_parent, 
                        deme_arr_next[wt_idxs]))
        individuals_post_migration = np.vstack((mut_types_next, deme_arr_next, ind_in_deme_arr_next)).T


    else: ##2-D migrations 
        '''
        Migration patterns from a given index i in 2-D:
        Top = (i-L)%(L*L)
        Bottom = (i+L)%(L*)
        Right = i + [-(L-1) if (i+1)%L ==0 else 1][0]
        Left = i + [(L-1) if (i)%L ==0 else -1][0]
    
        This takes care of edge cases also
        '''

        '''Creating the probablities of moving to the nearest neigbor if it is of same type. The probablity is weighted by the number of mutants/wiltypes in the neighborinf deme'''
        left_parent_prob_mut = m / 4 * np.take(rho_e_parent,
        [(deme_arr[i] + (L - 1)) if (deme_arr[i]) % L == 0 
        else (deme_arr[i] - 1) for i in range(len(deme_arr))])
        right_parent_prob_mut = m / 4 *np.take(rho_e_parent, 
        [(deme_arr[i] - (L - 1)) if (deme_arr[i] + 1) % L == 0 
        else (deme_arr[i] + 1) for i in range(len(deme_arr))])    
        top_parent_prob_mut = m / 4 * np.take(rho_e_parent, ((deme_arr - L) % (L * L)))
        bottom_parent_prob_mut = m / 4 * np.take(rho_e_parent, ((deme_arr + L) % (L * L)))
        mid_parent_prob_mut = (1 - m) * np.take(rho_e_parent, deme_arr)
        total_prob_mut = (left_parent_prob_mut + right_parent_prob_mut 
        + top_parent_prob_mut + bottom_parent_prob_mut + mid_parent_prob_mut)
    
        '''Set the cumulative probability and normalise'''
        left_parent_prob_mut_cumulative = safe_divide(left_parent_prob_mut, total_prob_mut)
        right_parent_prob_mut_cumulative = safe_divide(left_parent_prob_mut + right_parent_prob_mut, total_prob_mut)
        top_parent_prob_mut_cumulative = safe_divide(left_parent_prob_mut 
        + right_parent_prob_mut + top_parent_prob_mut, total_prob_mut)
        bottom_parent_prob_mut_cumulative = safe_divide(left_parent_prob_mut + right_parent_prob_mut
         + top_parent_prob_mut + bottom_parent_prob_mut, total_prob_mut)

        '''The location of individuals where probablty to move left is found is in left_idx. For those indivduals, move the location. Repeat for others'''
        left_parent_idxs_mut = np.where(np.logical_and(
            which_parent_rand < left_parent_prob_mut_cumulative,
            mut_types_next == 1))[0]
        deme_arr_next[left_parent_idxs_mut] = ([(deme_arr[i] + L - 1) if deme_arr[i] % L == 0 
        else (deme_arr[i] - 1) for i in left_parent_idxs_mut])

        right_parent_idxs_mut = np.where(np.logical_and.reduce((
            which_parent_rand > left_parent_prob_mut_cumulative,
            which_parent_rand < right_parent_prob_mut_cumulative,
            mut_types_next == 1)))[0]
        deme_arr_next[right_parent_idxs_mut] = ([(deme_arr[i] - (L - 1)) if (deme_arr[i] + 1) % L == 0 
        else (deme_arr[i] + 1) for i in right_parent_idxs_mut])
    
        top_parent_idxs_mut = np.where(np.logical_and.reduce((
            which_parent_rand > right_parent_prob_mut_cumulative,
            which_parent_rand < top_parent_prob_mut_cumulative,
            mut_types_next == 1)))[0]
        deme_arr_next[top_parent_idxs_mut] = (((deme_arr[top_parent_idxs_mut] - L) % (L * L))).astype(np.int64)
        
        bottom_parent_idxs_mut = np.where(np.logical_and.reduce((
            which_parent_rand > top_parent_prob_mut_cumulative,
            which_parent_rand < bottom_parent_prob_mut_cumulative,
            mut_types_next == 1)))[0]
        deme_arr_next[bottom_parent_idxs_mut] = (((deme_arr[bottom_parent_idxs_mut] + L) % (L * L))).astype(np.int64)
    
        mid_parent_idxs_mut = np.where(np.logical_and(
            which_parent_rand > bottom_parent_prob_mut_cumulative,
            mut_types_next == 1))[0]
        deme_arr_next[mid_parent_idxs_mut] = (deme_arr[mid_parent_idxs_mut]).astype(np.int64)
    
        mut_idxs = np.concatenate((left_parent_idxs_mut, 
        right_parent_idxs_mut, top_parent_idxs_mut, bottom_parent_idxs_mut, mid_parent_idxs_mut))
        ind_in_deme_arr_next[mut_idxs] = random.randint(0, np.take(rho_e_parent, 
                        deme_arr_next[mut_idxs]))


        # repeat for the WT

        left_parent_prob_wt = m / 4 * np.take(rho_wt_parent,
        [(deme_arr[i] + (L - 1)) if (deme_arr[i]) % L == 0 
        else (deme_arr[i] - 1) for i in range(len(deme_arr))])
        right_parent_prob_wt = m / 4 *np.take(rho_wt_parent, 
        [(deme_arr[i] - (L - 1)) if (deme_arr[i] + 1) % L == 0 
        else (deme_arr[i] + 1) for i in range(len(deme_arr))])    
        top_parent_prob_wt = m / 4 * np.take(rho_wt_parent, ((deme_arr - L) % (L * L)))
        bottom_parent_prob_wt = m / 4 * np.take(rho_wt_parent, ((deme_arr + L) % (L * L)))
        mid_parent_prob_wt = (1 - m) * np.take(rho_wt_parent, deme_arr)
        total_prob_wt = (left_parent_prob_wt + right_parent_prob_wt 
        + top_parent_prob_wt + bottom_parent_prob_wt + mid_parent_prob_wt)
    
        '''Set the cumulative probability and normalise'''
        left_parent_prob_wt_cumulative = safe_divide(left_parent_prob_wt, total_prob_wt)
        right_parent_prob_wt_cumulative = safe_divide(left_parent_prob_wt + right_parent_prob_wt, total_prob_wt)
        top_parent_prob_wt_cumulative = safe_divide(left_parent_prob_wt 
        + right_parent_prob_wt + top_parent_prob_wt, total_prob_wt)
        bottom_parent_prob_wt_cumulative = safe_divide(left_parent_prob_wt + right_parent_prob_wt
         + top_parent_prob_wt + bottom_parent_prob_wt, total_prob_wt)

        left_parent_idxs_wt = np.where(np.logical_and(
            which_parent_rand < left_parent_prob_wt_cumulative,
            mut_types_next == 0))[0]
        deme_arr_next[left_parent_idxs_wt] = ([(deme_arr[i] + L - 1) if deme_arr[i] % L == 0 
        else (deme_arr[i] - 1) for i in left_parent_idxs_wt])

        right_parent_idxs_wt = np.where(np.logical_and.reduce((
            which_parent_rand > left_parent_prob_wt_cumulative,
            which_parent_rand < right_parent_prob_wt_cumulative,
            mut_types_next == 0)))[0]
        deme_arr_next[right_parent_idxs_wt] = ([(deme_arr[i] - (L - 1)) if (deme_arr[i] + 1) % L == 0 
        else (deme_arr[i] + 1) for i in right_parent_idxs_wt])
    
        top_parent_idxs_wt = np.where(np.logical_and.reduce((
            which_parent_rand > right_parent_prob_wt_cumulative,
            which_parent_rand < top_parent_prob_wt_cumulative,
            mut_types_next == 0)))[0]
        deme_arr_next[top_parent_idxs_wt] = (((deme_arr[top_parent_idxs_wt] - L) % (L * L))).astype(np.int64)
        
        bottom_parent_idxs_wt = np.where(np.logical_and.reduce((
            which_parent_rand > top_parent_prob_wt_cumulative,
            which_parent_rand < bottom_parent_prob_wt_cumulative,
            mut_types_next == 0)))[0]
        deme_arr_next[bottom_parent_idxs_wt] = (((deme_arr[bottom_parent_idxs_wt] + L) % (L * L))).astype(np.int64)
    
        mid_parent_idxs_wt = np.where(np.logical_and(
            which_parent_rand > bottom_parent_prob_wt_cumulative,
            mut_types_next == 0))[0]
        deme_arr_next[mid_parent_idxs_wt] = (deme_arr[mid_parent_idxs_wt]).astype(np.int64)
    
        wt_idxs = np.concatenate((left_parent_idxs_wt, 
        right_parent_idxs_wt, top_parent_idxs_wt, bottom_parent_idxs_wt, mid_parent_idxs_wt))
        ind_in_deme_arr_next[wt_idxs] = random.randint(0, np.take(rho_wt_parent, 
                        deme_arr_next[wt_idxs]))
        individuals_post_migration = np.vstack((mut_types_next, deme_arr_next, ind_in_deme_arr_next)).T

    del mut_types_next
    del deme_arr_next
    del ind_in_deme_arr_next
    del which_parent_rand
    del deme_arr
    del mut_types
    del ind_in_deme_arr
    return individuals_post_migration


def sample_data(rho_e, n, rho):
    '''
    Pick random deme locations for as many individuals as we want to sample and what the index within a deme is 
    Currently, uniform probablity distributions
    '''
    individuals = [] # format will be [mut_type, deme index, individual index (inside the deme)]
    if n < round(sum(rho_e)):   ##If sample size is less than the total number of mutants 
        individuals_location = random.choice(np.arange(0, round(sum(rho_e))), size = n, replace = False)

    else:   ##If sample size is equal or more than mutants
        individuals_location = np.arange(0, round(sum(rho_e)))
        n = sum(rho_e)
    
    ind_inside_deme = np.mod(individuals_location, rho)
    deme_ind = (individuals_location - ind_inside_deme) // rho

    ###Adding mutants to the data structure
    for k in range(n):  
        # 1 is for having beneficial mutation
        individuals.append([1., deme_ind[k], ind_inside_deme[k]])
        
    individuals = np.array(individuals)
    return individuals

def coalescent(parents, leaf_counts_offspring):
    '''
    Input - list of parents (there could be duplicates), leaf counts for the offsprings. 
    Each element of the two list should correspond to the same offspring individual. 

    Find whether there are any common parents going one generation backward in time. 
    (i.e. check coalescent events in 1 gen)
    Based on the information, update the list of leaf counts, which always sums up to nsample.
    If two branches coalesce, the length of the leaf counts list will be reduced by 1 because the 
    common parent's leaf counts = sum of two offspring's leaf counts.  
    
    Output - the set of parents (i.e. no duplicate), and leaf counts for the parents
    '''
    num_parents = len(parents)
    parents_set = np.array([])
    leaf_counts_parents = np.array([])
    for i in range(num_parents):
        parent = parents[i]
        leaf_count_offspring = leaf_counts_offspring[i]
        if len(parents_set) == 0:
            parents_set = np.array([parent])
            leaf_counts_parents = np.array([leaf_count_offspring])
        elif (parent ==  parents_set).all(1).any():
            idx = np.where((parents_set == parent).all(1))[0][0]
            leaf_counts_parents[idx] += leaf_count_offspring
        else:
            parents_set = np.append(parents_set, [parent], axis = 0)
            leaf_counts_parents = np.append(leaf_counts_parents, leaf_count_offspring)
    
    return parents_set, leaf_counts_parents
