"""
update.py: code containing the main update function for the species matrix

CONTAINS:
    - update: function using the equilibrium concentrations of resources to calculate growth rates
              and perform the updates
    - run: running function

"""

from R_space_dynamics import *
from N_space_dynamics import *
from time import time

import numpy as np

from visualization import R_ongrid

#-----------------------------------------------------------------------------------------
# define update(R_curr,N_curr,param,mat) function to perform one step of update of the grid
# first solves for equilibrium R_eq, then calculates growth rates gr, and finally updates N

def update(R_curr,N_curr,param,mat):
    """
    R_curr: current state of resources vectors matrix
    N_curr: current state of species matrix
    param: parameters dictionary
    mat: matrices dictionary

    returns the updated species matrix (both encoded and decoded), the equilibrium concentrations

    """
    R_eq = SOR_multigrid(N_curr.shape[0],param['ref'],R_curr,N_curr,param,mat)     
    gr   = growth_rates(R_eq,N_curr,param,mat)
    N_curr_dec = decode(N_curr)
    N_new_dec  = death_birth(N_curr_dec,gr)
    N_new      = encode(N_new_dec)

    return R_eq,N_new,N_new_dec

#-----------------------------------------------------------------------------------------
# define run(times) function to run the simulation a desired number of times

def run(times,R0,N0,param,mat):
    """
    times: number of steps we want to run the simulation for
    R0: initial state of resources vectors matrix
    N0: initial state of species matrix
    param: parameters dictionary
    mat: matrices dictionary

    returns steps to produce a movie of the simulation and final R and N states

    """

    # start time
    t = time()


    # first time run on multigrid
    R_eq,N_new,N_new_dec = update(R0,N0,param,mat)
    # empty steps list
    steps = [N_new_dec]

    # then always solve on finest grid
    for i in range(times):
        print('iteration number: ',i)

        # adapt to finest grid
        r = param['ref']+1
        R = change_grid(R_eq,N_new.shape[0]*r)
        N = np.zeros((N_new.shape[0]*r,N_new.shape[0]*r,N_new.shape[2]))              # new grid for individuals

        # inverse restriction of individuals grid
        for k in range(N_new.shape[0]):
            for j in range(N_new.shape[0]):
                N[r*k:r*k+r,r*j:r*j+r,:]=N_new[k,j,:] 
        print('current grid size: ', N_new.shape[0]*r)

        # run on finest grid
        R_eq = SOR(R,N,param,mat)
        R_eq = change_grid(R,N0.shape[0])
        N    = change_grid(N,N0.shape[0])     
        gr   = growth_rates(R_eq,N,param,mat)
        N_dec = decode(N) 
        N_new_dec  = death_birth(N_dec,gr)
        N_new      = encode(N_new_dec)
        
        steps.append(N_new_dec)

    tf = time()
    print("\nTotal running time: %.2f min" % ((tf-t)/60))

    return steps, R_eq, N_new_dec


