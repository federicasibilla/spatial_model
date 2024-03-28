"""
    __init__.py: a python program to integrate the resource dynamics on a grid

"""

import numpy as np
import matplotlib.pyplot as plt

from time import time
from mpl_toolkits import mplot3d
from space_dynamics import *


# initialize R0
n_r = 3
n_s = 2
n   = 100
R0  = np.random.uniform(8, 12, size=(n, n, n_r))

# define parameters
param = {
    # model parameters
    'R0' : R0,                                         # initial conc. nxnxn_r [mass/vol]
    'w'  : np.random.uniform(0, 1, size=n_r),          # energy conversion     [energy/mass]
    'l'  : np.random.uniform(0, 1, size=n_r),          # leakage               [adim]
    'tau': 10,                                         # reinsertion rate inv. [time]
    # sor algorithm parameters
    'n'  : n,                                          # grid points in each dim
    'sor': 1.38,                                       # relaxation parameter
    'L0' : 10                                          # grid true size        [length]
}

# initialize species 
N = np.zeros((n, n, n_s))
# random species disposition
for i in range(n):
    for j in range(n):
        k = np.random.randint(0, n_s)  
        N[i, j, k] = 1

# make matrices
up_mat  = np.random.uniform(0, 1, size=(n_s,n_r))
met_mat = np.random.uniform(0, 1, size=(n_r,n_r))
print(up_mat,met_mat)

mat = {
    'uptake' : up_mat,
    'met'    : met_mat
}

# plot initial states
#R_ongrid(R0,param)
#R_ongrid(N,param)

# run SOR algorithm
R_eq = SOR(N,param,mat,1e-2)
R_ongrid(R_eq,param)


