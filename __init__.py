"""
    __init__.py: a python program to integrate the resource dynamics on a grid

"""

import numpy as np
import matplotlib.pyplot as plt

from time import time
from mpl_toolkits import mplot3d
from space_dynamics import *
from adhoc import *


# initialize R0
n_r = 10
n_s = 3
n   = 101
R0  = np.random.uniform(80, 120, size=(n, n, n_r))
#R0=R0_adhoc()

# define parameters
param = {
    # model parameters
    'R0' : R0,                                         # initial conc. nxnxn_r [mass/vol]
    'w'  : np.random.uniform(0, 1, size=n_r),          # energy conversion     [energy/mass]
    'l'  : np.random.uniform(0, 1, size=n_r),          # leakage               [adim]
    'tau': 10,                                         # reinsertion rate inv. [time]
    # sor algorithm parameters
    'n'  : n,                                          # grid points in each dim
    'sor': 1.55,                                       # relaxation parameter
    'L0' : 5                                           # grid true size        [length]
}

# initialize species 
#N = N0_adhoc()
N = np.random.uniform(80, 120, size=(n, n, n_s))

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
R_eq = SOR(N,param,mat,1e-1)
R_ongrid(R_eq,param)


