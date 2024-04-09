"""
    __init__.py: a python program to integrate the resource dynamics on a grid

"""

import numpy as np

from R_space_dynamics import *
from N_space_dynamics import *
from visualization  import *

# initialize R0
n_r = 2
n_s = 2
n   = 40
R0  = np.ones((n, n, n_r))
R0[:,:20,0]=10
g   = np.array([0.5,0.5]) # growth convertion factors

# define parameters
param = {
    # model parameters
    'R0' : R0,                                         # initial conc. nxnxn_r [monod constants]
    'w'  : np.ones((n_r)),                             # energy conversion     [energy/mass]
    'l'  : np.ones((n_r)),                             # leakage               [adim]
    'tau': 1,                                          # reinsertion rate inv. [time] 
    'g'  : g,                                          # growth conv. factors  [1/energy]
    'k'  : 1,                                          # monod constant        [mass/vol]
    # sor algorithm parameters
    'n'  : n,                                          # grid points in each dim
    'sor': 1.50,                                       # relaxation parameter
    'L0' : 40,                                         # grid true size        [length]
    'D'  : 1e3                                         # diffusion constant    [area/time]   
}

# initialize species 
N = np.zeros((n, n, n_s))
# random species disposition
#for i in range(n):
#    for j in range(n):
#        k = np.random.randint(0, n_s)  
#        N[i, j, k] = 1
N[:,:20,0]=1
N[:,20:,1]=1

# make matrices
up_mat  = np.array([[0.8,0.],[0.,1]])
met_mat = np.array([[0.,0.],[1.,0.]])
print(up_mat)
print(met_mat)

mat = {
    'uptake' : up_mat,
    'met'    : met_mat
}

# plot initial states
#R_ongrid(R0,param)
#N_ongrid(N)
makenet(met_mat)
vispreferences(up_mat)

# run SOR algorithm
R_eq = SOR_multigrid(n,0,param['R0'],N,param,mat,0.001)
R_ongrid(R_eq)

G = growth_rates(R_eq,N,param,mat)
print(np.min(G),np.max(G))
G_ongrid(G)

N_dec = decode(N)
N_dec = death_birth(N,G)
N_enc = encode(N_dec)
N_ongrid(N)

