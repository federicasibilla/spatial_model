
"""
__initCA__.py: code to run a simulation of random two species-four resources cellular automaton

"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation

from R_space_dynamics import *
from N_space_dynamics import *
from visualization  import *
from update import *

# initialize R0
n_r = 4
n_s = 2
n   = 40
R0  = np.ones((n, n, n_r))
R0[:,:,0]=10 # start with saturated supplied resource
g   = np.ones((n_s))*0.5                               # growth convertion factors

# initialize species 
N = np.zeros((n, n, n_s))
# random species disposition
for i in range(n):
    for j in range(n):
        k = np.random.randint(0, n_s)  
        N[i, j, k] = 1

# define parameters
param = {
    # model parameters
    'R0' : R0.copy(),                                  # initial conc. nxnxn_r [monod constants]
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

# make matrices
up_mat  = np.random.rand(n_s,n_r)
met_mat = np.random.rand(n_r,n_r)
met_mat[0,:] = 0 # fix supplied resource
print(up_mat)
print(met_mat)

mat = {
    'uptake' : up_mat,
    'met'    : met_mat
}

# list of arrays to store time steps 
steps = [decode(N)]

for i in range(1000):

    R_eq,N_new,N_new_dec = update(R0,N,param,mat)
    steps.append(N_new_dec)

    R0 = R_eq
    N  = N_new

R_ongrid(R_eq)

def update(frame):
    plt.clf()
    plt.imshow(steps[frame], cmap='cool')  # Usa 'binary' per colori bianco e nero
    plt.title(f'Timestep {frame}')
    plt.axis('off')

# Crea l'animazione
fig = plt.figure()
ani = FuncAnimation(fig, update, frames=len(steps), interval=50)

plt.show()
