"""
    __init__.py: a python program to integrate the resource dynamics on a grid

"""

import numpy as np

from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter


from R_space_dynamics import *
from N_space_dynamics import *
from visualization  import *
from update import *

# initialize R0
n_r = 2
n_s = 2
n   = 40
R0  = np.zeros((n, n, n_r))
R0[:,:,0]=10
g   = np.array([0.5,0.5]) # growth convertion factors

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
    'sor': 1.80,                                       # relaxation parameter
    'L0' : 40,                                         # grid true size        [length]
    'D'  : 1e3,                                        # diffusion constant    [area/time] 
    'acc': 1e-7,                                       # maximum accepted stopping criterion   
    'ref': 1                                           # number of grid refinements to perform 
}

# initialize species 
N = np.random.randint(2, size=(n, n, n_s))
#N = np.zeros((n,n,n_s))
# random species disposition
#for i in range(n):
#    for j in range(n):
#        k = np.random.randint(0, n_s)  
#        N[i, j, k] = 1
#N[:,:20,0]=1
#N[:,20:,1]=1

# make matrices
up_mat  = np.array([[0.8,0.],[0.,1]])
met_mat = np.array([[0.,0.],[1.,0.]])
print(up_mat)
print(met_mat)

mat = {
    'uptake' : up_mat,
    'met'    : met_mat
}

steps,R_fin,N_fin = run(3,R0,N,param,mat)

R_ongrid(R_fin)
N_ongrid(encode(N_fin))

def update(frame):
    plt.clf()
    plt.imshow(steps[frame], cmap='cool')  # Usa 'binary' per colori bianco e nero
    plt.title(f'Timestep {frame}')
    plt.axis('off')

# Crea l'animazione
fig = plt.figure()
ani = FuncAnimation(fig, update, frames=len(steps), interval=10)

#writer = PillowWriter(fps=1000/10)  # Specifica il numero di frame al secondo
#ani.save('animation.gif', writer=writer)

plt.show()


