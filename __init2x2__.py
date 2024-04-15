"""
__init2x2__.py: initialization file to run simulation of a community with two species and two resources
                one species consumes the externally supplied resource and produces a metabolite that is
                consumed by the second species: c->A->x->B. The grid is 40x40 and it is initialized to 
                have species A on one side and species B on the other; resources initial conditions are 
                saturation of resource c everywhere and absence of resource x everywhere.

                accuracy 10-6: 
                TEMPO DI ESECUZIONE SIMULAZIONE (D=10^5): 8.91 minuti
                TEMPO DI ESECUZIONE SIMULAZIONE (D=10^3): 38.10 minuti
                TEMPO PRODUZIONE VIDEO: about 5 minutes
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
# saturate resource c everywhere
R0[:,:,0]=10
g   = np.array([0.5,0.5]) 
m   = np.array([0.,0.])

# initialize species grid: A|B
N = np.zeros((n,n,n_s))
N[:,:20,0]=1
N[:,20:,1]=1

# define parameters
param = {
    # model parameters
    'R0' : R0.copy(),                                  # initial conc. nxnxn_r [monod constants]
    'w'  : np.ones((n_r))*20,                          # energy conversion     [energy/mass]
    'l'  : np.ones((n_r)),                             # leakage               [adim]
    'tau': np.array([1,1000000]),                      # reinsertion rate inv. [time], basically no replenishment for byprod and no dil.
    'g'  : g,                                          # growth conv. factors  [1/energy]
    'm'  : m,                                          # maintainance requ.    [energy/time]
    
    # sor algorithm parameters
    'n'  : n,                                          # grid points in each dim
    'sor': 1.80,                                       # relaxation parameter
    'L0' : 40,                                         # grid true size        [length]
    'D'  : 1e5,                                        # diffusion constant    [area/time] 
    'acc': 1e-6,                                       # maximum accepted stopping criterion   
    'ref': 1                                           # number of grid refinements to perform 
}

# make matrices
up_mat   = np.array([[0.8,0.],[0.,0.8]])
met_mat  = np.array([[0.,0.],[1.,0.]])
sign_mat = np.ones((n_s, n_r))
mat_ess  = np.array([[0.,0.],[0.,0.]])
print(up_mat)
print(met_mat)
print(sign_mat)
print(mat_ess)

mat = {
    'uptake' : up_mat,
    'met'    : met_mat,
    'sign'   : sign_mat,
    'ess'    : mat_ess 
}

# visualize matrices
vispreferences(mat)
makenet(met_mat)

#-------------------------------------------------------------------------------------------------------------
# SIMULATION

# run 1000 steps 
steps,R_fin,N_fin = run(1000,R0,N,param,mat)

# plot final R and N grids
R_ongrid(R_fin)
N_ongrid(encode(N_fin))

# produce a movie of the simulation
def next(frame):
    plt.clf()
    plt.imshow(steps[frame], cmap='cool')  
    plt.title(f'Timestep {frame}')
    plt.axis('off')

fig = plt.figure()
ani = FuncAnimation(fig, next, frames=len(steps), interval=5)

writer = PillowWriter(fps=1000/5)  
ani.save('animation.gif', writer=writer)




