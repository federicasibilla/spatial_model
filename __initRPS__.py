"""
__initRPS__.py: 

RUNNING TIME FOR 10^3 D AND 10^-6 STOP: 22.88 minutes (1000 steps)
                                        74.68 minutes (2000 steps)
            
"""
import numpy as np

from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

from R_space_dynamics import *
from N_space_dynamics import *
from visualization  import *
from update import *

# initialize R0
n_r = 6
n_s = 3
n   = 40

R0  = np.zeros((n, n, n_r))
# saturate resource c everywhere, zero toxins
R0[:,:,0:3]=10
g   = np.array([0.5,0.5,0.5]) 
m   = np.array([0.,0.,0.])

# initialize species grid: random
N = np.zeros((n,n,n_s))
for i in range(n):
    for j in range(n):
        idx = np.random.randint(3)
        N[i,j,idx]=1

# define parameters
param = {
    # model parameters
    'R0' : R0.copy(),                                   # initial conc. nxnxn_r [monod constants]
    'w'  : np.ones((n_r))*20,                           # energy conversion     [energy/mass]
    'l'  : np.array([1,1,1,0.,0.,0.]),                  # leakage               [adim]
    'tau': np.array([1,1,1,10000000,10000000,10000000]),# reinsertion rate inv. [time] 
    'g'  : g,                                           # growth conv. factors  [1/energy]
    'm'  : m,                                           # maintainance requ.    [energy/time]
    
    # sor algorithm parameters
    'n'  : n,                                           # grid points in each dim
    'sor': 1.80,                                        # relaxation parameter
    'L0' : 40,                                          # grid true size        [length]
    'D'  : 1e3,                                         # diffusion constant    [area/time] 
    'acc': 1e-6,                                        # maximum accepted stopping criterion   
    'ref': 1                                            # number of grid refinements to perform 
}

# make matrices
up_mat   = np.array([[1.,0.,0.,0.,0.,0.8],
                     [0.,1.,0.,0.8,0.,0.],
                     [0.,0.,1.,0.,0.8,0.]])
sign_mat = np.array([[1.,0.,0.,0.,0.,-1],
                     [0.,1.,0.,-1,0.,0.],
                     [0.,0.,1.,0.,-1,0.]])                   
met_mat  = np.array([[0.,0.,0.,0.,0.,0.],
                     [0.,0.,0.,0.,0.,0.],
                     [0.,0.,0.,0.,0.,0.],
                     [1.,0.,0.,0.,0.,0.],
                     [0.,1.,0.,0.,0.,0.],
                     [0.,0.,1.,0.,0.,0.]])
mat_ess   = np.array([[0.,0.,0.,0.,0.,0.],
                     [0.,0.,0.,0.,0.,0.],
                     [0.,0.,0.,0.,0.,0.]])
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
steps,R_fin,N_fin = run(2000,R0,N,param,mat)

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




