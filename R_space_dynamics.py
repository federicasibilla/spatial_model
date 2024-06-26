""" 
    R_space_dynamics.py: functions to simulate the steady state of nutrients diffusion in a 
    spatially structured consumer-resources model

    CONTAINS:
        - f: the temporal dynamics taking place at each spatial site
        - SOR: function to apply the successive over relaxation method to find 
               equilibrium distributions of resources, solves diffusion=f
        - SOR_multigrid: function to implement the SOR solution in a multigrid refinement context
        - change_grid: function needed for the multigrid method to go interpolate and restrict between 
                       grids at different refinement levels

"""

from audioop import mul
import numpy as np
import matplotlib.pyplot as plt

from time import time
from scipy.interpolate import RectBivariateSpline
from visualization import *

#-----------------------------------------------------------------------------------
# define f(R(x,y)) where R is the vector containing the concentrations, f the vector
# yielding the concentration dynamics f is the right hand side term of the equation
# ∇^2(R) = f(R(x,y))

def f(R,R0,N,param,mat):
    """
    R: resources vector (intended at point x,y)
    R0: initial concentration at specific point
    N: species vector (intended at point x,y)
    param: dictionary containing parameters
    mat: dictionary containing matrices

    returns a vector with value of dR/dt at a specific grid point

    """
    # check essential nutrients presence (at each site)
    up_eff = np.ones((mat['uptake'].shape))
    for i in range(len(N)):
        # calculate essential nutrients modulation for each species
        if (np.sum(mat['ess'][i]!=0)):
            mu = np.min(R[mat['ess'][i]==1]/(R[mat['ess'][i]==1]+1))
        else:
            mu = 1
        up_eff[i]=mat['uptake'][i]*mu
    # resource loss due to uptake
    out = np.dot((R*up_eff/(1+R)).T,N.T)
    # resource production due to metabolism
    inn = np.dot(((1/param['w']*(np.dot(param['w']*param['l']*R*up_eff/(1+R),mat['met'].T)))*mat['spec_met']).T,N.T)
    # resource replenishment
    ext = 1/param['tau']*(R0-R)

    return (ext+inn-out)

#-----------------------------------------------------------------------------------
# define SOR(param,mat,acc) function to implement the SOR algorithm for integration 
# on a grid; assumes fixed boundary conditions

def SOR(R,N,param,mat):
    """
    R: current resource concentration (on current grid) matrice n x n x n_R
    N: state vector for individuals (complete m4atrix on the current grid)
    param: dictionary containing parameters
    mat: dictionary containing matrices
    acc: desired accuracy, determins when to stop algorithm

    returns the equilibrium matrix n x n x n_r with concentrations and prints execution time

    """

    n   = N.shape[0]               # grid size
    n_r = param['R0'].shape[2]     # numer of nutrients
    dx  = param['L0']/n            # step size

    # initialization and boundary conditions setting
    R0 = change_grid(param['R0'],n)
    R  = R0.copy()

    # keep track of size of update
    delta = np.zeros((n+2,n+2,n_r))
    delta_max = np.max(R)
    delta_max_list = [2,1]

    # SOR algorithm and keep track of time
    t1 = time()

    # calculate value of f(R(x,y)) and store in an n x n x n_r matrix
    ff = np.zeros((n,n,n_r))
 
    # loop on each grid point and calculate vector there
    for i in range(n):
        for j in range(n):
            ff[i,j,:] = f(R[i,j,:],R0[i,j,:],N[i,j,:],param,mat)
    #for i in range(n_r):
    #    ff[:,:,i]=ff[:,:,i]-np.mean(ff[:,:,i])

    ff_temp = np.zeros((n+2,n+2,n_r))
    ff_temp[1:n+1,1:n+1,:]=ff
    ff = ff_temp

    # create ghost points initial guess
    R = np.zeros((n+2,n+2,n_r))
    #impose dirichelet BC through ghost points
    R[0, 1:-1,:]  = R0[0,:,:] #R[-2, 1:-1,:] 
    R[-1, 1:-1,:] = R0[-1,:,:] #R[1, 1:-1,:]  
    R[1:-1, 0,:]  = R0[:,0,:] #R[1:-1, -2,:] 
    R[1:-1, -1,:] = R0[:,-1,:] #R[1:-1, 1,:]
    R[0, -1,:]    = np.mean([R0[0,-2,:],R0[1,-1,:]],axis=0) #R[-2, 1,:]   
    R[-1, 0,:]    = np.mean([R0[-2,0,:],R0[-1,1,:]],axis=0) #R[1, -2,:]   
    R[-1, -1,:]   = np.mean([R0[-1,-1,:],R0[-2,-2,:]],axis=0) #R[1, 1,:]
    R[0,0,:]      = np.mean([R0[0,1,:],R0[1,0,:]],axis=0) #R[-2,-2,:]
    
    while (np.abs(delta_max_list[-1]-delta_max_list[-2])>param['acc']):

        # implement red/blue updates to make starting point ininfluential
        # loop on red
        for i in np.arange(1,n+1):
            start = 2 if i%2==0 else 1
            for j in np.arange(start,n+1,2):
                # next index for pbc (for the previous one pbc are automatic)
                delta[i,j] = (0.25+1e-13)*(R[i,j+1,:]+R[i,j-1,:]+R[i+1,j,:]+R[i-1,j,:]+(dx**2/param['D'])*ff[i,j,:])-R[i,j,:]
                R[i,j,:] += delta[i,j,:]*param['sor']
        # loop in blue
        for i in np.arange(1,n+1):
            start = 1 if i%2==0 else 2
            for j in np.arange(start,n+1,2):
                # next index for pbc (for the previous one pbc are automatic)
                delta[i,j] = (0.25+1e-13)*(R[i,j+1,:]+R[i,j-1,:]+R[i+1,j,:]+R[i-1,j,:]+(dx**2/param['D'])*ff[i,j,:])-R[i,j,:]
                R[i,j,:] += delta[i,j,:]*param['sor']
        
        # check updates
        delta_max = np.max(np.abs(delta))
        delta_max_list.append(delta_max)

        print("N_iter %d delta_max %e\r" % (len(delta_max_list), delta_max), end='')
   
    t2 = time()
    print("\nTotal running time: %.2f min" % ((t2-t1)/60))
    print("Code speed: %.1f iterations per second" %(len(delta_max_list)/(t2-t1)))

    plt.plot(delta_max_list[2:])
    plt.yscale('log') 
    plt.xlabel('Iteration')
    plt.ylabel('Delta max')
    plt.title('SOR algortihm convergence')
    plt.savefig('error.png')

    return R

#---------------------------------------------------------------------------------------------------
# define SOR_multigrid() function to implement the SOR algorithm on a grid that keeps getting refined
# the solution at the coarser grid is the initial guess of the finer grid

def SOR_multigrid(n_0,n_ref,R,N,param,mat):
    """
    n_0: initial grid size (corresponds to species grid size)
    n_ref: number of refinements
    R: current resource concentration (on current grid) matrice n x n x n_R
    N: state vector for individuals (complete matrix on the current grid)
    param: dictionary containing parameters
    mat: dictionary containing matrices
    acc: desired accuracy, determins when to stop algorithm

    returns R equilibrium concentration solved on the finer grid and averaged out at each
    coarse grid point

    """
    print('grid size is: ', R.shape[0])
    R = SOR(R,N,param,mat)                              # first run

    for i in range(2,n_ref+2):

        print('grid refinement number:', i-1)
        n = n_0*i                                       # current grid size
        print('grid size is: ', n)

        N_new = np.zeros((n,n,N.shape[2]))              # new grid for individuals
        # inverse restriction of individuals grid
        for k in range(n_0):
            for j in range(n_0):
                N_new[i*k:i*k+i,i*j:i*j+i,:]=N[k,j,:]
        R_new   = change_grid(R,n)                      # translate R on a refined grid 

        R_eq = SOR(R_new,N_new,param,mat)               # solve with SOR on this grid
        R = R_eq                                        # use this new solution as next initial guess

    # go back to initial grid size
    R_res = change_grid(R,n_0)

    return R_res

#-----------------------------------------------------------------------------------------------
# define change_grid(matrix,n_in,n_fin) function to translate a matrix defined on a coarser grid
# to a patrix defined on a finer grid (by interpolating) and viceversa (by restricting)

def change_grid(matrix,n_fin):
    """
    matrix: matrix that we wish to translate to a new grid
    n_fin: final grid size

    returns a matrix translated on the n_finxn_fin grid

    """
    # extract initial grid size
    n_in = matrix.shape[0]

    # create empty matrix of desired size
    interpolated_matrix = np.zeros((n_fin,n_fin, matrix.shape[2]))
    
    # go to finer grid: interpolate
    if n_fin>n_in:
        # creation of old and new grid
        x, y = np.linspace(0, n_in, n_in),np.linspace(0, n_in, n_in)
        x_fin, y_fin = np.meshgrid(np.linspace(0, n_in, n_fin), np.linspace(0, n_in, n_fin))
        for i in range(matrix.shape[2]):  
            # interpolate using RectBivariateSpline
            interpolator = RectBivariateSpline(y, x, matrix[:, :, i])
            interpolated_matrix[:, :, i] = interpolator.ev(y_fin, x_fin)
        return interpolated_matrix
    # go from finer to coarser: restrict
    else:
        scale = n_in // n_fin
        for i in range(n_fin):
            for j in range(n_fin):
                interpolated_matrix[i,j,:] = matrix[i*scale, j*scale, :]
        return interpolated_matrix


