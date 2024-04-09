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

import numpy as np

from time import time
from scipy.interpolate import RectBivariateSpline

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
    # resource loss due to uptake
    out = np.dot((-R*mat['uptake']/(1+R)).T,N.T)
    # resource production due to metabolism
    inn = np.dot((1/param['w']*(np.dot(param['w']*param['l']*R*mat['uptake']/(1+R),mat['met'].T))).T,N.T)
    # resource replenishment
    ext = 1/param['tau']*(R0-R)

    return (ext+inn+out)

#-----------------------------------------------------------------------------------
# define SOR(param,mat,acc) function to implement the SOR algorithm for integration 
# on a grid; assumes fixed boundary conditions

def SOR(R,N,param,mat,acc,n):
    """
    R: current resource concentration (on current grid) matrice n x n x n_R
    N: state vector for individuals (complete matrix on the current grid)
    param: dictionary containing parameters
    mat: dictionary containing matrices
    acc: desired accuracy, determins when to stop algorithm
    n: current grid size

    returns the equilibrium matrix n x n x n_r with concentrations and prints execution time

    """
    n_r = param['R0'].shape[2]     # numer of nutrients
    dx  = param['L0']/n            # step size

    # initialization and boundary conditions setting
    R0 = change_grid(param['R0'],param['R0'].shape[0],n)

    # keep track of size of update
    delta = R.copy()
    delta_max = np.max(R)
    delta_max_list = []

    # SOR algorithm and keep track of time
    t1 = time()

    # calculate value of f(R(x,y)) and store in an n x n x n_r matrix
    ff = np.zeros((n,n,n_r))
 
    # loop on each grid point and calculate vector there
    for i in range(n):
        for j in range(n):
            ff[i,j,:] += f(R[i,j,:],R0[i,j,:],N[i,j,:],param,mat)
    
    while ((delta>acc*R0).any()):

        # implement red/blue updates to make starting point ininfluential
        # loop on red
        for i in np.arange(0,n):
            start = 1 if i%2==0 else 0
            for j in np.arange(start,n,2):
              # next index for pbc (for the previous one pbc are automatic)
              nexti = (i+1) % n
              nextj = (j+1) % n
              delta[i,j] = 0.25*(R[i,nextj,:]+R[i,j-1,:]+R[nexti,j,:]+R[i-1,j,:]+(dx**2/param['D'])*ff[i,j,:])-R[i,j,:]
              R[i,j,:] += delta[i,j,:]*param['sor']
        # loop in blue
        for i in np.arange(0,n):
            start = 0 if i%2==0 else 1
            for j in np.arange(start,n,2):
              # next index for pbc (for the previous one pbc are automatic)
              nexti = (i+1) % n
              nextj = (j+1) % n
              delta[i,j] = 0.25*(R[i,nextj,:]+R[i,j-1,:]+R[nexti,j,:]+R[i-1,j,:]+(dx**2/param['D'])*ff[i,j,:])-R[i,j,:]
              R[i,j,:] += delta[i,j,:]*param['sor']
        # check updates
        delta_max = np.max(np.abs(delta))
        delta_max_list.append(delta_max)

        print("N_iter %d delta_max %e\r" % (len(delta_max_list), delta_max), end='')

    t2 = time()
    print("\nTotal running time: %.2f min" % ((t2-t1)/60))
    print("Code speed: %.1f iterations per second" %(len(delta_max_list)/(t2-t1)))

    return R

#---------------------------------------------------------------------------------------------------
# define SOR_multigrid() function to implement the SOR algorithm on a grid that keeps getting refined
# the solution at the coarser grid is the initial guess of the finer grid

def SOR_multigrid(n_0,n_ref,R,N,param,mat,acc):
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
    R = SOR(R,N,param,mat,acc,n_0)                # first run

    for i in range(1,n_ref+1):
        n = n_0*2*i                               # current grid size
        N_new = np.zeros((n,n,N.shape[2]))        # new grid for individuals
        # inverse restriction of individuals grid
        for k in range(n_0):
            for j in range(n_0):
                N_new[i*k:i*k+i,i*j:i*j+i,:]=N[k,j,:]
        R_new   = change_grid(R,R.shape[0],n)     # translate R on a refined grid 
        SOR(R_new,N_new,param,mat,acc/(3**i),n)   # solve with SOR on this grid
        R = R_new                                 # use this new solution as next initial guess

    # go back to initial grid size
    R_res = change_grid(R,R.shape[0],n_0)

    return R_res

#-----------------------------------------------------------------------------------------------
# define change_grid(matrix,n_in,n_fin) function to translate a matrix defined on a coarser grid
# to a patrix defined on a finer grid (by interpolating) and viceversa (by restricting)

def change_grid(matrix,n_in,n_fin):
    """
    matrix: matrix that we wish to translate to a new grid
    n_in: initial grid size
    n_fin: final grid size

    returns a matrix translated on the n_finxn_fin grid

    """

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


