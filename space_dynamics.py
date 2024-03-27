""" 
    space_dynamics.py: functions to simulate the steady state of nutrients diffusion in a 
    spatially structured consumer-resources model

    Author: Federica Sibilla, University of Lausanne

"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from time import time
from mpl_toolkits import mplot3d

sns.set(style='whitegrid')

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
    out = (-R*mat['uptake']/(1+R*mat['uptake'])).T
    # resource production due to metabolism
    inn = (1/param['w']*(np.dot(param['w']*(1-param['l'])*R*mat['uptake']/(1+R*mat['uptake']),mat['met'].T))).T
    # resource replenishment
    ext = 1/param['tau']*(R0-R)

    return -(ext+np.dot(inn+out,N))


# define SOR(param,mat,acc) function to implement the SOR algorithm for integration 
# on a grid; assumes fixed boundary conditions

def SOR(N,param,mat,acc):
    """
    N: state vector for individuals (complete matrix on the grid)
    param: dictionary containing parameters
    mat: dictionary containing matrices
    acc: desired accuracy, determins when to stop algorithm

    returns the equilibrium matrix n x n x n_r with concentrations and prints execution time

    """
    n   = param['n']               # grid size
    n_r = param['R0'].shape[2]     # numer of nutrients
    dx  = param['L0']/n            # step size

    # initialization and boundary conditions setting
    R  = param['R0'].copy()        # matrice n x n x n_R
    R0 = param['R0'].copy()        # store initial values to use in replenishment 
    for i in range(len(R)):
        R[i,0,:]  = 100              # Left boundary
        R[i,-1,:] = 100              # Right boundary
        R[0,i,:]  = 100              # Lower boudary
        R[-1,i,:] = 100              # Upeer boundary

    # keep track of size of update
    delta = R.copy()
    delta_max = 1
    delta_max_list = []

    # SOR algorithm and keep track of time
    t1 = time()

    # calculate value of f(R(x,y)) and store in an n x n x n_r matrix
    ff = np.zeros((n,n,n_r))
    # loop on each grid point and calculate vector there
    for i,point_r,point_r0,point_n in zip(range(n),R,R0,N):
        for j,r,r0,nn in zip(range(n),point_r,point_r0,point_n):
            ff[i,j,:] = f(r,r0,nn,param,mat)

    while delta_max > acc:

        # implement red/blue updates to make starting point ininfluential
        # loop on red
        for i in np.arange(1,n-1):
            start = 2 if i%2==0 else 1
            for j in np.arange(start,n-1,2):
              delta[i,j] = 0.25*(R[i,j+1,:]+R[i,j-1,:]+R[i+1,j,:]+R[i-1,j,:]+dx**2*ff[i,j,:])-R[i,j,:]
              R[i,j,:] += delta[i,j,:]*param['sor']
        # loop in blue
        for i in np.arange(1,n-1):
            start = 1 if i%2==0 else 2
            for j in np.arange(start,n-1,2):
              delta[i,j] = 0.25*(R[i,j+1,:]+R[i,j-1,:]+R[i+1,j,:]+R[i-1,j,:]+dx**2*ff[i,j,:])-R[i,j,:]
              R[i,j,:] += delta[i,j,:]*param['sor']
        # check updates
        delta_max = np.max(np.abs(delta[1:n-1,1:n-1,:]))
        delta_max_list.append(delta_max)

        print("N_iter %d delta_max %e\r" % (len(delta_max_list), delta_max), end='')

    t2 = time()
    print("\nTotal running time: %.2f min" % ((t2-t1)/60))
    print("Code speed: %.1f iterations per second" %(len(delta_max_list)/(t2-t1)))

    return R

# define R_ongrid(R) to visualize the equilibrium concentrations for each resource

def R_ongrid(R,param):
    """
    R: the matrix n x n x n_r with equilibrium concentrations
    param: dictionary containing parameters

    plots the equilibrium concentrations for each nutrient

    """

    # create the grid
    x = np.arange(param['n'])
    y = np.arange(param['n'])
    X, Y = np.meshgrid(x, y)

    # R matrix as function of x and y plot (one plot per nutrient)
    n_r = R.shape[2]
    fig = plt.figure(figsize=(18, 6))
    axs = [fig.add_subplot(1, n_r, i+1, projection='3d') for i in range(n_r)]
    cmaps = ['plasma', 'ocean','CMRmap'] # att.va bene solo per 2 nut

    for i in range(n_r):
        ax = axs[i]
        #cmap=cmaps[i]
        surf = ax.plot_surface(X, Y, R[:, :, i], cmap='ocean', edgecolor='none')
        fig.colorbar(surf, ax=ax, label='Concentration')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('Concentration')
        ax.set_title('Resource {}'.format(i+1))

    plt.show()

    return