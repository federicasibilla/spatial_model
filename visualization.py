"""
visualization.py: file containing the functions to make plots 

CONTAINS:
    - R_ongrid: takes the matrix nxnxn_R with concentrations and returns the plot of such distributions 
                in 3D, one surface for each nutrient on the grid it is defined on
    - N_ongrid: takes the matrix nxnxn_s with ones and zeros representing which species is present
                on each sites and returns the plot, with one color for each species
    - G_ongrid: takes the nxn matrix of growth rates and returns the plot on grid
    - makenet:  function to draw the metabolic network, takes the metabolic matrix as input
    - vispreferences: function to visualize the uptake preferences

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import networkx as nx
import seaborn as sns

from mpl_toolkits import mplot3d
from networkx.drawing.nx_agraph import to_agraph
sns.set(style='whitegrid')

#---------------------------------------------------------------------------------
# define R_ongrid(R) to visualize the equilibrium concentrations for each resource

def R_ongrid(R):
    """
    R: the matrix n x n x n_r with concentrations

    plots the equilibrium concentrations for each nutrient

    """

    # create the grid
    x = np.arange(R.shape[0])
    y = np.arange(R.shape[0])
    X, Y = np.meshgrid(x, y)

    # R matrix as function of x and y plot (one plot per nutrient)
    n_r = R.shape[2]
    fig = plt.figure(figsize=(18, 6))
    axs = [fig.add_subplot(1, n_r, i+1, projection='3d') for i in range(n_r)]

    for i in range(n_r):
        ax = axs[i]
        surf = ax.plot_surface(X, Y, R[:, :, i], cmap='ocean', edgecolor='none')
        fig.colorbar(surf, ax=ax, label='Concentration',ticks=np.linspace(np.min(R[:,:,i]),np.max(R[:,:,i]),10))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('Concentration')
        ax.set_title('Resource {}'.format(i+1))

    plt.savefig('/Users/federicasibilla/Documenti/Tesi/Codice/spatial_model/R.png')
    plt.close()

    return

#---------------------------------------------------------------------------------
# define N_ongrid(R) to visualize the disposition of species 

def N_ongrid(N):
    """
    N: the nxnxn_s matrix containing nxn elements with length n_s composed by all zeros and
       one corresponding to the species present in the grid point (1 species per grid point)

    plots the grid with current species disposition

    """
    # define colors for species distinction
    cmap = plt.cm.get_cmap('cool', N.shape[2])  
    norm = mc.Normalize(vmin=0, vmax=N.shape[2]-1)

    # plot gird
    colors = cmap(norm(np.argmax(N, axis=2)))
    plt.figure(figsize=(8, 8))

    plt.imshow(colors, interpolation='nearest')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm,ticks=np.arange(N.shape[2]),label='Species')

    plt.savefig('/Users/federicasibilla/Documenti/Tesi/Codice/spatial_model/N0.png')
    plt.close()

    return

#---------------------------------------------------------------------------------
# define G_ongrid(G) function to visualize growth rates

def G_ongrid(G):
    """
    G: growth rates vector

    returns grid with color gradient corresponding to growth rates

    """
    # define colors for species distinction
    cmap = plt.cm.get_cmap('hot') 

    # plot gird
    plt.figure(figsize=(8, 8))

    plt.imshow(G, cmap=cmap)

    plt.colorbar(label='Growth rate')

    plt.savefig('/Users/federicasibilla/Documenti/Tesi/Codice/spatial_model/G.png')
    plt.close()

    return

#---------------------------------------------------------------------------------
# define makenet(met_matrix) to visualize the metabolic processes network, with
# resources as nodes and allocation magnitude as edges thikness

def makenet(met_matrix):
    """
    met_matrix: metabolic matrix, with resources as rows and columns and allocation rates as
                entries

    returns the graph of metabolic allocations

    """
    G = nx.DiGraph()

    for i in range(met_matrix.shape[0]):
        for j in range(met_matrix.shape[1]):
            G.add_edge(f"Res{j+1}", f"Res{i+1}", weight=met_matrix[i, j])

    # draw graph
    agraph = to_agraph(G)
    agraph.layout(prog='dot', args='-GK=0.5 -Gsep=3 -Ncolor=lightblue -Nstyle=filled -Npenwidth=2 -Ecolor=gray -Nnodesep=0.1')
    for edge in agraph.edges():
        weight = G[edge[0]][edge[1]]['weight']
        agraph.get_edge(*edge).attr['penwidth'] = weight * 5
    img = agraph.draw(format='png')
    with open('met_net.png', 'wb') as f:
        f.write(img)

    return

#---------------------------------------------------------------------------------
# defining vispreferences(up_mat) function to visualize the uptake preferences 
# of the different species

def vispreferences(up_mat):
    """
    up_mat: uptake matrix of the different species and resources

    returns a graph to visualize uptake preferences 

    """
    plt.figure(figsize=(10, 6))

    colors = plt.cm.viridis(np.linspace(0, 1, up_mat.shape[1]))  

    legend = 0
    for i in range(up_mat.shape[0]):
        offset = 0 
    
        for j in range(up_mat.shape[1]):
            lunghezza_segmento = up_mat[i, j]  
            if legend<up_mat.shape[1]:
                plt.bar(i, lunghezza_segmento, bottom=offset, width=0.8, color=colors[j], label=f'Res {j+1}')
                offset += lunghezza_segmento
                legend +=1
            else:
                plt.bar(i, lunghezza_segmento, bottom=offset, width=0.8, color=colors[j])
                offset += lunghezza_segmento

    plt.xlabel('Species')
    plt.ylabel('Uptake')
    plt.title('Consumer preferences')
    plt.legend()
    plt.grid(True) 

    plt.savefig('/Users/federicasibilla/Documenti/Tesi/Codice/spatial_model/uptake_pref.png')
    plt.close()

    return