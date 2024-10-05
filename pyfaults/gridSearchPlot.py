##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# generates a scatter plot with grid search results ----------
def gridSearchScatterPlot(r2vals, sVec, prob, inclZ=False):
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection='3d')

    sx = []
    sy = []
    sz = []
    for i in range(len(sVec)):
        sx.append(sVec[i][0])
        sy.append(sVec[i][1])
        sz.append(sVec[i][2])
    
    if inclZ is False:
        for i in range(len(r2vals)):
            ax.set_zlabel('Probability')
            p = ax.scatter(sx, sy, prob, c=r2vals, cmap=cm.coolwarm)
    elif inclZ is True:
        for i in range(len(r2vals)):
            ax.set_zlabel('Sz')
            p = ax.scatter(sx, sy, sz, c=r2vals, cmap=cm.coolwarm)
    
    ax.set_xlabel('Sx')
    ax.set_ylabel('Sy')
    
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(0,1)
    
    fig.colorbar(p, shrink=0.5, aspect=5, location='left')
    
    plt.tight_layout()
    return p

# generates a surface plot with grid search results ----------
def gridSearchSurfacePlot(r2vals, sVec, prob, inclZ=False):
    
    fig = plt.figure()
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    sx = []
    sy = []
    sz = []
    for i in range(len(sVec)):
        sx.append(sVec[i][0])
        sy.append(sVec[i][1])
        sz.append(sVec[i][2])
    
    if inclZ is False:
        for i in range(len(r2vals)):
            ax.set_zlabel('Probability')
            z = np.array([prob, r2vals])
            p = ax.plot_surface(sx, sy, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    elif inclZ is True:
        for i in range(len(r2vals)):
            ax.set_zlabel('Sz')
            z = np.array([sz, r2vals])
            p = ax.plot_surface(sx, sy, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    
    ax.set_xlabel('Sx')
    ax.set_ylabel('Sy')
    
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(0,1)
    
    fig.colorbar(p, shrink=0.5, aspect=5, location='left')
    
    plt.tight_layout()
    return p
        