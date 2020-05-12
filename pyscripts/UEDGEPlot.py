"""
Created on Thu Feb 20 16:29:38 2020

@author: jguterl
"""


import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
from uedge import *
from mpldatacursor import datacursor
#decorator to set ax or new fig is request
# def AutoFig(f):
#     def wrapper_do_twice(*args, **kwargs):
#         if kwargs.get('ax') is None:
#             ax=plt.gca()
#         if kwargs.get('NewFig') is True:
#             fig,ax
#         func(*args, **kwargs)
#         return func(*args, **kwargs)
#     return wrapper_do_twice


#%%


#%%

    #plt.show(block=True)
        
    
def PlotMesh2(ax=None,iso=False,Cell=None):
    plt.ion()         
    if ax is None:
        fig,ax = plt.subplots()

    if (iso):
        ax.set_aspect('equal', 'datalim')
    else:
        ax.set_aspect('auto', 'datalim')

    
    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            ax.plot(com.rm[ix,iy,[1,2,4,3,1]],
                     com.zm[ix,iy,[1,2,4,3,1]], 
                     color="b")
    

    
    #ax.figure.show()     
    ax.xaxis.label('R [m]')
    ax.yaxis.ylabel('Z [m]')
    ax.figure.suptitle('UEDGE mesh')
    ax.grid(True)
    plt.show()
    
    return