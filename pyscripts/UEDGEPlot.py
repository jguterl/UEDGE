"""
Created on Thu Feb 20 16:29:38 2020

@author: jguterl
"""


import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import numpy as np
from uedge import *
from uedge.UEDGEMesh import UEDGEMesh
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

class UEDGEPlot1DBase():
    def __init__(self):
        self.eV=1.602176634e-19
    
    def PlotData1DBase(x,data,ax=None,DataLim=None,DataScale='linear',Verbose=False):
        
        if ax is None:
            ax=plt.gca()  
            
        ax.plot(x,data)
        
    
            
class UEDGEPlot2DBase():
    def __init__(self):
        self.eV=1.602176634e-19
        
    def PlotSeparatrix(self,rm,zm,iysptrx,ax=None,color='r',linewidth=1,**kwargs):
       sepx=np.concatenate((rm[:,iysptrx,3],np.array([rm[-1,iysptrx,4]])))
       sepy=np.concatenate((zm[:,iysptrx,3],np.array([zm[-1,iysptrx,4]])))
       
       if ax is None:
           ax=plt.gca()
       ax.plot(sepx,sepy,color=color,linewidth=linewidth,**kwargs)
       
    @staticmethod    
    def PlotData2DBase(r,z,data,ax=None,ColorMap='jet',DataLim=None,DataScale='linear',Verbose=False):
        """Plot UEDGE grid."""
        if ColorMap not in matplotlib.pyplot.colormaps():
            print('ColorMap {} not defined in matplotlib...')
            print('ColorMap must be chosen in the following list:')
            print(matplotlib.pyplot.colormaps())
            return
        
        if ax is None:
            ax=plt.gca()
        
        
        def onpick(evt):
            if isinstance(evt.artist,matplotlib.collections.PatchCollection):
                print(evt.artist.get_array()[evt.ind[0]])
            if evt.artist in Pos.keys():  
                annot.set_visible(False)
                annot.xy = Pos[evt.artist]
                evt.artist.set_facecolor='blue'
                evt.artist.set_fill=True
                annot.set_text(Dic[evt.artist])
                annot.set_visible(True)
            if evt.mouseevent.button == 3:
                annot.set_visible(False)    
            #ax.figure.canvas.draw_idle()
    
        Nx=len(r)
        Ny=len(r[0])
        data=data.reshape(Nx*Ny)
        if DataLim is None:
            DataLim=(min(data),max(data))
        if Verbose:
            print('DataLim:',DataLim)
            
        
        patches=[]
        if ax is None:
            ax=plt.gca()
            
        idx=[np.array([1,2,4,3,1])]
        Dic={}
        Pos={}
        Obj={}
        for i in range(Nx):
            for j in range(Ny):
                    Data=np.concatenate((r[i,j,idx],z[i,j,idx])).reshape(2,5).T
                    p=matplotlib.patches.Polygon(Data,closed=True,edgecolor=None,label='ix={},iy={}'.format(i,j),picker=5)
                    
                    patches.append(p)
                    
                    Obj[p]=p        
                    Dic[p]='ix={},iy={}'.format(i,j)
                    Pos[p]=(r[i,j,0],z[i,j,0])
        if Verbose:
            print('DataLim=',DataLim)
        Collec=matplotlib.collections.PatchCollection(patches)
        Collec.set_picker(True)
        if DataScale=='log':
            norm=matplotlib.colors.LogNorm(vmin=DataLim[0],vmax=DataLim[1])
        elif DataScale=='symlog':
            norm=matplotlib.colors.SymLogNorm(vmin=DataLim[0],vmax=DataLim[1])
        elif DataScale=='linear':
            norm=matplotlib.colors.Normalize(vmin=DataLim[0],vmax=DataLim[1])
        else:
            print('Unknow DataScale. Must be log|symlog|linear')
            return
        Collec.set_array(data)
        Collec.set_cmap(ColorMap)
        Collec.set_norm(norm)
        Collec.set_clim(vmin=DataLim[0],vmax=DataLim[1])
        #cmap = plt.get_cmap('jet')
        #colors = cmap(data)   
        
        
        ax.add_collection(Collec)
        aspect = 20
        pad_fraction = 0.5
        divider = make_axes_locatable(ax)
        width = axes_size.AxesY(ax, aspect=1./aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax = divider.append_axes("right", size=width, pad=pad)
        plt.colorbar(Collec,ax=cax,norm=norm)  
            
        ax.set_ylim(z.min(),z.max())
        ax.set_xlim(r.min(),r.max())
        annot = ax.annotate("", xy=(0,0), xytext=(-20,20),textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
         
        ax.figure.canvas.mpl_connect('pick_event', onpick)   
        ax.set_aspect('equal', 'box')
        

class UEDGEPlot(UEDGEPlot2DBase,UEDGEMesh):
    eV=1.602176634e-19
    PlasmaVars=['ni','up','te','ti','phi']
    def __init__(self,Verbose=False):
        self.Ncol=2
        self.Verbose=Verbose
        
    def Plot2D(self,DataField=[],Object=[],Grid=None,ColorMap='jet',DataLim=None,DataScale='linear',fig=None):
        Object=[self]+Object
        if DataField==[]:
            DataField=self.__class__.PlasmaVars
        self.__class__.PlotData2D(DataField,Object,Grid,ColorMap,DataLim,DataScale,fig,self.Verbose)
        
    def PlotSingle2D(self,DataField,Grid=None,ax=None,ColorMap='jet',DataLim=None,DataScale='linear',fig=None):
        Grid=self.GetGrid()
        self.__class__.PlotData2DBase(Grid['rm'],Grid['zm'],self.GetData(DataField),ax,ColorMap,DataLim,DataScale,self.Verbose)
                
    @staticmethod
    def PlotData2D(DataField=[],Object=[],Grid=None,ColorMap='jet',DataLim=None,DataScale='linear',fig=None,Verbose=False):
        if Verbose: print('Object:',Object)
        if type(Object)!=list:
            Object=list(Object)
        if len(Object)<1: 
            print('Nothing to plot...')
            return 
        if Verbose: print('Object:',Object)
        if DataField==[]: 
            print('At least one field is required for plotting....')
            return
        
        if not isinstance(DataField,list):
            DataField=list(DataField)
            
        Nrow=len(Object)
        Ncol=len(DataField)
        if Verbose: print('Ncol:',Ncol)
        if Verbose: print('Nrow:',Nrow)
        
        if len(Object)>1:
            FigTitle='Comparison '+'/'.join([o.GetCaseName() for o in Object])
        else:
            FigTitle=Object[0].GetCaseName()
            
        if fig is None:
            fig,ax=plt.subplots(Nrow,Ncol,sharex='all',sharey='all',num=FigTitle)
        else:
            ax=fig.subplots(Nrow,Ncol,sharex='all',sharey='all',num=FigTitle)
        
        if type(ax)!=np.ndarray:
            ax=[ax]
        i=0
        if Verbose: print('Object:',Object)
        for o in Object:
            CaseName=o.GetCaseName()
            if Grid is None:
                Grid=o.GetGrid()  
            for F in DataField:
                if Grid.get('rm') is None or Grid.get('zm') is None:
                    print('No grid available for {}. Skipping...'.format(CaseName))
                else:
                    Data=o.GetData(F)
                    if Verbose: print('Data.shape=',Data.shape)
                    if Verbose: print('rm.shape=',Grid['rm'].shape)
                    if Data is not None:
                        if F.lower()=='te' or F.lower()=='ti':
                            Data=Data/1.602176634e-19
                        if Verbose: print('Plotting ',F)
                        o.__class__.PlotData2DBase(Grid['rm'],Grid['zm'],Data,ax[i],ColorMap,DataLim,DataScale,Verbose)
                        ax[i].set_title('{}:{}'.format(CaseName,F))
                    else:
                        print('Cannot plot {F} for case {}'.format(F,CaseName))
                i=i+1
        plt.show()        
                
                
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