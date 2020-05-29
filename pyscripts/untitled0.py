#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 16:11:32 2020

@author: jguterl
"""
import numpy as np
import matplotlib.pyplot as plt
def DrawSeparatrix1D(self,ax=None,Target='outer'):
    if ax is None:
        ax=plt.gca()
    if Target=='outer':
        idx=self.Grid['rm'].shape[0]-1
    elif Target=='inner':
        idx=0
    else:
        print('Unknown target. Cannot plot separatrix in 1D...')
        return
    Rsep=self.Grid['rm'][idx,self.Grid['iysptrx'],3]
    ax.plot((Rsep,Rsep),ax.get_ylim(),'b--')
    
def PlotRadial(self,DataField,Data,Grid,xidx,ax=None,Label=None,DataScale='linear',iSp=0):
        if ax is None:
            ax=plt.gca()
    
        if Grid is None and hasattr(self,'GetGrid'):
            Grid=self.GetGrid()
        else:
            if XType!='idx':
                print('Cannot find a grid to plot against XType={}... Skipping plot'.format(XType))
                return
            
        if type(DataField)!=list:
            DataField=[DataField]
                
        if Data is None:
            self.SetData()
        Data=self.Data
        if Data is None:
            print('Cannot access data... Exiting ...')
            return
    
        for F in DataField:
            if self.Verbose: print('Plotting {} target profiles of {}'.format(Target,F))
            DataPlot=Data.get(F.lower())
            if DataPlot is None:
                print('Cannot find the field {} in available data. Skipping the plot...'.format(F))
                continue
            
            DataShape=DataPlot.shape
            if len(DataShape)<3:
                if iSp>0:
                    print('Cannot get data for species of index:{}. Shape of Data:{}. Skipping....'.format(iSp,DataShape))
                    continue
            elif iSp>=DataShape[3]:
                print('Cannot get data for species of index:{}. Shape of Data:{}. Skipping....'.format(iSp,DataShape))
                continue
            
            nx=DataPlot.shape[0]-1
            #if len(DataPlot.shape)>1:
            ny=DataPlot.shape[1]-1
            #else:
                #ny=None
                
            if Target=='outer':
                Idxg=(nx,np.arange(ny+1))
            elif Target=='inner':
                Idxg=(0,np.arange(ny+1))
            else:
                print('Target must be "inner" or "outer"... Exiting...')
                return 
            
            if len(DataShape)>2:
                Idx=(Idxg[0],Idxg[1],iSp)
            else:
                Idx=Idxg
            
            if XType=='r':
                X=self.Grid['rm'][:,:,0]
            elif XType=='psi':
                X=self.Grid['psi'][:,:,0]
            elif XType=='idx':
                X=np.arange(nx)
            else:
                print('XType must be "idxx","psi" or "r" ... Exiting...')
                return 
            
            LabelPlot="{}:{}".format(Target,F)
            ax.plot(X[Idxg],DataPlot[Idx],label=LabelPlot,marker='o')
            ax.legend()
            
            
def PlotTarget(self,DataField,Data=None,Grid=None,ax=None,Target='outer',Label=None,DataScale='linear',XType='r',iSp=0):
        if ax is None:
            ax=plt.gca()
    
        if Grid is None and hasattr(self,'GetGrid'):
            Grid=self.GetGrid()
        else:
            if XType!='idx':
                print('Cannot find a grid to plot against XType={}... Skipping plot'.format(XType))
                return
            
        if type(DataField)!=list:
            DataField=[DataField]
                
        if Data is None:
            self.SetData()
        Data=self.Data
        if Data is None:
            print('Cannot access data... Exiting ...')
            return
    
        for F in DataField:
            if self.Verbose: print('Plotting {} target profiles of {}'.format(Target,F))
            DataPlot=Data.get(F.lower())
            if DataPlot is None:
                print('Cannot find the field {} in available data. Skipping the plot...'.format(F))
                continue
            
            DataShape=DataPlot.shape
            if len(DataShape)<3:
                if iSp>0:
                    print('Cannot get data for species of index:{}. Shape of Data:{}. Skipping....'.format(iSp,DataShape))
                    continue
            elif iSp>=DataShape[3]:
                print('Cannot get data for species of index:{}. Shape of Data:{}. Skipping....'.format(iSp,DataShape))
                continue
            
            nx=DataPlot.shape[0]-1
            #if len(DataPlot.shape)>1:
            ny=DataPlot.shape[1]-1
            #else:
                #ny=None
                
            if Target=='outer':
                Idxg=(nx,np.arange(ny+1))
            elif Target=='inner':
                Idxg=(0,np.arange(ny+1))
            else:
                print('Target must be "inner" or "outer"... Exiting...')
                return 
            
            if len(DataShape)>2:
                Idx=(Idxg[0],Idxg[1],iSp)
            else:
                Idx=Idxg
            
            if XType=='r':
                X=self.Grid['rm'][:,:,0]
            elif XType=='psi':
                X=self.Grid['psi'][:,:,0]
            elif XType=='idx':
                X=np.arange(nx)
            else:
                print('XType must be "idxx","psi" or "r" ... Exiting...')
                return 
            
            LabelPlot="{}:{}".format(Target,F)
            ax.plot(X[Idxg],DataPlot[Idx],label=LabelPlot,marker='o')
            ax.legend()