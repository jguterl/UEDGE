#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:52:13 2020

@author: jguterl
"""
import uedge
import matplotlib.pyplot as plt
import os
import numpy
import os, sys, math, string, re
from pathlib import Path
import numpy
from . import UEDGEToolBox
from uedge import *

class UEDGESimulation(object):
       
    def __init__(self,  *args, **kwargs):
        # Import Uedge packages as attribute of the class instance 
        # WARNING: this is not a deep copy!!!
        self.ListPkg=UEDGEToolBox.GetListPackage()
        print('Loaded Packages:',self.ListPkg)
        for pkg in self.ListPkg:
            exec('self.' + pkg + '=' + pkg,globals(),locals())
        
            
    def ReadInput(self,FileName:str):
        for pkg in self.ListPkg:
            exec('from uedge import '+pkg)
        print('### Looking for input file:',FileName)
        if not os.path.exists(FileName):
            print('### Cannot find {}'.format(os.path.abspath(FileName)))
            FileName_=os.path.join(Settings.InputDir,FileName)
            if not os.path.exists(FileName_):
                raise IOError('Cannot find the requested file in InputDir either:{} '.format(FileName))
            else:
                print('### File found in InputDir:{}'.format(FileName_))
        else: 
            FileName_=os.path.abspath(FileName)
            print('### File found :{}'.format(FileName_))
            
        print('### Loading {} '.format(FileName_))    
        
        f=open(FileName_,'r')
        lines=f.read()
        f.close()
        print('### Input:{}'.format(lines))
        exec(lines)    
        
         
            # except Exception as e:
            #     print(repr(e))
            #     print('Could not import:', os.path.join(Settings.InputDir,FileName)) 
            
            
    def ReadDivertorPlateFile(self,FileName,Verbose=False):
    
        Dic={'rplate1':[],'zplate1':[],'rplate2':[],'zplate2':[],'Comment':'No comment','FileName':FileName,'InnerDivPlateFile':None,
        'OuterDivPlateFile':None}
        
        
        try:
            with  open(os.path.join(FileName)) as f:
                exec(f.read(),Dic)
        except:
            raise ValueError('Cannot Read File:'+FileName)
        
        
        
        if Dic['InnerDivPlateFile'] is not None:
             InnerData=numpy.loadtxt(Dic['InnerDivPlateFile'])
             if Verbose: print(InnerData)   
             Dic['rplate1']=list(InnerData[:,0])
             Dic['zplate1']=list(InnerData[:,1])
        
        if Dic['OuterDivPlateFile'] is not None:
            OuterData=numpy.loadtxt(Dic['OuterDivPlateFile'])  
            Dic['rplate2']=list(OuterData[:,0])
            Dic['zplate2']=list(OuterData[:,1])
        # Check if data are correct
        # (len(rplate2)<2 or len(zplate2)<2 or len(rplate1)<2 or len(zplate1)<2):
         #       raise ValueError('Wrong size of plates coordinates')
        return Dic
    
    # outer is 2, inner is 1


    def SetDivertorPlates(self,FileName):
        Dic=ReadDivertorPlateFile(self,FileName)
        if len(Dic['rplate1'])<2 or len(Dic['rplate1'])!=len(Dic['zplate1']):
            raise ValueError('wrong dimension of coordinates of plate #1')
        if len(Dic['rplate2'])<2 or len(Dic['rplate2'])!=len(Dic['zplate2']):
            raise ValueError('wrong dimension of coordinates of plate #2')
        print('FileName:'+Dic['FileName'])
        print('Comment:'+Dic['Comment'])
        self.grd.nplate1=len(Dic['rplate1'])
        self.grd.nplate2=len(Dic['rplate2'])
        self.grd.gchange("Mmod")
        self.grd.rplate2=Dic['rplate2']
        self.grd.zplate2=Dic['zplate2']
        self.grd.rplate2=Dic['rplate1']
        self.grd.zplate1=Dic['zplate1']
    
    def PlotDivertorPlates(self,FileName=None):
        if type(FileName)==str:
            Dic=ReadDivertorPlateFile(FileName)
            if len(Dic['rplate1'])<2 or len(Dic['rplate1'])!=len(Dic['zplate1']):
                print('rplate1=',Dic['rplate1'])
                print('zplate1=',Dic['zplate1'])
                raise ValueError('wrong dimension of coordinates of plate #1')
            if len(Dic['rplate2'])<2 or len(Dic['rplate2'])!=len(Dic['zplate2']):
                raise ValueError('wrong dimension of coordinates of plate #2')
        
            r1=Dic['rplate1']
            r2=Dic['rplate2']
            z1=Dic['zplate1']
            z2=Dic['zplate2']
            plt.suptitle=Dic['Comment']
            plt.xlabel('R[m]')
            plt.ylabel('Z[m]')
            plt.axis('equal')
            plt.grid(True)
        else:
            r1=grd.rplate1
            r2=grd.rplate2
            z1=grd.zplate1
            z2=grd.zplate2
        
        plt.plot(r1,z1,color='r')
        plt.plot(r2,z2,'g')
    
Sim=UEDGESimulation()

def ReadInput(FileName):
    Sim.ReadInput(FileName)
    
                   
        
        
     
