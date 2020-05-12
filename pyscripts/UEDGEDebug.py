#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 19:20:36 2020

@author: jguterl
"" 16.84704033613164
   16.84704064249413
"""

import numpy as np
import os
import uedge
from uedge.UEDGEDoc import *
def DumpData(FileName,DataName):
    com='np.savetxt("{}",{})'.format(FileName,DataName)
    exec(com,globals(),locals())
    
def BatchDumpData(Suffix:str,Folder:str,VarList:list,Verbose=False):
    import uedge
    from uedge import bbb,com,grd,flx,api,aph
    for V in VarList:
        Result=SearchSilent(V,exact=True)
        if len(Result)>0:
            Package=Result[0]['Package']
            FileName=os.path.join(Folder,'{}_{}.txt'.format(V,Suffix))
            DataName='{}.{}'.format(Package,V)
            if Verbose:
                print(FileName)
                print(DataName)
        
            comm='np.savetxt("{}",{}.flatten(),header=str({}.shape))'.format(FileName,DataName,DataName)
            exec(comm,globals(),locals())
        else:
            print('skipping {}'.format(V))

            
def GetList(FileName):
    with open(FileName) as f:
        out=f.read().splitlines()
    return list(dict.fromkeys(out))

def Compare(Suffixes:list,Folder:str,VarList:list,Verbose=False,Precision=1e-15):
    Dic={}
    if len(Suffixes)!=2: raise ValueError('Suffixes must be a list of two elements')
    for V in VarList:
        for S in Suffixes:
            FileName=os.path.join(Folder,'{}_{}.txt'.format(V,S))
            DataName='{}_{}'.format(V,S)
            com='{}=np.loadtxt("{}")'.format(DataName,FileName)
            exec(com,globals(),Dic)
            
    for V in VarList:
        Data=[]
        for S in Suffixes:
            DataName='{}_{}'.format(V,S)
            Data.append(Dic[DataName])
            
        diff=np.where(abs((Data[0]-Data[1])/Data[0])>Precision)[0].shape
        print('{}:{}'.format(V,diff))
    return Dic

def CompareOMPJacobian(Suffixes:list,Folder:str,VarList:list,Verbose=False,Precision=1e-15):
    Dic={}
    if len(Suffixes)!=2: raise ValueError('Suffixes must be a list of two elements')
    for V in VarList:
        for S in Suffixes:
            FileName=os.path.join(Folder,'{}_{}.txt'.format(V,S))
            DataName='{}_{}'.format(V,S)
            com='{}=np.loadtxt("{}")'.format(DataName,FileName)
            exec(com,globals(),Dic)
            
    for V in VarList:
        Data=[]
        for S in Suffixes:
            DataName='{}_{}'.format(V,S)
            Data.append(Dic[DataName])
            
        diff=np.where(abs((Data[0]-Data[1])/Data[0])>Precision)[0].shape
        print('{}:{}'.format(V,diff))
    return Dic
    
        
    
    