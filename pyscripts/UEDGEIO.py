#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 10:34:18 2020

@author: jguterl
"""
import types
# import uedge
# from uedge import *
from . import UEDGEToolBox
# from uedge import UEDGEDoc
from .UEDGEDoc import SearchSilent

import numpy as np


class UEDGEIOBase():
    """
    Base class for input/output of UEDGE data
    """
    
    
    def __init__(self):
        self.ListPkg=UEDGEToolBox.GetListPackage()
         
    def ProcessVars(self,VarList):
        Out=[]
        for V in VarList:
                Result=SearchSilent(V)
                if len(Result)>0:
                    Pkg=Result[0]['Package']
                    Out.append('{}.{}'.format(Pkg,V))
                else:
                    print('No package found for {}'.format(A))
        return list(dict.fromkeys(Out))
    
    def ShowVars(self):
        print(self.VarList)
        
    def GetDefaultVarList(self):
        Out=[]
        for k,v in self.__class__.GetClassAttr().items():
            Out.extend(v)
        return  Out
    
    @classmethod    
    def GetClassAttr(cls,Verbose=False):
        Attr = dict((k,v) for k,v in cls.__dict__.items() if '__' not in k and not isinstance(v,types.FunctionType) and not isinstance(v,classmethod))
        if Verbose: print(Attr)
        return Attr
        

        
    
            
    
class Numpy(UEDGEIOBase):
    """
    Class for load/save UEDGE data in numpy format
    TODO: implementation of full mode for dumping entire simulation
    """
    
    def Save(self,FileName,VarList=[],Tag:dict={},Verbose=False):
        for pkg in self.ListPkg:
            exec('from uedge import '+pkg)
        Header=VarList
        if len(Header)>1:
            HeaderTag=list(Tag.keys())
            TagData=list(Tag.values())
            if len(HeaderTag)<1:
                HeaderTag=['dummy']
                TagData=['0']
                
            D=locals()
            if Verbose: 
                print('Numpy:',Header)
                print('Numpy:',HeaderTag)
                print('Numpy:',TagData)
            exec('Data=tuple([{}])'.format(','.join(VarList)),globals(),D)
            exec('TagData=tuple([{}])'.format(','.join(TagData)),globals(),D)
            try: 
                np.save(FileName,(Header,HeaderTag,D['TagData'],D['Data']))
                if Verbose:
                    print('Numpy: Data saved in file:',FileName)
            except Exception as e:
                raise IOError('Could not save plasma state in numpy format in {}:'.format(FileName),repr(e))  
            
        else:
            print('No data to save... No file written')

    def Load(self,FileName,LoadList=[],ExcludeList=[],Enforce=True,Verbose=False,CheckSize=True):
        if not FileName.endswith('.npy'):
                FileName=FileName+'.npy'
        Header,HeaderTag,TagData,Data=self.LoadData(FileName,Verbose)
        if Verbose: print('Reading data in...')
        self.ReadData(Header,Data,LoadList,ExcludeList,Enforce=True,Verbose=Verbose,CheckSize=CheckSize)
        
    def LoadData(self,FileName,Verbose=False,Enforce=True):
        
        try: 
            return np.load(FileName,allow_pickle=True)
        except Exception as e:
            raise IOError('Could not load data from numpy file {}:{}'.format(FileName,repr(e)))
        if Verbose: print('Data loaded from {}'.format(FileName))
        
    def ReadData(self,VarList,Data,LoadList=[],ExcludeList=[],CheckSize=True,Verbose=False,Enforce=True):
        if Verbose: print('0: Reading data in')
        for pkg in self.ListPkg:
            exec('from uedge import '+pkg)
        if Verbose: print('1: Reading data in')
        if len(LoadList)==0:
            LoadAll=True
        else:
            LoadAll=False
        if Verbose: print(VarList)
        for V,D in zip(VarList,Data):
            if Verbose: print('V= {}'.format(V))
            if (LoadAll or V in LoadList) and V not in ExcludeList:
                if Verbose: print('Loading {}'.format(V))
                try:
                    #exec('=D'.format(D),globals(),globals())
                    if CheckSize:
                        if V.count('.')>0 :
                            Dic=locals()
                            exec("Condition=type({})==np.ndarray".format(V),globals(),Dic)
                            if Dic['Condition']:
                                exec("Condition2=({}.shape==D.shape)".format(V),globals(),Dic)
                                if not Dic['Condition2']:
                                    if Enforce:
                                        raise ValueError('Mismatch in dimension of {}'.format(V))
                                    else:
                                        print('Mismatch in dimension of {}'.format(V))
                                        
                    D=locals()        
                    exec('{}=D'.format(V),globals(),D)
                    if Verbose: print('-> Success')
                except Exception as e:
                    if Verbose: print('-> Failed')
                    if Enforce:
                        print(repr(e))    
                        raise ValueError('Cannot execute {}={}'.format(V,D))
                    else:
                        print('Skipping loading of data for {} '.format(V))
                    
class UEDGEIO(UEDGEIOBase):
    """
    Class to handle save/load of UEDGE data.
    Data dump in file 
    """
    PlasmaVars=['te','ti','phi','ng','up','tg','ni']
    PlasmaVarss=[V+'s' for V in PlasmaVars]
    GridVars=['nisp','ngsp','rm','zm','nx','ny']
    RunVars=['dtreal','dt_tot','ftol_dt','GridFileName']
    
    def __init__(self):
        self.VarList=[];
        
    
    def SelectVars(self,mode='regular',ExtraVars=[],VarGlobals=[],Verbose=False):
        cls=self.__class__
        if mode =='regular':
            self.VarList=cls.PlasmaVars+cls.PlasmaVarss+cls.RunVars
        elif mode=='plasma':
            self.VarList=cls.PlasmaVars+cls.PlasmaVarss
        elif mode=='grid':
            self.VarList=cls.GridVars
        elif mode=='full':
            pass
        elif mode=='run':
            self.VarList=cls.GridVars
        elif mode is None:
            pass
        else:
            raise KeyError('Unknown value for argument mode mode=regular|plasma|grid|full|run')
        self.VarList+=ExtraVars
        ProcessedVarList=self.ProcessVars(self.VarList)  
        
        ProcessedVarList+=VarGlobals
        if Verbose:
           print('Processed Variable List:',ProcessedVarList)
        return ProcessedVarList
           
    def Save(self,FileName,Mode='regular',ExtraVars=[],GlobalVars=[],Tag={},Format='numpy',Verbose=False):
        """ Method to save UEDGE data in a file
            Parameters:
                mode(str): define lists of variables set as attributes of the parent class UEDGEIO to build list of variables to be saved  
                         ='plasma': 
                         ='grid'
                         ='complete'
                         ='full'
        """
        ProcessedVarList=self.SelectVars(Mode,ExtraVars,GlobalVars,Verbose)
        if Format=='numpy' or Format=='npy':
            Worker=Numpy()
        else:
            raise KeyError('Unknown format. format=numpy|hdf5|json|txt')
        if Verbose: print('Format:',Format)
        Worker.Save(FileName,ProcessedVarList,Tag,Verbose)

    def Load(self,FileName,Format='numpy',LoadList=[],ExcludeList=[],CheckCompat=True,Verbose=False):
        if Format=='numpy' or Format=='npy':
            Worker=Numpy()
        else:
            raise KeyError('Unknown format. format=numpy|hdf5|json|txt')
        Worker.Load(FileName,LoadList,ExcludeList,CheckCompat,Verbose)


IO=UEDGEIO()        
        
    #Default lists of variables to be save
    
    
# N=Numpy()
                    
# S=UEDGEIO()
# bbb.dtreal=55
# S.Save('/home/jguterl/test2.npy',Verbose=True)
# print('bbb.dtreal=',bbb.dtreal)

# #N.Save('/home/jguterl/test2.npy',['bbb.dtreal'],{'v':'1'},Verbose=True)
# bbb.dtreal=54
# print('bbb.dtreal=',bbb.dtreal)
# S.Load('/home/jguterl/test2.npy')

# #S.LoadSave('/home/jguterl/test.npy',Verbose=True)
# print('bbb.dtreal=',bbb.dtreal)
