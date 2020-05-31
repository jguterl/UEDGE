#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 10:34:18 2020

@author: jguterl
"""
import types
try: 
    import uedge
    from uedge import UEDGEToolBox
except:
    pass
# from uedge import *

# from uedge import UEDGEDoc
from uedge.UEDGEDoc import SearchSilent
from uedge.UEDGEMesh import UEDGEMesh
from uedge.UEDGEPlot import UEDGEPlot
import numpy as np
import os
import scipy
from scipy import interpolate
class UEDGEIOBase():
    """
    Base class for input/output of UEDGE data
    """
    
    
    def __init__(self):
        self.ListPkg=UEDGEToolBox.GetListPackage()
    
    def RemovePkg(self,VarName:str):
            
        if len(VarName.split('.'))>1:
            return VarName.split('.')[1]
        else:
            return VarName
                        
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
    
    def GatherData(self,Mode='regular',ExtraVars=[],VarGlobals=[]):
        VarList=self.SelectVars(Mode,ExtraVars)+VarGlobals
        if self.Verbose: print("GatherData VarList:",VarList)
        ProcessedVarList=self.ProcessVars(VarList)
        if self.Verbose: print("GatherData Processed VarList:",VarList)
        
        
        Data={}    
        for pkg in self.ListPkg:
            exec('from uedge import '+pkg)
        Dic=locals() 
        for V in ProcessedVarList:
            if self.Verbose:print('Preparing data {}'.format(V))
            #remove package name 
            if len(V.split('.'))>0:
                VV=V.split('.')[1]
            else:
                VV=V
            try:    
                exec("Data['{}']={}".format(VV,V),globals(),Dic)
            except:
                print('Unable to process data {}'.format(V))
                
        return Dic['Data']
    
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
 
    def Interpolate2D(self,r,z,data,rnew,znew):
        if r.shape!=data.shape or z.shape!=data.shape:
            print('Mismatch in shape of data and grid:{}:{}/{}'.format(data.shape,r.shape,z.shape))
            return data
        else:
            points = np.array( (r.flatten(), z.flatten()) ).T
            values = data.flatten()
            return interpolate.griddata(points, values, (rnew,znew), method='nearest', rescale=False)
    
    def Interpolate2DData(self,Data,OldGrid,NewGrid,VarList=[]):
        r=OldGrid['rm'][:,:,0]
        z=OldGrid['zm'][:,:,0]
        rnew=NewGrid['rm'][:,:,0]
        znew=NewGrid['zm'][:,:,0]
        if self.Verbose: print('old shape:{}; new shape={}'.format(r,shape,rnew.shape))
        for (K,data) in Data.items():
            if self.Verbose: print('{}.shape={}'.format(K,data.shape))
            if K in VarList:
                if len(data.shape)>2:
                    dataout=np.zeros((rnew.shape[0],rnew.shape[1],data.shape[2]))
                    for i in range(data.shape[2]):
                        out=self.Interpolate2D(r,z,np.squeeze(data[:,:,i]),rnew,znew)
                        dataout[:,:,i]=out[:,:]
                    Data[K]=dataout
                else:
                    Data[K]=self.Interpolate2D(r,z,data,rnew,znew)
                        
        return Data
    
    def LoadDataPackage(self,Data,CheckDim=True,Enforce=True):
        print('Loading the following variables into UEDGE packages:')
        print('Variables:',list(Data.keys()))
        
        for pkg in self.ListPkg:
            exec('from uedge import '+pkg)
            
        for k,v in Data.items():
            Mismatch=False
            DicVar=SearchSilent(k)
            if self.Verbose: print('Loading variable {}'.format(k))
            if len(DicVar)<1:
                print('No package found for {}'.format(k))
            else:    
                if len(DicVar)>1:
                    raise ValueError('More than one package found for {}'.format(k))
                Pkg=DicVar[0]['Package']
                try:
                    if CheckDim:
                        if type(v)==np.ndarray:
                            Dic=locals()
                            exec("D={}.{}".format(Pkg,k),globals(),Dic)
                            if (Dic['D'].shape!=v.shape):
                                        if Enforce:
                                            raise ValueError('Mismatch in dimension of {}'.format(k))
                                        else:
                                            print('Mismatch in dimension of {}'.format(k))
                                            Mismatch=True
                            else:
                                if self.Verbose: print('Matching dimensions for {}'.format(k))
                                            
                    Dic=locals()
                    Dic['v']=v
                    if Mismatch and type(v)==np.ndarray:
                        Str='['+",".join(['0:{}'.format(v.shape[i]) for i in range(len(v.shape))])+"]"
                        if self.Verbose: print('exec to fix mismatch:','{}.{}{}=v'.format(Pkg,k,Str))
                        exec('{}.{}{}=v'.format(Pkg,k,Str),globals(),Dic)    
                    else:
                        if self.Verbose: print('exec:','{}.{}=v'.format(Pkg,k))
                        exec('{}.{}=v'.format(Pkg,k),globals(),Dic)    
                    #Dic['v']=v
                        
                    if self.Verbose: print('-> Success')
                except Exception as e:
                    if self.Verbose: print('-> Failed')
                    if self.Verbose: print(repr(e))
                    if Enforce:
                        if self.Verbose: print(repr(e))    
                        raise ValueError('Cannot set {}'.format(k))
                    else:
                        print('Skipping loading of data for {} '.format(k))
    
    def SelectVars(self,Mode='regular',ExtraVars=[]):
        cls=self.__class__
        if Mode =='regular':
            self.VarList=cls.PlasmaVars+cls.PlasmaVarss+cls.RunVars+cls.GridVars
        elif Mode=='plasma':
            self.VarList=cls.PlasmaVars+cls.PlasmaVarss
        elif Mode=='grid':
            self.VarList=cls.GridVars
        elif Mode=='full':
            pass
        elif Mode=='run':
            self.VarList=cls.GridVars
        elif Mode is None:
            pass
        else:
            raise KeyError('Unknown value for argument mode mode=regular|plasma|grid|full|run')
        self.VarList+=ExtraVars
        if self.Verbose:
           print('Variable List:',self.VarList)
        return self.VarList
    
    
    def SelectWorker(self,Format='numpy'):
       if Format=='numpy' or Format=='npy':
           return Numpy(self.Verbose)
       else:
           raise KeyError('Unknown format. format=numpy|hdf5|json|txt')
       
       if self.Verbose: print('Format:',Format)
        
        
class Numpy(UEDGEIOBase):
    """
    Class for load/save UEDGE data in numpy format
    TODO: implementation of full mode for dumping entire simulation
    """
    def __init__(self,Verbose=False):
        self.Verbose=Verbose
        
    # def _SaveData(self,FileName,VarList=[],Tag:dict={}):
    #     if not FileName.endswith('.npy'):
    #             FileName=FileName+'.npy'
    #     print('Saving data in {} ...'.format(FileName))
    #     for pkg in self.ListPkg:
    #         exec('from uedge import '+pkg)
    #     Header=VarList
    #     if len(Header)>1:
    #         HeaderTag=list(Tag.keys())
    #         DataTag=list(Tag.values())
    #         if len(HeaderTag)<1:
    #             HeaderTag=['dummy']
    #             DataTag=['0']
                
    #         Dic=locals()
    #         if Verbose: 
    #             print('- Header:',Header)
    #             print('- Header Tag:',HeaderTag)
    #             print('- Tag Data:',DataTag)
    #         exec('Data=tuple([{}])'.format(','.join(VarList)),globals(),Dic)
            
    #         exec('TagData=tuple(DataTag)',globals(),Dic)
    #         try: 
    #             np.save(FileName,(Header,HeaderTag,Dic['TagData'],Dic['Data']))
    #             if Verbose:
    #                 print('Numpy: Data saved in file:',FileName)
    #         except Exception as e:
    #             raise IOError('Could not save plasma state in numpy format in {}:'.format(FileName),repr(e))  
            
    #     else:
    #         print('No data to save... No file written')
       
    
    
    def SaveData(self,FileName,DataSave:dict={},Tag:dict={}):
        print('Saving data in {} ...'.format(FileName))
        if not DataSave:
            print('No data to save... Skipping')
            return False
        
        
        if not Tag: # check if dictionaary is empty
            HeaderTag=['dummy']
            DataTag=['0']
            
        HeaderTag=list(Tag.keys())
        DataTag=list(Tag.values())    
        HeaderData=list(DataSave.keys())
        Data=list(DataSave.values())

        if self.Verbose: 
            print('- Header:',HeaderData)
            print('- Header Tag:',HeaderTag)
            print('- Tag Data:',DataTag)
    
        try: 
            np.save(FileName,(HeaderData,HeaderTag,tuple(DataTag),tuple(Data)))
            if self.Verbose:
                print('Numpy: Data saved in file:',FileName)
            return True
        except Exception as e:
            raise IOError('Could not save plasma state in numpy format in {}:'.format(FileName),repr(e))  
        
        else:
            print('No data to save... No file written')
            return False


    # def _LoadData(self,FileName,LoadList=[],ExcludeList=[],Enforce=True,Verbose=False,CheckSize=True,LoadPackage=True):
    #     Out={}
    #     print('Loading data from {} ...'.format(FileName))
    #     Header,HeaderTag,TagData,Data=self.ReadData(FileName,Verbose)
    #     Out=self.ImportData(Header,Data,LoadList,ExcludeList,Enforce,Verbose,CheckSize,LoadPackage)
    #     Out['HeaderTag']=HeaderTag
    #     Out['TagData']=TagData
    #     return Out
    
    def LoadData(self,FileName):
        
        
        print('Loading data from {} ...'.format(FileName))
        
        try: 
            Header,HeaderTag,TagData,Data=np.load(FileName,allow_pickle=True)
        except Exception as e:
            raise IOError('Could not load data from numpy file {}:{}'.format(FileName,repr(e)))
        
        if self.Verbose: print('Data loaded from {}'.format(FileName))
        Dic={}
        Data=dict((K,V) for (K,V) in zip(Header,Data))
        Tag=dict((K,V) for (K,V) in zip(HeaderTag,TagData))
        return (Data,Tag)
        
        
        # Out=self.ImportData(Data,LoadList,ExcludeList,Enforce,Verbose,CheckSize,LoadPackage)
        # Out['HeaderTag']=HeaderTag
        # Out['TagData']=TagData
        # return Out
    
    
    
    
    
    # def LoadInterpolate(self,FileName,OldGrid,NewGrid=None,CaseName=None,Folder='SaveDir',Format='numpy',LoadList=[],ExcludeList=[]):
    #     print('Loading data from {} ...'.format(FileName))
    #     O=self.LoadData(FileName,CaseName,Folder,Format,LoadList,ExcludeList)
    #     if NewGrid=None:
    #         NewGrid=self.GetGrid()
    #     Dat=self.InterpolatePlasmaVars(O.Data)
    #     Header=Dat.keys()
    #     Data=Dat.values()
    #     Out=self.ReadData(Header,Data,LoadList,ExcludeList,Enforce=True,CheckSize,LoadPackage)
    #     Out['HeaderTag']=HeaderTag
    #     Out['TagData']=TagData
    
    # def _ReadData(self,FileName,Verbose=False,Enforce=True):
        
    #     try: 
    #         return np.load(FileName,allow_pickle=True)
    #     except Exception as e:
    #         raise IOError('Could not load data from numpy file {}:{}'.format(FileName,repr(e)))
    #     if Verbose: print('Data loaded from {}'.format(FileName))
        
            
                    
class UEDGEIO(UEDGEIOBase):
    """
    Class to handle save/load of UEDGE data.
    Data dump in file 
    """
    PlasmaVars=['te','ti','phi','ng','up','tg','ni']
    PlasmaVarss=[V+'s' for V in PlasmaVars]
    GridVars=['nisp','ngsp','rm','zm','nx','ny','iysptrx']
    RunVars=['dtreal','dt_tot','ftol_dt','GridFileName']
    
    def __init__(self,Verbose=False):
        self.VarList=[];
        self.Verbose=Verbose
        self.ListPkg=UEDGEToolBox.GetListPackage()
        self.InputLines=[]
        self.CaseName='None'
        
    def ImportData(self,Data:dict,LoadList=[],ExcludeList=[],CheckDim=True,Enforce=True,LoadPackage='plasma'):
        if self.Verbose:
            print('Importing:',Data.keys())
        Data=dict((self.RemovePkg(k),v) for k,v in Data.items())
        if self.Verbose:
            print('Filtered:',Data.keys())    
        if len(LoadList)==0:
            LoadList=list(Data.keys())
        
       
            
        if self.Verbose: print('Variables available:',list(Data.keys()))
        
        DataSelect=dict((k,v) for (k,v) in Data.items() if k in LoadList and k not in ExcludeList)
        if self.Verbose: print('Selected variables:',list(DataSelect.keys()))
        if self.Verbose: print('LoadPackage:',LoadPackage)
        if type(LoadPackage)==str:
            if LoadPackage=='all':
                DataPkg=dict((k,v) for (k,v) in DataSelect.items())
                self.LoadDataPackage(DataPkg,CheckDim,Enforce)
            elif LoadPackage.lower()=='none':
                return DataSelect
            else:    
                LoadPackage=self.SelectVars(LoadPackage)
                DataPkg=dict((k,v) for (k,v)in DataSelect.items() if k in LoadPackage)
                self.LoadDataPackage(DataPkg,CheckDim,Enforce)
        elif type(LoadPackage)==list:
            DataPkg=dict((k,v) for (k,v) in DataSelect.items() if k in LoadPackage)
            self.LoadDataPackage(DataPkg,CheckDim,Enforce)
        else:
            print('LoadPackage must a string or a list of variable names... Import of data aborted.')
            return None
        
        
        return DataSelect
    
    def SaveLog(self,FileName,Str,Tag={}):
        if not os.path.exists(FileName):
            with open(FileName, "a") as F:
                for (k,v) in Tag.items():
                    F.write('### {}:{}\n'.format(k,v))
        with open(FileName, "a") as F:
            F.write(Str+'\n')

    
   
           
    # def _Save(self,FileName,Mode='regular',ExtraVars=[],GlobalVars=[],Tag={},Format='numpy',Verbose=False):
    #     """ Method to save UEDGE data in a file
    #         Parameters:
    #             mode(str): define lists of variables set as attributes of the parent class UEDGEIO to build list of variables to be saved  
    #                      ='plasma': 
    #                      ='grid'
    #                      ='complete'
    #                      ='full'
    #     """
    #     Data=self.ProcessData(Mode,ExtraVars,GlobalVars,Verbose)
        
        
    #     if Format=='numpy' or Format=='npy':
    #         Worker=Numpy(self.Verbose)
    #     else:
    #         raise KeyError('Unknown format. format=numpy|hdf5|json|txt')
        
    #     if Verbose: print('Format:',Format)
        
    #     Worker.Save(FileName,Data,Tag,Verbose)
    
    def Save(self,FileName,Mode='regular',ExtraVars=[],GlobalVars=[],Tag={},Format='numpy'):
        """ Method to save UEDGE data in a file
            Parameters:
                mode(str): define lists of variables set as attributes of the parent class UEDGEIO to build list of variables to be saved  
                         ='plasma': 
                         ='grid'
                         ='complete'
                         ='full'
        """
        Data=self.GatherData(Mode,ExtraVars,GlobalVars)
        Data['_InputLines']=self.InputLines
        Data['_CaseName']=self.CaseName
        Worker=self.SelectWorker(Format)
        Worker.SaveData(FileName,Data,Tag)
    
   
        
    def Load(self,FileName,Format='numpy',LoadList=[],ExcludeList=[],CheckDim=True,Enforce=True,LoadPackage='plasma'):  
        Worker=self.SelectWorker(Format)
        (Data,Tag)=Worker.LoadData(FileName)
        Data=self.ImportData(Data,LoadList,ExcludeList,CheckDim,Enforce,LoadPackage)
        return (Data,Tag)

    def LoadInterpolate(self,FileName,OldGrid,NewGrid,Format='numpy',LoadList=[],ExcludeList=[],CheckDim=True,LoadPackage='plasma',InterpolateData='plasma'):
        D=UEDGEData(FileName,Format='numpy',LoadList=[],ExcludeList=[])
        VarList=self.SelectVars(Mode=InterpolateData)
        if self.Verbose:
            print('Variables to be interpolated:',VarList)
        
        if self.Verbose: 
            print('OldGrid',OldGrid)
        
        if OldGrid=='loaded':
            OldGrid=D.Grid
        elif type(OldGrid)==str:
            OldGrid=UEDGEMesh().ImportMesh(OldGrid)
        
        if type(OldGrid)!=dict:
            print('OldGrid must be of type grid')
        if not OldGrid:
            print('OldGrid is empty')
            return (D.Data,D.Tag)
        D.Data=self.Interpolate2DData(D.Data,OldGrid,NewGrid,VarList)
        Enforce=True
        D.Data=self.ImportData(D.Data,LoadList,ExcludeList,CheckDim,Enforce,LoadPackage)
        
        return (D.Data,D.Tag)    


IO=UEDGEIO()        
         
    #Default lists of variables to be save
    
class UEDGEData(UEDGEPlot,UEDGEIO):
    def __init__(self,FileName,CaseName=None,Folder='SaveDir',Format='numpy',LoadList=[],ExcludeList=[],Verbose=False):
        self.Verbose=Verbose
        self.LoadData(FileName,CaseName,Folder,Format,LoadList,ExcludeList)
        
    def LoadData(self,FileName,CaseName=None,Folder='SaveDir',Format='numpy',LoadList=[],ExcludeList=[]):
        self.FilePath=UEDGEToolBox.Source(FileName,Folder=Folder,Enforce=False,Verbose=True,CaseName=CaseName)
        if self.FilePath is not None:
            (Data,self.Tag)=IO.Load(self.FilePath,Format=Format,LoadList=[],ExcludeList=[],CheckDim=False,LoadPackage='none')
            self.ProcessData(Data)
        else:
            print('Cannot load {}'.format(self.FilePath))
    def ProcessData(self,Data):
        cls=self.__class__
        self.Data=dict((k,Data[k]) for k in Data.keys() if k in cls.PlasmaVars or k in cls.PlasmaVarss)
        self.Grid=dict((k,Data[k]) for k in Data.keys() if k in cls.GridVars)
        
    def GetData(self,Field):
        return self.Data.get(Field.lower())
    
    def SetData(self,Field): #to be compataible with UEDGEsimulation class
        pass
    
    def GetGrid(self):
        return self.Grid
    
    def SetGrid(self,Grid=None):
        if Grid is not None:
            self.Grid=Grid
    
        
    def GetCaseName(self):
        if self.Tag.get('CaseName') is None:

            return 'None'
        else:
            return self.Tag.get('CaseName')
        
    
    
        
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
