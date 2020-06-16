#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:52:13 2020

@author: jguterl
"""

import types
#import matplotlib.pyplot as plt
import numpy as np
import os, sys, math, string, re
from pathlib import Path
from uedge import UEDGEToolBox
from uedge.UEDGEToolBox import *
from uedge.UEDGEDoc import *

from uedge.UEDGEIO import UEDGEIO
from uedge import UEDGEPlot 
from uedge import UEDGEMesh
from uedge.UEDGEPlot import *
from uedge.UEDGEMesh import *
from colorama import Fore, Back, Style
from datetime import date,datetime

class UEDGESimBase(UEDGEMesh):    
    def __init__(self,  Verbose=False,*args, **kwargs):
        # Import Uedge packages as attribute of the class instance 
        # WARNING: this is not a deep copy!!!
        self.InputLines=[]
        self.ListPkg=UEDGEToolBox.GetListPackage()
        self.Verbose=Verbose
        print('Loaded Packages:',self.ListPkg)
        # for pkg in self.ListPkg:
        #     exec('self.' + pkg + '=' + pkg,globals(),locals())
        self.IO=UEDGEIO(self.Verbose)
        self.ExcludeList=['ExcludeList','ListPkg','IO']+self.ListPkg
        #self.SetVersion()
        
    def Updateftol(self):
        bbb.ylodt = bbb.yl
        bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
        fnrm_old=sqrt(sum((bbb.yldot[0:bbb.neq-1]*bbb.sfscal[0:bbb.neq-1])**2))   
        return max(min(bbb.ftol_dt, 0.01*fnrm_old),bbb.ftol_min)
    
    def GetData(self,Field):
        for pkg in self.ListPkg:
            exec('from uedge import '+pkg)
        Out=locals()
        comm='{}={}.{}'.format(Field.lower(),'bbb',Field.lower())
        if self.Verbose: print('Retrieving data:',comm)
        try: 
            exec(comm,globals(),Out)
        except:
            print('Cannot set {}={}.{} '.format(Field,'bbb',Field))
        if Out.get(Field) is None:
            print('Cound not get data for {}'.format(Field))

        return Out.get(Field.lower())
    
    def SetData(self):
        for pkg in self.ListPkg:
            exec('from uedge import '+pkg)
        Out=locals()
        comm='{}={}.{}'.format(Field.lower(),'bbb',Field.lower())
        if self.Verbose: print('Setting data:',comm)
        try: 
            exec(comm,globals(),Out)
        except:
            print('Cannot set {}={}.{} '.format(Field,'bbb',Field))
        if Out.get(Field) is None:
            print('Cound not get data for {}'.format(Field))

        self.Data=Out
    
    def GetGrid(self):
        Grid={'rm':com.rm,'zm':com.zm,'iysptrx':com.iysptrx,'nx':com.nx,'ny':com.ny,'psi':com.psi}
        return Grid
        
    def ReadInput(self,FileName:str=None,Folder:str='InputDir',Verbose:bool=True,Initialize:bool=False,Filter='*'):
        '''
        ReadInput(self,FileName:str,Folder:str=None,Verbose:bool=True,Initialize:bool=True)
        
        Parse and execute each line of FileName.
        FileName must a path toward a python script file (.py)
        The existence of FileName is first checked with its absolute path.
        If Folder is None, the path is the current working directory. 
        If the file is not found, then the file is search for in the 'InputDir' defined in Settings.
        
        Parameters
        ----------
        FileName : str
            Path to input file 
        Folder: str,optional [None]
            
        Verbose : bool, optional
            Verbose mode when parsing the input file. The default is True.
            

        Returns
        -------
        None.

        '''
        if FileName is None:
            FileName=LsFolder(Source(None,Folder),Filter=Filter,Ext='*.py',LoadMode=True)
        if FileName is None:
            print('Nothing to read... Exiting ...')
            return    
        # Looking for file 
        for pkg in self.ListPkg:
            exec('from uedge import '+pkg)
        FilePath=UEDGEToolBox.Source(FileName,Folder=Folder,Enforce=False,Verbose=True)
        if FilePath is None:
            FilePath=UEDGEToolBox.Source(FileName,Folder='InputDir',Enforce=True,Verbose=True)
        print('### Loading {} '.format(FilePath))    
        
        # parsing file    
        f=open(FilePath,'r')
        lines=f.read()
        f.close()
        Lines=lines.splitlines()
        count=1
        for L in Lines:
            if not L.strip().startswith('#'):
                if Verbose:
                    print('{} : {}'.format(count,L))
                self.InputLines.append(L)    
                exec(L)
            count=count+1
        self.CurrentInputFile=FilePath 
        if Initialize:
            self.Initialize()
    def Save(self,FileName,CaseName='auto',Folder='SaveDir',Mode='regular',ExtraVars=[],GlobalVars=[],Tag={},Format=None,ForceOverWrite=False,Verbose=False):
        '''
        Wrapper method to save UEDGE simulation data
        See Save method of UEDGEIO class
        
        Parameters
        ----------
        FileName : str
            Path to input file 
        Verbose : bool, optional
            Verbose mode when parsing the input file. The default is True.

        Returns
        -------
        None.

        '''
        if CaseName=='auto':
            CaseName=self.CaseName
        
            
        if Format is None:
            Format=self.__class__.Format
            
        Verbose=Verbose or self.Verbose
        
        self.Tag=GenerateTag(self)
        self.Tag.update(Tag)
        FilePath=UEDGEToolBox.Source(FileName,Folder=Folder,Enforce=False,Verbose=Verbose,CaseName=CaseName,CheckExistence=False,CreateFolder=True)
        if Verbose:
            print("Saving data in file:{}".format(FilePath))
        # Looking for file 
        if ForceOverWrite or UEDGEToolBox.CheckFileExist(FilePath):
            self.IO.Verbose=Verbose
            self.IO.Save(FilePath,Mode,ExtraVars,GlobalVars,self.Tag,Format)
        else:
            print("No data saved in file:{}".format(FilePath))        
    
    def SaveLog(self,FileName,Str,Tag={},Folder='SaveDir'):    
        self.Tag=GenerateTag(self)
        self.Tag.update(Tag)
        FilePath=UEDGEToolBox.Source(FileName,Folder=Folder,Enforce=False,Verbose=self.Verbose,CaseName=self.CaseName,CheckExistence=False,CreateFolder=True)
        self.IO.SaveLog(FilePath,Str,self.Tag)
            
    def LoadInterp(self,FileName,OldGrid='loaded',NewGrid=None,CaseName=None,Folder='SaveDir',LoadPackage='plasma',InterpolateData='plasma',LoadList=[],ExcludeList=[],Format=None,CheckDim=True,Verbose=False):
        
        if Format is None:
            Format=self.Format
        Verbose=Verbose or self.Verbose 
        
        if NewGrid is None:
            NewGrid=self.GetGrid()
        elif type(NewGrid)==str:
                NewGrid=UEDGEMesh().ImportMesh(NewGrid)
        elif type(NewGrid)!=dict:
            print("NewGrid must a Grid or a path toward a grid file")
             
        FilePath=UEDGEToolBox.Source(FileName,Folder=Folder,Enforce=False,Verbose=Verbose,CaseName=CaseName)
        if Verbose:
            print("Load data in file:{}".format(FilePath))
        # Looking for file 
        if os.path.isfile(FilePath):
            self.IO.Verbose=Verbose
            self.IO.LoadInterpolate(FilePath,OldGrid,NewGrid,Format,LoadList,ExcludeList,CheckDim,LoadPackage,InterpolateData)
            bbb.restart=1
        else:
            print("The file {} does not exist".format(FilePath))
            
            
        
                
    def Load(self,FileName=None,CaseName=None,Folder='SaveDir',LoadPackage='all',LoadList=[],ExcludeList=[],Format=None,CheckDim=True,Enforce=True,Verbose=False,Filter='*'):
        
        """
        Wrapper method to load UEDGE simulation data
        See Load method of UEDGEIO class
        Args:
            FileName (TYPE): DESCRIPTION.
            CaseName (TYPE, optional): DESCRIPTION. Defaults to None.
            Folder (TYPE, optional): DESCRIPTION. Defaults to 'SaveDir'.
            LoadList (TYPE, optional): DESCRIPTION. Defaults to [].
            ExcludeList (TYPE, optional): DESCRIPTION. Defaults to [].
            Format (TYPE, optional): DESCRIPTION. Defaults to 'numpy'.
            CheckCompat (TYPE, optional): DESCRIPTION. Defaults to True.
            Verbose (TYPE, optional): DESCRIPTION. Defaults to False.

        Returns:
            None.

        """
        if FileName is None:
            FileName=LsFolder(Source(None,Folder),Filter=Filter,Ext='*.npy',LoadMode=True)
        if FileName is None:
            print('Nothing to load... Exiting ...')
            return
        if Format is None:
            Format=self.Format
        Verbose=Verbose or self.Verbose 
        
        from uedge import bbb
        if Format=='npy' or Format=='numpy':
            if not FileName.endswith('.npy'):
                FileName=FileName+'.npy'
        
        FilePath=UEDGEToolBox.Source(FileName,Folder=Folder,Enforce=False,Verbose=Verbose,CaseName=CaseName)
        if Verbose:
            print("Load data in file:{}".format(FilePath))
        # Looking for file 
        if os.path.isfile(FilePath):
            self.IO.Verbose=Verbose
            self.IO.Load(FilePath,Format,LoadList,ExcludeList,CheckDim,Enforce,LoadPackage)
            bbb.restart=1
        else:
            print("The file {} does not exist".format(FilePath))
        
            
    @classmethod    
    def GetClassAttr(cls,Verbose=False):
        """
        

        Args:
            cls (TYPE): DESCRIPTION.
            Verbose (TYPE, optional): DESCRIPTION. Defaults to False.

        Returns:
            Attr (TYPE): DESCRIPTION.

        """
        
        Attr = dict((k,v) for k,v in cls.__dict__.items() if '__' not in k and not isinstance(v,types.FunctionType) and not isinstance(v,classmethod))
        if Verbose: print(Attr)
        return Attr
    
    def GetInstanceAttr(self,Verbose=False):
        Attr = dict((k,v) for k,v in self.__dict__.items() if '__' not in k and not isinstance(v,types.FunctionType) and not isinstance(v,classmethod))
        if Verbose: print(Attr)
        return Attr
    
    @classmethod    
    def ShowClass(cls,Verbose=False):
        Attr = dict((k,v) for k,v in cls.__dict__.items() if '__' not in k and not isinstance(v,types.FunctionType) and not isinstance(v,classmethod))
        if Verbose: print(Attr)
        for A,V in Attr.items():
            comm='print("{}","=",)'.format(A,V)
            print('{}={}'.format(A,V))
            
    def Show(self,Verbose=False):
        print('Internal UEDGE parameters:')
        for A,V in self.GetInstanceAttr().items(): 
            Result=SearchSilent(A)
            if len(Result)>0:
                Pkg=Result[0]['Package']
                print(' - {}.{}={}'.format(Pkg,A,V))
            else:
                print(' - XXX.{}={}'.format(A,V))
        print('\n Run parameters:')
        for A,V in self.__class__.GetClassAttr().items(): 
            print(' - {}={}'.format(A,V))
            
    def SetParams(self,**kwargs):
        """
        

        Args:
            Verbose (TYPE, optional): DESCRIPTION. Defaults to False.

        Returns:
            None.

        """
        self.OverrideObjectParams(**kwargs)
        self.OverrideClassParams(**kwargs)
        Params=dict((k,v) for (k,v) in self.GetInstanceAttr().items() if k not in self.ExcludeList)
        Params.update(**kwargs)
        if self.Verbose:
            print('Setting up following params in UEDGE:',list(Params.keys()))
        
        self.SetPackageParams(Params) 
                    
    def SetPackageParams(self,Dic:dict):
        """
        

        Args:
            Verbose (TYPE, optional): DESCRIPTION. Defaults to False.

        Returns:
            None.

        """
        
        for A,V in Dic.items():
            if V is not None:
                Result=SearchSilent(A.split('[')[0])
                if len(Result)>0:
                    Pkg=Result[0]['Package']
                    comm='{}.{}={}'.format(Pkg,A,V)
                    if self.Verbose: print(comm)
                    try: 
                        exec(comm,globals(),globals())
                    except:
                        print('Cannot set {}.{} = {} '.format(Pkg,A,V))
                else:
                    print('No package found for {}'.format(A))
    
    def OverrideObjectParams(self,**kwargs):
        """
        

        Args:
            **kwargs (TYPE): DESCRIPTION.

        Returns:
            None.

        """
        
                
        #Override object attribute 
        Dic=self.GetInstanceAttr()
        for k,v in kwargs.items():
            if k in Dic:
                self.__dict__[k]=v
                
    def OverrideClassParams(self,**kwargs):
        #Override class attribute 
        Dic=self.__class__.GetClassAttr()
        for k,v in kwargs.items():
            if k in Dic:
                self.__class__.__dict__[k]=v
                
    def PrintTimeStepModif(self,i):
        self.PrintInfo('New time-step = {:.4E}'.format(bbb.dtreal),color=Back.MAGENTA)
    
    def PrintCurrentIteration(self,i,j=None):
        if j is None:
            self.PrintInfo('{cn}: Main loop i={i}/{imax}       dtreal={dt:.4E}'.format(cn=self.CaseName,i=i,imax=self.Imax, dt=bbb.dtreal),color=Back.BLUE)
        else:
            self.PrintInfo('{cn}: Subloop   i={i}/{imax} j={j}/{jmax} dtreal={dt:.4E}'.format(cn=self.CaseName,i=i,imax=self.Imax,j=j,jmax=self.Jmax,dt=bbb.dtreal),color=Back.YELLOW)
        
    def TimeEvolution(self):
        """
        Allow the evolution of UEDGE settings in time. Set by SetTimeEvolution(). Print a log of values associated with files

        Returns:
            None.

        """
        pass
        #for k,v in self.TimeParameters:
                
            #bbb.dt_tot
    def Run(self,Verbose=False): 
        """
        

        Args:
            Verbose (TYPE, optional): DESCRIPTION. Defaults to False.

        Returns:
            Status='mainloop'|'tstop'|'dtkill'|'aborted'

        """
        
        bbb.exmain_aborted=0 
       
        while bbb.exmain_aborted==0:
            
            
            self.PrintInfo('----Starting Main Loop ----')
# Main loop----------------------------------------------- 
            for imain in range(self.Imax):
                self.Status='mainloop'
                
                bbb.icntnunk = 0
                
                bbb.ftol = self.Updateftol()
     
                if (bbb.initjac == 0): 
                    bbb.newgeo=0
                    
                self.PrintTimeStepModif(imain) 
                
                
                self.PrintCurrentIteration(imain)
                
                bbb.exmain() # take a single step at the present bbb.dtreal
                sys.stdout.flush()
                
                if bbb.exmain_aborted==1:
                    break
                
                if (bbb.iterm == 1 and bbb.exmain_aborted!=1):
                    self.AutoSave()
                    self.SaveLast() # Save data in file SaveDir/CaseName/last.npy
                    bbb.dt_tot += bbb.dtreal
                    self.dt_tot=bbb.dt_tot
                    self.TimeEvolution()

                    if bbb.dt_tot>=bbb.t_stop:
                            bbb.exmain_aborted=1
                            self.SaveFinalState()
                            self.Status='tstop'
                            return self.Status 
                        
             
# Second loop -----------------------------------------------------------------------------------
                    bbb.icntnunk = 1
                    bbb.isdtsfscal = 0
                    for ii2 in range(self.Jmax): #take ii2max steps at the present time-step
                        if bbb.exmain_aborted==1:
                            break
                        bbb.ftol = self.Updateftol()
                        self.PrintCurrentIteration(imain,ii2)
                        bbb.exmain()
                        sys.stdout.flush()
                        
                        if bbb.iterm == 1 and bbb.exmain_aborted!=1:
                            self.SaveLast() # Save data in file SaveDir/CaseName/last.npy
                            bbb.dt_tot += bbb.dtreal
                            self.dt_tot=bbb.dt_tot
                            self.TimeEvolution()
                        else:
                            break
                        
                        if bbb.dt_tot>=bbb.t_stop:
                            self.PrintInfo('SUCCESS: dt_tot >= t_stop')
                            self.SaveFinalState()
                            self.Status='tstop'
                            return self.Status 
# End Second loop ----------------------------------------------------------------------------------- 
                             
                if bbb.exmain_aborted==1:
                    break
                
# Handle success/error ------------------------------------------------------------------------------       
                if (bbb.iterm == 1):
                    bbb.dtreal *= self.mult_dt_fwd
                    self.dtreal=bbb.dtreal
                    
                else:    #print bad eqn, cut dtreal by 3
                    self.Itrouble()
             
                    self.PrintInfo('Converg. fails for bbb.dtreal; reduce time-step by 3',Back.RED) 
                    bbb.dtreal /= self.mult_dt_bwd
                    self.dtreal=bbb.dtreal
                    bbb.iterm = 1
                    
                    if (bbb.dtreal < bbb.dt_kill):
                        self.PrintInfo('FAILURE: time-step < dt_kill',Back.RED)
                        bbb.exmain_aborted=1
                        self.Status='dtkill'
                        return self.Status 
# End of main loop -------------------------------------------------------- --------------------------------            
        if bbb.exmain_aborted==1: self.Status='aborted'
        return self.Status
                                   
             
    def Initialize(self,ftol=1e20,restart=0,dtreal=1e10):
        bbb.dtreal=dtreal
        bbb.ftol=ftol
        bbb.restart=restart
        if (bbb.iterm == 1 and bbb.ijactot>1):
           self.PrintInfo("Initial successful time-step exists",Back.GREEN)
           return
        else:
           self.PrintInfo("Taking initial step with Jacobian:",Back.CYAN)
           bbb.icntnunk = 0
           bbb.exmain()
           sys.stdout.flush()
        if (bbb.iterm != 1):
            self.PrintInfo("Error: converge an initial time-step first",Back.RED)
            bbb.exmain_aborted=1
        else:
           self.PrintInfo("First initial time-step has converged",Back.GREEN) 
           return
       
    def PrintInfo(self,Str,color=Back.CYAN,Extra=False):
        if Extra: print("*---------------------------------------------------------*")
        print("{color}{}{reset}".format(Str,color=color,reset=Style.RESET_ALL))
        if Extra: print("*---------------------------------------------------------*")
        
    def Itrouble(self):
        ''' Function that displays information on the problematic equation '''
        from numpy import mod,argmax
        from uedge import bbb
        # Set scaling factor
        scalfac = bbb.sfscal
        if (bbb.svrpkg[0].decode('UTF-8').strip() != "nksol"): scalfac = 1/(bbb.yl + 1.e-30)  # for time-dep calc.
    
        # Find the fortran index of the troublemaking equation
        itrouble=argmax(abs(bbb.yldot[:bbb.neq]))+1
        print("** Fortran index of trouble making equation is:")
        print(itrouble)
    
        # Print equation information
        print("** Number of equations solved per cell:")
        print("numvar = {}".format(bbb.numvar))
        print(" ")
        iv_t = mod(itrouble-1,bbb.numvar) + 1 # Use basis indexing for equation number
        print("** Troublemaker equation is:")
        # Verbose troublemaker equation
        if abs(bbb.idxte-itrouble).min()==0:
            print('Electron energy equation: iv_t={}'.format(iv_t))           
        elif abs(bbb.idxti-itrouble).min()==0:
            print('Ion energy equation: iv_t={}'.format(iv_t))   
        elif abs(bbb.idxphi-itrouble).min()==0:
            print('Potential equation: iv_t={}'.format(iv_t))   
        elif abs(bbb.idxu-itrouble).min()==0:
            for species in range(bbb.idxu.shape[2]):
                if abs(bbb.idxu[:,:,species]-itrouble).min()==0:
                    print('Ion momentum equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxn-itrouble).min()==0:
            for species in range(bbb.idxn.shape[2]):
                if abs(bbb.idxn[:,:,species]-itrouble).min()==0:
                    print('Ion density equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxg-itrouble).min()==0:
            for species in range(bbb.idxg.shape[2]):
                if abs(bbb.idxg[:,:,species]-itrouble).min()==0:
                    print('Gas density equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxtg-itrouble).min()==0:
            for species in range(bbb.idxtg.shape[2]):
                if abs(bbb.idxtg[:,:,species]-itrouble).min()==0:
                    print('Gas temperature equation of species {}: iv_t={}'.format(species, iv_t))   
        # Display additional information about troublemaker cell
        print(" ")
        print("** Troublemaker cell (ix,iy) is:")
        print(bbb.igyl[itrouble-1,])
        print(" ")
        print("** Timestep for troublemaker equation:")
        print(bbb.dtuse[itrouble-1])
        print(" ")
        print("** yl for troublemaker equation:")
        print(bbb.yl[itrouble-1])
        print(" ")
        
    def WhichEq(self,itrouble):
        ''' Function that displays information on the problematic equation '''
        from numpy import mod,argmax
        from uedge import bbb
        # Set scaling factor
        scalfac = bbb.sfscal
        if (bbb.svrpkg[0].decode('UTF-8').strip() != "nksol"): scalfac = 1/(bbb.yl + 1.e-30)  # for time-dep calc.
    
        # Find the fortran index of the troublemaking equation
        print("** Fortran index of trouble making equation is:")
        print(itrouble)
    
        # Print equation information
        print("** Number of equations solved per cell:")
        print("numvar = {}".format(bbb.numvar))
        print(" ")
        iv_t = mod(itrouble-1,bbb.numvar) + 1 # Use basis indexing for equation number
        print("** Troublemaker equation is:")
        # Verbose troublemaker equation
        if abs(bbb.idxte-itrouble).min()==0:
            print('Electron energy equation: iv_t={}'.format(iv_t))           
        elif abs(bbb.idxti-itrouble).min()==0:
            print('Ion energy equation: iv_t={}'.format(iv_t))   
        elif abs(bbb.idxphi-itrouble).min()==0:
            print('Potential equation: iv_t={}'.format(iv_t))   
        elif abs(bbb.idxu-itrouble).min()==0:
            for species in range(bbb.idxu.shape[2]):
                if abs(bbb.idxu[:,:,species]-itrouble).min()==0:
                    print('Ion momentum equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxn-itrouble).min()==0:
            for species in range(bbb.idxn.shape[2]):
                if abs(bbb.idxn[:,:,species]-itrouble).min()==0:
                    print('Ion density equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxg-itrouble).min()==0:
            for species in range(bbb.idxg.shape[2]):
                if abs(bbb.idxg[:,:,species]-itrouble).min()==0:
                    print('Gas density equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxtg-itrouble).min()==0:
            for species in range(bbb.idxtg.shape[2]):
                if abs(bbb.idxtg[:,:,species]-itrouble).min()==0:
                    print('Gas temperature equation of species {}: iv_t={}'.format(species, iv_t))   
        # Display additional information about troublemaker cell
        print(" ")
        print("** Troublemaker cell (ix,iy) is:")
        print(bbb.igyl[itrouble-1,])
        print(" ")
        print("** Timestep for troublemaker equation:")
        print(bbb.dtuse[itrouble-1])
        print(" ")
        print("** yl for troublemaker equation:")
        print(bbb.yl[itrouble-1])
        print(" ")   
        
    def SetCaseName(self,CaseName):
        from uedge import bbb
        self.CaseName=CaseName
        try:
            bbb.CaseName=CaseName
            
        except:
            print('Cannot set CaseName in bbb')
            
    def GetCaseName(self):
        try:
            return bbb.CaseName[0].decode().rstrip()
        except:
            return None
        
    def SetDescription(self,Str):
        self.Description=Str
        
    def GetDescription(self):
        return self.Description
        
    def GetVersion(self):
        return uedge.__version__  
    def GetGitHash(self):
        return uedge.GitHash
    def SetVersion(self):
        try:
            bbb.uedge_ver=uedge.__version__
        except:
            print('Cannot set Uedge version')
            
    def SaveLast(self):
        self.Save('last.npy',CaseName=self.CaseName,Folder='SaveDir',Mode=self.__class__.Mode,ExtraVars=[],GlobalVars=[],Tag={},Format=self.__class__.Format,ForceOverWrite=True)
        
    def AutoSave(self):
        self.iSave+=1
        if self.iSave>self.ISave:
            self.iSave=1
            self.Save('save_{}'.format(GetTimeStamp()))
             
                
    def SaveFinalState(self):
        self.Save('final_state_{}'.format(GetTimeStamp()),CaseName=self.CaseName,Folder='SaveDir',Mode=self.__class__.Mode,Format=self.__class__.Format,ForceOverWrite=True)
        self.Save('final_state',CaseName=self.CaseName,Folder='SaveDir',Mode=self.__class__.Mode,Format=self.__class__.Format,ForceOverWrite=True)
        



#---------------------------------------------------------------------------------------------------------------     
class UEDGESim(UEDGESimBase,UEDGEPlot):
    """ Main class to run UEDGE simulation.

    This class contains methods to run, save and restore UEDGE simulation.
    An instance of this class "Sim" is automatically created when the module uedge is imported.
    This class is derived from UEDGESimBase, which containes the method Run() which is the equilent of rdrundt/rdcondt and Show  
    
    Methods:
        
    Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * For module TODOs
    * You have to also use ``sphinx.ext.todo`` extension

.. _Google Python Style Guide:
   http://google.github.io/styleguide/pyguide.html

"""
    mult_dt_fwd=3.4
    mult_dt_bwd=3
    iSave=1
    ISave=10        
    Imax=500
    Jmax=5
    Format='numpy'
    Mode='regular'
    
    
    
    def __init__(self,Verbose=False):
        self.DefaultSettings()
        self.Verbose=Verbose
        UEDGESimBase.__init__(self)
        
    def DefaultSettings(self):
        self.dtreal=1e-8
        self.t_stop=10
        self.ftol_min=1e-10
        self.ftol_dt=1e-10
        self.ftol=1e-8
        self.itermx=7
        self.rlx=0.9
        self.incpset=7
        self.dt_tot=None
        self.isbcwdt=1
        self.CaseName='None'
        self.Description='None'
        
        @staticmethod
        def Resetftol():
            bbb.ftol_min=1e-10
            bbb.ftol_dt=1e-10
            print("bbb.ftol_dt={};bbb.ftol_min={}".format(bbb.ftol_min,bbb.ftol_dt))
        @staticmethod
        def Changeftol(factor):
            bbb.ftol_min*=factor
            bbb.ftol_dt*=factor
            print("bbb.ftol_dt={};bbb.ftol_min={}".format(bbb.ftol_min,bbb.ftol_dt))    
         
    def Start(self,**kwargs):
        self.Initialize()
        self.Cont(**kwargs)
    
    def Cont(self,**kwargs):
        bbb.restart=1
        self.SetParams(**kwargs)
        return self.Run(Verbose=self.Verbose)
    
    def Restart(self,**kwargs):
        self.dtreal/=self.mult_dt_fwd
        self.SetParams(**kwargs)
        bbb.restart=1
        bbb.iterm=1
        return self.Run()
    
    def Help(self):
        help(self.__class__)
             

        
    def Restore(self,FileName=None,Folder='InputDir'):
        """
        
    
        Args:
            FileName (TYPE): DESCRIPTION.
            CaseName (TYPE, optional): DESCRIPTION. Defaults to Sim.CaseName.
            Folder (TYPE, optional): DESCRIPTION. Defaults to 'SaveDir'.
            LoadList (TYPE, optional): DESCRIPTION. Defaults to [].
            ExcludeList (TYPE, optional): DESCRIPTION. Defaults to [].
            Format (TYPE, optional): DESCRIPTION. Defaults to Sim.Format.
            CheckCompat (TYPE, optional): DESCRIPTION. Defaults to True.
            Verbose (TYPE, optional): DESCRIPTION. Defaults to False.
    
        Returns:
            None.
    
        """
        
        Sim.ReadInput(FileName,Folder=Folder,Initialize=True)
        Sim.Load('last.npy',CaseName=Sim.CaseName,LoadPackage='plasma')
        
        
    def RunRamp(self,Data:dict,Istart:int=0,dtreal_start:float=1e-8,tstop:float=10):
        """
        

        Args:
            Data (dict): ebbb.g. Dic={'ncore[0]':np.linspace(1e19,1e20,10),'pcoree':np.linspace(0.5e6,5e6,10)}
            Istart (int, optional): DESCRIPTION. Defaults to 0.
            dtreal_start (float, optional): DESCRIPTION. Defaults to 1e-8.
            tstop (float, optional): DESCRIPTION. Defaults to 10.

        Returns:
            None.

        """
        
        #Check if all data arrays have the same length
        List=[v.shape for (k,v) in Data.items()]
        if not all(L == List[0] for L in List):
            print('Arrays of different size... Cannot proceed...')
            return
        Istop=List[0][0]
        if Istart>=Istop:
            print('Istart={} >= Istop={}: cannot proceed...'.format(Istart,Istop))
        irun=Istart
        # Loop over data
        
        BaseCaseName=self.CaseName
        while irun <Istop:
          
            # 1) Set data in uedge packages
            Params=dict((k,v[irun]) for (k,v) in Data.items())
            ListParams=['{}:{}'.format(k,v) for k,v in Params.items()]
            StrParams=['{}_{:2.2e}'.format(k.split('[')[0],v) for k,v in Params.items()]
            
            self.CaseName=BaseCaseName+'_'.join(StrParams)
            self.SetPackageParams(Params)
            
            # 2) Run until completion
            self.PrintInfo('RAMP i={}/{} : '.format(irun,Istop)+','.join(ListParams),color=Back.MAGENTA)
            Status=self.Cont(dt_tot=0,dtreal=dtreal_start,t_stop=tstop)
            if Status=='tstop':
                ListValueParams=['{:2.2e}'.format(v) for k,v in Params.items()]
                self.Save('final_state_ramp_'+'_'.join(ListValueParams)) 
                self.SaveLog('logramp','{}:{}::'.format(self.Tag['Date'],self.Tag['Time'])+';'.join(ListParams))              
                irun+=1
            else:
                print('Exiting ramp... Need to add a routine to restart after dtkill')
                return

Sim=UEDGESim()

def TogglePlasma():
    bbb.isnion[0]=Toggle(bbb.isnion[0])
    bbb.isupon[0]=Toggle(bbb.isupon[0])
    bbb.istion=Toggle(bbb.istion)
    bbb.isteon=Toggle(bbb.isteon)
    
def TurnOffPlasma():
    bbb.isnion[0]=0
    bbb.isupon[0]=0
    bbb.istion=0
    bbb.isteon=0
    
def TurnOnPlasma():
    bbb.isnion[0]=1
    bbb.isupon[0]=1
    bbb.istion=1
    bbb.isteon=1    
    
def Toggle(State):
    if State ==0:
        return 1
    elif State == 1:
        return 0
    else:
        raise ValueError('Cannot toggle request variable in state:{}'.format(State))
        
def ReadInput(*args,**kwargs):
    Sim.ReadInput(*args,**kwargs)

def Read(*args,**kwargs):
    Sim.ReadInput(*args,**kwargs)     

    
def Load(*args,**kwargs):
    Sim.Load(*args,**kwargs)  
    
def Save(*args,**kwargs):
    Sim.Save(*args,**kwargs)
    
def Restore(*args,**kwargs):
    Sim.Restore(*args,**kwargs)    

def SetCaseName(CaseName:str):
    """
    Set the name of the simulation case.

    Args:
        CaseName (str): Name of the simulation.

    Returns:
        None.

    """
    
    Sim.SetCaseName(CaseName)
     


def GetCaseName():
    try:
        print("# CaseName:{}".format(Sim.GetCaseName()[0]))
    except:
        print('Cannot get CaseName')    
        pass
         
#---------------------------------------------------------------------------------------------------------------        
class UEDGEDivertorPlates():
    def ReadDivertorPlateFile(self,FileName,Verbose=False):
    
        Dic={'rplate1':[],'zplate1':[],'rplate2':[],'zplate2':[],'Comment':'No comment','FileName':FileName,'InnerDivPlateFile':None,
        'OuterDivPlateFile':None}
        
        
        try:
            with  open(os.path.join(FileName)) as f:
                exec(f.read(),Dic)
        except:
            raise ValueError('Cannot Read File:'+FileName)
        
        
        
        if Dic['InnerDivPlateFile'] is not None:
             InnerData=np.loadtxt(Dic['InnerDivPlateFile'])
             if Verbose: print(InnerData)   
             Dic['rplate1']=list(InnerData[:,0])
             Dic['zplate1']=list(InnerData[:,1])
        
        if Dic['OuterDivPlateFile'] is not None:
            OuterData=np.loadtxt(Dic['OuterDivPlateFile'])  
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
        

eV=1.602176634e-19     
#%% 
def GenerateTag(Sim=None):
    from uedge import __version__,GitHash
    today = date.today()
    now = datetime.now()
    Tag={}
    Tag['Date'] = today.strftime("%d%b%Y")
    Tag['Time'] = now.strftime("%H-%M-%S")
    Tag['User']=    Settings.UserName
    Tag['Version'] = __version__
    Tag['PlatForm'] = Settings.Platform
    Tag['GitHash']=GitHash
    Tag['InputDir']=Settings.InputDir
    Tag['RunDir']=Settings.RunDir
    Tag['SaveDir']=Settings.InputDir
    Tag['PythonVersion']=Settings.Platform['python_version']
    Tag['Machine']=Settings.Platform['machine']
    Tag['Processor']=Settings.Platform['processor']
    Tag['PlatForm']=Settings.Platform['platform']
    Tag['CaseName']='None'
    Tag['Description']='None'
    try:
        if Sim.CaseName is not None:
            Tag['CaseName']=Sim.CaseName
        if Sim.Description is not None: 
            Tag['Description']=Sim.Description
    except: 
        pass
    return Tag

def GetTimeStamp():
    today = date.today()
    now = datetime.now()
    return "{}_{}".format(today.strftime('%d%m%y'),now.strftime("%H%M%S"))

    def RunRamp(self,Data:dict,Istart=0,dtreal_start=1e-8,tstop=10):
        """
        
    
        Returns:
            None.
    
        """
        #Check if all data arrays have the same length
        List=[v.shape for (k,v) in Data.items()]
        if not all(L == List[0] for L in List):
            print('Arrays of different size... Cannot proceed...')
            return
        Istop=List[0][0]
        if Istart>=Istop:
            print('Istart={} >= Istop={}: cannot proceed...'.format(Istart,Istop))
        irun=Istart
        # Loop over data
        
        
        while irun <Istop:
    
            # 1) Set data in uedge packages
            Params=dict((k,v[irun]) for (k,v) in Data.items())
            ListParams=['{}:{}'.format(k,v) for k,v in Params.items()]
            self.SetPackageParams(Params)
            
            # 2) Run until completion
            self.PrintInfo('RAMP i={}/{} : '.format(irun,Istop)+','.join(ListParams),color=Back.MAGENTA)
            Status=self.Cont(dt_tot=0,dtreal=dtreal_start,t_stop=tstop)
            if Status=='tstop':
                ListValueParams=['{:2.2e}'.format(v) for k,v in Params.items()]
                self.Save('final_state_ramp_'+'_'.join(ListValueParams)) 
                self.SaveLog('logramp','{}:{}::'.format(self.Tag['Date'],self.Tag['Time'])+';'.join(ListParams))              
                irun+=1
            else:
                print('Exiting ramp... Need to add a routine to restart after dtkill')
                return