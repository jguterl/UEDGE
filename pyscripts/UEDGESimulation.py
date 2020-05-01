#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:52:13 2020

@author: jguterl
"""
import uedge
import types
#import matplotlib.pyplot as plt
import numpy as np
import os, sys, math, string, re
from pathlib import Path
from . import UEDGEToolBox
from .UEDGEToolBox import *
from .UEDGEDoc import *
from . import UEDGEIO
from uedge import *
from colorama import Fore, Back, Style

class UEDGESimBase():    
    def __init__(self,  *args, **kwargs):
        # Import Uedge packages as attribute of the class instance 
        # WARNING: this is not a deep copy!!!
        self.ListPkg=UEDGEToolBox.GetListPackage()
        print('Loaded Packages:',self.ListPkg)
        for pkg in self.ListPkg:
            exec('self.' + pkg + '=' + pkg,globals(),locals())
        self.IO=UEDGEIO.UEDGEIO()
        self.ExcludeList=['ExcludeList','ListPkg','IO']+self.ListPkg
        #self.SetVersion()

    def ReadInput(self,FileName:str,Folder:str=None,Verbose:bool=True):
        '''
        Parse and execute each line of FileName.
        FileName must a path toward a python script file (.py)
        The existence of FileName is first checked with its absolute path.
        If no path is specified, the path is the current working directory. 
        If the file is not found, then the file is search for in the 'InputDir' defined in Settings.
        
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
                exec(L)
            count=count+1
        self.CurrentInputFile=FilePath 
        
    def Save(self,FileName,CaseName=None,Folder='SaveDir',Mode='regular',ExtraVars=[],GlobalVars=[],Tag={},Format='numpy',ForceOverWrite=False,Verbose=False):
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
       
        FilePath=UEDGEToolBox.Source(FileName,Folder=Folder,Enforce=False,Verbose=Verbose,CaseFolder=CaseName,CheckExistence=False)
        if Verbose:
            print("Saving data in file:{}".format(FilePath))
        # Looking for file 
        if UEDGEToolBox.CheckFileExist(FilePath) or ForceOverWrite:
            self.IO.Save(FilePath,Mode,ExtraVars,GlobalVars,Tag,Format,Verbose)
        else:
            print("No data saved in file:{}".format(FilePath))
            
    def Load(self,FileName,CaseName=None,Folder='SaveDir',LoadList=[],ExcludeList=[],Format='numpy',CheckCompat=True,Verbose=False):
        '''
        Wrapper method to load UEDGE simulation data
        See Load method of UEDGEIO class
        
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
        #bbb.allocate() #< no allocation leads to segmentation fault
        if Format=='npy' or Format=='numpy':
            if not FileName.endswith('.npy'):
                FileName=FileName+'.npy'
        
        FilePath=UEDGEToolBox.Source(FileName,Folder=Folder,Enforce=False,Verbose=Verbose,CaseFolder=CaseName)
        if Verbose:
            print("Load data in file:{}".format(FilePath))
        # Looking for file 
        if os.path.isfile(FilePath):
            self.IO.Load(FilePath,Format,LoadList,ExcludeList,CheckCompat=CheckCompat,Verbose=Verbose)
        else:
            print("The file {} does not exist".format(FilePath))
        
            
    @classmethod    
    def GetClassAttr(cls,Verbose=False):
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
            
    def SetUEDGEParams(self,Verbose=False):
        
        for A,V in self.GetInstanceAttr().items():
            if V is not None and V in not ExcludeList:
                Result=SearchSilent(A)
                if len(Result)>0:
                    Pkg=Result[0]['Package']
                    comm='{}.{}={}'.format(Pkg,A,V)
                    if Verbose: print(comm)
                    try: 
                        exec(comm,globals(),globals())
                    except:
                        print('Cannot set {}.{} = {} '.format(Pkg,A,V))
                else:
                    print('No package found for {}'.format(A))
    
    def OverrideParams(self,**kwargs):
        #Override class attribute 
        Dic=self.__class__.GetClassAttr()
        for k,v in kwargs.items():
            if k in Dic:
                self.__class__.__dict__.[k]=v
                
        #Override object attribute 
        Dic=self.GetInstanceAttr()
        for k,v in kwargs.items():
            if k in Dic:
                self.__dict__[k]=v
    
    def PrintTimeStepModif(self,i):
        self.PrintInfo('New time-step = {:.4E}'.format(bbb.dtreal),color=Back.MAGENTA)
    
    def PrintCurrentIteration(self,i,j=None):
        if j is None:
            self.PrintInfo('Main loop i={i}/{imax}       dtreal={dt:.4E}'.format(i=i,imax=self.Imax, dt=bbb.dtreal),color=Back.BLUE)
        else:
            self.PrintInfo('Subloop   i={i}/{imax} j={j}/{jmax} dtreal={dt:.4E}'.format(i=i,imax=self.Imax,j=j,jmax=self.Jmax,dt=bbb.dtreal),color=Back.YELLOW)
        
    
    def Run(self,Verbose=True):        
        bbb.exmain_aborted=0 
        if Verbose: print('00:exmain_aborted:',bbb.exmain_aborted)
        while bbb.exmain_aborted==0:
            if Verbose: print('0:exmain_aborted:',bbb.exmain_aborted)
            if bbb.exmain_aborted==1:
                break
            self.PrintInfo('----Starting Main Loop ----')
            for imain in range(self.Imax):
                
                bbb.icntnunk = 0
                bbb.ylodt = bbb.yl #this is use in set_dt which provides option to select optimal timestep
                bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq]*bbb.sfscal[0:bbb.neq])**2))
                bbb.ftol = max(min(bbb.ftol_dt, 0.01*fnrm_old),bbb.ftol_min)
     
                if (bbb.initjac == 0): 
                    bbb.newgeo=0
                self.PrintTimeStepModif(imain) 
                #self.PrintInfo('Suggested timestep:{}'.format(self.GetMinTimeStep()))
                self.PrintCurrentIteration(imain)
                bbb.exmain() # take a single step at the present bbb.dtreal
                if bbb.exmain_aborted==1:
                    break
                
                if (bbb.iterm == 1):
                    bbb.ylodt = bbb.yl
                    bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                    fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq-1]*bbb.sfscal[0:bbb.neq-1])**2))
                    if bbb.dt_tot>=bbb.t_stop:
                            bbb.exmain_aborted=1
                            self.PrintInfo('SUCCESS: dt_tot >= t_stop')
                            if Verbose: print('11:exmain_aborted:',bbb.exmain_aborted)  
                            break
                    bbb.icntnunk = 1
                    bbb.isdtsfscal = 0
# Second loop -----------------------------------------------------------------------------------                    
                    for ii2 in range(self.Jmax): #take ii2max steps at the present time-step
                        if bbb.exmain_aborted==1:
                            break
                        bbb.ftol = max(min(bbb.ftol_dt, 0.01*fnrm_old),bbb.ftol_min)
                        self.PrintCurrentIteration(imain,ii2)
                        bbb.exmain()
                        if bbb.exmain_aborted==1:
                            break
                        if bbb.iterm == 1:
                            bbb.ylodt = bbb.yl
                            bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                            fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq-1]*bbb.sfscal[0:bbb.neq-1])**2))
                            bbb.dt_tot += bbb.dtreal
                            self.dt_tot=bbb.dt_tot
                            
                        if bbb.dt_tot>=bbb.t_stop:
                            bbb.exmain_aborted=1
                            self.PrintInfo('SUCCESS: dt_tot >= t_stop')
                            if Verbose: print('11:exmain_aborted:',bbb.exmain_aborted)  
                            break
# Second loop -----------------------------------------------------------------------------------                              
                if bbb.exmain_aborted==1:
                    break
# Handle success/error -----------------------------------------------------------------------------------        
                if (bbb.iterm == 1):
                    bbb.dtreal *= self.mult_dt_fwd
                    self.dtreal=bbb.dtreal
                
                else:    #print bad eqn, cut dtreal by 3
                    self.Itrouble()
                    
                    if (bbb.dtreal < bbb.dt_kill):
                        self.PrintInfo('FAILURE: time-step < dt_kill',Back.RED)
                        bbb.exmain_aborted=1
                        break
                    self.PrintInfo('Converg. fails for bbb.dtreal; reduce time-step by 3',Back.RED) 
                    bbb.dtreal /= self.mult_dt_bwd
                    self.dtreal=bbb.dtreal
                              
                bbb.iterm = 1
             
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
        try:
            bbb.CaseName=CaseName
            self.CaseName=CaseName
        except:
            print('Cannot set CaseName')
            
    def GetCaseName(self):
        try:
            return bbb.CaseName[0].decode().rstrip()
        except:
            return None
        
        
    def SetVersion(self):
        try:
            bbb.uedge_ver=uedge.__version__
        except:
            print('Cannot set Uedge version')
        
    


#---------------------------------------------------------------------------------------------------------------     
class UEDGESim(UEDGESimBase):
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
    Imax=500
    Jmax=5
    Format='numpy'
    Mode='regular'
    
    
    
    def __init__(self):
        self.DefaultSettings()
        UEDGESimBase.__init__(self)
        
    def DefaultSettings(self):
        self.dtreal=1e-10
        self.t_stop=10
        self.ftol_min=1e-10
        self.ftol_dt=1e-10
        self.ftol=1e-8
        self.itermx=7
        self.rlx=0.9
        self.incpset=7
        self.dt_tot=None
        self.iscontdt=1
        self.CaseName=None
        
    
    def Start(self,Verbose=False,**kwargs):
        self.Initialize()
        self.Cont(Verbose=Verbose,**kwargs)
    
    def Cont(self,Verbose=False,**kwargs):
        bbb.restart=1
        self.OverrideParams(**kwargs)
        self.SetUEDGEParams()
        self.Run(Verbose=Verbose)
    
    def Restart(self,**kwargs):
        bbb.restart=1
        self.OverrideParams(**kwargs)
        
        self.dtreal/=self.mult_dt_fwd
        self.SetUEDGEParams()
        bbb.iterm=1
        self.Run()
    
    def Help(self):
        help(self.__class__)
        

Sim=UEDGESim()
def ReadInput(FileName,Folder=None,Verbose=False):
     Sim.ReadInput(FileName,Folder,Verbose=Verbose)   
     
def Restore(FileName,CaseName=Sim.CaseName,Folder='SaveDir',LoadList=[],ExcludeList=[],Format=Sim.Format,CheckCompat=True,Verbose=False):
    Sim.Load(FileName,CaseName,Folder,LoadList,ExcludeList,Format,CheckCompat,Verbose)
    
def Save(FileName,CaseName=Sim.CaseName,Folder='SaveDir',Mode=Sim.Mode,ExtraVars=[],GlobalVars=[],Tag={},Format=Sim.Format,ForceOverWrite=False,Verbose=False):
    Sim.Save(FileName,CaseName,Folder,Mode,ExtraVars,GlobalVars,Tag,Format,ForceOverWrite,Verbose)
    
def SetCaseName(CaseName:str):
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
        
# def InitRun(self):
    #         if (bbb.iterm == 1 and bbb.ijactot>1):
    #            self.PrintInfo("Initial successful time-step exists",Back.GREEN)
    #            return
    #         else:
    #            self.PrintInfo("Need to take initial step with Jacobian; trying to do here",Back.CYAN)
    #            bbb.icntnunk = 0
    #            bbb.exmain()
               
    #         if (bbb.iterm != 1):
    #             self.PrintInfo("Error: converge an initial time-step first",Back.RED)
    #             bbb.exmain_aborted=1
    #         else:
    #            self.PrintInfo("First initial time-step has converged",Back.GREEN) 
    #            return
    
    # def Restart(self,dt=None,mul_dt=3):
    #     bbb.iterm=1
    #     bbb.dtreal/=mult_dt
    #     self.RunTime(dt,mul_dt)
              
    # def GetMinTimeStep(self):
        #evaluate yldot
        # bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
        # dt=np.zeros(bbb.neq) 
        # for i in range(bbb.neq):
        #     if bbb.iseqalg[i]==1:
        #         print('eq al:',i)
        #         dt[i]=1e20
        #     else:
        #         if bbb.yldot[i]==0.0:
        #             dt[i]=1e20
        #         else:    
        #             dt[i]=abs(bbb.yl[i]/(bbb.yldot[i]*bbb.sfscal[0:bbb.neq]))
        # return dt.min()
    #    return -1

    # def RunTime(self,dt=None,tstop=10,Imax=500,Jmax=5,ftol_dt=1e-10,itermx=7,rlx=0.9,incpset=7,method_dt=0,mult_dt_fwd=3.4,mult_dt_bwd=3):
    #     # this allow run restart 
    #     bbb.rlx=rlx
    #     bbb.ftol_dt=ftol_dt
    #     bbb.itermx=itermx
    #     bbb.incpset=incpset
    #     if dt is not None:
    #         bbb.dtreal=dt
    #     bbb.exmain_aborted=0    
    #     while bbb.exmain_aborted==0:
    #         self.InitRun()   
    #         if bbb.exmain_aborted==1:
    #             break
    #         self.PrintInfo('----Starting Main Loop ----')
    #         for imain in range(Imax):
    #             bbb.icntnunk = 0
                
    #             bbb.ylodt = bbb.yl #this is use in set_dt which provides option to select optimal timestep
    #             bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
    #             fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq]*bbb.sfscal[0:bbb.neq])**2))
    #             bbb.ftol = max(min(bbb.ftol_dt, 0.01*fnrm_old),bbb.ftol_min)
                
               
                
    #             if (bbb.initjac == 0): 
    #                 bbb.newgeo=0
    #             self.PrintInfo('Number time-step changes = {} New time-step = {colorT}{:.4E}{reset}'.format(imain, bbb.dtreal,reset=Style.RESET_ALL,colorT=Back.MAGENTA)) 
    #             self.PrintInfo('Suggested timestep:{}'.format(self.GetMinTimeStep()))
    #             bbb.exmain() # take a single step at the present bbb.dtreal
    #             if bbb.exmain_aborted==1:
    #                 break
                
    #             if (bbb.iterm == 1):
    #                 bbb.ylodt = bbb.yl
    #                 bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
    #                 fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq-1]*bbb.sfscal[0:bbb.neq-1])**2))
    #                 if (bbb.dt_tot>=bbb.t_stop):
    #                     self.PrintInfo('SUCCESS: dt_tot >= t_stop')
                        
    #                     break
    #                 bbb.icntnunk = 1
    #                 bbb.isdtsfscal = 0
    #                 for ii2 in range( 1, bbb.ii2max+1): #take ii2max steps at the present time-step
    #                     if bbb.exmain_aborted==1:
    #                         break
    #                     bbb.itermx = itermx
    #                     bbb.ftol = max(min(bbb.ftol_dt, 0.01*fnrm_old),bbb.ftol_min)
    #                     bbb.exmain()
    #                     if (bbb.dt_tot>=bbb.t_stop):
    #                         self.PrintInfo('SUCCESS: dt_tot >= t_stop')
    #                         bbb.exmain_aborted=1
    #                         break
    #                     if bbb.exmain_aborted==1:
    #                         break
    #                     if bbb.iterm == 1:
    #                         bbb.ylodt = bbb.yl
    #                         bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
    #                         fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq-1]*bbb.sfscal[0:bbb.neq-1])**2))
    #                         bbb.dt_tot += bbb.dtreal
    #                     else:
    #                         break
                        
                        
    #             if bbb.exmain_aborted==1:
    #                         break
                        
    #             if (bbb.iterm != 1):    #print bad eqn, cut dtreal by 3, set irev flag
    #                 #print('\n\n\n exmain_aborted:',bbb.exmain_aborted,'\n\n\n')
    #                 self.Itrouble()
                    
    #                 if (bbb.dtreal < bbb.dt_kill):
    #                     self.PrintInfo('FAILURE: time-step < dt_kill',Back.RED)
    #                     break
    #                 self.PrintInfo('Converg. fails for bbb.dtreal; reduce time-step by 3',Back.RED) 
    #                 bbb.dtreal /= mult_dt_bwd
                    
                    
    #             if (bbb.iterm == 1):
    #                 bbb.dtreal *= mult_dt_fwd
                    
    #             bbb.iterm = 1    
                   


# def RunTime(*arg,**kwargs):
#     Sim.RunTime(*arg,**kwargs)


    
# def ir():
#     Sim.InitRun()
        
     
