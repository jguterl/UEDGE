#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:46:09 2020

@author: jguterl
"""
import uedge
import os
from uedge import *
from .UEDGESettings import *

def GetListPackage()->list:
    import pkgutil
    import uedge
    ListPkgUEDGE=[]
    package = uedge
    PkgList=list(pkgutil.iter_modules(package.__path__))
    for pkg in PkgList:
        PkgName=pkg.name
        if PkgName.endswith('py'):
            ListPkgUEDGE.append(PkgName[:-2])
   
    return ListPkgUEDGE

def CheckFileExist(FilePath:str)->bool:
    import os
    if os.path.isfile(FilePath):
        out=' '
        while out!='y' and out!='n' and out!='': 
            out=input('The file {} already exists. Do you want to overwrite it? [y]/n'.format(FilePath))
            if out=='y' or out=='':
                return True
            elif out=='n':
                return False            
    else:
        return True        

    ListPkgUEDGE=[]
    package = uedge
    PkgList=list(pkgutil.iter_modules(package.__path__))
    for pkg in PkgList:
        PkgName=pkg.name
        if PkgName.endswith('py'):
            ListPkgUEDGE.append(PkgName[:-2])
   
    return 

def Source(ObjectName:str,Folder:str='InputDir',Enforce=True,Verbose:bool=False,CaseName=None,CheckExistence=True,CreateFolder=False):
    """
    

    Args:
        ObjectName (str): DESCRIPTION.
        Folder (str, optional): DESCRIPTION. Defaults to 'InputDir'.
        Enforce (TYPE, optional): DESCRIPTION. Defaults to True.
        Verbose (bool, optional): DESCRIPTION. Defaults to False.
        CaseName (TYPE, optional): DESCRIPTION. Defaults to None.
        CheckExistence (TYPE, optional): DESCRIPTION. Defaults to True.

    Raises:
        IOError: DESCRIPTION.

    Returns:
        ObjectPath (TYPE): DESCRIPTION.

    """

    if Verbose:
        print('# Looking for input file {} in {}'.format(ObjectName,Folder))
    if Folder=='InputDir':
        try:
            ObjectDir=Settings.InputDir
        except: 
            print('# Settings object for UEDGE not find... Looking for InputDir in current directory')
            ObjectDir='InputDir'        
    elif Folder=='RunDir':
        try:
            ObjectDir=Settings.RunDir
        except: 
            print('# Settings object for UEDGE not find... Looking for RunDir in current directory')
            ObjectDir='RunDir'
    elif Folder=='SaveDir':
        try:
            ObjectDir=Settings.SaveDir
        except: 
            print('# Settings object for UEDGE not find... Looking for SaveDir in current directory')
            ObjectDir='SaveDir'
    elif Folder is None:
        ObjectDir=None
    else:    
        ObjectDir=Folder
        
    
    if ObjectDir is None:    
        ObjectPath=os.path.abspath(ObjectName)
    else:
        if CaseName is  not None:
            ObjectDir=os.path.join(ObjectDir,CaseName)
            if not os.path.isdir(ObjectDir) and CreateFolder:
                try:
                    os.mkdir(ObjectDir)
                except OSError:
                    pass
                    #print ("Creation of the directory {} failed".format(ObjectDir))
        ObjectPath=os.path.join(os.path.abspath(ObjectDir),ObjectName)
        
    if CheckExistence and not os.path.exists(ObjectPath):
        if Enforce:
            raise IOError('Cannot find {}:'.format(ObjectPath))
        else:
            print('### Cannot find {}'.format(ObjectPath))
        return None
    else:
        if Verbose:
            print('### Found {}'.format(ObjectPath))
        return ObjectPath 
        
        
def PrintSummary():
    print('******** Summary **********')
    print('* neq  :',bbb.neq)
    print('* nisp :',com.nisp)
    print('* nhsp :',com.nhsp)
    print('* ngsp :',com.ngsp)
    print('* nzspt:',com.nzspt)
    print('* nusp :',com.nusp)
    print('* iigsp:',bbb.iigsp)
    print('***************************')
    
def MakeTrajectories(Npoints,Data):
        for k,v in Data.items():
            Min=v['Min'] 
            Max=v['Max']
            x=np.linspace(0,1,Npoints)
            f=v['f']
            Data[k]=Min+f(x)*(Max-Min)
        return Data

    
        
# def GetCaseFolder(CaseName:str):
#     if CaseName=='Default':
#         CaseName=UEDGESimulation.GetTag()
    
#         elif         
# def
#  PathFile(FileName:str,SaveFolder:str='SaveDir',CaseName='Default',Enforce=True,Verbose:bool=False):
    
#         CaseFolder=GetCaseFolder(CaseFolder)
            
#          if os.path.dirname(FileName)!=''
#         if Verbose:
#             print('### Looking for input file {} in {}'.format(FileName,Folder))
#         if Folder=='SaveDir':
#             try:
#                 FileDir=Settings.InputDir
#             except: 
#                 print('Settings object for UEDGE not find... Looking for SaveDir in current directory')
#                 FileDir='SaveDir'
#         elif Folder is None:
#             ObjectDir=None
#         else:
#             ObjectDir=Folder
            
#         if ObjectDir is None:    
#             FilePath=os.path.abspath(FileName)
#         else:
#             FilePath=os.path.join(os.path.abspath(ObjectDir),FileName)
            
#         if not os.path.exists(FilePath):
#             if Enforce:
#                 raise IOError('Cannot find {}:'.format(FiletPath))
#             else:
#                 print('### Cannot find {}'.format(ObjectPath))
#             return None
#         else:
#             if Verbose:
#                 print('### Found {}'.format(ObjectPath))
#             return ObjectPath
