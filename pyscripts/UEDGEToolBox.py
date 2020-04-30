#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:46:09 2020

@author: jguterl
"""
import uedge
import os
from uedge import *

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

def Source(ObjectName:str,Folder:str='InputDir',Enforce=True,Verbose:bool=False):
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
        elif Folder is None:
            ObjectDir=None
        else:    
            ObjectDir=Folder
            
        
        if ObjectDir is None:    
            ObjectPath=os.path.abspath(ObjectName)
        else:
            ObjectPath=os.path.join(os.path.abspath(ObjectDir),ObjectName)
            
        if not os.path.exists(ObjectPath):
            if Enforce:
                raise IOError('Cannot find {}:'.format(ObjectPath))
            else:
                print('### Cannot find {}'.format(ObjectPath))
            return None
        else:
            if Verbose:
                print('### Found {}'.format(ObjectPath))
            return ObjectPath 
        
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
