#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 22:09:43 2020

@author: jguterl
"""

import sys,os
import easygui
import platform,inspect
import configparser
import getpass
#from .UEDGEToolBox import LsFolder

import numpy as np
class UEDGESettings():
    def __init__(self):
        
        self.Platform=self.GetPlatform()
        self.config={'UEDGE':{'UserName':None,'RunDir':None,'SaveDir':None,'InputDir':None}}
        self.config=self.CheckUEDGEConfig()
        if self.config is None:
            print('Using default settings. Run InitConfig() to set the configuration...')
            self.UserName=getpass.getuser()
            self.RunDir=os.getcwd()
            self.SaveDir=os.getcwd()
            self.InputDir=os.getcwd() 
        else:
            print('# Loading UEDGE settings from {}'.format(self.ConfigFileName))
            self.UserName=self.config['UEDGE'].get('UserName')
            self.RunDir=self.config['UEDGE'].get('RunDir')
            self.SaveDir=self.config['UEDGE'].get('SaveDir')
            self.InputDir=self.config['UEDGE'].get('InputDir') 

    def Config(self):
        '''Print current UEDGE configuration'''
        print('******** UEDGE configuration ********')
        print('Username:',self.UserName)
        print('RunDir:',self.RunDir)
        print('SaveDir:',self.SaveDir)
        print('InputDir:',self.InputDir)
    def GetPlatform(self):
        PF={}
        AllFunctions = inspect.getmembers(platform, inspect.isfunction)
        for (n,f) in AllFunctions:
            if not n.startswith('_'):
                try: 
                    PF[n]=f()
                except:
                    pass
        return PF
    
    def CheckUEDGEConfig(self):
        self.ConfigFileName=os.path.join(os.path.expanduser("~"),'.UedgeSettings')
        if not os.path.exists(self.ConfigFileName):
            print ('Cannot find the config file ' + self.ConfigFileName +'....')
            return None
        config = configparser.ConfigParser()
        try:
            config.read(self.ConfigFileName)
            return config
        except: 
            print('Cannot parse the config file:',self.ConfigFileName)
            return None

        
def CreateUEDGESettingsFile():
        config = configparser.ConfigParser()
        GetInput=True
        while GetInput:
            UserName=input('Enter the username:')
            
            RunDir=easygui.diropenbox(title='Select RunDir for pyUEDGE',default='.')
            SaveDir=easygui.diropenbox(title='Select SaveDir for pyUEDGE',default='.')
            InputDir=easygui.diropenbox(title='Select InputDir for pyUEDGE',default='.')
            if RunDir is None:
                RunDir=os.getcwd()
            if SaveDir is None:
                RunDir=os.getcwd()
            if InputDir is None:
                InputDir=os.getcwd()
            if UserName is None:
                UserName='unknown'
            RunDir=os.path.os.path.abspath(RunDir)
            SaveDir=os.path.os.path.abspath(SaveDir)
            InputDir=os.path.os.path.abspath(InputDir)
            print(' # Settings:')
            print('  - UserName:',UserName)
            print('  - RunDir:',RunDir)
            print('  - SaveDir:',SaveDir)
            print('  - InputDir:',InputDir)
            if QueryYesNo('Are these settings correct?'):
               config['UEDGE']={}
               config['UEDGE']['UserName']=UserName
               config['UEDGE']['RunDir']=RunDir
               config['UEDGE']['SaveDir']=SaveDir
               config['UEDGE']['InputDir']=InputDir
               ConfigFileName=os.path.join(os.path.expanduser("~"),'.UedgeSettings')
               with open(ConfigFileName, 'w') as f:
                   config.write(f)
               print(' # Creation of the config file:'+ConfigFileName)
               GetInput=False
               return
            elif not QueryYesNo('Enter the settings again?'):
               print(' # The config file .UEDGEInfo has not been created in the home folder...') 
               GetInput=False
               return 





def QueryYesNo(question, default="no"):
        """Ask a yes/no question via input() and return their answer.
    
        "question" is a string that is presented to the user.
        "default" is the presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).
    
        The "answer" return value is True for "yes" or False for "no".
        """
        valid = {"yes": True, "y": True, "ye": True,
                 "no": False, "n": False}
        if default is None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)
    
        while True:
            sys.stdout.write(question + prompt)
            choice = input().lower()
            if default is not None and choice == '':
                return valid[default]
            elif choice in valid:
                return valid[choice]
            else:
                sys.stdout.write("Please respond with 'yes' or 'no' "
                                 "(or 'y' or 'n').\n")  
                
Settings=UEDGESettings()
S=Settings

def CdRunDir():
    os.chdir(Settings.RunDir)
    print('Current working directory:',os.getcwd())
 
def CdSaveDir():
    os.chdir(Settings.SaveDir)
    print('Current working directory:',os.getcwd())
    
def CdInputDir():
    os.chdir(Settings.InputDir)
    print('Current working directory:',os.getcwd())
    
def LsFolder(Folder,Filter='*',Ext="*.py",LoadMode=False):
    import glob
    if '.' in Ext:
        ListFile = [f for f in glob.glob(os.path.join(Folder,Ext))] +[f for f in glob.glob(os.path.join(Folder,Filter)) if os.path.isdir(f)]
    else:
        ListFile = [f for f in glob.glob(os.path.join(Folder,Ext))]
    Listfile=list(dict.fromkeys(ListFile).keys())
    ListFile.sort(key=str.casefold)
    print('### Content matching "{}" in {}:'.format(Ext,Folder))
    if ListFile is not None:
        for i,F in enumerate(ListFile):
            print(' [{}]: {}'.format(i,os.path.basename(F))) 
    print('')
    Message='Enter a number to look into a folder or file or press r (return) or q (exit)\n >>>: '
    Input=input(Message)
    while Input!='q' and Input!='r':
        if Input.isnumeric() and ListFile is not None and int(Input) in range(len(ListFile)):
            if os.path.isfile(ListFile[int(Input)]):
                print('File:{}'.format(ListFile[int(Input)]))
                if LoadMode:
                    return ListFile[int(Input)]
                Input=input(Message)
            elif os.path.isdir(ListFile[int(Input)]):
                Input=LsFolder(os.path.join(Folder,ListFile[int(Input)]),Filter,Ext,LoadMode)
                if Input=='r':
                   for i,F in enumerate(ListFile):
                       print(' [{}]: {}'.format(i,os.path.basename(F))) 
                   print('') 
                   Input=input(Message) 
                elif LoadMode  and Input!='q':
                    return Input
        else:
            Input=input(Message)
    if Input=='q' and LoadMode:
        return None
    else:
        return Input
        
    
        
def LsInputDir(Folder=None,Ext='*.py'):
    if Folder is None:
        LsFolder(Settings.InputDir,Ext)
    else:
        LsFolder(os.path.join(Settings.InputDir,Folder),Ext)
        
def LsSaveDir(Folder=None,Ext='*'):
    if Folder is None:
        LsFolder(Settings.SaveDir,Ext)
    else:
        LsFolder(os.path.join(Settings.SaveDir,Folder),Ext)
    
def Config():
    Settings.Config() 
    
def InitConfig():
    CreateUEDGESettingsFile()  
# start into the run folder
CdRunDir()    