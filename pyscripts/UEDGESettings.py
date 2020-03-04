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


import numpy as np
class UEDGESettings():
    def __init__(self):
        print(' # Loading UEDGE settings...')
        self.Platform=self.GetPlatform()
        try: 
            self.config=self.CheckUEDGEConfig()
        except Exception as e:
            print(repr(e))
        finally:
            UserName=self.config['UEDGE'].get('UserName')
            if UserName is None:
               self.UserName='unknown'
            else:
                self.UserName=UserName
            
            RunDir=self.config['UEDGE'].get('RunDir')    
            if RunDir is None:
               self.RunDir=os.getcwd()
            else:
                self.RunDir=RunDir
            
            SaveDir=self.config['UEDGE'].get('SaveDir')    
            if SaveDir is None:
               self.SaveDir=os.getcwd()
            else:
                self.SaveDir=SaveDir  
                
            InputDir=self.config['UEDGE'].get('InputDir')    
            if InputDir is None:
               self.InputDir=os.getcwd()
            else:
                self.InputDir=InputDir  
   
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
        FileName=os.path.join(os.path.expanduser("~"),'.UedgeSettings')
        if not os.path.exists(FileName):
            raise IOError('Cannot find the file: ' + FileName +'. Run CreateUEDGESettingsFile first....')
        config = configparser.ConfigParser()
        try:
            config.read(FileName)
            return config
        except: 
            raise IOError('Cannot read the config file:',FileName)
            
    def CreateStamp(self): 
            from uedge import bbb
            Stamp={}
            Stamp['time'] = time.time()
            Stamp['ctime'] = time.ctime()
            Stamp['User']=self.UserName
            Stamp['Version'] = bbb.uedge_ver
            Stamp['PlatForm'] = self.Platform 
        
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
               FileName=os.path.join(os.path.expanduser("~"),'.UedgeSettings')
               with open(FileName, 'w') as f:
                   config.write(f)
               print(' # Creation of the config file:'+FileName)
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
# start into the run folder
CdRunDir()    