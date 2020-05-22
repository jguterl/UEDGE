#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 09:00:51 2020

@author: jguterl
"""
from colorama import Fore, Back, Style
import numpy as np
import os
Dic={'ncore[0]':np.linspace(1e19,1e20,10),'pcoree':np.linspace(0.5e6,5e6,10),'pcoree':np.linspace(0.5e6,5e6,10)}

Dic={'pcoree':np.linspace(0.5e6,5e6,10),'pcorei':np.linspace(0.5e6,5e6,10)}   
Npoints=10
Data={}
f=lambda x:x**0.1
Data['ncore[0]']={'Min':1e19,'Max':1e20,'f':f}
f=lambda x:x**3
Data['pcoree']={'Min':0.5e6,'Max':5e6,'f':f}
Data['pcorei']=Data['pcoree']
Data=MakeTrajectories(Npoints,Data)
fig,ax=plt.subplots(len(list(Data.keys())))
for (i,(k,v)) in enumerate(Data.items()):
    print(i,k,v)
    ax[i].plot(v,'o',label=k)
plt.show()    

        
            
        
        
    
    
    
    
    
