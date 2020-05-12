#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 10:45:03 2020

@author: jguterl
"""

import numpy as np
s=np.loadtxt('serialjac.dat')
p=np.loadtxt('paralleljac.dat')
for i in range(s.shape[0]):
    for j in range(4):
        if s[i,j]!=p[i,j]:
            print('>>>>',i, s[i,:])
            print('    ',i, p[i,:])
        
