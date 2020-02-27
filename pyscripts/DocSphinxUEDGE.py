#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 23:43:21 2020

@author: jguterl
"""
import os
from uedge import *
Folder='/home/jguterl/UEDGE/docs/source/'
#for pkg in Doc.DocPkg.keys():
pkg='bbb'
with open(os.path.join(Folder,pkg+'Doc.rst'),'w') as f:
          #f.write('class {}Doc(object):'.format(pkg))
          f.write('Package bbb \n')
          f.write('============\n')
          f.write('')

          for V,D in Doc.DocPkg[pkg].items():
              f.write('.. py:attribute:: '+V+'\n')
              f.write('   \n')
              S=D['Comment'].split('\n')
              
              for s in S:
                  f.write('   {}\n'.format(s))
              f.write('   \n')    
              f.write('   :Default: {}\n'.format(D['Default']))
              f.write('   :Dimension: {}\n'.format(D['Dimension']))
              f.write('   :Group: {}\n'.format(D['Group']))
              f.write('   :Type: {}\n'.format(D['Type']))
              f.write('   :Unit: {}\n'.format(D['Unit']))