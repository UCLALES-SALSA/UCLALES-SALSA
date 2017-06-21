#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:33:06 2017

@author: aholaj
"""
# hakee log -tiedoston perusteella emulattoorin kokonaissuoritusajan
import os
import sys
#f = open(, 'r')

from PythonMethods import Muunnos

if ( len(sys.argv) > 1):
    filu = sys.argv[1]
    
else:
    filu = '/home/aholaj/koodit_IL/bash/emulatortime/log'

folder = os.path.dirname(os.path.realpath( filu ))
f = open(filu, 'r+')

unixtimes=[]
for line in f:
    if (line[0] != 'L' and line[0] !=' ' and line[0] !='S'):
        unixtimes.append(int(line))

aloitusaika  = min( unixtimes )
lopetusaika  = max( unixtimes )
suoritusaika = lopetusaika - aloitusaika
print 'Suoritusaika ' +  Muunnos(suoritusaika)+'\n' 
f.write( 'Suoritusaika ' +  Muunnos(suoritusaika)+'\n' )




f.close()
