# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 15:54:25 2016

@author: aholaj
"""
##########################################
###                                    ###
### PURPOSE OF THE SCRIPT              ###
### plot profiles from sound_in files  ###
###                                    ###
##########################################

import numpy as np

import matplotlib.pyplot as plt

folder='/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/'
#file = 'sound_in_DYCOMSIIRF02'
file = 'sound_in_eimear'
filu = folder+file
f = open(filu, 'r')

z = []
t = []
q = []
u = []
v = []


for line in f:
    zA, tA, qA, uA, vA = line.split()
    z.append(float(zA))
    t.append(float(tA))
    q.append(float(qA))
    u.append(float(uA))
    v.append(float(vA))
    
f.close()
z[0]=0.

plt.figure()
plt.plot( t, z )

plt.figure()
plt.plot( q, z )

plt.show()

#print z