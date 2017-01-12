# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 16:48:41 2016

@author: aholaj
"""

import numpy as np

import matplotlib.pyplot as plt

folder='/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/designs/'
#file = 'sound_in_DYCOMSIIRF02'
file = 'corr_design_15d.csv'
filu = folder+file
f = open(filu, 'r')

case     = []
q_inv    = []
tpot_inv = []
q_pbl    = []
tpot_pbl = []
pblh     = []
num_pbl  = []

i=0
for line in f:
    #print line
        
    caseA, q_invA, tpot_invA, q_pblA, tpot_pblA, pblhA, num_pblA = line.split(',')
    if ( i > 0 ):
        case.append( caseA.replace('"',"").zfill(2) )
        q_inv.append( float(q_invA) )
        tpot_inv.append( float(tpot_invA) )
        q_pbl.append( float(q_pblA) )
        tpot_pbl.append( float(tpot_pblA) )
        pblh.append( float(pblhA) )
        num_pbl.append( float(num_pblA))
        
    i=i+1   
    
    
f.close()

print pblh