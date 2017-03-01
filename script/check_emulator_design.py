#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 16:48:41 2016

@author: aholaj
"""

import numpy as np

import matplotlib.pyplot as plt

import ECLAIR_calcs
import sys

if ( len(sys.argv) > 1):
    filu = sys.argv[1]
else:
    folder='/home/aholaj/mounttauskansiot/ibrixmount/DESIGN/'
    #file = 'sound_in_DYCOMSIIRF02'
    file = 'corr_design.csv'
    filu = folder+file
    
f = open(filu, 'r')

LWC = False

case     = []
q_inv    = []
tpot_inv = []
tpot_pbl = []
pblh     = []
num_pbl  = []

if LWC:
    lwc_max = []       
else:
    q_pbl   = []

i=0
for line in f:
    #print line

    if LWC:
        caseA, q_invA, tpot_invA, lwc_maxA, tpot_pblA, pblhA, num_pblA = line.split(',')
    else:
        caseA, q_invA, tpot_invA, q_pblA,   tpot_pblA, pblhA, num_pblA = line.split(',')
    if ( i > 0 ):
        case.append( caseA.replace('"',"").zfill(2) )
        q_inv.append( float(q_invA) )
        tpot_inv.append( float(tpot_invA) )


        tpot_pbl.append( float(tpot_pblA) )
        pblh.append( float(pblhA)*1000. )
        num_pbl.append( float(num_pblA))

        if LWC:
            lwc_max.append( float(lwc_maxA) )        
        else:
            q_pbl.append( float(q_pblA) )      
        
    i=i+1   
        
f.close()




design = len(pblh)
cloudbase   = np.zeros(design)
if LWC:
    q_pbl   = np.zeros(design)       
else:    
    lwc_max = np.zeros(design)       

p_surf = 101780.
checkoutALA = np.zeros((len(case), 6))
checkoutYLA = np.zeros((len(case), 6))

print ' '

for i in xrange(design):
    checkoutALA[i, 0] = i+1
    checkoutYLA[i, 0] = i+1
    if LWC:
        q_pbl[i]     = ECLAIR_calcs.solve_rw( p_surf, tpot_pbl[i], lwc_max[i]*0.001, pblh[i] )
        cloudbase[i] = ECLAIR_calcs.calc_cloud_base( p_surf, tpot_pbl[i], q_pbl[i]  )
    else:
        lwc_max[i]   = ECLAIR_calcs.calc_lwc_altitude( p_surf, tpot_pbl[i], q_pbl[i], pblh[i] )
        cloudbase[i] = ECLAIR_calcs.calc_cloud_base(   p_surf, tpot_pbl[i], q_pbl[i]*0.001  )
    
    
    c = 0    
    #################################
    ###                           ###
    ### CONSTRAINT 1)             ###
    ### q_pbl > q_inv> 1e-3 kg/kg ###
    ###                           ###
    #################################
    c += 1
    if ( q_inv[i]*0.001 < 1e-3 ):
        print 'ÄLÄHDYS', case[i], 'rajoite', c
        print 'q_inv liian pieni', q_inv[i]*0.001
        checkoutALA[i, c] += 1
    if ( q_pbl[i] < q_inv[i]*0.001  ):
        print 'ÄLÄHDYS', case[i], 'rajoite', c
        print 'q_inv liian suuri'
        checkoutYLA[i, c] += 1


    #################################
    ###                           ###
    ### CONSTRAINT 2)             ###
    ### 1 K < t_inv < 15 K        ###
    ###                           ###
    #################################
    c += 1
    if ( tpot_inv[i] < 1. ):
        print 'ÄLÄHDYS', case[i], 'rajoite', c
        print 'tpot_inv liian pieni', tpot_inv[i]
        checkoutALA[i, c] += 1
    if ( tpot_inv[i] > 15. ):
        print 'ÄLÄHDYS', case[i], 'rajoite', c
        print 'tpot_inv liian suuri', tpot_inv[i]
        checkoutYLA[i, c] += 1


    #################################
    ###                           ###
    ### CONSTRAINT 3)             ###
    ### 30m + 50m < pblh < 3000 m ###
    ###                           ###
    #################################
    c += 1
    if ( pblh[i]*1000. < 80. ):
        print 'ÄLÄHDYS', case[i], 'rajoite', c
        print 'pblh liian pieni', pblh[i]
        checkoutALA[i, c] += 1
    if ( pblh[i] > 3000. ):
        print 'ÄLÄHDYS', case[i], 'rajoite', c
        print 'pbl liian suuri', pblh[i]
        checkoutYLA[i, c] += 1


    ####################################
    ###                              ###
    ### CONSTRAINT 4)                ###
    ### 30m < cloud base< pblh - 50m ###
    ###                              ###
    ####################################
    c += 1
    if ( cloudbase[i] < 30. ):
        print 'ÄLÄHDYS', case[i], 'rajoite', c
        print 'cloudbase liian pieni', cloudbase[i]
        checkoutALA[i, c] += 1
    if ( cloudbase[i] > pblh[i]-50. ):
        print 'ÄLÄHDYS', case[i], 'rajoite', c
        print 'cloudbase liian suuri', cloudbase[i], 'ylaraja', pblh[i]-50.
        checkoutYLA[i, c] += 1

    ####################################
    ###                              ###
    ### CONSTRAINT UNNAMED 5)        ###
    ### lwc_max > 0.0                ###
    ###                              ###
    ####################################
    c += 1
    if ( lwc_max[i] < 0.0 ):
        print 'ÄLÄHDYS', case[i], 'rajoite', c
        print 'lwc_max liian pieni', lwc_max
        checkoutALA[i, c] += 1
#    if ( lwc_max > 1000000000.0 ):
#        print 'ÄLÄHDYS', case[i]
#        print c, 'lwc_max liian suuri'
#        k += 1
    
    if ( sum(checkoutALA[i,1:])> 0 or sum(checkoutYLA[i,1:])> 0 ):
        print ' '
print ' '
print 'alarajat'
print checkoutALA
print ' '
print 'ylarajat'
print checkoutYLA
ala=sum(checkoutALA[:,1:])
yla=sum(checkoutYLA[:,1:])
print 'rajoiterikkomukset ALA', ala
print 'rajoiterikkomukset YLA', yla
print 'väärien määrä', sum(ala+yla)