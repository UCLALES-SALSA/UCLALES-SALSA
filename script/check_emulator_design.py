#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 16:48:41 2016

@author: aholaj
"""

import numpy as np

#try:
#	PLOT = True
#	import matplotlib.pyplot as plt
#	print "matplotlib.pyplot imported"
#except ImportError:
#	PLOT = False

try:
	COLOR = True
	from termcolor import colored
	print 'colored imported'
except ImportError: 
	COLOR = False

import ECLAIR_calcs
import sys
import os
import glob
import subprocess
from itertools import cycle
from netCDF4 import Dataset

ibrix = os.environ["IBRIXMOUNT"]
if ( len(sys.argv) > 1):
    filu = sys.argv[1]
    folder = os.path.dirname(os.path.realpath( filu ))
else:
    folder=ibrix+'/DESIGN/'
    #file = 'sound_in_DYCOMSIIRF02'
    cwd = os.getcwd()
    os.chdir(folder)
    for file in glob.glob("*.csv"):
        designbasename=file
    filu = folder+designbasename
    os.chdir(cwd)
print filu

ncfolder = ibrix + '/DESIGNnetcdf/'
    
f = open(filu, 'r')

LWC = True
riviLKM = subprocess.check_output( "cat " + filu + " | wc -l ", shell=True)
nroCases = int( riviLKM )-1
etunolla = len(str(nroCases))
print nroCases
design     = np.zeros( ( nroCases,6 ) )
caselist   = np.chararray(nroCases, itemsize = etunolla)

q_inv    = np.zeros( nroCases )
tpot_inv = np.zeros( nroCases )
clw_max  = np.zeros( nroCases )
tpot_pbl = np.zeros( nroCases )
pblh     = np.zeros( nroCases )
num_pbl  = np.zeros( nroCases )
q_pbl    = np.zeros( nroCases )

printtaus = []
       #1 2 3 4 5 6 7  8    9 10
tabs = [7,6,9,8,10,8,9, 10, 8, 6]
listtabs = tabs
Ltabs = len(tabs)
Stabs = sum(tabs)
if LWC:
    tabs[9] = 6
else:
    tabs[9] = 8 
    
tabs = cycle(tabs)



## READ CSV
i=0
for line in f:
    
    if ( i > 0 ):
        
        A0, A1, A2, A3, A4, A5, A6 = line.split(',')
        if ( i > 0 ):
            index = i-1 
            caselist[index]    = A0.replace('"',"").zfill(etunolla) # case
            design[ index, 0 ] = float( A1 ) # q_inv
            design[ index, 1 ] = float( A2 ) # tpot_inv
            design[ index, 2 ] = float( A3 ) # clw_max
            design[ index, 3 ] = float( A4 ) # tpot_pbl
            design[ index, 4 ] = float( A5 ) # pblh
            design[ index, 5 ] = float( A6 ) # num_pbl
        
        
    else:
#         names = line.split(',')
#         temp = ''
#         for k in xrange(len(names)):
#             if k == 0:
#                 temp = temp + "case".rjust(next(tabs))
#             else:
#                 temp = temp + names[k].replace('"', "").rjust(next(tabs))
#         temp = temp + 'cloudbase'.rjust(next(tabs)) + 'thick'.rjust( next(tabs) )
#         if LWC:
#             temp = temp + 'q_pbl'.rjust(next(tabs))
#         else:
#             temp = temp + 'clw_max'.rjust(next(tabs))
#         temp = temp.replace('\n', '').replace('\r','')
         temp = "case".rjust(listtabs[0])+"q_inv".rjust(listtabs[1])+"tpot_inv".rjust(listtabs[2])+"clw_max".rjust(listtabs[3])+"tpot_pbl".rjust(listtabs[4])+"pblh".rjust(listtabs[5])+"num_pbl".rjust(listtabs[6])+"cloudbase".rjust(listtabs[7])+"thick".rjust(listtabs[8])+"q_pbl".rjust(listtabs[9])

         printtaus.append(temp)
        
    i=i+1
        
f.close()

os.chdir(folder)
tag =subprocess.check_output("git describe --tags | tr -dc '[:alnum:].'", shell=True)   


p_surf = 101780.
cloudbase   = np.zeros(nroCases)

 

for i in xrange(nroCases):
    q_inv[i]    = design[ i, 0 ]
    tpot_inv[i] = design[ i, 1 ]
    clw_max[i]  = design[ i, 2 ]
    tpot_pbl[i] = design[ i, 3 ]
    pblh[i]     = design[ i, 4 ]
    num_pbl[i]  = design[ i, 5 ]

    if LWC:
#                                               Pa    K            g/kg -> kg/kg      m         
        q_pbl[i]     = ECLAIR_calcs.solve_rw( p_surf, tpot_pbl[i], clw_max[i]*0.001, pblh[i] ) 
#                                                    Pa      K            kg/kg from previous line        
        cloudbase[i] = ECLAIR_calcs.calc_cloud_base( p_surf, tpot_pbl[i], q_pbl[i]  ) # m
        
        q_pbl[i] = q_pbl[i] * 1000. # kg/kg -> g/kg
        
    else:
#                                                       Pa      K           g/kg -> kg/kg   m        
        clw_max[i]   = ECLAIR_calcs.calc_lwc_altitude( p_surf, tpot_pbl[i], q_pbl[i]*0.001, pblh[i] )*1000. # kg/kg -> g/kg
#                                                       Pa      K            g/kg -> kg/kg        
        cloudbase[i] = ECLAIR_calcs.calc_cloud_base(   p_surf, tpot_pbl[i], q_pbl[i]*0.001  ) # m


ncfile = Dataset( ncfolder + 'design_'+tag + '.nc', 'w' )

ncfile.createDimension('case', nroCases )
q_inv_ncf     = ncfile.createVariable( 'q_inv',     np.dtype('float32').char, ('case') )
tpot_inv_ncf  = ncfile.createVariable( 'tpot_inv',  np.dtype('float32').char, ('case') )
clw_max_ncf   = ncfile.createVariable( 'clw_max',   np.dtype('float32').char, ('case') )
tpot_pbl_ncf  = ncfile.createVariable( 'tpot_pbl',  np.dtype('float32').char, ('case') )
pblh_ncf      = ncfile.createVariable( 'pblh',      np.dtype('float32').char, ('case') )
num_pbl_ncf   = ncfile.createVariable( 'num_pbl',   np.dtype('float32').char, ('case') )
q_pbl_ncf     = ncfile.createVariable( 'q_pbl',     np.dtype('float32').char, ('case') )
cloudbase_ncf = ncfile.createVariable( 'cloudbase', np.dtype('float32').char, ('case') )
thickness_ncf = ncfile.createVariable( 'thickness', np.dtype('float32').char, ('case') )

q_inv_ncf[:]     = q_inv
tpot_inv_ncf[:]  = tpot_inv
clw_max_ncf[:]   = clw_max
tpot_pbl_ncf[:]  = tpot_pbl
pblh_ncf[:]      = pblh
num_pbl_ncf[:]   = num_pbl
q_pbl_ncf[:]     = q_pbl
cloudbase_ncf[:] = cloudbase
thickness_ncf[:] = pblh - cloudbase

ncfile.close()

def check_constrain( variable, lowerbound, upperbound, variablename, lowerboundNAME, upperboundNAME, unit, dimensions = 90 ):
    if check_constrain.counter >= 1:
        check_constrain.checkoutALA = np.column_stack(( check_constrain.checkoutALA, np.zeros(( dimensions, 1)) ))
        check_constrain.checkoutYLA = np.column_stack(( check_constrain.checkoutYLA, np.zeros(( dimensions, 1)) ))
    else:
        check_constrain.checkoutALA = np.column_stack(( np.arange(1,dimensions+1),   np.zeros(( dimensions, 1)) ))
        check_constrain.checkoutYLA = np.column_stack(( np.arange(1,dimensions+1),   np.zeros(( dimensions, 1)) ))
    check_constrain.counter += 1
    wi = int( np.log(max( np.max(variable), np.max(np.abs(lowerbound)), np.max(np.abs(upperbound)) )) ) + 2
    if isinstance( lowerbound, float ):
        lowerbound = lowerbound * np.ones(dimensions)
    if isinstance( upperbound, float ):
        upperbound = upperbound * np.ones(dimensions)

    print """
#################################
###                           
### CONSTRAINT """ + str(check_constrain.counter) + ')'
    print '###', lowerboundNAME, '<', variablename, '<', upperboundNAME, unit
    print """###
#################################
"""

        
    for i in xrange(dimensions):        
        if ( variable[i] < lowerbound[i] ):
            print 'VIOLATION', i+1, 'constraint', check_constrain.counter
            print str(variablename) +  ' too small value:' + str(round( variable[i], 1)).rjust(wi) + ' lower bound:' + str(round( lowerbound[i], 1)).rjust(wi) + ' unit ', unit
            check_constrain.checkoutALA[i, check_constrain.counter] += 1
        if ( variable[i] > upperbound[i]   ):
            print 'VIOLATION', i+1, 'constraint', check_constrain.counter
            print str(variablename) +  '  too big value:' + str(round( variable[i], 1)).rjust(wi) + ' upper bound:' + str(round( upperbound[i], 1)).rjust(wi) + ' unit ', unit
            check_constrain.checkoutYLA[i, check_constrain.counter] += 1
        if ( check_constrain.checkoutALA[ i, check_constrain.counter ] > 0 or check_constrain.checkoutYLA[ i, check_constrain.counter ] > 0 ):
            print ' '        
    if ( sum(check_constrain.checkoutALA[ :, check_constrain.counter ]) == 0. and sum(check_constrain.checkoutYLA[ :, check_constrain.counter ]) == 0. ):
        print 'Constraint', check_constrain.counter, 'is OK'
print ' '
check_constrain.counter = 0
alvl=2500.
cp=1005.
check_constrain( q_inv,     1.,                     q_pbl,                      'q_inv',      '1',                     'q_pbl',     'g/kg', nroCases )
check_constrain( tpot_inv,  1. + (alvl/cp)*clw_max, 15.,                        't_inv',      '1 + (alvl/cp)*clw_max', '15',        'K'   , nroCases )
check_constrain( pblh,      80.,                    3000.,                      'pblh',       '80',                    '3000',      'm'   , nroCases )
check_constrain( cloudbase, 30.,                    pblh-50.*np.ones(nroCases), 'cloud base', '30',                    'pblh - 50', 'm'   , nroCases )
check_constrain( clw_max,   0.0,                    10000.,                     'clw_max',    '0.0',                   'INF',       'g/kg', nroCases )

# forming CSV
for i in xrange(nroCases):
    temp = caselist[i].rjust(next(tabs)) 
    for k in xrange(np.size(design,1)):
        temp = temp + str( round( float( design[i, k] ),2 ) ).rjust( next(tabs) )
    temp = temp + str( round(cloudbase[i],2)).rjust(next(tabs)) + str( round(pblh[i] - cloudbase[i], 2)).rjust(next(tabs))
    if LWC:
        temp = temp + str(round(q_pbl[i],2)).rjust(next(tabs))
    else:
        temp = temp + str(round(clw_max[i],2)).rjust(next(tabs))
    printtaus.append( temp )

## ARGUMENT MINIMUM
#temp = "argmin".rjust( next(tabs))
#for k in xrange(np.size(design,1)):
#    temp = temp + str( np.argmin( design[:,k] ) +1 ).rjust( next(tabs))
#temp = temp + str( np.argmin(cloudbase) +1 ).rjust(next(tabs)) + str( np.argmin( pblh - cloudbase ) +1 ).rjust(next(tabs))
#if LWC:
#    temp = temp + str( np.argmin( q_pbl ) +1 ).rjust(next(tabs))
#else:
#    temp = temp + str( np.argmin(clw_max)  +1).rjust(next(tabs))
#print temp    
#
## ARGUMENT MAXIMUM
#temp = "argmax".rjust( next(tabs))
#for k in xrange(np.size(design,1)):
#    temp = temp + str( np.argmax( design[:,k] ) +1 ).rjust(next(tabs))
#temp = temp + str( np.argmax( cloudbase ) +1 ).rjust(next(tabs)) + str( np.argmax( pblh - cloudbase ) +1 ).rjust(next(tabs))
#if LWC:
#    temp = temp + str( np.argmax(  q_pbl) +1 ).rjust(next(tabs))
#else:
#    temp = temp + str( np.argmax(clw_max) +1 ).rjust(next(tabs))
#print temp 
   
argminimums = np.zeros(Ltabs-1)
argmaximums = np.zeros(Ltabs-1)

for k in xrange(np.size(design,1)):
    argminimums[k] = np.argmin( design[:,k])
argminimums[k+1] = np.argmin(cloudbase)
argminimums[k+2] = np.argmin(pblh-cloudbase)
if LWC:
    argminimums[k+3] = np.argmin(q_pbl)
else:
    argminimums[k+3] = np.argmin(clw_max)

for k in xrange(np.size(design,1)):
    argmaximums[k] = np.argmax( design[:,k])
argmaximums[k+1] = np.argmax(cloudbase)
argmaximums[k+2] = np.argmax(pblh-cloudbase)
if LWC:
    argmaximums[k+3] = np.argmax(q_pbl)
else:
    argmaximums[k+3] = np.argmax(clw_max)
    

#printing
 
kokonaispituus = sum(listtabs)
sys.stdout.write(printtaus[0])
sys.stdout.write('\n')
for printindeksi in xrange( 1, np.size(printtaus) ): #
    i = printindeksi - 1
    argu = np.column_stack(( argminimums, argmaximums))
    if i in argu: #i in argminimums:            
        summat = [0]
        MINI = np.where( argminimums == i)[0]
        MAXI = np.where( argmaximums == i)[0]
        molemmat = np.sort( np.concatenate((MINI,MAXI )) )
        for k in np.where( argu == i)[0]:
            
            summat.append( sum(listtabs[:min(k+1,len(listtabs))]) )
            summat.append( sum(listtabs[:min(k+2,len(listtabs))]) )
        summat.append( sum(listtabs) )
        
        itersu=iter(summat)
        alku = itersu.next()
        loppu = itersu.next()
        sys.stdout.write( printtaus[printindeksi][alku : loppu ] )
        for k in xrange(len(molemmat)): #(len(summat)-2)/2
            alku = loppu
            loppu = itersu.next()
            if molemmat[int(k)] in MINI and molemmat[int(k)] not in MAXI:
                sys.stdout.write( colored( printtaus[printindeksi][ alku : loppu ], "blue" ) ) if COLOR else sys.stdout.write( printtaus[printindeksi][ alku : loppu ] )
            elif molemmat[int(k)] in MAXI and molemmat[int(k)] not in MINI:
                sys.stdout.write( colored( printtaus[printindeksi][ alku : loppu ], "red" ) ) if COLOR else sys.stdout.write( printtaus[printindeksi][ alku : loppu ] )
            elif molemmat[int(k)] in MAXI and molemmat[int(k)] in MINI:
                sys.stdout.write( colored( printtaus[printindeksi][ alku : loppu ], "green" ) ) if COLOR else sys.stdout.write( printtaus[printindeksi][ alku : loppu ] )
            else:
                sys.stdout.write( printtaus[printindeksi][alku : loppu ] )                
                
            alku = loppu
            loppu = itersu.next()
            sys.stdout.write(printtaus[printindeksi][alku : loppu ])
        sys.stdout.write('\n')

    else:
        sys.stdout.write( printtaus[printindeksi] )
        sys.stdout.write('\n')

    
print '-'*Stabs
print printtaus[0]
print '-'*Stabs
# MINIMUM
temp = "min".rjust( next(tabs))
for k in xrange(np.size(design,1)):
    temp = temp + str( round(np.min( design[:,k] ),2) ).rjust( next(tabs))
temp = temp + str( round( np.min(cloudbase),2 )).rjust(next(tabs))+ str( round( np.min( pblh - cloudbase ),2 )).rjust(next(tabs))
if LWC:
    temp = temp + str(round( np.min(  q_pbl),2 )).rjust(next(tabs))
else:
    temp = temp + str(round( np.min(clw_max),2 )).rjust(next(tabs))
print temp    

# MAXIMUM
temp = "max".rjust( next(tabs))
for k in xrange(np.size(design,1)):
    temp = temp + str( round(np.max( design[:,k] ),2) ).rjust( next(tabs))
temp = temp + str( round( np.max(cloudbase),2 )).rjust(next(tabs))+ str( round( np.max( pblh - cloudbase ),2 )).rjust(next(tabs))
if LWC:
    temp = temp + str(round( np.max(  q_pbl),2 )).rjust(next(tabs))
else:
    temp = temp + str(round( np.max(clw_max),2 )).rjust(next(tabs))
print temp

# MEAN
temp = "mean".rjust( next(tabs))
for k in xrange(np.size(design,1)):
    temp = temp + str( round(np.average( design[:,k] ),2) ).rjust( next(tabs))
temp = temp + str( round( np.average(cloudbase),2 )).rjust(next(tabs))+ str( round( np.average( pblh - cloudbase ),2 )).rjust(next(tabs))
if LWC:
    temp = temp + str(round( np.average(  q_pbl),2 )).rjust(next(tabs))
else:
    temp = temp + str(round( np.average(clw_max),2 )).rjust(next(tabs))
print temp     

# ARGUMENT MINIMUM
temp = "argmin".rjust( next(tabs))
for k in xrange(np.size(design,1)):
    temp = temp + str( np.argmin( design[:,k] ) +1 ).rjust( next(tabs))
temp = temp + str( np.argmin(cloudbase) +1 ).rjust(next(tabs)) + str( np.argmin( pblh - cloudbase ) +1 ).rjust(next(tabs))
if LWC:
    temp = temp + str( np.argmin( q_pbl ) +1 ).rjust(next(tabs))
else:
    temp = temp + str( np.argmin(clw_max)  +1).rjust(next(tabs))
print temp    

# ARGUMENT MAXIMUM
temp = "argmax".rjust( next(tabs))
for k in xrange(np.size(design,1)):
    temp = temp + str( np.argmax( design[:,k] ) +1 ).rjust(next(tabs))
temp = temp + str( np.argmax( cloudbase ) +1 ).rjust(next(tabs)) + str( np.argmax( pblh - cloudbase ) +1 ).rjust(next(tabs))
if LWC:
    temp = temp + str( np.argmax(  q_pbl) +1 ).rjust(next(tabs))
else:
    temp = temp + str( np.argmax(clw_max) +1 ).rjust(next(tabs))
print temp 

print '-'*Stabs
print printtaus[0]


#print ' '
#print 'abs'
#for i in xrange(nroCases):
#    if ( pblh[i] - cloudbase[i]  < 50.  ):
#        print str(i+1).rjust(2), str(round(cloudbase[i], 2)).rjust(9), str(round(pblh[i], 2)).rjust(9), str(round(pblh[i]- cloudbase[i], 2)).rjust(9)

        
print ' '
print 'lower boundary violations'
for i in xrange(nroCases):
    a = str(int(check_constrain.checkoutALA[i,0])).rjust(2) + ' |'
    if ( sum( check_constrain.checkoutALA[i, 1:] ) > 0.0 ):
        for k in check_constrain.checkoutALA[i,1:]: 
            if int(k) == 0:
                p = '.'
            else:
                p = 'X'
            a = a + ' ' + p
        print a
            
print ' '
print 'upper boundary violations'
for i in xrange(nroCases):
    a = str(int(check_constrain.checkoutYLA[i,0])).rjust(2) + ' |'
    if ( sum( check_constrain.checkoutYLA[i, 1:] ) > 0.0 ):
        for k in check_constrain.checkoutYLA[i,1:]: 
            if int(k) == 0:
                p = '.'
            else:
                p = 'X'
            a = a + ' ' + p
        print a
        
ala = sum( check_constrain.checkoutALA[ :, 1: ] )
yla = sum( check_constrain.checkoutYLA[ :, 1: ] )
print ' '

lbv = 'lower boundaries violations'
vv  = len(lbv)
for i in ala:
    lbv = lbv + ' ' + str(int(i)).rjust(2)
print lbv
    
ubv = 'upper boundaries violations'
for i in yla:
    ubv = ubv + ' ' + str(int(i)).rjust(2)
print ubv

ww = 3*check_constrain.counter
vv = vv + ww
uu = '-'*vv
print uu
print 'total number of violations ' + str(int(sum( ala + yla ))).rjust(ww)
print ' '
print 'use LWC from csv', LWC


print 'version', tag
