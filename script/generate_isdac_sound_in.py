#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:50:49 2018

@author: aholaj
"""

import os
import numpy as np
from math import modf

def liq_wat_pot_temp( z ):
    
    if z<400 :
        thetaL = 265. + 0.004*(z-400.)
        
    elif z >= 400. and z < 825. :
        thetaL = 265.
    
    elif z >= 825. and z < 2045. :
        thetaL = 266+np.power( z-825., 0.3)
        
    elif z>= 2045. :
        thetaL = 271 + np.power(z-2000., 0.3)
    
    return thetaL

def tot_wat_mix_rat( z ):
    
    if z<400 :
        qt = 1.5-0.00075*(z-400.)
        
    elif z >= 400. and z < 825. :
        qt = 1.5
    
    elif z >= 825. and z < 2045. :
        qt = 1.2
        
    elif z>= 2045. :
        qt = 0.5-0.000075*(z-2045.)
    
    return qt        

def wind_u( z ):
    
    return -7.

def wind_v( z ):
    
    return -2 + 0.003*z

def main( dz = 3.5 ):
#    input_vector = case[0], q_inv[1], tpot_inv[2], q_pbl[3], tpot_pbl[4], pblh[5], num_pbl[6], dz[7], nzp[8], windprofile[9], pres0[10]
    les = os.environ["LES"]
    
    
    
    rootfolder = les + '/bin/case_isdac/'

    if modf(dz)[0] == 0.0 :
        postfix = int(dz)
    else:
        postfix = dz
    
    
    filename =   'sound_in' + str(postfix) # _2045_
#    print filename
    filu = rootfolder + filename
    print filu    
    
    f = open(filu,'w')

    pres0 = 1020.
    
    z = np.arange(0.,1508.+dz, dz)
    
    liqPotTemp = []
    wc = []
    u = []
    v = []

    for zz in z:

        liqPotTemp.append( liq_wat_pot_temp( zz ) ) 
        wc.append( tot_wat_mix_rat( zz ) )          
        u.append( wind_u( zz ) )
        v.append( wind_v( zz ) )

    z[0] = pres0
            


    for k in xrange(len(z)):
        
        row= ' {0:6.1f} {1:6.6f} {2:4.6f} {3:4.6f} {4:4.6f}\n'.format(            \
                                        z[k],                           \
                                        liqPotTemp[k],                  \
                                        wc[k],                          \
                                        u[k],                           \
                                        v[k]                            )
        f.write(row)                                        
    
    f.close()

if __name__ == "__main__":
    main()