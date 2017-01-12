# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 14:20:00 2016

@author: aholaj
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 13:40:45 2016

@author: aholaj
"""

import numpy as np

import matplotlib.pyplot as plt

from subprocess import call
#import shlex

from ModDataPros import plottaa
from ModDataPros import plot_alustus
from ModDataPros import plot_lopetus
import sys

global rootfolder
rootfolder='/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/bin/case_emulator/'

filu ='/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/designs/corr_design_15d_withconstraints.csv'

if (len(sys.argv) > 1):
    for i in sys.argv[1:]:
        rootfolder = sys.argv[1]
        filu = sys.argv[2]

def deltaz( pblh, nzp=200 ):
    deltaz = max( 1., min( 20. , round( 1.3333 * pblh / float(nzp) ) ) )
    
    
    return deltaz
    
def define_z( pblh, deltaz, nzp=200 ):
#z = [ 100., 200., 500.,790.,795.,800.,850.,900.,950.,1000.,1050.,1100.,1150.,1200.,1250.,1300.,1350.,1400.,1450.,1500.,1550.,1600.,1650.,1700.,1750.,1800.,1850.,1900.,1950.,2000.]
    
    zt = []    
    
    zt.append(-deltaz/2.)
    
    for i in xrange(nzp+4):
        zt.append(zt[i]+deltaz)
    
    print 'deltaz ' + str(deltaz)
    print 'korkeus ' + str(zt[-1])
    print 'pblh ' + str(pblh)
    
    return zt

def tot_wat_mix_rat( z, pblh, q_pbl, q_inv, q_toa = 2. ):
    if ( z < pblh ):
        q = q_pbl # g/kg
    else:
        q = (q_pbl - q_inv) - ( q_pbl - q_inv - q_toa) * (1. - np.exp( -( z-pblh ) / 500. ) ) # g/kg
       
    return q

def liq_pot_temp( z, pblh, tpot_pbl, tpot_inv ):
    if (z < pblh ):
        theta = tpot_pbl
    else:
        theta = tpot_pbl + tpot_inv +np.power( z-pblh , 1./3. )
    
    return theta
        


# westerly winds
def u_zonal(z):
#    u = 3. + 4.2*z/1000.
    u = 0.
    return u

# southerly winds
def v_meridional(z):
#    v = -9. + 5.6*z/1000.
    v= 0.
    return v

def read_design( filu  ):
    
    
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
            pblh.append( 1000.*float(pblhA) )           
            num_pbl.append( float(num_pblA) )
            
        i=i+1   
    
        
        
    f.close()
    
    return case, q_inv, tpot_inv, q_pbl, tpot_pbl, pblh, num_pbl


def write_sound_in( case, q_inv, tpot_inv, q_pbl, tpot_pbl, pblh, dz, nzp=200, pres0 = 1017.8, u0 = 0., v0 = 0. ):

    
    #pres0 = 1017.8
    t0 = tpot_pbl
    q0 = q_pbl
    #u0 = 3.
    #v0 = -9.
        
    folder= rootfolder+'emul' + case +'/'
    call(['mkdir','-p', folder])
    filename='sound_in' #+case
#    print filename
    filu=folder+filename
    
    
    f=open(filu,'w')
    
    #z = [ 100., 200., 500.,790.,795.,800.,850.,900.,950.,1000.,1050.,1100.,1150.,1200.,1250.,1300.,1350.,1400.,1450.,1500.,1550.,1600.,1650.,1700.,1750.,1800.,1850.,1900.,1950.,2000.]
    
#    print row

    
    z = define_z( pblh, dz, nzp )

    z[0] = pres0
    liqPotTemp = [t0]
    wc = [q0]
    u = [u0]
    v = [v0]


    for k in xrange(1,len(z)):
        
        liqPotTemp.append( liq_pot_temp( z[k], pblh, tpot_pbl, tpot_inv ) )
        wc.append(      tot_wat_mix_rat( z[k], pblh,    q_pbl,    q_inv ) )
        u.append(               u_zonal( z[k] ) )
        v.append(          v_meridional( z[k] ) )
    #    print 'lev '+ str(lev)+' lpt '+ str(liqPotTemp(lev))+ ' wc ' + str(wc(lev)) + ' u ' + str(u(lev)) + ' v ' + str(v(lev))

    for k in xrange(len(z)):
        
        row= ' {0:4f} {1:4f} {2:4f} {3:4f} {4:4f}\n'.format(            \
                                        z[k],                           \
                                        liqPotTemp[k],                  \
                                        wc[k],                          \
                                        u[k],                           \
                                        v[k]                            )
        f.write(row)                                        
    
    plot_alustus()
    plottaa( liqPotTemp[1:], z[1:], case+' liquid potential temperature', 'liquid potential temperature', 'height' )
    plt.savefig( folder + case + '_'+ 'liquid potential temperature'  + '.png', bbox_inches='tight')    
    
    plot_alustus()
    plottaa( wc[1:], z[1:], case+' water mixing ratio', 'water mixing ratio', 'height' )
    plt.savefig( folder + case + '_'+ 'water mixing ratio'  + '.png', bbox_inches='tight')    
    

       
    
    f.close()
    
def write_namelist(case, dz, num_pbl, nzp=200):
    import os
    folder = rootfolder +'emul' + case +'/'
    call(['mkdir','-p', folder])
    #filename = 'NAMELIST'
    #filu = folder+filename
    
#    print num_pbl
#    print len(num_pbl)
    command = 'dir='+folder                                         + \
              ' nzp='+str(nzp)                                      + \
              ' deltaz=' + str(dz)                                    + \
              ' dtlong=2.'                                         + \
              ' level=3'                                            + \
              ' filprf=' + '"' + "'emul" + case + "'"    +'"'           + \
              ' hfilin=' + '"' + "'emul" + case + ".rst'" +'"'           + \
              ' CCN=' + str(num_pbl)                                  + \
              ' n='+'"'+'125., ' + str(num_pbl*1e-6)+' , 0., 0., 0., 0., 0.' + '"'   + \
              ' /home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/script/generate_namelist.bash'

#               ' Tspinup=10.'     
 #               ' timmax=20.'

    #print command
    os.system(command)
#    args = shlex.split(command)
#    print  args  
#    args = ['timmax=20.','/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/script/generate_namelist.bash']
#    call(args)
    
def main( nzp=200, filu ='/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/designs/corr_design_15d_withconstraints.csv' ):
    args=['rm','-rf', rootfolder+'*']
    call(args)
    case, q_inv, tpot_inv, q_pbl, tpot_pbl, pblh, num_pbl = read_design( filu )
    

    for k in xrange(len(case)): #
        print ' '
        print case[k]
        dz = deltaz( pblh[k] )    
        write_sound_in( case[k], q_inv[k], tpot_inv[k], q_pbl[k], tpot_pbl[k], pblh[k], dz, nzp )
        write_namelist( case[k], dz, num_pbl[k], nzp )

def dycoms():
    call(['rm','-rf', rootfolder+'*'])
    case = 'dycoms'
    q_inv = 4.45
    tpot_inv = 6.7
    q_pbl = 9.45
    tpot_pbl = 288.3
    pblh = 795.
    write_sound_in( case, q_inv, tpot_inv, q_pbl, tpot_pbl, pblh)
    write_namelist( case, 20., 660. )
    
main(200, filu)
#dycoms()
#plot_lopetus()