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
from FindCloudBase import calc_rh_profile
from FindCloudBase import calc_cloud_droplet_diam
from plot_profiles import PlotProfiles
import sys

global rootfolder
rootfolder='/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/bin/case_emulator/'

filu ='/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/designs/corr_design_15d.csv'
#filu='/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/designs/vanhat/corr_design_15d_withconstraints.csv'

if (len(sys.argv) > 1):
    for i in sys.argv[1:]:
        rootfolder = sys.argv[1]
        filu = sys.argv[2]


def deltaz( pblh, nzp ):
    deltaz = max( 1., min( 20. , round( 1.3333 * pblh / float(nzp) ) ) )
   
   
    return deltaz, nzp

def deltaz2( pblh, nzp ):
    dz = max( 1., min( 20. , round( 1.3333 * pblh / float(nzp) ) ) )

    deltaz = 10.0    
    
    if ( dz < 10.0 ):
        nzp = int(nzp*dz/deltaz)
        deltaz = 10.0
    else:
        deltaz = dz
    
    return deltaz, nzp
    
def define_z( pblh, deltaz, nzp ):
#z = [ 100., 200., 500.,790.,795.,800.,850.,900.,950.,1000.,1050.,1100.,1150.,1200.,1250.,1300.,1350.,1400.,1450.,1500.,1550.,1600.,1650.,1700.,1750.,1800.,1850.,1900.,1950.,2000.]
    
    zt = []    
    
    zt.append(-deltaz/2.)
    
    for i in xrange(nzp+4):
        zt.append(zt[i]+deltaz)
    
    print 'deltaz '  + str(deltaz)
    print 'korkeus ' + str(zt[-1])
    print 'pblh '    + str(pblh)
    print 'nzp '     + str(nzp)
    if (pblh > zt[-1]):
        sys.exit("MODEL HEIGHT NOT ABOVE PBLH")
    
    return zt

def tot_wat_mix_rat( z, pblh, q_pbl, q_inv, invThi = 0., q_toa = 2. ):
    if ( z < pblh ):
        q = q_pbl # g/kg
    else:
        q = (q_pbl - q_inv) - ( q_pbl - q_inv - q_toa) * (1. - np.exp( -( z-pblh ) / 500. ) ) # g/kg
       
    return q

def pot_temp( z, pblh, tpot_pbl, tpot_inv, invThi=0.  ):
    if (z < pblh ):
        theta = tpot_pbl
    else:
        theta = tpot_pbl + tpot_inv +np.power( z-pblh , 1./3. )
    
    return theta
    
    
    
        
def tot_wat_mix_rat_IT( z, pblh, q_pbl, q_inv, invThi = 50., q_toa = 2.  ):
    if ( z < pblh ):
        q = q_pbl # g/kg
    elif (z > pblh + invThi ):
        q = (q_pbl - q_inv) - ( q_pbl - q_inv - q_toa) * (1. - np.exp( -( z-pblh ) / 500. ) ) # g/kg
    else:
        q = -q_inv/invThi * ( z - pblh ) + q_pbl
        
       
    return q

def pot_temp_IT( z, pblh, tpot_pbl, tpot_inv, invThi = 50. ):
    if   (z < pblh ):
        theta = tpot_pbl
    elif (z > pblh + invThi ):
        theta = tpot_pbl + tpot_inv +np.power( z-pblh , 1./3. )
    else:
        theta = tpot_inv/invThi * ( z - pblh )  + tpot_pbl
    
    return theta

    

def wind_0(z):
    u = 0.
    v = 0.
    return u,v
    
def wind_dycoms(z):
    u = 3. + 4.2*z/1000.
    v = -9. + 5.6*z/1000.
    return u,v


def wind_ascos(z):
    ascos=PlotProfiles('sound_in', "bin/case_ascos/")
    u = ascos.returnUAppr( z )
    v = ascos.returnVAppr( z )

    return u,v


# t_grad [ K / m]
def thickness( tpot_inv, t_grad = 0.2   ):

    invThi = tpot_inv / t_grad
    return invThi


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


def write_sound_in( case, q_inv, tpot_inv, q_pbl, tpot_pbl, pblh, dz, nzp, num_pbl, pres0 = 1017.8, u0 = 0., v0 = 0. ):

    
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
    print 'tpot_pbl ' + str(tpot_pbl)
    print 'tpot_inv ' + str(tpot_inv)
    print 'q_pbl ' + str(q_pbl)
    print 'q_inv ' + str(q_inv)
    print 'num_pbl ' + str(num_pbl)
    z[0] = pres0
    potTemp = [t0]
    wc = [q0]
    u = [u0]
    v = [v0]

    invThi = 50.    
#    invThi = thickness( tpot_inv, t_grad = 0.2 )
    print 'inversion thickness ' + str(invThi)

    for k in xrange(1,len(z)):
        
        

        potTemp.append( pot_temp_IT( z[k], pblh, tpot_pbl, tpot_inv, invThi ) )
        wc.append(      tot_wat_mix_rat_IT( z[k], pblh, q_pbl, q_inv, invThi ) )
        u_apu, v_apu = wind_0( z[k] )        
        u.append( u_apu )
        v.append( v_apu )
    
#    for k in xrange(len(z)):
#        print z[k], potTemp[k]
        
    #    print 'lev '+ str(lev)+' theta '+ str(potTemp(lev))+ ' wc ' + str(wc(lev)) + ' u ' + str(u(lev)) + ' v ' + str(v(lev))

    rh, pressL = calc_rh_profile(  potTemp, np.multiply(np.asarray(wc), 0.001),  z )
    drop, cloudwater = calc_cloud_droplet_diam( potTemp, np.multiply(np.asarray(wc), 0.001),  pressL, num_pbl)
    cloudwater = np.multiply(cloudwater,1000.)
#    if you want to modify water content use the two following commands
#    rh, wc2 = calc_rh_profile(  potTemp, np.multiply(np.asarray(wc), 0.001),  z, True )
#    wc = np.multiply(wc,1000.)
    
    # write to sound_in    

    for k in xrange(len(z)):
        
        row= ' {0:15f} {1:15.6f} {2:15.6f} {3:15.6f} {4:15.6f}\n'.format(            \
                                        z[k],                           \
                                        potTemp[k],                  \
                                        wc[k],                          \
                                        u[k],                           \
                                        v[k]                            )
        f.write(row)                                        
    
#     plotting
#    print 'max liq pot temp', np.max(potTemp)
#    print 'm
    z[0] = 0.
    plot_alustus()
    plottaa( potTemp, z, case+' liquid potential temperature', 'liquid potential temperature K', 'height m' )
    plt.savefig( folder + case + '_'+ 'liquid_potential_temperature'  + '.png', bbox_inches='tight')    
    
    plot_alustus()
    plottaa( wc, z, case+' water mixing ratio', 'water mixing ratio g/kg', 'height m' )
    plt.savefig( folder + case + '_'+ 'water_mixing_ratio'  + '.png', bbox_inches='tight')    
#
    plot_alustus()
    plottaa( rh, z, case+' relative humidity', 'relative humidity %', 'height m' )
    plt.savefig( folder + case + '_'+ 'relative_humidity'  + '.png', bbox_inches='tight')    

    plot_alustus()
    plottaa( drop, z, case+' cloud droplet diameter', r'cloud droplet diameter $ \mu m$', 'height m' )
    plt.savefig( folder + case + '_'+ 'cloud_droplet_diameter'  + '.png', bbox_inches='tight')   

    plot_alustus()
    plottaa( cloudwater, z, case+' cloud water mixing ratio', 'cloud water mixing ratio g/kg', 'height m' )
    plt.savefig( folder + case + '_'+ 'cloud_water_mixing_ratio'  + '.png', bbox_inches='tight')         
    
    f.close()
    
def write_namelist(case, dz, num_pbl, nzp):
    import os
    folder = rootfolder +'emul' + case +'/'
    call(['mkdir','-p', folder])
    #filename = 'NAMELIST'
    #filu = folder+filename
    
#    print num_pbl
#    print len(num_pbl)
#              ' timmax=100.'                                        + \
#              ' Tspinup=0.'                                         + \
    command = 'dir='+folder                                         + \
              ' nzp='+str(nzp)                                      + \
              ' deltaz=' + str(dz)                                  + \
              ' dtlong=2.'                                          + \
              ' level=3'                                            + \
              ' filprf=' + '"' + "'emul" + case + "'"    +'"'       + \
              ' hfilin=' + '"' + "'emul" + case + ".rst'" +'"'      + \
              ' CCN=' + str(num_pbl)                                + \
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
    
def main( nzp_orig=200, filu ='/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/designs/corr_design_15d_withconstraints.csv' ):
    args=['rm','-rf', rootfolder+'*']
    call(args)
    case, q_inv, tpot_inv, q_pbl, tpot_pbl, pblh, num_pbl = read_design( filu )
    
#    print 'nzp '     + str(nzp)

    for k in xrange(len(case)): #len(case)
        print ' '
        print case[k]
        
        dz, nzp = deltaz2( pblh[k], nzp_orig )    
#        print 'nzp '     + str(nzp)
#        print 'dz ' + str(dz)
#        print 'write_sound_in'
        write_sound_in( case[k], q_inv[k], tpot_inv[k], q_pbl[k], tpot_pbl[k], pblh[k], dz, nzp, num_pbl[k] )
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