#!/usr/bin/python
# -*- coding: utf-8 -*-
# #                    Mounted nzp          file    windprofile  pres0   par/serial   runNroBegin  runNroEnd
# ./emulator_inputs.py True    200         $DESIGN  ideal        1017.8  parallel     1           90
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



from subprocess import call
#import shlex


from FindCloudBase import calc_rh_profile
from FindCloudBase import calc_cloud_droplet_diam
from FindCloudBase import rslf
from ECLAIR_calcs import solve_rw

import sys
import os
import subprocess
import glob

global mounted # if run on local computer
global rootfolder
global designfilu
global tag

def bool_convert(s):
    if s=="True":
        r = True
    elif s=="False":
        r = False
    else:
        r = True
    return r

mounted = True

if __name__ == "__main__":

    if ( len(sys.argv) > 2):
        print "Mounted:", sys.argv[1], "nzp:",sys.argv[2], "filu:", sys.argv[3], "windprofile:",  sys.argv[4], "pres0:", sys.argv[5], "runmode:", sys.argv[6], "runNroBegin:", sys.argv[7], "runNroEnd:",  sys.argv[8], "level:",  sys.argv[9]
        mounted = bool_convert( sys.argv[1] )

      
if mounted:
    import matplotlib.pyplot as plt
    from ModDataPros import plottaa
    from plot_profiles import PlotProfiles
    from ModDataPros import initializeColors

home= os.environ["HOME"]
ibrix = os.environ["IBRIXMOUNT"]
les = os.environ["LES"]

lesroot    = les + '/'
designroot = ibrix + '/DESIGN/'

cwd = os.getcwd()

os.chdir(designroot)

tag =subprocess.check_output("git describe --tags | tr -dc '[:alnum:].'", shell=True)
rootfolder = lesroot + 'bin/case_emulator_DESIGN_' + tag + '/'
for file in glob.glob("*.csv"):
    designbasename=file

os.chdir(cwd)

designfilu = designroot + designbasename

#if (len(sys.argv) > 1):
#    for i in sys.argv[1:]:
#        rootfolder = sys.argv[1]
#        filu = sys.argv[2]


def deltazvanha( pblh, nzp ):
    deltaz = max( 1., min( 20. , round( 1.3333 * pblh / float(nzp) ) ) )
   
   
    return deltaz, nzp

def deltaz( pblh, nzp_orig ):
    dz = max( 1., min( 20. , round( max( 1.3333 * pblh, pblh + 500.) / float(nzp_orig) ) ) )

    if ( dz < 10.0 ):
        deltaz = 10.0
        nzp = int(( nzp_orig-1.5 )*dz/deltaz )
    else:
        deltaz = dz
        nzp = nzp_orig

    return deltaz, nzp
    
def define_z( pblh, deltaz, nzp ):
#z = [ 100., 200., 500.,790.,795.,800.,850.,900.,950.,1000.,1050.,1100.,1150.,1200.,1250.,1300.,1350.,1400.,1450.,1500.,1550.,1600.,1650.,1700.,1750.,1800.,1850.,1900.,1950.,2000.]
    
    zt = []    
    
    zt.append(-deltaz/2.)
    
    for i in xrange(nzp-1):
        zt.append(zt[i]+deltaz)
    
    if (pblh > zt[-1]):
        sys.exit("MODEL HEIGHT NOT ABOVE PBLH")
    
    return zt

def tot_wat_mix_rat_LIN(  z, pblh, q_pbl, q_inv, invThi = 50. ):
    if ( z < pblh ):
        q = q_pbl # g/kg
    elif (z > pblh + invThi ):
        q = q_pbl - q_inv - ( z - ( pblh + invThi ))*tot_wat_mix_rat_LIN( pblh + invThi, pblh, q_pbl, q_inv, invThi )/2000.
    else:
        q = q_pbl - q_inv/invThi * ( z - pblh )
       
    return q

def pot_temp_LIN( z, pblh, tpot_pbl, tpot_inv, t_grad, invThi = 50., dtdzFT = 0.003  ):
    if (z < pblh ):
        theta = tpot_pbl
    elif( z > pblh + invThi):
        theta = tpot_pbl + tpot_inv + dtdzFT*( z - ( pblh + invThi ))
    else:
        theta = tpot_pbl + t_grad * ( z - pblh )
    
    return theta


def tot_wat_mix_rat_IT( z, pblh, q_pbl, q_inv, invThi = 50., q_toa = 2.  ):
    if ( z < pblh ):
        q = q_pbl # g/kg
    elif (z > pblh + invThi ):
        q = (q_pbl - q_inv) - ( q_pbl - q_inv - q_toa) * (1. - np.exp( -( z-pblh-invThi ) / 500. ) ) # g/kg
    else:
        q = -q_inv/invThi * ( z - pblh ) + q_pbl
        
       
    return q

def pot_temp_IT( z, pblh, tpot_pbl, tpot_inv, invThi = 50. ):
    if   (z < pblh ):
        theta = tpot_pbl
    elif (z > pblh + invThi ):
        theta = tpot_pbl + tpot_inv +np.power( z-pblh-invThi , 1./3. )
    else:
        theta = tpot_inv/invThi * ( z - pblh )  + tpot_pbl
    
    return theta

def dthcon(tpot_pbl, sst, pres0 ):
    # CONSTANTS    
    R = 287.04    # Specific gas constant for dry air (R_specific=R/M), J/kg/K
    Rm = 461.5    # -||- for water
    ep = R/Rm
    ep2 = Rm/R-1.0 #M_air/M_water-1
    cp = 1005.0    # Specific heat for a constant pressure
    rcp = R/cp
    cpr = cp/R
    g = 9.8
    p00 =1.0e+05
    p00i = 1./p00
    alvl   = 2.5e+06
    ############
    dthcon = tpot_pbl - sst*(p00/pres0)**rcp
    
    return dthcon
    
def drtcon( q_pbl, sst, pres0 ):
    drtcon = q_pbl - rslf(pres0,sst)

    return drtcon    

def absT( theta, p, conversion = 100. ):
    p = p*conversion # conversion from hPa to Pa by default
    R = 287.04    # Specific gas constant for dry air (R_specific=R/M), J/kg/K
    p00 =1.0e+05
    p00i = 1./p00
    cp = 1005.0    # Specific heat for a constant pressure
    rcp = R/cp

    absolut = theta*(p*p00i)**rcp
    
    return absolut

def potT( t, p, conversion = 100.):
    p = p*conversion # conversion from hPa to Pa by default
    R = 287.04    # Specific gas constant for dry air (R_specific=R/M), J/kg/K
    p00 =1.0e+05
#    p00i = 1./p00
    cp = 1005.0    # Specific heat for a constant pressure
    rcp = R/cp

    theta = t*(p00/p)**rcp
    
    return theta

def wind_0(z):
    u = 0.
    v = 0.
    return u,v

def wind_ideal(z):
    u = 10.
    v = 0.
    return u,v    
    
def wind_dycoms(z):
    u = 3. + 4.2*z/1000.
    v = -9. + 5.6*z/1000.
    return u,v


def wind_ascos(z):
    ascos = PlotProfiles('sound_in', "bin/case_ascos/")
    u = ascos.returnUAppr( z )
    v = ascos.returnVAppr( z )

    return u,v

def calcWind(u,v):
    x = np.asarray(u)
    y = np.asarray(v)    
    t = np.column_stack(( x, y))
    wind = np.zeros(len(u))
    for k in xrange(len(u)):
        wind[k] = np.linalg.norm(t[k])
    
    return wind
    
def calcWindShear(u,v,z):
    size = len(z)-1
    shear = np.zeros(size)
    for k in xrange(size):
        shear[k] = np.sqrt( np.power( u[k+1] - u[k], 2 ) + np.power( v[k+1] - v[k], 2 ) ) / ( z[k+1] - z[k] )
    
    return shear
        

# t_grad [ K / m]
def thickness( tpot_inv, t_grad = 0.3   ):

    invThi = tpot_inv / t_grad
    return invThi


def read_design( filu ):
    
    riviLKM = subprocess.check_output( "cat " + filu + " | wc -l ", shell=True)
    nroCases = int( riviLKM )-1
    etunolla = 3 #int(len(str(nroCases)))
    
    f = open( filu, 'r' )
    
    design     = np.zeros( (   nroCases,6 ) )
    caselist   = np.chararray( nroCases, itemsize = etunolla)
    
    i=0
    for line in f:
        #print line
            
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
            
        i=i+1   
    
        
        
    f.close()
    
    #      case,        q_inv,       tpot_inv,    clw_max,     tpot_pbl,    pblh,        num_pbl
    return caselist,    design[:,0], design[:,1], design[:,2], design[:,3], design[:,4], design[:,5]
    


def write_sound_in( input_vector ):
#    input_vector = case[0], q_inv[1], tpot_inv[2], q_pbl[3], tpot_pbl[4], pblh[5], num_pbl[6], dz[7], nzp[8], windprofile[9], pres0[10]
    case        =        input_vector[0]
    q_inv       = float( input_vector[1] )
    tpot_inv    = float( input_vector[2] )
    q_pbl       = float( input_vector[3] )
    tpot_pbl    = float( input_vector[4] )
    pblh        = float( input_vector[5] )
    num_pbl     = float( input_vector[6] )
    dz          = float( input_vector[7] )
    nzp         = int(   input_vector[8] )
    windprofile =        input_vector[9]
    pres0       = float( input_vector[10] )
    
    folder = rootfolder+'emul' + case +'/'
    call(['mkdir','-p', folder])
    filename = 'sound_in' #+case
#    print filename
    filu = folder + filename
    
    
    f = open(filu,'w')
    
    #z = [ 100., 200., 500.,790.,795.,800.,850.,900.,950.,1000.,1050.,1100.,1150.,1200.,1250.,1300.,1350.,1400.,1450.,1500.,1550.,1600.,1650.,1700.,1750.,1800.,1850.,1900.,1950.,2000.]
    
#    print row

    
    z = define_z( pblh, dz, nzp )
    print ' '
    print 'case '    + str(case)
    print 'deltaz '  + str(dz)
    print 'korkeus ' + str(z[-1])
    print 'pblh '    + str(round(pblh,2))
    print 'nzp '     + str(nzp)

    print 'tpot_pbl ' + str(round(tpot_pbl,2))
    print 'tpot_inv ' + str(round(tpot_inv,2))
    print 'q_pbl '    + str(round(q_pbl,2))
    print 'q_inv '    + str(round(q_inv,2))
    print 'num_pbl '  + str(round(num_pbl,2))
    z[0] = pres0
    potTemp = [tpot_pbl]
    wc = [q_pbl]
    u = [0.]
    v = [0.]

#    invThi = 50.
    t_grad = 0.3   
    invThi = thickness( tpot_inv, t_grad )
    print 'inversion thickness ' + str(round(invThi,2))

    for k in xrange(1,len(z)):
        
        

#        potTemp.append( pot_temp_IT( z[k], pblh, tpot_pbl, tpot_inv, invThi ) )
#        wc.append(      tot_wat_mix_rat_IT( z[k], pblh, q_pbl, q_inv, invThi, q_toa = min( max( q_pbl-q_inv-0.01, 0.01) , 2. ) ) )

        potTemp.append( pot_temp_LIN(        z[k], pblh, tpot_pbl, tpot_inv, t_grad, invThi, dtdzFT = 0.003 ) ) 
        wc.append(      tot_wat_mix_rat_LIN( z[k], pblh,    q_pbl,    q_inv,  invThi  ) )

        if ( windprofile == 'zero' ):
                u_apu, v_apu = wind_0( z[k] )
        elif ( windprofile == 'ascos' ):
                u_apu, v_apu = wind_ascos( z[k] )
        elif ( windprofile == 'dycoms' ):
                u_apu, v_apu = wind_dycoms( z[k] )
        elif ( windprofile == 'ideal' ):
                u_apu, v_apu = wind_ideal( z[k] )                  
        u.append( u_apu )
        v.append( v_apu )
    if ( windprofile == 'zero' ):
        u_apu, v_apu = wind_0( 0. )
    elif ( windprofile == 'ascos' ):
        u_apu, v_apu = wind_ascos(  0. )
    elif ( windprofile == 'dycoms' ):
        u_apu, v_apu = wind_dycoms(  0. )
    elif ( windprofile == 'ideal' ):
        u_apu, v_apu = wind_ideal(  0. )
    u[0] = u_apu
    v[0] = v_apu        
#    for k in xrange(len(z)):
#        print z[k], potTemp[k]
        
    #    print 'lev '+ str(lev)+' theta '+ str(potTemp(lev))+ ' wc ' + str(wc(lev)) + ' u ' + str(u(lev)) + ' v ' + str(v(lev))

    rh, pressL = calc_rh_profile(  potTemp, np.multiply(np.asarray(wc), 0.001),  z )
    drop, cloudwater = calc_cloud_droplet_diam( potTemp, np.multiply(np.asarray(wc), 0.001),  pressL, num_pbl*1e6)
    cloudwater = np.multiply(cloudwater,1000.)
    wind = calcWind(u,v)
    windshear = calcWindShear(u,v,z)
#    if you want to modify water content use the two following commands
#    rh, wc2 = calc_rh_profile(  potTemp, np.multiply(np.asarray(wc), 0.001),  z, True )
#    wc = np.multiply(wc,1000.)[7]

    # write to sound_in    

    for k in xrange(len(z)):
        
        row= ' {0:15f} {1:15.6f} {2:15.6f} {3:15.6f} {4:15.6f}\n'.format(            \
                                        z[k],                           \
                                        potTemp[k],                  \
                                        wc[k],                          \
                                        u[k],                           \
                                        v[k]                            )
        f.write(row)                                        
    
    f.close()
    
#     plotting if mounted
    if ( mounted ):
                    
        markers=True
        z[0] = 0.
        initializeColors(7)

        plottaa( potTemp, z, tit = case+' liq. pot. temp., pblh: ' + str( round(pblh,2) ) + ' ilt.: ' +str( round(invThi,2) ), xl = 'liquid potential temperature K', yl = 'height m', markers=markers, uusikuva = True )
        plt.axhline( y = pblh )
#        plt.axhline( y = pblh + invThi )
#        plt.plot([tpot_pbl,tpot_pbl+tpot_inv], [pblh,pblh+invThi], color='r', marker='o')
#        plt.ylim( [pblh-1.*dz, pblh + invThi+1*dz])
#    #    ax.set_yticks(z)
#    #    ax.set_yticks([pblh, pblh+invThi], minor=True)
        plt.savefig( folder + case + '_0_'+ 'liquid_potential_temperature'  + '.png', bbox_inches='tight')    
        plt.close()
        
        plottaa( wc, z, tit = case+' water mix. rat., pblh: ' + str( round(pblh,2) ) + ' ilt.: ' +str( round(invThi,2) ), xl =  'water mixing ratio g/kg', yl = 'height m', markers=markers, uusikuva = True )
        plt.axhline( y = pblh )
#        plt.axhline( y = pblh + invThi )
#        plt.plot([q_pbl,q_pbl-q_inv], [pblh,pblh+invThi], color='r', marker='o')
#        plt.ylim( [pblh-1.*dz, pblh + invThi+1*dz])
        plt.savefig( folder + case + '_0_'+ 'water_mixing_ratio'  + '.png', bbox_inches='tight')    
        plt.close()

        plottaa( rh, z, tit = case+' relative humidity', xl = 'relative humidity %', yl = 'height m', markers=markers, uusikuva = True )
        plt.savefig( folder + case + '_'+ 'relative_humidity'  + '.png', bbox_inches='tight')    
        plt.close()
        

        plottaa( drop, z, tit = case+' cloud droplet diameter', xl = r'cloud droplet diameter $ \mu m$', yl = 'height m', markers=markers, uusikuva = True )
        plt.savefig( folder + case + '_'+ 'cloud_droplet_diameter'  + '.png', bbox_inches='tight')   
        plt.close()

        plottaa( cloudwater, z, tit = case+' cloud water mixing ratio', xl = 'cloud water mixing ratio g/kg', yl = 'height m', markers=markers, uusikuva = True )
        plt.savefig( folder + case + '_'+ 'cloud_water_mixing_ratio'  + '.png', bbox_inches='tight')         
        plt.close()

        plottaa( wind, z, tit = case+' wind '+ windprofile, xl = 'wind m/s', yl = 'height m', markers=markers, uusikuva = True )
        plt.savefig( folder + case + '_'+ 'wind'  + '.png', bbox_inches='tight')
        plt.close()

        plottaa( windshear, z[:-1], tit = case+' wind shear '+ windprofile, xl = 'wind shear '+ r'$s^{-1}$', yl = 'height m', markers=markers, uusikuva = True )
        plt.savefig( folder + case + '_'+ 'windshear'  + '.png', bbox_inches='tight')
        plt.close()

    else:
        print 'Not plotting initial conditions'
    
    return True
    
def write_namelist( input_vector ):
    case     =        input_vector[0]
    nzp      = int(   input_vector[1] )
    dz       = float( input_vector[2] )
    q_inv    = float( input_vector[3] )
    tpot_inv = float( input_vector[4] )
    clw_max  = float( input_vector[5] )
    tpot_pbl = float( input_vector[6] )    
    pblh     = float( input_vector[7] )
    num_pbl  = float( input_vector[8] )    
    pres0    = float( input_vector[9] )
    level    = input_vector[10]




    folder = rootfolder +'emul' + case 
    call(['mkdir','-p', folder])
    sst = absT( tpot_pbl, pres0 )
    nudge_zmax = pblh - 100.
    print 'generating NAMELIST ' + str(case)
#    print 'design tag', tag, type(tag)
#    dth = dthcon( tpot_pbl, sst, pres0 )
#    drt = drtcon( q_pbl, sst, pres0 )
    #filename = 'NAMELIST'
    #filu = folder+filename
    
#    print num_pbl
#    print len(num_pbl)
#              ' timmax=100.'                                        + \
#              ' Tspinup=0.'                                         + \
#              ' th00=' + str(tpot_pbl)                              +\
#              ' dthcon=' + str(dth)                              +\
#              ' drtcon=' + str(drt)                              +\
#                  ' dtlong=2.'                                          + \


    command = 'dir='       + folder                    +\
              ' design='   + '"' + tag +'"'        +\
              ' case='     + case                      +\
              ' q_inv='    + '"' + str(round(q_inv,2))    + '    [g/kg]' + '"' +\
              ' tpot_inv=' + '"' + str(round(tpot_inv,2)) + '    [g/kg]' + '"' +\
              ' clw_max='  + '"' + str(round(clw_max,2))  + '    [g/kg]' + '"' +\
              ' tpot_pbl=' + '"' + str(round(tpot_pbl,2)) +      '  [K]' + '"' +\
              ' pblh='     + '"' + str(round(pblh,2))     +       ' [m]' + '"' +\
              ' num_pbl='  + '"' + str(round(num_pbl,2))  +  '   [#/mg]' + '"' +\
              ' level='    + level               +\
              ' nzp='      + str(nzp)            +\
              ' deltaz='   + str(dz)             +\
              ' CCN='      + str(num_pbl*1e6)    +\
              ' sst='      + str(sst)            +\
              ' nudge_zmax=' + str(nudge_zmax)   +\
              ' filprf=' + '"' + "'emul" + case + "'"    +'"'                        +\
              ' hfilin=' + '"' + "'emul" + case + ".rst'" +'"'                       +\
              ' n='+'"'+ str(num_pbl)+', 0., 0., 0., 0., 0., 0.' + '"'   +\
              ' ' + cwd + '/generate_namelist.bash'
#
#    2D settings
#
#              ' nxp='      + str(5)            +\
#              ' nyp='      + str(204)            +\
#              ' ssam_intvl=120.' +\
#              ' savg_intvl=120.' +\
#              ' frqanl=1800.' +\
#              ' frqhis=1800.' +\
#              ' nxpart=.false.' +\
#

#               ' Tspinup=10.'     
 #               ' timmax=20.'

#    print command
    os.system(command)
#    args = shlex.split(command)
#    print  args  
#    args = ['timmax=20.',home+'/mounttauskansiot/voimahomemount/UCLALES-SALSA/script/generate_namelist.bash']
#    call(args)
    
    return True
    
def main( nzp_orig=200, filu = designfilu, windprofile = 'ideal', pres0 = 1017.8, runmode='parallel', runNroBegin = 1, runNroEnd = 90, level = '3' ):
    nzp_orig = int(nzp_orig)
    pres0 = float(pres0)
    runNroBegin = int(runNroBegin)
    runNroEnd   = int(runNroEnd)
    A = runNroBegin - 1
    B = runNroEnd    
#    args=['rm','-rf', rootfolder+'*']
#    call(args)
    args ='rm -rf '+ rootfolder+'*'
    os.system(args)
    cwd = os.getcwd()
    designfolder = os.path.dirname( os.path.realpath( filu ) )
    os.chdir( designfolder )
    os.chdir( cwd )

    case, q_inv, tpot_inv, clw_max, tpot_pbl, pblh, num_pbl = read_design( filu )

    if runNroEnd > len(case):
        runNroEnd = len(case)        
        sys.exit('runNroEnd too big: '+ str( runNroEnd ) + ' max value: ' + str( len(case) ) )
    
    
    # define resolutions and nzp's    
    dz = np.zeros(len(case))
    nzp = np.zeros(len(case))
    q_pbl = np.zeros(len(case))    
    for k in xrange( A, B ):
        dz[k], nzp[k] = deltaz( pblh[k], nzp_orig )
        q_pbl[k]      = solve_rw( pres0*100., tpot_pbl[k], clw_max[k]*0.001, pblh[k] )*1000.
        if (q_pbl[k] < 0. ):
            sys.exit('q_pbl NEGATIVE VALUE')
    nzp = nzp.astype(int)
    

#    case     =        input_vector[0]
#    nzp      = int(   input_vector[1] )
#    dz       = float( input_vector[2] )
#    q_inv    = float( input_vector[3] )
#    tpot_inv = float( input_vector[4] )
#    clw_max  = float( input_vector[5] )
#    tpot_pbl = float( input_vector[6] )    
#    pblh     = float( input_vector[7] )
#    num_pbl  = float( input_vector[8] )    
#    pres0    = float( input_vector[9] )
    
    if ( runmode == 'serial' ) :
        print 'serial mode'
        for k in xrange( A, B ): 
            write_sound_in( [ case[k], q_inv[k], tpot_inv[k], q_pbl[k], tpot_pbl[k], pblh[k], num_pbl[k], dz[k], nzp[k], windprofile, pres0 ] )
            write_namelist( [ case[k], nzp[k], dz[k], q_inv[k], tpot_inv[k], clw_max[k], tpot_pbl[k], pblh[k], num_pbl[k], pres0, level] )
            
    elif ( runmode == 'parallel' ):
        print 'parallel mode'
        koko = len(case)
        windprofile = [windprofile]*koko
        pres0       = [pres0]*koko
        level       = [level]*koko
        from multiprocessing import Pool
        pool = Pool(processes= 4)
        sound_in_iter = iter( np.column_stack( ( case[A:B], q_inv[A:B], tpot_inv[A:B], q_pbl[A:B], tpot_pbl[A:B], pblh[A:B], num_pbl[A:B], dz[A:B], nzp[A:B], windprofile[A:B], pres0[A:B] ) ) )
        namelist_iter = iter( np.column_stack( ( case[A:B], nzp[A:B], dz[A:B], q_inv[A:B], tpot_inv[A:B], clw_max[A:B], tpot_pbl[A:B], pblh[A:B], num_pbl[A:B], pres0[A:B], level[A:B] ) ) )
        for i in pool.imap_unordered( write_sound_in, sound_in_iter ):
            print i
        for k in pool.imap_unordered( write_namelist, namelist_iter ):
            print k
        #pool.imap_unordered( write_namelist, namelist_iter )
            

#def dycoms():
#    call(['rm','-rf', rootfolder+'*'])
#    case = 'dycoms'
#    q_inv = 4.45
#    tpot_inv = 6.7
#    q_pbl = 9.45
#    tpot_pbl = 288.3
#    pblh = 795.
#    write_sound_in( case, q_inv, tpot_inv, q_pbl, tpot_pbl, pblh)
#    write_namelist( case, 20., 660. )
#    
if __name__ == "__main__":
    
    if ( len(sys.argv) > 2):
        print sys.argv[2]
        main( nzp_orig = sys.argv[2], filu =  sys.argv[3], windprofile = sys.argv[4], pres0 = sys.argv[5], runmode = sys.argv[6], runNroBegin = sys.argv[7], runNroEnd =  sys.argv[8], level = sys.argv[9] )
        print 'generated by using Command Line arguments'
    else:
        main()
        print 'generated by using default values'


#main(200, filu, 'ascos')
#dycoms()
#plot_lopetus()
