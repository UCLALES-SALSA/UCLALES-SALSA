import sys
import numpy as np
import copy
def calc_psat_w(T):
    # Function calculates the saturation vapor pressure (Pa) of liquid water as a function of temperature (K)
    #
    # thrm.f90:  real function rslf(p,t)
    c0=0.6105851e+03
    c1=0.4440316e+02
    c2=0.1430341e+01
    c3=0.2641412e-01
    c4=0.2995057e-03
    c5=0.2031998e-05
    c6=0.6936113e-08
    c7=0.2564861e-11
    c8=-.3704404e-13
    #
    x=max(-80.,T-273.16)
    return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

def rslf(p,t):
#! This function calculates the liquid saturation vapor mixing ratio as
#! a function of temperature and pressure
    # thrm.f90:  real function rslf(p,t)

    c0=0.6105851e+03
    c1=0.4440316e+02
    c2=0.1430341e+01
    c3=0.2641412e-01
    c4=0.2995057e-03
    c5=0.2031998e-05
    c6=0.6936113e-08
    c7=0.2564861e-11
    c8=-.3704404e-13
    #
    x = min(max(-80.,t-273.16),50.)
    esl = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    rslf = .622*esl/(p-esl)
    
    return rslf

def calc_rh(rw,T,press):
    # Calculate RH (%) from water vapor mixing ratio rw (r=m_w/m_air [kg/kg]), temperature (K) and pressure (Pa)
    #
    # r=m_w//m_air=pw/Rm/(pair/R)=pw/(p-pw)*R/Rm => pw=p*r/(R/Rm+r)
    #
    R=287.04    # Specific gas constant for dry air (R_specific=R/M), J/kg/K
    Rm=461.5    # Specific gas constant for water
    ep=R/Rm
    #
    psat=calc_psat_w(T)
    return press*rw/(ep+rw)/psat*100.
    # When ep>>rw => RH=press*rw/(ep*psat)*100


def calc_cloud_base(p_surf,theta,rw):
    # Calulate cloud base heigh when potential temperature  (theta [k]) and water vapor
    # mixing ratio (rw [kg/kg]) are constants. Surface pressure p_surf is given in Pa.
    # For more information, see "lifted condensation level" (LCL).
    #
    # Constants
    R=287.04    # Specific gas constant for dry air (R_specific=R/M), J/kg/K
    Rm=461.5    # -||- for water
    ep2=Rm/R-1.0 #M_air/M_water-1
    cp=1005.0    # Specific heat for a constant pressure
    rcp=R/cp
    cpr=cp/R
    g=9.8
    p00=1.0e+05
    #
    # Iterate pressure for each altitude level (P & T => RH)
    dz=1.            #     1 m resolution
    z=-dz            # The first altitude will be z=0
    press=p_surf    # Start from surface
    RH=0
    while RH<100 and z<10000:
        # From z to z+dz
        z+=dz
        # Virtual temperature: T_virtual=T*(1+ep2*rl)
        xsi=(1+ep2*rw)
        # Temperature (K) and pressure (Pa)
        tavg=theta*(press/p00)**rcp
        press-=g*dz*press/(R*tavg*xsi)
        #
        # Current RH (%)
        RH=calc_rh(rw,tavg,press)
    #
    if z>=10000: return -1
    return z

def calc_rh_profile( theta, wc, z, modify = False, rh_max = 101. ):
    ts = copy.deepcopy(theta)
    rts = copy.deepcopy(wc)
    ps = copy.deepcopy(z)    
    if (rh_max <= 100. and modify ):
        sys.exit('rh_max must be over 100.')
        # Calulate cloud base heigh when potential temperature  (ts [k]) and water vapor
    # mixing ratio (rts [kg/kg]) are constants. Surface pressure p_surf is given in Pa.
    # For more information, see "lifted condensation level" (LCL).
    #
    # Constants
    
    R = 287.04    # Specific gas constant for dry air (R_specific=R/M), J/kg/K
    Rm = 461.5    # -||- for water
    ep2 = Rm/R-1.0 #M_air/M_water-1
    cp = 1005.0    # Specific heat for a constant pressure
    rcp = R/cp
    cpr = cp/R
    g = 9.8
    p00 =1.0e+05
    p00i = 1./p00
    #

    
    # Iterate pressure for each altitude level (P & T => RH)
#    dz_orig = dz
    k=len(ps)
    
    xs = np.zeros(k)
    tks = np.zeros(k)
    thds = np.zeros(k)
    hs = np.zeros(k)
    zold1 = 0.
    zold2 = 0.
    

    for ns in xrange(k):
        
        xs[ns] = (1.+ep2*rts[ns]) 
        if (ns == 0):
            ps[ns] = ps[ns]*100.
            zold2  = 0.
            hs[ns] = 0.
#             dz = dz_orig/2.
            
        else:
            hs[ns] = ps[ns]
            zold1  = zold2
            zold2  = ps[ns]
            tavg   =(ts[ns]*xs[ns]+ts[ns-1]*xs[ns-1]*(p00**rcp)/ps[ns-1]**rcp)*.5       
                 #tavg = ts[i]*( ps/p00 )**rcp
            
            ps[ns] = (ps[ns-1]**rcp-g*(zold2-zold1)*(p00**rcp)/(cp*tavg))**cpr
          
         #
         # Current RH (%)
         #rh_temp = calc_rh( rts[ns] , ts[ns], ps[ns])
#         while (rh_temp > rh_max):
#             
#             h2o = 100./rh_max*rts[i]
#             rts[i] = h2o*1000.
#             rh_temp = calc_rh( h2o , tavg, ps)
#         RH.append( rh_temp)        
#         print z, ps, tavg, RH[-1]
        tks[ns] = ts[ns]*(ps[ns]*p00i)**rcp # absolute temperature 
    for ns in xrange(k):
        if (ns == 0 and xs[ns] > 100.):
            print "WARNING: RH>100% AT GROUND LEVEL"
        
        xs[ns] = 100.*rts[ns] / rslf( ps[ns],tks[ns] )# relative humidity
        thds[ns] = tks[ns]*(p00/ps[ns])**rcp
        row= ' {0:10.1f}{1:10.1f}{2:10.2f}{3:10.2f}{4:10.5f}{5:10.1f}'.format(            \
                                        ps[ns],                       \
                                        hs[ns],                           \
                                        tks[ns],                        \
                                        thds[ns],                  \
                                        rts[ns],                          \
                                        xs[ns]                            )
         
#print row

        print row.expandtabs(10)
#         press -= g*dz*press/(R*tavg*xsi)
    
    if modify:
        return xs, rts
    else:
        return xs, ps

def calc_cloud_droplet_diam( theta, wc, press, num_pbl ): 
    ts = copy.deepcopy(theta)
    rts = copy.deepcopy(wc)
    ps = copy.deepcopy(press)    
    ccn = copy.deepcopy(num_pbl)
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
    #

    
    # Iterate pressure for each altitude level (P & T => RH)
#    dz_orig = dz
    k=len(ps)
    
#    xs = np.zeros(k)
    tks = np.zeros(k)
    diam = np.zeros(k)
    wcd = np.zeros(k)
#    thds = np.zeros(k)
#    hs = np.zeros(k)
    
    for ns in xrange(k):
        til = ts[ns]*( ps[ns]*p00i )**rcp
        xx = til
        yy = rslf( ps[ns],xx )
        zz = max( rts[ns]-yy, 0. )
        erotus = 10.
        
#        RH = calc_rh( rw, til, ps[ns] )
        i = 0
        if (zz > 0.):
            #for i in xrange(3):
            while ( erotus > 0.01 and i < 25 ): #abs( til -xx) > 0.1 or i == 0
                x1 = alvl / ( cp*xx )
                xx = xx - ( xx - til*( 1.+x1*zz ) )/( 1. + x1*til*( zz/xx + (1.+yy*ep)*yy*alvl/( Rm*xx*xx )))
                yy = rslf( ps[ns],xx )
                zz = max( rts[ns]-yy, 0. )
#                print 'erotus '+ str(abs( til -xx))
                erotus = abs( til -xx)
                i += 1
#            print 'iteraatiokierrosten lkm '+str(i)
#            print 'erotus '+ str(abs( til -xx))
        diam[ns] = np.power( (6.*zz/ccn)/(1000.*np.pi), 1./3. )*1e6      
        wcd[ns] = zz
        tks[ns]=xx
    
    return diam, wcd
#######################################################################################    
# Test
# =====
# Take relavant values from "corr_design_15d.csv"
#q_inv=9.04913
#t_inv=9.06275
q_pbl=7.10374    # Boundary layer H2O mixing ratio (g/kg)
t_pbl=293.081        # Boundar layer potential temperature (K) 
#pblh=1.81307
#CCN=6.04246e+06
#
# From Eimear's inputs
p_surf=101780. # Surface pressure (Pa)
#
#print "Cloud base (LCL) is at ",calc_cloud_base(p_surf,t_pbl,q_pbl*1e-3),"m"

