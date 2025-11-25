MODULE mo_salsa_cloud_ice_DET
  USE classSection, ONLY : Section
  USE mo_submctl, ONLY : in2a,fn2b, ira,fra, iia, fia, ncld, nprc, nice, pi, pi6,    & 
       nliquid, nfrozen,                        &
       nlim, prlim, lsicehom, lsiceimm, lsicedep,  &
       dinscheme, &
       boltz, planck, rg, rd, avog,             &
       fixinc, spec,                            &
       lsicenucl
  USE mo_salsa_types, ONLY : aero, cloud, ice, precp, liquid, frozen, rateDiag
  USE mo_particle_external_properties, ONLY : calcSweq
  USE math_functions, ONLY : cumlognorm

  IMPLICIT NONE

  CONTAINS

  !***********************************************
  ! Ice nucleation deterministic
  ! Sources
  !   Mor05   Morrison et al. (Journal of Aerosol Science, 62, 1665-1677, 2005)
  !   U17   Ullrich et al. (2017)
  !   KS98  Khvorostyanov and Sassen (Geophysical Research letters, 25, 3155-3158, 1998)
  !   PK97  Pruppacher and Klett, Microphysics of Clouds and Precipitation, 1997
  !***********************************************
  
  SUBROUTINE ice_nucl_driver(kproma,kbdim,klev,   &
                             ptemp,prv,prs,prsi,ptstep )
        
    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: kproma,kbdim,klev
    REAL, INTENT(in) :: ptstep
    REAL, INTENT(in) :: ptemp(kbdim,klev),  &
                        prv(kbdim,klev),    &
                        prs(kbdim,klev),    &
                        prsi(kbdim,klev)
    
    ! Which species are allowed to freeze
    LOGICAL, PARAMETER :: aerosol_ice = .TRUE., cloud_ice = .TRUE., precip_ice = .TRUE.
    ! Maximum temperature (K) for homogeneous nucleation
    REAL, PARAMETER :: tmax_homog=235.15 ! Kasper: Fixed so that it is 235.15 K Juha: Isnt this a bit high?

    REAL, PARAMETER :: dmin = 1.e-10 ! Diameter limit

    REAL :: frac

    INTEGER :: ii,jj,kk,ss,bb
    REAL :: pf_imm, pf_dep, pf_hom, jf, ns
    REAL :: nnum, omega
    REAL :: Sw_eq, Si, zvol, ra, rb
    REAL :: ddry, dwet, dins

    LOGICAL :: isdry
    
    INTEGER :: ibc, idu, iwa, irim, nspec, ndry, phase, in2b, fn2b
    REAL :: zinsol
    
    ! Ice nucleation modes
    ! 1) Homogeneous freezing: possible for any aqueous droplet with or without insoluble species (DU or BC)
    ! 2) Immersion freezing: possible for aqueous droplets that have an insoluble core (DU or BC)
    ! 3) Deposition freezing: possible for dry insoluble aerosol at sub-saturated conditions (RH < 100%)
    ! 4) Contact freezing: not implemented, because most INPs are immersed in liquid droplets;
    !    Juha: We could (should?) just extend the range of immersion freezing temperatures to replace this, 
    !          since at least in the boundary layer, there's practically never clean ice nuclei, which are 
    !          needed for this and higher up other processes likely dominate.
    
    nspec = spec%getNSpec(type="total")
    ndry = spec%getNSpec(type="dry")
    
    ! Mass (volume) array indices of BC, DU and water
    ibc = spec%getIndex("BC",notFoundValue=0)
    idu = spec%getIndex("DU",notFoundValue=0)
    iwa = spec%getIndex("H2O")
    irim = spec%getIndex("rime")

    in2b = 14
    fn2b = 17

    ! Loop over liquid phase bins
    DO kk = 1,nliquid 
       DO jj = 1,klev
          DO ii = 1,kproma
             IF (ptemp(ii,jj) > 273.15) CYCLE
             IF (liquid(ii,jj,kk)%numc < liquid(ii,jj,kk)%nlim) CYCLE

             IF (dinscheme == 1 .AND. lsicedep) THEN
               ! Part of Phillips 2013 schem to calculate the number of IN per unit surface area of insoluble material
               nnum = SUM(liquid(ii,jj,in2b:fn2b)%numc/(1-liquid(ii,jj,in2b:fn2b)%indef))
               omega = SUM((liquid(ii,jj,in2b:fn2b)%ddry)**2*pi*liquid(ii,jj,in2b:fn2b)%numc/(1-liquid(ii,jj,in2b:fn2b)%indef))
             END IF

             phase = liquid(ii,jj,kk)%phase
             
             ! Get the insoluble volume concentration
             zinsol = 0.
             IF ( ibc > 0 ) zinsol = zinsol + liquid(ii,jj,kk)%volc(ibc)
             IF ( idu > 0 ) zinsol = zinsol + liquid(ii,jj,kk)%volc(idu)

             ! Wet diameter and the diameter of the insoluble part
             CALL liquid(ii,jj,kk)%updateDiameter(limit=.TRUE.,type="all")
             dins = liquid(ii,jj,kk)%dins
             ddry = liquid(ii,jj,kk)%ddry
             dwet = liquid(ii,jj,kk)%dwet
             
             ! Equilibrium saturation ratio
             Sw_eq = calcSweq(liquid(ii,jj,kk),ptemp(ii,jj))

             ! Supersaturation over ice
             Si = prv(ii,jj)/prsi(ii,jj)
             
             ! Immersion freezing (not directly from aerosol)
             pf_imm = 0.
             !IF ( dins > dmin .AND. lsiceimm .AND. ANY(phase == [2,3]) ) THEN
             IF ( dins > dmin .AND. lsiceimm .AND. ANY(phase == [2,3]) .AND. Si>1. .AND. prv(ii,jj)/prs(ii,jj)>0.97) THEN
                
                jf = calc_Jhet(dins,ptemp(ii,jj),Sw_eq)
                pf_imm = 1. - EXP( -pi*dins**2*jf )
             END IF

             ! Deposition freezing
             pf_dep = 0.
             IF ( dins > dmin .AND. prv(ii,jj)/prs(ii,jj)<1.0 .AND. &
                  phase == 1 .AND. lsicedep                           ) THEN
                !write(*,*) 'Depositon nucleation!'
                !Si = prv(ii,jj)/prsi(ii,jj) ! Water vapor saturation ratio over ice
                IF (dinscheme == 0) THEN
                  jf = calc_Jdep(dins,ptemp(ii,jj),Si)
                  pf_dep = 1. - EXP( -pi*dins**2*jf )
                ELSE IF (dinscheme == 1) THEN
                  jf = calc_Jdep_Phi13(dins,ptemp(ii,jj),Si,nnum,omega)
                  pf_dep = 1. - EXP(- jf)
                END IF
                
             END IF
             
             ! Total fraction of particles nucleating ice
             frac = MAX(0., MIN(0.99,pf_imm+pf_dep))
             ! Determine the target ice bin
             bb = getIceBin(dwet)


             CALL iceNucleation(ii,jj,bb,ndry,iwa,irim,liquid(ii,jj,kk),frac)

             
             ! Homogeneous freezing
             pf_hom = 0.
             IF (dwet-dins > dmin .AND. ptemp(ii,jj) < tmax_homog .AND. lsicehom) THEN
                !Si = prv(ii,jj)/prsi(ii,jj) ! Water vapor saturation ratio over ice
                jf = calc_Jhf(ptemp(ii,jj),Sw_eq,Si)
                pf_hom = 1. - EXP( -jf*pi6*(dwet**3 - dins**3)*ptstep )
             END IF

             frac = MAX(0., MIN(0.99,pf_hom))

             CALL iceNucleation_hom(ii,jj,bb,ndry,iwa,irim,liquid(ii,jj,kk),frac)
 
             CALL iceDiagnostics(liquid(ii,jj,kk),pf_imm,pf_dep,pf_hom)
             
             CALL ice(ii,jj,bb)%updateRhomean()

          END DO !ii
       END DO !jj
    END DO

  END SUBROUTINE ice_nucl_driver
  
  !
  ! Ice nucleation based on Khvorostyanov and Curry, Geophys. Res. Lett., 27, 4081-4084, 2000 [KC00] and
  ! Khvorostyanov and Sassen, Geophys. Res. Lett., 25, 3155-3158, 1998 [KS98]
  !
  REAL FUNCTION calc_Jhet(rn,temp,Sw)
    ! The rate of germ formation (#/s) through heterogeneous freezing following
    ! Khvorostyanov and Curry, Geophys. Res. Lett., 27, 4081-4084, 2000 [KC00]
    ! Additional parameters from
    ! Khvorostyanov and Curry, J. Atmos. Sci., 61, 2676-2691, 2004 [KC04]
    ! Khvorostyanov and Sassen, Geophys. Res. Lett., 25, 3155-3158, 1998 [KS98]
    ! Jeffery and Austin, J. Geophys. Res., 102, 25269-25279, 1997
    ! Li et al., J. Geophys. Res., 118, 11213-11227, 2013
    ! - Here heterogeneos freezing includes just immersion freezing mode

    ! Changed all radiuses to diameters!!

    IMPLICIT NONE
    REAL, INTENT(in) :: rn,temp,Sw
    
    REAL :: Tc, act_energy, Lefm, GG, sigma_is, r_g, sf, crit_energy
    REAL, PARAMETER :: & ! Constants
         C = 1.7e10, & ! Constant (1.7e11 dyn cm^-2 = 1.7e11*1e-5/1e-4 N m^-2 = 1.7e10 N m^-2)
         rho_ice = 900., & ! Density of ice (kg/m^3)
         c_1s = 1e32, & ! The concentration of water molecules adsorbed on 1 m^2 a surface (#/m^2)
         T0 = 273.15 ! 0 C in Kelvins
    REAL, PARAMETER :: & ! Case-dependent parameters
         epsi = 0., & ! Elastic strain produced in ice embryo by the insoluble substrate
         alpha = 0.0, & ! Relative area of active sites
         mis = 0.5 ! Cosine of the contact angle
    
    calc_Jhet = 0.
    
    ! Must have a core
    IF (rn<1e-10) RETURN
    
    Tc = temp-T0 ! Temperature in Celsius
    
    calc_Jhet= exp(150.577-0.517*temp)
  END FUNCTION calc_Jhet
  
  
  REAL FUNCTION calc_Jdep(rn,temp,Si)
    ! The rate of germ formation (#/s) through deposition freezing following
    ! Khvorostyanov and Curry, Geophys. Res. Lett., 27, 4081-4084, 2000 [KC00]
    ! Additional parameters from
    ! Hoose et al., J. Atmos. Sci., 67, 2483-2503, 2010 [Ho10]

    ! Changed all radiuses to diameters
    
    IMPLICIT NONE
    REAL, INTENT(in) :: rn,temp,Si
    
    REAL :: Tc
    REAL, PARAMETER :: & ! Constants
         T0 = 273.15 ! 0 C in Kelvins
    REAL, PARAMETER :: & ! Case-dependent parameters
         alpha_d = 285.692, &
         beta_d = 0.017, &
         gamma_d = 256.692, &
         kappa_d = 0.08, &
         lamda_d = 200.745
    
    calc_Jdep = 0.
    
    ! Must have a core and supersaturation over ice
    IF (rn<1e-10 .OR. Si<1.0001) RETURN
    
    Tc = temp-T0 ! Temperature in Celsius
    
    calc_Jdep = exp(alpha_d*(Si-1)**(1./4.)*cos(beta_d*(temp-gamma_d))**2*(pi/2.-atan(kappa_d*(temp-lamda_d)))/pi)
 
  END FUNCTION calc_Jdep

  REAL FUNCTION calc_Jdep_Phi13(rn,temp,Si,nnum,omega)
    ! Deposition freezing
    
    IMPLICIT NONE
    REAL, INTENT(in) :: temp,            & ! Temperature in K
                        rn,           & ! Diameter of the insoluble aerosol
                        Si,           & ! Equilibrium saturation ratio
                        nnum,          & ! B-bin number concentration (insoluble)
                        omega            ! B-bin surface area (insoluble)
    
    REAL :: Tc, act_energy, sigma_iv, d_g, sf, crit_energy

    REAL, PARAMETER :: & ! Constants
         C = 1.7e10, & ! Constant (1.7e11 dyn cm^-2 = 1.7e11*1e-5/1e-4 N m^-2 = 1.7e10 N m^-2)
         rho_ice = 900., & ! Density of ice (kg/m^3)
         T0 = 273.15 ! 0 C in Kelvins
    REAL, PARAMETER :: & ! Case-dependent parameters
         epsi = 0., &    ! Elastic strain produced in ice embryo by the insoluble substrate
         b0 = -1.0261, & ! Parameters for Phillips 2013...
         b1 = 3.1656e-3, &
         b2 = 5.3938e-4, &
         b3 = 8.2584e-6, &
         T0_DU = -40., &
         T1_DU = -30., &
         T2_DU = -10., &
         dT = 5., &
         dSi = 0.1, &
         Sw_0 = 0.97, &
         h = 0.15, &
         gamma = 2., &
         alpha = 0.666667
    REAL :: thrad, costh   ! Contact angle in radians, cosine of contact angle

    REAL :: fc_comp1, fc_comp2, fc, wvp, Sw_r, H_comp1, H_x, Xi, Mu, x_cub, Si_0
    
    calc_Jdep_Phi13 = 0.

    ! Must have a core and supersaturation over ice
    IF (rn<1e-7 .OR. Si<1.0001) RETURN

    Tc = temp-T0 ! Temperature in Celsius

    x_cub = b0 + b1*Tc + b2*Tc**2 + b3*Tc**3

    Si_0 = 1+10**x_cub

    IF (Tc < T0_DU) THEN
      fc_comp1 = 1
    ELSE IF (Tc > T0_DU+dT) THEN
      fc_comp1 = h
    ELSE
      fc_comp1 = Cub_inter(1.,h,Tc,T0_DU,T0_DU+dT)
    END IF

    IF (Si < Si_0) THEN
       fc_comp2 = 0.
    ELSE IF (Si > Si_0+dSi) THEN
       fc_comp2 = 1.
    ELSE
       fc_comp2 = Cub_inter(0.,1.,Si,Si_0,Si_0+dSi) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF

    fc = fc_comp1*fc_comp2/gamma

    !Saturation ratio of vapor with respect to water
    wvp = Si*(EXP(9.550426-5723.265/(Tc+273.15)+3.53068*LOG(Tc+273.15)-0.00728332*(Tc+273.15)))/100
    Sw_r = wvp/(611.2*EXP(17.67*Tc/(Tc+243.5))*0.01)

    !Fractions reducing IN activity at low Si, warm T

    IF (Sw_r < Sw_0) THEN
       H_comp1 = 0.
    ELSE IF (Sw_r > 1) THEN
       H_comp1 = 1.
    ELSE
       H_comp1 = Cub_inter(0.,1.,Sw_r,Sw_0,1.)
    END IF
    
    H_x = MIN(fc + (1-fc)*H_comp1,1.)

    IF (Tc < T1_DU) THEN
       Xi = 1.
    ELSE IF (Tc > T2_DU) THEN
       Xi = 0.
    ELSE
       Xi = Cub_inter(1.,0.,Tc,T1_DU,T2_DU)
    END IF
    
    Mu = H_x*Xi*(nnum/omega)*pi*rn**2

    calc_Jdep_Phi13 = 1-EXP(-Mu)
 
  END FUNCTION calc_Jdep_Phi13
  
  REAL FUNCTION calc_Jhf(temp,Sw,Si)
    ! Homogeneous freezing based on Khvorostyanov and Sassen, Geophys. Res. Lett., 25, 3155-3158, 1998 [KS98]

    IMPLICIT NONE
    REAL, intent(in) :: temp, Sw, Si ! Temperature (K) and water vapor saturation ratio
    
    REAL :: Delta_aw, P, Tc, act_energy, Lefm, GG, sigma_is, r_g, crit_energy
    REAL, PARAMETER :: & ! Constants
         Nc=5.85e16, & ! The number of water molecules contacting unit area of ice germ (#/m^2) [from KC00]
         rho_ice = 900., & ! Density of ice (kg/m^3)
         T0 = 273.15, & ! 0 C in Kelvins
         a = -906.7, &
         b = 8502., &
         c = -26924., &
         d = 29180.
    
    calc_Jhf = 0.
    
    Tc = temp-T0 ! Temperature in Celsius
    ! Critical energy of germ formation
    ! a) Ice germ radius (eq. 9a)
    Lefm = (79.7+0.708*Tc-2.5e-3*Tc**2)*4.1868e3 ! Effective latent heat of melting (eq. 6)
    GG = rg*temp/Lefm/spec%mwa
    IF ( (T0/temp)*Sw**GG<1.0001 ) RETURN
    sigma_is = 28.e-3+0.25e-3*Tc ! Surface tension between ice and solution
    r_g = 2.*sigma_is/( rho_ice*Lefm*log((T0/temp)*Sw**GG) )
    IF (r_g<=1e-10) RETURN

    Delta_aw = (Si-1)*Sw/Si ! Activity of water

    IF ((Delta_aw <= 0.26) .AND. (Delta_aw >= 0.34)) RETURN

    P = a+b*Delta_aw+c*(Delta_aw)**2+d*(Delta_aw)**3
    
    calc_Jhf = 10**P*1e6
    
  END FUNCTION calc_Jhf


  
  REAL function calc_het_sf(x,mis)
    ! Calculate shape factor for heterogeneous ice nucleation.
    !   Khvorostyanov and Curry, J. Atmos. Sci., 61, 2676-2691, 2004
    REAL :: x ! x=r_core/r_crit
    REAL :: mis ! Cosine of the contact angle
    REAL :: fii, fpsi
    !
    IF (x>100.) THEN
       ! Problems with numerical accuracy when x>>1; x=100 seems to be good limit for using the limiting value
       calc_het_sf = (mis**3-3*mis+2)/4
    ELSE
       fii = sqrt(1.-2.*mis*x+x**2)
       fpsi = (x-mis)/fii
       calc_het_sf = 0.5* (1. + ( ( 1.-mis*x)/fii)**3 + (2.-3.*fpsi+fpsi**3)*x**3 + 3.*mis*(fpsi-1.)*x**2)
    ENDIF
  END function calc_het_sf
  
  ! ---------------------------------------------------
  ! Cubic interpolation function by Phillips et al. (2008)
  FUNCTION Cub_inter(a_in,b_in,y,y1,y2) result(inter)
    IMPLICIT NONE
    REAL, intent(in) :: a_in, b_in, y, y1, y2
    REAL :: inter     !The result of interpolation

    REAL :: A, B, a0, a1, a2, a3

    A = 6*(a_in-b_in)/(y2-y1)**3
    B = a_in+A*(y1**3)/6-A*y1**2*y2/2

    a0 = B
    a1 = A*y1*y2
    a2 = -A*(y1+y2)/2
    a3 = A/3
    inter = a0 + a1*y + a2*y**2 + a3*y**3
  END FUNCTION Cub_inter

  !***********************************************
  !
  ! Prescribed ice number concentration for cloudy regions (rc>0.001 g/kg) where ice
  ! supersaturation is at least 5%. Ice number concentration is increased towards the
  ! target concentration (fixinc, #/kg) by converting the largest cloud droplets to ice.
  !
  !***********************************************
  SUBROUTINE ice_fixed_NC(kproma, kbdim,  klev,   &
       ptemp, ppres, prv,  prsi     )
    

 !!!!!!!!!!!!!!!!!!!!!!!!! DOES NOT WORK ANYMORE AS IT IS, NEED TO UPDATE ARRAY INDICES IF THIS IS STILL REQUIRED

    
    IMPLICIT NONE
    INTEGER, INTENT(in) :: kproma,kbdim,klev
    
    REAL, INTENT(in) :: &
         ptemp(kbdim,klev),    &
         ppres(kbdim,klev),    &
         prv(kbdim,klev),    &
         prsi(kbdim,klev)
    
    INTEGER :: ii,jj,kk,ss
    
    REAL :: pdn, iceSupSat, rc_tot, Ni0,  &
         sumICE, iceTendecyNumber,    &
         liqToIceFrac
    
    INTEGER :: iwa, nspec,ndry
    
    iwa = spec%getIndex("H2O")
    nspec = spec%getNSpec(type="total")
    ndry = spec%getNSpec(type="dry")
    
    DO ii = 1,kbdim
       DO jj = 1,klev
          pdn=ppres(ii,jj)/(rd*ptemp(ii,jj)) ! Air density (kg/m^3)
          
          iceSupSat = prv(ii,jj) / prsi(ii,jj)  - 1.0 ! ice supersaturation
          rc_tot = 0.
          DO kk = 1,ncld
             rc_tot = rc_tot + cloud(ii,jj,kk)%volc(iwa) * spec%rhowa/pdn ! cloud water mixing ratio (kg/kg)
          END DO
          
          ! conditions for ice nucleation
          IF ( icesupsat < 0.05 .OR. rc_tot < 0.001e-3  ) CYCLE
          
          ! target number concentration of ice, converted to #/m^3
          Ni0     = fixinc * pdn
          
          ! current ice number concentration (#/m^3)
          sumICE    = SUM(   ice(ii,jj,:)%numc )
          
          IF ( sumICE > Ni0 ) CYCLE
          
          DO kk = nice,1,-1 ! Assuming nice=ncld
             IF ( sumICE < Ni0 .AND. cloud(ii,jj,kk)%numc > nlim) THEN
                
                iceTendecyNumber = max( 0.0, min( Ni0 - ice(ii,jj,kk)%numc , cloud(ii,jj,kk)%numc )  )
                
                ice(ii,jj,kk)%numc   = ice(ii,jj,kk)%numc   + iceTendecyNumber
                sumICE = sumICE + iceTendecyNumber
                
                liqToIceFrac   = MAX( 0.0, MIN( 1.0, iceTendecyNumber/cloud(ii,jj,kk)%numc ) )
                cloud(ii,jj,kk)%numc = cloud(ii,jj,kk)%numc - iceTendecyNumber
                
                DO ss = 1,ndry
                   ice(ii,jj,kk)%volc(ss) =   ice(ii,jj,kk)%volc(ss) + max(0., cloud(ii,jj,kk)%volc(ss)*liqToIceFrac )
                   cloud(ii,jj,kk)%volc(ss) = cloud(ii,jj,kk)%volc(ss) - max(0., cloud(ii,jj,kk)%volc(ss)*liqToIceFrac )              
                END DO
                ice(ii,jj,kk)%volc(iwa) = ice(ii,jj,kk)%volc(iwa) +   &
                     MAX(0., cloud(ii,jj,kk)%volc(iwa)*liqToIceFrac*spec%rhowa/spec%rhoic)
                cloud(ii,jj,kk)%volc(iwa) = cloud(ii,jj,kk)%volc(iwa) - MAX(0., cloud(ii,jj,kk)%volc(iwa)*liqToIceFrac)
             END IF
          END DO
       END DO
    END DO
    
  END SUBROUTINE ice_fixed_NC
  
  ! ------------------------------------------------------------
  
  SUBROUTINE ice_melt(kproma,kbdim,klev,   &
       ptemp )

    ! UPDATE INDICES
    
    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: kproma,kbdim,klev
    REAL, INTENT(in) :: ptemp(kbdim,klev)

    REAL :: dwet,rhomean

    REAL :: liqvol

    INTEGER :: ii,jj,kk,ss,bb
    INTEGER :: nspec, ndry, iwa, irim

    REAL, PARAMETER :: maxfrac = 0.999
    
    nspec = spec%getNSpec(type="total")
    ndry = spec%getNSpec(type="dry")
    iwa = spec%getIndex("H2O")
    irim = spec%getIndex("rime")
    
    ! JUHA: Should add some real parameterization for the freezing of ice? Instantaneous
    !       at 0 C probably not good at least for larger particles.
    
    DO ii = 1,kbdim
       DO jj = 1,klev
          ! Ice and snow melt when temperature above 273.15 K
          ! => should add the effect of freezing point depression
          IF (ptemp(ii,jj) < 273.15 ) CYCLE
          
          DO kk = 1,nice
             ! Ice => precipitating water?; Juha - not really sure what to assume here...
             ! For now, use precip as target, since likely ice particles that survive to the
             ! point where they are melted, they are probably fairly large....
             IF (ice(ii,jj,kk)%numc < ice(ii,jj,kk)%nlim) CYCLE
             
             ! Find The corresponding precip bin
             CALL ice(ii,jj,kk)%updateDiameter(limit=.TRUE.,type="wet")
             dwet = ice(ii,jj,kk)%dwet
             CALL ice(ii,jj,kk)%updateRhomean()
             rhomean = ice(ii,jj,kk)%rhomean
             bb = getPrecipBin(dwet,rhomean)
             
             DO ss = 1,ndry
                precp(ii,jj,bb)%volc(ss) = precp(ii,jj,bb)%volc(ss) + maxfrac*ice(ii,jj,kk)%volc(ss)
                ice(ii,jj,kk)%volc(ss) = (1.-maxfrac)*ice(ii,jj,kk)%volc(ss)
             END DO

             liqvol = ( ice(ii,jj,kk)%volc(iwa)*spec%rhoic + ice(ii,jj,kk)%volc(irim)*spec%rhori )/spec%rhowa

             precp(ii,jj,bb)%numc = precp(ii,jj,bb)%numc + maxfrac*ice(ii,jj,kk)%numc
             ice(ii,jj,kk)%numc = (1.-maxfrac)*ice(ii,jj,kk)%numc
             
             precp(ii,jj,bb)%volc(iwa) = precp(ii,jj,bb)%volc(iwa) + maxfrac*liqvol                       
             ice(ii,jj,kk)%volc(iwa) = (1.-maxfrac)*ice(ii,jj,kk)%volc(iwa)
             ice(ii,jj,kk)%volc(irim) = (1.-maxfrac)*ice(ii,jj,kk)%volc(irim)                

             ! Secondary ice diagnostics
             ice(ii,jj,kk)%SIP_drfr = (1.-maxfrac)*ice(ii,jj,kk)%SIP_drfr
             ice(ii,jj,kk)%SIP_rmspl = (1.-maxfrac)*ice(ii,jj,kk)%SIP_rmspl             
          END DO
             
       END DO
    END DO
    
  END SUBROUTINE ice_melt
  
  ! ------------------------------------------
  
  INTEGER FUNCTION getIceBin(ldwet)
    REAL, INTENT(in)    :: ldwet
    REAL :: vol
    INTEGER :: ii

    ! Find the ice bin where a freezing liquid droplet of given diameter will presumably belong
    ! -----------------------------------------------------------------------------------------
    vol = pi6*ldwet**3
    ! Correct for the the change in density (AS a first guess, just assume pristine ice density)
    vol = vol*spec%rhowa/spec%rhoic
    ii = COUNT( (ice(1,1,iia:fia)%vlolim < vol) )               
    getIceBin = MIN(MAX(1,ii),nice)

  END FUNCTION getIceBin

  INTEGER FUNCTION getPrecipBin(idwet,idens)
    REAL, INTENT(in) :: idwet   ! This should be given as a spherical effective diameter
    REAL, INTENT(in) :: idens   ! this should be the particle mean density (contribution by pristine and rimed ice)
    REAL :: vol
    INTEGER :: ii

    ! Find the precip bin where a melting ice particle of given diameter and mean density will presumably belong
    ! ------------------------------------------------------------------------------------------------------------
    vol = pi6*idwet**3
    ! Correct for the change in density
    vol = vol*idens/spec%rhowa
    ii = COUNT( (precp(1,1,ira:fra)%vlolim < vol) )
    getPrecipBin = MIN(MAX(1,ii),nprc)
    
  END FUNCTION getPrecipBin


  SUBROUTINE iceNucleation(ii,jj,iice,ndry,iwa,irim,pliq,frac) !Det
    INTEGER, INTENT(in) :: ii,jj,iice
    INTEGER, INTENT(in) :: iwa, irim, ndry
    TYPE(Section), INTENT(inout) :: pliq  ! Liquid particle properties
    REAL, INTENT(in)          :: frac  ! Fraction of nucleated particles from liquid phase
    REAL :: f0, frac2, V_tot, frac_DU, f1
    INTEGER :: ss, idu

   
    idu = spec%getIndex("DU",notFoundValue=0)

    ! DUST VOLUME FRACTION
    frac_DU=pliq%volc(idu)/SUM( pliq%volc(1:ndry))

    ! The old IN nucleated fraction
    f0 = pliq%indef

    frac2=frac
    if(frac_DU <= 0.1) then                         ! JUST TO CHANGE ACTIVATED AMOUNT IN A-BINS
       frac2=frac*frac_DU
    endif
    
    IF (frac2<f0) RETURN

    ! Update frozen fraction
    pliq%indef = frac2

    ! Modify frac2(ii,jj,kk) to get Delta phi, and to correct new additional fraction of frozen particles
    frac2 = (frac2-f0)/(1-f0) 

    ! TOTAL FROZEN VOLUME WITHOUT WATER
    V_tot=frac2*SUM( pliq%volc(1:ndry))

    f1 = frac2
    IF (f1 < -1.e-3 .OR. f1 > 2.) WRITE(*,*) 'update indef wrong', f1
    f1 = MIN( MAX(f1,0.),1.-1.e-20 )  ! the fraction can be a very small number, but still meaningful           
    
    ! Dry aerosol             
    IF(frac_DU <= 0.1 ) then
      IF (lsicenucl%state)  &    ! If mode=2 and state=false, do not produce new ice but just remove the aerosol
            ice(ii,jj,iice)%volc(idu) =    &
            MAX(0., ice(ii,jj,iice)%volc(idu) + V_tot*frac2) ! was frac_DU
      
         pliq%volc(idu) =   &
            MAX(0., pliq%volc(idu)-V_tot*frac2) ! was frac_DU
    ELSE
      IF (lsicenucl%state) &   !  If mode=2 and state=FALSE, do not produce new ice, just remove aerosol/droplets
            ice(ii,jj,iice)%volc(idu) = MAX(0., ice(ii,jj,iice)%volc(idu) + pliq%volc(idu)*frac2)
      pliq%volc(idu) = MAX(0., pliq%volc(idu)*(1.-frac2))
      !DO ss = 1,ndry
      !   IF (lsicenucl%state) &   !  If mode=2 and state=FALSE, do not produce new ice, just remove aerosol/droplets
      !         ice(ii,jj,iice)%volc(ss) = MAX(0., ice(ii,jj,iice)%volc(ss) + pliq%volc(ss)*frac2)
      !   pliq%volc(ss) = MAX(0., pliq%volc(ss)*(1.-frac2))
      !END DO
    END IF
    
    ! Water (total ice)
    IF (ANY(pliq%phase == [1,2])) THEN
       ! Aerosol or cloud droplets -> only pristine ice production
       IF (lsicenucl%state) &   ! If mode=2 and state=FALSE, do not produce new ice, just remove aerosol/droplets
            ice(ii,jj,iice)%volc(iwa) = MAX(0.,ice(ii,jj,iice)%volc(iwa) + pliq%volc(iwa)*frac2*spec%rhowa/spec%rhoic)
       pliq%volc(iwa) = MAX(0., pliq%volc(iwa)*(1.-frac2))
    ELSE IF (pliq%phase == 3) THEN
       ! Precip -> rimed ice
       IF (lsicenucl%state) &   ! Same as above
            ice(ii,jj,iice)%volc(irim) = MAX(0.,ice(ii,jj,iice)%volc(irim) + pliq%volc(iwa)*frac2*spec%rhowa/spec%rhori)
       pliq%volc(iwa) = MAX(0., pliq%volc(iwa)*(1.-frac2))
    END IF
    
    ! Number concentration
    IF (lsicenucl%state) &   ! Same as above
         ice(ii,jj,iice)%numc = MAX(0.,ice(ii,jj,iice)%numc + pliq%numc*frac2)
    pliq%numc = MAX(0.,pliq%numc*(1.-frac2))
    
  END SUBROUTINE iceNucleation

  SUBROUTINE iceNucleation_hom(ii,jj,iice,ndry,iwa,irim,pliq,frac) !Stochastic freezing
    INTEGER, INTENT(in) :: ii,jj,iice
    INTEGER, INTENT(in) :: iwa, irim, ndry
    TYPE(Section), INTENT(inout) :: pliq  ! Liquid particle properties
    REAL, INTENT(in)          :: frac  ! Fraction of nucleated particles from liquid phase
    !REAL :: f0, frac2
    REAL :: f0, f1, nnuc_new, ncur0
    INTEGER :: ss
    
    ! The old IN nucleated fraction

    !IF (frac<f0) CYCLE

    ! Modify frac2(ii,jj,kk) to get Delta phi, and to correct new additional fraction of frozen particles
    !frac2 = (frac-f0)/(1-f0) 

    ncur0 = pliq%numc
    nnuc_new = frac*pliq%numc
   
   ! The old IN nucleated fraction
    f0 = pliq%indef
   
   ! Get the new nucleated fraction
    f1 = nnuc_new + (ncur0 - nnuc_new)*f0
    f1 = f1/ncur0
    IF (f1 < -1.e-3 .OR. f1 > 2.) WRITE(*,*) 'update indef wrong', f1
    f1 = MIN( MAX(f1,0.),1.-1.e-20 )  ! the fraction can be a very small number, but still meaningful           
    pliq%indef = f1
    
    ! Dry aerosol
    DO ss = 1,ndry
       IF (lsicenucl%state) &   !  If mode=2 and state=FALSE, do not produce new ice, just remove aerosol/droplets
            ice(ii,jj,iice)%volc(ss) = MAX(0., ice(ii,jj,iice)%volc(ss) + pliq%volc(ss)*frac)
       pliq%volc(ss) = MAX(0., pliq%volc(ss)*(1.-frac))
    END DO
    
    ! Water (total ice)
    IF (ANY(pliq%phase == [1,2])) THEN
       ! Aerosol or cloud droplets -> only pristine ice production
       IF (lsicenucl%state) &   ! If mode=2 and state=FALSE, do not produce new ice, just remove aerosol/droplets
            ice(ii,jj,iice)%volc(iwa) = MAX(0.,ice(ii,jj,iice)%volc(iwa) + pliq%volc(iwa)*frac*spec%rhowa/spec%rhoic)
       pliq%volc(iwa) = MAX(0., pliq%volc(iwa)*(1.-frac))
    ELSE IF (pliq%phase == 3) THEN
       ! Precip -> rimed ice
       IF (lsicenucl%state) &   ! Same as above
            ice(ii,jj,iice)%volc(irim) = MAX(0.,ice(ii,jj,iice)%volc(irim) + pliq%volc(iwa)*frac*spec%rhowa/spec%rhori)
       pliq%volc(iwa) = MAX(0., pliq%volc(iwa)*(1.-frac))
    END IF
    
    ! Number concentration
    IF (lsicenucl%state) &   ! Same as above
         ice(ii,jj,iice)%numc = MAX(0.,ice(ii,jj,iice)%numc + pliq%numc*frac)
    pliq%numc = MAX(0.,pliq%numc*(1.-frac))
    
  END SUBROUTINE iceNucleation_hom


  SUBROUTINE iceDiagnostics(pliq,fimm,fdep,fhom)
    TYPE(Section), INTENT(in) :: pliq
    REAL, INTENT(in)          :: fimm,fdep,fhom
    
    CALL rateDiag%Ice_hom%accumulate(n=pliq%numc*fhom,    &
                                     v=pliq%volc(:)*fhom)
    CALL rateDiag%Ice_dep%accumulate(n=pliq%numc*fdep,    &
                                     v=pliq%volc(:)*fdep)
    CALL rateDiag%Ice_imm%accumulate(n=pliq%numc*fimm,    &
                                     v=pliq%volc(:)*fimm)
    
  END SUBROUTINE iceDiagnostics
  

END MODULE mo_salsa_cloud_ice_DET
