MODULE mo_salsa_cloud_ice
  USE classSection, ONLY : Section
  USE mo_submctl, ONLY : in2a,fn2b, iia, fia, ncld, nprc, nice, nsnw, pi, pi6,    &
       nliquid, nfrozen,                              &
       nlim, prlim, ice_hom, ice_imm, ice_dep, rhowa, &
       rhoic, rhosn, boltz, planck, rg, rd, avog,     &
       fixinc, spec
  USE mo_salsa_types, ONLY : aero, cloud, ice, precp, snow, liquid, frozen
  USE mo_particle_external_properties, ONLY : calcSweq
  USE util, ONLY : calc_correlation, cumlognorm, closest

  IMPLICIT NONE

  CONTAINS

  !***********************************************
  ! Ice nucleation
  ! Sources
  !   Mor05   Morrison et al. (Journal of Aerosol Science, 62, 1665-1677, 2005)
  !   KC00    Khvorostyanov and Curry (Geophysical Research letters, 27, 4081-4084, 2000)
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
    REAL, PARAMETER :: tmax_homog=243.15  ! Juha: Isnt this a bit high?

    REAL, PARAMETER :: dmin = 1.e-10 ! Diameter limit

    REAL :: frac

    INTEGER :: ii,jj,kk,ss,bb
    REAL :: pf_imm, pf_dep, pf_hom, jf
    REAL :: Sw_eq, Si, zvol, ra, rb
    REAL :: ddry, dwet, dins

    LOGICAL :: isdry
    
    INTEGER :: ibc, idu, iwa, nspec, ndry
    REAL :: zinsol
    
    ! Ice nucleation modes
    ! 1) Homogeneous freezing: possible for any aqueous droplet with or without insoluble species (DU or BC)
    ! 2) Immersion freezing: possible for aqueous droplets that have an insoluble core (DU or BC)
    ! 3) Deposition freezing: possible for dry insoluble aerosol at sub-saturated conditions (RH < 100%)
    ! 4) Contact freezing: not implemented, because most INPs are immersed in liquid droplets;
    !    Juha: We could (should?) just extend the range of immersion freezing temperatures to replace this, 
    !          since at least in the boundary layer, there's practically never clean ice nuclei, which are 
    !          needed for this and higher up other processes likely dominate.
    
    nspec = spec%getNSpec(type="wet")
    ndry = spec%getNSpec(type="dry")
    
    ! Mass (volume) array indices of BC, DU and water
    ibc = spec%getIndex("BC",notFoundValue=0)
    idu = spec%getIndex("DU",notFoundValue=0)
    iwa = spec%getIndex("H2O")


    DO kk = 1,nice
       DO jj = 1,klev
          DO ii = 1,kproma
             IF ( ANY([ice(ii,jj,kk)%volc(iwa),ice(ii,jj,kk)%vrime] > 1.e-30) .AND.  &
                  ice(ii,jj,kk)%numc > ice(ii,jj,kk)%nlim                     .AND.  &
                  ice(ii,jj,kk)%volc(iwa) < 0.5*ice(ii,jj,kk)%vrime                      )  &
                  WRITE(*,*) 'ENNEN',ice(ii,jj,kk)%volc(iwa), ice(ii,jj,kk)%vrime
          END DO
       END DO
    END DO

    ! Loop over liquid phase bins
    DO kk = 1,nliquid 
       DO ii = 1,kbdim
          DO jj = 1,klev
             IF (ptemp(ii,jj) > 273.15) CYCLE
             IF (liquid(ii,jj,kk)%numc < liquid(ii,jj,kk)%nlim) CYCLE

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
             
             ! Immersion freezing (similar for all categories)
             pf_imm = 0.
             IF (dins > dmin .AND. ice_imm) THEN
                jf = calc_Jhet(dins,ptemp(ii,jj),Sw_eq)
                pf_imm = 1. - EXP( -jf*ptstep )
             END IF

             ! Deposition freezing
             pf_dep = 0.
             IF (dins > dmin .AND. dwet-dins < dmin .AND. prv(ii,jj)/prs(ii,jj)<1.0 .AND. &
                 liquid(ii,jj,kk)%phase == 1 .AND. ice_dep                      ) THEN
                Si = prv(ii,jj)/prsi(ii,jj) ! Water vapor saturation ratio over ice
                jf = calc_Jdep(dins,ptemp(ii,jj),Si)
                pf_dep = 1. - exp( -jf*ptstep )

             END IF
             
             ! Homogeneous freezing
             pf_hom = 0.
             IF (dwet-dins > dmin .AND. ptemp(ii,jj) < tmax_homog .AND. ice_hom) THEN
                jf = calc_Jhf(ptemp(ii,jj),Sw_eq)
                pf_hom = 1. - exp( -jf*pi6*(dwet**3 - dins**3)*ptstep )
             END IF
             
             frac = MAX(0., MIN(1.,pf_imm+pf_hom+pf_dep-(pf_imm+pf_dep)*pf_hom))

             ! Determine the parallel ice bin
             bb = getIceBin(kk,ddry,liquid(ii,jj,kk)%phase)

             CALL iceNucleation(ii,jj,bb,ndry,iwa,liquid(ii,jj,kk),frac)

             IF (ice(ii,jj,bb)%numc > ice(ii,jj,bb)%nlim) &
                  CALL ice(ii,jj,bb)%updateRhomean()

          END DO !ii
       END DO !jj
    END DO

    DO kk = 1,nice
       DO jj = 1,klev
          DO ii = 1,kproma
             IF ( ANY([ice(ii,jj,kk)%volc(iwa),ice(ii,jj,kk)%vrime] > 1.e-30) .AND.  &
                  ice(ii,jj,kk)%numc > ice(ii,jj,kk)%nlim                     .AND.  &
                  ice(ii,jj,kk)%volc(iwa) < 0.5*ice(ii,jj,kk)%vrime                      )  &
                  WRITE(*,*) 'JALKEEN',ice(ii,jj,kk)%volc(iwa), ice(ii,jj,kk)%vrime
          END DO
       END DO
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
         c_1s = 1e19, & ! The concentration of water molecules adsorbed on 1 m^2 a surface (#/m^2)
         T0 = 273.15 ! 0 C in Kelvins
    REAL, PARAMETER :: & ! Case-dependent parameters
         epsi = 0., & ! Elastic strain produced in ice embryo by the insoluble substrate
         alpha = 0.0, & ! Relative area of active sites
         mis = 0.5 ! Cosine of the contact angle
    
    calc_Jhet = 0.
    
    ! Must have a core
    IF (rn<1e-10) RETURN
    
    Tc = temp-T0 ! Temperature in Celsius
    
    ! Activation energy (case-dependent)
    !   Eq. 15 in Jeffery and Austin (1997) and parameters for p=1 bar (Table 2) - used in KC04
    act_energy = rg*temp*(347./(temp-177.)-log(4.14/349.))/avog
    !   Khvorostyanov and Sassen (1998) for T < -30 C - used in KC00
    !act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))
    !   Eq. 2 in Li et al. (2013)
    !IF (Tc<=-30.0) THEN
    !    act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))
    !ELSE
    !    act_energy = 5.55*exp(-8.423e-3*Tc+6.384e-4*Tc**2+7.891e-6*Tc**3)/avog*4.1868e3
    !END IF
    
    ! Critical energy of germ formation
    ! a) Ice germ radius (eq. 2.6 in KC00)
    Lefm = (79.7+0.708*Tc-2.5e-3*Tc**2)*4.1868e3 ! Effective latent heat of melting (eq. 6 in KS98)
    GG = rg*temp/spec%mwa/Lefm ! Eq 2.7 in KC00
    IF ( (T0/temp)*Sw**GG<1.0001 ) RETURN
    sigma_is = 28.e-3+0.25e-3*Tc ! Surface tension between ice and solution (from KS98)
    r_g = 4.*sigma_is/( rho_ice*Lefm*log((T0/temp)*Sw**GG)-C*epsi**2) !*2 for diameter
    IF (r_g<=1e-10) RETURN
    ! b) Shape factor (eq. 2.9 in KC00)
    sf=calc_het_sf(rn/r_g,mis)
    ! c) Critical energy (eq. 2.10 in KC00)
    crit_energy = (pi/3.)*sigma_is*(r_g**2)*sf-(1./4.)*alpha*(1.-mis)*rn**2   ! /4 for diameter
    
    ! Eq 2.1 in KC00
    calc_Jhet= boltz*temp/planck*c_1s*pi*rn**2*exp((-act_energy-crit_energy)/(boltz*temp)) ! /4 for diameter
    
  END FUNCTION calc_Jhet
  
  
  REAL FUNCTION calc_Jdep(rn,temp,Si)
    ! The rate of germ formation (#/s) through deposition freezing following
    ! Khvorostyanov and Curry, Geophys. Res. Lett., 27, 4081-4084, 2000 [KC00]
    ! Additional parameters from
    ! Hoose et al., J. Atmos. Sci., 67, 2483-2503, 2010 [Ho10]

    ! Changed all radiuses to diameters
    
    IMPLICIT NONE
    REAL, INTENT(in) :: rn,temp,Si
    
    REAL :: Tc, act_energy, sigma_iv, r_g, sf, crit_energy
    REAL, PARAMETER :: & ! Constants
         C = 1.7e10, & ! Constant (1.7e11 dyn cm^-2 = 1.7e11*1e-5/1e-4 N m^-2 = 1.7e10 N m^-2)
         rho_ice = 900., & ! Density of ice (kg/m^3)
         T0 = 273.15 ! 0 C in Kelvins
    REAL, PARAMETER :: & ! Case-dependent parameters
         epsi = 0., & ! Elastic strain produced in ice embryo by the insoluble substrate
         mis = 0.5 ! Cosine of the contact angle
    
    calc_Jdep = 0.
    
    ! Must have a core and supersaturation over ice
    IF (rn<1e-10 .OR. Si<1.0001) RETURN
    
    Tc = temp-T0 ! Temperature in Celsius
    
    ! Activation energy (case-dependent)
    !   Set to zero in KC00
    act_energy = 0.
    
    ! Critical energy of germ formation
    ! a) Ice germ radius (based on eq. 2.12 in KC00)
    sigma_iv = ( (76.1-0.155*Tc) + (28.5+0.25*Tc) )*1e-3 ! Surface tension between ice and vapor (from Ho10)
    r_g = 4.*sigma_iv/( rg*rho_ice/spec%mwa*temp*log(Si)-C*epsi**2)   ! R_v=R*rho_ice/M_ice    *2 for diameter
    IF (r_g<=1e-10) RETURN
    ! b) Shape factor (eq. 2.9 in KC00)
    sf=calc_het_sf(rn/r_g,mis)
    ! c) Critical energy (eq. 2.12 in KC00)
    crit_energy = (pi/3.)*sigma_iv*r_g**2*sf ! / 4 for diameter
    
    ! Eq 2.13 in KC00
    !   The pre-exponential factor (kineticc oefficient) is about (1e26 cm^-2)*rn**2
    calc_Jdep = (1./4.)*1e30*rn**2*exp((-act_energy-crit_energy)/(boltz*temp))  ! / 4 for diameter
    
  END FUNCTION calc_Jdep
  
  REAL FUNCTION calc_Jhf(temp,Sw)
    ! Homogeneous freezing based on Khvorostyanov and Sassen, Geophys. Res. Lett., 25, 3155-3158, 1998 [KS98]

    IMPLICIT NONE
    REAL, intent(in) :: temp, Sw ! Temperature (K) and water vapor saturation ratio
    
    REAL :: Tc, act_energy, Lefm, GG, sigma_is, r_g, crit_energy
    REAL, PARAMETER :: & ! Constants
         Nc=5.85e16, & ! The number of water molecules contacting unit area of ice germ (#/m^2) [from KC00]
         rho_ice = 900., & ! Density of ice (kg/m^3)
         T0 = 273.15 ! 0 C in Kelvins
    
    calc_Jhf = 0.
    
    Tc = temp-T0 ! Temperature in Celsius
    
    ! Activation energy (case-dependent)
    !   Khvorostyanov and Sassen (1998) for T < -30 C
    act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))
    
    ! Critical energy of germ formation
    ! a) Ice germ radius (eq. 9a)
    Lefm = (79.7+0.708*Tc-2.5e-3*Tc**2)*4.1868e3 ! Effective latent heat of melting (eq. 6)
    GG = rg*temp/Lefm/spec%mwa
    IF ( (T0/temp)*Sw**GG<1.0001 ) RETURN
    sigma_is = 28e-3+0.25e-3*Tc ! Surface tension between ice and solution
    r_g = 2.*sigma_is/( rho_ice*Lefm*log((T0/temp)*Sw**GG) )
    IF (r_g<=1e-10) RETURN
    ! c) Critical energy (eq. 9b)
    crit_energy = 4.*pi/3.*sigma_is*r_g**2
    
    ! Eq. 1
    calc_Jhf = 2.0*Nc*(spec%rhowa*boltz*temp/rho_ice/planck)*sqrt(sigma_is/boltz/temp)*exp((-crit_energy-act_energy)/(boltz*temp))
    
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
  
  !***********************************************
  !
  ! Prescribed ice number concentration for cloudy regions (rc>0.001 g/kg) where ice
  ! supersaturation is at least 5%. Ice number concentration is increased towards the
  ! target concentration (fixinc, #/kg) by converting the largest cloud droplets to ice.
  !
  !***********************************************
  SUBROUTINE ice_fixed_NC(kproma, kbdim,  klev,   &
       ptemp, ppres, prv,  prsi     )
    

    
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
    
    INTEGER :: iwa, nspec
    
    iwa = spec%getIndex("H2O")
    nspec = spec%getNSpec()
    
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
                
                DO ss = 1,nspec-1
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
        
    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: kproma,kbdim,klev
    REAL, INTENT(in) :: ptemp(kbdim,klev)
    
    INTEGER :: ii,jj,kk,ss
    
    INTEGER :: nspec, iwa
    
    nspec = spec%getNSpec()
    iwa = spec%getIndex("H2O")
    
    ! JUHA: Should add some real parameterization for the freezing of ice? Instantaneous
    !       at 0 C probably not good at least for larger particles.
    
    DO ii = 1,kbdim
       DO jj = 1,klev
          ! Ice and snow melt when temperature above 273.15 K
          ! => should add the effect of freezing point depression
          IF (ptemp(ii,jj) <= 273.15 ) CYCLE
          
          DO kk = 1,nice
             ! Ice => cloud water
             IF (ice(ii,jj,kk)%numc<prlim) CYCLE
             DO ss = 1,nspec-1
                cloud(ii,jj,kk)%volc(ss) = cloud(ii,jj,kk)%volc(ss) + 0.999*ice(ii,jj,kk)%volc(ss)
                ice(ii,jj,kk)%volc(ss) = 0.001*ice(ii,jj,kk)%volc(ss)
             END DO
             ! Water
             cloud(ii,jj,kk)%volc(iwa) = cloud(ii,jj,kk)%volc(iwa) + 0.999*ice(ii,jj,kk)%volc(iwa)*ice(ii,jj,kk)%rhomean/spec%rhowa
             ice(ii,jj,kk)%volc(iwa) = 0.001*ice(ii,jj,kk)%volc(iwa)
             ice(ii,jj,kk)%vrime = 0.001*ice(ii,jj,kk)%vrime
             
             cloud(ii,jj,kk)%numc = cloud(ii,jj,kk)%numc + 0.999*ice(ii,jj,kk)%numc
             ice(ii,jj,kk)%numc = 0.001*ice(ii,jj,kk)%numc
          END DO
          
          DO kk =1,nsnw
             ! Snow => precipitation (bin 1)
             IF (snow(ii,jj,kk)%numc<prlim) CYCLE
             DO ss = 1,nspec-1
                precp(ii,jj,1)%volc(ss) = precp(ii,jj,1)%volc(ss) + snow(ii,jj,kk)%volc(ss)
                snow(ii,jj,kk)%volc(ss) = 0.
             END DO
             ! Water
             precp(ii,jj,1)%volc(iwa) = precp(ii,jj,1)%volc(iwa) + snow(ii,jj,kk)%volc(iwa)*spec%rhosn/spec%rhowa
             snow(ii,jj,kk)%volc(iwa) = 0.
             
             precp(ii,jj,1)%numc = precp(ii,jj,1)%numc + snow(ii,jj,kk)%numc
             snow(ii,jj,kk)%numc = 0.
          END DO
       END DO
    END DO
    
  END SUBROUTINE ice_melt
  
  SUBROUTINE autosnow(kproma,kbdim,klev)
    !
    ! Uses a more straightforward method for converting cloud droplets to drizzle.
    ! Assume a lognormal cloud droplet distribution for each bin. Sigma_g is an adjustable
    ! parameter and is set to 1.2 by default
    !
    

    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: kproma,kbdim,klev
    
    REAL :: Vrem, Nrem, Vtot, Ntot
    REAL :: dvg,dg
    
    REAL, PARAMETER :: zd0 = 250.e-6  ! Adjustable
    REAL, PARAMETER :: sigmag = 1.2   ! Adjustable
    !REAL, PARAMETER :: max_rate_autoc=1.0e10 ! Maximum autoconversion rate (#/m^3/s)
    
    INTEGER :: ii,jj,cc,ss
    
    INTEGER :: nspec, iwa
    
    nspec = spec%getNSpec()
    iwa = spec%getIndex("H2O")
    
    ! Find the ice particle bins where the mean droplet diameter is above 250 um
    ! Do some fitting...
    DO jj = 1,klev
       DO ii = 1,kbdim
          DO cc = 1,nice
             
             Ntot = ice(ii,jj,cc)%numc
             Vtot = SUM(ice(ii,jj,cc)%volc(1:nspec))
             
             IF ( Ntot > prlim .AND. Vtot > 0. ) THEN
                ! Volume geometric mean diameter
                dvg = ((Vtot/Ntot/pi6)**(1./3.))*EXP( (3.*LOG(sigmag)**2)/2. )
                dg = dvg*EXP( -3.*LOG(sigmag)**2 )
                
                Vrem = Max(0., Vtot*( 1. - cumlognorm(dvg,sigmag,zd0) ) )
                Nrem = Max(0., Ntot*( 1. - cumlognorm(dg,sigmag,zd0) )  )
                
                IF ( Vrem > 0. .AND. Nrem > prlim) THEN
                   
                   ! Put the mass and number to the first snow bin and remover from ice
                   DO ss = 1,nspec-1
                      snow(ii,jj,1)%volc(ss) = snow(ii,jj,1)%volc(ss) + ice(ii,jj,cc)%volc(ss)*(Nrem/Ntot)
                      ice(ii,jj,cc)%volc(ss) = ice(ii,jj,cc)%volc(ss)*(1. - (Nrem/Ntot))
                   END DO
                   ! From ice to snow volume
                   snow(ii,jj,1)%volc(iwa) = max(0., snow(ii,jj,1)%volc(iwa) + &
                        ice(ii,jj,cc)%volc(iwa)*(Vrem/Vtot)*spec%rhoic/spec%rhosn)
                   ice(ii,jj,cc)%volc(iwa) = max(0., ice(ii,jj,cc)%volc(iwa)*(1. - (Vrem/Vtot)))
                   
                   snow(ii,jj,1)%numc = snow(ii,jj,1)%numc + Nrem
                   ice(ii,jj,cc)%numc = ice(ii,jj,cc)%numc - Nrem
                   
                END IF ! Nrem Vrem
                
             END IF ! Ntot Vtot
             
          END DO ! cc
       END DO ! ii
    END DO ! jj
    
  END SUBROUTINE autosnow
  
  INTEGER FUNCTION getIceBin(lbin,lddry,phase)
    INTEGER, INTENT(in) :: lbin
    REAL, INTENT(in)    :: lddry
    INTEGER, INTENT(in) :: phase
    
    REAL :: vol
    INTEGER :: ii
    
    vol = pi6*lddry**3
    
    IF ( ANY(phase == [1,2]) ) THEN
       ASSOCIATE( hilim => ice(1,1,iia%cur:fia%cur)%vhilim )
         
         ii = COUNT( (vol < hilim) )
         
       END ASSOCIATE
       
    ELSE IF (phase == 3) THEN
       
       ! Move to ice bin with closest matching dry volume. Ain't perfect but
       ! the bin update subroutine in SALSA will take care of the rest.
       
       ASSOCIATE(dmid => ice(1,1,iia%cur:fia%cur)%dmid)
         ! 1) Find the closest matching bin dry diameter        
         ii = closest(dmid,lddry)
         
         ! Keeping it simple, put everythin just to a-bins. In the near future, the
         ! a/b bins for ice are going to be deprecated anyway.
       END ASSOCIATE
       
    END IF
    
    getIceBin = MAX(1,ii)

  END FUNCTION getIceBin
    
  SUBROUTINE iceNucleation(ii,jj,iice,ndry,iwa,pliq,frac)
    INTEGER, INTENT(in) :: ii,jj,iice
    INTEGER, INTENT(in) :: iwa, ndry
    TYPE(Section), INTENT(inout) :: pliq  ! Liquid particle properties
    REAL, INTENT(in)          :: frac  ! Fraction of nucleated particles from liquid phase
    
    INTEGER :: ss
    
    ! Dry aerosol
    DO ss = 1,ndry
       ice(ii,jj,iice)%volc(ss) = MAX(0., ice(ii,jj,iice)%volc(ss) + pliq%volc(ss)*frac)
       pliq%volc(ss) = MAX(0., pliq%volc(ss)*(1.-frac))
    END DO
    
    ! Water (total ice)
    IF (ANY(pliq%phase == [1,2])) THEN
       ! Aerosol or cloud droplets -> only pristine ice production
       ice(ii,jj,iice)%volc(iwa) = MAX(0.,ice(ii,jj,iice)%volc(iwa) + pliq%volc(iwa)*frac*spec%rhowa/spec%rhoic)
       pliq%volc(iwa) = MAX(0., pliq%volc(iwa)*(1.-frac))
    ELSE IF (pliq%phase == 3) THEN
       ! Precip -> rimed ice
       ice(ii,jj,iice)%volc(iwa) = MAX(0.,ice(ii,jj,iice)%volc(iwa) + pliq%volc(iwa)*frac*spec%rhowa/spec%rhori)
       ice(ii,jj,iice)%vrime = MAX(0.,ice(ii,jj,iice)%vrime + pliq%volc(iwa)*frac*spec%rhowa/spec%rhori)
       pliq%volc(iwa) = MAX(0., pliq%volc(iwa)*(1.-frac))
    END IF
    
    ! Number concentration
    ice(ii,jj,iice)%numc = MAX(0.,ice(ii,jj,iice)%numc + pliq%numc*frac)
    pliq%numc = MAX(0.,pliq%numc*(1.-frac))
    
  END SUBROUTINE iceNucleation
  

END MODULE mo_salsa_cloud_ice
