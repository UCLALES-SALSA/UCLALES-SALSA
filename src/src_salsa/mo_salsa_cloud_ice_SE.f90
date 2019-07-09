MODULE mo_salsa_cloud_ice_SE
  USE mo_salsa_types, ONLY : liquid, ice, precp, rateDiag
  USE mo_submctl, ONLY : nliquid, ira,fra, nprc, iia, fia, nice, pi6, ice_hom, ice_dep, ice_imm, spec,  &
                         boltz, pi, planck, rg, avog, lsFreeTheta, initMinTheta
  USE util, ONLY : erfm1, f_gauss
  USE mo_particle_external_properties, ONLY : calcSweq
  USE classSection, ONLY : Section
  IMPLICIT NONE

  REAL, PARAMETER, PRIVATE :: sqrt2 = SQRT(2.)

  ! Mean and standard deviation for the contact angle distribution
  ! For simplicity just use one set, even though they can be very different
  ! different IN constituents and the fits to measurements in different
  ! nucleation modes can also be quite different. Below default values for DU
  ! in immersion mode (in degrees) from Savre and Ekman 2015. Adjust these from the SALSA
  ! namelist to suit different setups.
  !
  REAL, SAVE :: mean_theta = 132. ! mean contact angle
  REAL, SAVE :: sigma_theta = 20.  ! STD of the contact angle distribution.

  ! This contains the ice nucleation paramterization procedures according
  ! to Savre and Ekman 2015

  
  CONTAINS

  SUBROUTINE ice_nucl_driver(kproma,kbdim,klev,ptemp,prv,prs,prsi,tstep)
    INTEGER, INTENT(in) :: kproma, kbdim, klev
    REAL, INTENT(in) :: ptemp(kbdim,klev), prv(kbdim,klev), prs(kbdim,klev), prsi(kbdim,klev)
    REAL, INTENT(in) :: tstep
    
    REAL :: th00(kbdim,klev,nliquid) 
    REAL :: Seq(kbdim,klev,nliquid), Si(kbdim,klev,nliquid)
    REAL :: f_dep(kbdim,klev,nliquid), f_imm(kbdim,klev,nliquid), f_hom(kbdim,klev,nliquid)
    REAL :: frac(kbdim,klev,nliquid)
    REAL :: Jhom
    
    INTEGER :: ii,jj,kk,is, ibc, idu
    INTEGER :: thind
    REAL :: dwet, dins
    REAL, PARAMETER :: dmin = 1.e-10
    REAL, PARAMETER :: Vmin = pi6*(1.e-10)**3
    REAL, PARAMETER :: tmax_homog = 243.15  ! Juha: Isnt this a bit high?   
    LOGICAL :: nuc_mask(kbdim,klev,nliquid) ! Mask for conditions on different ice nucleation processes    

    
    ! Use IN deficit fraction as a "prognostic" variable. From that, first update theta0. Then calculate nucleation.
    ! Finally, estimate the new total nucleation deficit fraction. This wont be exact, but hopefully close enough.

    ! If no insolubles present, no need to do this (except homogeneous nucleation if it will be implemented)
    IF (spec%Ninsoluble > 0) THEN

       f_dep = 0.; f_imm = 0.; f_hom = 0.

       ibc = spec%getIndex("BC",notFoundValue=0)
       idu = spec%getIndex("DU",notFoundValue=0)
       thind = MERGE(2,1, idu > 0 )

       Seq = 0.
       Si = 0.
       th00 = 0.
       
       DO kk = 1,nliquid
          DO jj = 1,klev
             DO ii = 1,kproma
                Seq(ii,jj,kk) = calcSweq(liquid(ii,jj,kk),ptemp(ii,jj))
                Si(ii,jj,kk) = prv(ii,jj)/prsi(ii,jj)

                ! Update particle diameters for later use
                CALL liquid(ii,jj,kk)%updateDiameter(type="all", limit=.TRUE.)                
                dwet = liquid(ii,jj,kk)%dwet
                dins = liquid(ii,jj,kk)%dins
                
                ! Calculate homogeneous freezing here separately since it does not need contact angle integration
                ! Homogeneous freezing
                IF (dwet-dins > dmin .AND. ptemp(ii,jj) < tmax_homog .AND. ice_hom) THEN
                   CALL J_hf(ptemp(ii,jj),Seq(ii,jj,kk),Jhom)
                   f_hom(ii,jj,kk) = 1. - EXP( -Jhom*pi6*(dwet**3-dins**3)*tstep )
                END IF
                
             END DO
          END DO
       END DO

       ! Update the low limit contact angle for the distribution integration
       CALL low_theta( kproma,kbdim,klev,th00,mean_theta,sigma_theta )

       ! Immersion freezing
       nuc_mask(:,:,:) = ( ( liquid(:,:,:)%phase == 2 .OR. liquid(:,:,:)%phase == 3 ) .AND. &
                           ( liquid(:,:,:)%dins > dmin) .AND. ice_imm )

       CALL gauss_legendre( kproma, kbdim, klev, ptemp, Seq, th00,        &
                            mean_theta, sigma_theta, tstep, nuc_mask, 1, f_imm   )

       ! Deposition freezing
       DO kk = 1,nliquid
          nuc_mask(:,:,kk) = ( ( liquid(:,:,kk)%phase == 1 .AND. prv(:,:)/prs(:,:) < 1.0 ) .AND.  &
                               ( liquid(:,:,kk)%dwet-liquid(:,:,kk)%dins < dmin ) .AND.           &
                               ( liquid(:,:,kk)%dins > dmin ) .AND. ice_dep )
       END DO

       CALL gauss_legendre( kproma, kbdim, klev, ptemp, Si, th00,       &
                            mean_theta, sigma_theta, tstep, nuc_mask, 2, f_dep  )

       CALL iceDiagnostics(kproma,kbdim,klev,liquid,f_imm,f_dep,f_hom)
             
       frac = MAX(0., MIN(0.99,f_imm+f_hom+f_dep-(f_imm+f_dep)*f_hom))
       
       CALL iceNucleation( kproma,kbdim,klev,frac )
       
    END IF ! insoluble
       
  END SUBROUTINE Ice_Nucl_Driver

  ! ----------------------------------------
  
  SUBROUTINE gauss_legendre( kproma, kbdim, klev, Tk, Seq, th00,     &
                             thm, ths, tstep, nuc_mask, mode, frac   )
    
    ! Integrate the ice nucleation across contact angles using the Gauss-Legenre quadrature
    INTEGER, INTENT(in) :: kproma,kbdim,klev
    REAL, INTENT(in) :: Tk(kbdim,klev), Seq(kbdim,klev,nliquid), th00(kbdim,klev,nliquid)  !!! Note for deposition mode, Seq should be saturation w.r.t ice
    REAL, INTENT(in) :: thm, ths, tstep
    LOGICAL, INTENT(in) :: nuc_mask(kbdim,klev,nliquid)
    INTEGER, INTENT(in) :: mode   ! 1 = immersion, 2=deposition
    REAL, INTENT(out) :: frac(kbdim,klev,nliquid)   ! Fraction nucleated
    
    REAL :: thmax, thmin, th1
    
    ! Points and weights for 10th order method
    REAL, DIMENSION(10), PARAMETER  :: x = [ -0.973908239224106,-0.865060566454542,   &
                                             -0.679410688379483,-0.433395397408358,   &
                                             -0.148874339007507, 0.148874339007507,   &
                                             0.433395397408358,  0.679410688379483,   &
                                             0.865060566454542,  0.973908239224106    ]
    REAL, DIMENSION(10), PARAMETER  :: w = [ 0.066667031481273176,0.14945422671662059,  &
                                             0.21908574324601154,0.26926671836765848,   &
                                             0.29552422471242434,0.29552422471242434,   &
                                             0.26926671836765848,0.21908574324601154,   &
                                             0.14945422671662059,0.066667031481273176   ]    
    REAL :: J, ftheta, Jacc, dins

    INTEGER :: i, ii,jj,kk
    INTEGER :: ins
    
    thmax = 180.
    
    frac = 0.
    
    DO kk = 1,nliquid
       DO jj = 1,klev
          DO ii = 1,kproma

             IF ( .NOT. nuc_mask(ii,jj,kk) .OR.                 &
                  Tk(ii,jj) > 273.15 .OR.                       &
                  liquid(ii,jj,kk)%numc < liquid(ii,jj,kk)%nlim ) CYCLE
             
             thmin = th00(ii,jj,kk)            
             Jacc = 0.

             dins = liquid(ii,jj,kk)%dins
                          
             DO i = 1,10
                
                th1 = 0.5*(thmax - thmin)*x(i) + 0.5*(thmax + thmin)

                IF (mode == 1) THEN
                   CALL J_imm(th1, Tk(ii,jj), dins, Seq(ii,jj,kk), J)
                ELSE IF (mode == 2) THEN
                   CALL J_dep(th1, Tk(ii,jj), dins, Seq(ii,jj,kk), J)
                END IF
                   
                IF (th1 > thmin) THEN
                   ftheta = f_gauss(th1,ths,thm)
                ELSE
                   ftheta = 0.
                END IF

                frac(ii,jj,kk) = frac(ii,jj,kk) + ftheta*w(i)*( 1. - EXP( -J*tstep ) )

             END DO

             frac(ii,jj,kk) = 0.5 * (thmax - thmin) * frac(ii,jj,kk) 
                          
          END DO
       END DO
    END DO

  END SUBROUTINE gauss_legendre

  ! ----------------------------------------
  
  SUBROUTINE J_imm(theta,Tk,Din,Seq,J)
    ! Immersion freezing according to Khvorostyanov and Curry 2000    

    IMPLICIT NONE
    REAL, INTENT(in) :: theta,    & ! Contact angle
                        Tk,       & ! temperature in K
                        Din,      & ! Diameter of the IN
                        Seq         ! Equilibrium saturation ratio 
    REAL, INTENT(out) :: J          ! Nucleation rate (per particle per sec)
    
    REAL, PARAMETER :: & 
         C = 1.7e10, &           ! Constant (1.7e11 dyn cm^-2 = 1.7e11*1e-5/1e-4 N m^-2 = 1.7e10 N m^-2)
         rho_ice = 900., &       ! Density of ice (kg/m^3)
         c_1s = 1.e19, &          ! The concentration of water molecules adsorbed on 1 m^2 a surface (#/m^2)
         T0 = 273.15             ! 0 C in Kelvins
    REAL, PARAMETER :: &         ! Case-dependent parameters
         epsi = 0., &            ! Elastic strain produced in ice embryo by the insoluble substrate
         alpha = 0.0             ! Relative area of active sites

    REAL :: costh                ! Cosine of the contact angle    
    REAL :: Tc, act_energy, Lefm, GG, sigma_is, d_g, sf, crit_energy
    REAL :: thrad
    
    J = 0.
    
    ! Must have a core
    IF (Din<1e-10) RETURN
    
    Tc = Tk-T0 ! Temperature in Celsius
    
    ! Activation energy (case-dependent)
    !   Eq. 15 in Jeffery and Austin (1997) and parameters for p=1 bar (Table 2) - used in KC04
    act_energy = rg*Tk*(347./(Tk-177.)-log(4.14/349.))/avog
    !   Khvorostyanov and Sassen (1998) for T < -30 C - used in KC00
    !act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))
    !   Eq. 2 in Li et al. (2013)
    !IF (Tc<=-30.0) THEN
    !    act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))
    !ELSE
    !    act_energy = 5.55*exp(-8.423e-3*Tc+6.384e-4*Tc**2+7.891e-6*Tc**3)/avog*4.1868e3
    !END IF

    ! Surface tension between ice and solution (from KS98)
    sigma_is = 28.e-3+0.25e-3*Tc 

    ! Critical energy of germ formation
    ! a) Ice germ Diameter (eq. 2.6 in KC00)
    Lefm = (79.7+0.708*Tc-2.5e-3*Tc**2)*4.1868e3 
    GG = rg*Tk/spec%mwa/Lefm
    
    IF ( (T0/Tk)*Seq**GG<1.0001 ) RETURN

    d_g = 4.*sigma_is/( rho_ice*Lefm*log((T0/Tk)*Seq**GG)-C*epsi**2) 

    IF (d_g<=1e-10) RETURN
    
    ! b) Shape factor (eq. 2.9 in KC00)
    thrad = theta * pi/180.
    costh = COS(thrad)
    
    sf=het_sf(Din/d_g,costh)
    
    ! c) Critical energy (eq. 2.10 in KC00)
    crit_energy = (pi/3.)*sigma_is*(d_g**2)*sf - (1./4.)*alpha*(1.-costh)*Din**2
    
    ! Eq 2.1 in KC00
    J = (boltz*Tk/planck) * c_1s*pi*(Din**2)*exp( -(act_energy+crit_energy)/(boltz*Tk) ) 
    
  END SUBROUTINE J_imm
 
  ! -----------------------------------------------------------

  SUBROUTINE J_dep(theta,Tk,Din,Seq,J)
    ! Deposition freezing
    
    IMPLICIT NONE
    REAL, INTENT(in) :: theta,         & ! Contact angle
                        Tk,            & ! Temperature in K
                        Din,           & ! Diameter of the insoluble aerosol
                        Seq              ! Equilibrium saturation ratio
    REAL, INTENT(out) :: J
    
    REAL :: Tc, act_energy, sigma_iv, d_g, sf, crit_energy

    REAL, PARAMETER :: & ! Constants
         C = 1.7e10, & ! Constant (1.7e11 dyn cm^-2 = 1.7e11*1e-5/1e-4 N m^-2 = 1.7e10 N m^-2)
         rho_ice = 900., & ! Density of ice (kg/m^3)
         T0 = 273.15 ! 0 C in Kelvins
    REAL, PARAMETER :: & ! Case-dependent parameters
         epsi = 0.    ! Elastic strain produced in ice embryo by the insoluble substrate

    REAL :: thrad, costh   ! Contact angle in radians, cosine of contact angle
    
    ! Must have a core and supersaturation over ice
    IF (Din<1e-10 .OR. Seq<1.0001) RETURN
    
    Tc = Tk-T0 ! Temperature in Celsius
    
    ! Activation energy (case-dependent)
    !   Set to zero in KC00
    act_energy = 0.
    
    ! Critical energy of germ formation
    ! a) Ice germ radius (based on eq. 2.12 in KC00)
    ! Surface tension between ice and vapor (from Ho10)
    sigma_iv = ( (76.1-0.155*Tc) + (28.5+0.25*Tc) )*1.e-3 
    d_g = 4.*sigma_iv/( rg*spec%rhoic/spec%mwa*Tk*log(Seq)-C*epsi**2)  
    IF (d_g<=1e-10) RETURN
    
    ! b) Shape factor (eq. 2.9 in KC00)
    thrad = theta * pi/180.
    costh = COS(thrad)
    sf=het_sf(Din/d_g,costh)
    
    ! c) Critical energy (eq. 2.12 in KC00)
    crit_energy = (pi/3.)*sigma_iv*(d_g**2)*sf ! / 4 for diameter
    
    ! Eq 2.13 in KC00
    !   The pre-exponential factor (kineticc oefficient) is about (1e26 cm^-2)*rn**2
    J = (1./4.)*1.e30*(Din**2)*exp( -(act_energy+crit_energy)/(boltz*Tk)) 
    
  END SUBROUTINE J_dep
  
  ! ---------------------------------------------------

  SUBROUTINE J_hf(Tk,Sw,J)
    ! Homogeneous freezing based on Khvorostyanov and Sassen, Geophys. Res. Lett., 25, 3155-3158, 1998
    
    IMPLICIT NONE
    REAL, intent(in) :: Tk, Sw ! Temperature (K) and water vapor saturation ratio
    REAL, INTENT(out) :: J
    
    REAL :: Tc, act_energy, Lefm, GG, sigma_is, d_g, crit_energy

    REAL, PARAMETER :: & ! Constants
         Nc=5.85e16, & ! The number of water molecules contacting unit area of ice germ (#/m^2) [from KC00]
         rho_ice = 900., & ! Density of ice (kg/m^3)
         T0 = 273.15 ! 0 C in Kelvins
        
    Tc = Tk-T0 ! Temperature in Celsius
    
    ! Activation energy (case-dependent)
    !   Khvorostyanov and Sassen (1998) for T < -30 C
    act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))
    
    ! Critical energy of germ formation
    ! a) Ice germ radius (eq. 9a)
    Lefm = (79.7+0.708*Tc-2.5e-3*Tc**2)*4.1868e3 ! Effective latent heat of melting (eq. 6)
    GG = rg*Tk/Lefm/spec%mwa
    IF ( (T0/Tk)*Sw**GG<1.0001 ) RETURN
    sigma_is = 28.e-3+0.25e-3*Tc ! Surface tension between ice and solution
    d_g = 4.*sigma_is/( spec%rhoic*Lefm*log((T0/Tk)*Sw**GG) )
    IF (d_g<=1e-10) RETURN
    
    ! c) Critical energy (eq. 9b)
    crit_energy = (pi/3.)*sigma_is*d_g**2
    
    ! Eq. 1
    J = 2.0*Nc*(spec%rhowa*boltz*Tk/spec%rhoic/planck)*sqrt(sigma_is/boltz/Tk) * &
        exp( -(crit_energy+act_energy)/(boltz*Tk) )
    
  END SUBROUTINE J_hf
  
  ! -----------------------------------------------------

  REAL FUNCTION het_sf(x,mis)
    ! Calculate shape factor for heterogeneous ice nucleation.
    !   Khvorostyanov and Curry, J. Atmos. Sci., 61, 2676-2691, 2004
    REAL :: x ! x=r_core/r_crit
    REAL :: mis ! Cosine of the contact angle
    REAL :: fii, fpsi
    !
    IF (x>100.) THEN
       ! Problems with numerical accuracy when x>>1; x=100 seems to be good limit for using the limiting value
       het_sf = (mis**3-3*mis+2)/4
    ELSE
       fii = sqrt(1.-2.*mis*x+x**2)
       fpsi = (x-mis)/fii
       het_sf = 0.5* (1. + ( ( 1.-mis*x)/fii)**3 + (2.-3.*fpsi+fpsi**3)*x**3 + 3.*mis*(fpsi-1.)*x**2)
    ENDIF
  END FUNCTION het_sf
  
  ! ------------------------------------

  SUBROUTINE low_theta(kproma,kbdim,klev,th00,thmean,thstd)
    
    INTEGER, INTENT(in) :: kproma,kbdim,klev
    REAL, INTENT(in) :: thmean, thstd
    REAL, INTENT(out) :: th00(kbdim,klev,nliquid) 
    INTEGER :: nb, ii, jj
    REAL :: indef        ! IN deficit ratio from the previous timestep

    th00 = 0.

    DO nb = 1,nliquid
       DO jj = 1,klev
          DO ii = 1,kproma
             indef = liquid(ii,jj,nb)%indef
             th00(ii,jj,nb) = thmean + sqrt2*thstd*erfm1( 2.*indef - 1. )
             th00(ii,jj,nb) = MIN( MAX( th00(ii,jj,nb),0. ),180. )
             
             IF (.NOT. lsFreeTheta%state ) th00(ii,jj,nb) = MAX(initMinTheta,th00(ii,jj,nb))
             
          END DO
       END DO
    END DO
    
  END SUBROUTINE low_theta

  ! ----------------------------------------------------

  SUBROUTINE iceNucleation(kproma,kbdim,klev,frac)
    ! Update the IN deficit fraction
    ! Update the number of ice particles
    INTEGER, INTENT(in) :: kproma, kbdim,klev
    REAL, INTENT(in) :: frac(kbdim,klev,nliquid)

    REAL, PARAMETER :: minVin = pi6*(1.e-8)**3  ! Minimum vol (D=10nm) for the insoluble core to count as IN
    
    REAL :: f0, f1, nnuc_new, ncur0
       
    INTEGER :: ii,jj,kk,bb,ss
    INTEGER :: ndry, iwa, irim
    REAL :: dwet

    ndry = spec%getNSpec(type="dry")
    iwa = spec%getIndex("H2O")
    irim = spec%getIndex("rime")
    
    DO kk = 1,nliquid
       DO jj = 1,klev
          DO ii = 1,kproma

             IF (liquid(ii,jj,kk)%numc < liquid(ii,jj,kk)%nlim) CYCLE  ! Not enough droplets
             IF ( SUM(liquid(ii,jj,kk)%volc(spec%ind_insoluble)) < minVin ) CYCLE ! Not enough IN material
             
             ! The old IN "deficit" fraction
             f0 = liquid(ii,jj,kk)%indef

             ! Current number of potential IN particles
             ncur0 = liquid(ii,jj,kk)%numc             
             nnuc_new = frac(ii,jj,kk)*ncur0
             
             f1 = nnuc_new + (ncur0 + nnuc_new)*f0
             f1 = f1/ncur0
             IF (f1 < 0. .OR. f1 > 1.) WRITE(*,*) 'update indef wrong', f1
             f1 = MIN( MAX(f1,1.e-6),1.-1.e-6 )
             
             liquid(ii,jj,kk)%indef = f1
             
             ! Determine the target ice bin
             CALL liquid(ii,jj,kk)%updateDiameter(type="wet",limit=.TRUE.)
             dwet = liquid(ii,jj,kk)%dwet             
             bb = getIceBin(dwet)

             ! Update the ice bins
             ! Dry aerosol
             DO ss = 1,ndry
                ice(ii,jj,bb)%volc(ss) =    &
                     MAX(0., ice(ii,jj,bb)%volc(ss) + liquid(ii,jj,kk)%volc(ss)*frac(ii,jj,kk))
                liquid(ii,jj,kk)%volc(ss) =   &
                     MAX(0., liquid(ii,jj,kk)%volc(ss)*(1.-frac(ii,jj,kk)))
             END DO
             
             ! Water (total ice)
             IF (ANY(liquid(ii,jj,kk)%phase == [1,2])) THEN
                ! Aerosol or cloud droplets -> only pristine ice production
                ice(ii,jj,bb)%volc(iwa) =   &
                     MAX(0.,ice(ii,jj,bb)%volc(iwa) + liquid(ii,jj,kk)%volc(iwa)*frac(ii,jj,kk)*spec%rhowa/spec%rhoic)
                liquid(ii,jj,kk)%volc(iwa) = MAX(0., liquid(ii,jj,kk)%volc(iwa)*(1.-frac(ii,jj,kk)))
             ELSE IF (liquid(ii,jj,kk)%phase == 3) THEN
                ! Precip -> rimed ice
                ice(ii,jj,bb)%volc(irim) =   &
                     MAX(0.,ice(ii,jj,bb)%volc(irim) + liquid(ii,jj,kk)%volc(iwa)*frac(ii,jj,kk)*spec%rhowa/spec%rhori)
                liquid(ii,jj,kk)%volc(iwa) = MAX(0., liquid(ii,jj,kk)%volc(iwa)*(1.-frac(ii,jj,kk)))
             END IF
             
             ! Number concentration
             ice(ii,jj,bb)%numc = MAX(0.,ice(ii,jj,bb)%numc + liquid(ii,jj,kk)%numc*frac(ii,jj,kk))
             liquid(ii,jj,kk)%numc = MAX(0.,liquid(ii,jj,kk)%numc*(1.-frac(ii,jj,kk)))

             CALL ice(ii,jj,bb)%updateRhomean()
             
          END DO
       END DO
    END DO
             
  END SUBROUTINE iceNucleation

  ! --------------------------------------
  
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

  ! ------------------------------------------------------------
  
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

  ! ---------------------------

  SUBROUTINE iceDiagnostics(kproma,kbdim,klev,pliq,fimm,fdep,fhom)
    INTEGER, INTENT(in) :: kproma,kbdim,klev
    TYPE(Section), INTENT(in) :: pliq(kbdim,klev,nliquid)
    REAL, INTENT(in)          :: fimm(kbdim,klev,nliquid),     &
                                 fdep(kbdim,klev,nliquid),     &
                                 fhom(kbdim,klev,nliquid)
    INTEGER :: kk

    DO kk = 1,nliquid
       CALL rateDiag%Ice_hom%accumulate(n=pliq(1,1,kk)%numc*fhom(1,1,kk),    &
                                        v=pliq(1,1,kk)%volc(:)*fhom(1,1,kk))
       CALL rateDiag%Ice_dep%accumulate(n=pliq(1,1,kk)%numc*fdep(1,1,kk),    &
                                        v=pliq(1,1,kk)%volc(:)*fdep(1,1,kk))
       CALL rateDiag%Ice_imm%accumulate(n=pliq(1,1,kk)%numc*fimm(1,1,kk),    &
                                        v=pliq(1,1,kk)%volc(:)*fimm(1,1,kk))
    END DO
       
  END SUBROUTINE iceDiagnostics


  
  
END MODULE mo_salsa_cloud_ice_SE
