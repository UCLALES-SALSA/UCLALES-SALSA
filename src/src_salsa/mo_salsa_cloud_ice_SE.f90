MODULE mo_salsa_cloud_ice_SE
  USE mo_salsa_types, ONLY : liquid, ice
  IMPLICIT NONE

  ! This contains the ice nucleation paramterization procedures according
  ! to Savre and Ekman 2015

  
  CONTAINS

  SUBROUTINE IceNucleation
    
    ! Get/update the low theta angle, needed for each ice nucleating species in each bin
    ! low_theta is a tracer in aerosol and cloud bins for each insoluble species
    ! It is initialized as the theta_mean + sigma*sqrt(2)*erf(-1). This is used to
    ! compute ice nucleation in from each bin. low_theta is then updated to be
    ! MAX(low_theta_old, new value) where new_value is obtained using the Nice_new/Ntot for the current bin
    

    
    ! Iterate the immersion and deposition nucleation modes in the legendre polynomial integration

    
  END SUBROUTINE IceNucleation



    
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
    
    calc_Jhet = 0.
    
    ! Must have a core
    IF (rn<1e-10) RETURN
    
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
    
    IF ( (T0/temp)*Sw**GG<1.0001 ) RETURN

    d_g = 4.*sigma_is/( rho_ice*Lefm*log((T0/temp)*Sw**GG)-C*epsi**2) 

    IF (r_g<=1e-10) RETURN
    
    ! b) Shape factor (eq. 2.9 in KC00)
    thrad = theta * pi/180.
    costh = COS(thrad)

    sf=calc_het_sf(Din/d_g,costh)
    
    ! c) Critical energy (eq. 2.10 in KC00)
    crit_energy = (pi/3.)*sigma_is*(d_g**2)*sf - (1./4.)*alpha*(1.-costh)*Din**2
    
    ! Eq 2.1 in KC00
    J = (boltz*Tk/planck) * c_1s*pi*(Din**2)*exp( -(act_energy+crit_energy)/(boltz*temp) ) 
    
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
    sf=calc_het_sf(Din/d_g,costh)
    
    ! c) Critical energy (eq. 2.12 in KC00)
    crit_energy = (pi/3.)*sigma_iv*(d_g**2)*sf ! / 4 for diameter
    
    ! Eq 2.13 in KC00
    !   The pre-exponential factor (kineticc oefficient) is about (1e26 cm^-2)*rn**2
    calc_Jdep = (1./4.)*1.e30*(Din**2)*exp( -(act_energy+crit_energy)/(boltz*Tk)) 
    
  END FUNCTION calc_Jdep
  
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
    d_g = 4.*sigma_is/( spec%rhoic*Lefm*log((T0/Tk)*Seq**GG) )
    IF (d_g<=1e-10) RETURN
    
    ! c) Critical energy (eq. 9b)
    crit_energy = (pi/3.)*sigma_is*d_g**2
    
    ! Eq. 1
    calc_Jhf = 2.0*Nc*(spec%rhowa*boltz*Tk/spec%rhoic/planck)*sqrt(sigma_is/boltz/Tk)*exp( -(crit_energy+act_energy)/(boltz*Tk) )
    
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

  SUBROUTINE low_theta()

    INTEGER :: ni, nins, nb, ii, jj
    REAL :: th00_old

    ! Number of warm phase bins and number of IN species
    nb = SIZE(liquid,DIM=3) 
    nins = spec%Ninsoluble

    ! Loop over insoluble species
    IF (nins > 0) RETURN

    DO jj = 1,klev
       DO ii = 1,kproma
          DO ni = 1,nins
             th00_old = liquid(ii,jj,ni)%thetcn
             
          END DO
       END DO
    END DO
    

  END SUBROUTINE low_theta

  

END MODULE mo_salsa_cloud_ice_SE
