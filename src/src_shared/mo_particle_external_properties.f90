MODULE mo_particle_external_properties
  USE mo_submctl, ONLY : pi6, eps, rg, surfw0, grav, spec
  USE classSection, ONLY : Section
  USE mo_ice_shape, ONLY : getDiameter, t_shape_coeffs
  IMPLICIT NONE

  ! This module contains a collection of function to calculate physical and thermodynamical particle properties,
  ! such as diameters, fall velocities, equilibirium saturation ratios at a droplet surface etc.
    
  CONTAINS
    
    !
    ! This function calculates the terminal velocity of liquid droplets falling in
    ! atmospheric moist air at temperature T, pressure P and density dn
    ! using the method of Beard (1976) employed by Pinsky (1998-2008) in his
    ! paper series about turbulence-enhanced collision using particle trajectory analysis
    ! This function takes into account the nonsphericity of falling raindrops 
    ! Beard, K. v. (1976). Terminal Velocity and Shape of Cloud and Precipitation Drops Aloft.
    ! Journal of Atmospheric Sciences, 33(5), 851–864. 
    ! https://doi.org/https://doi.org/10.1175/1520-0469(1976)033<0851:TVASOC>2.0.CO;2
    ! Pinsky, M. B., Khain, A. P., & Shapiro, M. (2007). Collisions of Cloud Droplets in a Turbulent Flow.
    ! Part IV: Droplet Hydrodynamic Interaction. Journal of the Atmospheric Sciences, 64(7), 2462–2482.
    ! https://doi.org/https://doi.org/10.1175/JAS3952.1
    ! Silvia: 27.05.2023
    ! 
    ! Juha : 
    ! Settling velocities of ice particles are calculated combining a modified version of the P3 scheme
    ! of Morrison and Milbrandt(2015) and Kh  and Curry (2005)
    ! See Ahola et al. (2020) for description of primary ice formation and ice microphysics in
    ! Morrison, H., & Milbrandt, J. A. (2015). 
    ! Parameterization of Cloud Microphysics Based on the Prediction of Bulk Ice Particle Properties. 
    ! Part I: Scheme Description and Idealized Tests. Journal of the Atmospheric Sciences, 72(1), 287–311.
    ! https://doi.org/https://doi.org/10.1175/JAS-D-14-0065.1
    ! Khvorostyanov, V. I.,  Curry, J. A. (2002).
    ! Terminal Velocities of Droplets and Crystals: Power Laws with Continuous Parameters over the Size Spectrum.
    ! Journal of the Atmospheric Sciences, 59(11), 1872–1884.
    ! https://doi.org/10.1175/1520-0469(2002)059<1872:TVODAC>2.0.CO;2
    !
    ! The function requires the surface tension of pure water to account for changes
    ! in droplet shape when large droplets fall
    ! Functions and references are included in the correspondent section
    ! 
    REAL FUNCTION terminal_vel(diam,rhop,rhoa,visc,beta,flag,shape,dnsp)
      IMPLICIT NONE
      REAL, INTENT(in) :: diam,  &      ! Particle diameter; for ice this should be the spherical equivalent diameter
                          rhop          ! Bulk density of particle
      REAL, INTENT(in) :: rhoa, visc, beta ! Air density, viscocity and Cunningham correction factor
      INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud droplets (2), precip (3), ice (4)
      TYPE(t_shape_coeffs), INTENT(in), OPTIONAL :: shape ! Shape coefficients needed for ice
      REAL, INTENT(in), OPTIONAL :: dnsp                  ! Maximum diameter of non-spherical ice particle
      
      ! Constants
      REAL, PARAMETER :: rhoa_ref = 1.225 ! reference air density (kg/m^3)

      REAL :: Vb, Ap   ! Bulk volume, cross sectional area (should revise Ap for nonspherical ice!!!)
      REAL :: Re2,X,Y,NDa,Np,Bo,X3,Y3,Re3, temp
           
      terminal_vel = 0.
      IF( ANY(flag == [1,2,3])) THEN
         ! Aerosol and cloud and rain droplets
         ! Modification to mo_salsa_dynamics 27.05.2023
         ! beta: Cc     = 1+2.51.*knud; knud = 2*mfpair/D 
         
         ! Selection of the flow regime using drop diameter
         IF (diam.le.19E-6) THEN
            !Small cloud droplets 0.5 um <= D < 19 um
            terminal_vel = (diam**2)*(rhop-rhoa)*grav*beta/(18.*visc)  ![m s-1] eq. 3.21
         ELSE IF (diam.gt.19E-6 .AND. diam.lt.1.07E-3) THEN
            ! Large cloud droplets and small raindrops 19 um <= D < 1.07 mm
            NDa    =  4./3.*diam**3.*grav*(rhop-rhoa)*rhoa/visc**2.; ! Davies number 
            X      =  log(NDa)!  log_e
            Y      =  -0.318657E1 + 0.992696*X -0.153193E-2*X**2 - 0.987059E-3*X**3 & 
                      -0.578878E-3*X**4 +0.855176E-4*X**5 -0.327815E-5*X**6
            Re2    =  beta*exp(Y)
            terminal_vel = visc/rhoa/diam*Re2 
         ELSE
            ! Raindrops
            ! sigma =  calcSurfW(temp) You could use temperature dependent surface tension
            ! To keep consistency with other routines a constant value of surfw0=72 mN/m is used
            Np    =  surfw0**3.*rhoa**2./(visc**4.*(rhop-rhoa)*grav) 
            Bo    =  4./3.*(rhop-rhoa)*grav/surfw0*diam**2. !Bond number 
            X3    =  log(Bo*Np**(1./6.)) !log_e
            Y3    =  -0.500015E1 + 0.523778E1*X3 -0.204914E1*X3**2 +0.475294*X3**3  &
                     -0.542819E-1*X3**4 +0.238449E-2*X3**5 
            Re3    =  Np**(1./6.)*exp(Y3) 
            terminal_vel = MIN(visc/rhoa/diam*Re3,9.2) 
            ! R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
            ! Note: this is valid up to 2 mm or 9 m/s (at 1000 mbar), where droplets start to break
            ! Silvia_Note: This equation agrees with experiments up to 4 mm when the velocity reaches 9.2 m/s at 293 K rho_air = 1.2 kg m-3
         END IF       
         ! Previous code overestimates the settling velocity for droplets smaller than 1.2mm and introduces a sharp change at 80um
         !IF (diam<80.0e-6) THEN
            ! Stokes law with Cunningham slip correction factor
          !  terminal_vel = (diam**2)*(rhop-rhoa)*grav*beta/(18.*visc)  ![m s-1]
         !ELSE IF (diam<1.2e-3) THEN
            ! Droplets from 40 um to 0.6 mm: linear dependence on particle radius and a correction for reduced pressure
            !   R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
          !  terminal_vel = 4.e3*diam*SQRT(rhoa_ref/rhoa)
         !ELSE
            ! Droplets larger than 0.6 mm: square root dependence on particle radius and a correction for reduced pressure
            !   R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
            ! Note: this is valid up to 2 mm or 9 m/s (at 1000 mbar), where droplets start to break
          !  terminal_vel = 2.01e2*SQRT( MIN(diam/2.,2.0e-3)*rhoa_ref/rhoa )
          !END IF
      ELSE IF (flag==4) THEN   ! Ice
         ! Khvorostyanov and Curry 2002
         Vb = pi6*diam**3     ! Bulk volume of the particle obtained from spherical equivalent diameter
         Ap = shape%gamma*dnsp**shape%sigma
         X = ( 2. * Vb * (rhop - rhoa) * grav * dnsp**2 ) /  &
              ( Ap * rhoa * visc**2 )
         terminal_vel = kcVt(shape,dnsp,X,visc,rhoa)                 
      END IF
      
    END FUNCTION terminal_vel
    
    !--
    REAL FUNCTION kc1213(X)
      ! Calculate the term needed in 2.12 and 2.13 in Khvorostyanov and Curry 2002
      REAL, INTENT(in) :: X
      REAL, PARAMETER :: c1 = 0.0902 !(KC2002)
      kc1213 = SQRT(1. + c1*SQRT(X))      
    END FUNCTION kc1213
    !--
    REAL FUNCTION kcbre(X)
      ! b_Re from 2.12 in KC2002
      REAL, INTENT(in) :: X
      REAL, PARAMETER :: c1 = 0.0902 !(KC2002)
      kcbre = 0.5 * c1 * SQRT(X)
      kcbre = kcbre / ( kc1213(X) - 1. )
      kcbre = kcbre / kc1213(X)
    END FUNCTION kcbre
    !--
    REAL FUNCTION kcare(X)
      ! a_Re from 2.13 in KC2002
      REAL, INTENT(in) :: X
      REAL, PARAMETER :: delta0 = 9.06, c1 = 0.0902 !(KC2002)
      REAL :: bre
      bre = kcbre(X)
      kcare = 0.25*delta0**2
      kcare = kcare * (kc1213(X) - 1.)**2
      kcare = kcare / (X**bre)
    END FUNCTION kcare
    !--
    REAL FUNCTION kcAv(shape,X,visc,rhoa)
      ! 2.24 from KC2002
      TYPE(t_shape_coeffs), INTENT(in) :: shape
      REAL, INTENT(in) :: X,visc, rhoa
      REAL :: are, bre
      are = kcare(X)
      bre = kcbre(X)      
      kcAv = are * visc**(1.-2.*bre)
      kcAv = kcAv * ( (2.*shape%alpha*grav)/(rhoa*shape%gamma) )**bre      
    END FUNCTION kcAv
    !--
    REAL FUNCTION kcBv(shape,X)
      ! 2.25 from KC2002
      TYPE(t_shape_coeffs), INTENT(in) :: shape
      REAL, INTENT(in) :: X
      REAL :: bre
      bre = kcbre(X)
      kcBv = bre * (shape%beta - shape%sigma + 2.) - 1.
    END FUNCTION kcBv
    !--
    REAL FUNCTION kcVt(shape,D,X,visc,rhoa)
      ! 2.23 from KC2002
      TYPE(t_shape_coeffs), INTENT(in) :: shape
      REAL, INTENT(in) :: D, X, visc, rhoa     
      REAL :: Av, Bv      
      Av = kcAv(shape,X,visc,rhoa)
      Bv = kcBv(shape,X)
      kcVt = Av * D**Bv
    END FUNCTION kcVt
      
    
    !
    ! Function for calculating effective (wet) radius for any particle type
    ! - Aerosol, cloud and rain are spherical
    ! - Snow and ice can be irregular and their densities can be size-dependent
    !
    ! Correct dimension is needed for irregular particles (e.g. ice) for calculating fall speed (deposition and coagulation)
    ! and capacitance (condensation). Otherwise spherical assumed. 
    !
    FUNCTION calcDiamLES(ns,numc,mass,flag,sph)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ns ! Number of species
      INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud droplets (2), precip (3) and ice (4) particle phases
      REAL, INTENT(IN) :: numc, mass(ns)
      LOGICAL, OPTIONAL, INTENT(in) :: sph
      REAL :: calcDiamLES
      
      LOGICAL :: l_sph
      
      ! By default, calculate diameter assuming spherical particles (relevant for ice)
      l_sph = .TRUE.
      IF (PRESENT(sph)) l_sph = sph
      
      calcDiamLES=0.

      IF (numc < 1.e-15) RETURN
            
      IF (flag==4) THEN   ! Ice
         IF (l_sph) THEN
            ! Spherical equivalent for ice
            calcDiamLES = ( SUM(mass(1:ns)/spec%rhoice(1:ns))/numc/pi6 )**(1./3.)
         ELSE
            ! non-spherical ice
            ! Get the effective ice diameter, i.e. the max diameter for non-spherical ice            
            calcDiamLES = getDiameter( SUM(mass(1:ns-1)),mass(ns),numc )
         END IF
      ELSE
         ! Radius from total volume of a spherical particle or aqueous droplet
         calcDiamLES = ( SUM(mass(1:ns)/spec%rholiq(1:ns))/numc/pi6 )**(1./3.)
      ENDIF

    END FUNCTION calcDiamLES

    ! -------------------------------------------------

    !
    ! Function for calculating equilibrium water saturation ratio at droplet surface based on Köhler theory
    !
    REAL FUNCTION calcSweq(part,T)
      TYPE(Section), INTENT(in) :: part ! Any particle
      REAL, INTENT(IN) :: T ! Absolute temperature (K)
      REAL :: dwet
      REAL :: sigma
      REAL :: znw,zns ! Moles of water and soluble material
      REAL :: zvw, zvs, zvtot ! Volume concentrations of water and soluble material and total dry
      INTEGER :: iwa, ndry ! Index for water, number of "dry" species
      INTEGER :: i
      
      iwa = spec%getIndex("H2O")
      ndry = spec%getNSpec(type="dry")

      calcSweq = 0.
      IF (part%numc < part%nlim) RETURN
      
      ! Wet diameter  !! USE THE FUNCTIONS PROVIDED FOR THIS??
      dwet = (SUM(part%volc(:))/part%numc/pi6)**(1./3.)
      
      ! Equilibrium saturation ratio = xw*exp(4*sigma*v_w/(R*T*Dwet))
      
      znw = part%volc(iwa)*spec%rhowa/spec%mwa
      zvw = part%volc(iwa)
      zns = 0.
      zvs = 0.
      zvtot = 0.
      DO i = 1,ndry
         zns = zns + spec%diss(i)*part%volc(i)*spec%rholiq(i)/spec%MM(i)
         zvs = zvs + MIN(1.,spec%diss(i)) * part%volc(i) ! Use "diss" here just to select the soluble species
         zvtot = zvtot + part%volc(i)
      END DO

            
       ! Combine the two cases from original code since they're exactly the same??
      IF (zvw > 1.e-28*part%numc .OR. zvs > 1.e-28*part%numc) THEN
         ! Aqueous droplet OR dry partially soluble particle
         calcSweq = (znw/(zns+znw)) * exp(4.*surfw0*spec%mwa/(rg*T*spec%rhowa*dwet))
      ELSE IF (zvtot-zvs > 1.e-28*part%numc) THEN
         ! Dry insoluble particle
         calcSweq = exp(4.*surfw0*spec%mwa/(rg*T*spec%rhowa*dwet))
      ELSE
         ! Just add eps to avoid divide by zero
         calcSweq = (znw/(eps+zns+znw)) * exp(4.*surfw0*spec%mwa/(rg*T*spec%rhowa*dwet))
      END IF
      
    END FUNCTION calcSweq

    ! -------------------------------------------------
    !
    ! Function for calculating surface tension of pure water as a function temperature
    ! It will allow to consider changes in shape of falling raindrops
    FUNCTION calcSurfW(T) RESULT(sigma)
      REAL, INTENT(IN) :: T ! Absolute temperature (K)
      REAL :: sigma
          ! Raindrops
            ! sigma: surface tension in mN/m
            ! tau= 1-T/Tc; 
            ! T: Absolute Temperature in K
            ! Tc: Critical Temperature in K
            ! B, mu, b: model's parameters
            ! Tc=647.096; %K
            ! B=235.8;    % mN/m
            ! mu=1.256; 
            ! b=-0.625; 
            ! tau=1-T./Tc;
            ! sigma=B.*tau.^mu.*(1+b.*tau); % mN/m
            ! sigma=sigma.*1E-3;  % N/m
            ! The International Association for the Properties of Water and Steam. 2014. 
            ! Revised Release on Surface Tension of Ordinary Water Substance: IAPWS R1-76 (2014).
            ! http://www.iapws.org
            sigma = (235.8*(1-T/647.096)**1.256*(1-0.625*(1-T/647.096)))/1000.
      
    END FUNCTION calcSurfW
    
    
  END MODULE mo_particle_external_properties
  
 
