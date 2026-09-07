MODULE mo_particle_external_properties
  USE mo_submctl, ONLY : pi6, eps, rg, surfw0, grav, spec, pi,mwa,pstand,als,alv,rv
  USE classSection, ONLY : Section
  USE mo_ice_shape, ONLY : getDiameter, t_shape_coeffs
  IMPLICIT NONE

  TYPE t_fit_par
     REAL :: slopelg,aresid,cresid,uresid,nup1,nup2, &
	     ndwn1,ndwn2  
  END TYPE t_fit_par
  
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
    REAL FUNCTION terminal_vel(diam,rhop,rhoa,visc,Cc,flag,shape,dnsp)
      IMPLICIT NONE
      REAL, INTENT(in) :: diam,  &      ! Particle diameter; for ice this should be the spherical equivalent diameter
                          rhop          ! Bulk density of particle
      REAL, INTENT(in) :: rhoa, visc, Cc! Air density, viscosity and Cunningham correction factor
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
            terminal_vel = (diam**2)*(rhop-rhoa)*grav*Cc/(18.*visc)  ![m s-1] eq. 3.21
         ELSE IF (diam.gt.19E-6 .AND. diam.lt.1.07E-3) THEN
            ! Large cloud droplets and small raindrops 19 um <= D < 1.07 mm
            NDa    =  4./3.*diam**3.*grav*(rhop-rhoa)*rhoa/visc**2.; ! Davies number 
            X      =  log(NDa)!  log_e
            Y      =  -0.318657E1 + 0.992696*X -0.153193E-2*X**2 - 0.987059E-3*X**3 & 
                      -0.578878E-3*X**4 +0.855176E-4*X**5 -0.327815E-5*X**6
            Re2    =  Cc*exp(Y)
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
            ! Silvia_Note: This equation agrees with experiments up to 7 mm when the velocity reaches 9.2 m/s at 293 K rho_air = 1.2 kg m-3
         END IF
      ELSE IF (flag==4) THEN   ! Ice
         ! Khvorostyanov and Curry 2002
         ! mp = rhoeff*Vsphere --> Vb*rhop with rhop = rhoeff
         ! rhoeff: mass divided by the volume of a circumscribed sphere whose 
         ! diameter is equal to the particle maximum dimension D=dnsp
         Vb = pi6*dnsp**3     
         Ap = shape%gamma*dnsp**shape%sigma
         X = ( 2. * Vb * (rhop - rhoa) * grav * dnsp**2 ) /  &
              ( Ap * rhoa * visc**2 )
         terminal_vel = kcVt(shape,dnsp,X,visc,rhoa)                 
      END IF
      
    END FUNCTION terminal_vel
    
    
    !-----------------------------------------------------------------------------------------------
    ! This function calculates the cross-sectional area using shape parameters derived from the implemented
    ! Morrison, H., & Milbrandt, J. A. (2015). src/src_shared/mp_ice_shape.f90
    ! Parameterization of Cloud Microphysics Based on the Prediction of Bulk Ice Particle Properties. 
    ! Part I: Scheme Description and Idealized Tests. Journal of the Atmospheric Sciences, 72(1), 287–311.
    ! https://doi.org/https://doi.org/10.1175/JAS-D-14-0065.1
    
    ! cross_sec_area = gamma*D**sigma 
    ! For ice particles gamma and sigma changes linearly between those of 
    ! pristine crystals (given in the runles) 
    ! and spherical ones when the ice rime fraction increases from 0 to 1
    
    REAL FUNCTION cross_sec_area(D,flag,shape) 
	IMPLICIT NONE
	
	REAL, INTENT(in) :: D          ! Particle diameter
        INTEGER, INTENT(in) :: flag    ! Parameter for identifying aerosol (1), cloud droplets (2), precip (3), ice (4)
        TYPE(t_shape_coeffs), INTENT(in), OPTIONAL :: shape ! Shape coefficients needed for ice
        
        IF( ANY(flag == [1,2,3])) THEN    
           ! Aerosol and cloud and rain droplets
           ! diam is dwet and we assume spherical droplets
           cross_sec_area = pi/4 * D**2
        ELSE IF (flag==4) THEN   
           ! Ice   
           ! D for non-spherical ice is defined as the maximum particle length or dimension
           cross_sec_area = shape%gamma*D**shape%sigma
        END IF          
    
    END FUNCTION cross_sec_area
    
    !-----------------------------------------------------------------------------------------------
    ! This function calculates the capacitance of droplets and ice crystals using information in 
    ! Pruppacher, H., & Klett, J. (1997). Microphysics of clouds and precipitation. Springer. 
    ! Equations 13-77 to 13-79
    REAL FUNCTION capacitance(D,flag, aspect_ratio) 
        IMPLICIT NONE
	
	REAL, INTENT(in) :: D          ! Particle diameter for phase <4, nonspherical diameter for ice
        INTEGER, INTENT(in) :: flag    ! Parameter for identifying aerosol (1), cloud droplets (2), precip (3), ice (4)
        REAL, INTENT(in) :: aspect_ratio
        REAL :: eccentricity, c
        
        ! capacitance in the units of D meters 
        
        IF( ANY(flag == [1,2,3])) THEN    
           ! Aerosol and cloud and rain droplets
           ! diam is dwet and we assume spherical droplets
           ! 
           capacitance = D/2
        ELSE IF (flag==4) THEN   
          ! Ice   
          ! D for non-spherical ice is defined as the maximum particle length or dimension
          ! Ice particle aspect ratio  
          ! D for non-spherical ice is defined as the maximum particle length 
          ! D is related to particle projected area then dnsp~a
	  ! aspect_ratio = c/a = polar radius /equatorial radius
	  ! c should correspond to the volume of the spheroid with eq.radius a
	  ! c = Vspheroid / (4/3*pi*a**2)  Vspheroid = (massice/numc)/rhoice 
          ! c = aspect_ratio * D/2  ! a= D/2
           IF (aspect_ratio < 0.999) THEN ! oblate
              eccentricity = SQRT(1-aspect_ratio**2)
              capacitance = D/2 * eccentricity / ASIN(eccentricity)
           ELSE IF (aspect_ratio >=0.999 .AND. aspect_ratio < 1.001) THEN ! sphere
              capacitance = D/2
           ELSE IF (aspect_ratio > 1.001) THEN ! prolate
              eccentricity = SQRT(1-aspect_ratio**(-2))
              capacitance = (D/2 * aspect_ratio) * eccentricity / LOG(aspect_ratio) / (1+eccentricity)
           END IF
        END IF        
           
    END FUNCTION capacitance
    
    !-----------------------------------------------------------------------------------------------
    REAL FUNCTION ventilation_factor(diam,rhop,rhoa,visc,beta,flag,shape,dnsp,zdfh2o,aspect_ratio)
       IMPLICIT NONE
       REAL, INTENT(in) :: diam,  &      ! Particle diameter; for ice this should be the spherical equivalent diameter
                          rhop          ! Bulk density of particle
       REAL, INTENT(in) :: rhoa, visc, beta ! Air density, viscosity and Cunningham correction factor
       INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud droplets (2), precip (3), ice (4)
       TYPE(t_shape_coeffs), INTENT(in), OPTIONAL :: shape ! Shape coefficients needed for ice
       REAL, INTENT(in), OPTIONAL :: dnsp                  ! Maximum diameter of non-spherical ice particle
       REAL, INTENT(in) :: zdfh2o           ! Diffusion coefficient of water vapor in air (m2/s)
       REAL, INTENT(in) :: aspect_ratio
      
      
       REAL :: kvisc, velocity, reynolds_number, schmidt_number, xqi,fv
      
       velocity = terminal_vel(diam,rhop,rhoa,visc,beta,flag,shape,dnsp)
       
       kvisc = visc/rhoa
       schmidt_number = kvisc / zdfh2o
       
       fv = 1.0
       ventilation_factor = 1.0 
       
       IF (ANY(flag == [2,3])) THEN 
           ! Aerosol and cloud and rain droplets
           ! diam is dwet and we assume spherical droplets
           ! ventilation effects are just important for raindrops
           reynolds_number = velocity*diam / kvisc
           xqi = schmidt_number**(1./3.) * SQRT(reynolds_number)  
           ! Pruppacher, H., & Klett, J. (1997). Microphysics of clouds and precipitation. Springer.
           IF (diam <= 60.E-6 .AND. xqi<1.4) THEN 
           	ventilation_factor = 1.00 + 0.108*xqi**2 ! eq.13-60
           ELSE IF ((diam > 60.E-6 .AND. diam <= 1500.E-6) .AND. &
                    (xqi>=1.4 .AND. xqi <=51.4)) THEN
                 ventilation_factor = 0.78 + 0.308*xqi   ! eq. 13-61           
           ELSE
           	 ventilation_factor = 1.0           
           END IF  
           ventilation_factor = MIN(ventilation_factor,16.)           
       ELSE IF (flag==4) THEN ! Ice particles
           reynolds_number = velocity*dnsp / kvisc
           xqi = schmidt_number**(1./3.) * SQRT(reynolds_number)  
           ! Pruppacher, H., & Klett, J. (1997). Microphysics of clouds and precipitation. Springer.
           ! fv ventilation coefficient for idealized snow crystals ~ hexagonal plates
           ! Reynolds number must be defined for dnsp, related to the particle projected_area A=gamma*dnsp**sigma
           ! Welss, J.‐N., Siewert, C., & Seifert, A.(2024). Explicit habit‐prediction in the Lagrangian 
           ! super‐particle ice microphysics model McSnow. 
           ! Journal of Advances in Modeling Earth Systems, 16, e2023MS003805. https://doi.org/10.1029/2023MS003805  
           IF (xqi < 1.0) THEN
              fv = 1. + 0.14*xqi**2 ! PK eq-13-88
              ventilation_factor = fv + 2.8E-3*xqi**1.5 / aspect_ratio ! Welss eq.34
           ELSE
              fv = 0.86 + 0.28*xqi ! PK eq-13-89
              ventilation_factor = fv + 2.8E-2*xqi * aspect_ratio    ! Welss eq.34
           END IF   
           ventilation_factor = MIN(ventilation_factor,7.)       
       END IF
                 
    END FUNCTION ventilation_factor
    
        !-----------------------------------------------------------------------------------------------
    REAL FUNCTION thermal_ventilation_factor(diam,rhop,rhoa,visc,beta,flag,shape,dnsp,zthcond,aspect_ratio,cpm)
       IMPLICIT NONE
       REAL, INTENT(in) :: diam,  &      ! Particle diameter; for ice this should be the spherical equivalent diameter
                          rhop           ! Bulk density of particle
       REAL, INTENT(in) :: rhoa, visc, beta ! Air density, viscosity and Cunningham correction factor
       INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud droplets (2), precip (3), ice (4)
       TYPE(t_shape_coeffs), INTENT(in), OPTIONAL :: shape ! Shape coefficients needed for ice
       REAL, INTENT(in), OPTIONAL :: dnsp                  ! Maximum diameter of non-spherical ice particle
       REAL, INTENT(in) :: zthcond        ! Thermal conductivity of dry air (J/m/K/s)
       REAL, INTENT(in) :: aspect_ratio
       REAL, INTENT(in) :: cpm            ! Heat capacity of moist air (J/kg/K)      
      
       REAL :: kvisc, velocity, reynolds_number, prandtl_number, xqi,fv
      
       velocity = terminal_vel(diam,rhop,rhoa,visc,beta,flag,shape,dnsp)
       
       kvisc = visc/rhoa
       prandtl_number = kvisc*cpm / zthcond
       
       fv = 1.0
       thermal_ventilation_factor = 1.0        
            
       IF (ANY(flag == [2,3])) THEN 
           ! Aerosol and cloud and rain droplets
           ! diam is dwet and we assume spherical droplets
           ! ventilation effects are just important for raindrops
           reynolds_number = velocity*diam / kvisc
           xqi = prandtl_number**(1./3.) * SQRT(reynolds_number)  
           ! Jacobson FAM eq.17.31
           IF (xqi.LE.1.4) THEN 
           	thermal_ventilation_factor = 1.00 + 0.108*xqi**2 
           ELSE IF (xqi.GT.1.4) THEN
                 thermal_ventilation_factor = 0.78 + 0.308*xqi             
           END IF              
           thermal_ventilation_factor = MIN(thermal_ventilation_factor,16.)
       ELSE IF (flag==4) THEN ! Ice particles
           ! Assumption thermal_ventilation_factor follows the same
           ! function type than ventilation_factor for ice crystals 
           ! Jacobson FAM also Ji, W., and P. K. Wang, 1999: J. Atmos. Sci., 56, 829–836
           ! 
           reynolds_number = velocity*dnsp / kvisc
           xqi = prandtl_number**(1./3.) * SQRT(reynolds_number)  
           ! Pruppacher, H., & Klett, J. (1997). Microphysics of clouds and precipitation. Springer.
           ! fv ventilation coefficient for idealized snow crystals ~ nearly spherical (oblate or prolate)
           ! Reynolds number must be defined for dnsp, related to the particle projected_area A=gamma*dnsp**sigma
           ! Welss, J.‐N., Siewert, C., & Seifert, A.(2024). Explicit habit‐prediction in the Lagrangian 
           ! super‐particle ice microphysics model McSnow. 
           ! Journal of Advances in Modeling Earth Systems, 16, e2023MS003805. https://doi.org/10.1029/2023MS003805  
           IF (xqi < 1.0) THEN
              fv = 1. + 0.14*xqi**2 ! PK eq-13-88
              thermal_ventilation_factor = fv + 2.8E-3*xqi**1.5 / aspect_ratio ! Welss eq.34
           ELSE
              fv = 0.86 + 0.28*xqi ! PK eq-13-89
              thermal_ventilation_factor = fv + 2.8E-2*xqi * aspect_ratio    ! Welss eq.34
           END IF  
           thermal_ventilation_factor = MIN(thermal_ventilation_factor,7.)         
       END IF
       
          
    END FUNCTION thermal_ventilation_factor
    
    !-----------------------------------------------------------------------------------------------
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
      REAL :: calcDiamLES, mass_p, mass_r
      
      LOGICAL :: l_sph
      
      ! By default, calculate diameter assuming spherical particles (relevant for ice)
      ! To keep consistency these are the indexes of bc and du 
      ! cnstr%rhobc  => allRho(3)
      ! cnstr%rhodu  => allRho(4)
      
      l_sph = .TRUE.
      IF (PRESENT(sph)) l_sph = sph
      
      calcDiamLES=0.

      IF (numc < 1.e-6) RETURN
            
      IF (flag==4) THEN   ! Ice
         mass_p = SUM(mass(1:ns-1))
         mass_r = mass(ns)
         IF (l_sph) THEN
            ! Spherical equivalent for ice
            ! rhoice = vector containing all dry species plus ice and rimed ice
            !       e.g. [rhooc rhodu rhoic rhori]   
            calcDiamLES = ( SUM(mass(1:ns)/spec%rhoice(1:ns))/numc/pi6 )**(1./3.)
         ELSE
            ! non-spherical ice
            ! Get the effective ice diameter, i.e. the max diameter for non-spherical ice 
            ! getDiameter(mpri,mrim,numc)           
            calcDiamLES = getDiameter( mass_p, mass_r, numc)
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
    
    
    REAL FUNCTION Dveff(temp,press,ssi,dnsp,aspect_ratio,mfp,Dv,xk,cap)
    
    ! This function is built using the following references: 
    ! ------------- INFLUENCES OF GAS-PHASE VAPOR DIFFUSION AND ATTACHMENT KINETICS----------------
    ! Fortran codes to calculate the characteristic and surface supersaturation, 
    ! and the deposition coefficients is also available through Data Commons at
    !https://doi.org/10.26208/s7de-et44. by 
    ! 1. Harrington, J. Y., G. A. Sokolowsky, and H. Morrison, 2021: 
    ! Semianalytic Functions to Calculate the Deposition Coefficients for Ice Crystal Vapor Growth 
    ! in Bin and Bulk Microphysical Models. J. Atmos. Sci., 78, 1735–1752, 
    ! 2. Kärcher, B., Jensen, E. J., Pokrifka, G. F., & Harrington, J. Y. (2023). 
    ! Ice supersaturation variability in cirrus clouds: Role of vertical wind ...
    ! Journal of Geophysical Research: Atmospheres, 128, e2023JD039324. 
    ! 3. Harrington, J. Y., A. Moyle, L. E. Hanson, and H. Morrison, 2019: 
    ! On Calculating Deposition Coefficients and Aspect-Ratio Evolution in 
    ! Approximate Models of Ice Crystal Vapor Growth. J. Atmos. Sci., 76, 1609–1625, 
    !------------------------------------------------------------------------
    	REAL, INTENT(IN) :: temp  ! Absolute temperature in K
    	REAL, INTENT(IN) :: press ! Pressure in Pa
    	REAL, INTENT(IN) :: ssi  ! Supersaturation over ice 0<ssi<1 ssi= rhi-1
    	REAL, INTENT(IN) :: dnsp ! non-spherical diameter
    	REAL, INTENT(IN) :: aspect_ratio  ! aspect ratio of ice crystal
    	REAL, INTENT(IN) :: mfp  ! mean free path of vapour molecules  
    	REAL, INTENT(IN) :: Dv   ! water vapor diffusivity
    	REAL, INTENT(IN) :: xk   ! air thermal conductivity 
    	REAL, INTENT(IN) :: cap  ! ice crystal capacitance
    	
    	TYPE(t_fit_par) :: param
        REAL, PARAMETER :: Runiv = 8.314
        REAL, PARAMETER :: scrat_low=0.01,scrat_hi=100. ! Fig2.in 1.
        INTEGER :: flag
        
        REAL :: ei0,gtp1,gtp2,gtp_stand,vel,deltaT,x,sc,sa!,Dv,xk,
        REAL :: scrit_a,scrit_c,capnew,alen,clen,amn,cmn,dmn,asp_rat_mfp,C3a,C3c,m
        REAL :: scratio_a,slrata_interp,slrata_interp_lim,sldiffa_upratio
        REAL :: scratio_c,slratc_interp,slratc_interp_lim,sldiffc_upratio
        REAL :: sisloc_a, slocal_ap,slocal_diff_a, alpha_a,caplengtha,num1
        REAL :: sisloc_c, slocal_cp,slocal_diff_c, alpha_c,caplengthc, num2
        
          ei0 = esi(press,temp)
          !Dv = 2.11e-5*(temp/273.15)**1.94 * (pstand/press) ! vapor diff. (m^2/s)
	  !xk = 2.3823e-2 +7.1177e-5*(temp-273.15) ! thermal conductivity in W/m/k
      	  gtp1 = rv*temp/(Dv*ei0)
      	  gtp2 = (alv**2.0/(xk*rv*temp**2.0)) - &
                 (alv/(xk*temp))
      	  gtp_stand = 1.0/(gtp1+gtp2)                       ! effective diffusivity        
          vel = sqrt(8.0*Runiv*temp/(pi*mwa))               ! molecular speed(m/s)
          
          ! getting scrit along crystal axis
          ! sa, sc: Characteristic supersaturation as a function of supercooling
          ! describes the supersaturation dependence of surface-kinetic mediated growth
          ! Fig.1/Table 1 in 3.  
	  deltaT = temp-273.15
	  x = max(-70.,min(deltaT,-1.0))
	  IF (deltaT.LT.-30.0) THEN
	     ! with bailey and hallett
	     sc = 3.7955 + 0.10614*x + 0.0075309*x**2
	     sa = sc
	  ELSEIF (deltaT.GE.-30.0.AND.deltaT.LE.-22.0) THEN
	     sc = 753.63 + 105.97 * x + 5.5532 * x**2 + 0.12809 * x**3  &
		  + 0.001103 * x**4
	  ELSEIF (deltaT.GT.-22.0.AND.deltaT.LE.-1.0) THEN
	     sc = 1.1217 + 0.038098 * x - 0.083749 * x**2  &
		  - 0.015734 * x**3 - 0.0010108 * x**4  &
		  - 2.9148e-05 * x**5 - 3.1823e-07 * x**6
	  END IF
	  IF(deltaT.GE.-30.0.AND.deltaT.LE.-22.0) THEN
	     sa =  -0.71057 - 0.14775*x + 0.0042304*x**2
	  ELSEIF (deltaT.gt.-22.0.and.deltaT.le.-15.0) THEN
	     sa = -5.2367 - 1.3184*x - 0.11066*x**2 - 0.0032303*x**3
	  ELSEIF (deltaT.GT.-15.0.AND.deltaT.LE.-10.0) THEN
	     sa = 0.34572 - 0.0093029*x + 0.00030832*x**2 ! fit to Nelson & Night
	  ELSEIF (deltaT.GT.-10.0.AND.deltaT.LE.-1.0) THEN
	     sa = 0.34572 - 0.0093029*x + 0.00030832*x**2 ! Nelson & Knight
	  END IF	     
	  scrit_a = sa/100.
	  scrit_c = sc/100.
	  
	  ! capacitance evaluated capacitance evaluated one mean free path from the surface.
	  flag = 4 
	  !cap = capacitance(dnsp,flag, aspect_ratio) 
	  ! Ice particle aspect ratio  
          ! D for non-spherical ice is defined as the maximum particle length 
          ! D is related to particle projected area then dnsp~a
	  ! aspect_ratio = c/a = polar radius /equatorial radius
	  ! c = Vspheroid / (4/3*pi*a**2)  Vspheroid = (massice/numc)/rhoice 
          ! c = aspect_ratio * D/2  ! a= D/2
          alen = dnsp/2 
      	  clen = aspect_ratio*alen
      	  amn = alen + mfp
          cmn = clen + mfp
          dmn = dnsp + 2*mfp
	  asp_rat_mfp = cmn/amn
	  capnew= capacitance(dmn,flag, asp_rat_mfp) 
	  
	  caplengtha = alen*clen/cap
      	  caplengthc = alen**2/cap
          C3a = (vel*caplengtha*cap/capnew)/(4.*Dv)       
          C3c = (vel*caplengthc*cap/capnew)/(4.*Dv)
          !_____________________________________________________
      	  ! Get M coefficient following 
      	  ! Pokrifka, G. F., et al., 2020: Estimating Surface Attachment Kinetic 
      	  ! and Growth Transition Influences on Vapor-Grown Ice Crystals. 
      	  ! J. Atmos. Sci., 77, 2393–2410, https://doi.org/10.1175/JAS-D-19-0303.1. 
          IF (dnsp.LE.20E-6) THEN
	  	m = 1  !dislocation growth
  	  ELSEIF ((dnsp.GT.20E-6).AND.(dnsp.LT.140E-6)) THEN
  	  	m = 1+10*(dnsp/2-10E-6)/(70E-6-10E-6)!1+14*(dnsp/2-10E-6)/(70E-6-10E-6)
  	  ELSEIF (dnsp.GE.140E-6) THEN
  	        m = 10 !15 ! strong step nucleation
  	  END IF 
      	  !_____________________________________________________
      	  ! Get alpha for a-axis next from parameterization
      
	  sisloc_a = Dv/(gtp_stand*rv*temp/ei0) /(1./(1.+C3a))
	  slocal_diff_a = ssi/sisloc_a
	  
  	  CALL getfittingparams(param,m)
	  scratio_a = min(max(slocal_diff_a/scrit_a,scrat_low),scrat_hi)
	  slrata_interp = 10.**(param%slopelg*alog10(slocal_diff_a/scrit_a))
	  slrata_interp_lim = max(min(slrata_interp,1.0) , &
	                      1.0/(sisloc_a))
	  
	  sldiffa_upratio = max(0.0,slrata_interp-1.) + 1.
	  
	  slrata_interp = max(slrata_interp , &
	                      1.0/(sisloc_a)) &
	             + param%aresid*slrata_interp * &
	   (slrata_interp/(1./sisloc_a))**(param%ndwn1) * &
	   min((slrata_interp/(1./sisloc_a))**(param%ndwn2),1.0) &
	   - param%uresid*slrata_interp * slrata_interp_lim**param%nup1 * &
	   sldiffa_upratio**(param%nup2)
	   
	  slocal_ap = slocal_diff_a/slrata_interp    
	  alpha_a=(slocal_ap/scrit_a)**m * tanh((scrit_a/slocal_ap)**m)
	  alpha_a = max(min(alpha_a,1.0),1.e-2) ! min limit was 1e-6
	  
	  !_____________________________________________________
      	  ! Get alpha for c-axis next from parameterization
      	  sisloc_c = Dv/(gtp_stand*rv*temp/ei0) /(1./(1.+C3c))
          slocal_diff_c = ssi/sisloc_c
          scratio_c = min(max(slocal_diff_c/scrit_c,scrat_low),scrat_hi)
          slratc_interp = 10.**(param%slopelg*alog10(slocal_diff_c/scrit_c))
          slratc_interp_lim = max(min(slratc_interp,1.0) , &
              1.0/(sisloc_c))
          sldiffc_upratio = max(0.0,slratc_interp-1.) + 1.
          slratc_interp = max(slratc_interp, &
              1./sisloc_c)     &
              + param%cresid*slratc_interp * &
              (slratc_interp/(1./sisloc_c))**(param%ndwn1) * &
              min((slratc_interp/(1./sisloc_c))**(param%ndwn2),1.0) &
              - param%uresid*slratc_interp * slratc_interp_lim**param%nup1 * &
              sldiffc_upratio**(param%nup2)
          slocal_cp = slocal_diff_c/slratc_interp
          alpha_c=(slocal_cp/scrit_c)**m * tanh((scrit_c/slocal_cp)**m)
          alpha_c = max(min(alpha_c,1.0),1.e-2)
          !_____________________________________________________      
          ! Effective diffusivity 
          num1 = 4*Dv*cap/(alpha_a*vel*alen*clen)+cap/capnew
          num2 = 4*Dv*cap/(alpha_c*vel*alen**2)+cap/capnew
          Dveff = 2./3.*(Dv/num1)+ 1./3.*(Dv/num2)      	
	END FUNCTION Dveff
   
   !----------------------------------------------------------------------
	! SUBROUTINE GETFITTINGPARAMS
	! This function is built using the following references: 
	    ! Fortran codes to calculate the characteristic and surface supersaturation, 
	    ! and the deposition coefficients is also available through Data Commons at
	    !  https://doi.org/10.26208/s7de-et44. by 
	    ! Harrington, J. Y., G. A. Sokolowsky, and H. Morrison, 2021: 
	    ! Semianalytic Functions to Calculate the Deposition Coefficients for Ice Crystal Vapor Growth 
	    ! in Bin and Bulk Microphysical Models. J. Atmos. Sci., 78, 1735–1752, 
	    ! https://doi.org/10.1175/JAS-D-20-0307.1. 
	    !------------------------------------------------------------------------
	!
	! This subroutine gets the various parameters used in the fits for the
	! surface supersaturation based on the chosen value of M
	!     JYH, Penn State University, May 15 2019
	!
	!  Input: m (parameter in Nelson and Baker (1996) approximation for
	!            the deposition coefficient)
	!
	!  Output: slopelg,aresid,cresid,uresid,nup1,nup2,ndwn1,ndwn2
	!          (Coefficients in surface supersaturation polynomial fit)
	!
	!----------------------------------------------------------------------
	SUBROUTINE getfittingparams(param,m)

	  IMPLICIT NONE
	  REAL, INTENT(IN) :: m ! exponent M in equation for deposition coefficient
	  TYPE(t_fit_par), INTENT(out) :: param
	 
	  param%slopelg = 0.1532 + 0.49078 * m - 0.18002 * m**2 + 0.04165 * m**3 - &
	       0.0061541 * m**4 + 0.00058009 * m**5 - 3.3791e-05 * m**6 + &
	       1.1092e-06 * m**7 - 1.5691e-08 * m**8
	  param%aresid = 0.46192 - 0.093575 * m + 0.012798 * m**2 - &
	       0.0008756 * m**3 + 2.2466e-05 * m**4
	  param%cresid = param%aresid
	  param%uresid = 0.25872 - 0.057957 * m + 0.0079602 * m**2 - &
	       0.00052905 * m**3 + 1.3145e-05 * m**4
	  
	  param%nup1 = -0.69061 + 10.765 * m - 9.1666 * m**2 + 3.8888 * m**3 - &
	       0.91779 * m**4 + 0.12544 * m**5 - 0.0098291 * m**6 +  &
	       0.00040868 * m**7 - 6.9744e-06 * m**8
	  
	  param%nup2 = 5.1239 - 19.204 * m + 15.236 * m**2 - 6.4802 * m**3 +  &
	       1.5151 * m**4 - 0.20317 * m**5 + 0.015601 * m**6 - &
	       0.00063752 * m**7 + 1.0736e-05 * m**8
	  
	  param%ndwn1 = -0.86527 - 0.03202 * m + 0.027342 * m**2 - 0.04378 * m**3 &
	       + 0.014734 * m**4 - 0.0021036 * m**5 + 0.0001431 * m**6 &
	       - 4.317e-06 * m**7 + 3.9242e-08 * m**8
	  
	  param%ndwn2 = -4.692 + 12.47 * m - 9.8936 * m**2 + 4.1881 * m**3 &
	       - 1.0003 * m**4 + 0.14011 * m**5 - 0.011359 * m**6 &
	       + 0.00049122 * m**7 - 8.7303e-06 * m**8      
	  
	END SUBROUTINE getfittingparams	
	
	! From thrm.f90 adapted for vapor pressure only
	! ---------------------------------------------------------------------
	! This function calculates the vapor pressure of water over ice as
	! function of temperature and pressure
	!
	REAL FUNCTION esi(p,t)
	  REAL, INTENT (in) :: p, t
	  REAL, PARAMETER   :: c0 = 0.6114327e+03, c1 = 0.5027041e+02,    &
			       c2 = 0.1875982e+01, c3 = 0.4158303e-01,    &
			       c4 = 0.5992408e-03, c5 = 0.5743775e-05,    &
			       c6 = 0.3566847e-07, c7 = 0.1306802e-09,    &
			       c8 = 0.2152144e-12

	  REAL  :: x
	  x = MIN(MAX(-80.,t-273.16),0.)
	  esi = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

	END FUNCTION esi
    
  END MODULE mo_particle_external_properties
 
