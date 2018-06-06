MODULE mo_particle_external_properties
  USE mo_submctl, ONLY : pi6, eps, rg, surfw0, grav, spec
  USE classSection, ONLY : Section
  IMPLICIT NONE

  ! This module contains a collection of function to calculate physical and thermodynamical particle properties,
  ! such as diameters, fall velocities, equilibirium saturation ratios at a droplet surface etc.
  
  
  CONTAINS
  
    !
    ! Function for calculating terminal velocities for different particle types and size ranges.
    !     Tomi Raatikainen (2.5.2017)
    !     - Changed from radius to diameter since ~the rest of the model as well as the calculations below take diameter anyway! -Juha
    REAL FUNCTION terminal_vel(diam,rhop,rhoa,visc,beta,flag)
      IMPLICIT NONE
      REAL, INTENT(in) :: diam, rhop ! Particle diameter and density
      REAL, INTENT(in) :: rhoa, visc, beta ! Air density, viscocity and Cunningham correction factor
      INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud droplets (2), precip (3), ice (4) and snow(5)
      ! Constants
      REAL, PARAMETER :: rhoa_ref = 1.225 ! reference air density (kg/m^3)
      

      terminal_vel = 0.
      IF( ANY(flag == [1,2,3])) THEN
         ! Aerosol and cloud and rain droplets
         IF (diam<80.0e-6) THEN
            ! Stokes law with Cunningham slip correction factor
            terminal_vel = (diam**2)*(rhop-rhoa)*grav*beta/(18.*visc) !(4.*radius**2)*(rhop-rhoa)*grav*beta/(18.*visc) ![m s-1]
         ELSE IF (diam<1.2e-3) THEN
            ! Droplets from 40 um to 0.6 mm: linear dependence on particle radius and a correction for reduced pressure
            !   R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
            terminal_vel = 4.e3*diam*SQRT(rhoa_ref/rhoa)!8.e3*radius*sqrt(rhoa_ref/rhoa)
         ELSE
            ! Droplets larger than 0.6 mm: square root dependence on particle radius and a correction for reduced pressure
            !   R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
            ! Note: this is valid up to 2 mm or 9 m/s (at 1000 mbar), where droplets start to break
            terminal_vel = 2.01e2*SQRT( MIN(diam/2.,2.0e-3)*rhoa_ref/rhoa ) !2.01e2*SQRT( MIN(radius,2.0e-3)*rhoa_ref/rhoa)
         END IF
      ELSE IF (flag==4) THEN   ! Ice
         ! Ice crystal terminal fall speed from Ovchinnikov et al. (2014)
         terminal_vel = 12.0*SQRT(diam)
      ELSE IF (flag==5) THEN   ! Snow
         ! The same for snow
         terminal_vel = 12.0*SQRT(diam)
      END IF

    END FUNCTION terminal_vel
    
    !!!!!!!!!!! THEN DIAMETER CALCULATIONS BELOW COULD BE PUT UNDER A COMMON INTERFACE TO FACILITATE DIFFERENT DATATYPES!!!!!
    ! 
    ! Function for calculating dimension (or wet diameter) for any particle type
    ! - Aerosol, cloud and rain are spherical
    ! - Snow and ice can be irregular and their densities can be size-dependent
    !
    ! Correct dimension is needed for irregular particles (e.g. ice and snow) for calculating fall speed (deposition and coagulation)
    ! and capacitance (condensation). Otherwise compact spherical structure can be expected,
    !
    SUBROUTINE calcDiamSALSA(n,ppart,dia)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: n
      TYPE(Section), INTENT(in) :: ppart(n)
      REAL, INTENT(OUT) :: dia(n)
      INTEGER :: i,nspec

      nspec = spec%getNSpec()
      
      dia(:) = 10.e-9
      DO i=1,n
         IF (ppart(i)%numc > ppart(i)%nlim) &
              dia(i)=(SUM(ppart(i)%volc(1:nspec))/ppart(i)%numc/pi6)**(1./3.)
      END DO
      
    END SUBROUTINE calcDiamSALSA
    
    !
    ! Function for calculating effective (wet) radius for any particle type
    ! - Aerosol, cloud and rain are spherical
    ! - Snow and ice can be irregular and their densities can be size-dependent
    !
    ! Correct dimension is needed for irregular particles (e.g. ice and snow) for calculating fall speed (deposition and coagulation)
    ! and capacitance (condensation). Otherwise compact spherical structure can be expected,
    !
    REAL FUNCTION calcDiamLES(n,numc,mass,flag)
      USE mo_submctl, ONLY : pi6
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n ! Number of species
      INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud droplets (2), precip (3), ice (4) and snow(5) particle phases
      REAL, INTENT(IN) :: numc, mass(n)
      
      calcDiamLES=0.
      
      ! Don't calculate if very low number concentration
      IF (numc<1e-15) RETURN
      
      IF (flag==4) THEN   ! Ice
         ! Spherical ice
         calcDiamLES=( SUM(mass(:)/spec%rhoice(:))/numc/pi6)**(1./3.)
      ELSE IF (flag==5) THEN   ! Snow
         ! Spherical snow
         calcDiamLES=( SUM(mass(:)/spec%rhosnow(:))/numc/pi6)**(1./3.)
      ELSE
         ! Radius from total volume of a spherical particle or aqueous droplet
         calcDiamLES=( SUM(mass(:)/spec%rholiq(:))/numc/pi6)**(1./3.)
      ENDIF
      
    END FUNCTION calcDiamLES
    
    !
    ! Function for calculating equilibrium water saturation ratio at droplet surface based on Köhler theory
    !
    REAL FUNCTION calcSweq(part,T)
      TYPE(Section), INTENT(in) :: part ! Any particle
      REAL, INTENT(IN) :: T ! Absolute temperature (K)
      REAL :: dwet
      
      REAL :: znw,zns ! Moles of water and soluble material
      REAL :: zvw, zvs, zvtot ! Volume concentrations of water and soluble material and total dry
      INTEGER :: iwa, ndry ! Index for water, number of "dry" species
      INTEGER :: i
      
      iwa = spec%getIndex("H2O")
      ndry = spec%getNSpec(type="dry")
      
      ! Wet diameter  !! USE THE FUNCTIONS PROVIDED FOR THIS??
      dwet = (SUM(part%volc(:))/part%numc/pi6)**(1./3.)
      
      ! Equilibrium saturation ratio = xw*exp(4*sigma*v_w/(R*T*Dwet))
      
      znw = part%volc(iwa)*spec%rhowa/spec%mwa
      zvw = part%volc(iwa)*spec%rhowa/spec%mwa
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
      
    END FUNCTIOn calcSweq
    
  END MODULE mo_particle_external_properties
  
