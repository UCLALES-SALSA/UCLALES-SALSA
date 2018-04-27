MODULE mo_salsa_coagulation
  IMPLICIT NONE
  
  ! ---------------------------------------------------
  ! Contains subroutines and functions used for 
  ! coagulation calculations in mo_salsa_dynamics.f90
  !

  CONTAINS

    SUBROUTINE build_coagulation_kernels(temp, press,                & 
                                         nb1,nlim1,npart1,mpart1,dpart1,   &
                                         nb2,nlim2,npart2,mpart2,dpart2    )
      
      REAL, INTENT(in) :: temp, press
      INTEGER, INTENT(in) :: nb1
      REAL, INTENT(in) :: nlim1, npart1(nb1), mpart1(nb1), dpart1(nb1) ! Number, mass and diameter
      INTEGER, OPTIONAL, INTENT(in) :: nb2
      REAL, OPTIONAL, INTENT(in) :: nlim2, npart2(nb2), mpart2(nb2), dpart2(nb2) ! If not present, do the coagulation kernels with the single particle type      
      
      INTEGER :: mm,nn

      IF ( ALL( [PRESENT(nb2),PRESENT(npart2),PRESENT(mpart2),PRESENT(dpart2)] ) ) THEN

         DO mm = 1, nb1         ! smaller colliding particle
            IF (npart1 < nlim) CYCLE
            DO nn = mm, nb2            ! larger colliding particle
               IF (aero(ii,jj,nn)%numc < nlim) CYCLE
               zcc(mm,nn) = coagc(zdpart(mm),zdpart(nn),zmpart(mm),zmpart(nn),temppi,pressi,1,1,1)
               zcc(nn,mm) = zcc(mm,nn)
            END DO
         END DO
         
      ELSE


   END SUBROUTINE build_coagulation_kernels



   !------------------------------------------------
   !
   ! ***************
   ! Function coagc
   ! ***************
   !
   ! Calculation of coagulation coefficients.
   ! Extended version of the function originally
   ! found in mo_salsa_init. This is now placed
   ! here to avoid cyclic dependencies between
   ! MODULEs upon coupling with UCLALES.
   !
   ! J. Tonttila, FMI, 05/2014
   !
   !-------------------------------------------------
   REAL FUNCTION coagc(diam1,diam2,mass1,mass2,temp,pres,kernel,flag1,flag2)

      USE mo_submctl, ONLY : pi, pi6, boltz, pstand, grav, rd 
      USE mo_particle_external_properties, ONLY : terminal_vel

      IMPLICIT NONE

      !-- Input variables ----------
      REAL, INTENT(IN) :: &
         diam1,  &   ! diameters of colliding particles [m]
         diam2,  &   !
         mass1,  &   ! masses -"- [kg]
         mass2,  &
         temp,   &   ! ambient temperature [K]
         pres        ! ambient pressure [fxm]

      INTEGER, INTENT(in) :: kernel ! select the type of kernel: 1 - aerosol-aerosol coagulation (the original version)
                                        !                            2 - hydrometeor-aerosol or hydrometeor-hydrometeor coagulation
      INTEGER, INTENT(in) :: flag1,flag2 ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)

      !-- Output variables ---------

      !-- Local variables ----------
      REAL ::  &
         visc,     &   ! viscosity of air [kg/(m s)]
         vkin,     &   ! Kinematic viscosity of air [m2 s-1]
         zrhoa,    &   ! Density of air [kg m-3]
         mfp,      &   ! mean free path of air molecules [m]
         mdiam,    &   ! mean diameter of colliding particles [m]
         fmdist,   &   ! distance of flux matching [m]
         eddy_dis, &   ! Eddy dissipation time
         zecoll,   &   ! Collision efficiency for graviational collection
         zev,      &   !
         zea,      &
         zbrown,   &   ! Components for coagulation kernel; Brownian
         zbrconv,  &   !                                    Convective diffusion enhancement
         zgrav,    &   !                                    Gravitational collection
         ztshear,  &   ! turbulent shear
         zturbinert    ! turbulent inertia

      REAL, DIMENSION (2) :: &
         diam,   &   ! diameters of particles [m]
         mpart,  &   ! masses of particles [kg]
         knud,   &   ! particle knudsen number [1]
         beta,   &   ! Cunningham correction factor [1]
         zrhop,  &   ! Particle density [kg m-3]
         dfpart, &   ! particle diffusion coefficient [m2/s]
         mtvel,  &   ! particle mean thermal velocity [m/s]
         termv,  &   ! Particle terminal velocity
         omega,  &   !
         tva,    &   ! temporary variable [m]
         flux        ! flux in continuum and free molec. regime [m/s]

      REAL ::  &
         schm(2), &   ! Schmidt nubmer
         reyn(2), &    ! Reynolds number
         stok             ! Stokes number
      INTEGER :: lrg,sml

      zbrown = 0.
      zbrconv = 0.
      zgrav = 0.
      zev = 0.
      coagc = 0.

      !-------------------------------------------------------------------------------

      !-- 0) Initializing particle and ambient air variables --------------------
      diam  = (/ diam1, diam2 /)       ! particle diameters [m]
      mpart = (/ mass1, mass2 /)       ! particle masses [kg]

      visc = (7.44523e-3*SQRT(temp**3))/(5093.*(temp+110.4)) ! viscosity of air [kg/(m s)]

      mfp = (1.656e-10*temp+1.828e-8)*pstand/pres ! mean free path of air [m]

      !-- 2) Slip correction factor for small particles -------------------------

      knud = 2.*mfp/diam                                    ! Knudsen number
      beta = 1.+knud*(1.142+0.558*exp(-0.999/knud))! Cunningham correction factor
      ! (Allen and Raabe, Aerosol Sci. Tech. 4, 269)

      !-- 3) Particle properties ------------------------------------------------

      dfpart = beta*boltz*temp/(3.*pi*visc*diam)  ! diffusion coefficient [m2/s]
      mtvel  = sqrt((8.*boltz*temp)/(pi*mpart))    ! mean thermal velocity [m/s]
      omega  = 8.*dfpart/(pi*mtvel)

      mdiam = 0.5*(diam(1)+diam(2))               ! mean diameter [m]

      !-- 4) Calculation of fluxes and flux matching ----------------------------

      flux(1) = 4.*pi*mdiam*( dfpart(1)+dfpart(2) )    ! flux in continuum regime [m3/s]
      flux(2) = pi*sqrt((mtvel(1)**2)+(mtvel(2)**2))*(mdiam**2) !  -"- in free molec. regime [m3/s]

      tva(1) = ((mdiam+omega(1))**3 - &              ! temporary variable [m]
                (mdiam**2+omega(1)**2)* &
                sqrt((mdiam**2+omega(1)**2)))/ &
               (3.*mdiam*omega(1)) - mdiam

      tva(2) = ((mdiam+omega(2))**3 - &              ! temporary variable [m]
                (mdiam**2+omega(2)**2)* &
                sqrt((mdiam**2+omega(2)**2)))/ &
               (3.*mdiam*omega(2)) - mdiam

      fmdist = sqrt(tva(1)**2+tva(2)**2)             ! flux matching distance [m]

      SELECT CASE(kernel)
         CASE(1)

            ! Aerosol-Aerosol coagulation - like the f version
            !-- 5) Coagulation coefficient [m3/s] -------------------------------------
            coagc = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2))

         CASE(2)

            ! Which particle is larger?
            sml = 1; lrg = 2
            IF (diam(1) >= diam(2)) THEN
               lrg = 1; sml = 2
            END IF

            zrhoa = pres/(rd*temp)   ! Density of air
            zrhop = mpart/(pi6*diam**3)             ! Density of particles
            vkin  = visc/zrhoa   ! Kinematic viscosity of air [m2 s-1]

            termv(1) = terminal_vel(diam(1)/2.,zrhop(1),zrhoa,visc,beta(1),flag1)
            termv(2) = terminal_vel(diam(2)/2.,zrhop(2),zrhoa,visc,beta(2),flag2)

            ! Reynolds number
            reyn = diam*termv/vkin
            ! Schmidt number for the smaller particle
            schm = vkin/dfpart
            ! Stokes number
            stok = 2.*termv(sml)*ABS(termv(1) - termv(2))/( diam(lrg)*grav )

            !Brownian component
            zbrown = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2))

            ! Convective enhancement
            IF (reyn(lrg) <= 1.) THEN
               zbrconv = 0.45*zbrown*( (reyn(lrg)*schm(sml))**(1./3.) )
            ELSE IF (reyn(lrg) > 1.) THEN
               zbrconv = 0.45*zbrown*SQRT(reyn(lrg))*( schm(sml)**(1./3.) )
            END IF

            ! Turbulent Shear
            eddy_dis = 10.e-4 ! Values suggested by Sami - could be taken from the LES model?
            ztshear = SQRT(8.*pi*eddy_dis/(15.*vkin))*(0.5*(diam(1)+diam(2)))**3
            zturbinert = pi*eddy_dis**(0.75) /(grav* SQRT(SQRT( vkin )))  &
                         *(0.5*(diam(1)+diam(2)))**2* ABS(termv(1)-termv(2))

            ! gravitational collection
            zea = stok**2/( stok + 0.5 )**2
            IF (stok > 1.214) THEN
               zev = 0.75*LOG(2.*stok)/(stok - 1.214)
               zev = (1. + zev)**(-2)
            ELSE IF (stok <= 1.214) THEN
               zev = 0.
            END IF

            zecoll = (60.*zev + zea*reyn(lrg))/(60. + reyn(lrg))
            zgrav = zecoll * pi * mdiam**2
            zgrav = zgrav * ABS(termv(1)-termv(2))

            ! Total coagulation kernel
            coagc = zbrown  + zbrconv + SQRT(zgrav**2+ ztshear**2+ zturbinert**2)

      END SELECT

   END FUNCTION coagc




END MODULE mo_salsa_coagulation 
