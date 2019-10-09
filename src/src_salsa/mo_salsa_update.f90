!****************************************************************
!*                                                              *
!*   MODULE MO_SALSA_UPDATE                                 *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate aerosol dynamics                              *
!*                                                              *
!****************************************************************
!
! -- Added update for cloud bins (05/2014 J Tonttila, FMI)
!
!****************************************************************

MODULE mo_salsa_update
   USE classSection
   IMPLICIT NONE

CONTAINS

   SUBROUTINE distr_update(kproma, kbdim, klev, allSALSA, level )
      USE mo_submctl
    

      IMPLICIT NONE

      !-- Input and output variables ----------
      INTEGER, INTENT(IN) ::      &
         kproma,                    & ! number of horiz. grid points 
         kbdim,                     & ! dimension for arrays
         klev                         ! number of vertical levels
      
      TYPE(Section), INTENT(inout) :: allSALSA(kbdim,klev,ntotal)

      INTEGER, INTENT(in) :: level                         ! thermodynamical level

      !-- Local variables ----------------------
      INTEGER :: ii, jj, kk, mm
      REAL    :: zvpart, znfrac, zvfrac, zVrat, zVilo, zVihi, zVexc, zvdec
      LOGICAL :: within_bins
      INTEGER :: count
      INTEGER :: zndry,znwet

      zndry = spec%getNSpec(type="dry")
      znwet = spec%getNSpec(type="wet")

      zvpart = 0.
      zvfrac = 0.

      DO jj = 1, klev
         DO ii = 1, kbdim
          
            ! ------------------------------------------------------------------------
            ! ************* AEROSOLS **************
            ! ------------------------------------------------------------------------

            within_bins = .FALSE.
            !-- Check if the volume of the bin is within bin limits after update
            count = 0
            DO WHILE(.NOT. within_bins)
               within_bins = .TRUE.

               DO kk = fn2b-1, in1a, -1
                  mm = 0
                  IF (aero(ii,jj,kk)%numc > nlim) THEN

                     zvpart = 0.
                     zvfrac = 0.

                     IF (kk == fn2a) CYCLE

                     ! Dry volume
                     zvpart = sum(aero(ii,jj,kk)%volc(1:zndry))/aero(ii,jj,kk)%numc

                     ! Smallest bin cannot decrease
                     IF (aero(ii,jj,kk)%vlolim > zvpart .AND. kk == in1a) CYCLE

                     ! Decreasing bins
                     IF(aero(ii,jj,kk)%vlolim > zvpart) THEN
                        mm = kk - 1
                        IF(kk == in2b) mm = fn1a ! 2b goes to 1a
                      
                        aero(ii,jj,mm)%numc = aero(ii,jj,mm)%numc + aero(ii,jj,kk)%numc
                        aero(ii,jj,kk)%numc = 0.
                        aero(ii,jj,mm)%volc(1:znwet) = aero(ii,jj,mm)%volc(1:znwet) + aero(ii,jj,kk)%volc(1:znwet)
                        aero(ii,jj,kk)%volc(1:znwet) = 0.
                        CYCLE
                     END IF

                     !-- If size bin has not grown, CYCLE
                     IF (zvpart <= pi6*aero(ii,jj,kk)%dmid**3) CYCLE

                     !-- volume ratio of the size bin
                     zVrat = aero(ii,jj,kk)%vhilim/aero(ii,jj,kk)%vlolim
                
                     !-- particle volume at the low end of the bin
                     zVilo = 2.*zvpart/(1. + zVrat)

                     !-- particle volume at the high end of the bin
                     zVihi = zVrat * zVilo
                   
                     !-- volume in the grown bin which exceeds
                     !   the bin upper limit
                     zVexc = 0.5*(zVihi + aero(ii,jj,kk)%vhilim)

                     !-- number fraction to be moved to the larger bin
                     znfrac = min(1., (zVihi-aero(ii,jj,kk)%vhilim) / (zVihi - zVilo))
          
                     !-- volume fraction to be moved to the larger bin
                     !zvfrac = znfrac * zVexc / (1./2.*(zVihi+zVilo))
                     zvfrac = MIN(0.99, znfrac*zVexc/zvpart)

                     IF(zvfrac < 0.) STOP 'Error aerosol 0'
                     !-- update bin
                     mm = kk+1
                     !-- volume
                     aero(ii,jj,mm)%volc(1:znwet) = aero(ii,jj,mm)%volc(1:znwet) &
                                               + znfrac * aero(ii,jj,kk)%numc * zVexc * aero(ii,jj,kk)%volc(1:znwet) / &
                                               sum(aero(ii,jj,kk)%volc(1:zndry))

                     aero(ii,jj,kk)%volc(1:znwet) = aero(ii,jj,kk)%volc(1:znwet) &
                                               - znfrac * aero(ii,jj,kk)%numc * zVexc * aero(ii,jj,kk)%volc(1:znwet) / &
                                               sum(aero(ii,jj,kk)%volc(1:zndry))

                     !-- number
                     aero(ii,jj,mm)%numc = aero(ii,jj,mm)%numc + znfrac * aero(ii,jj,kk)%numc

                     aero(ii,jj,kk)%numc = aero(ii,jj,kk)%numc * (1. - znfrac)

                  END IF ! nlim

                  IF ( aero(ii,jj,kk)%numc > nlim ) THEN
                     zvpart = sum(aero(ii,jj,kk)%volc(1:zndry))/aero(ii,jj,kk)%numc  ! Note: dry volume
                     within_bins = (aero(ii,jj,kk)%vlolim < zvpart .AND. zvpart < aero(ii,jj,kk)%vhilim)
                  END IF

               END DO ! - kk

               count = count + 1
               IF (count > 100)  STOP  "Aerosol bin update not converged"


            END DO ! - within bins

            ! ------------------------------------------------------------------------
            ! ************* CLOUD DROPLETS  **************
            ! ------------------------------------------------------------------------
          
            within_bins = .FALSE.
            count = 0
            DO WHILE (.NOT. within_bins)
               within_bins = .TRUE.

               DO kk = ncld, ica%cur, -1
                  mm = 0

                  IF ( cloud(ii,jj,kk)%numc > nlim .AND. sum(cloud(ii,jj,kk)%volc(1:zndry)) > 1.e-30 ) THEN
                
                     ! Don't convert cloud or rain droplets to anything else here.
                     zvpart = sum(cloud(ii,jj,kk)%volc(1:zndry))/cloud(ii,jj,kk)%numc
                
                     !-- volume ratio of the size bin
                     zVrat = cloud(ii,jj,kk)%vhilim/cloud(ii,jj,kk)%vlolim
                
                     !-- particle volume at the low end of the bin
                     zVilo = 2.*zvpart/(1. + zVrat)

                     !-- particle volume at the high end of the bin
                     zVihi = zVrat * zVilo

                     !-- Decreasing droplets
                     IF ( zvpart < pi6*cloud(ii,jj,kk)%vlolim .AND.  &
                         (kk /= ica%cur .AND. kk /= icb%cur)    ) THEN

                        !-- Volume in the decreased bin which is below the bin lower limit
                        zVexc = 0.5*(zVilo + cloud(ii,jj,kk)%vlolim)
                      
                        !-- Number fraction to be moved to the smaller bin
                        znfrac = min(1.,(cloud(ii,jj,kk)%vlolim-zVilo) / (zVihi-zVilo))
                      
                        !-- Index for the smaller bin
                        mm = kk - 1

                     !-- Increasing droplets !! #mergemod
                     ELSE IF ( zvpart > pi6*cloud(ii,jj,kk)%dmid**3 .AND.  &
                              (kk /= fca%cur .AND. kk /= fcb%cur)     )  THEN  ! Increasing droplets

                        !-- volume in the grown bin which exceeds the bin upper limit
                        zVexc = 0.5*(zVihi + cloud(ii,jj,kk)%vhilim)

                        !-- number fraction to be moved to the larger bin
                        znfrac = min(1.,(zVihi-cloud(ii,jj,kk)%vhilim) / (zVihi-zVilo))
          
                        !-- Index for the larger bin
                        mm = kk + 1

                     ELSE  ! Particle size unchanged
                        CYCLE
                     END IF

                     !-- volume fraction to be moved
                     !zvfrac = znfrac * zVexc / (0.5*(zVihi+zVilo))
                     !zvfrac = min(zvfrac,1.)
                     zvfrac = MIN(0.99,znfrac*zVexc/zvpart)

                     IF(zvfrac < 0.) STOP 'Error cloud 0'

                     !-- volume
                     cloud(ii,jj,mm)%volc(1:znwet) = cloud(ii,jj,mm)%volc(1:znwet)     &
                                                + zvfrac*cloud(ii,jj,kk)%volc(1:znwet)
                   
                     cloud(ii,jj,kk)%volc(1:znwet) = cloud(ii,jj,kk)%volc(1:znwet)     &
                                                - zvfrac*cloud(ii,jj,kk)%volc(1:znwet)
 
                     !-- number
                     cloud(ii,jj,mm)%numc = cloud(ii,jj,mm)%numc    &
                                             + znfrac*cloud(ii,jj,kk)%numc
                   
                     cloud(ii,jj,kk)%numc = cloud(ii,jj,kk)%numc    &
                                             - znfrac*cloud(ii,jj,kk)%numc
              
                  END IF !nlim
                
                  IF ( cloud(ii,jj,kk)%numc > nlim .AND.  sum(cloud(ii,jj,kk)%volc(1:zndry)) > 1.e-30 ) THEN
                     zvpart = sum(cloud(ii,jj,kk)%volc(1:zndry))/cloud(ii,jj,kk)%numc ! Note: dry volume
                     within_bins = (cloud(ii,jj,kk)%vlolim < zvpart .AND. zvpart < cloud(ii,jj,kk)%vhilim)
                  END IF

               END DO !kk

               count = count + 1

               IF (count > 100)  STOP "Cloud bin update not converged"


            END DO !within_bins


            ! ------------------------------------------------------------------------
            ! ************* RAIN DROPS **************
            ! Everything else the same as with cloud
            ! droplets & aerosols, except that the rain
            ! bins are organized according to the wet
            ! diameter.
            ! ------------------------------------------------------------------------
          
            within_bins = .FALSE.
            count = 0
            DO WHILE (.NOT. within_bins)
               within_bins = .TRUE.

               ! -- Juha: now the same for precp bins
               DO kk = nprc, ira, -1
                  mm = 0
                  IF ( precp(ii,jj,kk)%numc > prlim ) THEN
                
                     zvpart = sum(precp(ii,jj,kk)%volc(1:znwet))/precp(ii,jj,kk)%numc
                
                     !-- volume ratio of the size bin
                     zVrat = precp(ii,jj,kk)%vhilim/precp(ii,jj,kk)%vlolim
                
                     !-- particle volume at the low end of the bin
                     zVilo = 2.*zvpart/(1. + zVrat)

                     !-- particle volume at the high end of the bin
                     zVihi = zVrat * zVilo

                     ! Calculate the threshold particle volume for decreasing
                     zvdec = (pi6*precp(ii,jj,kk)%dmid**3) - &
                             0.2*((pi6*precp(ii,jj,kk)%dmid**3) - precp(ii,jj,kk)%vlolim)

                     !-- Decreasing droplets - This is now more critical since we are following the wet diameter!!!
                     IF ( zvpart < zvdec .AND. kk /= ira ) THEN

                        !-- Volume in the decreased bin which is below the bin lower limit
                        zVexc = 0.5*(zVilo + precp(ii,jj,kk)%vlolim)
                     
                        !-- Number fraction to be moved to the smaller bin
                        znfrac = min(1.,(precp(ii,jj,kk)%vlolim-zVilo) / (zVihi-zVilo))
                        IF (znfrac < 0.) STOP 'Error, numc precp 0'

                        !-- Index for the smaller bin
                        mm = kk - 1

                     !-- Increasing droplets
                     ELSE IF ( zvpart > pi6*precp(ii,jj,kk)%dmid**3 .AND. kk /= fra )  THEN  ! Increasing droplets

                        !-- volume in the grown bin which exceeds the bin upper limit
                        zVexc = 0.5*(zVihi + precp(ii,jj,kk)%vhilim)

                        !-- number fraction to be moved to the larger bin
                        znfrac = min(.99,(zVihi-precp(ii,jj,kk)%vhilim) / (zVihi-zVilo))

                        !-- Index for the larger bin
                        mm = kk + 1

                     ELSE  ! Particle size unchanged
                        CYCLE
                     END IF

                     !-- volume fraction to be moved
                     !zvfrac = min(0.99,znfrac * zVexc / (0.5*(zVihi+zVilo)))
                     zvfrac = MIN(0.99,znfrac*zVexc/zvpart)
                     IF(zvfrac < 0.) STOP 'Error, volc precp 0'

                     !-- volume
                     precp(ii,jj,mm)%volc(1:znwet) = precp(ii,jj,mm)%volc(1:znwet)     &
                                                + zvfrac*precp(ii,jj,kk)%volc(1:znwet)
                   
                     precp(ii,jj,kk)%volc(1:znwet) = precp(ii,jj,kk)%volc(1:znwet)     &
                                                - zvfrac*precp(ii,jj,kk)%volc(1:znwet)
 
                     !-- number
                     precp(ii,jj,mm)%numc = precp(ii,jj,mm)%numc    &
                                             + znfrac*precp(ii,jj,kk)%numc
                   
                     precp(ii,jj,kk)%numc = precp(ii,jj,kk)%numc    &
                                             - znfrac*precp(ii,jj,kk)%numc
              
                  END IF !nlim
                
                  IF ( precp(ii,jj,kk)%numc > prlim ) THEN
                     zvpart = sum(precp(ii,jj,kk)%volc(1:znwet))/precp(ii,jj,kk)%numc ! Note: droplet volume
                     within_bins = (precp(ii,jj,kk)%vlolim < zvpart .AND. zvpart < precp(ii,jj,kk)%vhilim )
                  END IF

               END DO !kk

               count = count + 1

               IF (count > 100)  STOP "Precip bin update not converged"

            END DO !within_bins

            IF(level < 5 ) CYCLE ! skip ice and snow distr. updates if thermodynamical level doesn't include ice microphysics

            ! ------------------------------------------------------------------------
            ! ************* ICE PARTICLES  **************
            ! ------------------------------------------------------------------------

            within_bins = .FALSE.
            count = 0
            DO WHILE (.NOT. within_bins)
               within_bins = .TRUE.

               DO kk = nice, iia%cur, -1
                  mm = 0

                  IF ( ice(ii,jj,kk)%numc > prlim .AND. sum(ice(ii,jj,kk)%volc(1:zndry)) > 1.e-30 ) THEN

                     ! Don't convert ice to anything else here. 
                     zvpart = sum(ice(ii,jj,kk)%volc(1:zndry))/ice(ii,jj,kk)%numc

                     !-- volume ratio of the size bin
                     zVrat = ice(ii,jj,kk)%vhilim/ice(ii,jj,kk)%vlolim

                     !-- particle volume at the low end of the bin
                     zVilo = 2.*zvpart/(1. + zVrat)

                     !-- particle volume at the high end of the bin
                     zVihi = zVrat * zVilo

                     !-- Decreasing size
                     IF ( zvpart < pi6*ice(ii,jj,kk)%vlolim .AND.  &
                         (kk /= iia%cur .AND. kk /= iib%cur)    ) THEN

                        !-- Volume in the decreased bin which is below the bin lower limit
                        zVexc = 0.5*(zVilo + ice(ii,jj,kk)%vlolim)

                        !-- Number fraction to be moved to the smaller bin
                        znfrac = min(1.,(ice(ii,jj,kk)%vlolim-zVilo) / (zVihi-zVilo))

                        !-- Index for the smaller bin
                        mm = kk - 1

                     !-- Increasing size
                     ELSE IF ( zvpart > pi6*ice(ii,jj,kk)%dmid**3 .AND.  &
                              (kk /= fia%cur .AND. kk /= fib%cur)    ) THEN !#mod

                        !-- volume in the grown bin which exceeds the bin upper limit
                        zVexc = 0.5*(zVihi + ice(ii,jj,kk)%vhilim)

                        !-- number fraction to be moved to the larger bin
                        znfrac = min(1.,(zVihi-ice(ii,jj,kk)%vhilim) / (zVihi-zVilo))

                        !-- Index for the larger bin
                        mm = kk + 1

                     ELSE  ! Particle size unchanged
                        CYCLE
                     END IF

                     !-- volume fraction to be moved

                     zvfrac = MIN(0.99,znfrac*zVexc/zvpart)
                     IF(zvfrac < 0.) STOP 'Error ice volc 0'

                     !-- volume
                     ice(ii,jj,mm)%volc(1:znwet) = ice(ii,jj,mm)%volc(1:znwet)     &
                                              + zvfrac*ice(ii,jj,kk)%volc(1:znwet)

                     ice(ii,jj,kk)%volc(1:znwet) = ice(ii,jj,kk)%volc(1:znwet)     &
                                              - zvfrac*ice(ii,jj,kk)%volc(1:znwet)

                     !-- number
                     ice(ii,jj,mm)%numc = ice(ii,jj,mm)%numc    &
                                           + znfrac*ice(ii,jj,kk)%numc

                     ice(ii,jj,kk)%numc = ice(ii,jj,kk)%numc    &
                                           - znfrac*ice(ii,jj,kk)%numc

                  END IF !nlim

                  IF ( ice(ii,jj,kk)%numc > prlim .AND. sum(ice(ii,jj,kk)%volc(1:zndry)) > 1.e-30 ) THEN
                     zvpart = sum(ice(ii,jj,kk)%volc(1:zndry))/ice(ii,jj,kk)%numc ! Note: dry volume
                     within_bins = (ice(ii,jj,kk)%vlolim < zvpart .AND. zvpart < ice(ii,jj,kk)%vhilim)
                  END IF

               END DO !kk

               count = count + 1
               IF (count > 100)  STOP "Ice bin update not converged"

            END DO !within_bins


            ! ------------------------------------------------------------------------
            ! ************* SNOW DROPS **************
            ! Everything else the same as with cloud
            ! droplets & aerosols, except that the snow
            ! bins are organized according to the wet
            ! diameter.
            ! ------------------------------------------------------------------------

            within_bins = .FALSE.
            count = 0
            DO WHILE (.NOT. within_bins)
               within_bins = .TRUE.


               DO kk = nsnw, isa, -1
                  mm = 0
                  IF ( snow(ii,jj,kk)%numc > prlim ) THEN

                     zvpart = sum(snow(ii,jj,kk)%volc(1:znwet))/snow(ii,jj,kk)%numc

                     !-- volume ratio of the size bin
                     zVrat = snow(ii,jj,kk)%vhilim/snow(ii,jj,kk)%vlolim

                     !-- particle volume at the low end of the bin
                     zVilo = 2.*zvpart/(1. + zVrat)

                     !-- particle volume at the high end of the bin
                     zVihi = zVrat * zVilo

                     ! Calculate the threshold particle volume for decreasing
                     zvdec = (pi6*snow(ii,jj,kk)%dmid**3) - &
                             0.2*((pi6*snow(ii,jj,kk)%dmid**3) - snow(ii,jj,kk)%vlolim)

                     !-- Decreasing droplets - This is now more critical since we are following the wet diameter!!!
                     IF ( zvpart < zvdec .AND. kk /= ira ) THEN

                        !-- Volume in the decreased bin which is below the bin lower limit
                        zVexc = 0.5*(zVilo + snow(ii,jj,kk)%vlolim)

                        !-- Number fraction to be moved to the smaller bin
                        znfrac = min(1.,(snow(ii,jj,kk)%vlolim-zVilo) / (zVihi-zVilo))
                        IF (znfrac < 0.) STOP 'Error, numc snow 0'

                        !-- Index for the smaller bin
                        mm = kk - 1

                     !-- Increasing droplets
                     ELSE IF ( zvpart > pi6*snow(ii,jj,kk)%dmid**3 .AND. kk /= fsa )  THEN  ! Increasing droplets

                        !-- volume in the grown bin which exceeds the bin upper limit
                        zVexc = 0.5*(zVihi + snow(ii,jj,kk)%vhilim)

                        !-- number fraction to be moved to the larger bin
                        znfrac = min(.99,(zVihi-snow(ii,jj,kk)%vhilim) / (zVihi-zVilo))
                        IF (znfrac < 0.) STOP 'Error, snow numc 0'

                        !-- Index for the larger bin
                        mm = kk + 1

                     ELSE  ! Particle size unchanged
                        CYCLE
                     END IF

                     !-- volume fraction to be moved
                     zvfrac = MIN(0.99,znfrac*zVexc/zvpart)
                     IF(zvfrac < 0.) STOP 'Error: snow volc 0'

                     !-- volume
                     snow(ii,jj,mm)%volc(1:znwet) = snow(ii,jj,mm)%volc(1:znwet)     &
                                               + zvfrac*snow(ii,jj,kk)%volc(1:znwet)

                     snow(ii,jj,kk)%volc(1:znwet) = snow(ii,jj,kk)%volc(1:znwet)     &
                                               - zvfrac*snow(ii,jj,kk)%volc(1:znwet)

                     !-- number
                     snow(ii,jj,mm)%numc = snow(ii,jj,mm)%numc    &
                                            + znfrac*snow(ii,jj,kk)%numc

                     snow(ii,jj,kk)%numc = snow(ii,jj,kk)%numc    &
                                            - znfrac*snow(ii,jj,kk)%numc

                  END IF !nlim

                  IF ( snow(ii,jj,kk)%numc > prlim ) THEN
                     zvpart = sum(snow(ii,jj,kk)%volc(1:znwet))/snow(ii,jj,kk)%numc
                     within_bins = (snow(ii,jj,kk)%vlolim < zvpart .AND. zvpart < snow(ii,jj,kk)%vhilim)
                  END IF

               END DO !kk

               count = count + 1

               IF (count > 100)  STOP  "Snow bin update not converged"


            END DO !within_bins

         END DO    ! - ii
      END DO       ! - jj


   END SUBROUTINE distr_update

END MODULE mo_salsa_update
