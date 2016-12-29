!****************************************************************
!*                                                              *
!*   module MO_SALSA_UPDATE                                 *
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

CONTAINS

  SUBROUTINE distr_update(kproma, kbdim, klev, &
                          paero, pcloud, pprecp, &
                          pice, psnow, level )

    USE mo_submctl
    

    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid points 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical levels

    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), & ! Cloud size distribution and properties
                                      pprecp(kbdim,klev,nprc), & ! Rain drop size distribution and properties
                                      paero(kbdim,klev,fn2b),  & ! Aerosols particle size distribution and properties
                                      pice(kbdim,klev,nice),   & ! Ice particle size distribution and properties
                                      psnow(kbdim,klev,nsnw)     ! Snow flake size distribution and properties
    INTEGER, INTENT(in) :: level                         ! thermodynamical level

    !-- Local variables ----------------------
    INTEGER :: ii, jj, kk, mm
    REAL :: zvpart, znfrac, zvfrac, zVrat, zVilo, zVihi, zVexc, zvdec
    LOGICAL  :: within_bins
    INTEGER :: count

    zvpart = 0.
    zvfrac = 0.

    DO jj = 1,klev
       DO ii = 1,kbdim
          
          ! ------------------------------------------------------------------------
          ! ************* AEROSOLS **************
          ! ------------------------------------------------------------------------

          within_bins = .FALSE.
          !-- Check if the volume of the bin is within bin limits after update
          count = 0
          DO WHILE(.NOT.within_bins)
             within_bins = .TRUE.

             DO kk = fn2b-1,in1a,-1
                mm = 0
                IF (paero(ii,jj,kk)%numc > nlim) THEN

                   zvpart = 0.
                   zvfrac = 0.

                   IF (kk == fn2a) CYCLE 

                   ! Dry volume
                   zvpart = sum(paero(ii,jj,kk)%volc(1:7))/paero(ii,jj,kk)%numc

                   ! Smallest bin cannot decrease
                   IF (paero(ii,jj,kk)%vlolim > zvpart .AND. kk == in1a) CYCLE

                   ! Decreasing bins
                   IF(paero(ii,jj,kk)%vlolim > zvpart) THEN
                      mm = kk - 1
                      IF(kk == in2b) mm = fn1a ! 2b goes to 1a
                      
                      paero(ii,jj,mm)%numc = paero(ii,jj,mm)%numc + paero(ii,jj,kk)%numc 
                      paero(ii,jj,kk)%numc = 0.
                      paero(ii,jj,mm)%volc(:) = paero(ii,jj,mm)%volc(:) + paero(ii,jj,kk)%volc(:)
                      paero(ii,jj,kk)%volc(:) = 0.
                      CYCLE
                   END IF

                   !-- If size bin has not grown, cycle
                   IF(zvpart <= pi6*paero(ii,jj,kk)%dmid**3) CYCLE

                   !-- volume ratio of the size bin
                   zVrat = paero(ii,jj,kk)%vhilim/paero(ii,jj,kk)%vlolim
                
                   !-- particle volume at the low end of the bin
                   zVilo = 2.*zvpart/(1. + zVrat)

                   !-- particle volume at the high end of the bin
                   zVihi = zVrat * zVilo
                   
                   !-- volume in the grown bin which exceeds 
                   !   the bin upper limit
                   zVexc = 0.5*(zVihi + paero(ii,jj,kk)%vhilim)

                   !-- number fraction to be moved to the larger bin
                   znfrac = min(1.,(zVihi-paero(ii,jj,kk)%vhilim) / (zVihi - zVilo))
          
                   !-- volume fraction to be moved to the larger bin
                   !zvfrac = znfrac * zVexc / (1./2.*(zVihi+zVilo))
                   zvfrac = MIN(0.99,znfrac*zVexc/zvpart)

                   IF(zvfrac < 0.) STOP 'Error aerosol 0'
                   !-- update bin
                   mm = kk+1
                   !-- volume
                   paero(ii,jj,mm)%volc(:) = paero(ii,jj,mm)%volc(:) &
                        + znfrac * paero(ii,jj,kk)%numc * zVexc * paero(ii,jj,kk)%volc(:) / &
                        sum(paero(ii,jj,kk)%volc(1:7))

                   paero(ii,jj,kk)%volc(:) = paero(ii,jj,kk)%volc(:) &
                        - znfrac * paero(ii,jj,kk)%numc * zVexc * paero(ii,jj,kk)%volc(:) / &
                        sum(paero(ii,jj,kk)%volc(1:7))

                   !-- number
                   paero(ii,jj,mm)%numc = paero(ii,jj,mm)%numc + znfrac * paero(ii,jj,kk)%numc

                   paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc * (1. - znfrac)

                END IF ! nlim

                IF ( paero(ii,jj,kk)%numc > nlim ) THEN
                   zvpart = sum(paero(ii,jj,kk)%volc(1:7))/paero(ii,jj,kk)%numc  ! Note: dry volume
                   within_bins = (paero(ii,jj,kk)%vlolim<zvpart .AND. zvpart<paero(ii,jj,kk)%vhilim)
                END IF

             END DO ! - kk

             count = count + 1
             IF (count > 100) STOP 'Error: Aerosol bin update not converged'

          END DO ! - within bins

          ! ------------------------------------------------------------------------
          ! ************* CLOUD DROPLETS  **************
          ! ------------------------------------------------------------------------
          
          within_bins = .FALSE.
          count = 0
          DO WHILE (.NOT. within_bins) 
             within_bins = .TRUE.

             DO kk = ncld,ica%cur,-1
                mm = 0

                IF ( pcloud(ii,jj,kk)%numc > nlim .AND. sum(pcloud(ii,jj,kk)%volc(1:7)) > 1.e-30 ) THEN
                
                   ! Don't convert cloud or rain droplets to anything else here.
                   zvpart = sum(pcloud(ii,jj,kk)%volc(1:7))/pcloud(ii,jj,kk)%numc
                
                   !-- volume ratio of the size bin
                   zVrat = pcloud(ii,jj,kk)%vhilim/pcloud(ii,jj,kk)%vlolim
                
                   !-- particle volume at the low end of the bin
                   zVilo = 2.*zvpart/(1. + zVrat)

                   !-- particle volume at the high end of the bin
                   zVihi = zVrat * zVilo

                   !-- Decreasing droplets
                   IF ( zvpart < pi6*pcloud(ii,jj,kk)%vlolim .AND.  &
                        (kk /= ica%cur .AND. kk /= icb%cur)    ) THEN 

                      !-- Volume in the decreased bin which is below the bin lower limit
                      zVexc = 0.5*(zVilo + pcloud(ii,jj,kk)%vlolim)
                      
                      !-- Number fraction to be moved to the smaller bin
                      znfrac = min(1.,(pcloud(ii,jj,kk)%vlolim-zVilo) / (zVihi-zVilo))
                      
                      !-- Index for the smaller bin
                      mm = kk - 1

                   !-- Increasing droplets !! #mergemod
                   ELSE IF ( zvpart > pi6*pcloud(ii,jj,kk)%dmid**3 .AND.  &
                             (kk /= fca%cur .AND. kk /= fcb%cur)     )  THEN  ! Increasing droplets  

                      !-- volume in the grown bin which exceeds the bin upper limit
                      zVexc = 0.5*(zVihi + pcloud(ii,jj,kk)%vhilim)

                      !-- number fraction to be moved to the larger bin
                      znfrac = min(1.,(zVihi-pcloud(ii,jj,kk)%vhilim) / (zVihi-zVilo))
          
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
                   pcloud(ii,jj,mm)%volc(:) = pcloud(ii,jj,mm)%volc(:)     &
                        + zvfrac*pcloud(ii,jj,kk)%volc(:)
                   
                   pcloud(ii,jj,kk)%volc(:) = pcloud(ii,jj,kk)%volc(:)     &
                        - zvfrac*pcloud(ii,jj,kk)%volc(:)
 
                   !-- number
                   pcloud(ii,jj,mm)%numc = pcloud(ii,jj,mm)%numc    &
                        + znfrac*pcloud(ii,jj,kk)%numc
                   
                   pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc    &
                        - znfrac*pcloud(ii,jj,kk)%numc
              
                END IF !nlim
                
                IF ( pcloud(ii,jj,kk)%numc > nlim .AND.  sum(pcloud(ii,jj,kk)%volc(1:7)) > 1.e-30 ) THEN
                   zvpart = sum(pcloud(ii,jj,kk)%volc(1:7))/pcloud(ii,jj,kk)%numc ! Note: dry volume
                   within_bins = (pcloud(ii,jj,kk)%vlolim<zvpart .AND. zvpart<pcloud(ii,jj,kk)%vhilim)
                END IF

             END DO !kk

             count = count + 1
             IF (count > 100) STOP 'Error: Cloud bin update not converged'

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

             ! -- Juha: now the same for cloud bins
             DO kk = nprc,ira,-1
                mm = 0
                IF ( pprecp(ii,jj,kk)%numc > prlim ) THEN
                
                   zvpart = sum(pprecp(ii,jj,kk)%volc(1:8))/pprecp(ii,jj,kk)%numc
                
                   !-- volume ratio of the size bin
                   zVrat = pprecp(ii,jj,kk)%vhilim/pprecp(ii,jj,kk)%vlolim
                
                   !-- particle volume at the low end of the bin
                   zVilo = 2.*zvpart/(1. + zVrat)

                   !-- particle volume at the high end of the bin
                   zVihi = zVrat * zVilo

                   ! Calculate the threshold particle volume for decreasing
                   zvdec = (pi6*pprecp(ii,jj,kk)%dmid**3) - &
                           0.2*((pi6*pprecp(ii,jj,kk)%dmid**3) - pprecp(ii,jj,kk)%vlolim)

                   !-- Decreasing droplets - This is now more critical since we are following the wet diameter!!!
                   IF ( zvpart < zvdec .AND. kk /= ira ) THEN 

                      !-- Volume in the decreased bin which is below the bin lower limit
                      zVexc = 0.5*(zVilo + pprecp(ii,jj,kk)%vlolim)
                      
                      !-- Number fraction to be moved to the smaller bin
                      znfrac = min(1.,(pprecp(ii,jj,kk)%vlolim-zVilo) / (zVihi-zVilo))
                      IF (znfrac < 0.) STOP 'Error, numc precp 0'

                      !-- Index for the smaller bin
                      mm = kk - 1

                   !-- Increasing droplets
                   ELSE IF ( zvpart > pi6*pprecp(ii,jj,kk)%dmid**3 .AND. kk /= fra )  THEN  ! Increasing droplets  

                      !-- volume in the grown bin which exceeds the bin upper limit
                      zVexc = 0.5*(zVihi + pprecp(ii,jj,kk)%vhilim)

                      !-- number fraction to be moved to the larger bin
                      znfrac = min(.99,(zVihi-pprecp(ii,jj,kk)%vhilim) / (zVihi-zVilo))
                      IF (znfrac < 0.) STOP 'Error, numc precp 0, INC'

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
                   pprecp(ii,jj,mm)%volc(:) = pprecp(ii,jj,mm)%volc(:)     &
                        + zvfrac*pprecp(ii,jj,kk)%volc(:)
                   
                   pprecp(ii,jj,kk)%volc(:) = pprecp(ii,jj,kk)%volc(:)     &
                        - zvfrac*pprecp(ii,jj,kk)%volc(:)
 
                   !-- number
                   pprecp(ii,jj,mm)%numc = pprecp(ii,jj,mm)%numc    &
                        + znfrac*pprecp(ii,jj,kk)%numc
                   
                   pprecp(ii,jj,kk)%numc = pprecp(ii,jj,kk)%numc    &
                        - znfrac*pprecp(ii,jj,kk)%numc
              
                END IF !nlim
                
                IF ( pprecp(ii,jj,kk)%numc > prlim ) THEN
                   zvpart = sum(pprecp(ii,jj,kk)%volc(1:8))/pprecp(ii,jj,kk)%numc ! Note: droplet volume
                   within_bins = (pprecp(ii,jj,kk)%vlolim<zvpart .AND. zvpart<pprecp(ii,jj,kk)%vhilim )
                END IF

             END DO !kk

             count = count + 1
             IF (count > 100) STOP 'Error: precipitation bin update not converged'

          END DO !within_bins


        IF(level < 5 ) CYCLE ! skip ice and snow distr. updates if thermodynamical level doesn't include ice microphysics

          ! ------------------------------------------------------------------------
          ! ************* ICE PARTICLES  **************
          ! ------------------------------------------------------------------------

          within_bins = .FALSE.
          count = 0
          DO WHILE (.NOT. within_bins)
             within_bins = .TRUE.

             DO kk = nice,iia%cur,-1
                mm = 0

                IF ( pice(ii,jj,kk)%numc > prlim .AND. sum(pice(ii,jj,kk)%volc(1:7)) > 1.e-30 ) THEN

                   ! Don't convert cloud or rain droplets to anything else here. !!huomhuom
                   zvpart = sum(pice(ii,jj,kk)%volc(1:7))/pice(ii,jj,kk)%numc

                   !-- volume ratio of the size bin
                   zVrat = pice(ii,jj,kk)%vhilim/pice(ii,jj,kk)%vlolim

                   !-- particle volume at the low end of the bin
                   zVilo = 2.*zvpart/(1. + zVrat)

                   !-- particle volume at the high end of the bin
                   zVihi = zVrat * zVilo

                   !-- Decreasing droplets
                   IF ( zvpart < pi6*pice(ii,jj,kk)%vlolim .AND.  &
                        (kk /= iia%cur .AND. kk /= iib%cur)    ) THEN

                      !-- Volume in the decreased bin which is below the bin lower limit
                      zVexc = 0.5*(zVilo + pice(ii,jj,kk)%vlolim)

                      !-- Number fraction to be moved to the smaller bin
                      znfrac = min(1.,(pice(ii,jj,kk)%vlolim-zVilo) / (zVihi-zVilo))

                      !-- Index for the smaller bin
                      mm = kk - 1

                   !-- Increasing droplets
                   ELSE IF ( zvpart > pi6*pice(ii,jj,kk)%dmid**3 .AND.  &
                        (kk /= fia%cur .AND. kk /= fib%cur)    ) THEN !#mod
                    ! Increasing droplets

                      !-- volume in the grown bin which exceeds the bin upper limit
                      zVexc = 0.5*(zVihi + pice(ii,jj,kk)%vhilim)

                      !-- number fraction to be moved to the larger bin
                      znfrac = min(1.,(zVihi-pice(ii,jj,kk)%vhilim) / (zVihi-zVilo))

                      !-- Index for the larger bin
                      mm = kk + 1

                   ELSE  ! Particle size unchanged
                      CYCLE
                   END IF

                   !-- volume fraction to be moved

                   zvfrac = MIN(0.99,znfrac*zVexc/zvpart)
                   IF(zvfrac < 0.) STOP 'Error cloud volc 0'

                   !-- volume
                   pice(ii,jj,mm)%volc(:) = pice(ii,jj,mm)%volc(:)     &
                        + zvfrac*pice(ii,jj,kk)%volc(:)

                   pice(ii,jj,kk)%volc(:) = pice(ii,jj,kk)%volc(:)     &
                        - zvfrac*pice(ii,jj,kk)%volc(:)

                   !-- number
                   pice(ii,jj,mm)%numc = pice(ii,jj,mm)%numc    &
                        + znfrac*pice(ii,jj,kk)%numc

                   pice(ii,jj,kk)%numc = pice(ii,jj,kk)%numc    &
                        - znfrac*pice(ii,jj,kk)%numc

                END IF !nlim

                IF ( pice(ii,jj,kk)%numc > prlim .AND.  sum(pice(ii,jj,kk)%volc(1:7)) > 1.e-30 ) THEN
                   zvpart = sum(pice(ii,jj,kk)%volc(1:7))/pice(ii,jj,kk)%numc ! Note: dry volume
                   within_bins = (pice(ii,jj,kk)%vlolim<zvpart .AND. zvpart<pice(ii,jj,kk)%vhilim)
                END IF

             END DO !kk

             count = count + 1
             IF (count > 100) STOP 'Error: Ice bin update not converged'

          END DO !within_bins


          ! ------------------------------------------------------------------------
          ! ************* SNOW DROPS **************
          ! Everything else the same as with cloud
          ! droplets & aerosols, except that the snow
          ! bins are organized according to the wet
          ! diameter. !!huomhuom onko näin?
          ! ------------------------------------------------------------------------

          within_bins = .FALSE.
          count = 0
          DO WHILE (.NOT. within_bins)
             within_bins = .TRUE.


             DO kk = nsnw,isa,-1
                mm = 0
                IF ( psnow(ii,jj,kk)%numc > prlim ) THEN

                   zvpart = sum(psnow(ii,jj,kk)%volc(1:8))/psnow(ii,jj,kk)%numc

                   !-- volume ratio of the size bin
                   zVrat = psnow(ii,jj,kk)%vhilim/psnow(ii,jj,kk)%vlolim

                   !-- particle volume at the low end of the bin
                   zVilo = 2.*zvpart/(1. + zVrat)

                   !-- particle volume at the high end of the bin
                   zVihi = zVrat * zVilo

                   ! Calculate the threshold particle volume for decreasing
                   zvdec = (pi6*psnow(ii,jj,kk)%dmid**3) - &
                           0.2*((pi6*psnow(ii,jj,kk)%dmid**3) - psnow(ii,jj,kk)%vlolim)

                   !-- Decreasing droplets - This is now more critical since we are following the wet diameter!!!
                   IF ( zvpart < zvdec .AND. kk /= ira ) THEN

                      !-- Volume in the decreased bin which is below the bin lower limit
                      zVexc = 0.5*(zVilo + pprecp(ii,jj,kk)%vlolim)

                      !-- Number fraction to be moved to the smaller bin
                      znfrac = min(1.,(psnow(ii,jj,kk)%vlolim-zVilo) / (zVihi-zVilo))
                      IF (znfrac < 0.) STOP 'Error, numc precp 0'

                      !-- Index for the smaller bin
                      mm = kk - 1

                   !-- Increasing droplets
                   ELSE IF ( zvpart > pi6*psnow(ii,jj,kk)%dmid**3 .AND. kk /= fsa )  THEN  ! Increasing droplets

                      !-- volume in the grown bin which exceeds the bin upper limit
                      zVexc = 0.5*(zVihi + psnow(ii,jj,kk)%vhilim)

                      !-- number fraction to be moved to the larger bin
                      znfrac = min(.99,(zVihi-psnow(ii,jj,kk)%vhilim) / (zVihi-zVilo))
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
                   psnow(ii,jj,mm)%volc(:) = psnow(ii,jj,mm)%volc(:)     &
                        + zvfrac*psnow(ii,jj,kk)%volc(:)

                   psnow(ii,jj,kk)%volc(:) = psnow(ii,jj,kk)%volc(:)     &
                        - zvfrac*psnow(ii,jj,kk)%volc(:)

                   !-- number
                   psnow(ii,jj,mm)%numc = psnow(ii,jj,mm)%numc    &
                        + znfrac*psnow(ii,jj,kk)%numc

                   psnow(ii,jj,kk)%numc = psnow(ii,jj,kk)%numc    &
                        - znfrac*psnow(ii,jj,kk)%numc

                END IF !nlim

                IF ( psnow(ii,jj,kk)%numc > prlim ) THEN
                   zvpart = sum(psnow(ii,jj,kk)%volc(:))/psnow(ii,jj,kk)%numc
                   within_bins = (psnow(ii,jj,kk)%vlolim<zvpart .AND. zvpart<psnow(ii,jj,kk)%vhilim)
                END IF

             END DO !kk

             count = count + 1
             IF (count > 100) STOP 'Error: snow bin update not converged'

          END DO !within_bins

       END DO    ! - ii
    END DO       ! - jj


  END SUBROUTINE distr_update

END MODULE mo_salsa_update
