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

  ! ------------------------------------------------------------------------
  ! Size distribution is updated when particle or droplet volume exceeds certain limits.
  ! The limits are relative to bin midpoint volume and bin width:
  !     Vmid-f_low*(Vmid-Vlo) < V < Vmid+f_high*(Vhi-Vmid)
  ! Equally:
  !     -f_low < (V-Vmid)/((Vhi-Vlo)/2) < f_high
  !     -f_low < (2*V-Vhi-Vlo)/(Vhi-Vlo) < f_high
  ! Update moves particles that exceed the high or low limit to next of previous bin, respectively,
  ! so this is an iterative approach.
  !
  SUBROUTINE distr_update(kbdim, klev, paero, pcloud, pprecp, pice, psnow, level)
    USE mo_submctl, ONLY : t_section, nspec, nbins, ncld, nprc, nice, nsnw, &
        fn2a, fnp2a, nlim, prlim
    IMPLICIT NONE
    !-- Parameters --------------------------
    ! Asymmetric limits (default) - distributions adjuted especially after cloud activation
    !REAL :: dry_lo=1.0, dry_hi=0.0, wet_lo=0.2, wet_hi=0.0
    ! Strict limits - number concentrations start to drift due to numerical diffusion, not good
    !REAL :: dry_lo=0.1, dry_hi=0.1, wet_lo=0.1, wet_hi=0.1
    ! Symmetric and relaxed limits - less numerical diffusion, looks good
    REAL :: dry_lo=1.0, dry_hi=1.0, wet_lo=0.2, wet_hi=0.2

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) :: kbdim, klev
    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,nbins), &
        pcloud(kbdim,klev,ncld), pprecp(kbdim,klev,nprc), &
        pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    INTEGER, INTENT(in) :: level

    !-- Local variables ----------------------
    INTEGER :: ii, jj, nn

    nn = nspec + 1 ! Aerosol + water
    DO jj = 1,klev
        DO ii = 1,kbdim
            ! Aerosol (dry size)
            CALL bin_adjust(paero(ii,jj,:),nbins,fn2a,2,nn,nlim,dry_lo,dry_hi)

            ! Cloud droplets (dry size)
            CALL bin_adjust(pcloud(ii,jj,:),ncld,fnp2a,2,nn,nlim,dry_lo,dry_hi)

            ! Rain drops (wet size)
            CALL bin_adjust(pprecp(ii,jj,:),nprc,nprc,1,nn,prlim,wet_lo,wet_hi)

            IF(level > 4) THEN
                ! Ice (dry size)
                CALL bin_adjust(pice(ii,jj,:),nice,fnp2a,2,nn,prlim,dry_lo,dry_hi)

                ! Snow (wet size)
                CALL bin_adjust(psnow(ii,jj,:),nsnw,nsnw,1,nn,prlim,wet_lo,wet_hi)
            ENDIF
        END DO
    END DO

  END SUBROUTINE distr_update
  !
  SUBROUTINE bin_adjust(distr,n,na,n1,nn,lim,frac_lo,frac_hi)
    USE mo_submctl, ONLY : t_section
    IMPLICIT NONE
    !-- Input and output variables ----------
    TYPE(t_section), INTENT(inout) :: distr(n)
    INTEGER, INTENT(in) :: n, na  ! n is the total number of bins; na is the number of a-bins (na=n for rain and snow)
    INTEGER, INTENT(in) :: n1, nn ! n1=1 for wet size and 2 for dry size; nn is the total number of species
    REAL, INTENT(in) :: lim       ! Concentration limit
    REAL, INTENT(in) :: frac_lo, frac_hi ! Low and high limits for distribution update
    !-- Local variables ----------------------
    INTEGER :: kk, mm, iter, b1
    REAL :: zvpart, znfrac, zvfrac, zVrat, zVilo, zVihi, zVexc, dist
    LOGICAL  :: within_bins

    ! The first b-bin
    b1=na+1

    !-- Check if the volume of the bin is within bin limits after update
    within_bins = .FALSE.
    iter = 0
    DO WHILE(.NOT.within_bins)
        within_bins = .TRUE.
        DO kk = n,1,-1
            IF (distr(kk)%numc > lim) THEN
                !-- particle volume (dry or wet)
                zvpart = sum(distr(kk)%volc(n1:nn))/distr(kk)%numc

                !-- relative distance (-1...+1) from the bin midpoint volume
                dist = (2.*zvpart-distr(kk)%vhilim-distr(kk)%vlolim)/(distr(kk)%vhilim-distr(kk)%vlolim)

                !-- volume ratio of the size bin
                zVrat = distr(kk)%vhilim/distr(kk)%vlolim
                !-- particle volume at the low end of the bin
                zVilo = 2.*zvpart/(1. + zVrat)
                !-- particle volume at the high end of the bin
                zVihi = zVrat * zVilo

                IF(-1.*dist>frac_lo .AND. kk/=1 .AND. kk/=b1) THEN
                    !-- volume in the bin which exceeds the bin lower limit
                    zVexc = 0.5*(zVilo + distr(kk)%vlolim)
                    !-- number fraction to be moved to the smaller bin
                    znfrac = min(1.,(distr(kk)%vlolim-zVilo) / (zVihi - zVilo))
                    !-- index for the smaller bin
                    mm = kk - 1
                ELSEIF (dist>frac_hi .AND. kk/=na .AND. kk/=n) THEN
                    !-- volume in the bin which exceeds the bin upper limit
                    zVexc = 0.5*(zVihi + distr(kk)%vhilim)
                    !-- number fraction to be moved to the larger bin
                    znfrac = min(1.,(zVihi-distr(kk)%vhilim) / (zVihi - zVilo))
                    !-- index for the larger bin
                    mm = kk+1
                ELSE
                    CYCLE
                END IF

                !-- ignore irrelevant changes
                IF (ABS(znfrac * distr(kk)%numc) < lim) CYCLE

                !-- volume fraction to be moved
                zvfrac = znfrac*zVexc/zvpart

                IF (znfrac>1. .OR. zvfrac>1.) THEN
                    ! Move complete bin
                    znfrac=1.
                    zvfrac=1.
                ELSEIF (znfrac<0. .OR. zvfrac<0) THEN
                    WRITE(*,*) znfrac, zvfrac
                    WRITE(*,*) zvpart,distr(kk)%vlolim,0.5*(distr(kk)%vlolim+distr(kk)%vhilim),distr(kk)%vhilim
                    WRITE(*,*) dist, frac_lo, frac_hi
                    STOP 'Error: Bin update with invalid numbers'
                ENDIF

                !-- update bin volume
                distr(mm)%volc(:) = distr(mm)%volc(:) + zvfrac * distr(kk)%volc(:)
                distr(kk)%volc(:) = distr(kk)%volc(:) - zvfrac * distr(kk)%volc(:)
                !-- and number
                distr(mm)%numc = distr(mm)%numc + znfrac * distr(kk)%numc
                distr(kk)%numc = distr(kk)%numc - znfrac * distr(kk)%numc

            END IF ! lim

            IF ( within_bins .AND. distr(kk)%numc > lim ) THEN
                zvpart = sum(distr(kk)%volc(n1:nn))/distr(kk)%numc
                within_bins = distr(kk)%vlolim<zvpart .AND. zvpart<distr(kk)%vhilim
            END IF

        END DO ! - kk

        iter = iter + 1
        IF (iter > 100) STOP 'Error: Bin update not converged'

    ENDDO

  END SUBROUTINE bin_adjust

END MODULE mo_salsa_update
