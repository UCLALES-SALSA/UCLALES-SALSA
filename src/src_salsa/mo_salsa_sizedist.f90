MODULE mo_salsa_sizedist

CONTAINS

  SUBROUTINE size_distribution(kproma, kbdim, klev, &
       n, dpg, sigmag, naero)

    USE mo_submctl, ONLY :      &
         nreg,                      &
         pi6,                       &
         pi,                        &
         in1a,                      &
         in2b,                      &
         fn2a,                      &
         fn2b, listspec

    USE mo_salsa_driver, ONLY : aero ! This is needed for size bins spacings
!    USE grid, ONLY : prtcl !! this should be avoided ! #mod
    
!    USE class_componentIndex, ONLY : IsUsed !#mod
    IMPLICIT NONE

    INTEGER, PARAMETER :: nmod = 7

    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid points
         kbdim,                     & ! dimension for arrays
         klev                         ! number of vertical levels

    REAL, INTENT(IN) ::         &
         n(nmod)                  , & ! total concentration of a mode
         dpg(nmod)                , & ! geometric-mean diameter of a mode
         sigmag(nmod)                 ! standard deviation of a mode

    REAL, INTENT(OUT) ::        &
         naero(kbdim,klev,fn2b)      ! number concentration  [#/m3]

    !-- local variables
    REAL ::                     &
         deltadp                      ! bin width [m]

    REAL :: d1,d2,delta_d,dmid

    INTEGER :: ii, jj, kk,ib

    naero = 0.

    DO jj = 1,klev    ! vertical grid
       DO ii = 1,kproma ! horizontal grid

          DO kk = in1a, fn2a
             naero(ii,jj,kk) = 0.0

             d1 = (aero(ii,jj,kk)%vlolim/pi6)**(1./3.)
             d2 = (aero(ii,jj,kk)%vhilim/pi6)**(1./3.)
             delta_d=(d2-d1)/10
             DO ib = 1,10
                d1=(aero(ii,jj,kk)%vlolim/pi6)**(1./3.)+(ib-1.)*delta_d
                d2=d1+delta_d
                dmid=(d1+d2)/2
                deltadp = log(d2/d1)
                
                !deltadp = (aero(ii,jj,kk)%vhilim**(1./3.)-aero(ii,jj,kk)%vlolim**(1./3.))/   &
                !     pi6**(1./3.)
                
                !-- size distribution
                !   ntot = total number, total area, or total volume concentration
                !   dpg = geometric-mean number, area, or volume diameter
                !   n(kk) = number, area, or volume concentration in a bin
!                naero(ii,jj,kk) = naero(ii,jj,kk)+sum(n*deltadp/                        &
!                     (sqrt(2.*pi)*log(sigmag))*                 &
!                     exp(-log(aero(ii,jj,kk)%dmid/dpg)**2/(2.*log(sigmag)**2)))
                  naero(ii,jj,kk) = naero(ii,jj,kk)+sum(n*deltadp/                        &
                     (sqrt(2.*pi)*log(sigmag))*                 &
                     exp(-log(dmid/dpg)**2/(2.*log(sigmag)**2)))
                

 !                write(*,*) naero(ii,jj,kk)*1e-6,d1,d2,dmid
             END DO
 !            write(*,*) dpg(1:2),n(1:2)*1e-6
 !            write(*,*) kk,naero(ii,jj,kk)*1e-6,aero(ii,jj,kk)%dmid, (aero(ii,jj,kk)%vlolim/pi6)**(1./3.), &
 !                 (aero(ii,jj,kk)%vhilim/pi6)**(1./3.)
 !            write(*,*) sum(naero(ii,jj,:))
 !            write(*,*) 'here'
 !            pause
          END DO

          IF ( ANY( listspec=='BC' ) .OR. ANY( listspec=='DU' ) ) THEN
          ! insolubles
            DO kk = in2b, fn2b !!! #mod

             d1 = (aero(ii,jj,kk)%vlolim/pi6)**(1./3.)
             d2 = (aero(ii,jj,kk)%vhilim/pi6)**(1./3.)

             deltadp = log(d2/d1)

             !deltadp = (aero(ii,jj,kk)%vhilim**(1./3.)-aero(ii,jj,kk)%vlolim**(1./3.))/   &
             !     pi6**(1./3.)

             !-- size distribution
             !   ntot = total number, total area, or total volume concentration
             !   dpg = geometric-mean number, area, or volume diameter
             !   n(kk) = number, area, or volume concentration in a bin
             naero(ii,jj,kk) = sum(n*deltadp/                        &
                  (sqrt(2.*pi)*log(sigmag))*                 &
                  exp(-log(aero(ii,jj,kk)%dmid/dpg)**2/(2.*log(sigmag)**2)))

            END DO
          END IF
       END DO

    END DO

  END SUBROUTINE size_distribution

END MODULE mo_salsa_sizedist
