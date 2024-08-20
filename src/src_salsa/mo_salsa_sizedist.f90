MODULE mo_salsa_sizedist
   IMPLICIT NONE

CONTAINS

   SUBROUTINE size_distribution( kproma, kbdim, klev, nmod, nstr, nend,  &
                                 n, dpg, sigmag, naero)

      USE mo_submctl, ONLY :      &
         pi6,                       &
         pi,                        &
         fn2b

      USE mo_salsa_types, ONLY : aero ! This is needed for size bins spacings
      IMPLICIT NONE

      !INTEGER, PARAMETER :: nmod = 7
      INTEGER, INTENT(in) :: nmod
      INTEGER, INTENT(in) :: nstr, nend  ! Start and end bin indices
      
      INTEGER, INTENT(IN) ::      &          
         kproma,                    & ! number of horiz. grid points
         kbdim,                     & ! dimension for arrays
         klev                         ! number of vertical levels

      REAL, INTENT(IN) ::         &
         n(nmod),                   & ! total concentration of a mode
         dpg(nmod),                 & ! geometric-mean diameter of a mode
         sigmag(nmod)                 ! standard deviation of a mode

      REAL, INTENT(OUT) ::        &
         naero(kbdim,klev,fn2b)      ! number concentration  [#/m3]

      !-- local variables
      REAL ::                     &
         deltadp                      ! bin width [m]

      REAL :: d1, d2, delta_d, dmid

      INTEGER :: ii, jj, kk, ib

      naero = 0.
      

      
      DO jj = 1, klev    ! vertical grid
         DO ii = 1, kbdim ! horizontal grid

            DO kk = nstr, nend
               naero(ii,jj,kk) = 0.0

               d1 = (aero(ii,jj,kk)%vlolim/pi6)**(1./3.)
               d2 = (aero(ii,jj,kk)%vhilim/pi6)**(1./3.)
               delta_d = (d2-d1)/10
               DO ib = 1, 10
                  d1 = (aero(ii,jj,kk)%vlolim/pi6)**(1./3.)+(REAL(ib)-1.)*delta_d
                  d2 = d1+delta_d
                  dmid = (d1+d2)/2
                  deltadp = log(d2/d1)
                
                  !-- size distribution
                  !   ntot = total number, total area, or total volume concentration
                  !   dpg = geometric-mean number, area, or volume diameter
                  !   n(kk) = number, area, or volume concentration in a bin
                  naero(ii,jj,kk) = naero(ii,jj,kk)+sum(n*deltadp/                        &
                                                        (sqrt(2.*pi)*log(MAX(sigmag,1.e-9)))*        &
                                                        exp(-log(dmid/MAX(dpg,1.e-9))**2/(2.*log(MAX(sigmag,1.e-9))**2)))
               END DO

            END DO

         END DO
      END DO

   END SUBROUTINE size_distribution

END MODULE mo_salsa_sizedist
