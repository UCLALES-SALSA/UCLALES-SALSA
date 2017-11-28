!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
! MODIFICATIONS:
!
! Added new thermodynamics subroutine to be used with SALSA *SALSAthrm*.
! It basically calculates the water vapour pressure, allowing supersaturation.
! SALSA provides liquid water content for cloud droplets and rain
! by calculating activation, condensation, autoconversion etc., providing
! the corresponding sink/source term for ambient vapour mixing ratio.
!
! Juha Tonttila, FMI, 2014
!
!
!
MODULE thrm

  IMPLICIT NONE

CONTAINS

!
! -------------------------------------------------------------------------
! THERMO: calculates thermodynamics quantities according to level.  Level
! is passed in to allow level of diagnosis to be determined by call rather
! than by runtype
!
  SUBROUTINE thermo (level)

    USE grid, ONLY : a_rc, a_rv, a_rh, a_theta, a_pexnr, a_press, a_temp,  &
                     a_rsl, a_rp, a_tp, nxp, nyp, nzp, th00, pi0, pi1,a_rpp,   &
                     a_srp, a_ri, a_rsi, a_rhi, a_srs

    INTEGER, INTENT (in) :: level

    SELECT CASE (level)
    CASE DEFAULT
       CALL drythrm(nzp,nxp,nyp,a_pexnr,a_press,a_tp,a_theta,a_temp,pi0,   &
                    pi1,th00,a_rp,a_rv)
    CASE (2)
       CALL satadjst(nzp,nxp,nyp,a_pexnr,a_press,a_tp,a_theta,a_temp,pi0,  &
                     pi1,th00,a_rp,a_rv,a_rc,a_rsl)
    CASE (3)
       CALL satadjst3(nzp,nxp,nyp,a_pexnr,a_press,a_tp,a_theta,a_temp,pi0, &
                      pi1,th00,a_rp,a_rv,a_rc,a_rsl,a_rpp)
    CASE (4:5)
       CALL SALSAthrm(level,nzp,nxp,nyp,a_pexnr,pi0,pi1,th00,a_rp,a_tp,a_theta, &
                      a_temp,a_press,a_rsl,a_rh,a_rc,a_srp,a_ri,a_rsi,a_rhi,a_srs)
    END SELECT

  END SUBROUTINE thermo
!
! -------------------------------------------------------------------------
! update_pi1:  this routine updates a pressure associated with the
! subtraction of a mean acceleration, only incrementing it for dynamic and
! thermal effects for layers above the surface
!
  SUBROUTINE update_pi1(n1,awtbar,pi1)

    USE grid, ONLY : th00, zt

    INTEGER, INTENT (in) :: n1
    REAL, INTENT (in) , DIMENSION (n1)   :: awtbar
    REAL, INTENT (inout), DIMENSION (n1) :: pi1

    INTEGER :: k

    DO k = 2, n1
       pi1(k) = pi1(k-1) + awtbar(k-1)*(zt(k)-zt(k-1))/th00
    END DO

  END SUBROUTINE update_pi1

  !
  !----------------------------------------------------------------------
  ! SALSAthrm: Calculates potential and absolute temperatures, pressure,
  !            and total cloud/rain water mixing ratios with microphysics
  !            provided by the SALSA model. NOTE, no saturation adjustment
  !            takes place -> the resulting water vapour mixing ratio
  !            can be supersaturated, allowing the microphysical calculations
  !            in SALSA.
  !

  SUBROUTINE SALSAthrm(level,n1,n2,n3,pp,pi0,pi1,th00,rv,tl,th,tk,p,rs,rh,rc,srp,ri,rsi,rhi,srs)
    USE defs, ONLY : R, cp, cpr, p00, alvl, alvi
    USE grid, ONLY : a_dn
    IMPLICIT NONE

    INTEGER, INTENT(in) :: level,n1,n2,n3
    REAL, INTENT(in) :: rv(n1,n2,n3),   &     ! Water vapour mixing ratio
                        pp(n1,n2,n3),   &     ! Exner function
                        tl(n1,n2,n3)          ! liquid potential temp
    REAL, INTENT(in) :: th00
    REAL, INTENT(in) :: pi0(n1),pi1(n1)
    REAL, INTENT(IN) :: rc(n1,n2,n3),   &  ! Total cloud condensate mix rat
                        srp(n1,n2,n3)            ! Precipitation mix rat
    REAL, INTENT(OUT) :: rs(n1,n2,n3),  &   ! Saturation mix rat
                         rh(n1,n2,n3),  &     ! Relative humidity
                         th(n1,n2,n3),  &     ! Potential temperature
                         tk(n1,n2,n3),  &     ! Absolute temperature
                         p(n1,n2,n3)          ! Air pressure
    REAL, INTENT(IN)  :: ri(n1,n2,n3),   &  ! Total ice condensate mix rat
                         srs(n1,n2,n3)      ! Snow mix rat
    REAL, INTENT(OUT) :: rsi(n1,n2,n3), &  ! Saturation mix rat over ice
                         rhi(n1,n2,n3)     ! relative humidity over ice
    REAL    :: exner
    INTEGER :: k,i,j
    REAL    :: thil

     DO j = 3, n3-2
       DO i = 3, n2-2
          DO k = 1, n1

             ! Pressure
             exner = (pi0(k) + pi1(k) + pp(k,i,j))/cp
             p(k,i,j) = p00*exner**cpr
             thil = tl(k,i,j)+th00

             ! Potential and absolute temperatures
             th(k,i,j) = thil + (alvl*( rc(k,i,j) + srp(k,i,j) ))/cp/exner
             IF(level == 5) th(k,i,j) = th(k,i,j) + (alvi*( ri(k,i,j)+ srs(k,i,j) ))/cp/exner

             tk(k,i,j) = th(k,i,j)*exner

             ! Saturation mixing ratio
             rs(k,i,j) = rslf(p(k,i,j),tk(k,i,j))
             rh(k,i,j) = rv(k,i,j)/rs(k,i,j)

             IF (level==5) THEN
                 rsi(k,i,j) = rsif(p(k,i,j),tk(k,i,j))
                 rhi(k,i,j) = rv(k,i,j)/rsi(k,i,j)
             END IF

             ! True air density
             a_dn(k,i,j) = p(k,i,j)/(R*tk(k,i,j))

          END DO
       END DO
    END DO

  END SUBROUTINE SALSAthrm
!
! -------------------------------------------------------------------------
! DRYTHRM:  this routine calculates theta, and pressure for
! the case when no moisture is present
!
  SUBROUTINE drythrm(n1,n2,n3,pp,p,thil,theta,t,pi0,pi1,th00,rt,rv)

  USE defs, ONLY : cp, cpr, p00

  INTEGER, INTENT (in) :: n1,n2,n3
  REAL, INTENT (in)    :: pi0(n1),pi1(n1),th00
  REAL, INTENT (in)    :: pp(n1,n2,n3),thil(n1,n2,n3),rt(n1,n2,n3)
  REAL, INTENT (out)   :: p(n1,n2,n3),theta(n1,n2,n3),rv(n1,n2,n3),t(n1,n2,n3)

  INTEGER :: i,j,k
  REAL    :: exner

  DO j = 3, n3-2
    DO i = 3, n2-2
      DO k = 1, n1
        exner  = (pi0(k)+pi1(k)+pp(k,i,j))/cp
        p(k,i,j) = p00 * (exner)**cpr
        theta(k,i,j) = thil(k,i,j)+th00
        t(k,i,j) = theta(k,i,j)*exner
        rv(k,i,j) = rt(k,i,j)
      END DO
    END DO
  END DO

  END SUBROUTINE drythrm
!
! -------------------------------------------------------------------------
! SATADJST:  this routine calculates theta, and pressure and diagnoses
! liquid water using a saturation adjustment for warm-phase systems
!
  SUBROUTINE satadjst(n1,n2,n3,pp,p,tl,th,tk,pi0,pi1,th00,rt,rv,rc,rs)

    USE defs, ONLY : cp, cpr, alvl, ep, Rm, p00

    INTEGER, INTENT (in) ::  n1,n2,n3

    REAL, INTENT (in), DIMENSION (n1,n2,n3)  :: pp, tl, rt
    REAL, INTENT (in), DIMENSION (n1)        :: pi0, pi1
    REAL, INTENT (in)                        :: th00
    REAL, INTENT (out), DIMENSION (n1,n2,n3) :: rc,rv,rs,th,tk,p

    INTEGER :: k, i, j, iterate
    REAL    :: exner,til,x1,xx,yy,zz

    DO j = 3, n3-2
       DO i = 3, n2-2
          DO k = 1, n1
             exner = (pi0(k)+pi1(k)+pp(k,i,j))/cp
             p(k,i,j) = p00 * (exner)**cpr
             til = (tl(k,i,j)+th00)*exner
             xx = til
             yy = rslf(p(k,i,j),xx)
             zz = max(rt(k,i,j)-yy,0.)
             IF (zz > 0.) THEN
                DO iterate = 1, 3
                   x1 = alvl/(cp*xx)
                   xx = xx - (xx - til*(1.+x1*zz))/(1. + x1*til                &
                        *(zz/xx+(1.+yy*ep)*yy*alvl/(Rm*xx*xx)))
                   yy = rslf(p(k,i,j),xx)
                   zz = max(rt(k,i,j)-yy,0.)
                END DO
             END IF
             rc(k,i,j) = zz
             rv(k,i,j) = rt(k,i,j)-rc(k,i,j)
             rs(k,i,j) = yy
             tk(k,i,j) = xx
             th(k,i,j) = tk(k,i,j)/exner
          END DO
       END DO
    END DO

  END SUBROUTINE satadjst
!
! -------------------------------------------------------------------------
! SATADJST3:  this routine calculates theta, and pressure and diagnoses
! liquid water using a saturation adjustment for warm-phase systems; in
! addition, takes in the account the precipitable water when present
!
  SUBROUTINE satadjst3(n1,n2,n3,pp,p,tl,th,tk,pi0,pi1,th00,rt,rv,rc,rs,rp)

    USE defs, ONLY : cp, cpr, alvl, ep, Rm, p00
    USE mpi_interface, ONLY : myid, appl_abort

    INTEGER, INTENT (in) ::  n1,n2,n3

    REAL, INTENT (in), DIMENSION (n1,n2,n3)  :: pp, tl, rt, rp
    REAL, INTENT (in), DIMENSION (n1)        :: pi0, pi1
    REAL, INTENT (in)                        :: th00
    REAL, INTENT (out), DIMENSION (n1,n2,n3) :: rc, rv, rs, th, tk, p

    INTEGER :: k, i, j, iterate
    REAL    :: exner, tli, tx, txi, rsx, rcx, rpc, tx1, dtx
    REAL, PARAMETER :: epsln = 1.e-4

    DO j = 3, n3-2
       DO i = 3, n2-2
          DO k = 1, n1
             exner = (pi0(k)+pi1(k)+pp(k,i,j))/cp
             p(k,i,j) = p00 * (exner)**cpr
             tli = (tl(k,i,j)+th00)*exner
             tx = tli
             rsx = rslf(p(k,i,j),tx)
             rcx = max(rt(k,i,j)-rsx,0.)
             rpc = rp(k,i,j)
             IF (rcx > 0. .OR. rpc > 0.) THEN
                iterate = 1
                dtx = 1.
                IF (rcx < rpc) THEN
                   DO WHILE (dtx > epsln .AND. iterate < 10)
                      txi = alvl*rpc/(cp*tx)
                      tx1 = tx - (tx - tli*(1+txi)) / (1+txi*tli/tx)
                      dtx = abs(tx1-tx)
                      tx  = tx1
                      iterate = iterate+1
                   END DO
                   rsx = rslf(p(k,i,j),tx)
                   rcx = max(rt(k,i,j)-rsx,0.)
                ELSE
                   DO WHILE(dtx > epsln .AND. iterate < 10)
                      txi = alvl/(cp*tx)
                      tx1 = tx - (tx - tli*(1.+txi*rcx))/(1. + txi*tli         &
                             *(rcx/tx+(1.+rsx*ep)*rsx*alvl/(Rm*tx*tx)))
                      dtx = abs(tx1-tx)
                      tx  = tx1
                      rsx = rslf(p(k,i,j),tx)
                      rcx = max(rt(k,i,j)-rsx,0.)
                      iterate = iterate+1
                   END DO
                END IF

                IF (dtx > epsln) THEN
                   IF (myid == 0) PRINT *, '  ABORTING: thrm', dtx, epsln
                   IF (myid == 0) WRITE(*,*) pp(k,i,j),p(k,i,j),tl(k,i,j),th(k,i,j), &
                                             tk(k,i,j),rt(k,i,j),rv(k,i,j),rc(k,i,j),rs(k,i,j),rp(k,i,j)
                   CALL appl_abort(0)
                END IF
             END IF
             rc(k,i,j) = rcx
             rv(k,i,j) = rt(k,i,j)-rc(k,i,j)
             rs(k,i,j) = rsx
             tk(k,i,j) = tx
             th(k,i,j) = tk(k,i,j)/exner
          END DO
       END DO
    END DO

  END SUBROUTINE satadjst3
!
! ---------------------------------------------------------------------
! This function calculates the liquid saturation vapor mixing ratio as
! a function of temperature and pressure
!
  REAL FUNCTION rslf(p,t)

  REAL, INTENT (in) :: p, t
  REAL, PARAMETER   :: c0 = 0.6105851e+03, c1 = 0.4440316e+02,    &
                       c2 = 0.1430341e+01, c3 = 0.2641412e-01,    &
                       c4 = 0.2995057e-03, c5 = 0.2031998e-05,    &
                       c6 = 0.6936113e-08, c7 = 0.2564861e-11,    &
                       c8 = -.3704404e-13

  REAL ::  esl, x

  x = min(max(-80.,t-273.16),50.)
  esl = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
  rslf = .622*esl/(p-esl)

  END FUNCTION rslf
!
! ---------------------------------------------------------------------
! This function calculates the ice saturation vapor mixing ratio as a
! function of temperature and pressure
!
  REAL FUNCTION rsif(p,t)

  REAL, INTENT (in) :: p, t
  REAL, PARAMETER   :: c0 = 0.6114327e+03, c1 = 0.5027041e+02,    &
                       c2 = 0.1875982e+01, c3 = 0.4158303e-01,    &
                       c4 = 0.5992408e-03, c5 = 0.5743775e-05,    &
                       c6 = 0.3566847e-07, c7 = 0.1306802e-09,    &
                       c8 = 0.2152144e-12

  REAL  :: esi, x

  x = max(-80.,t-273.16)
  esi = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
  rsif = .622*esi/(p-esi)

  END FUNCTION rsif
!
! -------------------------------------------------------------------------
! FLL_TKRS: Updates scratch arrays with temperature and saturation mixing
! ratio
!
  SUBROUTINE fll_tkrs(n1,n2,n3,th,pp,pi0,pi1,tk,rs)

  USE defs, ONLY : cp, cpr, p00

  INTEGER, INTENT (in) :: n1,n2,n3
  REAL, INTENT (in)    :: th(n1,n2,n3), pp(n1,n2,n3)
  REAL, INTENT (in)    :: pi0(n1), pi1(n1)
  REAL, INTENT (out)   :: tk(n1,n2,n3)
  REAL, OPTIONAL, INTENT (out) :: rs(n1,n2,n3)

  INTEGER :: i, j, k
  REAL    :: exner

  DO j = 3, n3-2
    DO i = 3, n2-2
      DO k = 1, n1
        exner = (pi0(k)+pi1(k)+pp(k,i,j))/cp
        tk(k,i,j) = th(k,i,j)*exner
        IF (present(rs)) rs(k,i,j) = rslf(p00*exner**cpr,tk(k,i,j))
      END DO
    END DO
  END DO

  END SUBROUTINE fll_tkrs
!
! -------------------------------------------------------------------------
! BRUVAIS:  Calculates the brunt-vaisaila frequency in accordance with the
! thermodynamic level
!
! Modified for level 4,
! Juha Tonttila, FMI, 2014
!
  SUBROUTINE bruvais(n1,n2,n3,level,th,tl,rt,rs,en2,dzm,th00)

  USE defs, ONLY : g, R, cp, alvl, ep, ep2

  INTEGER, INTENT (in) :: n1, n2, n3, level
  REAL, INTENT (in)    :: th(n1,n2,n3), tl(n1,n2,n3), rt(n1,n2,n3),         &
                          rs(n1,n2,n3), dzm(n1), th00
  REAL, INTENT (out)   :: en2(n1,n2,n3)

  INTEGER :: i, k, j, kp1
  REAL    :: c1, c2, c3, tvk, tvkp1, rtbar, rsbar, aa, bb

  DO j = 3, n3-2
     DO i = 3, n2-2
        DO k = 1, n1-1
           SELECT CASE(level)
           CASE (0)
              en2(k,i,j) = g*dzm(k)*((th(k+1,i,j)-th(k,i,j))/th00)
           CASE (1)
              tvk = th(k,i,j)*(1.+ep2*rt(k,i,j))
              tvkp1 = th(k+1,i,j)*(1.+ep2*rt(k+1,i,j))
              en2(k,i,j) = g*dzm(k)*(tvkp1-tvk)/th00
           CASE (2)
              rtbar = 0.5*(rt(k,i,j)+rt(k+1,i,j))
              rsbar = 0.5*(rs(k,i,j)+rs(k+1,i,j))
              kp1 = min(n1-1,k+1)
              IF (rt(k,i,j) > rs(k,i,j) .AND. rt(kp1,i,j) > rs(kp1,i,j)) THEN
                 c1=(1.+ep*alvl/R/th(k,i,j))/ep
                 c2=ep*alvl*alvl/(R*cp*th(k,i,j)*th(k,i,j))
                 c3=alvl/(cp*th(k,i,j))
                 aa = (1. - rtbar + rsbar*c1)/(1. + c2*rsbar)
                 bb = (c3*aa - 1.)
              ELSE
                 aa = (1.00 + ep2*rtbar)
                 bb = ep2
              END IF
              en2(k,i,j) = g*dzm(k)*(aa*(tl(k+1,i,j)-tl(k,i,j))/th00        &
                           + bb*(rt(k+1,i,j)-rt(k,i,j)))
           CASE (3,4,5)
              rtbar = 0.5*(rt(k,i,j)+rt(k+1,i,j))
              rsbar = 0.5*(rs(k,i,j)+rs(k+1,i,j))
              kp1 = min(n1-1,k+2)
              IF (rt(k,i,j) > rs(k,i,j) .AND. rt(kp1,i,j) > rs(kp1,i,j)) THEN
                 c1=(1.+ep*alvl/R/th(k,i,j))/ep
                 c2=ep*alvl*alvl/(R*cp*th(k,i,j)*th(k,i,j))
                 c3=alvl/(cp*th(k,i,j))
                 aa = (1. - rtbar + rsbar*c1)/(1. + c2*rsbar)
                 bb = (c3*aa - 1.)
              ELSE
                 aa = (1.00 + ep2*rtbar)
                 bb = ep2
              END IF
              en2(k,i,j) = g*dzm(k)*(aa*(tl(k+1,i,j)-tl(k,i,j))/th00        &
                           + bb*(rt(k+1,i,j)-rt(k,i,j)))
           CASE DEFAULT
              STOP 'level not supported in bruvais'
           END SELECT
        END DO
        en2(n1,i,j) = en2(n1-1,i,j)
     END DO
  END DO

  END SUBROUTINE bruvais

END MODULE thrm
