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
module thrm

  implicit none

contains

!
! -------------------------------------------------------------------------
! THERMO: calculates thermodynamics quantities according to level.  Level
! is passed in to allow level of diagnosis to be determined by call rather
! than by runtype
!
  subroutine thermo (level)

    use grid, only : a_rc, a_rv, a_theta, a_pexnr, a_press, a_temp,  &
         a_rsl, a_rp, a_tp, nxp, nyp, nzp, th00, pi0, pi1,a_rpp,   &
         a_maerop, a_mcloudp, a_mprecpp, a_micep, a_msnowp, &
         nbins, ncld, nprc, nice, nsnw, &
         a_ri, a_rsi, a_dn, a_rip, a_rsp, a_rgp, a_rhp
    USE defs, ONLY : Rd

    integer, intent (in) :: level

    select case (level)
    case (0)
       a_ri = a_rip + a_rsp + a_rgp + a_rhp ! Total ice+snow+graupel+hail
       call satadjst(nzp,nxp,nyp,a_pexnr,a_press,a_tp,a_theta,a_temp,pi0,  &
                     pi1,th00,a_rp,a_rv,a_rc,a_rsl,a_rpp,a_ri,a_rsi)
    case (1)
       call drythrm(nzp,nxp,nyp,a_pexnr,a_press,a_tp,a_theta,a_temp,pi0,   &
                    pi1,th00,a_rp,a_rv)
    case (2)
       call satadjst(nzp,nxp,nyp,a_pexnr,a_press,a_tp,a_theta,a_temp,pi0,  &
                     pi1,th00,a_rp,a_rv,a_rc,a_rsl)
    case (3)
       call satadjst(nzp,nxp,nyp,a_pexnr,a_press,a_tp,a_theta,a_temp,pi0,  &
                     pi1,th00,a_rp,a_rv,a_rc,a_rsl,a_rpp)
    case (4:5)
       ! Update diagnostic variables: total liquid and ice
       a_rc(:,:,:) = SUM(a_maerop(:,:,:,1:nbins),DIM=4) + &
                     SUM(a_mcloudp(:,:,:,1:ncld),DIM=4) + &
                     SUM(a_mprecpp(:,:,:,1:nprc),DIM=4)
       a_ri(:,:,:) = SUM(a_micep(:,:,:,1:nice),DIM=4) + &
                     SUM(a_msnowp(:,:,:,1:nsnw),DIM=4)

       CALL SALSAthrm(level,nzp,nxp,nyp,a_pexnr,pi0,pi1,th00,a_tp,a_theta, &
                      a_temp,a_press,a_rsl,a_rc,a_ri,a_rsi)
    end select

    ! Air density
    a_dn(:,3:nxp-2,3:nyp-2) = a_press(:,3:nxp-2,3:nyp-2)/(Rd*a_temp(:,3:nxp-2,3:nyp-2))

  end subroutine thermo
!
! -------------------------------------------------------------------------
! update_pi1:  this routine updates a pressure associated with the
! subtraction of a mean acceleration, only incrementing it for dynamic and
! thermal effects for layers above the surface
!
  subroutine update_pi1(n1,awtbar,pi1)

    use grid, only : th00, zt

    integer, intent (in) :: n1
    real, intent (in) , dimension (n1) :: awtbar
    real, intent (inout), dimension (n1) :: pi1

    integer :: k

    do k=2,n1
       pi1(k) = pi1(k-1) + awtbar(k-1)*(zt(k)-zt(k-1))/th00
    end do

  end subroutine update_pi1

!
!----------------------------------------------------------------------
! SALSAthrm: Calculates potential and absolute temperatures, pressure,
!            and total cloud/rain water mixing ratios with microphysics
!            provided by the SALSA model. NOTE, no saturation adjustment
!            takes place -> the resulting water vapour mixing ratio
!            can be supersaturated, allowing the microphysical calculations
!            in SALSA.
!

  SUBROUTINE SALSAthrm(level,n1,n2,n3,pp,pi0,pi1,th00,tl,th,tk,p,rs,rc,ri,rsi)
    USE defs, ONLY : cp, cpr, p00, alvl, alvi
    IMPLICIT NONE

    INTEGER, INTENT(in) :: level,n1,n2,n3
    REAL, INTENT(in) :: pp(n1,n2,n3),pi0(n1),pi1(n1)
    REAL, INTENT(in) :: th00,tl(n1,n2,n3)
    REAL, INTENT(IN) :: rc(n1,n2,n3), ri(n1,n2,n3) ! Total liquid and ice mixing ratios
    REAL, INTENT(OUT) :: rs(n1,n2,n3),  &   ! Saturation mix rat
                         rsi(n1,n2,n3), &   ! Saturation mixing rat over ice
                         th(n1,n2,n3),  &     ! Potential temperature
                         tk(n1,n2,n3),  &     ! Absolute temperature
                         p(n1,n2,n3)           ! Air pressure
    REAL :: exner
    INTEGER :: k,i,j
    REAL :: thil

     DO j = 3,n3-2
       DO i = 3,n2-2
          DO k = 1,n1

             ! Pressure
             exner = (pi0(k) + pi1(k) + pp(k,i,j))/cp
             p(k,i,j) = p00*exner**cpr
             thil = tl(k,i,j)+th00

             ! Potential and absolute temperatures

             th(k,i,j) = thil + rc(k,i,j)*alvl/(cp*exner)

             if(level==5) th(k,i,j) = th(k,i,j) + ri(k,i,j)*alvi/(cp*exner)

             tk(k,i,j) = th(k,i,j)*exner

             ! Saturation mixing ratio
             rs(k,i,j) = rslf(p(k,i,j),tk(k,i,j))
             rsi(k,i,j) = rsif(p(k,i,j),tk(k,i,j))

          END DO
       END DO
    END DO

  END SUBROUTINE SALSAthrm
!
! -------------------------------------------------------------------------
! DRYTHRM:  this routine calculates theta, and pressure for
! the case when no moisture is present
!
  subroutine drythrm(n1,n2,n3,pp,p,thil,theta,t,pi0,pi1,th00,rt,rv)

  use defs, only : cp, cpr, p00

  integer, intent (in) :: n1,n2,n3
  real, intent (in)    :: pi0(n1),pi1(n1),th00
  real, intent (in)    :: pp(n1,n2,n3),thil(n1,n2,n3),rt(n1,n2,n3)
  real, intent (out)   :: p(n1,n2,n3),theta(n1,n2,n3),rv(n1,n2,n3),t(n1,n2,n3)

  integer :: i,j,k
  real    :: exner

  do j=3,n3-2
    do i=3,n2-2
      do k=1,n1
        exner  = (pi0(k)+pi1(k)+pp(k,i,j))/cp
        p(k,i,j) = p00 * (exner)**cpr
        theta(k,i,j)=thil(k,i,j)+th00
        t(k,i,j)=theta(k,i,j)*exner
        rv(k,i,j)=rt(k,i,j)
      enddo
    enddo
  enddo

  end subroutine drythrm
!
! -------------------------------------------------------------------------
! SATADJST:  this routine calculates theta, and pressure and diagnoses
! liquid water using a saturation adjustment for warm-phase systems; in
! addition, takes in the account the precipitable water when present
!
  subroutine satadjst(n1,n2,n3,pp,p,tl,th,tk,pi0,pi1,th00,rt,rv,rc,rs,rp,ri,rsi)

    use defs, only : cp, cpr, alvl, alvi, ep, Rm, p00
    use mpi_interface, only : appl_abort

    integer, intent (in) ::  n1,n2,n3

    real, intent (in), dimension (n1,n2,n3)  :: pp, tl, rt
    real, intent (in), dimension (n1)        :: pi0, pi1
    real, intent (in)                        :: th00
    real, intent (out), dimension (n1,n2,n3) :: rc, rv, rs, th, tk, p
    real, intent (in), optional, dimension (n1,n2,n3) :: rp, ri
    real, intent (out), optional, dimension (n1,n2,n3) :: rsi

    integer :: k, i, j, iterate
    real    :: exner, tli, tx, txi, rsx, rcx, rpc, tx1, dtx, ric
    real, parameter :: epsln = 1.e-4

    rpc = 0.
    ric = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             IF (PRESENT(rp)) rpc = MAX(0.,rp(k,i,j))
             IF (PRESENT(ri)) ric = MAX(0.,ri(k,i,j))
             exner=(pi0(k)+pi1(k)+pp(k,i,j))/cp
             p(k,i,j) = p00 * (exner)**cpr
             ! Adjust tl and rt so that they are for cloud water only
             tli=(tl(k,i,j)+th00)*exner+alvl/cp*rpc+alvi/cp*ric
             tx=tli
             rsx=rslf(p(k,i,j),tx) ! Saturation mixing ratio
             rcx=max(rt(k,i,j)-rpc-ric-rsx,0.) ! Cloud condensate mixing ratio
             if (rcx > 0.) then
                dtx = 2.*epsln
                iterate = 1
                do while(dtx > epsln .and. iterate < 20)
                   txi=alvl/(cp*tx)
                   tx1=tx - (tx - tli*(1.+txi*rcx))/(1. + txi*tli         &
                        *(rcx/tx+(1.+rsx*ep)*rsx*alvl/(Rm*tx*tx)))
                   dtx = abs(tx1-tx)
                   tx  = tx1
                   rsx=rslf(p(k,i,j),tx)
                   rcx=max(rt(k,i,j)-rpc-ric-rsx,0.)
                   iterate = iterate+1
                enddo
                if (dtx > epsln) then
                    print *, '  ABORTING: thrm', dtx, epsln
                    WRITE(*,*) pp(k,i,j),p(k,i,j),tl(k,i,j),th(k,i,j), &
                               tk(k,i,j),rt(k,i,j),rv(k,i,j),rc(k,i,j),rs(k,i,j)
                    call appl_abort(0)
                endif
             endif
             rc(k,i,j)=rcx
             rv(k,i,j)=rt(k,i,j)-rpc-ric-rc(k,i,j)
             rs(k,i,j)=rsx
             tk(k,i,j)=tx
             th(k,i,j)=tk(k,i,j)/exner
             IF (present(rsi)) rsi(k,i,j) = rsif(p(k,i,j),tk(k,i,j))
          enddo
       enddo
    enddo

  end subroutine satadjst
!
! ---------------------------------------------------------------------
! This function calculates the water saturation vapor mixing ratio as
! a function of temperature and pressure
!
  real function rslf(p,t)
  real, intent (in) :: p, t
  real ::  e
  e=esl(t)
  rslf=.622*e/(p-e)
  end function rslf

  real elemental function esl(t)
  real, intent (in) :: t
  real, parameter :: c0=0.6105851e+03, c1=0.4440316e+02,    &
                     c2=0.1430341e+01, c3=0.2641412e-01,    &
                     c4=0.2995057e-03, c5=0.2031998e-05,    &
                     c6=0.6936113e-08, c7=0.2564861e-11,    &
                     c8=-.3704404e-13
  real  :: x
  x=min(max(-80.,t-273.16),50.)
  esl=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
  end function esl
!
! ---------------------------------------------------------------------
! This function calculates the ice saturation vapor mixing ratio as a
! function of temperature and pressure
!
  real function rsif(p,t)
  real, intent (in) :: p, t
  real  :: e
  e=esi(t)
  rsif=.622*e/(p-e)
  end function rsif

  real elemental function esi(t)
  real, intent (in) :: t
  real, parameter :: c0=0.6114327e+03, c1=0.5027041e+02,    &
                     c2=0.1875982e+01, c3=0.4158303e-01,    &
                     c4=0.5992408e-03, c5=0.5743775e-05,    &
                     c6=0.3566847e-07, c7=0.1306802e-09,    &
                     c8=0.2152144e-12
  real  :: x
  x=max(-80.,t-273.16)
  esi=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
  end function esi
!
! -------------------------------------------------------------------------
! FLL_TKRS: Updates scratch arrays with temperature and saturation mixing
! ratio
!
  subroutine fll_tkrs(n1,n2,n3,th,pp,pi0,pi1,tk,rs)

  use defs, only : cp, cpr, p00

  integer, intent (in) :: n1,n2,n3
  real, intent (in)    :: th(n1,n2,n3), pp(n1,n2,n3)
  real, intent (in)    :: pi0(n1), pi1(n1)
  real, intent (out)   :: tk(n1,n2,n3)
  real, optional, intent (out)   :: rs(n1,n2,n3)

  integer :: i, j, k
  real    :: exner

  do j=3,n3-2
    do i=3,n2-2
      do k=1,n1
        exner=(pi0(k)+pi1(k)+pp(k,i,j))/cp
        tk(k,i,j)=th(k,i,j)*exner
        if (present(rs)) rs(k,i,j)=rslf( p00*exner**cpr ,tk(k,i,j))
      end do
    end do
  end do

  end subroutine fll_tkrs
!
! -------------------------------------------------------------------------
! BRUVAIS:  Calculates the brunt-vaisaila frequency in accordance with the
! thermodynamic level
!
! Modified for level 4,
! Juha Tonttila, FMI, 2014
!
  subroutine bruvais(n1,n2,n3,level,th,tl,rt,rs,en2,dzm,th00)

  use defs, only : g, Rm, cp, alvl, ep, ep2

  integer, intent (in) ::  n1, n2, n3, level
  real, intent (in)    ::  th(n1,n2,n3), tl(n1,n2,n3), rt(n1,n2,n3),         &
                           rs(n1,n2,n3), dzm(n1), th00
  real, intent (out)   ::  en2(n1,n2,n3)

  integer :: i, k, j, kp1
  real    :: c1, c2, c3, tvk, tvkp1, rtbar, rsbar, aa, bb

  select case(level)
  case (1)
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1-1
              tvk=th(k,i,j)*(1.+ep2*rt(k,i,j))
              tvkp1=th(k+1,i,j)*(1.+ep2*rt(k+1,i,j))
              en2(k,i,j)=g*dzm(k)*(tvkp1-tvk)/th00
          end do
          en2(n1,i,j)=en2(n1-1,i,j)
       end do
    end do
  case (2)
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1-1
              rtbar=0.5*(rt(k,i,j)+rt(k+1,i,j))
              rsbar=0.5*(rs(k,i,j)+rs(k+1,i,j))
              kp1=min(n1-1,k+1)
              if (rt(k,i,j) > rs(k,i,j) .and. rt(kp1,i,j) > rs(kp1,i,j)) then
                 c1=(1.+alvl/Rm/th(k,i,j))/ep
                 c2=alvl*alvl/(Rm*cp*th(k,i,j)*th(k,i,j))
                 c3=alvl/(cp*th(k,i,j))
                 aa=(1. - rtbar + rsbar*c1)/(1. + c2*rsbar)
                 bb=(c3*aa - 1.)
              else
                 aa=(1.00 + ep2*rtbar)
                 bb=ep2
              end if
              en2(k,i,j)=g*dzm(k)*(aa*(tl(k+1,i,j)-tl(k,i,j))/th00        &
                   + bb*(rt(k+1,i,j)-rt(k,i,j)))
          end do
          en2(n1,i,j)=en2(n1-1,i,j)
       end do
    end do
  case (0,3,4,5)
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1-1
              rtbar=0.5*(rt(k,i,j)+rt(k+1,i,j))
              rsbar=0.5*(rs(k,i,j)+rs(k+1,i,j))
              kp1=min(n1-1,k+2)
              if (rt(k,i,j) > rs(k,i,j) .and. rt(kp1,i,j) > rs(kp1,i,j)) then
                 c1=(1.+alvl/Rm/th(k,i,j))/ep
                 c2=alvl*alvl/(Rm*cp*th(k,i,j)*th(k,i,j))
                 c3=alvl/(cp*th(k,i,j))
                 aa=(1. - rtbar + rsbar*c1)/(1. + c2*rsbar)
                 bb=(c3*aa - 1.)
              else
                 aa=(1.00 + ep2*rtbar)
                 bb=ep2
              end if
              en2(k,i,j)=g*dzm(k)*(aa*(tl(k+1,i,j)-tl(k,i,j))/th00        &
                   + bb*(rt(k+1,i,j)-rt(k,i,j)))
          end do
          en2(n1,i,j)=en2(n1-1,i,j)
       end do
    end do
  case default
     stop 'level not supported in bruvais'
  end select

  end subroutine bruvais

end module thrm
