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
module mcrp

  use defs, only : alvl, alvi, rowt, pi, Rm, cp, kb, g, vonk
  use grid, only : dtlt, dzt, nxp, nyp, nzp,a_pexnr, a_rp, a_tp, th00, CCN,     &
       dn0, pi0, a_rt, a_tt, a_rpp, a_rpt, a_npp, a_npt, a_rv, a_rc, a_theta,   &
       a_press, a_temp, a_rsl, precip, a_dn, a_ustar,                  &
       a_naerop,  a_naerot,  a_maerop,  a_maerot,                               &
       a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,                              &
       a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt, mc_ApVdom,         &
       a_nicep,   a_nicet,   a_micep,   a_micet,                                &
       a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,                               &
       snowin,    prtcl, calc_eff_radius
  use thrm, only : thermo
  use stat, only : sflg, updtst, acc_removal, mcflg, acc_massbudged, cs_rem_set
  USE mo_submctl, ONLY : terminal_vel
  implicit none

  logical, parameter :: droplet_sedim = .False., khairoutdinov = .False.

  LOGICAL :: sed_aero = .TRUE.,  &
             sed_cloud = .TRUE., &
             sed_precp = .TRUE., &
             sed_ice = .TRUE., &
             sed_snow = .TRUE.
  !
  ! drop sizes definition is based on vanZanten (2005)
  ! cloud droplets' diameter: 2-50 e-6 m
  ! drizzle drops' diameter: 50-1000 e-6 m
  !
  real, parameter :: eps0 = 1e-20       ! small number
  real, parameter :: eps1 = 1e-9        ! small number
  real, parameter :: rho_0 = 1.21       ! air density at surface

  real, parameter :: D_min = 2.e-6      ! minimum diameter of cloud droplets
  real, parameter :: D_bnd = 80.e-6     ! precip/cloud diameter boundary
  real, parameter :: D_max = 1.e-3      ! maximum diameter of prcp drops

  real, parameter :: X_min = (D_min**3)*rowt*pi/6. ! min cld mass
  real, parameter :: X_bnd = (D_bnd**3)*rowt*pi/6. ! precip/cld bound mass
  real, parameter :: X_max = (D_max**3)*rowt*pi/6. ! max prcp mass

  real, parameter :: prw = pi * rowt / 6.

contains

  !
  ! ---------------------------------------------------------------------
  ! MICRO: sets up call to microphysics
  !
  subroutine micro(level)
    USE class_componentIndex, ONLY : GetNcomp
    integer, intent (in) :: level
    INTEGER :: nn

    select case (level)
    case(2)
       if (droplet_sedim)  &
            call sedim_cd(nzp,nxp,nyp,a_theta,a_temp,a_rc,precip,a_rt,a_tt)
    case(3)
       call mcrph(nzp,nxp,nyp,dn0,a_theta,a_temp,a_rv,a_rsl,a_rc,a_rpp,   &
                  a_npp,precip,a_rt,a_tt,a_rpt,a_npt)
    case(4,5)
       IF (level < 5) THEN
            sed_ice = .FALSE.; sed_snow = .FALSE.
       ENDIF
       nn = GetNcomp(prtcl)+1
       CALL sedim_SALSA(nzp,nxp,nyp,nn, dtlt, a_temp, a_theta,               &
                        a_naerop,  a_naerot,  a_maerop,  a_maerot,           &
                        a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,          &
                        a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,          &
                        a_nicep,   a_nicet,   a_micep,   a_micet,            &
                        a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,           &
                        a_ustar, precip, snowin, a_tt                )
    end select

  end subroutine micro

  !
  ! ---------------------------------------------------------------------
  ! MCRPH: calls microphysical parameterization
  !
  subroutine mcrph(n1,n2,n3,dn0,th,tk,rv,rs,rc,rp,np,rrate,         &
       rtt,tlt,rpt,npt)

    integer, intent (in) :: n1,n2,n3
    real, dimension(n1,n2,n3), intent (in)    :: th, tk, rv, rs
    real, dimension(n1)      , intent (in)    :: dn0
    real, dimension(n1,n2,n3), intent (inout) :: rc, rtt, tlt, rpt, npt, np, rp
    real, intent (out)                        :: rrate(n1,n2,n3)

    integer :: i, j, k

    !
    ! Microphysics following Seifert Beheng (2001, 2005)
    ! note that the order below is important as the rc array is
    ! redefined in cld_dgn below and is assumed to be cloud water
    ! after that and total condensate priort to that
    !

    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             rp(k,i,j) = max(0., rp(k,i,j))
             np(k,i,j) = max(min(rp(k,i,j)/X_bnd,np(k,i,j)),rp(k,i,j)/X_max)
          end do
       end do
    end do

    call wtr_dff_SB(n1,n2,n3,dn0,rp,np,rc,rs,rv,tk,rpt,npt)

    call auto_SB(n1,n2,n3,dn0,rc,rp,rpt,npt)

    call accr_SB(n1,n2,n3,dn0,rc,rp,np,rpt,npt)

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             rp(k,i,j) = rp(k,i,j) + max(-rp(k,i,j)/dtlt,rpt(k,i,j))*dtlt
             np(k,i,j) = np(k,i,j) + max(-np(k,i,j)/dtlt,npt(k,i,j))*dtlt
             rpt(k,i,j)= 0.
             npt(k,i,j)= 0.
             rp(k,i,j) = max(0., rp(k,i,j))
             np(k,i,j) = max(min(rp(k,i,j)/X_bnd,np(k,i,j)),rp(k,i,j)/X_max)
          end do
       end do
    end do

    call sedim_rd(n1,n2,n3,dtlt,dn0,rp,np,tk,th,rrate,rtt,tlt,rpt,npt)

    if (droplet_sedim) call sedim_cd(n1,n2,n3,th,tk,rc,rrate,rtt,tlt)

  end subroutine mcrph
  !
  ! ---------------------------------------------------------------------
  ! WTR_DFF_SB: calculates the evolution of the both number- and
  ! mass mixing ratio large drops due to evaporation in the absence of
  ! cloud water.
  !
  subroutine wtr_dff_SB(n1,n2,n3,dn0,rp,np,rl,rs,rv,tk,rpt,npt)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: tk(n1,n2,n3), np(n1,n2,n3), rp(n1,n2,n3),        &
         rs(n1,n2,n3),rv(n1,n2,n3), dn0(n1)
    real, intent (inout) :: rpt(n1,n2,n3), npt(n1,n2,n3), rl(n1,n2,n3)

    real, parameter    :: Kt = 2.5e-2    ! conductivity of heat [J/(sKm)]
    real, parameter    :: Dv = 3.e-5     ! diffusivity of water vapor [m2/s]

    integer             :: i, j, k
    real                :: Xp, Dp, G, S, cerpt, cenpt, xnpts
    real, dimension(n1) :: v1

    if(sflg) then
       xnpts = 1./((n3-4)*(n2-4))
       do k=1,n1
          v1(k) = 0.
       end do
    end if

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1
             if (rp(k,i,j) > rl(k,i,j)) then
                Xp = rp(k,i,j)/ (np(k,i,j)+eps0)
                Xp = MIN(MAX(Xp,X_bnd),X_max)
                Dp = ( Xp / prw )**(1./3.)

                G = 1. / (1. / (dn0(k)*rs(k,i,j)*Dv) + &
                     alvl*(alvl/(Rm*tk(k,i,j))-1.) / (Kt*tk(k,i,j)))
                S = rv(k,i,j)/rs(k,i,j) - 1.

                if (S < 0) then
                   cerpt = 2. * pi * Dp * G * S * np(k,i,j)
                   cenpt = cerpt * np(k,i,j) / rp(k,i,j)
                   rpt(k,i,j)=rpt(k,i,j) + cerpt
                   npt(k,i,j)=npt(k,i,j) + cenpt
                   if (sflg) v1(k) = v1(k) + cerpt * xnpts
                end if
             end if
             rl(k,i,j) = max(0.,rl(k,i,j) - rp(k,i,j))
          end do
       end do
    end do

    if (sflg) call updtst(n1,'prc',2,v1,1)

  end subroutine wtr_dff_SB
  !
  ! ---------------------------------------------------------------------
  ! AUTO_SB:  calculates the evolution of mass- and number mxg-ratio for
  ! drizzle drops due to autoconversion. The autoconversion rate assumes
  ! f(x)=A*x**(nu_c)*exp(-Bx), an exponential in drop MASS x. It can
  ! be reformulated for f(x)=A*x**(nu_c)*exp(-Bx**(mu)), where formu=1/3
  ! one would get a gamma dist in drop diam -> faster rain formation.
  !
  subroutine auto_SB(n1,n2,n3,dn0,rc,rp,rpt,npt)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dn0(n1), rc(n1,n2,n3), rp(n1,n2,n3)
    real, intent (inout) :: rpt(n1,n2,n3), npt(n1,n2,n3)

    real, parameter :: nu_c  = 0           ! width parameter of cloud DSD
    real, parameter :: k_c  = 9.44e+9      ! Long-Kernel
    real, parameter :: k_1  = 6.e+2        ! Parameter for phi function
    real, parameter :: k_2  = 0.68         ! Parameter for phi function

    real, parameter :: Cau = 4.1e-15 ! autoconv. coefficient in KK param.
    real, parameter :: Eau = 5.67    ! autoconv. exponent in KK param.
    real, parameter :: mmt = 1.e+6   ! transformation from m to \mu m

    integer :: i, j, k
    real    :: k_au, Xc, Dc, au, tau, phi

    k_au  = k_c / (20.*X_bnd) * (nu_c+2.)*(nu_c+4.)/(nu_c+1.)**2

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             Xc = rc(k,i,j)/(CCN+eps0)
             if (Xc > 0.) then
                Xc = MIN(MAX(Xc,X_min),X_bnd)
                au = k_au * dn0(k) * rc(k,i,j)**2 * Xc**2
                !
                ! small threshold that should not influence the result
                !
                if (rc(k,i,j) > 1.e-6) then
                   tau = 1.0-rc(k,i,j)/(rc(k,i,j)+rp(k,i,j)+eps0)
                   tau = MIN(MAX(tau,eps0),0.9)
                   phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
                   au  = au * (1.0 + phi/(1.0 - tau)**2)
                endif
                !
                ! Khairoutdinov and Kogan
                !
                if (khairoutdinov) then
                   Dc = ( Xc / prw )**(1./3.)
                   au = Cau * (Dc * mmt / 2.)**Eau
                end if

                rpt(k,i,j) = rpt(k,i,j) + au
                npt(k,i,j) = npt(k,i,j) + au/X_bnd
                !
             end if
          end do
       end do
    end do

  end subroutine auto_SB
  !
  ! ---------------------------------------------------------------------
  ! ACCR_SB calculates the evolution of mass mxng-ratio due to accretion
  ! and self collection following Seifert & Beheng (2001).  Included is
  ! an alternative formulation for accretion only, following
  ! Khairoutdinov and Kogan
  !
  subroutine accr_SB(n1,n2,n3,dn0,rc,rp,np,rpt,npt)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: rc(n1,n2,n3), rp(n1,n2,n3), np(n1,n2,n3), dn0(n1)
    real, intent (inout) :: rpt(n1,n2,n3),npt(n1,n2,n3)

    real, parameter :: k_r = 5.78
    real, parameter :: k_1 = 5.e-4
    real, parameter :: Cac = 67.     ! accretion coefficient in KK param.
    real, parameter :: Eac = 1.15    ! accretion exponent in KK param.

    integer :: i, j, k
    real    :: tau, phi, ac, sc

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             if (rc(k,i,j) > 0. .and. rp(k,i,j) > 0.) then
                tau = 1.0-rc(k,i,j)/(rc(k,i,j)+rp(k,i,j)+eps0)
                tau = MIN(MAX(tau,eps0),1.)
                phi = (tau/(tau+k_1))**4
                ac  = k_r * rc(k,i,j) * rp(k,i,j) * phi * sqrt(rho_0*dn0(k))
                !
                ! Khairoutdinov and Kogan
                !
                !ac = Cac * (rc(k,i,j) * rp(k,i,j))**Eac
                !
                rpt(k,i,j) = rpt(k,i,j) + ac

             end if
             sc = k_r * np(k,i,j) * rp(k,i,j) * sqrt(rho_0*dn0(k))
             npt(k,i,j) = npt(k,i,j) - sc
          end do
       end do
    end do

  end subroutine accr_SB
   !
   ! ---------------------------------------------------------------------
   ! SEDIM_RD: calculates the sedimentation of the rain drops and its
   ! effect on the evolution of theta_l and r_t.  This is expressed in
   ! terms of Dp the mean diameter, not the mass weighted mean diameter
   ! as is used elsewhere.  This is just 1/lambda in the exponential
   ! distribution
   !
   subroutine sedim_rd(n1,n2,n3,dt,dn0,rp,np,tk,th,rrate,rtt,tlt,rpt,npt)

     integer, intent (in)                      :: n1,n2,n3
     real, intent (in)                         :: dt
     real, intent (in),    dimension(n1)       :: dn0
     real, intent (in),    dimension(n1,n2,n3) :: rp, np, th, tk
     real, intent (out),   dimension(n1,n2,n3) :: rrate
     real, intent (inout), dimension(n1,n2,n3) :: rtt, tlt, rpt, npt

     real, parameter :: a2 = 9.65       ! in SI [m/s]
     real, parameter :: c2 = 6e2        ! in SI [1/m]
     real, parameter :: Dv = 25.0e-6    ! in SI [m/s]
     real, parameter :: cmur1 = 10.0    ! mu-Dm-relation for rain following
     real, parameter :: cmur2 = 1.20e+3 ! Milbrandt&Yau 2005, JAS, but with
     real, parameter :: cmur3 = 1.5e-3  ! revised constants
     real, parameter :: aq = 6.0e3
     real, parameter :: bq = -0.2
     real, parameter :: an = 3.5e3
     real, parameter :: bn = -0.1


     integer :: i, j, k, kp1, kk, km1
     real    :: b2, Xp, Dp, Dm, mu, flxdiv, tot,sk, mini, maxi, cc, zz
     real, dimension(n1) :: nslope,rslope,dn,dr, rfl, nfl, vn, vr, cn, cr

     b2 = a2*exp(c2*Dv)

     do j=3,n3-2
        do i=3,n2-2

           nfl(n1) = 0.
           rfl(n1) = 0.
           do k=n1-1,2,-1
              Xp = rp(k,i,j) / (np(k,i,j)+eps0)
              Xp = MIN(MAX(Xp,X_bnd),X_max)
              !
              ! Adjust Dm and mu-Dm and Dp=1/lambda following Milbrandt & Yau
              !
              Dm = ( 6. / (rowt*pi) * Xp )**(1./3.)
              mu = cmur1*(1.+tanh(cmur2*(Dm-cmur3)))
              Dp = (Dm**3/((mu+3.)*(mu+2.)*(mu+1.)))**(1./3.)

              vn(k) = sqrt(dn0(k)/1.2)*(a2 - b2*(1.+c2*Dp)**(-(1.+mu)))
              vr(k) = sqrt(dn0(k)/1.2)*(a2 - b2*(1.+c2*Dp)**(-(4.+mu)))
              !
              ! Set fall speeds following Khairoutdinov and Kogan

              if (khairoutdinov) then
                 vn(k) = max(0.,an * Dp + bn)
                 vr(k) = max(0.,aq * Dp + bq)
              end if

           end do

           do k=2,n1-1
              kp1 = min(k+1,n1-1)
              km1 = max(k,2)
              cn(k) = 0.25*(vn(kp1)+2.*vn(k)+vn(km1))*dzt(k)*dt
              cr(k) = 0.25*(vr(kp1)+2.*vr(k)+vr(km1))*dzt(k)*dt
           end do

           !...piecewise linear method: get slopes
           do k=n1-1,2,-1
              dn(k) = np(k+1,i,j)-np(k,i,j)
              dr(k) = rp(k+1,i,j)-rp(k,i,j)
           enddo
           dn(1)  = dn(2)
           dn(n1) = dn(n1-1)
           dr(1)  = dr(2)
           dr(n1) = dr(n1-1)
           do k=n1-1,2,-1
              !...slope with monotone limiter for np
              sk = 0.5 * (dn(k-1) + dn(k))
              mini = min(np(k-1,i,j),np(k,i,j),np(k+1,i,j))
              maxi = max(np(k-1,i,j),np(k,i,j),np(k+1,i,j))
              nslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(np(k,i,j)-mini), &
                   &                                     2.*(maxi-np(k,i,j)))
              !...slope with monotone limiter for rp
              sk = 0.5 * (dr(k-1) + dr(k))
              mini = min(rp(k-1,i,j),rp(k,i,j),rp(k+1,i,j))
              maxi = max(rp(k-1,i,j),rp(k,i,j),rp(k+1,i,j))
              rslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(rp(k,i,j)-mini), &
                   &                                     2.*(maxi-rp(k,i,j)))
           enddo

           rfl(n1-1) = 0.
           nfl(n1-1) = 0.
           do k=n1-2,2,-1

              kk = k
              tot = 0.0
              zz  = 0.0
              cc  = min(1.,cn(k))
              do while (cc > 0 .and. kk <= n1-1)
                 tot = tot + dn0(kk)*(np(kk,i,j)+nslope(kk)*(1.-cc))*cc/dzt(kk)
                 zz  = zz + 1./dzt(kk)
                 kk  = kk + 1
                 cc  = min(1.,cn(kk) - zz*dzt(kk))
              enddo
              nfl(k) = -tot /dt

              kk = k
              tot = 0.0
              zz  = 0.0
              cc  = min(1.,cr(k))
              do while (cc > 0 .and. kk <= n1-1)
                 tot = tot + dn0(kk)*(rp(kk,i,j)+rslope(kk)*(1.-cc))*cc/dzt(kk)
                 zz  = zz + 1./dzt(kk)
                 kk  = kk + 1
                 cc  = min(1.,cr(kk) - zz*dzt(kk))
              enddo
              rfl(k) = -tot /dt

              kp1=k+1
              flxdiv = (rfl(kp1)-rfl(k))*dzt(k)/dn0(k)
              rpt(k,i,j) =rpt(k,i,j)-flxdiv
              rtt(k,i,j) =rtt(k,i,j)-flxdiv
              tlt(k,i,j) =tlt(k,i,j)+flxdiv*(alvl/cp)*th(k,i,j)/tk(k,i,j)

              npt(k,i,j) = npt(k,i,j)-(nfl(kp1)-nfl(k))*dzt(k)/dn0(k)

              rrate(k,i,j)    = -rfl(k)/dn0(k) * alvl*0.5*(dn0(k)+dn0(kp1))

           end do
        end do
     end do

   end subroutine sedim_rd
  !
  ! ---------------------------------------------------------------------
  ! SEDIM_CD: calculates the cloud-droplet sedimentation flux and its effect
  ! on the evolution of r_t and theta_l assuming a log-normal distribution
  !
  subroutine sedim_cd(n1,n2,n3,th,tk,rc,rrate,rtt,tlt)


    integer, intent (in):: n1,n2,n3
    real, intent (in),   dimension(n1,n2,n3) :: th,tk,rc
    real, intent (out),  dimension(n1,n2,n3) :: rrate
    real, intent (inout),dimension(n1,n2,n3) :: rtt,tlt

    real, parameter :: c = 1.19e8 ! Stokes fall velocity coef [m^-1 s^-1]
    real, parameter :: sgg = 1.2  ! geometric standard dev of cloud droplets

    integer :: i, j, k, kp1
    real    :: Dc, Xc, vc, flxdiv
    real    :: rfl(n1)

    !
    ! calculate the precipitation flux and its effect on r_t and theta_l
    !
    do j=3,n3-2
       do i=3,n2-2
          rfl(n1) = 0.
          do k=n1-1,2,-1
             Xc = rc(k,i,j) / (CCN+eps0)
             Dc = ( Xc / prw )**(1./3.)
             Dc = MIN(MAX(Dc,D_min),D_bnd)
             vc = min(c*(Dc*0.5)**2 * exp(4.5*(log(sgg))**2),1./(dzt(k)*dtlt))
             rfl(k) = - rc(k,i,j) * vc
             !
             kp1=k+1
             flxdiv = (rfl(kp1)-rfl(k))*dzt(k)
             rtt(k,i,j) = rtt(k,i,j)-flxdiv
             tlt(k,i,j) = tlt(k,i,j)+flxdiv*(alvl/cp)*th(k,i,j)/tk(k,i,j)
             rrate(k,i,j) = -rfl(k)
          end do
       end do
    end do

  end subroutine sedim_cd



  ! ---------------------------------------------------------------------
  ! SEDIM_AERO: calculates the salsa particles sedimentation and dry deposition flux  (.. Zubair) !
  !
  ! Juha: The code below is a modified version of the original one by Zubair
  !
  ! Juha: Rain is now treated completely separately (20151013)
  !
  ! Jaakko: Modified for the use of ice and snow bins

  SUBROUTINE sedim_SALSA(n1,n2,n3,n4,tstep,tk,th,          &
                         naerop, naerot, maerop, maerot,   &
                         ncloudp,ncloudt,mcloudp,mcloudt,  &
                         nprecpp,nprecpt,mprecpp,mprecpt,  &
                         nicep,  nicet,  micep,  micet,    &
                         nsnowp, nsnowt, msnowp, msnowt,   &
                         ustar,rrate, srate, tlt          )

    USE mo_submctl, ONLY : nbins, ncld, nprc,           &
                               nice,  nsnw,                 &
                               nlim,prlim,  &
                               rhowa,rhoic
    USE class_ComponentIndex, ONLY : GetIndex, GetNcomp
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3,n4
    REAL, INTENT(in) :: tstep,                    &
                        tk(n1,n2,n3),             &
                        th(n1,n2,n3),             &
                        ustar(n2,n3),             &
                        naerop(n1,n2,n3,nbins),   &
                        maerop(n1,n2,n3,n4*nbins), &
                        ncloudp(n1,n2,n3,ncld),    &
                        mcloudp(n1,n2,n3,n4*ncld), &
                        nprecpp(n1,n2,n3,nprc),    &
                        mprecpp(n1,n2,n3,n4*nprc), &
                        nicep(n1,n2,n3,nice),    &
                        micep(n1,n2,n3,n4*nice), &
                        nsnowp(n1,n2,n3,nsnw),    &
                        msnowp(n1,n2,n3,n4*nsnw)

    REAL, INTENT(inout) :: naerot(n1,n2,n3,nbins),    &
                           maerot(n1,n2,n3,n4*nbins), &
                           ncloudt(n1,n2,n3,ncld),    &
                           mcloudt(n1,n2,n3,n4*ncld), &
                           nprecpt(n1,n2,n3,nprc),    &
                           mprecpt(n1,n2,n3,n4*nprc), &
                           nicet(n1,n2,n3,nice),      &
                           micet(n1,n2,n3,n4*nice),   &
                           nsnowt(n1,n2,n3,nsnw),    &
                           msnowt(n1,n2,n3,n4*nsnw), &
                           tlt(n1,n2,n3),             & ! Liquid water pot temp tendency
                           rrate(n1,n2,n3),           &
                           srate(n1,n2,n3)              ! snowing rate

    INTEGER :: ss
    INTEGER :: i,j,k,nc,istr,iend

    REAL :: prnt(n1,n2,n3,nprc), prvt(n1,n2,n3,n4*nprc)  ! Rain number and mass tendencies due to fallout
    REAL :: srnt(n1,n2,n3,nsnw), srvt(n1,n2,n3,n4*nsnw)  ! Snow number and mass tendencies due to fallout

    ! Particle mass removal arrays, given in kg/(m2 s)
    REAL :: remaer(n2,n3,n4*nbins),   &
            remcld(n2,n3,n4*ncld),    &
            remprc(n2,n3,n4*nprc),    &
            remice(n2,n3,n4*nice),    &
            remsnw(n2,n3,n4*nsnw)

    ! Particle number removal arrays
    REAL :: andep(n2,n3,nbins),     &
            cndep(n2,n3,ncld),      &
            indep(n2,n3,nice)

    REAL :: mctmp(n2,n3) ! Helper for mass conservation calculations

    ! Divergence fields
    REAL :: amdiv(n1,n2,n3,n4*nbins),    &
            cmdiv(n1,n2,n3,n4*ncld),     &
            imdiv(n1,n2,n3,n4*nice)
    REAL :: andiv(n1,n2,n3,nbins),       &
            cndiv(n1,n2,n3,ncld),        &
            indiv(n1,n2,n3,nice)


    remaer = 0.; remcld = 0.; remprc = 0.; remice = 0.; remsnw = 0.

    ! Sedimentation for slow (non-precipitating) particles
    !-------------------------------------------------------
    IF (sed_aero) THEN

       CALL DepositionSlow(n1,n2,n3,n4,nbins,tk,a_dn,1500.,ustar,naerop,maerop,dzt,nlim,andiv,amdiv,andep,remaer,1)

       naerot = naerot - andiv
       maerot = maerot - amdiv

       ! Account for changes in liquid water pot temperature
       nc = GetIndex(prtcl,'H2O')
       istr = (nc-1)*nbins+1
       iend = nc*nbins
       DO j = 3,n3-2
          DO i = 3,n2-2
             DO k = 2,n1
                tlt(k,i,j) = tlt(k,i,j) + SUM(amdiv(k,i,j,istr:iend))*(alvl/cp)*th(k,i,j)/tk(k,i,j)
             END DO
          END DO
       END DO

    END IF ! sed_aero

    IF (sed_cloud) THEN

       CALL DepositionSlow(n1,n2,n3,n4,ncld,tk,a_dn,rhowa,ustar,ncloudp,mcloudp,dzt,nlim,cndiv,cmdiv,cndep,remcld,2)

       ncloudt = ncloudt - cndiv
       mcloudt = mcloudt - cmdiv

       ! Account for changes in liquid water pot temperature
       nc = GetIndex(prtcl,'H2O')
       istr = (nc-1)*ncld+1
       iend = nc*ncld
       DO j = 3,n3-2
          DO i = 3,n2-2
             DO k = 2,n1
                tlt(k,i,j) = tlt(k,i,j) + SUM(cmdiv(k,i,j,istr:iend))*(alvl/cp)*th(k,i,j)/tk(k,i,j)
             END DO
          END DO
       END DO

    END IF ! sed_cloud

    IF (sed_ice) THEN

       CALL DepositionSlow(n1,n2,n3,n4,nice,tk,a_dn,rhoic,ustar,nicep,micep,dzt,prlim,indiv,imdiv,indep,remice,4)

       nicet = nicet - indiv 
       micet = micet - imdiv 

       ! Account for changes in liquid water pot temperature
       nc = GetIndex(prtcl,'H2O')
       istr = (nc-1)*nice+1
       iend = nc*nice
       DO j = 3,n3-2
          DO i = 3,n2-2
             DO k = 2,n1
                tlt(k,i,j) = tlt(k,i,j) + SUM(imdiv(k,i,j,istr:iend))*(alvi/cp)*th(k,i,j)/tk(k,i,j)
             END DO
          END DO
       END DO

    END IF ! sed_ice


    ! ---------------------------------------------------------
    ! SEDIMENTATION/DEPOSITION OF FAST PRECIPITATING PARTICLES
    IF (sed_precp) THEN
       CALL DepositionFast(n1,n2,n3,n4,nprc,tk,a_dn,rowt,nprecpp,mprecpp,tstep,dzt,prnt,prvt,remprc,prlim,rrate,3)

       nprecpt(:,:,:,:) = nprecpt(:,:,:,:) + prnt(:,:,:,:)/tstep
       mprecpt(:,:,:,:) = mprecpt(:,:,:,:) + prvt(:,:,:,:)/tstep

       ! Convert mass flux to heat flux (W/m^2)
       rrate(:,:,:)=rrate(:,:,:)*alvl

       nc = GetIndex(prtcl,'H2O')
       istr = (nc-1)*nprc + 1
       iend = nc*nprc
       DO j = 3,n3-2
          DO i = 3,n2-2
             DO k = 1,n1-1
                tlt(k,i,j) = tlt(k,i,j) - SUM(prvt(k,i,j,istr:iend))/tstep*(alvl/cp)*th(k,i,j)/tk(k,i,j)
             END DO
          END DO
       END DO
    END IF

    IF (sed_snow) THEN
       CALL DepositionFast(n1,n2,n3,n4,nsnw,tk,a_dn,rowt,nsnowp,msnowp,tstep,dzt,srnt,srvt,remsnw,prlim,srate,5)

       nsnowt(:,:,:,:) = nsnowt(:,:,:,:) + srnt(:,:,:,:)/tstep
       msnowt(:,:,:,:) = msnowt(:,:,:,:) + srvt(:,:,:,:)/tstep

       ! Convert mass flux to heat flux (W/m^2)
       srate(:,:,:)=srate(:,:,:)*alvi

       nc = GetIndex(prtcl,'H2O')
       istr = (nc-1)*nsnw + 1
       iend = nc*nsnw
       DO j = 3,n3-2
          DO i = 3,n2-2
             DO k = 1,n1-1
                tlt(k,i,j) = tlt(k,i,j) - SUM(srvt(k,i,j,istr:iend))/tstep*(alvi/cp)*th(k,i,j)/tk(k,i,j)
             END DO
          END DO
       END DO
    END IF

    IF (mcflg) THEN
       ! For mass conservation statistics
       mctmp(:,:) = 0.
       ss = getIndex(prtcl,'H2O')
       istr = (ss-1)*nbins; iend = ss*nbins
       mctmp(:,:) = mctmp(:,:) + SUM(remaer(:,:,istr:iend),dim=3)
       istr = (ss-1)*ncld; iend = ss*ncld
       mctmp(:,:) = mctmp(:,:) + SUM(remcld(:,:,istr:iend),dim=3)
       istr = (ss-1)*nprc; iend = ss*nprc
       mctmp(:,:) = mctmp(:,:) + SUM(remprc(:,:,istr:iend),dim=3)
       istr = (ss-1)*nice; iend = ss*nice
       mctmp(:,:) = mctmp(:,:) + SUM(remice(:,:,istr:iend),dim=3)
       istr = (ss-1)*nsnw; iend = ss*nsnw
       mctmp(:,:) = mctmp(:,:) + SUM(remsnw(:,:,istr:iend),dim=3)
       CALL acc_massbudged(n1,n2,n3,3,tstep,dzt,a_dn,rdep=mctmp,ApVdom=mc_ApVdom)
    END IF !mcflg
    ! Aerosol removal statistics
    IF (sflg) CALL acc_removal(n2,n3,n4,remaer,remcld,remprc,remice,remsnw)
    IF (sflg) CALL cs_rem_set(n2,n3,n4,remaer,remcld,remprc,remice,remsnw)

  END SUBROUTINE !sedim_SALSA


  ! -----------------------------------------------------------------



  SUBROUTINE DepositionSlow(n1,n2,n3,n4,nn,tk,adn,pdn,ustar,numc,mass,dzt,clim,flxdivn,flxdivm,depflxn,depflxm,flag)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3,n4       ! Grid numbers, number of chemical species
    INTEGER, INTENT(in) :: nn                ! Number of bins
    REAL, INTENT(in) :: tk(n1,n2,n3)         ! Absolute temprature
    REAL, INTENT(in) :: adn(n1,n2,n3)        ! Air density
    REAL, INTENT(in) :: pdn                  ! Particle density
    REAL, INTENT(in) :: ustar(n2,n3)         ! Friction velocity
    REAL, INTENT(in) :: numc(n1,n2,n3,nn)    ! Particle number concentration
    REAL, INTENT(in) :: mass(n1,n2,n3,nn*n4) ! Particle mass mixing ratio
    REAL, INTENT(in) :: dzt(n1)              ! Inverse of grid level thickness
    REAL, INTENT(IN) :: clim                ! Concentration limit (#/m^3)
    INTEGER, INTENT(IN) :: flag         ! An option for identifying aerosol, cloud, precipitation, ice and snow
    REAL, INTENT(OUT) :: flxdivm(n1,n2,n3,nn*n4), flxdivn(n1,n2,n3,nn) ! Mass and number divergence
    REAL, INTENT(OUT) :: depflxn(n2,n3,nn), depflxm(n2,n3,nn*n4) ! Mass and number deposition fluxes to the surface

    INTEGER :: i,j,k,kp1
    INTEGER :: bin,bs

    real, parameter :: M = 4.8096e-26 ! average mass of one air molecule, eq2.3 fundamentals of atm.
                                      ! modelling [kg molec-1]
    real, parameter :: A = 1.249      ! fundamentals of atm. modelling pg509
    real, parameter :: B = 0.42
    real, parameter :: C = 0.87

    REAL :: avis, kvis           ! Viscosity of air, kinematic viscosity
    REAL :: lambda              ! Mean free path
    REAL :: va                    ! Thermal speed of air molecule
    REAL :: Kn, GG               ! Knudsen number, slip correction
    REAL :: vc                    ! Critical fall speed
    REAL :: mdiff                ! Particle diffusivity
    REAL :: rt, Sc, St

    REAL :: rflm(n1,nn*n4), rfln(n1,nn), prvolc(n4), rwet
    flxdivm = 0.
    flxdivn = 0.
    depflxm = 0.
    depflxn = 0.

    DO j = 3,n3-2
       DO i = 3,n2-2

          rflm = 0.
          rfln = 0.
          
          DO k=n1-1,2,-1
             kp1 = k+1

             ! atm modelling Eq.4.54
             avis=1.8325e-5*(416.16/(tk(k,i,j)+120.0))*(tk(k,i,j)/296.16)**1.5
             kvis =  avis/adn(k,i,j)
             va = sqrt(8*kb*tk(k,i,j)/(pi*M)) ! thermal speed of air molecule
             lambda = 2*avis/(adn(k,i,j)*va) !mean free path

             ! Fluxes
             !------------------
             ! -- Calculate the *corrections* for small particles
             DO bin = 1,nn
                IF (numc(k,i,j,bin)*adn(k,i,j)<clim) CYCLE

                ! Calculate wet size
                !   n4 = number of active species
                !   bin = size bin
                prvolc(:)=mass(k,i,j,bin:(n4-1)*nn+bin:nn)
                rwet=calc_eff_radius(n4,numc(k,i,j,bin),prvolc,flag)

                ! Terminal velocity
                Kn = lambda/rwet
                GG = 1.+ Kn*(A+B*exp(-C/Kn))
                vc = terminal_vel(rwet,pdn,adn(k,i,j),avis,GG,flag)

                IF (k==2) THEN ! The level just above surface
                    ! Particle diffusitivity  (15.29) in jacobson book
                    mdiff = (kb*tk(k,i,j)*GG)/(6.0*pi*avis*rwet)
                    Sc = kvis/mdiff
                    St = vc*ustar(i,j)**2.0/g*kvis
                    if (St<0.01) St=0.01
                    rt = 1.0/MAX(epsilon(1.0),(ustar(i,j)*(Sc**(-2.0/3.0)+10**(-3.0/St)))) ! atm chem&phy eq19.18
                    vc = (1./rt) + vc
                ENDIF

                ! Flux for the particle mass
                DO bs = bin, (n4-1)*nn + bin, nn
                    rflm(k,bs) = -mass(k,i,j,bs)*vc
                END DO

                ! Flux for the particle number
                 rfln(k,bin) = -numc(k,i,j,bin)*vc
             END DO

             flxdivm(k,i,j,:) = (rflm(kp1,:)-rflm(k,:))*dzt(k)
             flxdivn(k,i,j,:) = (rfln(kp1,:)-rfln(k,:))*dzt(k)

          END DO ! k

          ! Deposition flux to surface
          k=2
          depflxm(i,j,:) = -rflm(k,:)*dzt(k)
          depflxn(i,j,:) = -rfln(k,:)*dzt(k)

       END DO ! i
    END DO ! j

  END SUBROUTINE DepositionSlow


  !------------------------------------------------------------------
  SUBROUTINE DepositionFast(n1,n2,n3,n4,nn,tk,adn,pdn,numc,mass,tstep,dzt,prnt,prvt,remprc,clim,rate,flag)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3,n4,nn
    REAL, INTENT(in) :: tk(n1,n2,n3)
    REAL, INTENT(in) :: adn(n1,n2,n3)
    REAL, INTENT(in) :: pdn
    REAL, INTENT(in) :: numc(n1,n2,n3,nn)
    REAL, INTENT(in) :: mass(n1,n2,n3,nn*n4)
    REAL, INTENT(in) :: tstep
    REAL, INTENT(in) :: dzt(n1)
    REAL, INTENT(IN) :: clim                ! Concentration limit (#/m^3)
    INTEGER, INTENT(IN) :: flag         ! An option for identifying aerosol, cloud, precipitation, ice and snow
    REAL, INTENT(out) :: prnt(n1,n2,n3,nn), prvt(n1,n2,n3,nn*n4)     ! Number and mass tendencies due to fallout
    REAL, INTENT(out) :: remprc(n2,n3,nn*n4)
    REAL, INTENT(out) :: rate(n1,n2,n3) ! Rain rate (kg/s/m^2)

    INTEGER :: k,i,j,bin
    INTEGER :: istr,iend

    real, parameter :: A = 1.249 ! fundamentals of atm. modelling pg509
    real, parameter :: B = 0.42
    real, parameter :: C = 0.87
    real, parameter :: M = 4.8096e-26 ! average mass of one air molecule, eq2.3 fundamentals of atm.
                                      ! modelling [kg molec-1]

    REAL :: GG
    REAL :: Kn
    REAL :: lambda      ! Mean free path
    REAL :: avis,kvis   ! Air viscosity, kinematic viscosity
    REAL :: va          ! Thermal speed of air molecule
    REAL :: vc

    ! For precipitation:
    REAL :: fd,fdmax,fdos ! Fall distance for rain drops, max fall distance, overshoot from nearest grid level
    REAL :: prnchg(n1,nn), prvchg(n1,nn,n4) ! Instantaneous changes in precipitation number and mass (volume)
    REAL :: rwet
 
    REAL :: prnumc, prvolc(n4)  ! Instantaneous source number and mass
    INTEGER :: kf, ni,fi
    LOGICAL :: prcdep  ! Deposition flag

    remprc(:,:,:) = 0.
    rate(:,:,:) = 0.
    prnt(:,:,:,:) = 0.
    prvt(:,:,:,:) = 0.

    DO j = 3,n3-2
       
       DO i = 3,n2-2

          prnchg = 0.
          prvchg = 0.
          
          DO k=n1-1,2,-1
          
             ! atm modelling Eq.4.54
             avis = 1.8325e-5*(416.16/(tk(k,i,j)+120.0))*(tk(k,i,j)/296.16)**1.5
             kvis = avis/adn(k,i,j) !actual density ???
             va = sqrt(8.*kb*tk(k,i,j)/(pi*M)) ! thermal speed of air molecule
             lambda = 2.*avis/(adn(k,i,j)*va) !mean free path

             ! Precipitation bin loop
             DO bin = 1,nn
                IF (numc(k,i,j,bin)*adn(k,i,j) < clim) CYCLE

                ! Calculate wet size
                !   n4 = number of active species
                !   bin = size bin
                prvolc(:)=mass(k,i,j,bin:(n4-1)*nn+bin:nn)
                rwet=calc_eff_radius(n4,numc(k,i,j,bin),prvolc,flag)

                ! Terminal velocity
                Kn = lambda/rwet
                GG = 1.+ Kn*(A+B*exp(-C/Kn))
                vc = terminal_vel(rwet,pdn,adn(k,i,j),avis,GG,flag)

                ! Rain rate statistics: removal of water from the current bin is accounted for
                ! Water is the last (n4) species and rain rate is given here kg/s/m^2
                rate(k,i,j)=rate(k,i,j)+mass(k,i,j,(n4-1)*nn+bin)*adn(k,i,j)*vc

                ! Determine output flux for current level: Find the closest level to which the
                ! current drop parcel can fall within 1 timestep. If the lowest atmospheric level
                ! is reached, the drops are sedimented.
                
                ! Maximum fall distance:
                fdmax = tstep*vc
             
                fd = 0.
                fi = 0
                prcdep = .FALSE. ! deposition flag
                DO WHILE ( fd < fdmax )
                   fd = fd + ( 1./dzt(k-fi) )
                   fi = fi + 1
                   ! Check if sedimentation occurs for current parcel
                   IF (k-fi <= 1) THEN
                      prcdep = .TRUE.
                      EXIT
                   END IF
                END DO
                ! Back up one level
                fi = fi-1
                kf = k - fi
                fd = fd - ( 1./dzt(kf) )
             
                ! How much the actual fall distance overshoots below the layer kf
                fdos = MIN(MAX(fdmax-fd,0.),1./dzt(kf))
             
                ! Remove the drops from the original level
                prnumc = numc(k,i,j,bin)
                prnchg(k,bin) = prnchg(k,bin) - prnumc
                DO ni = 1,n4
                   prvolc(ni) = mass(k,i,j,(ni-1)*nn+bin)
                   prvchg(k,bin,ni) = prvchg(k,bin,ni) - prvolc(ni)
                END DO ! ni

                ! Removal statistics
                IF (prcdep) THEN
                   DO ni=1,n4
                      remprc(i,j,(ni-1)*nn+bin) = remprc(i,j,(ni-1)*nn+bin) +    &
                           prvolc(ni)*adn(k,i,j)*vc
                   END DO
                ENDIF ! prcdep

                ! Put the drops to new positions (may be partially the original grid cell as well!)
                IF (fdos*dzt(kf) > 0.5) THEN  ! Reduce numerical diffusion
                   prnchg(kf-1,bin) = prnchg(kf-1,bin) + prnumc
                   DO ni = 1,n4
                      prvchg(kf-1,bin,ni) = prvchg(kf-1,bin,ni) + prvolc(ni)
                   END DO
                ELSE
                   prnchg(kf-1,bin) = prnchg(kf-1,bin) + ( fdos*dzt(kf) )*prnumc
                   prnchg(kf,bin) = prnchg(kf,bin) + ( 1. - fdos*dzt(kf) )*prnumc
                   DO ni = 1,n4
                      prvchg(kf-1,bin,ni) = prvchg(kf-1,bin,ni) + ( fdos*dzt(kf) )*prvolc(ni)
                      prvchg(kf,bin,ni) = prvchg(kf,bin,ni) + ( 1. - fdos*dzt(kf) )*prvolc(ni)
                   END DO
                END IF ! diffusion
             
             END DO !bin
          
          END DO ! k

          prnt(:,i,j,:) = prnt(:,i,j,:) + prnchg(:,:)
          DO ni = 1,n4
             istr = (ni-1)*nn+1
             iend = ni*nn
             prvt(:,i,j,istr:iend) = prvt(:,i,j,istr:iend) + prvchg(:,:,ni)
          END DO !ni

       END DO ! i

    END DO ! j

  END SUBROUTINE DepositionFast


end module mcrp
