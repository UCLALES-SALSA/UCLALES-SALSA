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
       a_Radry, a_Rawet, a_Rcdry, a_Rcwet, a_Rpdry, a_Rpwet,                    &
       a_Riwet, a_Rswet,                                                        &
       a_naerop,  a_naerot,  a_maerop,  a_maerot,                               &
       a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,                              &
       a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt, mc_ApVdom,         &
       a_nicep,   a_nicet,   a_micep,   a_micet,                                &
       a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,                               &
       snowin,    prtcl
  use thrm, only : thermo
  use stat, only : sflg, updtst, acc_removal, mcflg, acc_massbudged, cs_rem_set
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
       CALL sedim_SALSA(nzp,nxp,nyp,nn, level,dtlt, a_temp, a_theta,               &
                        a_Rawet,   a_Rcwet,   a_Rpwet,                       &
                        a_Riwet,   a_Rswet,                                  &
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
     real    :: b2, Xp, Dp, Dm, mu, flxdiv, tot,sk, mini, maxi, cc, zz, xnpts
     real, dimension(n1) :: nslope,rslope,dn,dr, rfl, nfl, vn, vr, cn, cr, v1

    if(sflg) then
       xnpts = 1./((n3-4)*(n2-4))
       do k=1,n1
          v1(k) = 0.
       end do
    end if

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

              rrate(k,i,j)    = -rfl(k) * alvl*0.5*(dn0(k)+dn0(kp1))
              if (sflg) v1(k) = v1(k) + rrate(k,i,j)*xnpts

           end do
        end do
     end do
     if (sflg) call updtst(n1,'prc',1,v1,1)

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

  SUBROUTINE sedim_SALSA(n1,n2,n3,n4,level,tstep,tk,th,          &
                         Rawet, Rcwet, Rpwet,              &
                         Riwet, Rswet,                     &
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
    INTEGER, INTENT(in) :: level 
    REAL, INTENT(in) :: tstep,                    &
                        tk(n1,n2,n3),             &
                        th(n1,n2,n3),             &

                        Rawet(n1,n2,n3,nbins),    &
                        Rcwet(n1,n2,n3,ncld),     &
                        Rpwet(n1,n2,n3,nprc),     &
                        Riwet(n1,n2,n3,nice),     &
                        Rswet(n1,n2,n3,nsnw),     &

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

    INTEGER :: ss, kp1
    INTEGER :: i,j,k,nc,istr,iend
    REAL :: xnpts
    REAL :: v1(n1)

    REAL :: prnt(n1,n2,n3,nprc), prvt(n1,n2,n3,n4*nprc)     ! Number and mass tendencies due to fallout

    REAL :: srnt(n1,n2,n3,nsnw), srvt(n1,n2,n3,n4*nsnw)     ! Number and mass tendencies due to fallout

    ! PArticle removal arrays, given in kg/(m2 s)
    REAL :: remaer(n2,n3,n4*nbins),   &
            remcld(n2,n3,n4*ncld),    &
            remprc(n2,n3,n4*nprc),    &
            remice(n2,n3,n4*nice),    &
            remsnw(n2,n3,n4*nsnw)

    REAL :: mctmp(n2,n3) ! Helper for mass conservation calculations

    ! Divergence fields
    REAL :: amdiv(n1,n2,n3,n4*nbins),    &
            cmdiv(n1,n2,n3,n4*ncld),     &
            imdiv(n1,n2,n3,n4*nice)
    REAL :: andiv(n1,n2,n3,nbins),       &
            cndiv(n1,n2,n3,ncld),        &
            indiv(n1,n2,n3,nice)
    
    ! Deposition fields (given as particle divergence)
    REAL :: amdep(n2,n3,n4*nbins),  &
            cmdep(n2,n3,n4*ncld),   &
            imdep(n2,n3,n4*nice)
    REAL :: andep(n2,n3,nbins),     &
            cndep(n2,n3,ncld),      &
            indep(n2,n3,nice)


    remaer = 0.; remcld = 0.; remprc = 0.; remice = 0.; remsnw = 0.
    amdiv = 0.; amdep = 0.; cmdiv = 0.; cmdep = 0.; imdiv = 0.; imdep = 0.
    andiv = 0.; andep = 0.; cndiv = 0.; cndep = 0.; indiv = 0.; indep = 0.
    prnt = 0.; prvt = 0.; srnt = 0.; srvt = 0.

    if(sflg) then
       xnpts = 1./((n3-4)*(n2-4))
       do k=1,n1
          v1(k) = 0.
       end do
    end if

    ! Sedimentation for slow (non-precipitating) particles
    !-------------------------------------------------------
    IF (sed_aero) THEN

       andiv = NumDivergence(n1,n2,n3,nbins,tk,a_dn,1500.,Rawet,naerop,dzt,nlim)
       amdiv = MassDivergence(n1,n2,n3,n4,nbins,tk,a_dn,1500.,Rawet,naerop,maerop,dzt,nlim)
       
       andep = NumDepositionSlow(n1,n2,n3,nbins,a_dn,1500.,tk,ustar,naerop,Rawet,nlim)
       amdep = MassDepositionSlow(n1,n2,n3,n4,nbins,a_dn,1500.,tk,ustar,naerop,maerop,Rawet,nlim)
       
       naerot = naerot - andiv
       naerot(2,:,:,:) = naerot(2,:,:,:) - andep
       maerot = maerot - amdiv
       maerot(2,:,:,:) = maerot(2,:,:,:) - amdep

       remaer(:,:,:) = amdep(:,:,:)

       ! Account for changes in in liquid water pot temperature
       nc = GetIndex(prtcl,'H2O')
       istr = (nc-1)*nbins+1
       iend = nc*nbins
       DO j = 3,n3-2
          DO i = 3,n2-2
             DO k = 1,n1
                tlt(k,i,j) = tlt(k,i,j) + SUM(amdiv(k,i,j,istr:iend))*(alvl/cp)*th(k,i,j)/tk(k,i,j)
             END DO
          END DO
       END DO

    END IF ! sed_aero

    IF (sed_cloud) THEN

       cndiv = NumDivergence(n1,n2,n3,ncld,tk,a_dn,rhowa,Rcwet,ncloudp,dzt,nlim)
       cmdiv = MassDivergence(n1,n2,n3,n4,ncld,tk,a_dn,rhowa,Rcwet,ncloudp,mcloudp,dzt,nlim)
    
       cndep = NumDepositionSlow(n1,n2,n3,ncld,a_dn,rhowa,tk,ustar,ncloudp,Rcwet,nlim)
       cmdep = MassDepositionSlow(n1,n2,n3,n4,ncld,a_dn,rhowa,tk,ustar,ncloudp,mcloudp,Rcwet,nlim)

       ncloudt = ncloudt - cndiv
       ncloudt(2,:,:,:) = ncloudt(2,:,:,:) - cndep
       mcloudt = mcloudt - cmdiv
       mcloudt(2,:,:,:) = mcloudt(2,:,:,:) - cmdep

       remcld(:,:,:) = cmdep(:,:,:)

       ! Account for changes in in liquid water pot temperature
       nc = GetIndex(prtcl,'H2O')
       istr = (nc-1)*ncld+1
       iend = nc*ncld
       DO j = 3,n3-2
          DO i = 3,n2-2
             DO k = 1,n1
                tlt(k,i,j) = tlt(k,i,j) + SUM(cmdiv(k,i,j,istr:iend))*(alvl/cp)*th(k,i,j)/tk(k,i,j)
             END DO
          END DO
       END DO

    END IF ! sed_cloud

    IF (sed_ice) THEN

       indiv = NumDivergence(n1,n2,n3,nice,tk,a_dn,rhoic,Riwet,nicep,dzt,prlim)
       imdiv = MassDivergence(n1,n2,n3,n4,nice,tk,a_dn,rhoic,Riwet,nicep,micep,dzt,prlim)

       indep = NumDepositionSlow(n1,n2,n3,nice,a_dn,rhoic,tk,ustar,nicep,Riwet,prlim)
       imdep = MassDepositionSlow(n1,n2,n3,n4,nice,a_dn,rhoic,tk,ustar,nicep,micep,Riwet,prlim)

       nicet = nicet - indiv 
       nicet(2,:,:,:) = nicet(2,:,:,:) - indep
       micet = micet - imdiv 
       micet(2,:,:,:) = micet(2,:,:,:) - imdep

       remice(:,:,:) = imdep(:,:,:)

       ! Account for changes in in liquid water pot temperature
       nc = GetIndex(prtcl,'H2O')
       istr = (nc-1)*nice+1
       iend = nc*nice
       DO j = 3,n3-2
          DO i = 3,n2-2
             DO k = 1,n1
                tlt(k,i,j) = tlt(k,i,j) + SUM(imdiv(k,i,j,istr:iend))*(alvi/cp)*th(k,i,j)/tk(k,i,j)
             END DO
          END DO
       END DO

    END IF ! sed_ice


    ! ---------------------------------------------------------
    ! SEDIMENTATION/DEPOSITION OF FAST PRECIPITATING PARTICLES
    IF (sed_precp) THEN
       CALL DepositionFast(n1,n2,n3,n4,nprc,tk,a_dn,rowt,Rpwet,nprecpp,mprecpp,tstep,dzt,prnt,prvt,remprc,prlim)
       
       ! Rain rate
       nc = GetIndex(prtcl,'H2O')
       istr = (nc-1)*nprc + 1
       iend = nc*nprc
       DO j = 3,n3-2
          DO i = 3,n2-2
             DO k = n1-2,2,-1

                nprecpt(k,i,j,:) = nprecpt(k,i,j,:) + prnt(k,i,j,:)/tstep
                mprecpt(k,i,j,:) = mprecpt(k,i,j,:) + prvt(k,i,j,:)/tstep

                kp1 = k + 1
                rrate(k,i,j) = -SUM(prvt(k,i,j,istr:iend)*2.) * alvl*0.5*(a_dn(k,i,j)+a_dn(kp1,i,j))

                tlt(k,i,j) = tlt(k,i,j) - SUM(prvt(k,i,j,istr:iend)/tstep)*(alvl/cp)*th(k,i,j)/tk(k,i,j)
                
                IF (sflg) v1(k) = v1(k) + rrate(k,i,j)*xnpts
             END DO
          END DO
       END DO
    END IF
    
    IF (sed_snow) THEN
       CALL DepositionFast(n1,n2,n3,n4,nsnw,tk,a_dn,rowt,Rswet,nsnowp,msnowp,tstep,dzt,srnt,srvt,remsnw,prlim)
           
       ! Snow rate
       nc = GetIndex(prtcl,'H2O')
       istr = (nc-1)*nsnw + 1
       iend = nc*nsnw
       DO j = 3,n3-2
          DO i = 3,n2-2
             DO k = n1-2,2,-1

                nsnowt(k,i,j,:) = nsnowt(k,i,j,:) + srnt(k,i,j,:)/tstep
                msnowt(k,i,j,:) = msnowt(k,i,j,:) + srvt(k,i,j,:)/tstep

                kp1 = k + 1
                srate(k,i,j) = -SUM(srvt(k,i,j,istr:iend)*1.) * alvi*0.5*(a_dn(k,i,j)+a_dn(kp1,i,j)) 

                tlt(k,i,j) = tlt(k,i,j) - SUM(srvt(k,i,j,istr:iend)/tstep)*(alvi/cp)*th(k,i,j)/tk(k,i,j)
                
                IF (sflg) v1(k) = v1(k) + srate(k,i,j)*xnpts
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
    if (sflg) call updtst(n1,'prc',1,v1,1)
    ! Aerosol removal statistics
    IF (sflg) CALL acc_removal(n2,n3,n4,remaer,remcld,remprc,remice,remsnw)
    IF (sflg) CALL cs_rem_set(n2,n3,n4,remaer,remcld,remprc,remice,remsnw)

  END SUBROUTINE !sedim_SALSA


  ! -----------------------------------------------------------------

  FUNCTION NumDivergence(n1,n2,n3,nn,tk,adn,pdn,rwet,numc,dzt,clim) RESULT(flxdiv)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3       ! Grid numbers
    INTEGER, INTENT(in) :: nn             ! Number of bins
    REAL, INTENT(in) :: tk(n1,n2,n3)      ! Absolute temprature
    REAL, INTENT(in) :: adn(n1,n2,n3)     ! Air density
    REAL, INTENT(in) :: pdn  ! Particle density
    REAL, INTENT(in) :: rwet(n1,n2,n3,nn) ! Particle wet radius
    REAL, INTENT(in) :: numc(n1,n2,n3,nn) ! Particle number concentrations
    REAL, INTENT(in) :: dzt(n1)           ! Inverse of grid level thickness
    REAL, INTENT(IN) :: clim                ! Concentration limit

    INTEGER :: i,j,k,kp1
    INTEGER :: bin

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
    REAL :: rhoref               ! Reference air density in STP conditions

    REAL :: rfl(n1,nn)
    REAL :: flxdiv(n1,n2,n3,nn)

    flxdiv = 0.

    rhoref = 1.01325e5/(287.*273.15)

    DO j = 3,n3-2

       DO i = 3,n2-2

          rfl = 0.

          DO k=n1-1,2,-1
             kp1 = k+1

             ! atm modelling Eq.4.54
             avis=1.8325e-5*(416.16/(tk(k,i,j)+120.0))*(tk(k,i,j)/296.16)**1.5
             kvis =  avis/adn(k,i,j) !actual density ???
             va = sqrt(8*kb*tk(k,i,j)/(pi*M)) ! thermal speed of air molecule
             lambda = 2*avis/(adn(k,i,j)*va) !mean free path

             ! Fluxes
             !------------------
             ! -- Calculate the *corrections* for small particles
             DO bin = 1,nn
                IF (numc(k,i,j,bin)<clim) CYCLE

                IF (rwet(k,i,j,bin) < 20.e-6 ) THEN
                   Kn = lambda/rwet(k,i,j,bin)
                   GG = 1.+ Kn*(A+B*exp(-C/Kn))
                   vc = (2.*(rwet(k,i,j,bin)**2)*(pdn-adn(k,i,j))*g/(9.*avis))*GG
                   vc = MIN(vc,1.)
                ELSE
                   vc = 2.e3*(2.*rwet(k,i,j,bin))*(rhoref/adn(k,i,j))**2
                   vc = MIN(vc,5.) 
                END IF

                ! Flux for the particle number
                IF ( k > 2 ) rfl(k,bin) = -numc(k,i,j,bin)*vc

             END DO

             flxdiv(k,i,j,:) = (rfl(kp1,:)-rfl(k,:))*dzt(k)

          END DO ! k

       END DO ! i
       
    END DO ! j
   
  END FUNCTION NumDivergence
  
  ! ----------------------------------------------------

  FUNCTION MassDivergence(n1,n2,n3,n4,nn,tk,adn,pdn,rwet,numc,mass,dzt,clim) RESULT(flxdiv)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3,n4       ! Grid numbers, number of chemical species
    INTEGER, INTENT(in) :: nn                ! Number of bins
    REAL, INTENT(in) :: tk(n1,n2,n3)         ! Absolute temprature
    REAL, INTENT(in) :: adn(n1,n2,n3)        ! Air density
    REAL, INTENT(in) :: pdn                  ! Particle density
    REAL, INTENT(in) :: rwet(n1,n2,n3,nn)    ! Particle wet radius
    REAL, INTENT(in) :: numc(n1,n2,n3,nn)    ! Particle number concentration
    REAL, INTENT(in) :: mass(n1,n2,n3,nn*n4) ! Particle mass mixing ratio
    REAL, INTENT(in) :: dzt(n1)              ! Inverse of grid level thickness
    REAL, INTENT(IN) :: clim                ! Concentration limit

    INTEGER :: i,j,k,kp1
    INTEGER :: bin,ss,bs

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
    REAL :: rhoref               ! Reference air density in STP conditions

    REAL :: rfl(n1,nn,n4)
    REAL :: flxdiv(n1,n2,n3,nn*n4)

    flxdiv = 0.

    rhoref = 1.01325e5/(287.*273.15)

    DO j = 3,n3-2

       DO i = 3,n2-2

          rfl = 0.
          
          DO k=n1-1,2,-1
             kp1 = k+1

             ! atm modelling Eq.4.54
             avis=1.8325e-5*(416.16/(tk(k,i,j)+120.0))*(tk(k,i,j)/296.16)**1.5
             kvis =  avis/adn(k,i,j) !actual density ???
             va = sqrt(8*kb*tk(k,i,j)/(pi*M)) ! thermal speed of air molecule
             lambda = 2*avis/(adn(k,i,j)*va) !mean free path

             ! Fluxes
             !------------------
             ! -- Calculate the *corrections* for small particles
             DO bin = 1,nn
                IF (numc(k,i,j,bin)<clim) CYCLE
        
                IF (rwet(k,i,j,bin) < 20.e-6 ) THEN
                   Kn = lambda/rwet(k,i,j,bin)
                   GG = 1.+ Kn*(A+B*exp(-C/Kn))
                   vc = (2.*(rwet(k,i,j,bin)**2)*(pdn-adn(k,i,j))*g/(9.*avis))*GG
                   vc = MIN(vc,1.)
                ELSE
                   vc = 2.e3*(2.*rwet(k,i,j,bin))*(rhoref/adn(k,i,j))**2
                   vc = MIN(vc,5.) 
                END IF

                ! Flux for the particle mass
                IF ( k > 2 ) THEN
                   DO ss = 1,n4
                      bs = (ss-1)*nn + bin
                      rfl(k,bin,ss) = -mass(k,i,j,bs)*vc
                      flxdiv(k,i,j,bs) = (rfl(kp1,bin,ss)-rfl(k,bin,ss))*dzt(k)
                   END DO
                ELSE
                   rfl(k,bin,:) = 0.
                   DO ss = 1,n4
                      bs = (ss-1)*nn + bin
                      flxdiv(k,i,j,bs) = (rfl(kp1,bin,ss)-rfl(k,bin,ss))*dzt(k)
                   END DO
                END IF

             END DO ! bin

          END DO ! k

       END DO ! i
  
    END DO ! j

  END FUNCTION MassDivergence

  ! ------------------------------------------

  FUNCTION NumDepositionSlow(n1,n2,n3,nn,adn,pdn,tk,ustar,numc,rwet,clim) RESULT(depflx)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3,nn
    REAL, INTENT(in) :: adn(n1,n2,n3) 
    REAL, INTENT(in) :: pdn
    REAL, INTENT(in) :: tk(n1,n2,n3)
    REAL, INTENT(in) :: ustar(n2,n3)
    REAL, INTENT(in) :: numc(n1,n2,n3,nn)
    REAL, INTENT(in) :: rwet(n1,n2,n3,nn)
    REAL, INTENT(IN) :: clim                ! Concentration limit

    INTEGER :: i,j,k,bin

    real, parameter :: A = 1.249 ! fundamentals of atm. modelling pg509
    real, parameter :: B = 0.42
    real, parameter :: C = 0.87
    real, parameter :: M = 4.8096e-26 ! average mass of one air molecule, eq2.3 fundamentals of atm.
                                      ! modelling [kg molec-1]

    REAL :: GG
    REAL :: St, Sc, Kn
    REAL :: lambda      ! Mean free path
    REAL :: mdiff       ! Particle diffusivity
    REAL :: avis,kvis   ! Air viscosity, kinematic viscosity
    REAL :: va          ! Thermal speed of air molecule
    REAL :: vc,vd       ! Particle fall speed, deposition velocity
    REAL :: rt
    REAL :: rhoref      ! Reference air density in STP conditions

    REAL :: rfl(nn)

    REAL :: depflx(n2,n3,nn)

    depflx = 0.

    rhoref = 1.01325e5/(287.*273.15)

    ! Fix level to lowest above ground level
    k = 2

    DO j = 3,n3-2

       DO i = 3,n2-2

          ! atm modelling Eq.4.54
          avis=1.8325e-5*(416.16/(tk(k,i,j)+120.0))*(tk(k,i,j)/296.16)**1.5
          kvis = avis/adn(k,i,j) !actual density ???
          va = sqrt(8.*kb*tk(k,i,j)/(pi*M)) ! thermal speed of air molecule
          lambda = 2.*avis/(adn(k,i,j)*va) !mean free path

          rfl = 0.

          DO bin = 1,nn
             IF (numc(k,i,j,bin) < clim) CYCLE

             Kn = lambda/rwet(k,i,j,bin)
             GG = 1.+ Kn*(A+B*exp(-C/Kn))

             IF (rwet(k,i,j,bin) < 20.e-6 ) THEN
                vc = (2.*(rwet(k,i,j,bin)**2)*(pdn-adn(k,i,j))*g/(9.*avis))*GG
                vc = MIN(vc,1.)
             ELSE
                vc = 2.e3*(2.*rwet(k,i,j,bin))*(rhoref/adn(k,i,j))**2
                vc = MIN(vc,5.) 
             END IF
             
             ! Particle diffusitivity  (15.29) in jacobson book
             mdiff = (kb*tk(k,i,j)*GG)/(6.0*pi*avis*rwet(k,i,j,bin))
             
             Sc = kvis/mdiff
             St = vc*ustar(i,j)**2.0/g*kvis
             if (St<0.01) St=0.01
             rt = 1.0/MAX(epsilon(1.0),(ustar(i,j)*(Sc**(-2.0/3.0)+10**(-3.0/St)))) ! atm chem&phy eq19.18

             vd = (1./rt) + vc

             rfl(bin) = -numc(k,i,j,bin)*vd

          END DO ! bin

          depflx(i,j,:) = -rfl(:)*dzt(k) 

       END DO ! i

    END DO ! j

  END FUNCTION NumDepositionSlow

  ! ------------------------------------------

  FUNCTION MassDepositionSlow(n1,n2,n3,n4,nn,adn,pdn,tk,ustar,numc,mass,rwet,clim) RESULT(depflx)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3,n4,nn
    REAL, INTENT(in) :: adn(n1,n2,n3) 
    REAL, INTENT(in) :: pdn
    REAL, INTENT(in) :: tk(n1,n2,n3)
    REAL, INTENT(in) :: ustar(n2,n3)
    REAL, INTENT(in) :: numc(n1,n2,n3,nn)
    REAL, INTENT(in) :: mass(n1,n2,n3,nn*n4)
    REAL, INTENT(in) :: rwet(n1,n2,n3,nn)
    REAL, INTENT(IN) :: clim                ! Concentration limit

    INTEGER :: i,j,k,bin
    INTEGER :: ss,bs

    real, parameter :: A = 1.249 ! fundamentals of atm. modelling pg509
    real, parameter :: B = 0.42
    real, parameter :: C = 0.87
    real, parameter :: M = 4.8096e-26 ! average mass of one air molecule, eq2.3 fundamentals of atm.
                                      ! modelling [kg molec-1]

    REAL :: GG
    REAL :: St, Sc, Kn
    REAL :: lambda      ! Mean free path
    REAL :: mdiff       ! Particle diffusivity
    REAL :: avis,kvis   ! Air viscosity, kinematic viscosity
    REAL :: va          ! Thermal speed of air molecule
    REAL :: vc,vd       ! Particle fall speed, deposition velocity
    REAL :: rt
    REAL :: rhoref      ! Reference air density in STP conditions

    REAL :: rfl(nn,n4)

    REAL :: depflx(n2,n3,nn*n4)

    depflx = 0.

    rhoref = 1.01325e5/(287.*273.15)

    ! Fix level to lowest above ground level
    k = 2

    DO j = 3,n3-2

       DO i = 3,n2-2

          ! atm modelling Eq.4.54
          avis=1.8325e-5*(416.16/(tk(k,i,j)+120.0))*(tk(k,i,j)/296.16)**1.5
          kvis = avis/adn(k,i,j) !actual density ???
          va = sqrt(8*kb*tk(k,i,j)/(pi*M)) ! thermal speed of air molecule
          lambda = 2*avis/(adn(k,i,j)*va) !mean free path

          rfl = 0.

          DO bin = 1,nn
             IF (numc(k,i,j,bin) < clim) CYCLE

             Kn = lambda/rwet(k,i,j,bin)
             GG = 1.+ Kn*(A+B*exp(-C/Kn))

             IF (rwet(k,i,j,bin) < 20.e-6 ) THEN
                vc = (2.*(rwet(k,i,j,bin)**2)*(pdn-adn(k,i,j))*g/(9.*avis))*GG
                vc = MIN(vc,1.)
             ELSE
                vc = 2.e3*(2.*rwet(k,i,j,bin))*(rhoref/adn(k,i,j))**2
                vc = MIN(vc,5.) 
             END IF

             ! Particle diffusitivity  (15.29) in jacobson book
             mdiff = (kb*tk(k,i,j)*GG)/(6.0*pi*avis*rwet(k,i,j,bin))
             
             Sc = kvis/mdiff
             St = vc*ustar(i,j)**2.0/g*kvis
             if (St<0.01) St=0.01
             rt = 1.0/MAX(epsilon(1.0),(ustar(i,j)*(Sc**(-2.0/3.0)+10**(-3.0/St)))) ! atm chem&phy eq19.18

             vd = (1./rt) + vc
    
             DO ss = 1,n4
                bs = (ss-1)*nn + bin
                rfl(bin,ss) = -mass(k,i,j,bs)*vd
                depflx(i,j,bs) = -rfl(bin,ss)*dzt(k)
             END DO

          END DO ! bin

       END DO ! i

    END DO ! j


  END FUNCTION MassDepositionSlow
  
  !------------------------------------------------------------------
  SUBROUTINE DepositionFast(n1,n2,n3,n4,nn,tk,adn,pdn,rwet,numc,mass,tstep,dzt,prnt,prvt,remprc,clim)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3,n4,nn
    REAL, INTENT(in) :: tk(n1,n2,n3)
    REAL, INTENT(in) :: adn(n1,n2,n3)
    REAL, INTENT(in) :: pdn
    REAL, INTENT(in) :: rwet(n1,n2,n3,nn)
    REAL, INTENT(in) :: numc(n1,n2,n3,nn)
    REAL, INTENT(in) :: mass(n1,n2,n3,nn*n4)
    REAL, INTENT(in) :: tstep
    REAL, INTENT(in) :: dzt(n1)
    REAL, INTENT(IN) :: clim                ! Concentration limit

    REAL, INTENT(out) :: prnt(n1,n2,n3,nn), prvt(n1,n2,n3,nn*n4)     ! Number and mass tendencies due to fallout
    REAL, INTENT(out) :: remprc(n2,n3,nn*n4)

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
    REAL :: rhoref      ! Reference air density in STP conditions

    ! For precipitation:
    REAL :: fd,fdmax,fdos ! Fall distance for rain drops, max fall distance, overshoot from nearest grid level
    REAL :: prnchg(n1,nn), prvchg(n1,nn,n4) ! Instantaneous changes in precipitation number and mass (volume)
 
    REAL :: prnumc, prvolc(n4)  ! Instantaneous source number and volume
    INTEGER :: kf, ni,fi
    LOGICAL :: prcdep  ! Deposition flag

    remprc(:,:,:) = 0.
    prnt(:,:,:,:) = 0.
    prvt(:,:,:,:) = 0.

    rhoref = 1.01325e5/(287.*273.15)

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
          
             prnumc = 0.
             prvolc = 0.

             ! Precipitation bin loop
             DO bin = 1,nn
                IF (numc(k,i,j,bin) < clim) CYCLE
             
                ! Terminal velocity
                IF (rwet(k,i,j,bin) < 20.e-6 ) THEN
                   Kn = lambda/rwet(k,i,j,bin)
                   GG = 1.+ Kn*(A+B*exp(-C/Kn))
                   vc = (2.*(rwet(k,i,j,bin)**2)*(pdn-adn(k,i,j))*g/(9.*avis))*GG
                ELSE
                   vc = 2.e3*(2.*rwet(k,i,j,bin))*(rhoref/adn(k,i,j))**2
                END IF
                vc = MIN(vc,10.)

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
             
                ! Put the drops to new positions (may be partially the original grid cell as well!)
                IF (prcdep) THEN
                   DO ni=1,n4
                      remprc(i,j,(ni-1)*nn+bin) = remprc(i,j,(ni-1)*nn+bin) +    &
                           prvolc(ni)*adn(k,i,j)*vc
                   END DO
                ELSE
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
                END IF ! prcdep
             
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
