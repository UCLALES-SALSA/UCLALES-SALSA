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

  use defs, only : alvl, alvi, rowt, pi, Rm, cp, kb, g

  implicit none

  logical, parameter :: khairoutdinov = .False.

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
    use mcrp_ice, only : micro_ice
    use stat, only : sflg, out_mcrp_data
    use grid, only : lev_sb,nzp,nxp,nyp,dtl,dzt,a_dn,a_theta,a_temp,a_rv, &
                  a_rsl, a_rc,CCN,a_rpp,a_npp,a_rp,a_rt,a_tt,a_rpt,a_npt, &
                  cldin,precip,sed_cloud,sed_precp
    integer, intent (in) :: level

    select case (level)
    case(0)
       ! Seifert and Beheng microphysics
       CALL micro_ice(lev_sb)
    case(2)
       IF (sflg) out_mcrp_data(:,:,:,:) = 0.
       if (sed_cloud)  &
            call sedim_cd(nzp,nxp,nyp,dtl,dzt,a_dn,a_theta,a_temp,a_rc,CCN,cldin,a_rt,a_tt)
    case(3)
       IF (sflg) out_mcrp_data(:,:,:,:) = 0.
       call mcrph(nzp,nxp,nyp,dtl,dzt,a_dn,a_theta,a_temp,a_rv,a_rsl,a_rc,CCN,a_rpp, &
                  a_npp,cldin,precip,a_rp,a_rt,a_tt,a_rpt,a_npt,sed_cloud,sed_precp)
    end select

  end subroutine micro

  !
  ! ---------------------------------------------------------------------
  ! MCRPH: calls microphysical parameterization
  !
  subroutine mcrph(n1,n2,n3,dtl,dzt,dn,th,tk,rv,rs,rc,ccn,rp,np,crate,rrate,  &
       rtp,rtt,tlt,rpt,npt,sed_cloud,sed_precp)
    use stat, only : sflg, out_mcrp_nout, out_mcrp_data, out_mcrp_list
    integer, intent (in) :: n1,n2,n3
    real, intent (in)                         :: dtl, dzt(n1), ccn
    real, dimension(n1,n2,n3), intent (in)    :: dn, th, tk, rv, rs
    real, dimension(n1,n2,n3), intent (inout) :: rc, rp, np, rpt, npt, rtp, rtt, tlt
    real, intent (out)                        :: crate(n1,n2,n3), rrate(n1,n2,n3)
    logical, intent (in)                      :: sed_cloud, sed_precp

    integer :: i, j, k
    REAL :: tmp_nr(n1,n2,n3), tmp_rr(n1,n2,n3), tmp_rt(n1,n2,n3)
    !
    ! Microphysics following Seifert Beheng (2001, 2005)

    ! Diagnostics (ignoring/avoiding negative values, and changes in number due to minimum and maximum volumes)
    if(sflg) CALL sb_var_stat('diag',2) ! Use rain concentrations
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             rp(k,i,j) = max(0., rp(k,i,j))
             np(k,i,j) = max(min(rp(k,i,j)/X_bnd,np(k,i,j)),rp(k,i,j)/X_max)
          end do
       end do
    end do
    if(sflg) CALL sb_var_stat('diag',3) ! ... simple difference

    ! Condensation/evaporation
    if(sflg) CALL sb_var_stat('cond',0) ! Use rain tendencies
    call wtr_dff_SB(n1,n2,n3,dn,rp,np,rs,rv,tk,rpt,npt)
    if(sflg) CALL sb_var_stat('cond',1) ! ... simple difference

    ! Autoconversion
    if(sflg) CALL sb_var_stat('auto',0) ! Use rain tendencies
    call auto_SB(n1,n2,n3,dn,rc,ccn,rp,rpt,npt)
    if(sflg) CALL sb_var_stat('auto',1) ! ... simple difference

    ! Accretion - coagulation
    if(sflg) CALL sb_var_stat('coag',0) ! Use rain tendenciess
    call accr_SB(n1,n2,n3,dn,rc,rp,np,rpt,npt)
    if(sflg) CALL sb_var_stat('coag',1) ! ... simple difference

    ! Apply tendencies
    if(sflg) CALL sb_var_stat('diag',2) ! Use rain concentrations
    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             rp(k,i,j) = rp(k,i,j) + max(-rp(k,i,j)/dtl,rpt(k,i,j))*dtl
             np(k,i,j) = np(k,i,j) + max(-np(k,i,j)/dtl,npt(k,i,j))*dtl
             rp(k,i,j) = max(0., rp(k,i,j))
             np(k,i,j) = max(min(rp(k,i,j)/X_bnd,np(k,i,j)),rp(k,i,j)/X_max)
          end do
       end do
    end do
    if(sflg) CALL sb_var_stat('diag',4) ! ... the real change compared with the tendency
    rpt(:,:,:)= 0.
    npt(:,:,:)= 0.

    ! Sedimentation
    if(sflg) CALL sb_var_stat('sedi',0) ! Use rain and total water tendencies
    rrate(:,:,:)=0.
    if (sed_precp) call sedim_rd(n1,n2,n3,dtl,dzt,dn,rp,np,tk,th,rrate,rtt,tlt,rpt,npt)

    ! Note: rc is not updated after autoconversion and accretion!
    crate(:,:,:)=0.
    if (sed_cloud) call sedim_cd(n1,n2,n3,dtl,dzt,dn,th,tk,rc,ccn,crate,rtt,tlt)
    if(sflg) CALL sb_var_stat('sedi',1) ! ... simple difference

  CONTAINS

      ! Collecting statistical outputs from the Seifert and Beheng micophysics
      !     Main program: out_mcrp_nout, out_mcrp_data, out_mcrp_list
      !     Calling subroutine: tmp_nr, tmp_rr and tmp_rt
      ! The first call (flag is 0 or 2) saves the current state and the second
      ! call (flag is 1, 3 or 4) is for calculating the change.
      SUBROUTINE sb_var_stat(prefix,flag)
        IMPLICIT NONE
        ! Input
        character (len=4), intent (in) :: prefix ! Process name
        INTEGER, INTENT(IN) :: flag ! Option flag
        ! Local
        INTEGER :: i
        !
        IF (out_mcrp_nout==0) RETURN
        !
        ! a) The first call
        IF (flag==0) THEN
            ! Save current tendency
            tmp_nr(:,:,:) = npt(:,:,:) ! Rain drop number
            tmp_rr(:,:,:) = rpt(:,:,:) ! Rain water mixing ratio
            tmp_rt(:,:,:) = rtt(:,:,:) ! Total water mixing ratio
            RETURN
        ELSEIF (flag==2) THEN
            ! Save current absolute concentration
            tmp_nr(:,:,:) = np(:,:,:)
            tmp_rr(:,:,:) = rp(:,:,:)
            tmp_rt(:,:,:) = rtp(:,:,:)
            RETURN
        ENDIF
        !
        ! b) The second call
        ! Find the requested ouput
        DO i=1,out_mcrp_nout
            IF ( prefix//'_nr' == out_mcrp_list(i) ) THEN
                ! Rain number (a_npt or a_npp)
                IF (flag==1) THEN
                    ! Calculate the change in tendency
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (npt(:,:,:) - tmp_nr(:,:,:))
                ELSEIF (flag==3) THEN
                    ! Calculate the change in absolute concentrations (divide by time step)
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (np(:,:,:) - tmp_nr(:,:,:))/dtl
                ELSEIF (flag==4) THEN
                    ! Compare the actual change in absolute concetration to the expected change
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (np(:,:,:) - tmp_nr(:,:,:))/dtl - npt(:,:,:)
                ENDIF
            ELSEIF ( prefix//'_rr' == out_mcrp_list(i) ) THEN
                ! Rain mixing ratio (a_rpt and a_rpp)
                IF (flag==1) THEN
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (rpt(:,:,:) - tmp_rr(:,:,:))
                ELSEIF (flag==3) THEN
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (rp(:,:,:) - tmp_rr(:,:,:))/dtl
                ELSEIF (flag==4) THEN
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (rp(:,:,:) - tmp_rr(:,:,:))/dtl - rpt(:,:,:)
                ENDIF
            ELSEIF ( prefix//'_rt' == out_mcrp_list(i) ) THEN
                ! Total water (a_rt and a_rp)
                IF (flag==0) THEN
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (rtt(:,:,:) - tmp_rt(:,:,:))
                ELSEIF (flag==3) THEN
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (rtp(:,:,:) - tmp_rt(:,:,:))/dtl
                ELSEIF (flag==4) THEN
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (rtp(:,:,:) - tmp_rt(:,:,:))/dtl - rtt(:,:,:)
                ENDIF
            ELSEIF ( prefix//'_rc' == out_mcrp_list(i) ) THEN
                ! Cloud water - diagnostic: when water vapor mixing ratio is constant, rc = const + rt - rr
                IF (flag==1) THEN
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (rtt(:,:,:)-rpt(:,:,:) - (tmp_rt(:,:,:)-tmp_rr(:,:,:)))
                ELSEIF (flag==3) THEN
                    out_mcrp_data(:,:,:,i) = out_mcrp_data(:,:,:,i) + (rtp(:,:,:)-rp(:,:,:) - (tmp_rt(:,:,:)-tmp_rr(:,:,:)))/dtl
                ELSEIF (flag==4) THEN
                ENDIF
            ENDIF
        ENDDO
      END SUBROUTINE sb_var_stat

  end subroutine mcrph
  !
  ! ---------------------------------------------------------------------
  ! WTR_DFF_SB: calculates the evolution of the both number- and
  ! mass mixing ratio large drops due to evaporation in the absence of
  ! cloud water.
  !
  subroutine wtr_dff_SB(n1,n2,n3,dn,rp,np,rs,rv,tk,rpt,npt)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: tk(n1,n2,n3), np(n1,n2,n3), rp(n1,n2,n3),        &
         rs(n1,n2,n3),rv(n1,n2,n3), dn(n1,n2,n3)
    real, intent (inout) :: rpt(n1,n2,n3), npt(n1,n2,n3)

    real, parameter    :: Kt = 2.5e-2    ! conductivity of heat [J/(sKm)]
    real, parameter    :: Dv = 3.e-5     ! diffusivity of water vapor [m2/s]

    integer             :: i, j, k
    real                :: Xp, Dp, G, S, cerpt, cenpt

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1
             if (rp(k,i,j) > 0. .AND. rv(k,i,j)<rs(k,i,j)) then
                Xp = rp(k,i,j)/ (np(k,i,j)+eps0)
                Xp = MIN(MAX(Xp,X_bnd),X_max)
                Dp = ( Xp / prw )**(1./3.)

                G = 1. / (1. / (dn(k,i,j)*rs(k,i,j)*Dv) + &
                     alvl*(alvl/(Rm*tk(k,i,j))-1.) / (Kt*tk(k,i,j)))
                S = rv(k,i,j)/rs(k,i,j) - 1.

                cerpt = 2. * pi * Dp * G * S * np(k,i,j)
                cenpt = cerpt * np(k,i,j) / rp(k,i,j)
                rpt(k,i,j)=rpt(k,i,j) + cerpt
                npt(k,i,j)=npt(k,i,j) + cenpt
             end if
          end do
       end do
    end do

  end subroutine wtr_dff_SB
  !
  ! ---------------------------------------------------------------------
  ! AUTO_SB:  calculates the evolution of mass- and number mxg-ratio for
  ! drizzle drops due to autoconversion. The autoconversion rate assumes
  ! f(x)=A*x**(nu_c)*exp(-Bx), an exponential in drop MASS x. It can
  ! be reformulated for f(x)=A*x**(nu_c)*exp(-Bx**(mu)), where formu=1/3
  ! one would get a gamma dist in drop diam -> faster rain formation.
  !
  subroutine auto_SB(n1,n2,n3,dn,rc,CCN,rp,rpt,npt)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dn(n1,n2,n3), rc(n1,n2,n3), CCN, rp(n1,n2,n3)
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
                au = k_au * dn(k,i,j) * rc(k,i,j)**2 * Xc**2
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
  subroutine accr_SB(n1,n2,n3,dn,rc,rp,np,rpt,npt)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: rc(n1,n2,n3), rp(n1,n2,n3), np(n1,n2,n3), dn(n1,n2,n3)
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
             ! accretion
             if (rc(k,i,j) > 0. .and. rp(k,i,j) > 0.) then
                tau = 1.0-rc(k,i,j)/(rc(k,i,j)+rp(k,i,j)+eps0)
                tau = MIN(MAX(tau,eps0),1.)
                phi = (tau/(tau+k_1))**4
                ac  = k_r * rc(k,i,j) * rp(k,i,j) * phi * sqrt(rho_0*dn(k,i,j))
                !
                ! Khairoutdinov and Kogan
                if (khairoutdinov) then
                   ac = Cac * (rc(k,i,j) * rp(k,i,j))**Eac
                end if
                !
                rpt(k,i,j) = rpt(k,i,j) + ac
             end if

             ! self-collection
             sc = k_r * np(k,i,j) * rp(k,i,j) * sqrt(rho_0*dn(k,i,j))
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
   subroutine sedim_rd(n1,n2,n3,dt,dzt,dens,rp,np,tk,th,rrate,rtt,tlt,rpt,npt)

     integer, intent (in)                      :: n1,n2,n3
     real, intent (in)                         :: dt, dzt(n1)
     real, intent (in),    dimension(n1,n2,n3) :: dens, rp, np, th, tk
     real, intent (inout), dimension(n1,n2,n3) :: rrate
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
              Dm = ( 6. / (rowt*pi) * Xp )**(1./3.)
              mu = cmur1*(1.+tanh(cmur2*(Dm-cmur3)))
              Dp = (Dm**3/((mu+3.)*(mu+2.)*(mu+1.)))**(1./3.)
              ! fall speeds
              vn(k) = sqrt(rho_0/dens(k,i,j))*(a2 - b2*(1.+c2*Dp)**(-(1.+mu)))
              vr(k) = sqrt(rho_0/dens(k,i,j))*(a2 - b2*(1.+c2*Dp)**(-(4.+mu)))
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
                 tot = tot + dens(kk,i,j)*(np(kk,i,j)+nslope(kk)*(1.-cc))*cc/dzt(kk)
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
                 tot = tot + dens(kk,i,j)*(rp(kk,i,j)+rslope(kk)*(1.-cc))*cc/dzt(kk)
                 zz  = zz + 1./dzt(kk)
                 kk  = kk + 1
                 cc  = min(1.,cr(kk) - zz*dzt(kk))
              enddo
              rfl(k) = -tot /dt

              kp1=k+1
              flxdiv = (rfl(kp1)-rfl(k))*dzt(k)/dens(k,i,j)
              rpt(k,i,j) =rpt(k,i,j)-flxdiv
              rtt(k,i,j) =rtt(k,i,j)-flxdiv
              tlt(k,i,j) =tlt(k,i,j)+flxdiv*(alvl/cp)*th(k,i,j)/tk(k,i,j)
              npt(k,i,j) = npt(k,i,j)-(nfl(kp1)-nfl(k))*dzt(k)/dens(k,i,j)
              rrate(k,i,j) = rrate(k,i,j) -rfl(k) * alvl
           end do
        end do
     end do

   end subroutine sedim_rd
  !
  ! ---------------------------------------------------------------------
  ! SEDIM_CD: calculates the cloud-droplet sedimentation flux and its effect
  ! on the evolution of r_t and theta_l assuming a log-normal distribution
  !
  subroutine sedim_cd(n1,n2,n3,dt,dzt,dens,th,tk,rc,CCN,rrate,rtt,tlt)

    integer, intent (in):: n1,n2,n3
    real, intent (in)   :: dt, dzt(n1), CCN
    real, intent (in),   dimension(n1,n2,n3) :: dens,th,tk,rc
    real, intent (inout),dimension(n1,n2,n3) :: rrate
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
          rfl(:) = 0.
          do k=n1-1,2,-1
             if (rc(k,i,j) > 0.) then
                Xc = rc(k,i,j) / CCN
                Dc = ( Xc / prw )**(1./3.)
                Dc = MIN(MAX(Dc,D_min),D_bnd)
                vc = min(c*(Dc*0.5)**2 * exp(5.0*(log(sgg))**2),1./(dzt(k)*dt))
                rfl(k) = - dens(k,i,j) * rc(k,i,j) * vc
             end if
             !
             kp1=k+1
             flxdiv = (rfl(kp1)-rfl(k))*dzt(k)/dens(k,i,j)
             rtt(k,i,j) = rtt(k,i,j)-flxdiv
             tlt(k,i,j) = tlt(k,i,j)+flxdiv*(alvl/cp)*th(k,i,j)/tk(k,i,j)
             rrate(k,i,j) = rrate(k,i,j) -rfl(k) * alvl
          end do
       end do
    end do

  end subroutine sedim_cd

end module mcrp
