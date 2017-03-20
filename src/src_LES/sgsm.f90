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
module sgsm

  use stat, only : sflg, updtst, acc_tend, sgsflxs, sgs_vel
  use util, only : tridiff, vel_bc
  implicit none
!
! setting the prandtl number to a value less than zero enforces an exponential
! decay in its value as a function of height from the surface with a scale
! height of 100m.  This is the default
!

  real, parameter     :: tkemin=1.e-20
  real :: csx = 0.23
  real :: prndtl = 0.3333333333

  real, allocatable, dimension (:,:) :: sxy1, sxy2, sxy3, sxz1, sxz2, sxz3    &
       , sxz4, sxz5, sxz6  ,szx1, szx2, szx3, szx4, szx5
  real, allocatable, dimension (:)   :: sz1, sz2, sz3, sz4, sz5, sz6, sz7, sz8

  integer :: k, i, j, indh, req(16)
  real    :: dti
  logical, save :: Initialized = .false.

contains

  subroutine diffuse_init(n1,n2,n3)

    integer, intent (in) :: n1, n2, n3

    integer :: nm

    nm = n1-1

    allocate(sxy1(n2,n3), sxy2(n2,n3), sxy3(n2,n3), sxz1(n2,nm), sxz2(n2,nm))
    allocate(sxz3(n2,nm), sxz4(n2,nm), sxz5(n2,nm), sxz6(n2,nm))
    allocate(szx1(n1,n2), szx2(n1,n2), szx3(n1,n2), szx4(n1,n2), szx5(n1,n2))
    allocate(sz1(n1),sz2(n1),sz3(n1),sz4(n1),sz5(n1),sz6(n1),sz7(n1),sz8(n1))

    initialized = .true.
  end subroutine

  !
  ! ---------------------------------------------------------------------
  ! SUBROUTINE DIFFUSE: Driver for calculating sub-grid fluxes (thus it
  ! includes call to surface routines)  Depending on value of ISGSTYP,
  ! the model computes diffusivities based on smaorinsky equilibrium model
  ! or a subgrid tke model.
  !
  ! Modified for Level 4:
  ! For water vapour: a_rv replaced by rx, which contains values
  !                   a_rv if level < 4 and a_rp if level == 4.
  !
  ! For total water amount: a_rp is replaced by rxt which is
  !                         similarly either a_rp (level < 4)
  !                         or a_rp + a_rc (level == 4).
  !
  ! Juha Tonttila, FMI, 2014
  !
  subroutine diffuse

    use grid, only : a_up, a_uc, a_ut, a_vp, a_vc, a_vt, a_wp, a_wc, a_wt    &
         , a_rv, a_rc, a_rp, a_tp, a_tt, a_sp, a_st, a_qt, a_qp, a_pexnr, a_theta  &
         , a_temp, a_rsl, nscl, nxp, nyp    &
         , nzp, zm, dxi, dyi, dzt, dzm, dtlt, dtlv , th00, dn0  &
         , pi0, pi1, newsclr, level, isgstyp, uw_sfc, vw_sfc, ww_sfc, wt_sfc &
         , wq_sfc

    use util, only         : get_avg3
    use mpi_interface, only: cyclics, cyclicc
    use thrm, only         : bruvais, fll_tkrs

    integer :: n
    REAL :: rx(nzp,nxp,nyp), rxt(nzp,nxp,nyp), a_tmp1(nzp,nxp,nyp), &
        a_tmp2(nzp,nxp,nyp), a_tmp3(nzp,nxp,nyp), a_tmp4(nzp,nxp,nyp), &
        a_tmp5(nzp,nxp,nyp), a_tmp6(nzp,nxp,nyp)

    SELECT CASE(level)
       CASE(1,2,3)
          rx = a_rv
          rxt = a_rp
       CASE(4,5)
          rx = a_rp
          rxt = a_rp + a_rc
    END SELECT


    if (.not.Initialized) call diffuse_init(nzp, nxp, nyp)
    !
    ! ----------
    ! Calculate Deformation and stability for SGS calculations
    !
    call fll_tkrs(nzp,nxp,nyp,a_theta,a_pexnr,pi0,pi1,dn0,th00,a_temp,rs=a_rsl)

    call bruvais(nzp,nxp,nyp,level,a_theta,a_tp,rxt,a_rsl,a_tmp3,dzm,th00)

    !
    ! the a_ut, a_wt, a_ut arrays are first used when the diffusive tendencies
    ! are calculated and applied.  Until then use them as scratch by
    ! associating them with scratch pointers (a-c)
    !
    call deform(nzp,nxp,nyp,dzm,dzt,dxi,dyi,a_up,a_vp,a_wp,a_tmp5,a_tmp6     &
         ,a_tmp4,a_tmp2)

    ! ----------
    ! Calculate Eddy Viscosity/Diffusivity according to specified SGS model
    !
    select case (isgstyp)
    case (1)
       call smagor(nzp,nxp,nyp,sflg,dxi,dn0,a_tmp3,a_tmp2,a_tmp1,zm)
    case (2)
       call deardf(nzp,nxp,nyp,sflg,dxi,zm,dn0,a_qp,a_qt,a_tmp3,a_tmp2,a_tmp1)

       call solv_tke(nzp,nxp,nyp,a_tmp3,a_tmp1,a_qp,a_qt,dn0,dzm,dzt,dxi,dyi  &
            ,dtlt)
    end select
    !
    ! Diffuse momentum
    !
    if (sflg) call acc_tend(nzp,nxp,nyp,a_uc,a_vc,a_wc,a_ut,a_vt,a_wt         &
         ,sz4,sz5,sz6,1,'sgs')

    call diff_prep(nzp,nxp,nyp,a_tmp5,a_tmp6,a_tmp4,a_tmp1)
    sxy1=0.; sxy2=0.

    call diff_vpt(nzp,nxp,nyp,dn0,dzm,dzt,dxi,dyi,dtlv,vw_sfc,sxy2,a_tmp6     &
         ,a_tmp5,a_tmp1,a_vp,a_wp,a_vt,sz2)

    call diff_upt(nzp,nxp,nyp,dn0,dzm,dzt,dxi,dyi,dtlv,uw_sfc,sxy1,a_tmp5     &
         ,a_tmp1,a_up,a_wp,a_ut,sz1)

    call diff_wpt(nzp,nxp,nyp,dn0,dzm,dzt,dyi,dxi,dtlv,ww_sfc,sxy1,a_tmp4     &
         ,a_tmp1,a_wp,a_up,a_wt,sz3)

    call cyclics(nzp,nxp,nyp,a_wt,req)
    call cyclicc(nzp,nxp,nyp,a_wt,req)
    call cyclics(nzp,nxp,nyp,a_vt,req)
    call cyclicc(nzp,nxp,nyp,a_vt,req)
    call cyclics(nzp,nxp,nyp,a_ut,req)
    call cyclicc(nzp,nxp,nyp,a_ut,req)

    if (sflg) then
       call sgs_vel(nzp,nxp,nyp,sz1,sz2,sz3)
       call acc_tend(nzp,nxp,nyp,a_uc,a_vc,a_wc,a_ut,a_vt,a_wt,sz4,sz5,sz6    &
            ,2,'sgs')
    end if
    !
    ! Diffuse scalars
    !
    a_tt=0.
    do n=1,nscl
       call newsclr(n)
       sxy1=0.
       sxy2=0.
       if ( associated(a_tp,a_sp) ) sxy1=wt_sfc
       if ( associated(a_tp,a_sp) ) sxy2=wt_sfc
       if ( associated(a_rp,a_sp) ) sxy1=wq_sfc

       if (sflg) a_tmp1=0.
       if ( isgstyp <= 1) then
          call diffsclr(nzp,nxp,nyp,dtlt,dxi,dyi,dzm,dzt,dn0,sxy1,sxy2   &
               ,a_sp,a_tmp2,a_st,a_tmp1)
       else if ( .not.associated(a_qp,a_sp) ) then
          call diffsclr(nzp,nxp,nyp,dtlt,dxi,dyi,dzm,dzt,dn0,sxy1,sxy2   &
               ,a_sp,a_tmp2,a_st,a_tmp1)
       end if
       if (sflg) then
          call get_avg3(nzp,nxp,nyp,a_tmp1,sz1)
          call updtst(nzp,'sgs',n,sz1,1)
          if (associated(a_sp,a_tp))                                          &
             call sgsflxs(nzp,nxp,nyp,level,a_tmp3,rx,a_theta,a_tmp1,'tl')
          if (associated(a_sp,a_rp))                          &
             call sgsflxs(nzp,nxp,nyp,level,a_tmp3,rx,a_theta,a_tmp1,'rt')
       endif
       call cyclics(nzp,nxp,nyp,a_st,req)
       call cyclicc(nzp,nxp,nyp,a_st,req)
    enddo

  end subroutine diffuse
  !
  ! ---------------------------------------------------------------------
  ! subroutine deform: computes the components of the deviatoric strain
  ! tensor, then computes s = du_i/dx_j(du_i/dx_j + du_j/dx_i) at a
  ! thermo point, the dummy arrays are respectively:
  !
  ! sxz1 -> div;
  ! szx2 -> s11;
  ! szx3 -> s33,
  ! szx4 -> s13
  !
  subroutine deform(n1,n2,n3,dzm,dzt,dx,dy,u,v,w,s12,s22,s23,s)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
    real, intent(in)    :: dzm(n1),dzt(n1),dx,dy

    real, intent(out)   :: s(n1,n2,n3), s22(n1,n2,n3)
    real, intent(out)   :: s12(n1,n2,n3),s23(n1,n2,n3)

    integer :: ip,im,jp,jm,kp
    real    :: y1a,y2a,y3a,y1b,y2b,y3b
    real    :: s11_wpt,s22_wpt,s33_wpt,s12_wpt,s13_wpt,s23_wpt
    !
    ! calculate components of the stress tensor at their natural locations
    !
    do j=1,n3
       jm=max(j-1,1)
       jp=min(j+1,n3)
       do i=1,n2
          im=max(i-1,1)
          ip=min(i+1,n2)
          do k=1,n1
             szx2(k,i) = 2.*(u(k,i,j)-u(k,im,j))*dx
             s22(k,i,j)= 2.*(v(k,i,j)-v(k,i,jm))*dy
             szx3(k,i) = 2.*(w(k,i,j)-w(max(1,k-1),i,j))*dzt(k)
             s12(k,i,j)= (u(k,i,jp)-u(k,i,j))*dy + (v(k,ip,j)-v(k,i,j))*dx
          enddo
          szx3(1,i)  = 0.
          szx3(n1,i) = 0.

          do k=1,n1-1
             szx1(k,i) = 0.333333*(szx2(k,i)+s22(k,i,j)+szx3(k,i))
             szx4(k,i) =(u(k+1,i,j)-u(k,i,j))*dzm(k)+(w(k,ip,j)-w(k,i,j))*dx
             s23(k,i,j)=(v(k+1,i,j)-v(k,i,j))*dzm(k)+(w(k,i,jp)-w(k,i,j))*dy
          end do
          szx1(n1,i) = 0.333333*(szx2(n1,i)+s22(n1,i,j)+szx3(n1,i))
       end do
       !
       ! average to a w-point
       !
       do i=1,n2
          im=max(i-1,1)
          do k=1,n1-1
             kp=k+1
             y1a=(szx2(k,i)-szx1(k,i))
             y2a=(s22(k,i,j)-szx1(k,i))
             y3a=(szx3(k,i)-szx1(k,i))
             y1b=(szx2(kp,i)-szx1(kp,i))
             y2b=(s22(kp,i,j)-szx1(kp,i))
             y3b=(szx3(kp,i)-szx1(kp,i))

             s11_wpt=0.5*(y1a*y1a+y1b*y1b)
             s22_wpt=0.5*(y2a*y2a+y2b*y2b)
             s33_wpt=0.5*(y3a*y3a+y3b*y3b)
             s12_wpt=0.125*(s12(k,i,j)*s12(k,i,j) + s12(kp,i,j)*s12(kp,i,j) &
                  +s12(k,im,j)*s12(k,im,j) + s12(k,i,jm)*s12(k,i,jm)        &
                  +s12(kp,im,j)*s12(kp,im,j)+s12(k,im,jm)*s12(k,im,jm)      &
                  +s12(kp,i,jm)*s12(kp,i,jm)+s12(kp,im,jm)*s12(kp,im,jm))
             s13_wpt=0.5*(szx4(k,i)*szx4(k,i)+szx4(k,im)*szx4(k,im))
             s23_wpt=0.5*(s23(k,i,j)*s23(k,i,j)+s23(k,i,jm)*s23(k,i,jm))
             s(k,i,j)= 0.5*(s11_wpt+s22_wpt+s33_wpt)                        &
                  + s12_wpt + s13_wpt + s23_wpt
          end do
       end do
    end do

  end subroutine deform
  !
  ! ----------------------------------------------------------------------
  ! Subroutine smagor:  computes visc/diff, upon entering the routine
  ! kh is filled with s2 and ri is filled with N^2.  On statistical
  ! timsteps, SGS energy, dissipation, viscosity, diffusivity and
  ! lengthscales are stored.
  !
  subroutine smagor(n1,n2,n3,sflg,dxi,dn0,ri,kh,km,zm)

    use defs, only          : pi, vonk
    use stat, only          : tke_sgs
    use util, only          : get_avg3, get_cor3
    use mpi_interface, only : cyclics, cyclicc

    implicit none

    logical, intent(in) :: sflg
    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: dxi,zm(n1),dn0(n1)
    real, intent(inout) :: ri(n1,n2,n3),kh(n1,n2,n3)
    real, intent(out)   :: km(n1,n2,n3)

    real    :: delta,pr,yplus

    pr    = abs(prndtl)

    delta = 1./dxi
    delta = (zm(2)/dxi/dxi)**0.333333333

    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1-1
             yplus = max(zm(2),min(zm(k),(zm(n1-1)-zm(k))))
             ri(k,i,j) = max( -1., ri(k,i,j)/(kh(k,i,j) + 1.e-12) )
             km(k,i,j) = sqrt(max(0.,kh(k,i,j)*(1.-ri(k,i,j)/pr))) &
                  *0.5*(dn0(k)+dn0(k+1))/(1./(delta*csx)**2+1./(yplus*vonk)**2)
          enddo
          kh(n1,i,j)   = kh(n1-1,i,j)
          km(n1,i,j)   = km(n1-1,i,j)
       enddo
    enddo

    call cyclics(n1,n2,n3,km,req)
    call cyclicc(n1,n2,n3,km,req)


    if (sflg) then
       call get_cor3(n1,n2,n3,km,km,sz1)
       call get_cor3(n1,n2,n3,km,kh,sz2)
       call updtst(n1,'sgs',-2,sz2,1)      ! dissipation
       do k=1,n1
          tke_sgs(k) = sz1(k)/(delta*pi*(csx*0.18))**2
          sz1(k) = 1./sqrt(1./(delta*csx)**2+1./(zm(k)*vonk+0.001)**2)
       end do
       call updtst(n1,'sgs',-1,tke_sgs,1) ! sgs tke
       call updtst(n1,'sgs',-5,sz1,1)      ! mixing length
       call updtst(n1,'sgs',-6,sz1,1)      ! dissipation lengthscale
    end if

    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
            kh(k,i,j) = km(k,i,j)/pr
            if (prndtl < 0.) then
               kh(k,i,j) = kh(k,i,j) * exp(zm(k)/(-100.))
            end if
          enddo
       enddo
    enddo
    call cyclics(n1,n2,n3,kh,req)
    call cyclicc(n1,n2,n3,kh,req)

    if (sflg) then
       call get_avg3(n1,n2,n3,km,sz3)
       call updtst(n1,'sgs',-3,sz3,1)  ! eddy viscosity
       call get_avg3(n1,n2,n3,kh,sz2)
       call updtst(n1,'sgs',-4,sz2,1)  ! eddy diffusivity
    end if

  end subroutine smagor
  !
  ! ----------------------------------------------------------------------
  ! Subroutine deardf:  computes visc/diff according to the prognostic
  ! TKE formulation of Deardorff (1980). Also, forcing of TKE is
  ! calculated.  Lastly, on statistical timsteps, SGS energy, dissipation,
  ! viscosity, diffusivity and lengthscales are stored.
  !
  ! upon entering:
  !   kh = deform
  !   xx = en2
  !   yy =  ---
  !   zz = deform
  !
  subroutine deardf(n1,n2,n3,sflg,dxi,zm,dn0,tke,tket,xx,zz,yy)

    use defs, only : vonk
    use stat, only : tke_sgs
    use util, only : get_avg3

    implicit none

    real, parameter :: cm=0.1,eps=1.e-12,cs=0.82,c_e1=0.225,c_e2=0.705

    logical, intent(in) :: sflg
    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: dxi,dn0(n1),zm(n1)
    real, intent(in)    :: tke(n1,n2,n3)
    real, intent(inout) :: xx(n1,n2,n3),zz(n1,n2,n3)
    real, intent(out)   :: yy(n1,n2,n3),tket(n1,n2,n3)

    real    :: ln,lm,ch
    !
    ! ----------
    ! Get neutral and effective lengthscales (ln,lm). Update tke tendency
    ! due to shear and buoyancy, fill:
    !   xx = dissipation lengthscale for later use
    !   yy = eddy viscosity (density weighted)
    !   zz = eddy diffusivity (density weighted)
    !
    do k=1,n1
       sz6(k) = 0.
       sz5(k) = 0.
       sz4(k) = 0.
    end do

    do j=1,n3
       do i=1,n2
          do k=2,n1-2
             ln   = sqrt(1./(dxi**2 + (0.23/(zm(k)*vonk))**2) )
             lm   = min(cs*sqrt(tke(k,i,j)/max(eps,xx(k,i,j))),ln)
             ch   = (1.+2.*lm/ln)
             tket(k,i,j) = cm*lm*sqrt(tke(k,i,j))*(zz(k,i,j) - ch*xx(k,i,j))
             xx(k,i,j) = lm/(c_e1 + c_e2*lm/ln)
             yy(k,i,j) = cm*lm*sqrt(tke(k,i,j))*.5*(dn0(k)+dn0(k+1))
             zz(k,i,j) = ch*yy(k,i,j)
             if (sflg) then
                sz6(k)  = sz6(k) + lm
                sz5(k)  = sz5(k) + lm*ch
                sz4(k)  = sz4(k) + (tke(k,i,j)**1.5)/xx(k,i,j)
             end if
          end do
       end do
    end do
    !
    ! --- accumululate SGS model stats
    !
    if (sflg) then
       call get_avg3(n1,n2,n3,tke,sz1)
       call get_avg3(n1,n2,n3,yy,sz2)
       call get_avg3(n1,n2,n3,zz,sz3)

       do k=1,n1
          sz6(k) = sz6(k)/real(n2*n3)
          sz5(k) = sz5(k)/real(n2*n3)
          sz4(k) = sz4(k)/real(n2*n3)
          tke_sgs(k) = sz1(k)
       end do
       call updtst(n1,'sgs',-1,sz1,1)
       call updtst(n1,'sgs',-2,sz4,1)
       call updtst(n1,'sgs',-3,sz2,1)
       call updtst(n1,'sgs',-4,sz3,1)
       call updtst(n1,'sgs',-5,sz5,1)
       call updtst(n1,'sgs',-6,sz6,1)
    end if

  end subroutine deardf
  !
  ! ----------------------------------------------------------------------
  ! subroutine solv_tke: solves for dissipation and diffusion of tke,
  ! using an implicit solver for vertical diffusion and dissipation
  !
  subroutine  solv_tke(n1,n2,n3,le,km,ep,et,dn0,dzm,dzt,dxi,dyi,dt)

    implicit none

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dn0(n1),dzm(n1),dzt(n1),dxi,dyi,dt
    real, intent (in)    :: km(n1,n2,n3),ep(n1,n2,n3),le(n1,n2,n3)
    real, intent (inout) :: et(n1,n2,n3)

    !
    ! Add horizontal diffusivities to tendency arrays
    !
    do j=2,n3-1
       do i=2,n2-1
          do k=2,n1-2
             et(k,i,j)=et(k,i,j)+(                                           &
                  ((km(k,i+1,j)+km(k,i,j))*(ep(k,i+1,j)-ep(k,i,j))           &
                  -(km(k,i,j)+km(k,i-1,j))*(ep(k,i,j)-ep(k,i-1,j)))*dxi*dxi  &
                  +((km(k,i,j+1)+km(k,i,j))*(ep(k,i,j+1)-ep(k,i,j))          &
                  -(km(k,i,j)+km(k,i,j-1))*(ep(k,i,j)-ep(k,i,j-1)))*dyi*dyi  &
                  )/((dn0(k)+dn0(k+1))*.5)
          enddo
       enddo
    enddo
    !
    ! Set up and solve implicit system
    !
    dti  = 1.0/dt
    do k=1,n1
       sz7(k) = 0.
    end do

    do j=2,n3-1
       indh=0
       do i=2,n2-1
          indh=indh+1
          sz7(1)=dzt(2)*(km(2,i,j)+km(1,i,j))
          do k=2,n1-2
             sz7(k)       = dzt(k+1)*(km(k+1,i,j)+km(k,i,j))
             sxz4(indh,k) = 0.5*(dn0(k)+dn0(k+1))*(ep(k,i,j)+dt*et(k,i,j))
             sxz3(indh,k) = -dt*dzm(k)*sz7(k)
             sxz2(indh,k) = -dt*dzm(k)*sz7(k-1)
             sxz1(indh,k) = (dn0(k)+dn0(k+1))*0.5                             &
                  *(1.+dt*sqrt(ep(k,i,j))/le(k,i,j))-sxz2(indh,k)-sxz3(indh,k)

          enddo
       enddo

       call tridiff(n2,n1-2,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
       !
       ! Back out tendencies from implicit solver
       !
       indh=0
       do i=2,n2-1
          indh=indh+1
          do k=2,n1-2
             et(k,i,j)=dti*(sxz1(indh,k)-ep(k,i,j))
          enddo
       enddo
    enddo

  end subroutine solv_tke
  !
  ! ----------------------------------------------------------------------
  ! subroutine diff_prep: multiplies the strain components computed in
  ! "deform" by the appropriately averaged value of km, in preperation
  ! for use by the diffusion routines
  !
  subroutine  diff_prep(n1,n2,n3,s12,s22,s23,km)

    integer, intent(in) :: n1,n2,n3
    real, intent(inout) :: s22(n1,n2,n3)
    real, intent(inout) :: s12(n1,n2,n3),s23(n1,n2,n3)
    real, intent(inout) :: km(n1,n2,n3)

    integer :: ip,jp

    do j=2,n3-1
       jp=min(j+1,n3)
       do i=2,n2-1
          ip=min(i+1,n2)
          do k=2,n1-1
             s22(k,i,j)=-s22(k,i,j)*0.5*(km(k,i,j)+km(k-1,i,j))
             s12(k,i,j)= -s12(k,i,j)*0.125*(km(k,i,j)+km(k,ip,j)+km(k,i,jp) &
                  +km(k,ip,jp)+km(k-1,i,j)+km(k-1,ip,j)+km(k-1,i,jp)        &
                  +km(k-1,ip,jp))
          enddo

          do k=1,n1-1
             s23(k,i,j)=-s23(k,i,j)*0.5*(km(k,i,j)+km(k,i,jp))
          enddo
       enddo
    enddo

  end subroutine diff_prep
  !
  ! ----------------------------------------------------------------------
  ! subroutine diff_upt: computes the diffusivity of velocities using
  ! a tri-diagnonal solver in the vertical at a u point, the deformation
  ! tensor component d31 is passed in via the tendency array
  !
  subroutine  diff_upt(n1,n2,n3,dn0,dzm,dzt,dxi,dyi,dt,sflx,tflx,sij,km, &
       u,w,tnd,flx)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: sij(n1,n2,n3)
    real, intent(in)    :: km(n1,n2,n3),u(n1,n2,n3),w(n1,n2,n3)
    real, intent(in)    :: sflx(n2,n3),tflx(n2,n3)
    real, intent(in)    :: dn0(n1),dzm(n1),dzt(n1),dxi,dyi,dt

    real, intent(inout) :: flx(n1)
    real, intent(out)   :: tnd(n1,n2,n3)

    dti   = 1.0/dt
    do k=1,n1
       sz7(k) = 0.
       sz8(k) = 0.
       flx(k) = 0.
       szx5(k,n2) = 0.
    end do

    do j=3,n3-2
       indh=0
       do i=3,n2-2
          indh=indh+1
          sz8(1)=dzm(1)*(km(1,i,j)+km(1,i+1,j))
          sz7(n1-1)  =.5*tflx(i,j)*(dn0(n1)+dn0(n1-1))
          sz7(1)     =.5*sflx(i,j)*(dn0(1)+dn0(2))
          do k=2,n1-1
             if (k < n1-1)  sz7(k)= (-(w(k,i+1,j)-w(k,i,j))*dxi) &
                                    *0.5*(km(k,i,j)+km(k,i+1,j))
             sz8(k)=dzm(k)*(km(k,i,j)+km(k,i+1,j))
             sxz4(indh,k)=u(k,i,j)*dn0(k) + dt*dzt(k)*(sz7(k-1)-sz7(k))
             sxz3(indh,k)=-0.5*dt*dzt(k)*sz8(k)
             sxz2(indh,k)=-0.5*dt*dzt(k)*sz8(k-1)
             sxz1(indh,k)=dn0(k)-sxz2(indh,k)-sxz3(indh,k)
          end do
          !
          ! Boundary conditions
          !
          if (vel_bc == 'noslip') then
             sxz1(indh,2)    = dn0(2)-2.*sxz2(indh,2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)-2.*sxz3(indh,n1-1)
          else
             sxz1(indh,2)    = dn0(2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)
          end if

          sxz3(indh,n1-1) = 0.
          sxz2(indh,2)    = 0.
          flx(1)=flx(1)+sz7(1)
          flx(n1-1)=flx(n1-1)+sz7(n1-1)

          do k=2,n1-1
             szx5(k,i) = (-2.*(u(k,i,j)-u(k,i-1,j))*dxi)*0.5*                 &
                  (km(k,i,j)+km(k-1,i,j))
          end do
       enddo

       do k=2,n1
          szx5(k,n2-1) = (-2.*(u(k,n2-1,j)-u(k,n2-2,j))*dxi)*0.5*             &
                  (km(k,n2-1,j)+km(k-1,n2-1,j))
       end do

       call tridiff(n2,n1-1,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
       !
       ! Back out diffusive tendencies from vertical implicit solution and
       ! solve for horizontal tendencies.  Accumulate vertical flux array
       !
       indh=0
       do i=3,n2-2
          indh=indh+1
          tnd(1,i,j) = 0.
          do k=2,n1-1
             tnd(k,i,j)=dti*(sxz1(indh,k)-u(k,i,j))-((szx5(k,i+1)-szx5(k,i))  &
                  *dxi + (sij(k,i,j)-sij(k,i,j-1))*dyi)/dn0(k)

             if (k < n1-1) flx(k)= flx(k)-dzm(k)*(km(k,i,j)+km(k,i+1,j))      &
                  *(sxz1(indh,k+1)-sxz1(indh,k))*.5
          enddo
          tnd(n1,i,j) = 0.
       enddo
    enddo

  end subroutine diff_upt
  !
  ! ----------------------------------------------------------------------
  ! subroutine diff_vpt: computes the diffusivity of velocities using
  ! a tri-diagnonal solver in the vertical at u or v pts depending on
  ! the values of ip and jp and the input arguments
  !
  subroutine  diff_vpt(n1,n2,n3,dn0,dzm,dzt,dxi,dyi,dt,sflx,tflx,sii,sij, &
       km,v,w,tnd,flx)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: dn0(n1),dzm(n1),dzt(n1),dxi,dyi,dt
    real, intent(in)    :: sflx(n2,n3),tflx(n2,n3)
    real, intent(in)    :: sii(n1,n2,n3),sij(n1,n2,n3)
    real, intent(in)    :: km(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)

    real, intent(inout) :: flx(n1)
    real, intent(out)   :: tnd(n1,n2,n3)

    dti = 1.0/dt
    do k=1,n1
       sz7(k)  = 0.
       sz8(k)  = 0.
       flx(k)  = 0.
    end do

    do j=3,n3-2
       indh  = 0
       do i=3,n2-2
          indh=indh+1
          sz8(1)=dzm(1)*(km(1,i,j)+km(1,i,j+1))
          sz7(n1-1)  =.5*(tflx(i,j))*(dn0(n1)+dn0(n1-1))
          sz7(1)     =.5*(sflx(i,j))*(dn0(1)+dn0(2))
          do k=2,n1-1
             if (k < n1-1) sz7(k) = (-(w(k,i,j+1)-w(k,i,j))*dyi) &
                                    *0.5*(km(k,i,j)+km(k,i,j+1))
             sz8(k)=dzm(k)*(km(k,i,j)+km(k,i,j+1))
             sxz4(indh,k)=v(k,i,j)*dn0(k) + dt*dzt(k)*(sz7(k-1)-sz7(k))
             sxz3(indh,k)=-0.5*dt*dzt(k)*sz8(k)
             sxz2(indh,k)=-0.5*dt*dzt(k)*sz8(k-1)
             sxz1(indh,k)=dn0(k)-sxz2(indh,k)-sxz3(indh,k)
          end do
          !
          ! Boundary conditions
          !
          if (vel_bc == 'noslip') then
             sxz1(indh,2)    = dn0(2)-2.*sxz2(indh,2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)-2.*sxz3(indh,n1-1)
          else
             sxz1(indh,2)    = dn0(2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)
          end if

          sxz3(indh,n1-1) = 0.
          sxz2(indh,2)    = 0.
          flx(1)=flx(1)+sz7(1)
          flx(n1-1)=flx(n1-1)+sz7(n1-1)
       enddo

       call tridiff(n2,n1-1,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
       !
       ! Back out diffusive tendencies from vertical implicit solution and
       ! solve for horizontal tendencies.  Accumulate vertical flux array
       !
       indh=0
       do i=3,n2-2
          indh=indh+1
          tnd(1,i,j) = 0.
          do k=2,n1-1
             tnd(k,i,j)=dti*(sxz1(indh,k)-v(k,i,j))-((sii(k,i,j+1)-sii(k,i,j))&
                  *dyi+(sij(k,i,j)-sij(k,i-1,j))*dxi)/dn0(k)
             if (k < n1-1) flx(k)= flx(k)-dzm(k)*(km(k,i,j)+km(k,i,j+1))      &
                  *(sxz1(indh,k+1)-sxz1(indh,k))*.5
          enddo
          tnd(n1,i,j) = 0.
       enddo
    enddo

  end subroutine diff_vpt
  !
  ! ----------------------------------------------------------------------
  ! subroutine diff_wpt: computes the diffusivity of velocities at a
  ! wpt
  !
  subroutine  diff_wpt(n1,n2,n3,dn0,dzm,dzt,dxi,dyi,dt,sflx,tflx,s23,km,w,u   &
       ,tnd,flx)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: s23(n1,n2,n3)
    real, intent(in)    :: km(n1,n2,n3),w(n1,n2,n3),u(n1,n2,n3)
    real, intent(in)    :: sflx(n2,n3),tflx(n2,n3)
    real, intent(in)    :: dn0(n1),dzm(n1),dzt(n1),dxi,dyi,dt

    real, intent(inout) :: flx(n1)
    real, intent(out)   :: tnd(n1,n2,n3)

    integer :: kp1,im1,jm1

    dti = 1.0/dt
    do k=1,n1
       sz8(k)  = 0.
       flx(k)  = 0.
    end do

    do j=1,n3
       do i=1,n2
          do k=1,n1
             tnd(k,i,j) = 0.
          end do
       end do
    end do

    sxz1(:,:) = 0.0

    do j=3,n3-2
       indh=0
       do i=3,n2-2
          indh=indh+1

          sz8(1)=dzt(2)*.5*(km(1,i,j)+km(2,i,j))
          do k=2,n1-2
             kp1 = k+1
             sz8(k)=dzt(kp1)*.5*(km(k,i,j)+km(kp1,i,j))
             sxz4(indh,k)=w(k,i,j)*(dn0(k)+dn0(kp1))*.5
             sxz3(indh,k)=-dt*dzm(k)*sz8(k)
             sxz2(indh,k)=-dt*dzm(k)*sz8(k-1)
             sxz1(indh,k)=(dn0(k)+dn0(kp1))*0.5 - sxz2(indh,k) - sxz3(indh,k)
          end do
          sxz2(indh,2)    = 0.
          sxz3(indh,n1-2) = 0.

          flx(1)=flx(1)+sflx(i,j)*dn0(2)
          flx(n1-1)=flx(n1-1)+tflx(i,j)*dn0(n1-1)

          do k=2,n1-1
             szx5(k,i) = ((u(k+1,i,j)-u(k,i,j))*dzm(k)+(w(k,i+1,j)-w(k,i,j))  &
                  *dxi)*(-0.5)*(km(k,i,j)+km(k,i+1,j))
          end do
          sxz4(indh,2)   =sxz4(indh,2) + dt*dzm(2)*sflx(i,j)*dn0(2)
          sxz4(indh,n1-2)=sxz4(indh,n1-2)-dt*dzm(n1-2)*tflx(i,j)*dn0(n1-1)
       end do

       do k=2,n1-1
          szx5(k,2) =  ((u(k+1,2,j)-u(k,2,j))*dzm(k) + (w(k,3,j)-w(k,2,j))    &
               *dxi)*(-0.5)*(km(k,2,j)+km(k,3,j))
       end do

       call tridiff(n2,n1-2,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
       !
       ! Back out diffusive tendencies from vertical implicit solution and
       ! solve for horizontal tendencies.  Accumulate vertical flux array
       !
       indh= 0
       do i=3,n2-2
          indh=indh+1
          im1 = max(i-1,2)
          jm1 = max(j-1,2)
          do k=2,n1-2
             tnd(k,i,j)=dti*(sxz1(indh,k)-w(k,i,j))-((szx5(k,i)-szx5(k,im1))  &
                  *dxi + (s23(k,i,j)-s23(k,i,jm1))*dyi)/((dn0(k)+dn0(k+1))*.5)
             flx(k) = flx(k)-dzt(k)*(km(k,i,j)+km(k+1,i,j))*0.5               &
                  *(sxz1(indh,k)-sxz1(indh,k-1))
          enddo
       enddo
    enddo

  end subroutine diff_wpt
  !
  ! -----------------------------------------------------------------------
  ! subroutine diffsclr: computes the diffusivity of a scalar using
  ! a tri-diagnonal solver in the vertical
  !
  subroutine diffsclr(n1,n2,n3,dtlt,dxi,dyi,dzm,dzt,dn0,sflx,tflx,scp,xkh,sct &
       ,flx)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: xkh(n1,n2,n3),scp(n1,n2,n3)
    real, intent(in)    :: sflx(n2,n3),tflx(n2,n3),dn0(n1)
    real, intent(in)    :: dxi,dyi,dzm(n1),dzt(n1),dtlt

    real, intent(out)   :: flx(n1,n2,n3),sct(n1,n2,n3)
    !
    ! compute vertical diffusion matrix coefficients for scalars,
    ! Coefficients need only be calculated once and can be used repeatedly
    ! for other scalars
    !
    dti       = 1.0/dtlt
    do k=1,n1
       sz7(k)   = 0.
    end do

    do j=3,n3-2
       do i=2,n2-2
          do k=2,n1-1
             szx1(k,i)=-(scp(k,i+1,j)-scp(k,i,j))*dxi*.25*(xkh(k,i,j)  +     &
                  xkh(k,i+1,j)+xkh(k-1,i,j)+xkh(k-1,i+1,j))
          enddo
       enddo
       !
       ! Set up Tri-diagonal Matrix
       !
       indh=0
       do i=3,n2-2
          indh=indh+1
          do k=2,n1-1
             if (k < n1-1) sz7(k)=dtlt*dzm(k)*xkh(k,i,j)
             sxz1(indh,k)=-dzt(k)*sz7(k-1)
             sxz2(indh,k)=-dzt(k)*sz7(k)
             sxz3(indh,k)=dn0(k)-sxz1(indh,k)-sxz2(indh,k)
             sxz4(indh,k)=scp(k,i,j)*dn0(k)
          enddo
          sxz4(indh,2)=scp(2,i,j)*dn0(2)                                     &
               + sflx(i,j)*(dn0(1)+dn0(2))     *.5 *dtlt*dzt(2)
          sxz4(indh,n1-1)=scp(n1-1,i,j)*dn0(n1-1)                            &
               - tflx(i,j)*(dn0(n1-1)+dn0(n1)) *.5*dtlt*dzt(n1-1)
       enddo

       call tridiff(n2,n1-1,indh,sxz1,sxz3,sxz2,sxz4,sxz5,sxz6)
       !
       ! compute scalar tendency in addition to vertical flux
       !
       indh=0
       do i=3,n2-2
          flx(1,i,j)   =sflx(i,j)*(dn0(1)+dn0(2))*.5
          flx(n1-1,i,j)=tflx(i,j)*(dn0(n1)+dn0(n1-1))*.5
          flx(n1,i,j)  =0.
          indh=indh+1
          do k=2,n1-1
             sct(k,i,j)=dti*(sxz5(indh,k)-scp(k,i,j))                       &
                  -((szx1(k,i)-szx1(k,i-1))                                   &
                  *dxi + (-(scp(k,i,j+1)-scp(k,i,j))*dyi*0.25*(xkh(k,i,j)     &
                  +xkh(k,i,j+1)+xkh(k-1,i,j)+xkh(k-1,i,j+1))+(scp(k,i,j)      &
                  -scp(k,i,j-1))*dyi*0.25*(xkh(k,i,j-1)+xkh(k,i,j)            &
                  +xkh(k-1,i,j-1)+xkh(k-1,i,j)))*dyi) /dn0(k)
             if (k<n1-1) flx(k,i,j)=-xkh(k,i,j)*(sxz5(indh,k+1)-sxz5(indh,k)) &
                  *dzm(k)
          end do
       enddo
    enddo

  end subroutine diffsclr
  !
  ! ---------------------------------------------------------------------
  ! Initializes sub-grid tke to a minimum value... also used to
  ! supress spuriously negative values
  !
  subroutine tkeinit(n1,tke)

    integer, intent(in) :: n1
    real, intent(inout) :: tke(n1)

    tke=max(tkemin,tke)

  end subroutine tkeinit

end module sgsm



