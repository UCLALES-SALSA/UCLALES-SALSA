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
MODULE sgsm

  USE stat, ONLY : sflg, updtst, acc_tend, sgsflxs, sgs_vel
  USE util, ONLY : tridiff, vel_bc
  IMPLICIT NONE
!
! setting the prandtl number to a value less than zero enforces an exponential
! decay in its value as a function of height from the surface with a scale
! height of 100m.  This is the default
!

  REAL, PARAMETER :: tkemin = 1.e-20
  REAL :: csx = 0.23
  REAL :: prndtl = 0.3333333333

  REAL, ALLOCATABLE, DIMENSION (:,:) :: sxy1, sxy2, sxy3, sxz1, sxz2, sxz3,    &
                                        sxz4, sxz5, sxz6  ,szx1, szx2, szx3, szx4, szx5
  REAL, ALLOCATABLE, DIMENSION (:)   :: sz1, sz2, sz3, sz4, sz5, sz6, sz7, sz8

  INTEGER :: k, i, j, indh, req(16)
  REAL    :: dti
  LOGICAL, SAVE :: Initialized = .FALSE.

CONTAINS

  SUBROUTINE diffuse_init(n1,n2,n3)

    INTEGER, INTENT (in) :: n1, n2, n3

    INTEGER :: nm

    nm = n1-1

    ALLOCATE(sxy1(n2,n3), sxy2(n2,n3), sxy3(n2,n3), sxz1(n2,nm), sxz2(n2,nm))
    ALLOCATE(sxz3(n2,nm), sxz4(n2,nm), sxz5(n2,nm), sxz6(n2,nm))
    ALLOCATE(szx1(n1,n2), szx2(n1,n2), szx3(n1,n2), szx4(n1,n2), szx5(n1,n2))
    ALLOCATE(sz1(n1),sz2(n1),sz3(n1),sz4(n1),sz5(n1),sz6(n1),sz7(n1),sz8(n1))

    initialized = .TRUE.
  END SUBROUTINE

  !
  ! ---------------------------------------------------------------------
  ! Subroutine DIFFUSE: Driver for calculating sub-grid fluxes (thus it
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
  SUBROUTINE diffuse

    USE grid, ONLY : a_up, a_uc, a_ut, a_vp, a_vc, a_vt, a_wp, a_wc, a_wt,    &
                     a_rv, a_rc, a_rp, a_ri, a_srp, a_tp, a_tt, a_sp, a_st, a_qt, a_qp, a_pexnr, a_theta,  &
                     a_temp, a_rsl, nscl, nxp, nyp,    &
                     nzp, zm, dxi, dyi, dzt, dzm, dtlt, dtlv , th00, dn0,  &
                     pi0, pi1, newsclr, level, isgstyp, uw_sfc, vw_sfc, ww_sfc, wt_sfc, &
                     wq_sfc

    USE util, ONLY          : get_avg3
    USE mpi_interface, ONLY : cyclics, cyclicc
    USE thrm, ONLY          : bruvais, fll_tkrs

    INTEGER :: n
    REAL    :: rx(nzp,nxp,nyp), rxt(nzp,nxp,nyp), a_tmp1(nzp,nxp,nyp), &
               a_tmp2(nzp,nxp,nyp), a_tmp3(nzp,nxp,nyp), a_tmp4(nzp,nxp,nyp), &
               a_tmp5(nzp,nxp,nyp), a_tmp6(nzp,nxp,nyp)

    SELECT CASE(level)
       CASE(1,2,3)
          rx = a_rv
          rxt = a_rp
       CASE(4,5)
          rx = a_rp
          rxt = a_rp + a_rc + a_ri + a_srp
    END SELECT


    IF (.NOT.Initialized) CALL diffuse_init(nzp, nxp, nyp)
    !
    ! ----------
    ! Calculate Deformation and stability for SGS calculations
    !
    CALL fll_tkrs(nzp,nxp,nyp,a_theta,a_pexnr,pi0,pi1,a_temp,rs=a_rsl)

    CALL bruvais(nzp,nxp,nyp,level,a_theta,a_tp,rxt,a_rsl,a_tmp3,dzm,th00)

    !
    ! the a_ut, a_wt, a_ut arrays are first used when the diffusive tendencies
    ! are calculated and applied.  Until then use them as scratch by
    ! associating them with scratch pointers (a-c)
    !
    CALL deform(nzp,nxp,nyp,dzm,dzt,dxi,dyi,a_up,a_vp,a_wp,a_tmp5,a_tmp6,     &
                a_tmp4,a_tmp2)

    ! ----------
    ! Calculate Eddy Viscosity/Diffusivity according to specified SGS model
    !
    SELECT CASE (isgstyp)
    CASE (1)
       CALL smagor(nzp,nxp,nyp,sflg,dxi,dn0,a_tmp3,a_tmp2,a_tmp1,zm)
    CASE (2)
       CALL deardf(nzp,nxp,nyp,sflg,dxi,zm,dn0,a_qp,a_qt,a_tmp3,a_tmp2,a_tmp1)

       CALL solv_tke(nzp,nxp,nyp,a_tmp3,a_tmp1,a_qp,a_qt,dn0,dzm,dzt,dxi,dyi,  &
                     dtlt)
    END SELECT
    !
    ! Diffuse momentum
    !
    IF (sflg) CALL acc_tend(nzp,nxp,nyp,a_uc,a_vc,a_wc,a_ut,a_vt,a_wt,         &
                            sz4,sz5,sz6,1,'sgs')

    CALL diff_prep(nzp,nxp,nyp,a_tmp5,a_tmp6,a_tmp4,a_tmp1)
    sxy1 = 0.; sxy2 = 0.

    CALL diff_vpt(nzp,nxp,nyp,dn0,dzm,dzt,dxi,dyi,dtlv,vw_sfc,sxy2,a_tmp6,     &
                  a_tmp5,a_tmp1,a_vp,a_wp,a_vt,sz2)

    CALL diff_upt(nzp,nxp,nyp,dn0,dzm,dzt,dxi,dyi,dtlv,uw_sfc,sxy1,a_tmp5,     &
                  a_tmp1,a_up,a_wp,a_ut,sz1)

    CALL diff_wpt(nzp,nxp,nyp,dn0,dzm,dzt,dyi,dxi,dtlv,ww_sfc,sxy1,a_tmp4,     &
                  a_tmp1,a_wp,a_up,a_wt,sz3)

    CALL cyclics(nzp,nxp,nyp,a_wt,req)
    CALL cyclicc(nzp,nxp,nyp,a_wt,req)
    CALL cyclics(nzp,nxp,nyp,a_vt,req)
    CALL cyclicc(nzp,nxp,nyp,a_vt,req)
    CALL cyclics(nzp,nxp,nyp,a_ut,req)
    CALL cyclicc(nzp,nxp,nyp,a_ut,req)

    IF (sflg) THEN
       CALL sgs_vel(nzp,nxp,nyp,sz1,sz2,sz3)
       CALL acc_tend(nzp,nxp,nyp,a_uc,a_vc,a_wc,a_ut,a_vt,a_wt,sz4,sz5,sz6,    &
                     2,'sgs')
    END IF
    !
    ! Diffuse scalars
    !
    a_tt = 0.
    DO n = 1, nscl
       CALL newsclr(n)
       sxy1 = 0.
       sxy2 = 0.
       IF ( associated(a_tp,a_sp) ) sxy1 = wt_sfc
       IF ( associated(a_tp,a_sp) ) sxy2 = wt_sfc
       IF ( associated(a_rp,a_sp) ) sxy1 = wq_sfc

       ! WHERE(abs(a_sp) < 1.e-40) a_sp=0. !stop denormal AZ
       IF (sflg) a_tmp1 = 0.
       IF ( isgstyp <= 1) THEN
          CALL diffsclr(nzp,nxp,nyp,dtlt,dxi,dyi,dzm,dzt,dn0,sxy1,sxy2,   &
                        a_sp,a_tmp2,a_st,a_tmp1)
       ELSE IF ( .NOT. associated(a_qp,a_sp) ) THEN
          CALL diffsclr(nzp,nxp,nyp,dtlt,dxi,dyi,dzm,dzt,dn0,sxy1,sxy2,   &
                        a_sp,a_tmp2,a_st,a_tmp1)
       END IF
       IF (sflg) THEN
          CALL get_avg3(nzp,nxp,nyp,a_tmp1,sz1)
          CALL updtst(nzp,'sgs',n,sz1,1)
          IF (associated(a_sp,a_tp))                                          &
             CALL sgsflxs(nzp,nxp,nyp,level,a_tmp3,rx,a_theta,a_tmp1,'tl')
          IF (associated(a_sp,a_rp))                          &
             CALL sgsflxs(nzp,nxp,nyp,level,a_tmp3,rx,a_theta,a_tmp1,'rt')
       END IF
       CALL cyclics(nzp,nxp,nyp,a_st,req)
       CALL cyclicc(nzp,nxp,nyp,a_st,req)
    END DO

  END SUBROUTINE diffuse
  !
  ! ---------------------------------------------------------------------
  ! Subroutine deform: computes the components of the deviatoric strain
  ! tensor, then computes s = du_i/dx_j(du_i/dx_j + du_j/dx_i) at a
  ! thermo point, the dummy arrays are respectively:
  !
  ! sxz1 -> div;
  ! szx2 -> s11;
  ! szx3 -> s33,
  ! szx4 -> s13
  !
  SUBROUTINE deform(n1,n2,n3,dzm,dzt,dx,dy,u,v,w,s12,s22,s23,s)

    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(in)    :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
    REAL, INTENT(in)    :: dzm(n1),dzt(n1),dx,dy

    REAL, INTENT(out)   :: s(n1,n2,n3), s22(n1,n2,n3)
    REAL, INTENT(out)   :: s12(n1,n2,n3),s23(n1,n2,n3)

    INTEGER :: ip,im,jp,jm,kp
    REAL    :: y1a,y2a,y3a,y1b,y2b,y3b
    REAL    :: s11_wpt,s22_wpt,s33_wpt,s12_wpt,s13_wpt,s23_wpt
    !
    ! calculate components of the stress tensor at their natural locations
    !
    DO j = 1, n3
       jm = max(j-1,1)
       jp = min(j+1,n3)
       DO i = 1, n2
          im = max(i-1,1)
          ip = min(i+1,n2)
          DO k = 1, n1
             szx2(k,i)  = 2.*(u(k,i,j)-u(k,im,j))*dx
             s22(k,i,j) = 2.*(v(k,i,j)-v(k,i,jm))*dy
             szx3(k,i)  = 2.*(w(k,i,j)-w(max(1,k-1),i,j))*dzt(k)
             s12(k,i,j) = (u(k,i,jp)-u(k,i,j))*dy + (v(k,ip,j)-v(k,i,j))*dx
          END DO
          szx3(1,i)  = 0.
          szx3(n1,i) = 0.

          DO k = 1, n1-1
             szx1(k,i)  = 0.333333*(szx2(k,i)+s22(k,i,j)+szx3(k,i))
             szx4(k,i)  = (u(k+1,i,j)-u(k,i,j))*dzm(k)+(w(k,ip,j)-w(k,i,j))*dx
             s23(k,i,j) = (v(k+1,i,j)-v(k,i,j))*dzm(k)+(w(k,i,jp)-w(k,i,j))*dy
          END DO
          szx1(n1,i) = 0.333333*(szx2(n1,i)+s22(n1,i,j)+szx3(n1,i))
       END DO
       !
       ! average to a w-point
       !
       DO i = 1, n2
          im = max(i-1,1)
          DO k = 1, n1-1
             kp = k+1
             y1a = (szx2(k,i)-szx1(k,i))
             y2a = (s22(k,i,j)-szx1(k,i))
             y3a = (szx3(k,i)-szx1(k,i))
             y1b = (szx2(kp,i)-szx1(kp,i))
             y2b = (s22(kp,i,j)-szx1(kp,i))
             y3b = (szx3(kp,i)-szx1(kp,i))

             s11_wpt = 0.5*(y1a*y1a+y1b*y1b)
             s22_wpt = 0.5*(y2a*y2a+y2b*y2b)
             s33_wpt = 0.5*(y3a*y3a+y3b*y3b)
             s12_wpt = 0.125*(s12(k,i,j)*s12(k,i,j) + s12(kp,i,j)*s12(kp,i,j)    &
                       +s12(k,im,j)*s12(k,im,j) + s12(k,i,jm)*s12(k,i,jm)        &
                       +s12(kp,im,j)*s12(kp,im,j)+s12(k,im,jm)*s12(k,im,jm)      &
                       +s12(kp,i,jm)*s12(kp,i,jm)+s12(kp,im,jm)*s12(kp,im,jm))
             s13_wpt = 0.5*(szx4(k,i)*szx4(k,i)+szx4(k,im)*szx4(k,im))
             s23_wpt = 0.5*(s23(k,i,j)*s23(k,i,j)+s23(k,i,jm)*s23(k,i,jm))
             s(k,i,j) = 0.5*(s11_wpt+s22_wpt+s33_wpt)                        &
                        + s12_wpt + s13_wpt + s23_wpt
          END DO
       END DO
    END DO

  END SUBROUTINE deform
  !
  ! ----------------------------------------------------------------------
  ! Subroutine smagor:  computes visc/diff, upon entering the routine
  ! kh is filled with s2 and ri is filled with N^2.  On statistical
  ! timsteps, SGS energy, dissipation, viscosity, diffusivity and
  ! lengthscales are stored.
  !
  SUBROUTINE smagor(n1,n2,n3,sflg,dxi,dn0,ri,kh,km,zm)

    USE defs, ONLY          : pi, vonk
    USE stat, ONLY          : tke_sgs
    USE util, ONLY          : get_avg3, get_cor3
    USE mpi_interface, ONLY : cyclics, cyclicc

    IMPLICIT NONE

    LOGICAL, INTENT(in) :: sflg
    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(in)    :: dxi,zm(n1),dn0(n1)
    REAL, INTENT(inout) :: ri(n1,n2,n3),kh(n1,n2,n3)
    REAL, INTENT(out)   :: km(n1,n2,n3)

    REAL :: delta,pr,yplus

    pr    = abs(prndtl)

    delta = 1./dxi
    delta = (zm(2)/dxi/dxi)**(1./3.)

    DO j = 3, n3-2
       DO i = 3, n2-2
          DO k = 1, n1-1
             yplus = max(zm(2),min(zm(k),(zm(n1-1)-zm(k))))
             ri(k,i,j) = max( -1., ri(k,i,j)/(kh(k,i,j) + 1.e-12) )
             km(k,i,j) = sqrt(max(0.,kh(k,i,j)*(1.-ri(k,i,j)/pr))) &
                         *0.5*(dn0(k)+dn0(k+1))/(1./(delta*csx)**2+1./(yplus*vonk)**2)
          END DO
          kh(n1,i,j) = kh(n1-1,i,j)
          km(n1,i,j) = km(n1-1,i,j)
       END DO
    END DO

    CALL cyclics(n1,n2,n3,km,req)
    CALL cyclicc(n1,n2,n3,km,req)


    IF (sflg) THEN
       CALL get_cor3(n1,n2,n3,km,km,sz1)
       CALL get_cor3(n1,n2,n3,km,kh,sz2)
       CALL updtst(n1,'sgs',-2,sz2,1)      ! dissipation

       DO k = 1, n1
          tke_sgs(k) = sz1(k)/(delta*pi*(csx*0.18))**2
          sz1(k) = 1./sqrt(1./(delta*csx)**2+1./(zm(k)*vonk+0.001)**2)
       END DO
       CALL updtst(n1,'sgs',-1,tke_sgs,1)  ! sgs tke
       CALL updtst(n1,'sgs',-5,sz1,1)      ! mixing length
       CALL updtst(n1,'sgs',-6,sz1,1)      ! dissipation lengthscale
    END IF

    DO j = 3, n3-2
       DO i = 3, n2-2
          DO k = 1, n1
            kh(k,i,j) = km(k,i,j)/pr
            IF (prndtl < 0.) THEN
               kh(k,i,j) = kh(k,i,j) * exp(zm(k)/(-100.))
            END IF
          END DO
       END DO
    END DO
    CALL cyclics(n1,n2,n3,kh,req)
    CALL cyclicc(n1,n2,n3,kh,req)

    IF (sflg) THEN
       CALL get_avg3(n1,n2,n3,km,sz3)
       CALL updtst(n1,'sgs',-3,sz3,1)  ! eddy viscosity
       CALL get_avg3(n1,n2,n3,kh,sz2)
       CALL updtst(n1,'sgs',-4,sz2,1)  ! eddy diffusivity
    END IF

  END SUBROUTINE smagor
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
  SUBROUTINE deardf(n1,n2,n3,sflg,dxi,zm,dn0,tke,tket,xx,zz,yy)

    USE defs, ONLY : vonk
    USE stat, ONLY : tke_sgs
    USE util, ONLY : get_avg3

    IMPLICIT NONE

    REAL, PARAMETER :: cm = 0.1, eps = 1.e-12, cs = 0.82, c_e1 = 0.225, c_e2 = 0.705

    LOGICAL, INTENT(in) :: sflg
    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(in)    :: dxi,dn0(n1),zm(n1)
    REAL, INTENT(in)    :: tke(n1,n2,n3)
    REAL, INTENT(inout) :: xx(n1,n2,n3),zz(n1,n2,n3)
    REAL, INTENT(out)   :: yy(n1,n2,n3),tket(n1,n2,n3)

    REAL :: ln,lm,ch
    !
    ! ----------
    ! Get neutral and effective lengthscales (ln,lm). Update tke tendency
    ! due to shear and buoyancy, fill:
    !   xx = dissipation lengthscale for later use
    !   yy = eddy viscosity (density weighted)
    !   zz = eddy diffusivity (density weighted)
    !
    DO k = 1, n1
       sz6(k) = 0.
       sz5(k) = 0.
       sz4(k) = 0.
    END DO

    DO j = 1, n3
       DO i = 1, n2
          DO k = 2, n1-2
             ln = sqrt(1./(dxi**2 + (0.23/(zm(k)*vonk))**2) )
             lm = min(cs*sqrt(tke(k,i,j)/max(eps,xx(k,i,j))),ln)
             ch = (1.+2.*lm/ln)
             tket(k,i,j) = cm*lm*sqrt(tke(k,i,j))*(zz(k,i,j) - ch*xx(k,i,j))
             xx(k,i,j) = lm/(c_e1 + c_e2*lm/ln)
             yy(k,i,j) = cm*lm*sqrt(tke(k,i,j))*.5*(dn0(k)+dn0(k+1))
             zz(k,i,j) = ch*yy(k,i,j)
             IF (sflg) THEN
                sz6(k) = sz6(k) + lm
                sz5(k) = sz5(k) + lm*ch
                sz4(k) = sz4(k) + (sqrt(tke(k,i,j)**3))/xx(k,i,j)
             END IF
          END DO
       END DO
    END DO
    !
    ! --- accumululate SGS model stats
    !
    IF (sflg) THEN
       CALL get_avg3(n1,n2,n3,tke,sz1)
       CALL get_avg3(n1,n2,n3,yy,sz2)
       CALL get_avg3(n1,n2,n3,zz,sz3)

       DO k = 1, n1
          sz6(k) = sz6(k)/REAL(n2*n3)
          sz5(k) = sz5(k)/REAL(n2*n3)
          sz4(k) = sz4(k)/REAL(n2*n3)
          tke_sgs(k) = sz1(k)
       END DO
       CALL updtst(n1,'sgs',-1,sz1,1)
       CALL updtst(n1,'sgs',-2,sz4,1)
       CALL updtst(n1,'sgs',-3,sz2,1)
       CALL updtst(n1,'sgs',-4,sz3,1)
       CALL updtst(n1,'sgs',-5,sz5,1)
       CALL updtst(n1,'sgs',-6,sz6,1)
    END IF

  END SUBROUTINE deardf
  !
  ! ----------------------------------------------------------------------
  ! Subroutine solv_tke: solves for dissipation and diffusion of tke,
  ! using an implicit solver for vertical diffusion and dissipation
  !
  SUBROUTINE  solv_tke(n1,n2,n3,le,km,ep,et,dn0,dzm,dzt,dxi,dyi,dt)

    IMPLICIT NONE

    INTEGER, INTENT (in) :: n1,n2,n3
    REAL, INTENT (in)    :: dn0(n1),dzm(n1),dzt(n1),dxi,dyi,dt
    REAL, INTENT (in)    :: km(n1,n2,n3),ep(n1,n2,n3),le(n1,n2,n3)
    REAL, INTENT (inout) :: et(n1,n2,n3)

    !
    ! Add horizontal diffusivities to tendency arrays
    !
    DO j = 2, n3-1
       DO i = 2, n2-1
          DO k = 2, n1-2
             et(k,i,j) = et(k,i,j)+(                                           &
                         ((km(k,i+1,j)+km(k,i,j))*(ep(k,i+1,j)-ep(k,i,j))           &
                         -(km(k,i,j)+km(k,i-1,j))*(ep(k,i,j)-ep(k,i-1,j)))*dxi*dxi  &
                         +((km(k,i,j+1)+km(k,i,j))*(ep(k,i,j+1)-ep(k,i,j))          &
                         -(km(k,i,j)+km(k,i,j-1))*(ep(k,i,j)-ep(k,i,j-1)))*dyi*dyi  &
                          )/((dn0(k)+dn0(k+1))*.5)
          END DO
       END DO
    END DO
    !
    ! Set up and solve implicit system
    !
    dti = 1.0/dt
    DO k = 1, n1
       sz7(k) = 0.
    END DO

    DO j = 2, n3-1
       indh = 0
       DO i = 2, n2-1
          indh = indh+1
          sz7(1) = dzt(2)*(km(2,i,j)+km(1,i,j))
          DO k = 2, n1-2
             sz7(k)       = dzt(k+1)*(km(k+1,i,j)+km(k,i,j))
             sxz4(indh,k) = 0.5*(dn0(k)+dn0(k+1))*(ep(k,i,j)+dt*et(k,i,j))
             sxz3(indh,k) = -dt*dzm(k)*sz7(k)
             sxz2(indh,k) = -dt*dzm(k)*sz7(k-1)
             sxz1(indh,k) = (dn0(k)+dn0(k+1))*0.5                             &
                            *(1.+dt*sqrt(ep(k,i,j))/le(k,i,j))-sxz2(indh,k)-sxz3(indh,k)

          END DO
       END DO

       CALL tridiff(n2,n1-2,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
       !
       ! Back out tendencies from implicit solver
       !
       indh = 0
       DO i = 2, n2-1
          indh = indh+1
          DO k = 2, n1-2
             et(k,i,j) = dti*(sxz1(indh,k)-ep(k,i,j))
          END DO
       END DO
    END DO

  END SUBROUTINE solv_tke
  !
  ! ----------------------------------------------------------------------
  ! Subroutine diff_prep: multiplies the strain components computed in
  ! "deform" by the appropriately averaged value of km, in preperation
  ! for use by the diffusion routines
  !
  SUBROUTINE  diff_prep(n1,n2,n3,s12,s22,s23,km)

    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(inout) :: s22(n1,n2,n3)
    REAL, INTENT(inout) :: s12(n1,n2,n3),s23(n1,n2,n3)
    REAL, INTENT(inout) :: km(n1,n2,n3)

    INTEGER :: ip,jp

    DO j = 2, n3-1
       jp = min(j+1,n3)
       DO i = 2, n2-1
          ip = min(i+1,n2)
          DO k = 2, n1-1
             s22(k,i,j) = -s22(k,i,j)*0.5*(km(k,i,j)+km(k-1,i,j))
             s12(k,i,j) = -s12(k,i,j)*0.125*(km(k,i,j)+km(k,ip,j)+km(k,i,jp)   &
                          +km(k,ip,jp)+km(k-1,i,j)+km(k-1,ip,j)+km(k-1,i,jp)   &
                          +km(k-1,ip,jp))
          END DO

          DO k = 1, n1-1
             s23(k,i,j) = -s23(k,i,j)*0.5*(km(k,i,j)+km(k,i,jp))
          END DO
       END DO
    END DO

  END SUBROUTINE diff_prep
  !
  ! ----------------------------------------------------------------------
  ! Subroutine diff_upt: computes the diffusivity of velocities using
  ! a tri-diagnonal solver in the vertical at a u point, the deformation
  ! tensor component d31 is passed in via the tendency array
  !
  SUBROUTINE  diff_upt(n1,n2,n3,dn0,dzm,dzt,dxi,dyi,dt,sflx,tflx,sij,km, &
                       u,w,tnd,flx)

    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(in)    :: sij(n1,n2,n3)
    REAL, INTENT(in)    :: km(n1,n2,n3),u(n1,n2,n3),w(n1,n2,n3)
    REAL, INTENT(in)    :: sflx(n2,n3),tflx(n2,n3)
    REAL, INTENT(in)    :: dn0(n1),dzm(n1),dzt(n1),dxi,dyi,dt

    REAL, INTENT(inout) :: flx(n1)
    REAL, INTENT(out)   :: tnd(n1,n2,n3)

    dti = 1.0/dt
    DO k = 1, n1
       sz7(k) = 0.
       sz8(k) = 0.
       flx(k) = 0.
       szx5(k,n2) = 0.
    END DO

    DO j = 3, n3-2
       indh = 0
       DO i = 3, n2-2
          indh = indh+1
          sz8(1) = dzm(1)*(km(1,i,j)+km(1,i+1,j))
          sz7(n1-1) =.5*tflx(i,j)*(dn0(n1)+dn0(n1-1))
          sz7(1)    =.5*sflx(i,j)*(dn0(1)+dn0(2))
          DO k = 2, n1-1
             IF (k < n1-1)  sz7(k) = (-(w(k,i+1,j)-w(k,i,j))*dxi) &
                                     *0.5*(km(k,i,j)+km(k,i+1,j))
             sz8(k) = dzm(k)*(km(k,i,j)+km(k,i+1,j))
             sxz4(indh,k) = u(k,i,j)*dn0(k) + dt*dzt(k)*(sz7(k-1)-sz7(k))
             sxz3(indh,k) = -0.5*dt*dzt(k)*sz8(k)
             sxz2(indh,k) = -0.5*dt*dzt(k)*sz8(k-1)
             sxz1(indh,k) = dn0(k)-sxz2(indh,k)-sxz3(indh,k)
          END DO
          !
          ! Boundary conditions
          !
          IF (vel_bc == 'noslip') THEN
             sxz1(indh,2)    = dn0(2)-2.*sxz2(indh,2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)-2.*sxz3(indh,n1-1)
          ELSE
             sxz1(indh,2)    = dn0(2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)
          END IF

          sxz3(indh,n1-1) = 0.
          sxz2(indh,2)    = 0.
          flx(1) = flx(1)+sz7(1)
          flx(n1-1) = flx(n1-1)+sz7(n1-1)

          DO k = 2, n1-1
             szx5(k,i) = (-2.*(u(k,i,j)-u(k,i-1,j))*dxi)*0.5*                 &
                         (km(k,i,j)+km(k-1,i,j))
          END DO
       END DO

       DO k = 2, n1
          szx5(k,n2-1) = (-2.*(u(k,n2-1,j)-u(k,n2-2,j))*dxi)*0.5*             &
                         (km(k,n2-1,j)+km(k-1,n2-1,j))
       END DO

       CALL tridiff(n2,n1-1,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
       !
       ! Back out diffusive tendencies from vertical implicit solution and
       ! solve for horizontal tendencies.  Accumulate vertical flux array
       !
       indh = 0
       DO i = 3, n2-2
          indh = indh+1
          tnd(1,i,j) = 0.
          DO k = 2, n1-1
             tnd(k,i,j) = dti*(sxz1(indh,k)-u(k,i,j))-((szx5(k,i+1)-szx5(k,i))  &
                          *dxi + (sij(k,i,j)-sij(k,i,j-1))*dyi)/dn0(k)

             IF (k < n1-1) flx(k) = flx(k)-dzm(k)*(km(k,i,j)+km(k,i+1,j))      &
                                    *(sxz1(indh,k+1)-sxz1(indh,k))*.5
          END DO
          tnd(n1,i,j) = 0.
       END DO
    END DO

  END SUBROUTINE diff_upt
  !
  ! ----------------------------------------------------------------------
  ! Subroutine diff_vpt: computes the diffusivity of velocities using
  ! a tri-diagnonal solver in the vertical at u or v pts depending on
  ! the values of ip and jp and the input arguments
  !
  SUBROUTINE  diff_vpt(n1,n2,n3,dn0,dzm,dzt,dxi,dyi,dt,sflx,tflx,sii,sij, &
                       km,v,w,tnd,flx)

    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(in)    :: dn0(n1),dzm(n1),dzt(n1),dxi,dyi,dt
    REAL, INTENT(in)    :: sflx(n2,n3),tflx(n2,n3)
    REAL, INTENT(in)    :: sii(n1,n2,n3),sij(n1,n2,n3)
    REAL, INTENT(in)    :: km(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)

    REAL, INTENT(inout) :: flx(n1)
    REAL, INTENT(out)   :: tnd(n1,n2,n3)

    dti = 1.0/dt
    DO k = 1, n1
       sz7(k) = 0.
       sz8(k) = 0.
       flx(k) = 0.
    END DO

    DO j = 3, n3-2
       indh = 0
       DO i = 3, n2-2
          indh = indh+1
          sz8(1) = dzm(1)*(km(1,i,j)+km(1,i,j+1))
          sz7(n1-1) =.5*(tflx(i,j))*(dn0(n1)+dn0(n1-1))
          sz7(1)    =.5*(sflx(i,j))*(dn0(1)+dn0(2))
          DO k = 2, n1-1
             IF (k < n1-1) sz7(k) = (-(w(k,i,j+1)-w(k,i,j))*dyi) &
                                    *0.5*(km(k,i,j)+km(k,i,j+1))
             sz8(k) = dzm(k)*(km(k,i,j)+km(k,i,j+1))
             sxz4(indh,k) = v(k,i,j)*dn0(k) + dt*dzt(k)*(sz7(k-1)-sz7(k))
             sxz3(indh,k) = -0.5*dt*dzt(k)*sz8(k)
             sxz2(indh,k) = -0.5*dt*dzt(k)*sz8(k-1)
             sxz1(indh,k) = dn0(k)-sxz2(indh,k)-sxz3(indh,k)
          END DO
          !
          ! Boundary conditions
          !
          IF (vel_bc == 'noslip') THEN
             sxz1(indh,2)    = dn0(2)-2.*sxz2(indh,2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)-2.*sxz3(indh,n1-1)
          ELSE
             sxz1(indh,2)    = dn0(2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)
          END IF

          sxz3(indh,n1-1) = 0.
          sxz2(indh,2)    = 0.
          flx(1) = flx(1)+sz7(1)
          flx(n1-1) = flx(n1-1)+sz7(n1-1)
       END DO

       CALL tridiff(n2,n1-1,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
       !
       ! Back out diffusive tendencies from vertical implicit solution and
       ! solve for horizontal tendencies.  Accumulate vertical flux array
       !
       indh = 0
       DO i = 3, n2-2
          indh = indh+1
          tnd(1,i,j) = 0.
          DO k = 2, n1-1
             tnd(k,i,j) = dti*(sxz1(indh,k)-v(k,i,j))-((sii(k,i,j+1)-sii(k,i,j))&
                          *dyi+(sij(k,i,j)-sij(k,i-1,j))*dxi)/dn0(k)
             IF (k < n1-1) flx(k) = flx(k)-dzm(k)*(km(k,i,j)+km(k,i,j+1))      &
                                    *(sxz1(indh,k+1)-sxz1(indh,k))*.5
          END DO
          tnd(n1,i,j) = 0.
       END DO
    END DO

  END SUBROUTINE diff_vpt
  !
  ! ----------------------------------------------------------------------
  ! Subroutine diff_wpt: computes the diffusivity of velocities at a
  ! wpt
  !
  SUBROUTINE  diff_wpt(n1,n2,n3,dn0,dzm,dzt,dxi,dyi,dt,sflx,tflx,s23,km,w,u,   &
                       tnd,flx)

    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(in)    :: s23(n1,n2,n3)
    REAL, INTENT(in)    :: km(n1,n2,n3),w(n1,n2,n3),u(n1,n2,n3)
    REAL, INTENT(in)    :: sflx(n2,n3),tflx(n2,n3)
    REAL, INTENT(in)    :: dn0(n1),dzm(n1),dzt(n1),dxi,dyi,dt

    REAL, INTENT(inout) :: flx(n1)
    REAL, INTENT(out)   :: tnd(n1,n2,n3)

    INTEGER :: kp1,im1,jm1

    dti = 1.0/dt
    DO k = 1, n1
       sz8(k) = 0.
       flx(k) = 0.
    END DO

    DO j = 1, n3
       DO i = 1, n2
          DO k = 1, n1
             tnd(k,i,j) = 0.
          END DO
       END DO
    END DO

    sxz1(:,:) = 0.0

    DO j = 3, n3-2
       indh = 0
       DO i = 3, n2-2
          indh = indh+1

          sz8(1) = dzt(2)*.5*(km(1,i,j)+km(2,i,j))
          DO k = 2, n1-2
             kp1 = k+1
             sz8(k) = dzt(kp1)*.5*(km(k,i,j)+km(kp1,i,j))
             sxz4(indh,k) = w(k,i,j)*(dn0(k)+dn0(kp1))*.5
             sxz3(indh,k) = -dt*dzm(k)*sz8(k)
             sxz2(indh,k) = -dt*dzm(k)*sz8(k-1)
             sxz1(indh,k) = (dn0(k)+dn0(kp1))*0.5 - sxz2(indh,k) - sxz3(indh,k)
          END DO
          sxz2(indh,2)    = 0.
          sxz3(indh,n1-2) = 0.

          flx(1) = flx(1)+sflx(i,j)*dn0(2)
          flx(n1-1) = flx(n1-1)+tflx(i,j)*dn0(n1-1)

          DO k = 2, n1-1
             szx5(k,i) = ((u(k+1,i,j)-u(k,i,j))*dzm(k)+(w(k,i+1,j)-w(k,i,j))  &
                         *dxi)*(-0.5)*(km(k,i,j)+km(k,i+1,j))
          END DO
          sxz4(indh,2)    = sxz4(indh,2) + dt*dzm(2)*sflx(i,j)*dn0(2)
          sxz4(indh,n1-2) = sxz4(indh,n1-2)-dt*dzm(n1-2)*tflx(i,j)*dn0(n1-1)
       END DO

       DO k = 2, n1-1
          szx5(k,2) = ((u(k+1,2,j)-u(k,2,j))*dzm(k) + (w(k,3,j)-w(k,2,j))    &
                      *dxi)*(-0.5)*(km(k,2,j)+km(k,3,j))
       END DO

       CALL tridiff(n2,n1-2,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
       !
       ! Back out diffusive tendencies from vertical implicit solution and
       ! solve for horizontal tendencies.  Accumulate vertical flux array
       !
       indh = 0
       DO i = 3, n2-2
          indh = indh+1
          im1  = max(i-1,2)
          jm1  = max(j-1,2)
          DO k = 2, n1-2
             tnd(k,i,j) = dti*(sxz1(indh,k)-w(k,i,j))-((szx5(k,i)-szx5(k,im1))  &
                          *dxi + (s23(k,i,j)-s23(k,i,jm1))*dyi)/((dn0(k)+dn0(k+1))*.5)
             flx(k) = flx(k)-dzt(k)*(km(k,i,j)+km(k+1,i,j))*0.5               &
                      *(sxz1(indh,k)-sxz1(indh,k-1))
          END DO
       END DO
    END DO

  END SUBROUTINE diff_wpt
  !
  ! -----------------------------------------------------------------------
  ! Subroutine diffsclr: computes the diffusivity of a scalar using
  ! a tri-diagnonal solver in the vertical
  !
  SUBROUTINE diffsclr(n1,n2,n3,dtlt,dxi,dyi,dzm,dzt,dn0,sflx,tflx,scp,xkh,sct, &
                      flx)

    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(in)    :: xkh(n1,n2,n3),scp(n1,n2,n3)
    REAL, INTENT(in)    :: sflx(n2,n3),tflx(n2,n3),dn0(n1)
    REAL, INTENT(in)    :: dxi,dyi,dzm(n1),dzt(n1),dtlt

    REAL, INTENT(out)   :: flx(n1,n2,n3),sct(n1,n2,n3)
    !
    ! compute vertical diffusion matrix coefficients for scalars,
    ! Coefficients need only be calculated once and can be used repeatedly
    ! for other scalars
    !
    dti = 1.0/dtlt
    DO k = 1, n1
       sz7(k)   = 0.
    END DO

    DO j = 3, n3-2
       DO i = 2, n2-2
          DO k = 2, n1-1

             szx1(k,i) = -(scp(k,i+1,j)-scp(k,i,j))*dxi*.25*(xkh(k,i,j)  +     &
                         xkh(k,i+1,j)+xkh(k-1,i,j)+xkh(k-1,i+1,j))
          END DO
       END DO
       !
       ! Set up Tri-diagonal Matrix
       !
       indh = 0
       DO i = 3, n2-2
          indh = indh+1
          DO k = 2, n1-1
             IF (k < n1-1) sz7(k) = dtlt*dzm(k)*xkh(k,i,j)
             sxz1(indh,k) = -dzt(k)*sz7(k-1)
             sxz2(indh,k) = -dzt(k)*sz7(k)
             sxz3(indh,k) = dn0(k)-sxz1(indh,k)-sxz2(indh,k)
             sxz4(indh,k) = scp(k,i,j)*dn0(k)
          END DO
          sxz4(indh,2) = scp(2,i,j)*dn0(2)                                     &
                         + sflx(i,j)*(dn0(1)+dn0(2))*.5   *dtlt*dzt(2)
          sxz4(indh,n1-1) = scp(n1-1,i,j)*dn0(n1-1)                            &
                            - tflx(i,j)*(dn0(n1-1)+dn0(n1))*.5   *dtlt*dzt(n1-1)
       END DO
       CALL tridiff(n2,n1-1,indh,sxz1,sxz3,sxz2,sxz4,sxz5,sxz6)
       !
       ! compute scalar tendency in addition to vertical flux
       !
       indh = 0
       DO i = 3, n2-2
          flx(1,i,j)    = sflx(i,j)*(dn0(1)+dn0(2))*.5
          flx(n1-1,i,j) = tflx(i,j)*(dn0(n1)+dn0(n1-1))*.5
          flx(n1,i,j)   = 0.
          indh = indh+1
          DO k = 2, n1-1
             sct(k,i,j) = dti*(sxz5(indh,k)-scp(k,i,j))                               &
                          -((szx1(k,i)-szx1(k,i-1))                                   &
                          *dxi + (-(scp(k,i,j+1)-scp(k,i,j))*dyi*0.25*(xkh(k,i,j)     &
                          +xkh(k,i,j+1)+xkh(k-1,i,j)+xkh(k-1,i,j+1))+(scp(k,i,j)      &
                          -scp(k,i,j-1))*dyi*0.25*(xkh(k,i,j-1)+xkh(k,i,j)            &
                          +xkh(k-1,i,j-1)+xkh(k-1,i,j)))*dyi) /dn0(k)
             IF (k < n1-1) flx(k,i,j) = -xkh(k,i,j)*(sxz5(indh,k+1)-sxz5(indh,k)) &
                                      *dzm(k)
          END DO
       END DO
    END DO

  END SUBROUTINE diffsclr
  !
  ! ---------------------------------------------------------------------
  ! Initializes sub-grid tke to a minimum value... also used to
  ! supress spuriously negative values
  !
  SUBROUTINE tkeinit(n1,tke)

    INTEGER, INTENT(in) :: n1
    REAL, INTENT(inout) :: tke(n1)

    tke = max(tkemin,tke)

  END SUBROUTINE tkeinit

END MODULE sgsm



