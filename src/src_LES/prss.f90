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
! PRSS: Pressure Solver:  Solves the anelastic or bousinessq system 
! for the pressure using a fractional step method, which is implemented
! using fft's and a tri-diagonal solver in the vertical
!
MODULE prss

   USE defs, ONLY : pi
   IMPLICIT NONE

CONTAINS
   !
   !----------------------------------------------------------------------
   ! Subroutine poisson: called by timesteping driver to invert the
   ! poisson equation for pressure and apply the velocity tendencies.
   !
   SUBROUTINE poisson

      USE grid, ONLY : nxp, nyp, nzp, dtlt, dxi, dyi, dzm, dzt, a_up,                     &
                       a_uc, a_ut, a_vp, a_vc, a_vt, a_wp, a_wc, a_wt, a_press, a_pexnr,  &
                       th00, dn0, wsavex, wsavey
      USE stat, ONLY : fill_scalar, sflg
      USE util, ONLY : ae1mm

      COMPLEX, ALLOCATABLE :: s1(:,:,:)
      REAL    :: mxdiv, awpbar(nzp)
      INTEGER :: ix,iy

      ix = max(1,nxp-4)
      iy = max(1,nyp-4)

      ALLOCATE (s1(ix,iy,nzp))
      s1 = 0.0
      !
      ! -------
      ! Do first step of asselin filter, first saving corrlations of
      ! tendencies for TKE budget on statistics timestep
      !

      CALL asselin(1)
      CALL apl_tnd(nzp,nxp,nyp,a_up,a_vp,a_wp,a_ut,a_vt,a_wt,dtlt)
      !
      ! ------
      ! Pressure Solve
      !
      CALL poiss(nzp,nxp,nyp,ix,iy,a_up,a_vp,a_wp,a_pexnr,a_press,dn0,th00,dzt, &
                 dzm,dxi,dyi,dtlt,s1,wsavex,wsavey)
      CALL ae1mm(nzp,nxp,nyp,a_wp,awpbar)
      !
      ! -------
      ! Do second step of asselin filter, first saving corrlations of
      ! tendencies for TKE budget on statistics timestep, note that the
      ! old centered velocity resides in up,vp,wp after the second step
      ! of the asselin filter, hence the pressure correlation terms in the
      ! tke budget include the effects of time filtering
      !
      CALL asselin(2)
      CALL velocity_bcs
      CALL get_diverg(nzp,nxp,nyp,ix,iy,s1,a_up,a_vp,a_wp,dn0,dzt,dxi,dyi,  &
                      dtlt,mxdiv)

      IF (sflg) THEN
         CALL get_diverg(nzp,nxp,nyp,ix,iy,s1,a_up,a_vp,a_wp,dn0,dzt,dxi,dyi,  &
                         dtlt,mxdiv)
         CALL fill_scalar(2,mxdiv)
         CALL prs_cor(nzp,nxp,nyp,a_pexnr,a_up,a_vp,a_wp,dzm,dxi,dyi,th00)
         CALL chk_tsplt(nzp,nxp,nyp,a_up,a_vp,a_wp,a_uc,a_vc,a_wc)
      END IF
      DEALLOCATE (s1)

   END SUBROUTINE poisson

   !
   ! --------------------------------------------------------------------
   ! Subroutine apl_tnd: applies tendencies to velocity field
   !
   SUBROUTINE apl_tnd(n1,n2,n3,u,v,w,ut,vt,wt,dtlt)

      USE util, ONLY : velset

      INTEGER :: n1,n2,n3,i,k,j
      REAL    :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
      REAL    :: ut(n1,n2,n3),vt(n1,n2,n3),wt(n1,n2,n3),dtlt,dt

      dt = dtlt*2.

      DO j = 1, n3
         DO i = 1, n2
            DO k = 2, n1-1
               w(k,i,j) = w(k,i,j)+wt(k,i,j)*dt
               u(k,i,j) = u(k,i,j)+ut(k,i,j)*dt
               v(k,i,j) = v(k,i,j)+vt(k,i,j)*dt
            END DO
         END DO
      END DO

      CALL velset(n1,n2,n3,u,v,w)

   END SUBROUTINE apl_tnd
   !
   ! --------------------------------------------------------------------
   ! Subroutine poiss: called each timestep to evaluate the pressure
   ! in accordance with the anelastic continuity equation, and then apply
   ! the pressure to the velocity terms for three dimensional flow,
   ! cyclic in x and y.  pp and pc are used as scratch arrays in the
   ! call to trdprs.  pp is filled with its diagnostic value in fll_prs
   !
   SUBROUTINE poiss(n1,n2,n3,ix,iy,u,v,w,pp,pc,dn0,th00,dzt,dzm,dx,dy,  &
                    dtlt,s1,wsvx,wsvy)

      USE util, ONLY  : velset, get_fft_twodim

      INTEGER :: n1,n2,n3,ix,iy
      REAL    :: pp(n1,n2,n3),pc(n1,n2,n3),dmy
      REAL    :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
      REAL    :: wsvx(1:),wsvy(1:),dn0(n1),dzt(n1),dzm(n1),dx,dy,dtlt,th00
      COMPLEX :: s1(ix,iy,n1)

      CALL get_diverg(n1,n2,n3,ix,iy,s1,u,v,w,dn0,dzt,dx,dy,dtlt,dmy)
      CALL get_fft_twodim(ix,iy,n1,s1,wsvx,wsvy,-1)
      CALL trdprs(n1,ix,iy,s1,dn0,dzt,dzm,dx,dy)
      CALL get_fft_twodim(ix,iy,n1,s1,wsvx,wsvy,+1)

      CALL fll_prs(n1,n2,n3,ix,iy,pp,s1)
      CALL prs_grd(n1,n2,n3,pp,u,v,w,dzm,dx,dy,dtlt)

      CALL velset(n1,n2,n3,u,v,w)

      pp(:,:,:) = pp(:,:,:)/th00
      pc(:,:,:) = pp(:,:,:)

   END SUBROUTINE poiss
   !
   ! --------------------------------------------------------------------
   ! Subroutine get_diverg: gets velocity divergence and puts it into
   ! a complexvalue array for use in pressure calculation
   !
   SUBROUTINE get_diverg(n1,n2,n3,ix,iy,s1,u,v,w,dn0,dz,dx,dy,dt,mxdiv)

      INTEGER, INTENT (in)  :: n1,n2,n3,ix,iy
      REAL, INTENT (in)     :: dz(n1),dn0(n1),dx,dy,dt
      REAL, INTENT (in)     :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
      REAL, INTENT (out)    :: mxdiv
      COMPLEX, INTENT (out) :: s1(ix,iy,n1)

      INTEGER :: k,i,j,l,m
      REAL    :: xf,yf,zf,wf1,wf2,dtv,dti

      s1(:,:,:) = (0.0,0.0)
      dtv = dt*2.
      dti = 1./dtv
      m = 0
      DO j = 3, n3-2
         m = m+1
         l = 0
         DO i = 3, n2-2
            l = l+1
            DO k = 2, n1-1
               wf1 = 0.5*(dn0(k+1)+dn0(k))
               wf2 = 0.5*(dn0(k)+dn0(k-1))
               IF (k == 2 )   wf2 = 0.
               IF (k == n1-1) wf1 = 0.
               xf = dn0(k)*dx*dti
               yf = dn0(k)*dy*dti
               zf = dti*dz(k)
               s1(l,m,k) = (wf1*w(k,i,j)-wf2*w(k-1,i,j))*zf                  &
                           +(v(k,i,j)-v(k,i,j-1))*yf+(u(k,i,j)-u(k,i-1,j))*xf
            END DO
         END DO
      END DO
      !
      ! save mxdiv to a statistics array, no reduction necessary as this is done
      ! in post processing
      !
      mxdiv = maxval(REAL(s1))

   END SUBROUTINE get_diverg
   !
   !----------------------------------------------------------------------
   ! Subroutine fll_prs: writes the pressure to the appropriate array
   !
   SUBROUTINE fll_prs(n1,n2,n3,ix,iy,pp,s1)

      USE mpi_interface, ONLY : cyclics,cyclicc

      INTEGER :: n1,n2,n3,ix,iy,k,i,j,l,m,req(16)
      REAL    :: pp(n1,n2,n3)
      COMPLEX :: s1(ix,iy,n1)

      pp(:,:,:) = 0.0
      DO k = 2, n1-1
         l = 0
         DO i = 3, n2-2
            l = l+1
            m = 0
            DO j = 3, n3-2
               m = m+1
               pp(k,i,j) = REAL(s1(l,m,k))
            END DO
         END DO
      END DO
      CALL cyclics(n1,n2,n3,pp,req)
      CALL cyclicc(n1,n2,n3,pp,req)

   END SUBROUTINE fll_prs
   !
   !---------------------------------------------------------------------
   ! TRDPRS: solves for the wave number (l,m) component of
   ! pressure in a vertical column using a tri-diagonal solver.
   !
   SUBROUTINE trdprs(n1,ix,iy,s1,dn0,dzt,dzm,dx,dy)

      USE mpi_interface, ONLY : yoffset, nypg, xoffset, wrxid, wryid, nxpg
      USE util, ONLY          : tridiff

      INTEGER, INTENT (in)    :: n1,ix,iy
      REAL, INTENT (in)       :: dn0(n1),dzt(n1),dzm(n1),dx,dy
      COMPLEX, INTENT (inout) :: s1(ix,iy,n1)

      REAL    :: ak(ix,n1),dk(ix,n1),bk(ix,n1),ck(ix,n1)
      REAL    :: xk(ix,n1),yk(ix,n1),wv(ix,iy)

      INTEGER :: k,l,m
      REAL    :: fctl,fctm,xl,xm,af,cf
      fctl = 2.*pi/float(nxpg-4)
      fctm = 2.*pi/float(nypg-4)

      DO l = 1, ix
         IF(l+xoffset(wrxid) <= (nxpg-4)/2+1) THEN
            xl = float(l-1+xoffset(wrxid))
         ELSE
            xl = float(l-(nxpg-4)-1+xoffset(wrxid))
         END IF
     
         DO m = 1, iy
            IF(m+yoffset(wryid) <= (nypg-4)/2+1) THEN
               xm = float(m-1+yoffset(wryid))
            ELSE
               xm = float(m-(nypg-4)-1+yoffset(wryid))
            END IF
            wv(l,m) = 2.*((cos(fctl*xl)-1.)*dx*dx + (cos(fctm*xm)-1.)*dy*dy)
         END DO
      END DO

      IF(wrxid == 0 .AND. wryid == 0 ) THEN
         wv(1,1) = 0.
      END IF
      !
      ! configure vectors for tri-diagonal solver
      !
      DO m = 1, iy
         DO k = 2, n1-1
            af = (dn0(k)+dn0(k-1))*.5
            cf = (dn0(k+1)+dn0(k))*.5
            IF (k == 2   ) af = 0.
            IF (k == n1-1) cf =0.
            DO l = 1, ix
               ak(l,k) = dzt(k)*dzm(k-1)*af
               bk(l,k) = s1(l,m,k)
               ck(l,k) = dzt(k)*dzm(k)*cf
               dk(l,k) = dn0(k)*wv(l,m)-(ak(l,k)+ck(l,k))
            END DO
         END DO
         !
         ! solve for fourier components, x_k, given a tri-diagonal matrix of the
         ! form a_k x_k-1 + d_k x_k + c_k x_k+1 = b_k.  y_k is a scratch array.
         !
         CALL tridiff(ix,n1-1,ix,ak,dk,ck,bk,xk,yk)

         DO k = 2, n1-1
            DO l = 1, ix
               IF (m+yoffset(wryid)+l+xoffset(wrxid) > 2) bk(l,k) = aimag(s1(l,m,k))
               IF (m+yoffset(wryid)+l+xoffset(wrxid) > 2) s1(l,m,k) = xk(l,k)
            END DO
         END DO
      
         CALL tridiff(ix,n1-1,ix,ak,dk,ck,bk,xk,yk)

         DO k = 2, n1-1
            DO l = 1, ix
               IF (m+yoffset(wryid)+l+xoffset(wrxid) > 2)                &
                  s1(l,m,k) = cmplx(REAL(s1(l,m,k)),xk(l,k))
            END DO
         END DO

      END DO

   END SUBROUTINE trdprs
   !
   !---------------------------------------------------------------------
   ! Subroutine Prs_grd: apply the pressure gradient term
   !
   SUBROUTINE prs_grd(n1,n2,n3,p,u,v,w,dz,dx,dy,dtlt)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: p(n1,n2,n3),dz(n1),dx,dy,dtlt
      REAL, INTENT (inout) :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)

      INTEGER :: i,j,k

      DO j = 1, n3-1
         DO i = 1, n2-1
            DO k = 2, n1-1
               IF(k /= n1-1)w(k,i,j) = w(k,i,j)-dz(k)*2.*dtlt*(p(k+1,i,j)-p(k,i,j))
               u(k,i,j) = u(k,i,j)-dx*2.*dtlt*(p(k,i+1,j)-p(k,i,j))
               v(k,i,j) = v(k,i,j)-dy*2.*dtlt*(p(k,i,j+1)-p(k,i,j))
            END DO
         END DO
      END DO

   END SUBROUTINE prs_grd
   !
   !---------------------------------------------------------------------
   ! Subroutine Prs_cor: correlate the pressure tendency with velocity
   ! field for TKE budget
   !
   SUBROUTINE prs_cor(n1,n2,n3,p,u,v,w,dz,dx,dy,th00)

      USE stat, ONLY : updtst
      USE util, ONLY : get_cor

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: p(n1,n2,n3),dz(n1),dx,dy,th00
      REAL, INTENT (in)    :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)

      REAL, DIMENSION (n2,n3) :: pgx, pgy, pgz, ufld, vfld, wfld
      REAL                    :: fx, fy, fz, v1da(n1), v1db(n1), v1dc(n1)
      INTEGER                 :: i,j,k,ip1,jp1

      v1da = 0.0
      v1db = 0.0
      v1dc = 0.0

      DO k = 2, n1-1
         fx = dx*th00
         fy = dy*th00
         fz = dz(k)*th00
         DO j = 1, n3
            DO i = 1, n2
               ip1 = min(n2,i+1)
               jp1 = min(n3,j+1)
               pgx(i,j) = -fx*(p(k,ip1,j)-p(k,i,j))
               pgy(i,j) = -fy*(p(k,i,jp1)-p(k,i,j))
               pgz(i,j) = -fz*(p(k+1,i,j)-p(k,i,j))
               ufld(i,j) = u(k,i,j)
               vfld(i,j) = v(k,i,j)
               wfld(i,j) = w(k,i,j)
            END DO
         END DO
         v1da(k) = get_cor(1,n2,n3,1,ufld,pgx)
         v1db(k) = get_cor(1,n2,n3,1,vfld,pgy)
         v1dc(k) = get_cor(1,n2,n3,1,wfld,pgz)
      END DO
      CALL updtst(n1,'prs',1,v1da,1)
      CALL updtst(n1,'prs',2,v1db,1)
      CALL updtst(n1,'prs',3,v1dc,1)

   END SUBROUTINE prs_cor
   !
   ! --------------------------------------------------------------------
   ! Subroutine chk_tsplt
   !
   SUBROUTINE chk_tsplt(n1,n2,n3,up,vp,wp,uc,vc,wc)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: up(n1,n2,n3),vp(n1,n2,n3),wp(n1,n2,n3)
      REAL, INTENT (in)    :: uc(n1,n2,n3),vc(n1,n2,n3),wc(n1,n2,n3)

      REAL :: wmx,umx,vmx

      wmx = maxval((wp-wc)/(wp+wc+1.e-5))
      umx = maxval((up-uc)/(up+uc+1.e-5))
      vmx = maxval((vp-vc)/(vp+vc+1.e-5))

   END SUBROUTINE chk_tsplt
   !
   !----------------------------------------------------------------------
   ! Subroutine Asselin:  Applies the asselin filter in two stages
   ! depending on the value of iac
   !
   SUBROUTINE asselin(iac)

      USE grid, ONLY : a_up,a_vp,a_wp,a_uc,a_vc,a_wc, nxyzp, runtype

      INTEGER :: iac
      INTEGER, SAVE :: ncall = 0

      IF (runtype == 'HISTORY') ncall = 1

      CALL predict(nxyzp,a_uc,a_up,iac,ncall)
      CALL predict(nxyzp,a_vc,a_vp,iac,ncall)
      CALL predict(nxyzp,a_wc,a_wp,iac,ncall)
      IF (iac == 2) ncall = ncall+1

   END SUBROUTINE asselin
   !
   !----------------------------------------------------------------------
   ! Subroutine predict:  This subroutine advances the leapfrog terms
   ! in two stages.  It applies the filter equation:
   !
   !         a(n) = a(n) + eps * (a(n-1) - 2*a(n) + a(n+1))
   !
   ! the first stage of the filter applies all but the a(n+2) term.
   ! the second stage renames the variables and applies this term, i.e.,
   ! a(n+1) -> a(n), a(n+2) -> a(n+1).  Note that for iac=2 ap=a(n+2)
   ! because the tendencies have been updated in pressure solver.  Durran,
   ! in his text cites values of eps of 0.2 for convective cloud models,
   ! we seem to get by with eps=0.1, perhaps because of the coupling
   ! provided by the staggered forward step.
   !
   SUBROUTINE predict(npts,ac,ap,iac,iflag)

      INTEGER :: m,npts,iac,iflag
      REAL    :: ac(npts),ap(npts),af(npts),epsu
      REAL, PARAMETER :: eps = 0.1

      epsu = eps
      IF (iflag == 0) epsu = 0.5
      IF (iac == 1) THEN
         DO m = 1, npts
            ac(m)=ac(m)+epsu*(ap(m)-2.*ac(m))
         END DO
      ELSE IF (iac == 2) THEN
         DO m = 1, npts
            af(m) = ap(m)
            ap(m) = ac(m)+epsu*af(m)
            ac(m) = af(m)
         END DO
      END IF

   END SUBROUTINE predict
   !
   !----------------------------------------------------------------------
   ! Subroutine Velocity_bcs: Applies boundary conditions on veolicities
   !
   SUBROUTINE velocity_bcs

      USE grid, ONLY : a_up,a_vp,a_wp,a_uc,a_vc,a_wc,a_pexnr,a_press,           &
                       nxp, nyp, nzp, dzm
      USE util, ONLY : velset, sclrset

      CALL velset(nzp,nxp,nyp,a_up,a_vp,a_wp)
      CALL velset(nzp,nxp,nyp,a_uc,a_vc,a_wc)
      CALL sclrset('grad',nzp,nxp,nyp,a_pexnr,dzm)
      CALL sclrset('grad',nzp,nxp,nyp,a_press,dzm)

   END SUBROUTINE velocity_bcs

END MODULE prss
