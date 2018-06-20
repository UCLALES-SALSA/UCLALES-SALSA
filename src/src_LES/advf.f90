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
MODULE advf

  IMPLICIT NONE

  INTEGER :: lmtr = 1

CONTAINS
  !
  !----------------------------------------------------------------------
  ! Subroutine fadvect: This is the driver for scalar advection.  It
  ! advects using the average of the velocities at the current and past
  ! times.
  !
  SUBROUTINE fadvect
    USE grid, ONLY : a_up, a_vp, a_wp, a_uc, a_vc, a_wc, a_rc, a_qp, newsclr,  &
                     nscl, a_sp, a_st, dn0 , nxp, nyp, nzp, dtlt,  &
                     dzt, dzm, zt, dxi, dyi, level, isgstyp
    USE stat, ONLY : sflg, updtst
    USE util, ONLY : get_avg3

    REAL    :: v1da(nzp), a_tmp1(nzp,nxp,nyp), a_tmp2(nzp,nxp,nyp)
    INTEGER :: n
    LOGICAL :: iw
    !
    ! diagnose liquid water flux
    !
    IF (sflg .AND. level > 1) THEN
       a_tmp1 = a_rc
       CALL add_vel(nzp,nxp,nyp,a_tmp2,a_wp,a_wc,.FALSE.)
       CALL mamaos(nzp,nxp,nyp,a_tmp2,a_rc,a_tmp1,zt,dzm,dn0,dtlt,.FALSE.)
       CALL get_avg3(nzp,nxp,nyp,a_tmp2,v1da)
       CALL updtst(nzp,'adv',0,v1da,1)
    END IF
    !
    ! loop through the scalar table, setting iscp and isct to the
    ! appropriate scalar pointer and do the advection, also add large
    ! scale subsidence.  Don't advect TKE here since it resides at a
    ! w-point
    !
    DO n = 1, nscl
      CALL newsclr(n)

      IF ( ANY(a_sp /= 0.0 ) ) THEN ! TR added: no need to calculate advection for zero arrays
         a_tmp1 = a_sp

         IF (isgstyp > 1 .AND. associated(a_qp,a_sp)) THEN
            iw = .TRUE.
         ELSE
            iw = .FALSE.
         END IF

         CALL add_vel(nzp,nxp,nyp,a_tmp2,a_vp,a_vc,iw)
         CALL mamaos_y(nzp,nxp,nyp,a_tmp2,a_sp,a_tmp1,dyi,dtlt)

         CALL add_vel(nzp,nxp,nyp,a_tmp2,a_up,a_uc,iw)
         CALL mamaos_x(nzp,nxp,nyp,a_tmp2,a_sp,a_tmp1,dxi,dtlt)

         CALL add_vel(nzp,nxp,nyp,a_tmp2,a_wp,a_wc,iw)
         CALL mamaos(nzp,nxp,nyp,a_tmp2,a_sp,a_tmp1,dzt,dzm,dn0,dtlt,iw)

         IF (sflg) THEN
            CALL get_avg3(nzp,nxp,nyp,a_tmp2,v1da)
            CALL updtst(nzp,'adv',n,v1da,1)
         END IF

         CALL advtnd(nzp,nxp,nyp,a_sp,a_tmp1,a_st,dtlt)

      ELSE IF (sflg) THEN
         ! Averages & statistics even for zeros (might be non-zero elsewhere)
         a_tmp2(:,:,:) = 0.
         CALL get_avg3(nzp,nxp,nyp,a_tmp2,v1da)
         CALL updtst(nzp,'adv',n,v1da,1)
      END IF

    END DO

  END SUBROUTINE fadvect
  !
  !----------------------------------------------------------------------
  ! Subroutine newdroplet: Calculates the actual tendency in cloud droplets
  ! and aerosols due to cloud activation.
  !
  SUBROUTINE newdroplet(pactmask)
    USE mo_submctl, ONLY : ncld,nbins,ica,fca,eps,spec
    USE grid, ONLY : nxp,nyp,nzp,dzt,            &
                     a_wp,a_wc,  &
                     a_naerop, a_naerot, a_maerop, a_maerot,  &
                     a_ncloudt, a_mcloudt,  &
                     a_nactd,  a_vactd,  a_rt
    IMPLICIT NONE

    LOGICAL, INTENT(in) :: pactmask(nzp,nxp,nyp)
    REAL    :: zw(nxp,nyp)
    INTEGER :: ii,jj,kk,bb,bbpar,ss,mm,mmpar,kp1,nc
    REAL    :: frac(nxp,nyp)
    REAL    :: fix_flux(nxp,nyp)
    REAL    :: dn(nxp,nyp),dv(nxp,nyp)

    ! Calculate the rate of change in cloud droplets due to cloud activation

    DO kk = 2, nzp-1

       kp1 = kk + 1 ! Tendency is applied in the level above where the potential number of activated is computed
       zw = 0.25*( a_wp(kk,:,:) + a_wc(kk,:,:) + a_wp(kp1,:,:) + a_wc(kp1,:,:) )
       frac = 0.

       DO bb = ica%cur, fca%cur

          bbpar = ica%par + (bb-ica%cur)

          dn(:,:) = zw(:,:)*a_nactd(kk,:,:,bb)*dzt(kk)

          DO jj = 3, nyp-2
             DO ii = 3, nxp-2
                fix_flux(ii,jj) = max(min(1.0,a_naerot(kp1,ii,jj,bbpar)/max(eps,dn(ii,jj))),0.8)
             END DO
          END DO

          ! Add to cloud droplets
          a_ncloudt(kp1,:,:,bb) = a_ncloudt(kp1,:,:,bb) +   &
                                  MERGE( dn(:,:)*fix_flux(:,:), 0., pactmask(kk,:,:) )

          ! Remove from aerosols
          a_naerot(kp1,:,:,bbpar) = a_naerot(kp1,:,:,bbpar) -   &
                                    MERGE( dn(:,:)*fix_flux(:,:), 0., pactmask(kk,:,:) )

          ! Change in dry ccn/aerosol mass
          DO ss = 1, spec%getNSpec()-1

             mm = (ss-1)*ncld + bb
             mmpar = (ss-1)*nbins + bbpar

             dv(:,:) = zw(:,:)*a_vactd(kk,:,:,mm)*dzt(kk)

             DO jj = 3, nyp-2
                DO ii = 3, nxp-2
                   fix_flux(ii,jj) = max(min(1.0,a_maerot(kp1,ii,jj,mmpar)/max(eps,dv(ii,jj))),0.8)
                END DO
             END DO

             ! Add to cloud droplets
             a_mcloudt(kp1,:,:,mm) = a_mcloudt(kp1,:,:,mm) +   &
                                     MERGE( dv(:,:)*fix_flux(:,:), 0., pactmask(kk,:,:) )

             ! Remove from aerosols
             a_maerot(kp1,:,:,mmpar) = a_maerot(kp1,:,:,mmpar) -   &
                                       MERGE( dv(:,:)*fix_flux(:,:), 0., pactmask(kk,:,:) )

          END DO ! ss

          ! Change in water content
          ! Assume that existing condensate from aerosols bins is taken in the same relation
          ! as the number concentration....
          frac(:,:) = a_nactd(kk,:,:,bb)/MAX(a_naerop(kk,:,:,bbpar),1.)

          nc = spec%getIndex('H2O')
          mm = (nc-1)*ncld + bb
          mmpar = (nc-1)*nbins + bbpar

          ! Amount of water condensed from vapor
          a_rt(kp1,:,:) = a_rt(kp1,:,:) -   &
                          MERGE( fix_flux(:,:)*zw(:,:)*(a_vactd(kk,:,:,mm)-frac(:,:)*a_maerop(kk,:,:,mmpar))*dzt(kk), 0., &
                                 pactmask(kk,:,:) )
          ! Add to cloud droplets
          a_mcloudt(kp1,:,:,mm) = a_mcloudt(kp1,:,:,mm) +   &
                                  MERGE( fix_flux(:,:)*zw(:,:)*a_vactd(kk,:,:,mm)*dzt(kk), 0., pactmask(kk,:,:) )

          ! Remove from aerosols
          a_maerot(kp1,:,:,mmpar) = a_maerot(kp1,:,:,mmpar) -   &
                                    MERGE( fix_flux(:,:)*zw(:,:)*frac(:,:)*a_maerop(kk,:,:,mmpar)*dzt(kk), 0., pactmask(kk,:,:) )

       END DO ! bb

    END DO ! kk

    ! Mask the number of activated for diagnostic output to be shown only in
    ! the points defined by the activation mask (0 elsewhere)
    DO bb = ica%cur, fca%cur
       a_nactd(:,:,:,bb) = MERGE(a_nactd(:,:,:,bb),0.,pactmask(:,:,:))
       DO ss = 1, spec%getNSpec()
          mm = (ss-1)*ncld + bb
          a_vactd(:,:,:,mm) = MERGE(a_vactd(:,:,:,mm),0.,pactmask(:,:,:))
       END DO
    END DO


  END SUBROUTINE newdroplet


  !
  !----------------------------------------------------------------------
  ! Subroutine add_vel: Adds current and past timelevels of velocities
  ! into scratch arrays for advection
  !
  SUBROUTINE add_vel(n1,n2,n3,su,up,uc,lwpt)

    INTEGER, INTENT(in)  ::  n1,n2,n3
    REAL, INTENT(in)     ::  up(n1,n2,n3),uc(n1,n2,n3)
    LOGICAL, INTENT (in) ::  lwpt
    REAL, INTENT(out)    ::  su(n1,n2,n3)

    INTEGER :: i,j,k

    IF (lwpt) THEN
       DO j= 1, n3
          DO i= 1, n2
             DO k= 1, n1-1
                su(k,i,j) = 0.25*(up(k,i,j)+uc(k,i,j)+up(k+1,i,j)+uc(k+1,i,j))
             END DO
             su(n1,i,j) = su(n1-1,i,j)
          END DO
       END DO
    ELSE
       DO j = 1, n3
          DO i = 1, n2
             DO k = 1, n1
                su(k,i,j) = 0.5*(up(k,i,j)+uc(k,i,j))
             END DO
          END DO
       END DO
    END IF

  END SUBROUTINE add_vel
  !
  ! ----------------------------------------------------------------------
  ! Subroutine advtnd: Backs out the advective tendencies
  !
  SUBROUTINE advtnd(n1,n2,n3,varo,varn,tnd,dt)

    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(in)    :: varo(n1,n2,n3),varn(n1,n2,n3),dt
    REAL, INTENT(inout) :: tnd(n1,n2,n3)

    REAL    :: dti
    INTEGER :: i,j,k

    dti = 1./dt

    DO j = 3, n3-2
       DO i = 3, n2-2
          tnd(1,i,j) = 0.
          DO k = 2, n1-1
             tnd(k,i,j) = tnd(k,i,j)+(varn(k,i,j)-varo(k,i,j))*dti
          END DO
          tnd(n1,i,j) = 0.
       END DO
    END DO

  END SUBROUTINE advtnd
  !
  !----------------------------------------------------------------------
  ! Subroutine mamaos: An alternative second order flux limited scheme
  ! written by Verica and Christiana as part of the MAMAOS program.
  !
  ! July 21, 2003
  !
  SUBROUTINE mamaos(n1,n2,n3,wp,scp0,scp,dzt,dzm,dn0,dt,lwpt)

    USE mpi_interface, ONLY : myid, appl_abort
    INTEGER, INTENT (in)    :: n1,n2,n3
    REAL, INTENT (in)       :: scp0(n1,n2,n3)
    REAL, INTENT (in)       :: dn0(n1),dzt(n1),dzm(n1)
    REAL, INTENT (in)       :: dt
    LOGICAL, INTENT (in)    :: lwpt
    REAL, INTENT (inout)    :: wp(n1,n2,n3),scp(n1,n2,n3)

    REAL    :: density(n1)   ! averaged density
    REAL    :: dzt_local(n1) ! grid spacing for scalars
    REAL    :: dzm_local(n1) ! grid spacing for velocity
    REAL    :: cfl(n1)       ! cfl numbers at the interface (staggered)
    REAL    :: C(n1)         ! limiter
    REAL    :: r(n1)         ! slope ratio
    REAL    :: wpdn(n1)      ! momentum: wp*density
    INTEGER :: i, j, k, kp1, k1, k2
    INTEGER :: gamma
    !
    ! initialize fields for use later
    !
    DO k = 1, n1
       kp1 = min(k+1,n1)
       density(k) = 0.5 * (dn0(k) + dn0(kp1))
       IF (lwpt) THEN
          dzt_local(k) = dzm(k)
          dzm_local(k) = dzt(kp1)
       ELSE
          dzt_local(k) = dzt(k)
          dzm_local(k) = dzm(k)
       END IF
    END DO

    DO j = 3, n3-2
       DO i = 3, n2-2
          !
          ! compute CFL and momentum
          !
          DO k = 1, n1-1
             cfl(k)  = wp(k,i,j) * dt * dzm_local(k)
             wpdn(k) = wp(k,i,j) * density(k)
             IF (abs(cfl(k)) > 1.0) THEN
                IF (myid == 0) WRITE(*,*) '  ABORTING: mamaos_z', cfl(k),wp(k,i,j),k,i,j
                CALL appl_abort (0)
             END IF
          END DO
          !
          ! calculate the ratio of slopes
          !
          DO k = 1, n1-1
             gamma = -INT(sign(1.,cfl(k)))
             IF (ABS(scp0(k+1,i,j)-scp0(k,i,j)) > SPACING(scp0(k,i,j)) ) THEN !.AND. scp0(k+1,i,j)-scp0(k,i,j) > 1.e-40) THEN
                k2 = max(1,k+gamma)
                k1 = min(n1,k+gamma+1)
                r(k) = (scp0(k1,i,j) - scp0(k2,i,j)) / (scp0(k+1,i,j) - scp0(k,i,j))
             ELSE
                r(k) = 0.
             END IF
          END DO
          !
          ! calculate the flux limiters
          !
          SELECT CASE (lmtr)
          CASE (1) ! minmod
             DO k = 1, n1-2
                C(k) = max(0., min(1., r(k)))
             END DO
          CASE(2)  ! superbee
             DO k = 1, n1-2
                C(k) = max(0., min(1., 2.*r(k)), min(2., r(k)))
             END DO
          CASE(3)  ! mc
             DO k = 1, n1-2
                C(k) = max(0., min(2.*r(k),(1.+r(k))/2., 2.))
             END DO
          CASE(4)  ! van Leer
             DO k = 1, n1-2
                C(k) = (r(k) + abs(r(k)))/(1. + abs(r(k)))
             END DO
          CASE DEFAULT ! no limiter
             DO k = 1, n1-2
                C(k) = 1.0
             END DO
          END SELECT

          wp(1,i,j) = 0.
          wp(n1-1,i,j) = 0.
          DO k = 2, n1-2
             wp(k,i,j) = 0.5 * wpdn(k) * (scp0(k+1,i,j)+scp0(k,i,j)) - &
                         0.5 * (scp0(k+1,i,j)-scp0(k,i,j)) *           &
                         ((1.-C(k))*abs(wpdn(k)) + wpdn(k)*cfl(k)*C(k))
          END DO
          DO k = 2, n1-1
             scp(k,i,j) = scp(k,i,j) - ((wp(k,i,j)-wp(k-1,i,j)) -      &
                          scp0(k,i,j)*(wpdn(k)-wpdn(k-1))) *           &
                          dt*dzt_local(k)/dn0(k)
          END DO

       END DO
    END DO

  END SUBROUTINE mamaos
 !
  !----------------------------------------------------------------------
  ! Subroutine mamaos_x: An alternative second order flux limited scheme
  ! for advection in the x direction.  (adapted from mamaos)
  !
  ! September 3, 2003
  !
  SUBROUTINE mamaos_x(n1,n2,n3,up,scp0,scp,dxi,dt)

    USE mpi_interface, ONLY : myid, appl_abort

    INTEGER, INTENT (in) :: n1,n2,n3
    REAL, INTENT (in)    :: dxi,dt
    REAL, INTENT (in)    :: scp0(n1,n2,n3),up(n1,n2,n3)
    REAL, INTENT (inout) :: scp(n1,n2,n3)

    REAL    :: cfl(n2,n1)       ! cfl numbers at the interface (staggered)
    REAL    :: C(n2,n1)         ! limiter
    REAL    :: r(n2,n1)         ! slope ratio
    REAL    :: scr(n2,n1)       ! flux scratch array
    INTEGER :: i, j, k, i1, i2
    INTEGER :: gamma
    !

    DO j = 3, n3-2
       !
       ! compute CFL and scr array for down-grid value of scalar
       !
       DO k = 2, n1-1
          DO i = 1, n2-1
             cfl(i,k) = up(k,i,j) * dt * dxi
             scr(i,k) = scp0(k,i+1,j)
             IF (abs(cfl(i,k)) > 1.0) THEN
                IF (myid == 0) PRINT *, '  ABORTING: mamaos_x'
                CALL appl_abort(0)
             END IF
          END DO
       END DO
          !
          ! calculate the ratio of slopes
          !
       DO k = 2, n1-1
          DO i = 2, n2-2
             gamma = INT(-sign(1.,cfl(i,k)))
             IF (abs(scr(i,k) - scp0(k,i,j)) > SPACING(scr(i,k))) THEN
                i2 = i+gamma
                i1 = i+gamma+1
                r(i,k) = (scp0(k,i1,j)-scp0(k,i2,j))/(scr(i,k)-scp0(k,i,j))
             ELSE
                r(i,k) = 0.
             END IF

             SELECT CASE (lmtr)
             CASE (1) ! minmod
                C(i,k) = max(0., min(1., r(i,k)))
             CASE(2)  ! superbee
                C(i,k) = max(0., min(1., 2.*r(i,k)), min(2., r(i,k)))
             CASE(3)  ! mc
                C(i,k) = max(0., min(2.*r(i,k),(1.+r(i,k))/2., 2.))
             CASE(4)  ! van Leer
                C(i,k) = (r(i,k) + abs(r(i,k)))/(1. + abs(r(i,k)))
             CASE DEFAULT ! no limiter
                C(i,k) = 1.0
             END SELECT

             scr(i,k) = 0.5 * up(k,i,j) * (scr(i,k)+scp0(k,i,j)) -      &
                        0.5 * (scr(i,k)-scp0(k,i,j)) *                  &
                        ((1.-C(i,k))*abs(up(k,i,j)) + up(k,i,j)*cfl(i,k)*C(i,k))
          END DO

          DO i = 3, n2-2
             scp(k,i,j) = scp(k,i,j) - ((scr(i,k)-scr(i-1,k)) -         &
                          scp0(k,i,j)*(up(k,i,j)-up(k,i-1,j)))*dt*dxi
          END DO
       END DO

    END DO

  END SUBROUTINE mamaos_x
  !
  !----------------------------------------------------------------------
  ! Subroutine mamaos_y: An alternative second order flux limited scheme
  ! for advection in the y direction.  (adapted from mamaos)
  !
  ! September 3, 2003
  !
  SUBROUTINE mamaos_y(n1,n2,n3,vp,scp0,scp,dyi,dt)

    USE mpi_interface, ONLY : myid, appl_abort

    INTEGER, INTENT (in) :: n1,n2,n3
    REAL, INTENT (in)    :: dyi,dt
    REAL, INTENT (in)    :: scp0(n1,n2,n3),vp(n1,n2,n3)
    REAL, INTENT (inout) :: scp(n1,n2,n3)

    REAL    :: cfl(n3,n1)       ! cfl numbers at the interface (staggered)
    REAL    :: C(n3,n1)         ! limiter
    REAL    :: r(n3,n1)         ! slope ratio
    REAL    :: scr(n3,n1)       ! flux scratch array
    INTEGER :: i, j, k, j1, j2
    INTEGER :: gamma
    !

    DO i = 1, n2
       !
       ! compute CFL and scr array for down-grid value of scalar
       !
       DO k = 2, n1-1
          DO j = 1, n3-1
             cfl(j,k) = vp(k,i,j) * dt * dyi
             scr(j,k) = scp0(k,i,j+1)
             IF (abs(cfl(j,k)) > 1.0) THEN
                IF (myid == 0) PRINT *, '  ABORTING: mamaos_y'
                CALL appl_abort(0)
             END IF
          END DO
       END DO
          !
          ! calculate the ratio of slopes
          !
       DO k = 2, n1-1
          DO j = 2, n3-2
             gamma = INT(-sign(1.,cfl(j,k)))
             IF (ABS(scr(j,k) - scp0(k,i,j)) > SPACING(scr(j,k)) ) THEN
                j2 = j+gamma
                j1 = j+gamma+1
                r(j,k) = (scp0(k,i,j1)-scp0(k,i,j2))/(scr(j,k)-scp0(k,i,j))
             ELSE
                r(j,k) = 0.
             END IF

             SELECT CASE (lmtr)
             CASE (1) ! minmod
                C(j,k) = max(0., min(1., r(j,k)))
             CASE(2)  ! superbee
                C(j,k) = max(0., min(1., 2.*r(j,k)), min(2., r(j,k)))
             CASE(3)  ! mc
                C(j,k) = max(0., min(2.*r(j,k),(1.+r(j,k))/2., 2.))
             CASE(4)  ! van Leer
                C(j,k) = (r(j,k) + abs(r(j,k)))/(1. + abs(r(j,k)))
             CASE DEFAULT ! no limiter
                C(j,k) = 1.0
             END SELECT

             scr(j,k) = 0.5 * vp(k,i,j) * (scr(j,k)+scp0(k,i,j)) -      &
                        0.5 * (scr(j,k)-scp0(k,i,j)) *                  &
                        ((1.-C(j,k))*abs(vp(k,i,j)) + vp(k,i,j)*cfl(j,k)*C(j,k))
          END DO

          DO j = 3, n3-2
             scp(k,i,j) = scp(k,i,j) - ((scr(j,k)-scr(j-1,k)) -         &
                          scp0(k,i,j)*(vp(k,i,j)-vp(k,i,j-1)))*dt*dyi
          END DO
       END DO

    END DO

  END SUBROUTINE mamaos_y

END MODULE advf
