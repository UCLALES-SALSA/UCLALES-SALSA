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
module advf

  implicit none

  integer :: lmtr = 1

contains
  !
  !----------------------------------------------------------------------
  ! subroutine fadvect: This is the driver for scalar advection.  It
  ! advects using the average of the velocities at the current and past
  ! times.
  !
  subroutine fadvect
    use grid, only : a_up, a_vp, a_wp, a_uc, a_vc, a_wc, a_rc, a_qp, newsclr  &
         , nscl, a_sp, a_st, dn0 , nxp, nyp, nzp, dtl  &
         , dzt, dzm, zt, dxi, dyi, isgstyp
    use stat, only      : sflg, updtst
    use util, only      : get_avg3

    real    :: v1da(nzp), a_tmp1(nzp,nxp,nyp), a_tmp2(nzp,nxp,nyp)
    integer :: n
    logical :: iw
    !
    ! diagnose liquid water flux
    !
    if (sflg) then
       a_tmp1=a_rc
       call add_vel(nzp,nxp,nyp,a_tmp2,a_wp,a_wc,.false.)
       call mamaos(nzp,nxp,nyp,a_tmp2,a_rc,a_tmp1,zt,dzm,dn0,dtl,.false.)
       call get_avg3(nzp,nxp,nyp,a_tmp2,v1da)
       call updtst(nzp,'adv',0,v1da,1)
    end if
    !
    ! loop through the scalar table, setting iscp and isct to the
    ! appropriate scalar pointer and do the advection, also add large
    ! scale subsidence.  Don't advect TKE here since it resides at a
    ! w-point
    !
    do n=1,nscl
       call newsclr(n)
      IF ( ANY(a_sp /= 0.0 ) ) THEN ! TR added: no need to calculate advection for zero arrays
       a_tmp1=a_sp

       if (isgstyp > 1 .and. associated(a_qp,a_sp)) then
          iw= .true.
       else
          iw= .false.
       end if

       call add_vel(nzp,nxp,nyp,a_tmp2,a_vp,a_vc,iw)
       call mamaos_y(nzp,nxp,nyp,a_tmp2,a_sp,a_tmp1,dyi,dtl)

       call add_vel(nzp,nxp,nyp,a_tmp2,a_up,a_uc,iw)
       call mamaos_x(nzp,nxp,nyp,a_tmp2,a_sp,a_tmp1,dxi,dtl)

       call add_vel(nzp,nxp,nyp,a_tmp2,a_wp,a_wc,iw)
       call mamaos(nzp,nxp,nyp,a_tmp2,a_sp,a_tmp1,dzt,dzm,dn0,dtl,iw)
       if (sflg) then
          call get_avg3(nzp,nxp,nyp,a_tmp2,v1da)
          call updtst(nzp,'adv',n,v1da,1)
       end if

       call advtnd(nzp,nxp,nyp,a_sp,a_tmp1,a_st,dtl)
      ELSEIF (sflg) THEN
       ! Averages & statistics even for zeros (might be non-zero elsewhere)
       a_tmp2(:,:,:)=0.
       call get_avg3(nzp,nxp,nyp,a_tmp2,v1da)
       call updtst(nzp,'adv',n,v1da,1)
      ENDIF
    end do

  end subroutine fadvect
  !
  !----------------------------------------------------------------------
  ! subroutine add_vel: Adds current and past timelevels of velocities
  ! into scratch arrays for advection
  !
  subroutine add_vel(n1,n2,n3,su,up,uc,lwpt)

    integer, intent(in) ::  n1,n2,n3
    real, intent(in)    ::  up(n1,n2,n3),uc(n1,n2,n3)
    logical, intent (in)::  lwpt
    real, intent(out)   ::  su(n1,n2,n3)

    integer i,j,k

    if (lwpt) then
       do j=1,n3
          do i=1,n2
             do k=1,n1-1
                su(k,i,j)=0.25*(up(k,i,j)+uc(k,i,j)+up(k+1,i,j)+uc(k+1,i,j))
             end do
             su(n1,i,j)=su(n1-1,i,j)
          enddo
       enddo
    else
       do j=1,n3
          do i=1,n2
             do k=1,n1
                su(k,i,j)=0.5*(up(k,i,j)+uc(k,i,j))
             end do
          enddo
       enddo
    end if

  end subroutine add_vel
  !
  ! ----------------------------------------------------------------------
  ! Subroutine advtnd: Backs out the advective tendencies
  !
  subroutine advtnd(n1,n2,n3,varo,varn,tnd,dt)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)   :: varo(n1,n2,n3),varn(n1,n2,n3),dt
    real, intent(inout)  :: tnd(n1,n2,n3)

    real :: dti
    integer :: i,j,k

    dti=1./dt

    do j=3,n3-2
       do i=3,n2-2
          tnd(1,i,j)  = 0.
          do k=2,n1-1
             tnd(k,i,j)=tnd(k,i,j)+(varn(k,i,j)-varo(k,i,j))*dti
          end do
          tnd(n1,i,j) = 0.
       enddo
    enddo

  end subroutine advtnd
  !
  !----------------------------------------------------------------------
  ! Subroutine mamaos: An alternative second order flux limited scheme
  ! written by Verica and Christiana as part of the MAMAOS program.
  !
  ! July 21, 2003
  !
  subroutine mamaos(n1,n2,n3,wp,scp0,scp,dzt,dzm,dn0,dt,lwpt)

    use mpi_interface, only : myid, appl_abort
    integer, intent (in)    :: n1,n2,n3
    real, intent (in)       :: scp0(n1,n2,n3)
    real, intent (in)       :: dn0(n1),dzt(n1),dzm(n1)
    real, intent (in)       :: dt
    logical, intent (in)    :: lwpt
    real, intent (inout)    :: wp(n1,n2,n3),scp(n1,n2,n3)

    real    :: density(n1)   ! averaged density
    real    :: dzt_local(n1) ! grid spacing for scalars
    real    :: dzm_local(n1) ! grid spacing for velocity
    real    :: cfl(n1)       ! cfl numbers at the interface (staggered)
    real    :: C(n1)         ! limiter
    real    :: r(n1)         ! slope ratio
    real    :: wpdn(n1)      ! momentum: wp*density
    integer :: i, j, k, kp1, k1, k2
    integer :: gamma
    !
    ! initialize fields for use later
    !
    do k = 1, n1
       kp1 = min(k+1,n1)
       density(k) = 0.5 * (dn0(k) + dn0(kp1))
       if (lwpt) then
          dzt_local(k) = dzm(k)
          dzm_local(k) = dzt(kp1)
       else
          dzt_local(k) = dzt(k)
          dzm_local(k) = dzm(k)
       endif
    enddo

    do j = 3, n3-2
       do i = 3, n2-2
          !
          ! compute CFL and momentum
          !
          do k = 1, n1-1
             cfl(k)  = wp(k,i,j) * dt * dzm_local(k)
             wpdn(k) = wp(k,i,j) * density(k)
             if (abs(cfl(k)) > 1.0) then
                if (myid == 0) WRITE(*,*) '  ABORTING: mamaos_z', cfl(k),wp(k,i,j),k,i,j
                call appl_abort (0)
             end if
          enddo
          !
          ! calculate the ratio of slopes
          !
          do k = 1, n1-1
             gamma = -sign(1.,cfl(k))
             if (abs(scp0(k+1,i,j)-scp0(k,i,j)) > spacing(scp0(k,i,j))) then
                k2 = max(1,k+gamma)
                k1 = min(n1,k+gamma+1)
                r(k) = (scp0(k1,i,j) - scp0(k2,i,j)) / (scp0(k+1,i,j) - scp0(k,i,j))
             else
                r(k) = 0.
             endif
          enddo
          !
          ! calculate the flux limiters
          !
          select case (lmtr)
          case (1) ! minmod
             do k = 1, n1-2
                C(k) = max(0., min(1., r(k)))
             enddo
          case(2)  ! superbee
             do k = 1, n1-2
                C(k) = max(0., min(1., 2.*r(k)), min(2., r(k)))
             enddo
          case(3)  ! mc
             do k = 1, n1-2
                C(k) = max(0., min(2.*r(k),(1.+r(k))/2., 2.))
             enddo
          case(4)  ! van Leer
             do k = 1, n1-2
                C(k) = (r(k) + abs(r(k)))/(1. + abs(r(k)))
             enddo
          case default ! no limiter
             do k = 1, n1-2
                C(k) = 1.0
             enddo
          end select

          wp(1,i,j) = 0.
          wp(n1-1,i,j) = 0.
          do k = 2, n1-2
             wp(k,i,j) = 0.5 * wpdn(k) * (scp0(k+1,i,j)+scp0(k,i,j)) - &
                  0.5 * (scp0(k+1,i,j)-scp0(k,i,j)) *                  &
                  ((1.-C(k))*abs(wpdn(k)) + wpdn(k)*cfl(k)*C(k))
          end do
          do k = 2,n1-1
             scp(k,i,j) = scp(k,i,j) - ((wp(k,i,j)-wp(k-1,i,j)) -      &
                  scp0(k,i,j)*(wpdn(k)-wpdn(k-1))) *                   &
                  dt*dzt_local(k)/dn0(k)
          enddo

       enddo
    enddo

  end subroutine mamaos
 !
  !----------------------------------------------------------------------
  ! Subroutine mamaos_x: An alternative second order flux limited scheme
  ! for advection in the x direction.  (adapted from mamaos)
  !
  ! September 3, 2003
  !
  subroutine mamaos_x(n1,n2,n3,up,scp0,scp,dxi,dt)

    use mpi_interface, only : myid, appl_abort

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dxi,dt
    real, intent (in)    :: scp0(n1,n2,n3),up(n1,n2,n3)
    real, intent (inout) :: scp(n1,n2,n3)

    real    :: cfl(n2,n1)       ! cfl numbers at the interface (staggered)
    real    :: C(n2,n1)         ! limiter
    real    :: r(n2,n1)         ! slope ratio
    real    :: scr(n2,n1)       ! flux scratch array
    integer :: i, j, k, i1, i2
    integer :: gamma
    !

    do j = 3,n3-2
       !
       ! compute CFL and scr array for down-grid value of scalar
       !
       do k = 2, n1-1
          do i = 1,n2-1
             cfl(i,k)  = up(k,i,j) * dt * dxi
             scr(i,k)  = scp0(k,i+1,j)
             if (abs(cfl(i,k)) > 1.0) then
                if (myid == 0) print *, '  ABORTING: mamaos_x'
                call appl_abort(0)
             end if
          end do
       end do
          !
          ! calculate the ratio of slopes
          !
       do k = 2, n1-1
          do i = 2,n2-2
             gamma = int(-sign(1.,cfl(i,k)))
             if (abs(scr(i,k) - scp0(k,i,j)) > spacing(scr(i,k))) then
                i2 = i+gamma
                i1 = i+gamma+1
                r(i,k) = (scp0(k,i1,j)-scp0(k,i2,j))/(scr(i,k)-scp0(k,i,j))
             else
                r(i,k) = 0.
             endif

             select case (lmtr)
             case (1) ! minmod
                C(i,k) = max(0., min(1., r(i,k)))
             case(2)  ! superbee
                C(i,k) = max(0., min(1., 2.*r(i,k)), min(2., r(i,k)))
             case(3)  ! mc
                C(i,k) = max(0., min(2.*r(i,k),(1.+r(i,k))/2., 2.))
             case(4)  ! van Leer
                C(i,k) = (r(i,k) + abs(r(i,k)))/(1. + abs(r(i,k)))
             case default ! no limiter
                C(i,k) = 1.0
             end select

             scr(i,k) = 0.5 * up(k,i,j) * (scr(i,k)+scp0(k,i,j)) -      &
                  0.5 * (scr(i,k)-scp0(k,i,j)) *                        &
                  ((1.-C(i,k))*abs(up(k,i,j)) + up(k,i,j)*cfl(i,k)*C(i,k))
          end do

          do i = 3,n2-2
             scp(k,i,j) = scp(k,i,j) - ((scr(i,k)-scr(i-1,k)) -         &
                  scp0(k,i,j)*(up(k,i,j)-up(k,i-1,j)))*dt*dxi
          enddo
       enddo

    enddo

  end subroutine mamaos_x
  !
  !----------------------------------------------------------------------
  ! Subroutine mamaos_y: An alternative second order flux limited scheme
  ! for advection in the y direction.  (adapted from mamaos)
  !
  ! September 3, 2003
  !
  subroutine mamaos_y(n1,n2,n3,vp,scp0,scp,dyi,dt)

    use mpi_interface, only : myid, appl_abort

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dyi,dt
    real, intent (in)    :: scp0(n1,n2,n3),vp(n1,n2,n3)
    real, intent (inout) :: scp(n1,n2,n3)

    real    :: cfl(n3,n1)       ! cfl numbers at the interface (staggered)
    real    :: C(n3,n1)         ! limiter
    real    :: r(n3,n1)         ! slope ratio
    real    :: scr(n3,n1)       ! flux scratch array
    integer :: i, j, k, j1, j2
    integer :: gamma
    !

    do i = 1, n2
       !
       ! compute CFL and scr array for down-grid value of scalar
       !
       do k = 2, n1-1
          do j = 1,n3-1
             cfl(j,k)  = vp(k,i,j) * dt * dyi
             scr(j,k)  = scp0(k,i,j+1)
             if (abs(cfl(j,k)) > 1.0) then
                if (myid == 0) print *, '  ABORTING: mamaos_y'
                call appl_abort(0)
             end if
          end do
       end do
          !
          ! calculate the ratio of slopes
          !
       do k = 2, n1-1
          do j = 2,n3-2
             gamma = int(-sign(1.,cfl(j,k)))
             if (abs(scr(j,k) - scp0(k,i,j)) > spacing(scr(j,k))) then
                j2 = j+gamma
                j1 = j+gamma+1
                r(j,k) = (scp0(k,i,j1)-scp0(k,i,j2))/(scr(j,k)-scp0(k,i,j))
             else
                r(j,k) = 0.
             endif

             select case (lmtr)
             case (1) ! minmod
                C(j,k) = max(0., min(1., r(j,k)))
             case(2)  ! superbee
                C(j,k) = max(0., min(1., 2.*r(j,k)), min(2., r(j,k)))
             case(3)  ! mc
                C(j,k) = max(0., min(2.*r(j,k),(1.+r(j,k))/2., 2.))
             case(4)  ! van Leer
                C(j,k) = (r(j,k) + abs(r(j,k)))/(1. + abs(r(j,k)))
             case default ! no limiter
                C(j,k) = 1.0
             end select

             scr(j,k) = 0.5 * vp(k,i,j) * (scr(j,k)+scp0(k,i,j)) -      &
                  0.5 * (scr(j,k)-scp0(k,i,j)) *                        &
                  ((1.-C(j,k))*abs(vp(k,i,j)) + vp(k,i,j)*cfl(j,k)*C(j,k))
          end do

          do j = 3,n3-2
             scp(k,i,j) = scp(k,i,j) - ((scr(j,k)-scr(j-1,k)) -         &
                  scp0(k,i,j)*(vp(k,i,j)-vp(k,i,j-1)))*dt*dyi
          enddo
       enddo

    enddo

  end subroutine mamaos_y

end module advf
