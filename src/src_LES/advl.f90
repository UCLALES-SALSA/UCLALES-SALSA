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
module advl

  implicit none

contains
  !
  ! ----------------------------------------------------------------------
  ! LADVECT:  This does the leap-frog advection using a fourth
  ! order centered in space algorithm.  Boundary point fluxes are not 
  ! computed... hence the routine is hardwired to a cyclic domain.  
  ! Vertical advection is density weighted consistent with the anelastic
  ! approximation
!
  subroutine ladvect

    use grid, only : a_uc, a_vc, a_wc, a_ut, a_vt, a_wt,      &
         nxp, nyp, nzp, dzt, dzm, dxi, dyi, dn0
    use stat, only : sflg, updtst, acc_tend
    use util, only : get_avg3

    real ::  v1da(nzp), v1db(nzp), v1dc(nzp), v1dd(nzp), v1de(nzp), &
         a_tmp1(nzp,nxp,nyp), a_tmp2(nzp,nxp,nyp)
    !
    ! prepare density weights for use vertical advection
    !
    if (sflg) call acc_tend(nzp,nxp,nyp,a_uc,a_vc,a_wc,a_ut,a_vt,a_wt,        &
         v1dc,v1dd,v1de,1,'adv')

    call advl_prep(nzp,nxp,nyp,a_wc,a_tmp1,dn0,dzt,dzm,v1da,v1db)
    !
    ! advection of u by (u,v,w) all at current timelevel.  also when flag
    ! is set updated statistical array with uw flux derived from ladvzu
    !
    call ladvxu(nzp,nxp,nyp,a_uc,a_ut,a_tmp2,dxi)
    call ladvyu(nzp,nxp,nyp,a_uc,a_ut,a_vc,a_tmp2,dyi)
    call ladvzu(nzp,nxp,nyp,a_uc,a_ut,a_tmp1,a_tmp2,v1da)
    if (sflg) then
       call get_avg3(nzp,nxp,nyp,a_tmp2,v1da)
       call updtst(nzp,'adv',-1,v1da,1)
    end if
    !
    ! advection of v by (u,v,w) all at current timelevel.  also when flag
    ! is set updated statistical array with uw flux derived from ladvzu
    !    
    call ladvxv(nzp,nxp,nyp,a_uc,a_vc,a_vt,a_tmp2,dxi)
    call ladvyv(nzp,nxp,nyp,a_vc,a_vt,a_tmp2,dyi)
    call ladvzv(nzp,nxp,nyp,a_vc,a_vt,a_tmp1,a_tmp2,v1da)
    if (sflg) then
       call get_avg3(nzp,nxp,nyp,a_tmp2,v1da)
       call updtst(nzp,'adv',-2,v1da,1)
    end if
    !
    ! advection of w by (u,v,w) all at current timelevel.  also when flag
    ! is set updated statistical array with uw flux derived from ladvzu
    !
    call ladvxw(nzp,nxp,nyp,a_uc,a_wc,a_wt,a_tmp2,dxi)
    call ladvyw(nzp,nxp,nyp,a_vc,a_wc,a_wt,a_tmp2,dyi)
    call ladvzw(nzp,nxp,nyp,a_wc,a_wt,a_tmp1,a_tmp2,v1db)

    if (sflg) then
       call get_avg3(nzp,nxp,nyp,a_tmp2,v1da)
       call updtst(nzp,'adv',-3,v1da,1)
    end if

    if (sflg) call acc_tend(nzp,nxp,nyp,a_uc,a_vc,a_wc,a_ut,a_vt,a_wt,        &
         v1dc,v1dd,v1de,2,'adv')

  end subroutine ladvect
  !
  ! ----------------------------------------------------------------------
  ! LADVXU: Advection of U, by U, div is a scratch array, 
  ! tendencies are accumulated if ut variable, dxi is the inverse of 
  ! delta-x the grid spacing in the horizontal direction.
  !
  subroutine ladvxu(n1,n2,n3,uc,ut,flx,dxi)

    integer, intent (in) ::  n1,n2,n3
    real, intent(in)     :: uc(n1,n2,n3),dxi
    real, intent(inout)  :: ut(n1,n2,n3)
    real, intent(out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=3,n2-1
          do k=2,n1-1
             flx(k,i,j)=(0.583333*(uc(k,i,j)+uc(k,i-1,j)) - 0.083333          &
                  *(uc(k,i+1,j)+uc(k,i-2,j)))*.5*(uc(k,i,j)+uc(k,i-1,j))
          end do
       end do

       do i=3,n2-2
          do k=2,n1-1
             ut(k,i,j)=ut(k,i,j)-(flx(k,i+1,j)-flx(k,i,j))*dxi
          end do
       end do
    end do

  end subroutine ladvxu
  !
  ! ----------------------------------------------------------------------
  ! LADVYU: Advection of U, by V, div is a scratch array, 
  ! tendencies are accumulated if fu variable, dyi is the inverse of 
  ! delta-y.
  !
  subroutine ladvyu(n1,n2,n3,uc,ut,vc,flx,dyi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: uc(n1,n2,n3),vc(n1,n2,n3),dyi
    real, intent (in out) :: ut(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=2,n3-2
       do i=1,n2-1
          do k=2,n1-1
             flx(k,i,j)=(0.583333*(uc(k,i,j)+uc(k,i,j+1)) - 0.083333          &
                  *(uc(k,i,j-1)+uc(k,i,j+2)))*.5*(vc(k,i,j)+vc(k,i+1,j))
          end do
       end do
    end do

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             ut(k,i,j)=ut(k,i,j)-(flx(k,i,j)-flx(k,i,j-1))*dyi
          end do
       end do
    end do

  end subroutine ladvyu
  !
  ! ----------------------------------------------------------------------
  ! Subroutine ladvzu: Advection of U, by W*rho, div is a scratch array, 
  ! tendencies are accumulated if ut variable, v1 is the inverse density
  ! array valid at thermo points
  !
  subroutine ladvzu(n1,n2,n3,uc,ut,wm,flx,v1)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: uc(n1,n2,n3),wm(n1,n2,n3),v1(n1)
    real, intent (in out) :: ut(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=1,n2-1
          flx(1,i,j)=0.
          do k=2,n1-2
             flx(k,i,j)=(0.583333*(uc(k,i,j)+uc(k+1,i,j)) - 0.083333          &
                  *(uc(k-1,i,j)+uc(k+2,i,j)))*.5*(wm(k,i,j)+wm(k,i+1,j))
          end do
          flx(n1-1,i,j)=0.
          flx(n1,i,j)=0.
       end do

       do i=1,n2-1
          do k=2,n1-1
             ut(k,i,j)=ut(k,i,j)-(flx(k,i,j)-flx(k-1,i,j))*v1(k)
          end do
       end do
    end do

  end subroutine ladvzu
  !
  ! ----------------------------------------------------------------------
  ! LADVXV: Advection of V, by U, div is a scratch array, 
  ! tendencies are accumulated if vt variable, dxi is the inverse of 
  ! delta-x the grid spacing in the horizontal direction.
  !
  subroutine ladvxv(n1,n2,n3,uc,vc,vt,flx,dxi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: uc(n1,n2,n3),vc(n1,n2,n3),dxi
    real, intent (in out) :: vt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3-1
       do i=2,n2-2
          do k=2,n1-1
             flx(k,i,j)=(0.583333*(vc(k,i,j)+vc(k,i+1,j)) - 0.083333          &
                  *(vc(k,i-1,j)+vc(k,i+2,j)))*.5*(uc(k,i,j)+uc(k,i,j+1))
          end do
       end do

       do i=3,n2-2
          do k=2,n1-1
             vt(k,i,j)=vt(k,i,j)-(flx(k,i,j)-flx(k,i-1,j))*dxi
          end do
       end do
    end do

  end subroutine ladvxv
  !
  ! ----------------------------------------------------------------------
  ! LADVYV: Advection of V, by V, div is a scratch array, 
  ! tendencies are accumulated if vt variable, dyi is the inverse of 
  ! delta-y.
  !
  subroutine ladvyv(n1,n2,n3,vc,vt,flx,dyi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: vc(n1,n2,n3),dyi
    real, intent (in out) :: vt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=3,n3-1
       do i=1,n2
          do k=2,n1-1
             flx(k,i,j)=(0.583333*(vc(k,i,j)+vc(k,i,j-1))-0.083333           &
                  *(vc(k,i,j+1)+vc(k,i,j-2)))*.5*(vc(k,i,j)+vc(k,i,j-1))
          end do
       end do
    end do

    do j=3,n3-2
       do i=1,n2
          do k=2,n1-1
             vt(k,i,j)=vt(k,i,j)-(flx(k,i,j+1)-flx(k,i,j))*dyi
          end do
       end do
    end do

  end subroutine ladvyv
  !
  ! ----------------------------------------------------------------------
  ! LADVZV: Advection of V, by W*rho, div is a scratch array, 
  ! tendencies are accumulated if vt variable, v1 is the inverse density
  ! array valid at thermo points
  !
  subroutine ladvzv(n1,n2,n3,vc,vt,wm,flx,v1)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: vc(n1,n2,n3),wm(n1,n2,n3),v1(n1)
    real, intent (in out) :: vt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3-1
       do i=1,n2
          flx(1,i,j)=0.
          do k=2,n1-2
             flx(k,i,j)=(0.583333*(vc(k,i,j)+vc(k+1,i,j))-0.083333            &
                  *(vc(k-1,i,j)+vc(k+2,i,j)))*.5*(wm(k,i,j)+wm(k,i,j+1))
          end do
          flx(n1-1,i,j)=0.
          flx(n1,i,j)=0.
       end do

       do i=1,n2
          do k=2,n1-1
             vt(k,i,j)=vt(k,i,j)-(flx(k,i,j)-flx(k-1,i,j))*v1(k)
          end do
       end do
    end do

  end subroutine ladvzv
  !
  ! ----------------------------------------------------------------------
  ! LADVXW: Advection of W, by U, div is a scratch array, 
  ! tendencies are accumulated if wt variable, dxi is the inverse of 
  ! delta-x the grid spacing in the horizontal direction.
  !
  subroutine ladvxw(n1,n2,n3,uc,wc,wt,flx,dxi)

    integer, intent (in)  :: n1,n2,n3
    real, intent (in)     :: uc(n1,n2,n3),wc(n1,n2,n3),dxi
    real, intent (in out) :: wt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=2,n2-2
          do k=2,n1-2
             flx(k,i,j)=(0.583333*(wc(k,i,j)+wc(k,i+1,j))-0.083333*          &
                  (wc(k,i-1,j)+wc(k,i+2,j)))*.5*(uc(k,i,j)+uc(k+1,i,j))
          end do
       end do

       do i=3,n2-2
          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)-(flx(k,i,j)-flx(k,i-1,j))*dxi
          end do
       end do
    end do

  end subroutine ladvxw
  !
  ! ----------------------------------------------------------------------
  ! LADVYW: Advection of W, by V, div is a scratch array, 
  ! tendencies are accumulated if wt variable, dyi is the inverse of 
  ! delta-y.
  !
  subroutine ladvyw(n1,n2,n3,vm,wc,wt,flx,dyi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: vm(n1,n2,n3),wc(n1,n2,n3),dyi
    real, intent (in out) :: wt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=2,n3-2
       do i=1,n2
          do k=2,n1-2
             flx(k,i,j)=(0.583333*(wc(k,i,j)+wc(k,i,j+1))-0.083333            &
                  *(wc(k,i,j-1)+wc(k,i,j+2)))*.5*(vm(k,i,j)+vm(k+1,i,j))
          end do
       end do
    end do

    do j=3,n3-2
       do i=1,n2
          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)-(flx(k,i,j)-flx(k,i,j-1))*dyi
          end do
       end do
    end do

  end subroutine ladvyw
  !
  ! ----------------------------------------------------------------------
  ! LADVZW: Advection of W, by W*rho, div is a scratch array, 
  ! tendencies are accumulated if wt variable, v2 is the inverse density
  ! array valid at w-points
  !
  subroutine ladvzw(n1,n2,n3,wc,wt,wm,flx,v2)

    integer, intent (in) ::  n1,n2,n3
    real, intent (in)    :: wm(n1,n2,n3),wc(n1,n2,n3),v2(n1)
    real, intent (inout) :: wt(n1,n2,n3)
    real, intent (out)   :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=1,n2
          flx(1,i,j)=0.
          flx(2,i,j)=wm(2,i,j)*.5*(0.583333*(wc(2,i,j))-0.083333*(wc(3,i,j)))
          do k=3,n1-1
             flx(k,i,j)=(0.583333*(wc(k,i,j)+wc(k-1,i,j))-0.083333            &
                  *(wc(k+1,i,j)+wc(k-2,i,j)))*.5*(wm(k,i,j)+wm(k-1,i,j))
          end do
          flx(n1,i,j)=0.

          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)-(flx(k+1,i,j)-flx(k,i,j))*v2(k)
          end do
       end do
    end do

  end subroutine ladvzw
  !
  ! ----------------------------------------------------------------------
  ! ADVL_PREP: prepares two scratch arrays with the inverse
  ! densities as they locate on thermo levels (v1) and w-levels (v2)
  !
  subroutine advl_prep(n1,n2,n3,w,wm,dn0,dzt,dzm,v1,v2)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    ::  w(n1,n2,n3),dn0(n1),dzt(n1),dzm(n1)
    real, intent (out)   ::  wm(n1,n2,n3),v1(n1),v2(n1)

    integer :: k

    do k=1,n1-1
       wm(k,:,:)=w(k,:,:)*(dn0(k)+dn0(k+1))*.5
       v1(k)=dzt(k)/dn0(k)
       v2(k)=2.*dzm(k)/(dn0(k)+dn0(k+1))
    end do

  end subroutine advl_prep

end module advl
