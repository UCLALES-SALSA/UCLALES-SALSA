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
! Copyright 199-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module lsvar

    ! Large scale forcing parameters
    integer, SAVE :: nt, nz
    REAL, ALLOCATABLE :: t_ls(:), sst_ls(:), ugeo_ls(:,:), vgeo_ls(:,:)

contains
  !----------------------------------------------------------------------
  ! Read the forcing data
  !
  subroutine lsvar_init(nzp,zg_in)
    implicit none
    INTEGER, INTENT(IN) :: nzp
    REAL, INTENT(IN) :: zg_in(nzp) ! Current grid
    INTEGER :: n, i, j, nt_old, nz_old
    REAL :: lf, hf
    REAL, ALLOCATABLE :: zz(:), uu(:,:), vv(:,:)

    ! SST file
    open (NEWUNIT=n,file='input_sst.dat',status='old',form='formatted')
    READ(n,*) nt ! The number of time values
    ALLOCATE(t_ls(nt),sst_ls(nt))
    DO i=1,nt
        READ(n,*) t_ls(i), sst_ls(i) ! Elapsed time (s) and SST (K)
    ENDDO
    CLOSE(n)
    nt_old=nt

    ! u and v files
    open (NEWUNIT=n,file='input_u_forcing.dat',status='old',form='formatted')
    READ(n,*) nz, nt ! The number of altitude and time values
    if (nt_old/=nt) STOP 'Incorrect number of time values!'
    ALLOCATE(zz(nz), uu(nz,nt), vv(nz,nt))
    DO j=1,nz
        read (n,*) zz(j), (uu(j,i),i=1,nt)
    ENDDO
    close (n)
    nz_old=nz

    open (NEWUNIT=n,file='input_v_forcing.dat',status='old',form='formatted')
    READ(n,*) nz, nt ! The number of altitude and time values
    if (nt_old/=nt .OR. nz_old/=nz) STOP 'Incorrect number of time or altitude values!'
    DO j=1,nz
        read (n,*) zz(j), (vv(j,i),i=1,nt)
    ENDDO
    close (n)

    ! Interpolate z-grid
    ALLOCATE(ugeo_ls(nzp,nt),vgeo_ls(nzp,nt))
    DO i=1,nzp
        IF (zg_in(i)<=zz(1)) THEN
            ugeo_ls(i,:)=uu(1,:)
            vgeo_ls(i,:)=vv(1,:)
        ELSEIF (zg_in(i)>=zz(nz)) THEN
            ugeo_ls(i,:)=uu(nz,:)
            vgeo_ls(i,:)=vv(nz,:)
        ELSE
            j=COUNT(zz<zg_in(i))
            hf=(zg_in(i)-zz(j))/(zz(j+1)-zz(j))
            lf=(zz(j+1)-zg_in(i))/(zz(j+1)-zz(j))
            ugeo_ls(i,:)=hf*uu(j+1,:)+lf*uu(j,:)
            vgeo_ls(i,:)=hf*vv(j+1,:)+lf*vv(j,:)
        ENDIF
    ENDDO
    ! Clean
    DEALLOCATE(zz, uu, vv)

    ! The final altitude count
    nz = nzp

  end subroutine lsvar_init
  !
  ! ----------------------------------------------------------------------
  ! subroutine varlscale : computes the time variation of the sst and geostrophic winds
  !
  subroutine varlscale(time_in,sst,u0,v0)
    implicit none
    real, intent (in) :: time_in ! Time after spinup (s)
    !integer, intent (in) :: nz
    real, intent(inout) :: sst
    real, intent(inout) :: u0(*),v0(*)
    real :: lf, hf
    integer :: k

    ! Find the time interval
    k=COUNT(t_ls<time_in)

    IF (k==0) THEN
        ! Use the first SST
        sst=sst_ls(1)
        ! ... and winds
         u0(1:nz) = ugeo_ls(1:nz,1)
         v0(1:nz) = vgeo_ls(1:nz,1)
    ELSEIF (k==nt) THEN
        ! Use the last SST
        sst=sst_ls(nt)
        ! ... and winds
         u0(1:nz) = ugeo_ls(1:nz,nt)
         v0(1:nz) = vgeo_ls(1:nz,nt)
    ELSE
        ! Interpolate between k and k+1
        hf=(time_in-t_ls(k))/(t_ls(k+1)-t_ls(k))
        lf=(t_ls(k+1)-time_in)/(t_ls(k+1)-t_ls(k))
        ! SST
        sst=lf*sst_ls(k)+hf*sst_ls(k+1)
        ! ...and winds
         u0(1:nz) = lf*ugeo_ls(1:nz,k)+hf*ugeo_ls(1:nz,k+1)
         v0(1:nz) = lf*vgeo_ls(1:nz,k)+hf*vgeo_ls(1:nz,k+1)
    ENDIF
   end subroutine varlscale

end module lsvar
