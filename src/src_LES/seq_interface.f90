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
! Copyright 199-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA,
! and TV Singh, Academic and Technology Services
!----------------------------------------------------------------------------
!
MODULE mpi_interface
   !
   !    nxg = nxpg-4
   !    nyg = nypg-4
   !    xcomm, commxid - communicator for x side processors and rank wrt it
   !    ycomm, commyid - communicator for y side processors and rank wrt it
   !    nxnzp = nx*nzp
   !    nynzp = ny*nzp
   !    wrxid, wryid, nxprocs,nyprocs:(wrxid,wryid)=myid
   !       in ranktable (nxprocs,nyprocs)
   !    nxpa,nypa: arrays containing nxp and nyp for all nxprocs and nyprocs resp.
   !    nynza, nxnza: arrays containing nynzp and nxnzp on nxprocs and nyprocs
   !    resp.
   !
   IMPLICIT NONE

   ! Juha: Some interfaces for communication subroutines; Makes it also easier to set corresponding
   !       dummy subroutines in seq_interface.f90 to avoid hairy looking tests in the main source code 
   !       whether MPI is used or not. Another option of course would be to use compiler directives,
   !       maye future work? However it will make the source code again a bit more messy.
   INTERFACE broadcast
      MODULE PROCEDURE broadcastRealArray1d, broadcastRealArray3d, broadcastInteger
   END INTERFACE broadcast
   
   INTEGER :: myid, pecount, nxpg, nypg, nxg, nyg, nbytes, intsize
   INTEGER :: xcomm, ycomm,commxid,commyid, MY_CMPLX, MY_SIZE
   INTEGER :: nxnzp,nynzp
   INTEGER :: wrxid, wryid, nxprocs, nyprocs
   INTEGER, ALLOCATABLE, DIMENSION(:) :: xoffset, yoffset, nxpa, nypa
   CHARACTER(len=80) :: ver='', author=''
   ! Additional, e.g. case specific, information
   CHARACTER(len=180), PARAMETER :: info=''

  ! these are the parameters used in the alltoallw call in the fft

CONTAINS
   !
   !----------------------------------------------------------------------
   ! INIT_MP: Initializes MPI
   !
   SUBROUTINE init_mpi

      CHARACTER (len=8) date

      myid = 0
      pecount = 1

      SELECT CASE (kind(0.0))
         CASE (4)
            nbytes = 4
            MY_SIZE = 4
            MY_CMPLX = 8
         CASE (8)
            nbytes = 8
            MY_SIZE = 8
            MY_CMPLX = 16
         CASE DEFAULT
            STOP "kind not supported"
      END SELECT

      SELECT CASE(kind(0))
         CASE (4)
            intsize = 4
         CASE (8)
            intsize = 8
         CASE DEFAULT
            STOP "int kind not supported"
      END SELECT
      !
      CALL date_and_time(date)
      IF (myid == 0) PRINT "(/1x,75('-'),/2x,A22,/2x,A15,I2,A15,I2,A14)", &
         'UCLALES-SALSA '//date, 'Computing using',nbytes,' byte REALs and', &
         intsize," byte INTEGERs"
      IF (myid == 0 .AND. len(info) > 0) PRINT *, ' '//trim(info)

   END SUBROUTINE init_mpi
   !
   !----------------------------------------------------------------------
   ! DEFINE_DECOMP: Defines MPI Decomposition
   !
   SUBROUTINE define_decomp(nxp, nyp, nxpart)

      INTEGER, INTENT(inout) :: nxp, nyp
      LOGICAL, INTENT(in)    :: nxpart

      nxprocs = 1
      nyprocs = 1

      !
      !   ranktable is the matrix having ranks of processes in x-y domain
      !
      wrxid = 0
      wryid = 0
      commxid = 0
      commyid = 0

      !
      ! there are two boundary points in each direction
      !
      nxpg = nxp
      nypg = nyp

      nxg = nxpg-4
      nyg = nypg-4

      ALLOCATE (nxpa(0:nxprocs-1), nypa(0:nyprocs-1))
      nxpa(0) = nxg
      nypa(0) = nyg

      !
      !  offsets for ecah processor in x and y, for a given grid (nxp x nyp)
      !
      ALLOCATE(xoffset(0:nxprocs-1),yoffset(0:nyprocs-1))

      xoffset = 0
      yoffset = 0

      IF(nxp < 5) THEN
         PRINT *, 'ABORT: X Horizontal domain size too small for ',nxprocs,    &
            ' processors.'
         PRINT *, '       Increase nyp to ',nxprocs*5, ' or run on ',nxpg/5,   &
            ' or fewer processors'
         CALL appl_abort(0)
      END IF
      IF(nyp < 5) THEN
         PRINT *, 'ABORT: Y Horizontal domain size too small for ',nyprocs,    &
            ' processors.'
         PRINT *, '       Increase nyp to ',nyprocs*5, ' or run on ',nypg/5,   &
            ' or fewer processors'
         CALL appl_abort(0)
      END IF

      IF (myid == 0) PRINT 61,'Sequential (shared memory) simulation'

61    FORMAT (/1x,49('-')/2x,A37)

   END SUBROUTINE define_decomp
   !
   !----------------------------------------------------------------------
   ! INIT_ALLTOALL_REORDERXY: Defines the mpi derived types to do a data
   ! movement of the form A(m,n/p,z) -> B(n,m/p,z) for data of type MY_CMPLX
   !
   SUBROUTINE init_alltoall_reorder(nxp,nyp,nzp)

      INTEGER, INTENT(in) :: nxp,nyp,nzp

      nxg = nxp-4
      nyg = nyp-4
      nxnzp = nxg*nzp
      nynzp = nyg*nzp

   END SUBROUTINE init_alltoall_reorder
   ! ---------------------------------------------------------------------
   ! Subroutine cyclics: commits exchange cyclic x boundary conditions
   !
   SUBROUTINE cyclics(n1,n2,n3,var,req)

      INTEGER, INTENT(in) :: n1,n2,n3,req(16)
      REAL, INTENT(inout) :: var(n1,n2,n3)

      IF (n3 == 5) THEN
         var(:,:,1) = var(:,:,3)
         var(:,:,2) = var(:,:,3)
         var(:,:,4) = var(:,:,3)
         var(:,:,5) = var(:,:,3)
      END IF
      IF (n2 == 5) THEN
         var(:,1,:) = var(:,3,:)
         var(:,2,:) = var(:,3,:)
         var(:,4,:) = var(:,3,:)
         var(:,5,:) = var(:,3,:)
      END IF

      var(:,:2,:) = var(:,n2-3:n2-2,:)
      var(:,n2-1:,:) = var(:,3:4,:)

      var(:,:,:2) = var(:,:,n3-3:n3-2)
      var(:,:,n3-1:) = var(:,:,3:4)
      var(:,:2,:2) = var(:,n2-3:n2-2,n3-3:n3-2)
      var(:,n2-1:,n3-1:) = var(:,3:4,3:4)
    
   END SUBROUTINE cyclics
   !
   ! ---------------------------------------------------------------------
   ! Subroutine cyclicc: comits excahnging cyclic boundary conditions
   SUBROUTINE cyclicc(n1,n2,n3,var,req)

      INTEGER :: n1,n2,n3,req(16)
      REAL    :: var(n1,n2,n3)

   END SUBROUTINE cyclicc
   !
   ! ---------------------------------------------------------------------
   SUBROUTINE appl_abort(ierr)

      INTEGER :: ierr
      STOP 'Program Aborted'

   END SUBROUTINE appl_abort
   !
   ! ---------------------------------------------------------------------
   SUBROUTINE appl_finalize(ierr)

      INTEGER :: ierr

   END SUBROUTINE appl_finalize
   !
   !---------------------------------------------------------------------------
   SUBROUTINE xshuffle(a,atmp,nx,ny,nz,isign)

      INTEGER, INTENT(in)    :: nx,ny,nz,isign
      COMPLEX, INTENT(inout) :: a(nx,ny,nz),atmp((nx+1)*(ny+1)*(nz+1))
      INTEGER ll,i,j,k

      IF(isign == 1) THEN
         ll = 0
         DO k = 1, nz
            DO j = 1, ny
               DO i = 1, nx
                  ll = ll+1
                  atmp(ll) = a(i,j,k)
               END DO
            END DO
         END DO

      ELSE
         ll = 0
         DO k = 1, nz
            DO j = 1, ny
               DO i = 1, nx
                  ll = ll+1
                  a(i,j,k) = atmp(ll)
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE xshuffle
   !
   !---------------------------------------------------------------------------
   SUBROUTINE yshuffle(a,atmp,nx,ny,nz,isign)

      INTEGER, INTENT(in)    :: nx,ny,nz,isign
      COMPLEX, INTENT(inout) :: a(ny,nx,nz),atmp((nx+1)*(ny+1)*(nz+1))
      INTEGER ll,i,j,k

      IF(isign == 1) THEN
         ll = 0
         DO k = 1, nz
            DO j = 1, nx
               DO i = 1, ny
                  ll = ll+1
                  atmp(ll) = a(i,j,k)
               END DO
            END DO
         END DO
      ELSE
         ll = 0
         DO k = 1, nz
            DO j = 1, nx
               DO i = 1, ny
                  ll = ll+1
                  a(i,j,k) = atmp(ll)
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE yshuffle
   !
   !---------------------------------------------------------------------------
   ! get maximum across processors
   !
   SUBROUTINE double_scalar_par_max(xxl,xxg)

      REAL(kind=8), INTENT(out) :: xxg
      REAL(kind=8), INTENT(in)  :: xxl

      xxg = xxl

   END SUBROUTINE double_scalar_par_max
   !
   !---------------------------------------------------------------------------
   SUBROUTINE double_scalar_par_sum(xxl,xxg)

      REAL(kind=8), INTENT(out) :: xxg
      REAL(kind=8), INTENT(in)  :: xxl

      xxg = xxl

   END SUBROUTINE double_scalar_par_sum
   !
   !---------------------------------------------------------------------------
   SUBROUTINE double_array_par_sum(xxl,xxg,n)

      INTEGER, INTENT(in) :: n
      REAL(kind=8), INTENT(out) :: xxg(n)
      REAL(kind=8), INTENT(in)  :: xxl(n)

      xxg = xxl

   END SUBROUTINE double_array_par_sum

 ! Juha added: Broadcast real arrays and stuff (Need interface to cover everything!)
 SUBROUTINE broadcastRealArray1d(NNdims,rootid,sendbuff)
   IMPLICIT NONE
   INTEGER, INTENT(in) :: NNdims(1)
   INTEGER, INTENT(in) :: rootid
   REAL, INTENT(inout) :: sendbuff(NNdims(1))
 END SUBROUTINE 
 
 SUBROUTINE broadcastRealArray3d(NNdims,rootid,sendbuff)
   IMPLICIT NONE
   INTEGER, INTENT(in) :: NNdims(3)
   INTEGER, INTENT(in) :: rootid
   REAL, INTENT(inout) :: sendbuff(NNdims(1),NNdims(2),NNdims(3))
 END SUBROUTINE broadcastRealArray3d

 SUBROUTINE broadcastInteger(sendbuff,rootid)
   IMPLICIT NONE
   INTEGER, INTENT(inout) :: sendbuff
   INTEGER, INTENT(in) :: rootid
 END SUBROUTINE broadcastInteger

END MODULE mpi_interface
