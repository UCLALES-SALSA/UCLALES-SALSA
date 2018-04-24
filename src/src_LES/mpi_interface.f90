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

   USE mpi
   IMPLICIT NONE
   
   ! Juha: Some interfaces for communication subroutines; Makes it also easier to set corresponding
   !       dummy subroutines in seq_interface.f90 to avoid hairy looking tests in the main source code 
   !       whether MPI is used or not. Another option of course would be to use compiler directives,
   !       maye future work? However it will make the source code again a bit more messy.
   INTERFACE broadcast
      MODULE PROCEDURE broadcastRealArray1d, broadcastRealArray3d, broadcastInteger
   END INTERFACE broadcast
   
   !
   !    nxg = nxpg-4
   !    nyg = nypg-4
   !    xcomm, commxid - communicator for x side processors and rank wrt it
   !    ycomm, commyid - communicator for y side processors and rank wrt it
   !    nxnzp = nx*nzp
   !    nynzp = ny*nzp
   !    wrxid, wryid, nxprocs,nyprocs: (wrxid,wryid)=myid in
   !        ranktable (nxprocs,nyprocs)
   !    nxpa,nypa: arrays containing nxp and nyp for all nxprocs and nyprocs
   !         respectively
   !    nynza, nxnza: arrays containing nynzp and nxnzp on nxprocs and nyprocs
   !         respectively
   !

   INTEGER :: myid, pecount, nxpg, nypg, nxg, nyg, nbytes, intsize, &
              MY_SIZE, MY_CMPLX
   INTEGER :: xcomm, ycomm,commxid,commyid
   INTEGER :: nxnzp,nynzp,fftinix,fftiniy
   INTEGER :: wrxid, wryid, nxprocs, nyprocs
   INTEGER, ALLOCATABLE, DIMENSION(:) :: xoffset, yoffset, nxpa, nypa, &
                                         nynza, nxnza

   ! these are the parameters used in the alltoallw call in the fft

   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ranktable,xtype,ytype,xdisp,&
                                           ydisp,xcount,ycount

   INTEGER :: stridetype,xstride,ystride,xystride,xylarry,xyzlarry,&
              fxytype,fxyztype

   CHARACTER(len=80) :: ver='', author=''
   ! Additional, e.g. case specific, information
   CHARACTER(len=180), PARAMETER :: info=''

CONTAINS
   !
   !----------------------------------------------------------------------
   ! INIT_MP: Initializes MPI
   !
   SUBROUTINE init_mpi

      INTEGER :: ierror
      CHARACTER (len=8) :: date

      CALL mpi_init(ierror)
      CALL mpi_comm_size(MPI_COMM_WORLD, pecount, ierror)
      CALL mpi_comm_rank(MPI_COMM_WORLD, myid, ierror)

      SELECT CASE (kind(0.0))
         CASE (4)
            nbytes = 4
            MY_SIZE = MPI_REAL
            MY_CMPLX = MPI_COMPLEX
         CASE (8)
            nbytes = 8
            MY_SIZE = MPI_DOUBLE_PRECISION
            MY_CMPLX = MPI_doUBLE_COMPLEX
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

      INTEGER :: ierror, i,j, modx,mody, nxpj, nypj, irank
      INTEGER :: worldgroup, xgroup, ygroup

      nxprocs = 1
      nyprocs = pecount
      IF (nyp > nxp) THEN
         nyprocs = pecount
         nxprocs = 1
         IF( (nyp-4)/nyprocs < 5 .OR. nxpart) THEN
            nyprocs = int(sqrt(REAL(pecount)))
            nxprocs = pecount/nyprocs
            DO WHILE (nyprocs*nxprocs /= pecount)
               nyprocs = nyprocs+1
               nxprocs = pecount/nyprocs
            END DO

            IF(nxprocs > nyprocs)then
               i = nxprocs
               nxprocs = nyprocs
               nyprocs = i
            END IF

            IF(nxp < 5) THEN
               PRINT *, 'ABORTING: NXP too small, increase to at least 5'
               CALL mpi_abort(MPI_COMM_WORLD,0,ierror)
            ELSE IF( (nxp-4)/nxprocs <= 5) THEN
               nxprocs = 1
               nyprocs = pecount
            END IF

            IF ( (nyp-4)/nyprocs < 5) THEN
               PRINT *, '  ABORTING: NYP too small for ',nyprocs,' processors.'
               PRINT *, '  Increase to ',nyprocs*9, ' or run on ',nypg/9,       &
                  ' or fewer processors'
               CALL mpi_abort(MPI_COMM_WORLD,0,ierror)
            END IF
         END IF

      ELSE

         nxprocs = pecount
         nyprocs = 1
         IF( (nxp-4)/nxprocs < 5 .OR. nxpart) THEN
            nxprocs = int(sqrt(REAL(pecount)))
            nyprocs = pecount/nxprocs

            DO WHILE (nyprocs*nxprocs /= pecount)
               nxprocs = nxprocs+1
               nyprocs = pecount/nxprocs
            END DO

            IF(nyprocs > nxprocs)then
               i = nyprocs
               nyprocs = nxprocs
               nxprocs = i
            END IF

            IF(nyp < 5) THEN
               PRINT *, 'ABORTING: NYP too small, increase to at least 5'
               CALL mpi_abort(MPI_COMM_WORLD,0,ierror)
            ELSE IF( (nyp-4)/nyprocs <= 5) THEN
               nyprocs = 1
               nxprocs = pecount
            END IF

            IF ( (nxp-4)/nxprocs < 5) THEN
               PRINT *, '  ABORTING: NXP too small for ',nxprocs,' processors.'
               PRINT *, '  Increase to ',nxprocs*9, ' or run on ',nxpg/9,       &
                  ' or fewer processors'
               CALL mpi_abort(MPI_COMM_WORLD,0,ierror)
            END IF
         END IF
      END IF

      IF( (nyp-4)/nyprocs < 5 .OR. nxpart) THEN
         nyprocs = int(sqrt(REAL(pecount)))
         nxprocs = pecount/nyprocs
         DO WHILE (nyprocs*nxprocs /= pecount)
            nyprocs = nyprocs+1
            nxprocs = pecount/nyprocs
         END DO
         IF(nxprocs > nyprocs)then
            i = nxprocs
            nxprocs = nyprocs
            nyprocs = i
         END IF
      END IF

      !ctvs    nxprocs=1
      !ctvs    nyprocs=2
      !
      !   ranktable is the matrix having ranks of processes in x-y domain
      !

      ALLOCATE(ranktable(-1:nxprocs,-1:nyprocs))

      DO i = -1, nxprocs
         DO j = -1, nyprocs
            ranktable(i,j) = 0
         END DO
      END DO
      irank = 0
      DO j = 0, nyprocs-1
         DO i = 0, nxprocs-1
            ranktable(i,j) = irank
            IF (myid == irank) THEN
               wrxid = i
               wryid = j
            END IF
            irank = irank + 1
         END DO
      END DO
      DO i = 0, nxprocs-1
         ranktable(i, -1) = ranktable(i, nyprocs - 1)
         ranktable(i, nyprocs) = ranktable(i, 0)
      END DO
      DO j = 0, nyprocs-1
         ranktable(-1, j) = ranktable(nxprocs - 1, j)
         ranktable(nxprocs, j) = ranktable(0, j)
      END DO
      ranktable(-1,-1) = ranktable(nxprocs-1,nyprocs-1)
      ranktable(nxprocs,nyprocs) = ranktable(0,0)
      ranktable(-1,nyprocs) = ranktable(nxprocs-1,0)
      ranktable(nxprocs,-1) = ranktable(0,nyprocs-1)

      CALL mpi_comm_group(mpi_comm_world, worldgroup, ierror)
      CALL mpi_group_incl(worldgroup, nxprocs,ranktable(0:nxprocs-1,wryid),&
                          xgroup,ierror)
      CALL mpi_comm_create(mpi_comm_world, xgroup,xcomm,ierror)
      CALL mpi_group_incl(worldgroup, nyprocs,ranktable(wrxid,0:nyprocs-1),&
                          ygroup,ierror)
      CALL mpi_comm_create(mpi_comm_world, ygroup,ycomm,ierror)

      CALL mpi_comm_rank(xcomm,commxid,ierror)
      CALL mpi_comm_rank(ycomm,commyid,ierror)

      !
      ! there are two boundary points in each direction
      !
      nxpg = nxp
      nypg = nyp
      nxp = (nxpg-4)/nxprocs + 4
      nyp = (nypg-4)/nyprocs + 4

      modx = modulo(nxpg-4,nxprocs)
      mody = modulo(nypg-4,nyprocs)
      !
      ! offsets for each processor in x and y directons(nxp x nyp)
      !
      ALLOCATE(xoffset(0:nxprocs-1),yoffset(0:nyprocs-1))

      xoffset = 0

      DO j = 1, nxprocs-1
         IF(j <= modx) THEN
            nxpj = nxp+1
         ELSE
            nxpj = nxp
         END IF
         xoffset(j) = xoffset(j-1)+nxpj-4
      END DO

      yoffset = 0

      DO j = 1, nyprocs-1
         IF(j <= mody) THEN
            nypj = nyp+1
         ELSE
            nypj = nyp
         END IF
         yoffset(j) = yoffset(j-1)+nypj-4
      END DO

      IF (nxpg > 5 .AND. nxp == 5) THEN
         PRINT *, 'ABORTING: Subdomain too finely discretized in x', nxpg, nxp
         CALL mpi_abort(MPI_COMM_WORLD,0,ierror)
      END IF

      IF (nypg > 5 .AND. nyp == 5) THEN
         PRINT *, 'ABORTING: Subdomain too finely discretized in y', nypg, nyp
         CALL mpi_abort(MPI_COMM_WORLD,0,ierror)
      END IF

      IF(nxp < 5) THEN
         PRINT *, 'ABORT: X Horizontal domain size too small for ',nxprocs,    &
            ' processors.'
         PRINT *, '       Increase nyp to ',nxprocs*5, ' or run on ',nxpg/5, &
            ' or fewer processors'
         CALL mpi_abort(MPI_COMM_WORLD,0,ierror)
      END IF
      IF(nyp < 5) THEN
         PRINT *, 'ABORT: Y Horizontal domain size too small for ',nyprocs,    &
            ' processors.'
         PRINT *, '       Increase nyp to ',nyprocs*5, ' or run on ',nypg/5, &
            ' or fewer processors'
         CALL mpi_abort(MPI_COMM_WORLD,0,ierror)
      END IF

      IF (myid == 0) THEN
         PRINT 61, 'Processor count', pecount,'nxpl =', nxp,' nypl = ',nyp

         DO i = 0, min(nxprocs,nyprocs)-1
            PRINT "(2x,A13,2I5)", 'x/y offset = ', xoffset(i), yoffset(i)
         END DO

         IF (nxprocs > nyprocs) PRINT "(15x,I5)", xoffset(nyprocs:nxprocs-1)
         IF (nxprocs < nyprocs) PRINT "(15x,I5)", yoffset(nxprocs:nyprocs-1)
      END IF

61    FORMAT (/1x,49('-')/2x,A15,I5,2(A6,I5))

   END SUBROUTINE define_decomp
   !
   !----------------------------------------------------------------------
   ! INIT_ALLTOALL_REORDERXY: Defines the mpi derived types to do a data
   ! movement of the form A(m,n/p,z) -> B(n,m/p,z) for data of type MY_CMPLX
   !
   SUBROUTINE init_alltoall_reorder(nxp,nyp,nzp)

      INTEGER, INTENT(in) :: nxp,nyp,nzp

      INTEGER :: nx, ny, i, j, k, ii, jj, ierr, cnt, typesize,nynzg, nxnzg

      nx  = max(1,nxp-4)
      ny  = max(1,nyp-4)
      nxg = nxpg-4
      nyg = nypg-4

      ALLOCATE(nxpa(0:nxprocs-1), nypa(0:nyprocs-1))
      ALLOCATE(nxnza(0:nyprocs-1), nynza(0:nxprocs-1))
      ALLOCATE(xcount(0:nxprocs-1,2),xtype(0:nxprocs-1,2),xdisp(0:nxprocs-1,2))
      ALLOCATE(ycount(0:nyprocs-1,2),ytype(0:nyprocs-1,2),ydisp(0:nyprocs-1,2))

      ii = nxg/nxprocs
      ii = nxg - nxprocs*ii

      DO i = 0, nxprocs-1
         nxpa(i) = nxg/nxprocs
         IF(i < ii) nxpa(i) = nxpa(i)+1
      END DO

      jj = nyg/nyprocs
      jj = nyg - nyprocs*jj

      DO i = 0, nyprocs-1
         nypa(i) = nyg/nyprocs
         IF(i < jj) nypa(i) = nypa(i)+1
      END DO

      nynzg = ny*nzp
      ii = nynzg/nxprocs
      ii = nynzg - nxprocs*ii
      DO i = 0, nxprocs-1
         nynza(i) = nynzg/nxprocs
         IF (i < ii) nynza(i) = nynza(i)+1
      END DO

      nxnzg = nx*nzp
      jj = nxnzg/nyprocs
      jj = nxnzg - nyprocs*jj
      DO i = 0, nyprocs-1
         nxnza(i) = nxnzg/nyprocs
         IF (i < jj) nxnza(i) = nxnza(i)+1
      END DO

      nxnzp = nxnza(wryid)
      nynzp = nynza(wrxid)

      CALL MPI_TYPE_SIZE(MY_CMPLX,typesize, ierr)

      xcount = 1
      ycount = 1
      xdisp = 0
      ydisp = 0

      DO i = 1, nxprocs-1
         xdisp(i,1) = xdisp(i-1,1)+nx*nynza(i-1)*typesize
         xdisp(i,2) = xdisp(i-1,2)+nxpa(i-1)*typesize
      END DO

      DO i = 1, nyprocs-1
         ydisp(i,1) = ydisp(i-1,1)+ny*nxnza(i-1)*typesize
         ydisp(i,2) = ydisp(i-1,2)+nypa(i-1)*typesize
      END DO

      DO i = 0, nxprocs-1
         CALL mpi_type_contiguous(nx*nynza(i),MY_CMPLX, xtype(i,1),ierr)
         CALL mpi_type_commit(xtype(i,1),ierr)
         CALL mpi_type_vector(nynzp, nxpa(i), nxg, MY_CMPLX, xtype(i,2),ierr)
         CALL mpi_type_commit(xtype(i,2),ierr)
      END DO

      DO i = 0, nyprocs-1
         CALL mpi_type_contiguous(ny*nxnza(i),MY_CMPLX, ytype(i,1),ierr)
         CALL mpi_type_commit(ytype(i,1),ierr)
         CALL mpi_type_vector(nxnzp, nypa(i), nyg, MY_CMPLX, ytype(i,2),ierr)
         CALL mpi_type_commit(ytype(i,2),ierr)
      END DO

      CALL MPI_TYPE_VECTOR(nyp-4,nzp*2,nxp*nzp,MY_SIZE,stridetype,ierr)
      CALL MPI_TYPE_COMMIT(stridetype,ierr)
      CALL MPI_TYPE_VECTOR(nyp-4,nzp*2,nxp*nzp,MY_SIZE,xstride,ierr)
      CALL MPI_TYPE_COMMIT(xstride,ierr)
      CALL MPI_TYPE_VECTOR(2,nzp*(nxp-4),nxp*nzp,MY_SIZE,ystride,ierr)
      CALL MPI_TYPE_COMMIT(ystride,ierr)
      CALL MPI_TYPE_VECTOR(2,2*nzp,nxp*nzp,MY_SIZE,xystride,ierr)
      CALL MPI_TYPE_COMMIT(xystride,ierr)

      CALL MPI_TYPE_VECTOR(nyp-4,nxp-4,nxpg-4,MY_SIZE,fxytype,ierr)
      CALL MPI_TYPE_COMMIT(fxytype,ierr)

      CALL MPI_TYPE_VECTOR(nyp-4,(nxp-4)*nzp,(nxpg-4)*nzp,MY_SIZE,fxyztype,ierr)
      CALL MPI_TYPE_COMMIT(fxyztype,ierr)

      CALL MPI_TYPE_VECTOR(nyp-4,nxp-4,nxp,MY_SIZE,xylarry,ierr)
      CALL MPI_TYPE_COMMIT(xylarry,ierr)
      CALL MPI_TYPE_VECTOR(nyp-4,(nxp-4)*nzp,nxp*nzp,MY_SIZE,xyzlarry,ierr)
      CALL MPI_TYPE_COMMIT(xyzlarry,ierr)

   END SUBROUTINE init_alltoall_reorder

   ! ---------------------------------------------------------------------
   ! Subroutine cyclics: commits exchange cyclic x boundary conditions
   !
   SUBROUTINE cyclics(n1,n2,n3,var,req)

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(inout) :: var(n1,n2,n3)
      INTEGER :: req(16)
      INTEGER :: ierror, stats(MPI_STATUS_SIZE,16), pxfwd, pxback, pyfwd, pyback
      INTEGER :: pxyne,pxyse,pxynw,pxysw

      IF (nypg == 5) THEN
         var(:,:,1) = var(:,:,3)
         var(:,:,2) = var(:,:,3)
         var(:,:,4) = var(:,:,3)
         var(:,:,5) = var(:,:,3)
      END IF
      IF (nxpg == 5) THEN
         var(:,1,:) = var(:,3,:)
         var(:,2,:) = var(:,3,:)
         var(:,4,:) = var(:,3,:)
         var(:,5,:) = var(:,3,:)
      END IF

      pxfwd  = ranktable(wrxid+1,wryid)
      pxback = ranktable(wrxid-1,wryid)
      pyfwd  = ranktable(wrxid,wryid+1)
      pyback = ranktable(wrxid,wryid-1)

      pxyne = ranktable(wrxid+1,wryid+1)
      pxyse = ranktable(wrxid+1,wryid-1)
      pxynw = ranktable(wrxid-1,wryid+1)
      pxysw = ranktable(wrxid-1,wryid-1)

      CALL mpi_isend(var(1,n2-3,3),1,xstride, pxfwd, 130, &
                     MPI_COMM_WORLD, req(1), ierror)
      CALL mpi_isend(var(1,3,3), 1, xstride, pxback, 140, &
                     MPI_COMM_WORLD, req(2), ierror)
      CALL mpi_irecv(var(1,n2-1,3), 1, xstride, pxfwd, 140, &
                     MPI_COMM_WORLD, req(3), ierror)
      CALL mpi_irecv(var(1,1,3), 1, xstride, pxback, 130, &
                     MPI_COMM_WORLD, req(4), ierror)

      CALL mpi_isend(var(1,3,3), 1, ystride, pyback, 110, &
                     MPI_COMM_WORLD, req(5), ierror)
      CALL mpi_irecv(var(1,3,n3-1), 1, ystride, pyfwd, 110, &
                     MPI_COMM_WORLD, req(6), ierror)
      CALL mpi_isend(var(1,3,n3-3), 1, ystride, pyfwd, 120, &
                     MPI_COMM_WORLD, req(7), ierror)
      CALL mpi_irecv(var(1,3,1), 1, ystride, pyback, 120, &
                     MPI_COMM_WORLD, req(8), ierror)

      CALL mpi_isend(var(1,n2-3,n3-3), 1, xystride, pxyne, 150, &
                     MPI_COMM_WORLD, req(9), ierror)
      CALL mpi_irecv(var(1,1,1), 1, xystride, pxysw, 150, &
                     MPI_COMM_WORLD, req(10), ierror)
      CALL mpi_isend(var(1,n2-3,3), 1, xystride, pxyse, 160, &
                     MPI_COMM_WORLD, req(11), ierror)
      CALL mpi_irecv(var(1,1,n3-1), 1, xystride, pxynw, 160, &
                     MPI_COMM_WORLD, req(12), ierror)

      CALL mpi_isend(var(1,3,n3-3), 1, xystride, pxynw, 170, &
                     MPI_COMM_WORLD, req(13), ierror)
      CALL mpi_irecv(var(1,n2-1,1), 1, xystride, pxyse, 170, &
                     MPI_COMM_WORLD, req(14), ierror)
      CALL mpi_isend(var(1,3,3), 1, xystride, pxysw, 180, &
                     MPI_COMM_WORLD, req(15), ierror)
      CALL mpi_irecv(var(1,n2-1,n3-1), 1, xystride, pxyne, 180, &
                     MPI_COMM_WORLD, req(16), ierror)

   END SUBROUTINE cyclics
   !
   !
   ! ---------------------------------------------------------------------
   ! Subroutine cyclicc: comits excahnging cyclic boundary conditions
   SUBROUTINE cyclicc(n1,n2,n3,var,req)

      INTEGER :: ierror, stats(MPI_STATUS_SIZE,16)
      INTEGER :: req(16),n1,n2,n3
      REAL    :: var(n1,n2,n3)

      CALL mpi_waitall(16,req,stats,ierror)

   END SUBROUTINE cyclicc

   SUBROUTINE appl_abort(apperr)

      INTEGER :: apperr,ierr
      CALL mpi_abort(MPI_COMM_WORLD,apperr,ierr)

   END SUBROUTINE appl_abort

   SUBROUTINE appl_finalize(ierr)

      INTEGER :: ierr
      CALL mpi_finalize(ierr)

   END SUBROUTINE appl_finalize

   SUBROUTINE xshuffle(a,atmp,nx,ny,nz,isign)

      INTEGER, INTENT(in)    :: nx,ny,nz,isign
      COMPLEX, INTENT(inout) :: a(nx,ny,nz),atmp((nx+1)*(ny+1)*(nz+1))
      INTEGER :: ierr,ll,i,j,k

      IF(isign == 1) THEN
         IF(nxprocs /= 1)then
            CALL mpi_alltoallw( a,xcount(0:,1) , xdisp(0:,1), xtype(0:,1), atmp, &
                                xcount(0:,2), xdisp(0:,2),xtype(0:,2),xcomm,ierr)
         ELSE
            ll = 0
            DO k = 1, nz
               DO j = 1, ny
                  DO i = 1, nx
                     ll = ll+1
                     atmp(ll) = a(i,j,k)
                  END DO
               END DO
            END DO

         END IF
      ELSE
         IF(nxprocs /= 1)then
            CALL mpi_alltoallw(atmp,xcount(0:,2),xdisp(0:,2),xtype(0:,2),a, &
                               xcount(0:,1), xdisp(0:,1),xtype(0:,1),xcomm,ierr)
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
      END IF

   END SUBROUTINE xshuffle

   SUBROUTINE yshuffle(a,atmp,nx,ny,nz,isign)

      INTEGER, INTENT(in)    :: nx,ny,nz,isign
      COMPLEX, INTENT(inout) :: a(ny,nx,nz),atmp((nx+1)*(ny+1)*(nz+1))
      INTEGER :: ierr,ll,i,j,k

      IF(isign == 1) THEN
         IF(nyprocs /= 1)then
            CALL mpi_alltoallw( a,ycount(0:,1),ydisp(0:,1),ytype(0:,1),atmp, &
                                ycount(0:,2),ydisp(0:,2),ytype(0:,2),ycomm,ierr)
         ELSE
            ll = 0
            DO k = 1, nz
               DO j = 1, ny
                  DO i = 1, nx
                     ll = ll+1
                     atmp(ll) = a(j,i,k)  ! Fixed i & j
                  END DO
               END DO
            END DO
         END IF
      ELSE
         IF(nyprocs /= 1)then
            CALL mpi_alltoallw(atmp,ycount(0:,2),ydisp(0:,2),ytype(0:,2),a, &
               ycount(0:,1),ydisp(0:,1),ytype(0:,1),ycomm,ierr)
         ELSE
            ll = 0
            DO k = 1, nz
               DO j = 1, ny
                  DO i = 1, nx
                     ll = ll+1
                     a(j,i,k) = atmp(ll)
                  END DO
               END DO
            END DO

         END IF
      END IF

   END SUBROUTINE yshuffle
   !
   !---------------------------------------------------------------------------
   ! get maximum across processors
   !
   SUBROUTINE double_scalar_par_max(xxl,xxg)

      REAL(kind=8), INTENT(out) :: xxg
      REAL(kind=8), INTENT(in)  :: xxl
      INTEGER :: mpiop,ierror


      CALL mpi_allreduce(xxl,xxg,1,MPI_DOUBLE_PRECISION, MPI_MAX, &
                         MPI_COMM_WORLD, ierror)

   END SUBROUTINE double_scalar_par_max


   SUBROUTINE double_scalar_par_sum(xxl,xxg)

      REAL(kind=8), INTENT(out) :: xxg
      REAL(kind=8), INTENT(in)  :: xxl
      INTEGER :: mpiop,ierror


      CALL mpi_allreduce(xxl,xxg,1,MPI_DOUBLE_PRECISION, MPI_SUM, &
                         MPI_COMM_WORLD, ierror)

   END SUBROUTINE double_scalar_par_sum

   SUBROUTINE double_array_par_sum(xxl,xxg,n)

      INTEGER, INTENT(in) :: n
      REAL(kind=8), INTENT(out) :: xxg(n)
      REAL(kind=8), INTENT(in)  :: xxl(n)
      INTEGER :: mpiop,ierror


      CALL mpi_allreduce(xxl,xxg,n,MPI_DOUBLE_PRECISION, MPI_SUM, &
                         MPI_COMM_WORLD, ierror)

   END SUBROUTINE double_array_par_sum

 ! Juha added: Broadcast real arrays and stuff (Need interface to cover everything!)
 SUBROUTINE broadcastRealArray1d(NNdims,rootid,sendbuff)
   IMPLICIT NONE
   
   INTEGER, INTENT(in) :: NNdims(1)
   INTEGER, INTENT(in) :: rootid
   REAL, INTENT(inout) :: sendbuff(NNdims(1))
   INTEGER :: ierror

   CALL MPI_BCAST(sendbuff,NNdims(1),MPI_DOUBLE_PRECISION,rootid,MPI_COMM_WORLD,ierror)

 END SUBROUTINE broadcastRealArray1d
 
 SUBROUTINE broadcastRealArray3d(NNdims,rootid,sendbuff)
   IMPLICIT NONE

   INTEGER, INTENT(in) :: NNdims(3)
   INTEGER, INTENT(in) :: rootid
   REAL, INTENT(inout) :: sendbuff(NNdims(1),NNdims(2),NNdims(3))
   INTEGER :: NNtot
   INTEGER :: ierror

   REAL, ALLOCATABLE :: buff1d(:)

   NNtot = NNdims(1)*NNdims(2)*NNdims(3)
   IF (myid /= rootid) sendbuff = 0.

   ALLOCATE(buff1d(NNtot))
   IF (myid == 0) THEN
      buff1d = RESHAPE(sendbuff,(/NNtot/))
   ELSE
      buff1d = 0.
   END IF

   CALL MPI_BCAST(buff1d,NNtot,MPI_DOUBLE_PRECISION,rootid,MPI_COMM_WORLD,ierror)
  
   sendbuff = RESHAPE(buff1d,NNdims)
   DEALLOCATE(buff1d)


 END SUBROUTINE broadcastRealArray3d

 SUBROUTINE broadcastInteger(sendbuff,rootid)
   IMPLICIT NONE

   INTEGER, INTENT(inout) :: sendbuff
   INTEGER, INTENT(in) :: rootid
   INTEGER :: ierror

   CALL MPI_BCAST(sendbuff,1,MPI_INTEGER,rootid,MPI_COMM_WORLD,ierror)


 END SUBROUTINE broadcastInteger




END MODULE mpi_interface
