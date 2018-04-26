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
MODULE util

   USE mpi_interface, ONLY : cyclics, cyclicc
   IMPLICIT NONE

   INTEGER, SAVE :: fftinix = 0, fftiniy = 0
   CHARACTER (len=6), PARAMETER :: vel_bc = 'frslip'

CONTAINS
   ! ----------------------------------------------------------------------
   ! Subroutine sclrset: Sets upper and lower boundaries to a constant
   ! gradient via extrapolation, or a zero-gradient condition depending
   ! on the flag, typically used to set boundary conditions for scalars
   !
   SUBROUTINE sclrset(type,n1,n2,n3,a,dz)

      USE mpi_interface, ONLY : myid, appl_abort

      INTEGER, INTENT(in)        :: n1,n2,n3
      REAL, INTENT(in), OPTIONAL :: dz(n1)
      REAL, INTENT(inout)        :: a(n1,n2,n3)
      CHARACTER (len=4)          :: type

      INTEGER :: i,j,req(16)
      REAL    :: dzf1,dzf2
      SELECT CASE (type)
         CASE ('grad')
            dzf1 = dz(2)/dz(1)
            dzf2 = dz(n1-2)/dz(n1-1)
         CASE ('cnst')
            dzf1 = 0.
            dzf2 = 0.
         CASE ('mixd')
            dzf1 = 0.
            dzf2 = dz(n1-2)/dz(n1-1)
         CASE DEFAULT
            IF (myid == 0) PRINT *, '  ABORTING:  BCs not supported'
            CALL appl_abort(0)
      END SELECT

      DO j = 1, n3
         DO i = 1, n2
            a(1,i,j)  = a(2,i,j)    - dzf1*(a(3,i,j)    - a(2,i,j))
            a(n1,i,j) = a(n1-1,i,j) + dzf2*(a(n1-1,i,j) - a(n1-2,i,j))
         END DO
      END DO

      CALL cyclics(n1,n2,n3,a,req)
      CALL cyclicc(n1,n2,n3,a,req)

   END SUBROUTINE sclrset
   !
   ! ----------------------------------------------------------------------
   ! VELSET:  Sets boundary conditions for velocity
   !
   SUBROUTINE velset(n1,n2,n3,u,v,w)

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(inout) :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)

      INTEGER :: i, j,req(16)

      CALL cyclics (n1,n2,n3,u,req)
      CALL cyclicc (n1,n2,n3,u,req)
      CALL cyclics (n1,n2,n3,v,req)
      CALL cyclicc (n1,n2,n3,v,req)
      CALL cyclics (n1,n2,n3,w,req)
      CALL cyclicc (n1,n2,n3,w,req)

      DO j = 1 ,n3
         DO i = 1, n2
            w(n1,i,j)    = 0.
            w(n1-1,i,j)  = 0.
            w(1,i,j)     = 0.
         END DO
      END DO
      IF (vel_bc == 'noslip') THEN
         DO j = 1, n3
            DO i = 1, n2
               u(1,i,j)  = -u(2,i,j)
               u(n1,i,j) = -u(n1-1,i,j)
               v(1,i,j)  = -v(2,i,j)
               v(n1,i,j) = -v(n1-1,i,j)
            END DO
         END DO
      ELSE
         DO j = 1, n3
            DO i = 1, n2
               u(1,i,j)  =  u(2,i,j)
               u(n1,i,j) = u(n1-1,i,j)
               v(1,i,j)  = v(2,i,j)
               v(n1,i,j) = v(n1-1,i,j)
            END DO
         END DO
      END IF

   END SUBROUTINE velset
   !
   !---------------------------------------------------------------------
   ! GET_AVG2dh: Get the average of a 2 dimensional (horizontal) input field
   !
   REAL FUNCTION get_avg2dh(n2,n3,a)

      INTEGER, INTENT(in) :: n2,n3
      REAL, INTENT(in)    :: a(n2,n3)

      INTEGER :: i,j
    
      get_avg2dh = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            get_avg2dh = get_avg2dh + a(i,j)
         END DO
      END DO

      get_avg2dh = get_avg2dh/REAL((n3-4)*(n2-4))

   END FUNCTION get_avg2dh
   !
   !---------------------------------------------------------------------
   ! GET_AVG_ts: gets average field value from the whole domain
   ! Implemented by Zubair Maalick
   ! Possibility for conditional sampling added by Juha Tonttila
   !
   ! Weighting by layer thickness implemented for non-uniform vertical resolution
   !
   REAL FUNCTION get_avg_ts(n1,n2,n3,a,dz,cond)

      INTEGER, INTENT (in) :: n1, n2, n3
      REAL, INTENT(in)     :: dz(n1)  ! Reciprocal of layer depth!
      REAL, INTENT (in)    :: a(n1,n2,n3)
      LOGICAL, OPTIONAL, INTENT(in) :: cond(n1,n2,n3)
    
      INTEGER :: i,j,k, npnt
      REAL    :: ztmp,ztot

      npnt = 0
      get_avg_ts = 0.
      IF (present(cond)) THEN
         DO j = 3, n3-2
            DO i = 3, n2-2
               ztot = 0.
               ztmp = 0.
               DO k = 2, n1
                  IF (cond(k,i,j)) THEN
                     ztmp = ztmp + a(k,i,j)*(1./dz(k))
                     ztot = ztot + (1./dz(k))
                  END IF
               END DO
               ! Grid weighted vertical average for columns with at least one available value
               IF (ztot /= 0.0 ) THEN
                  get_avg_ts = get_avg_ts + ztmp/ztot
                  npnt = npnt + 1
               END IF
            END DO
         END DO
       
      ELSE
         DO j = 3, n3-2
            DO i = 3, n2-2
               ztot = 0.
               ztmp = 0.
               DO k = 2, n1
                  ztmp = ztmp + a(k,i,j)*(1./dz(k))
                  ztot = ztot + (1./dz(k))
               END DO
               ! Grid wighted vertical average
               ztmp = ztmp/ztot
               npnt = npnt + 1
               get_avg_ts = get_avg_ts + ztmp
            END DO
         END DO

      END IF

      IF (npnt > 0) get_avg_ts = get_avg_ts/REAL(npnt)

   END FUNCTION get_avg_ts
   !
   !---------------------------------------------------------------------
   ! GET_AVG3: gets average across outer two dimensions at each
   ! point along inner dimension - calculated over all PUs
   !
   SUBROUTINE get_avg3(n1,n2,n3,a,avg,normalize,cond)

      USE mpi_interface, ONLY : nypg,nxpg,double_array_par_sum

      INTEGER,INTENT(in) :: n1,n2,n3
      REAL,INTENT(in)    :: a(n1,n2,n3)
      REAL,INTENT(out)   :: avg(n1)
      LOGICAL, OPTIONAL  :: normalize
      LOGICAL, OPTIONAL  :: cond(n1,n2,n3)

      INTEGER      :: k,i,j
      REAL(kind=8) :: lavg(n1),gavg(n1),counts(n1), x
      LOGICAL :: norm

      norm = .TRUE. ! Default
      IF (present(normalize)) norm = normalize

      IF (present(cond)) THEN
         gavg(:) = 0.
         counts(:) = 0.
         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1
                  IF (cond(k,i,j)) THEN
                     gavg(k) = gavg(k)+a(k,i,j)
                     counts(k) = counts(k)+1.
                  END IF
               END DO
            END DO
         END DO
         lavg = gavg
         CALL double_array_par_sum(lavg,gavg,n1)
         IF (norm) THEN
            lavg = counts
            CALL double_array_par_sum(lavg,counts,n1)
            avg(:) = 0.
            DO k = 1, n1
               IF (counts(k) > 0.) avg(k) = REAL(gavg(k)/counts(k))
            END DO
         ELSE
            avg(:) = REAL(gavg(:))
         END IF
      ELSE
         x = 1.
         IF (norm) x = 1./(REAL(nypg-4)*REAL(nxpg-4))

         gavg(:) = 0.
         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1
                  gavg(k) = gavg(k)+a(k,i,j)
               END DO
            END DO
         END DO
         lavg = gavg
         CALL double_array_par_sum(lavg,gavg,n1)
         avg(:) = REAL(gavg(:) * x)
      END IF

   END SUBROUTINE get_avg3
   !
   !---------------------------------------------------------------------
   ! Function get_cor: gets mean correlation between two fields at a
   ! given level
   !
   REAL FUNCTION get_cor(n1,n2,n3,k,a,b)

      INTEGER, INTENT (in) :: n1,n2,n3,k
      REAL, INTENT (inout) :: a(n1,n2,n3),b(n1,n2,n3)

      INTEGER :: i,j

      get_cor = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            get_cor = get_cor+a(k,i,j)*b(k,i,j)
         END DO
      END DO
      get_cor = get_cor/REAL((n3-4)*(n2-4))

   END FUNCTION get_cor
   !
   !---------------------------------------------------------------------
   ! Function get_cor3: gets mean correlation accross outer two dimensions
   ! at each point along inner dimension
   !
   SUBROUTINE get_cor3(n1,n2,n3,a,b,avg)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: a(n1,n2,n3),b(n1,n2,n3)
      REAL, INTENT (out)   :: avg(n1)

      INTEGER :: k,i,j

      avg(:) = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               avg(k) = avg(k)+a(k,i,j)*b(k,i,j)
            END DO
         END DO
      END DO
      avg(:) = avg(:)/REAL((n3-4)*(n2-4))

   END SUBROUTINE get_cor3
   !
   !---------------------------------------------------------------------
   ! Function get_var3: gets variance for a field whose mean is known
   !
   SUBROUTINE get_var3(n1,n2,n3,a,b,avg)

      INTEGER :: n1,n2,n3,k,i,j
      REAL    :: a(n1,n2,n3),b(n1),avg(n1)

      avg(:) = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               avg(k) = avg(k)+(a(k,i,j)-b(k))**2
            END DO
         END DO
      END DO
      avg(:) = avg(:)/REAL((n3-4)*(n2-4))

   END SUBROUTINE get_var3
   !
   !---------------------------------------------------------------------
   ! Function get_3rd3: gets the third moment for a field whose mean is known
   !
   SUBROUTINE get_3rd3(n1,n2,n3,a,b,avg)

      INTEGER :: n1,n2,n3,k,i,j
      REAL    :: a(n1,n2,n3),b(n1),avg(n1)

      avg(:) = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               avg(k) = avg(k)+(a(k,i,j)-b(k))**3
            END DO
         END DO
      END DO
      avg(:) = avg(:)/REAL((n3-4)*(n2-4))

   END SUBROUTINE get_3rd3
   !
   ! ----------------------------------------------------------------------
   ! Subroutine tridiff: standard tri-diagonal solver for nh columns
   ! uses the LU decomposition of the equation
   !
   !     cin1(k)*x(k-1)+ci(k)*x(k)+cip1(k)*x(k+1) = b(k)
   !
   ! such that
   !
   !         |  1   0   0   0   0 ... 0 | |u11 a12  0   0   0   0  ...  0  |
   !     L = | l21  1   0   0   0 ... 0 | | 0  u22 a23  0   0   0  ...  0  |
   !         |  0  l32  1   0   0 ... 0 | | 0   0  u33 a34  0   0  ...  0  |
   !         |                          | |                                |
   !         |  0   0   0   0   0 ... 1 | | 0   0   0   0   0   0  ... unn |
   !
   ! where aik =cip1(k) (i=k-1), =cin1(k) i=k+1
   !    u(1,1)   = a(1,1)
   !    l(i,i-1) = a(i,i-1)/u(i-1,i-1)
   !    u(i,i)   = a(i,i1)-l(i,i-1)*a(i-1,i)
   !
   ! and solves for x(k) = y(k)/u(k,k) - a(k,k+1)*x(k+1)/u(k,k) where
   ! y(k) = b(k) -l(k)y(k-1)
   !
   SUBROUTINE tridiff(nh,n1,nhdo,cin1,ci,cip1,rhs,cj,cjp1)

      INTEGER, INTENT(in) :: nh,n1,nhdo
      REAL, INTENT(in)    :: cin1(nh,n1),ci(nh,n1),cip1(nh,n1)
      REAL, INTENT(inout) :: rhs(nh,n1),cj(nh,n1),cjp1(nh,n1)

      INTEGER k,i
      REAL eps

      eps = sqrt(tiny(1.))

      DO i = 1, nhdo
         cjp1(i,2) = cip1(i,2)/ci(i,2)
         rhs(i,2) = rhs(i,2)/ci(i,2)
         DO k = 3, n1
            cj(i,k) = ci(i,k)-cin1(i,k)*cjp1(i,k-1)+eps
            cjp1(i,k) = cip1(i,k)/cj(i,k)
            rhs(i,k) = (rhs(i,k)-cin1(i,k)*rhs(i,k-1))/cj(i,k)
         END DO
         !
         ! here rhs = y(k)/u(k,k), cjp1=a(k,k+1)/u(k,k)
         !
         cj(i,n1) = rhs(i,n1)
         DO k = n1-1, 2, -1
            cj(i,k) = rhs(i,k)-cjp1(i,k)*cj(i,k+1)
         END DO
      END DO

   END SUBROUTINE tridiff
   !
   ! --------------------------------------------------------------------
   ! Subroutine ae1mm: subtracts mean value from given field (a=a-a_bar)
   !
   SUBROUTINE ae1mm(n1,n2,n3,a,abar)

      USE mpi_interface, ONLY : nypg,nxpg,double_array_par_sum

      INTEGER n1,n2,n3
      REAL, INTENT (inout), DIMENSION (n1,n2,n3) :: a(n1,n2,n3)
      REAL, INTENT (out), DIMENSION (n1)         :: abar(n1)

      INTEGER :: i,j,k

      REAL(kind=8) :: lavg(n1),gavg(n1)

      ! TR: this used to be for the whole domain ...
      !call get_avg3(n1,n2,n3,a,abar)

      gavg(:) = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               gavg(k) = gavg(k)+a(k,i,j)
            END DO
         END DO
      END DO
      lavg = gavg
      CALL double_array_par_sum(lavg,gavg,n1)
      abar(:) = REAL( gavg(:)/REAL((nypg-4)*(nxpg-4)) )

      DO j = 1, n3
         DO i = 1, n2
            DO k = 1, n1
               a(k,i,j) = a(k,i,j)-abar(k)
            END DO
         END DO
      END DO

   END SUBROUTINE ae1mm
   !
   !---------------------------------------------------------------------
   ! CRAYFFTUSE:  Uses the cray routines to do a 2D transform
   !
   SUBROUTINE get_fft_twodim(nx,ny,nz,a,wsavex,wsavey,isgn)

      USE mpi_interface, ONLY : xshuffle,yshuffle,nxg,nyg,nynzp,nxnzp

      INTEGER, INTENT(in)    :: nx,ny,nz,isgn
      COMPLEX, INTENT(inout) :: a(nx,ny,nz)
      REAL, INTENT(inout)    :: wsavex(4*nxg+100),wsavey(4*nyg+100)

      INTEGER :: k, j, i
      COMPLEX :: atmp((nx+1)*(ny+1)*(nz+1)),btmp(ny,nx,nz)

      CALL xshuffle(a,atmp,nx,ny,nz,1)
      CALL fft1dc(nxg,nynzp,atmp,wsavex,isgn,fftinix)
      CALL xshuffle(a,atmp,nx,ny,nz,-1)

      DO k = 1, nz
         DO j = 1, ny
            DO i = 1, nx
               btmp(j,i,k) = A(i,j,k)
            END DO
         END DO
      END DO

      CALL yshuffle(btmp,atmp,nx,ny,nz,1)
      CALL fft1dc(nyg,nxnzp,atmp,wsavey,isgn,fftiniy)
      CALL yshuffle(btmp,atmp,nx,ny,nz,-1)

      DO k = 1, nz
         DO j = 1, ny
            DO i = 1, nx
               A(i,j,k) = btmp(j,i,k)
            END DO
         END DO
      END DO

   END SUBROUTINE get_fft_twodim
   !
   !---------------------------------------------------------------
   ! MASKACTIV: Create a LOGICAL mask for grid points where cloud
   !            activation will be calculated.
   !
   ! Juha Tonttila, FMI, 2014
   !
   SUBROUTINE maskactiv(act_mask,nx,ny,nz,mode,rh,rc,w)
      IMPLICIT NONE

      INTEGER, INTENT(in) :: nx,ny,nz
      REAL, INTENT(in)    :: rh(nz,nx,ny)
      REAL, OPTIONAL, INTENT(in) :: rc(nz,nx,ny),w(nz,nx,ny)

      INTEGER, INTENT(in) :: mode ! 1 = Initialization; 2 = Normal timestepping

      LOGICAL, INTENT(out) :: act_mask(nz,nx,ny)

      LOGICAL :: actmask_newcloud(nz,nx,ny)
      LOGICAL :: actmask_oldcloud(nz,nx,ny)
      LOGICAL :: cldmask(nz,nx,ny)
      LOGICAL :: cldm1(nz,nx,ny),   & !inverse cloud mask offset downwards by one grid level
                 cldp1(nz,nx,ny),   & !inverse cloud mask offset upwards by one grid level
                 cldpm(nz,nx,ny)

      LOGICAL :: notused(nx,ny)

      INTEGER :: k

      actmask_newcloud = .FALSE.
      actmask_oldcloud = .FALSE.
      cldm1 = .TRUE.
      cldp1 = .TRUE.
      cldpm = .TRUE.
      act_mask = .FALSE.

      SELECT CASE(mode)
         CASE(1)
            ! Calculate cloud activation in all cloudy grid points for initialization
            ! i.e. points with RH >= threshold
            !zrh = rv/rs
            act_mask = ( rh >= 1.000 )

         CASE(2)

            ! Normal opperation
      
            IF ( .NOT. (present(rc) )) STOP 'maskactiv: invalid arguments'

            ! Mask grid points just below cloud or where new cloud is expected to form
            !
       
            cldmask(:,:,:) = ( rc(:,:,:) >= 1.e-5 )

      
            ! Check for the lowest point where RH > 1.
            cldpm(:,:,:) = .FALSE.
            cldpm(1:nz-1,:,:) = ( rh(2:nz,:,:) > 1.000 )

            cldp1(:,:,:) = .TRUE.
            cldp1(2:nz,:,:) = ( .NOT. cldpm(1:nz-1,:,:) )
            actmask_newcloud(:,:,:) = ( cldpm(:,:,:) .AND. cldp1(:,:,:) ) .AND. ( w(:,:,:) > 0. ) &
                                       .AND. (rc(:,:,:) < 5.e-5)

            ! Base of existing cloud
            cldpm(:,:,:) = .FALSE.
            cldpm(1:nz-1,:,:) = cldmask(2:nz,:,:)
       
            cldp1(:,:,:) = .TRUE.
            cldp1(2:nz,:,:) = ( .NOT. cldpm(1:nz-1,:,:) )

            ! Take the lowest level of the two cases
       
            notused = .TRUE.
            DO k = 2, nz-1
               ! New cloud + no old cloud
               act_mask(k,:,:) = MERGE( (actmask_oldcloud(k,:,:) .OR. actmask_newcloud(k,:,:))  &
                                 .AND. cldp1(k,:,:), .FALSE., notused(:,:))
               notused(:,:) = notused(:,:) .AND. .NOT. act_mask(k,:,:)
               notused(:,:) = (rh(k,:,:) < 0.99 .OR. notused(:,:))
            END DO

      END SELECT
     
   END SUBROUTINE maskactiv

   FUNCTION closest(array,val)
     ! Find the index of "array" with value closest to "val"
     IMPLICIT NONE
     INTEGER :: closest
     REAL, INTENT(in) :: array(:)
     REAL, INTENT(in) :: val

     INTEGER :: NN, N
     LOGICAL, ALLOCATABLE :: comp(:)

     NN = SIZE(array)
     N = smaller(array,val)
     
     IF ( N < NN .AND.                                 &
          ( ABS(array(N)-val) > ABS(array(N+1)-val) ) )  &
        N = N + 1
  
     closest = MAX(MIN(N,NN),1)

   END FUNCTION closest

   FUNCTION smaller(array,val)
     ! Find out how many elements of "array" have value smaller than "val"
     IMPLICIT NONE
     INTEGER :: smaller
     REAL, INTENT(in) :: array(:)
     REAL, INTENT(in) :: val

     INTEGER :: NN,N
     LOGICAL, ALLOCATABLE :: comp(:)

     NN = SIZE(array)

     ALLOCATE(comp(NN))
     
     comp = .FALSE.
     comp(:) = (array(:) < val)

     N = COUNT(comp)
     smaller = MAX(MIN(N,NN),1)

   END FUNCTION smaller

   ! ------------------------------
   
   !
   ! --------------------------------------------------------------------------- 
   ! For Level >= 4: Returns the index for mass mixing ratio in the prognostic
   ! tracer arrays for a given SALSA size bin and a given aerosol species (or 
   ! water)
   ! Input arguments: nbtot = total number of bins
   !                  nb    = number of the bin for which the index is fetched
   !                  nm    = index of the aerosol species (use classSpecies for this)
   !
   INTEGER FUNCTION getMassIndex(nbtot,nb,nm)
     IMPLICIT NONE
     INTEGER, INTENT(in) :: nbtot, nb, nm

     getMassIndex = (nm-1)*nbtot+nb

   END FUNCTION getMassIndex


   ! Function for calculating Pearson's correlation coefficient for two vectors
   REAL FUNCTION calc_correlation(x,y,n)
     INTEGER, INTENT(in) :: n
     REAL, INTENT(in) :: x(n), y(n)
     REAL :: sx, sy, sx2, sy2, sxy
     INTEGER :: i
     REAL, PARAMETER :: eps = EPSILON(1.0)
     sx=0.; sy=0.; sx2=0.; sy2=0.; sxy=0.
     DO i=1,n
        sx=sx+x(i)
        sy=sy+y(i)
        sx2=sx2+x(i)**2
        sy2=sy2+y(i)**2
        sxy=x(i)*y(i)
     END DO
     IF (sx2*n-sx**2<eps .OR. sy2*n-sy**2<eps) THEN
        calc_correlation = 0.
     ELSE
        calc_correlation = ( sxy*n-sx*sy )/( SQRT(sx2*n-sx**2)*SQRT(sy2*n-sy**2) )
     END IF
   END FUNCTION calc_correlation
   

END MODULE util
