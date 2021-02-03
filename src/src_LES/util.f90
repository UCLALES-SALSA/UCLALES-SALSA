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
module util

  use mpi_interface, only : cyclics, cyclicc
  implicit none

  integer, save :: fftinix=0, fftiniy=0
  logical, parameter  :: noslip = .false.

contains
  ! ----------------------------------------------------------------------
  ! Subroutine sclrset: Sets upper and lower boundaries to a constant
  ! gradient via extrapolation, or a zero-gradient condition depending
  ! on the flag, typically used to set boundary conditions for scalars
  !
  subroutine sclrset(type,n1,n2,n3,a,dz)

    use mpi_interface, only : myid, appl_abort

    integer, intent(in)           :: n1,n2,n3
    real, intent(in), optional    :: dz(n1)
    real, intent(inout)           :: a(n1,n2,n3)
    character (len=4)             :: type

    integer :: i,j,req(16)
    real :: dzf1,dzf2
    select case (type)
    case ('grad')
       dzf1  = dz(2)/dz(1)
       dzf2  = dz(n1-2)/dz(n1-1)
    case ('cnst')
       dzf1  = 0.
       dzf2  = 0.
    case ('mixd')
       dzf1  = 0.
       dzf2  = dz(n1-2)/dz(n1-1)
    case default
       if (myid == 0) print *, '  ABORTING:  BCs not supported'
       call appl_abort(0)
    end select

    do j=1,n3
       do i=1,n2
          a(1,i,j)  = a(2,i,j)    - dzf1*(a(3,i,j)    - a(2,i,j))
          a(n1,i,j) = a(n1-1,i,j) + dzf2*(a(n1-1,i,j) - a(n1-2,i,j))
       end do
    end do

    call cyclics(n1,n2,n3,a,req)
    call cyclicc(n1,n2,n3,a,req)

  end subroutine sclrset
  !
  ! ----------------------------------------------------------------------
  ! VELSET:  Sets boundary conditions for velocity
  !
  subroutine velset(n1,n2,n3,u,v,w)

    integer, intent(in) :: n1,n2,n3
    real, intent(inout) :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)

    integer :: i, j,req(16)

    call cyclics (n1,n2,n3,u,req)
    call cyclicc (n1,n2,n3,u,req)
    call cyclics (n1,n2,n3,v,req)
    call cyclicc (n1,n2,n3,v,req)
    call cyclics (n1,n2,n3,w,req)
    call cyclicc (n1,n2,n3,w,req)

    do j=1,n3
       do i=1,n2
          w(n1,i,j)    = 0.
          w(n1-1,i,j)  = 0.
          w(1,i,j)     = 0.
       enddo
    enddo
    if (noslip) then
       do j=1,n3
          do i=1,n2
             u(1,i,j)  = -u(2,i,j)
             u(n1,i,j) = -u(n1-1,i,j)
             v(1,i,j)  = -v(2,i,j)
             v(n1,i,j) = -v(n1-1,i,j)
          end do
       end do
    else
       do j=1,n3
          do i=1,n2
             u(1,i,j)  =  u(2,i,j)
             u(n1,i,j) = u(n1-1,i,j)
             v(1,i,j)  = v(2,i,j)
             v(n1,i,j) = v(n1-1,i,j)
          end do
       end do
    end if

  end subroutine velset
  !
  !---------------------------------------------------------------------
  ! GET_AVG2dh: Get the average of a 2 dimensional (horizontal) input field - calculated over all PUs
  !
  REAL FUNCTION get_avg2dh(n2,n3,a)
    use mpi_interface, only : nypg,nxpg,double_scalar_par_sum

    INTEGER, INTENT(in) :: n2,n3
    REAL, INTENT(in)    :: a(n2,n3)

    REAL(kind=8) :: lavg,gavg
    INTEGER :: i,j
    
    get_avg2dh = 0.
    DO j = 3,n3-2
       DO i = 3,n2-2
          get_avg2dh = get_avg2dh + a(i,j)
       END DO
    END DO

    lavg = get_avg2dh
    call double_scalar_par_sum(lavg,gavg)

    get_avg2dh = real(gavg)/real((nypg-4)*(nxpg-4))

  END FUNCTION get_avg2dh
  !
  !---------------------------------------------------------------------
  ! GET_AVG_ts: gets average field value from the whole domain 
  ! Implemented by Zubair Maalick
  ! Possibility for conditional sampling added by Juha Tonttila
  ! Possibility for vertical integral with air density as a weight
  ! Weighting by layer thickness implemented for non-uniform vertical resolution  - calculated over all PUs
  !
  real function get_avg_ts(n1,n2,n3,a,dz,cond,dens)
    use mpi_interface, only : double_scalar_par_sum

    integer, intent (in):: n1, n2, n3
    REAL, INTENT(in)    :: dz(n1)  ! Reciprocal of layer depth!
    real, intent (in)   :: a(n1,n2,n3)
    LOGICAL, OPTIONAL, INTENT(in) :: cond(n1,n2,n3)
    REAL, OPTIONAL, INTENT(in) :: dens(n1,n2,n3)
    
    integer :: i,j,k, npnt
    REAL :: ztmp,ztot
    REAL(kind=8) :: lavg,gavg,nn

    npnt = 0
    get_avg_ts=0.
    IF (PRESENT(cond) .AND. PRESENT(dens)) THEN
       ! Conditional vertical integral with density weights
       DO j=3,n3-2
          DO i=3,n2-2
             ztmp = 0.
             DO k = 2,n1
                IF (cond(k,i,j)) ztmp = ztmp + a(k,i,j)*(dens(k,i,j)/dz(k))
             END DO
             ! Vertical integral, so no need to normalize
             get_avg_ts = get_avg_ts + ztmp
             npnt = npnt + 1
          END DO
       END DO
    ELSEIF (PRESENT(cond)) THEN
       ! Conditional average
       DO j=3,n3-2
          DO i=3,n2-2
             ztot = 0.
             ztmp = 0.
             DO k = 2,n1
                IF (cond(k,i,j)) THEN
                   ztmp = ztmp + a(k,i,j)*(1./dz(k))
                   ztot = ztot + (1./dz(k))
                END IF
             END DO
             ! Grid weighted vertical average for columns with at least one available value
             if (ztot /=0.0 ) THEN
                get_avg_ts = get_avg_ts + ztmp/ztot
                npnt = npnt + 1
             END IF
          END DO
       END DO
    ELSEIF (PRESENT(dens)) THEN
       ! Vertical integral with density weights
       DO j=3,n3-2
          DO i=3,n2-2
             ztmp = 0.
             DO k = 2,n1
                ztmp = ztmp + a(k,i,j)*(dens(k,i,j)/dz(k))
             END DO
             ! Vertical integral, so no need to normalize
             get_avg_ts = get_avg_ts + ztmp
             npnt = npnt + 1
          END DO
       END DO
    ELSE
       ! Average
       DO j=3,n3-2
          DO i=3,n2-2
             ztot = 0.
             ztmp = 0.
             DO k = 2,n1
                ztmp = ztmp + a(k,i,j)*(1./dz(k))
                ztot = ztot + (1./dz(k))
             END DO
             ! Grid weighted vertical average
             get_avg_ts = get_avg_ts + ztmp/ztot
             npnt = npnt + 1
          END DO
       END DO
    END IF

    lavg = REAL(npnt)
    call double_scalar_par_sum(lavg,gavg)
    IF (gavg>0.) THEN
        nn = gavg ! npnt
        lavg = get_avg_ts
        call double_scalar_par_sum(lavg,gavg)
        get_avg_ts = real(gavg/nn)
    ELSE
        get_avg_ts = -999.
    ENDIF

  end function get_avg_ts
  !
  !---------------------------------------------------------------------
  ! GET_AVG3: gets average across outer two dimensions at each
  ! point along inner dimension - calculated over all PUs
  !
  subroutine get_avg3(n1,n2,n3,a,avg,normalize,cond)

    use mpi_interface, only : nypg,nxpg,double_array_par_sum

    integer,intent(in) :: n1,n2,n3
    real,intent(in) :: a(n1,n2,n3)
    real,intent(out) :: avg(n1)
    LOGICAL, OPTIONAL :: normalize
    LOGICAL, OPTIONAL :: cond(n1,n2,n3)

    integer      :: k,i,j
    real(kind=8) :: lavg(n1),gavg(n1),counts(n1), x
    LOGICAL :: norm

    norm = .TRUE. ! Default
    IF (present(normalize)) norm=normalize

    IF (present(cond)) THEN
        gavg(:) = 0.
        counts(:) = 0.
        do j=3,n3-2
           do i=3,n2-2
              do k=1,n1
                 IF (cond(k,i,j)) THEN
                    gavg(k)=gavg(k)+a(k,i,j)
                    counts(k)=counts(k)+1.
                 ENDIF
              end do
           end do
        end do
        lavg = gavg
        call double_array_par_sum(lavg,gavg,n1)
        IF (norm) THEN
            lavg = counts
            call double_array_par_sum(lavg,counts,n1)
            avg(:)=0.
            do k=1,n1
                IF (counts(k)>0.) avg(k) = real(gavg(k)/counts(k))
            END DO
        ELSE
            avg(:)=REAL(gavg(:))
        ENDIF
    ELSE
        x=1.
        IF (norm) x = 1./(real(nypg-4)*real(nxpg-4))

        gavg(:) = 0.
        do j=3,n3-2
           do i=3,n2-2
              do k=1,n1
                 gavg(k)=gavg(k)+a(k,i,j)
              end do
           end do
        end do
        lavg = gavg
        call double_array_par_sum(lavg,gavg,n1)
        avg(:) = real(gavg(:) * x)
    ENDIF

  end subroutine get_avg3
  !
  ! -------------------------------------------------------------------------
  !
  ! find the mean level at which sx=threshold - calculated over all PUs
  real function get_zi_val(n1, n2, n3, sx, z, threshold)
    use mpi_interface, only : double_scalar_par_sum

    integer, intent (in) :: n1, n2, n3
    real, intent (in)    :: z(n1), sx(n1,n2,n3), threshold

    integer :: i, j, k, n
    real    :: zibar
    REAL(kind=8) :: lavg,gavg,nn

    n = 0
    zibar = 0.
    do j=3,n3-2
        do i=3,n2-2
            k = 2
            do while (k < n1-2 .and. sx(k,i,j) > threshold)
                k = k+1
            end do
            IF (k < n1-2) THEN
                ! Level found
                n = n+1
                ! Interpolation
                zibar = zibar + z(k-1) +  &
                  (threshold - sx(k-1,i,j))/(z(k)-z(k-1))     /  &
                  (sx(k,i,j) - sx(k-1,i,j) + epsilon(1.))
            ENDIF
        end do
    end do

    ! Global mean
    lavg = REAL(n)
    call double_scalar_par_sum(lavg,gavg)
    IF (gavg>0.) THEN
        nn =gavg
        lavg = zibar
        call double_scalar_par_sum(lavg,gavg)
        get_zi_val = real(gavg/nn)
    ELSE
        get_zi_val = -999.
    ENDIF

  END function get_zi_val
  !
  ! -------------------------------------------------------------------------
  !
  ! Find the mean height where sx has its maximum gradient, positive or negative - calculated over all PUs
  real function get_zi_dmax(n1, n2, n3, sx, z)
    use mpi_interface, only : nypg,nxpg,double_scalar_par_sum

    integer, intent (in) :: n1, n2, n3
    real, intent (in)    :: z(n1), sx(n1,n2,n3)

    integer :: i, j, k
    real    :: sval, dmy, scr
    REAL(kind=8) :: lavg,gavg

    get_zi_dmax = 0.
    do j=3,n3-2
        do i=3,n2-2
            sval = 0. ! for k=1
            scr = z(1)
            do k=2,n1-5
                dmy = abs( (sx(k+1,i,j)-sx(k,i,j))/(z(k+1)-z(k)) )
                if (dmy > sval) then
                   sval = dmy
                   scr = z(k)
                end if
            end do
            get_zi_dmax = get_zi_dmax + scr
        end do
    end do

    lavg = get_zi_dmax
    call double_scalar_par_sum(lavg,gavg)
    get_zi_dmax = real(gavg)/real((nypg-4)*(nxpg-4))

  end function get_zi_dmax
  !
  ! -------------------------------------------------------------------------
  !
  ! Find the maximum of sx - calculated over all PUs
  real function get_max_val(n1, n2, n3, sx)
    use mpi_interface, only : double_scalar_par_max

    integer, intent (in) :: n1, n2, n3
    real, intent (in)    :: sx(n1,n2,n3)

    REAL(kind=8) :: lavg,gavg

    lavg = maxval(sx(2:n1,3:n2-2,3:n3-2))
    call double_scalar_par_max(lavg,gavg)
    get_max_val = REAL(gavg)

  end function get_max_val
  !
  ! -------------------------------------------------------------------------
  !
  ! Statistics of a scalar calculated over all PUs
  real function get_pustat_scalar(op, sx, wx)
    use mpi_interface, only : pecount, double_scalar_par_max, double_scalar_par_sum

    CHARACTER(LEN=3) :: op ! Operation
    real, intent (in)    :: sx ! Data
    real, OPTIONAL, intent (in) :: wx ! Weight (optional, for average only)

    REAL(kind=8) :: lavg,gavg,sw

    select case(op)
    CASE('avg')
        ! Average
        IF (PRESENT(wx)) THEN
            ! Weighted average: avg = sum(x(i)*w(i),i=1,n)/sum(w(i),i=1,n)
            lavg = wx
            call double_scalar_par_sum(lavg,gavg)
            IF (ABS(gavg)>1e-30) THEN
                sw = gavg
                lavg = sx
                call double_scalar_par_sum(lavg,gavg)
                get_pustat_scalar = REAL(gavg/sw)
            ELSE
                get_pustat_scalar = -999.
            ENDIF
        ELSE
            ! Average: avg = sum(x(i),i=1,n)/n
            lavg = sx
            call double_scalar_par_sum(lavg,gavg)
            get_pustat_scalar = REAL(gavg)/REAL(pecount)
        ENDIF
    CASE('sum')
        ! Sum
        lavg = sx
        call double_scalar_par_sum(lavg,gavg)
        get_pustat_scalar = REAL(gavg)
    CASE('max')
        ! Maximum
        lavg = sx
        call double_scalar_par_max(lavg,gavg)
        get_pustat_scalar = REAL(gavg)
    CASE('min')
        ! Minimum (<1e30)
        lavg = -sx
        call double_scalar_par_max(lavg,gavg)
        get_pustat_scalar = -REAL(gavg)
    case default
        WRITE(*,*) op
        STOP 'Bad option for get_pustat_scalar!'
    END SELECT

  end function get_pustat_scalar
  ! -------------------------------------------------------------------------
  !
  ! Statistics of a vector calculated over all PUs
  SUBROUTINE get_pustat_vector(op, n, sx, wx)
    use mpi_interface, only : pecount, double_array_par_sum

    CHARACTER(LEN=3) :: op        ! Operation
    integer, intent(in) :: n      ! Dimension
    real, intent (inout) :: sx(n) ! Data
    real, OPTIONAL, intent (in) :: wx(n) ! Weight (optional, for average only)

    INTEGER :: i
    REAL(kind=8) :: lavg(n),gavg(n),sw(n)

    select case(op)
    CASE('avg')
        ! Average
        IF (PRESENT(wx)) THEN
            ! Weighted average: avg = sum(x(i)*w(i),i=1,n)/sum(w(i),i=1,n)
            lavg = wx
            call double_array_par_sum(lavg,gavg,n)
            sw = gavg
            lavg = sx
            call double_array_par_sum(lavg,gavg,n)
            DO i=1,n
                IF (ABS(gavg(i))>1e-30) THEN
                    sx(i) = REAL(gavg(i)/sw(i))
                ELSE
                    sx(i) = -999.
                ENDIF
            ENDDO
        ELSE
            ! Average: avg = sum(x(i),i=1,n)/n
            lavg = sx
            call double_array_par_sum(lavg,gavg,n)
            sx(:) = REAL(gavg(:))/REAL(pecount)
        ENDIF
    CASE('sum')
        ! Sum
        lavg = sx
        call double_array_par_sum(lavg,gavg,n)
        sx(:) = REAL(gavg(:))
    case default
        WRITE(*,*) op
        STOP 'Bad option for get_pustat_vector!'
    END SELECT

  end SUBROUTINE get_pustat_vector
  !
  !---------------------------------------------------------------------
  ! Calculate histograms from concentration (num) and radius (rad) data - calculated over all PUs
  !
  SUBROUTINE HistDistr(n1,n2,n3,n4,rad,num,rbins,nout,hist)
    use mpi_interface, only : nypg,nxpg,double_array_par_sum
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1, n2, n3, n4, nout ! Dimensions
    REAL, INTENT(IN) :: rad(n1,n2,n3,n4), num(n1,n2,n3,n4) ! Size and number data
    REAL, INTENT(IN) :: rbins(nout+1) ! Size bin limits
    REAL, INTENT(OUT) :: hist(n1,nout) ! Output histogram
    REAL :: rmin, rmax
    INTEGER :: i, j, k, l, ii
    real(kind=8) :: lavg(n1*nout),gavg(n1*nout)
    !
    ! Include values that are within the bin limits
    rmin=MAX(rbins(1),1e-10) ! r=0 when there are no particles
    rmax=rbins(nout+1)
    !
    lavg=0.
    DO l=1,n4
      do j=3,n3-2
        do i=3,n2-2
          do k=2,n1
            IF (rmin<=rad(k,i,j,l) .AND. rad(k,i,j,l)<=rmax) THEN
                ii=nout
                DO WHILE (rad(k,i,j,l)<rbins(ii) .AND. ii>1)
                    ii=ii-1
                ENDDO
                ii=k+(ii-1)*n1
                lavg(ii)=lavg(ii)+num(k,i,j,l)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    ! Sum over PUs
    call double_array_par_sum(lavg,gavg,n1*nout)
    hist = RESHAPE( gavg, (/n1,nout/) )
    !
    ! Normalize by the number of columns
    hist(:,:)=hist(:,:)/(real(nypg-4)*real(nxpg-4))
    !
  END SUBROUTINE HistDistr
  !
  !---------------------------------------------------------------------
  ! function get_cor: gets mean correlation between two fields at a 
  ! given level - calculated over all PUs
  !
  real function get_cor(n1,n2,n3,k,a,b)
    use mpi_interface, only : nypg,nxpg,double_scalar_par_sum

    integer, intent (in) :: n1,n2,n3,k
    real, intent (in) :: a(n1,n2,n3),b(n1,n2,n3)

    REAL(kind=8) :: lavg,gavg
    integer :: i,j

    get_cor=0.
    do j=3,n3-2
       do i=3,n2-2
          get_cor=get_cor+a(k,i,j)*b(k,i,j)
       end do
    end do

    lavg = get_cor
    call double_scalar_par_sum(lavg,gavg)
    get_cor = real(gavg)/real((nypg-4)*(nxpg-4))

  end function get_cor
  !
  !---------------------------------------------------------------------
  ! function get_cor3: gets mean correlation accross outer two dimensions
  ! at each point along inner dimension - calculated over all PUs
  !
  subroutine get_cor3(n1,n2,n3,a,b,avg)
    use mpi_interface, only : nypg,nxpg,double_array_par_sum

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: a(n1,n2,n3),b(n1,n2,n3)
    real, intent (out)   :: avg(n1)

    real(kind=8) :: lavg(n1),gavg(n1)
    integer :: k,i,j

    avg(:) = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             avg(k)=avg(k)+a(k,i,j)*b(k,i,j)
          end do
       enddo
    enddo

    lavg(:) = avg(:)
    call double_array_par_sum(lavg,gavg,n1)
    avg(:) = real(gavg(:))/real((nypg-4)*(nxpg-4))

  end subroutine get_cor3
  !
  !---------------------------------------------------------------------
  ! function get_var3: gets variance for a field whose mean is known - calculated over all PUs
  !
  subroutine get_var3(n1,n2,n3,a,b,avg)
    use mpi_interface, only : nypg,nxpg,double_array_par_sum

    integer n1,n2,n3,k,i,j
    real a(n1,n2,n3),b(n1),avg(n1)
    real(kind=8) :: lavg(n1),gavg(n1)

    avg(:) = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             avg(k)=avg(k)+(a(k,i,j)-b(k))**2
          end do
       enddo
    enddo

    lavg = avg
    call double_array_par_sum(lavg,gavg,n1)
    avg(:) = real(gavg(:))/real((nypg-4)*(nxpg-4))

  end subroutine get_var3
  !
  !---------------------------------------------------------------------
  ! function get_3rd3: gets the third moment for a field whose mean is known - calculated over all PUs
  !
  subroutine get_3rd3(n1,n2,n3,a,b,avg)
    use mpi_interface, only : nypg,nxpg,double_array_par_sum

    integer n1,n2,n3,k,i,j
    real a(n1,n2,n3),b(n1),avg(n1)
    real(kind=8) :: lavg(n1),gavg(n1)

    avg(:) = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             avg(k)=avg(k)+(a(k,i,j)-b(k))**3
          end do
       enddo
    enddo

    lavg = avg
    call double_array_par_sum(lavg,gavg,n1)
    avg(:) = real(gavg(:))/real((nypg-4)*(nxpg-4))
  end subroutine get_3rd3
  !
  ! ----------------------------------------------------------------------
  ! subroutine tridiff: standard tri-diagonal solver for nh columns
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
  subroutine tridiff(nh,n1,nhdo,cin1,ci,cip1,rhs,cj,cjp1)

    integer, intent(in) :: nh,n1,nhdo
    real, intent(in)    :: cin1(nh,n1),ci(nh,n1),cip1(nh,n1)
    real, intent(inout) :: rhs(nh,n1),cj(nh,n1),cjp1(nh,n1)

    integer k,i
    real eps

    eps=sqrt(tiny(1.))

    do i=1,nhdo
       cjp1(i,2)=cip1(i,2)/ci(i,2)
       rhs(i,2)=rhs(i,2)/ci(i,2)
       do k=3,n1
          cj(i,k)=ci(i,k)-cin1(i,k)*cjp1(i,k-1)+eps
          cjp1(i,k)=cip1(i,k)/cj(i,k)
          rhs(i,k)=(rhs(i,k)-cin1(i,k)*rhs(i,k-1))/cj(i,k)
       enddo
       !
       ! here rhs = y(k)/u(k,k), cjp1=a(k,k+1)/u(k,k)
       !
       cj(i,n1)=rhs(i,n1)
       do k=n1-1,2,-1
          cj(i,k)=rhs(i,k)-cjp1(i,k)*cj(i,k+1)
       enddo
    end do

  end subroutine tridiff
  !
  ! --------------------------------------------------------------------
  ! subroutine ae1mm: subtracts mean value from given field (a=a-a_bar)
  !
  subroutine ae1mm(n1,n2,n3,a,abar)

    use mpi_interface, only : nypg,nxpg,double_array_par_sum

    integer n1,n2,n3
    real, intent (inout), dimension (n1,n2,n3) :: a(n1,n2,n3)
    real, intent (out), dimension (n1)         :: abar(n1)

    integer :: i,j,k

    real(kind=8) :: lavg(n1),gavg(n1)

    ! TR: this used to be for the whole domain ...
    !call get_avg3(n1,n2,n3,a,abar)

    gavg(:) = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             gavg(k)=gavg(k)+a(k,i,j)
          end do
       end do
    end do
    lavg = gavg
    call double_array_par_sum(lavg,gavg,n1)
    abar(:) = real( gavg(:)/real((nypg-4)*(nxpg-4)) )

    do j=1,n3
       do i=1,n2
          do k=1,n1
             a(k,i,j)=a(k,i,j)-abar(k)
          enddo
       enddo
    enddo

  end subroutine ae1mm
  !
  !---------------------------------------------------------------------
  ! CRAYFFTUSE:  Uses the cray routines to do a 2D transform
  !
  subroutine get_fft_twodim(nx,ny,nz,a,wsavex,wsavey,isgn)

    use mpi_interface, only : xshuffle,yshuffle,nxg,nyg,nynzp,nxnzp

    integer, intent(in)    :: nx,ny,nz,isgn
    complex, intent(inout) :: a(nx,ny,nz)
    real, intent(inout)    :: wsavex(4*nxg+100),wsavey(4*nyg+100)

    integer :: k, j, i
    complex :: atmp((nx+1)*(ny+1)*(nz+1)),btmp(ny,nx,nz)

    call xshuffle(a,atmp,nx,ny,nz,1)
    call fft1dc(nxg,nynzp,atmp,wsavex,isgn,fftinix)
    call xshuffle(a,atmp,nx,ny,nz,-1)

    do k=1,nz
       do j=1,ny
          do i=1,nx
             btmp(j,i,k)=A(i,j,k)
          enddo
       enddo
    enddo

    call yshuffle(btmp,atmp,nx,ny,nz,1)
    call fft1dc(nyg,nxnzp,atmp,wsavey,isgn,fftiniy)
    call yshuffle(btmp,atmp,nx,ny,nz,-1)

    do k=1,nz
       do j=1,ny
          do i=1,nx
             A(i,j,k)= btmp(j,i,k)
          enddo
       enddo
    enddo

  end subroutine get_fft_twodim

end module util
