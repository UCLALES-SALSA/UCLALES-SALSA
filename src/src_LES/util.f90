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
  character (len=6), parameter :: vel_bc = 'frslip'

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
    if (vel_bc == 'noslip') then
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
  ! GET_AVG: gets average field value from k'th level
  !
  real function get_avg(n1,n2,n3,k,a)

    integer, intent (in):: n1, n2, n3, k
    real, intent (in)   :: a(n1,n2,n3)

    integer :: i,j

    get_avg=0.
    do j=3,n3-2
       do i=3,n2-2
          get_avg = get_avg + a(k,i,j)
       enddo
    enddo
    get_avg = get_avg/real((n3-4)*(n2-4))

  end function get_avg
  !
  !---------------------------------------------------------------------
  ! GET_AVG2dh: Get the average of a 2 dimensional (horizontal) input field
  !
  REAL FUNCTION get_avg2dh(n2,n3,a)

    INTEGER, INTENT(in) :: n2,n3
    REAL, INTENT(in)    :: a(n2,n3)

    INTEGER :: i,j
    
    get_avg2dh = 0.
    DO j = 3,n3-2
       DO i = 3,n2-2
          get_avg2dh = get_avg2dh + a(i,j)
       END DO
    END DO

    get_avg2dh = get_avg2dh/real((n3-4)*(n2-4))

  END FUNCTION get_avg2dh
  !
  !---------------------------------------------------------------------
  ! GET_AVG_ts: gets average field value from the whole domain 
  ! Implemented by Zubair Maalick
  ! Possibility for conditional sampling added by Juha Tonttila
  ! 
  ! Weighting by layer thickness implemented for non-uniform vertical resolution
  !
  real function get_avg_ts(n1,n2,n3,a,dz,cond)

    integer, intent (in):: n1, n2, n3
    REAL, INTENT(in)    :: dz(n1)  ! Reciprocal of layer depth!
    real, intent (in)   :: a(n1,n2,n3)
    LOGICAL, OPTIONAL, INTENT(in) :: cond(n1,n2,n3)
    
    integer :: i,j,k, npnt
    REAL :: ztmp,ztot

    npnt = 0
    IF (PRESENT(cond)) THEN
       get_avg_ts=0.
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
             ! Grid weighted vertical average
             if (ztot /=0.0 ) ztmp = ztmp/ztot
             npnt = npnt + 1 ! Note: this should not be used (use ztot instead!)
             get_avg_ts = get_avg_ts + ztmp
          END DO
       END DO
       
    ELSE
       get_avg_ts=0.
       DO j=3,n3-2
          DO i=3,n2-2
             ztot = 0.
             ztmp = 0.
             DO k = 2,n1
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

    get_avg_ts = get_avg_ts/max(10e-5,real(npnt))

  end function get_avg_ts
  !
  !---------------------------------------------------------------------
  ! GET_AVG3: gets average across outer two dimensions at each
  ! point along inner dimension - calculated over all PUs
  !
  subroutine get_avg3(n1,n2,n3,a,avg)

    use mpi_interface, only : nypg,nxpg,double_array_par_sum

    integer,intent(in) :: n1,n2,n3
    real,intent(in) :: a(n1,n2,n3)
    real,intent(out) :: avg(n1)

    integer      :: k,i,j
    real(kind=8) :: lavg(n1),gavg(n1),x

    x = 1./(real(nypg-4)*real(nxpg-4))
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

  end subroutine get_avg3
  !
  !---------------------------------------------------------------------
  ! function get_cavg: gets conditionally average field from kth level
  !
  real function get_cavg(n1,n2,n3,k,a,x)

    integer, intent (in) :: n1,n2,n3,k
    real, intent (in)    :: a(n1,n2,n3),x(n2,n3)

    integer :: i,j
    real    :: cnt

    get_cavg=0.
    cnt=0.
    do j=3,n3-2
       do i=3,n2-2
          get_cavg=get_cavg+a(k,i,j)*x(i,j)
          cnt=cnt+x(i,j)
       end do
    end do
    IF (cnt /= 0.0) get_cavg=get_cavg/cnt

  end function get_cavg
  !
  !---------------------------------------------------------------------
  ! function get_csum: conditionally sum field a over indicator x
  !
  real function get_csum(n1,n2,n3,k,a,x)

    integer, intent (in) :: n1,n2,n3,k
    real, intent (in)    :: a(n1,n2,n3),x(n2,n3)

    integer :: i,j

    get_csum=0.
    do j=3,n3-2
       do i=3,n2-2
          get_csum=get_csum+a(k,i,j)*x(i,j)
       enddo
    enddo

  end function get_csum
  !
  !---------------------------------------------------------------------
  ! function get_cor: gets mean correlation between two fields at a 
  ! given level
  !
  real function get_cor(n1,n2,n3,k,a,b)

    integer, intent (in) :: n1,n2,n3,k
    real, intent (inout) :: a(n1,n2,n3),b(n1,n2,n3)

    integer :: i,j

    get_cor=0.
    do j=3,n3-2
       do i=3,n2-2
          get_cor=get_cor+a(k,i,j)*b(k,i,j)
       end do
    end do
    get_cor=get_cor/real((n3-4)*(n2-4))

  end function get_cor
  !
  !---------------------------------------------------------------------
  ! function get_cor3: gets mean correlation accross outer two dimensions
  ! at each point along inner dimension
  !
  subroutine get_cor3(n1,n2,n3,a,b,avg)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: a(n1,n2,n3),b(n1,n2,n3)
    real, intent (out)   :: avg(n1)

    integer :: k,i,j

    avg(:) = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             avg(k)=avg(k)+a(k,i,j)*b(k,i,j)
          end do
       enddo
    enddo
    avg(:)=avg(:)/real((n3-4)*(n2-4))

  end subroutine get_cor3
  !
  !---------------------------------------------------------------------
  ! function get_var3: gets variance for a field whose mean is known
  !
  subroutine get_var3(n1,n2,n3,a,b,avg)

    integer n1,n2,n3,k,i,j
    real a(n1,n2,n3),b(n1),avg(n1)

    avg(:) = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             avg(k)=avg(k)+(a(k,i,j)-b(k))**2
          end do
       enddo
    enddo
    avg(:)=avg(:)/real((n3-4)*(n2-4))

  end subroutine get_var3
  !
  !---------------------------------------------------------------------
  ! function get_3rd3: gets the third moment for a field whose mean is known
  !
  subroutine get_3rd3(n1,n2,n3,a,b,avg)

    integer n1,n2,n3,k,i,j
    real a(n1,n2,n3),b(n1),avg(n1)

    avg(:) = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             avg(k)=avg(k)+(a(k,i,j)-b(k))**3
          end do
       enddo
    enddo
    avg(:)=avg(:)/real((n3-4)*(n2-4))

  end subroutine get_3rd3
  ! 
  !-------------------------------------------------------------------
  ! function get_var: gets square of a field from k'th level 
  ! 
  real function get_sqr(n1,n2,n3,k,a)

    integer, intent (in):: n1, n2, n3, k
    real, intent (in)   :: a(n1,n2,n3)

    integer :: i,j
    real    :: avg

    avg=0.
    do j=3,n3-2
       do i=3,n2-2
          avg=avg+a(k,i,j)**2
       end do
    end do
    get_sqr = avg/real((n3-4)*(n2-4))

  end function get_sqr
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
  !
  !---------------------------------------------------------------
  ! MASKACTIV: Create a logical mask for grid points where cloud
  !            activation will be calculated.
  !
  ! Juha Tonttila, FMI, 2014
  !
  SUBROUTINE maskactiv(act_mask,nx,ny,nz,nbins,mode,prtcl,rh,    &
                       rc,pa_naerop, pa_maerop, pt, Rpwet, w, pa_ncloud  )
    USE mo_submctl, ONLY : rhowa, rhosu, rhooc, rhoss, mwa, msu, moc, mss, pi6, nlim
    USE class_ComponentIndex, ONLY : ComponentIndex,GetIndex,IsUsed
    IMPLICIT NONE

    INTEGER, INTENT(in) :: nx,ny,nz,nbins
    REAL, INTENT(in)    :: rh(nz,nx,ny)

    REAL, OPTIONAL, INTENT(in) :: pa_naerop(nz,nx,ny,nbins),    &
                                  pa_maerop(nz,nx,ny,8*nbins),  &
                                  pt(nz,nx,ny),                 &
                                  Rpwet(nz,nx,ny,nbins),        &
                                  w(nz,nx,ny),                  &
                                  rc(nz,nx,ny),                 &
                                  pa_ncloud(nz,nx,ny,nbins)
    
    INTEGER, INTENT(in) :: mode ! 1 = Initialization; 2 = Normal timestepping
    TYPE(ComponentIndex), INTENT(in) :: prtcl


    LOGICAL, INTENT(out) :: act_mask(nz,nx,ny)

    LOGICAL :: actmask_newcloud(nz,nx,ny)
    LOGICAL :: actmask_oldcloud(nz,nx,ny)
    LOGICAL :: cldmask(nz,nx,ny)
    LOGICAL :: cldm1(nz,nx,ny),   & !inverse cloud mask offset downwards by one grid level
               cldp1(nz,nx,ny),   & !inverse cloud mask offset upwards by one grid level 
               cldpm(nz,nx,ny)

    REAL :: nwpure ! number of moles for 1 micron droplet of pure water
    REAL :: nwkelvin ! kelvin effect for 1 micron droplet
    REAL :: nwtrue(nbins)  ! moles of water minus the core particle
    REAL :: nsaero(nbins) ! moles of solute
    REAL :: vsaero(nbins) ! Volume of solute
    REAL :: nwact(nz,nx,ny) ! Minimum activity from the aerosol bins assuming a 1 micron droplet for each gridpoint
    REAL :: nwactbin(nbins) ! Activities for each aerosol bin

    LOGICAL :: notused(nx,ny)

    INTEGER :: j, i, k, b, nc

    actmask_newcloud = .FALSE.
    actmask_oldcloud = .FALSE.
    cldm1 = .TRUE.
    cldp1 = .TRUE.
    cldpm = .TRUE.
    act_mask = .FALSE. 

    nwactbin = 0.

    nwpure = pi6*((1.e-6)**3)*rhowa ! MAss of pure water droplet


    SELECT CASE(mode)
    CASE(1)
       ! Calculate cloud activation in all cloudy grid points for initialization
       ! i.e. points with RH >= threshold
       !zrh = rv/rs
       act_mask = ( rh >= 1.000 )

    CASE(2)

       ! Normal opperation
       
       IF ( .NOT. (PRESENT(rc) .AND. PRESENT(Rpwet)) ) STOP 'maskactiv: invalid arguments'

       ! Mask grid points just below cloud or where new cloud is expected to form
       ! 
       
       ! Check critical RH for 1 micron droplets
       DO j = 3,ny-2
          DO i = 3,nx-2
             DO k = 1,nz
                ! -- Get the water content and the soluble material content for imaginary 1 micron droplets
                ! -- Kelvin effect for 1 micron droplets
                nwkelvin = EXP( 4.*0.073*mwa /           &  
                                (8.314*pt(k,i,j)*rhowa*(1.e-6)) )

                nsaero(:) = 0.
                vsaero(:) = 0.
                DO b = 1,nbins
                   IF (pa_naerop(k,i,j,b) > nlim) THEN

                      IF ( IsUsed(prtcl,'SO4') ) THEN
                         nc = GetIndex(prtcl,'SO4')
                         nc = (nc-1)*nbins+b
                         vsaero(b) = vsaero(b) + pa_maerop(k,i,j,nc)
                         nsaero(b) = nsaero(b) + 3.*(pa_maerop(k,i,j,nc)/msu)     
                      END IF
                      
                      IF ( IsUsed(prtcl,'OC') ) THEN
                         nc = GetIndex(prtcl,'OC')
                         nc = (nc-1)*nbins+b
                         vsaero(b) = vsaero(b) + pa_maerop(k,i,j,nc)
                         nsaero(b) = nsaero(b) + pa_maerop(k,i,j,nc)/moc
                      END IF

                      IF ( IsUsed(prtcl,'SS') ) THEN
                         nc = GetIndex(prtcl,'SS')
                         nc = (nc-1)*nbins+b
                         vsaero(b) = vsaero(b) + pa_maerop(k,i,j,nc)
                         nsaero(b) = nsaero(b) + 2.*(pa_maerop(k,i,j,nc)/mss)
                      END IF
                      
                      nsaero(b) = nsaero(b)/pa_naerop(k,i,j,b)
                      vsaero(b) = vsaero(b)/pa_naerop(k,i,j,b)

                   END IF
                END DO

                ! True mass of water in a droplet
                nwtrue(1:nbins) = nwpure - vsaero(1:nbins)
                ! Convert to number of moles of water in a droplet
                nwtrue(1:nbins) = nwtrue(1:nbins)/mwa

                nwactbin(:) = 999. ! undefined
                WHERE(pa_naerop(k,i,j,1:nbins) > nlim) &
                     nwactbin(1:nbins) = nwkelvin*nwtrue(1:nbins)/( nwtrue(1:nbins) + nsaero(1:nbins) )

                nwact(k,i,j) = MINVAL(nwactbin(:))

             END DO ! k
          END DO ! i
       END DO ! j
       
       cldmask(:,:,:) = ( rc(:,:,:) >= 1.e-5 )

      
       ! Check for the lowest point where RH > 1. 
       cldpm(:,:,:) = .FALSE.
       cldpm(1:nz-1,:,:) = ( rh(2:nz,:,:) > 1.000 )

       cldp1(:,:,:) = .TRUE.
       cldp1(2:nz,:,:) = ( .NOT. cldpm(1:nz-1,:,:) )
       actmask_newcloud(:,:,:) = ( cldpm(:,:,:) .AND. cldp1(:,:,:) ) .AND. ( w(:,:,:) > 0. ) &
            .AND. (rc(:,:,:)<5.e-5)

       ! Base of existing cloud
       cldpm(:,:,:) = .FALSE.
       cldpm(1:nz-1,:,:) = cldmask(2:nz,:,:)
       
       cldp1(:,:,:) = .TRUE.
       cldp1(2:nz,:,:) = ( .NOT. cldpm(1:nz-1,:,:) ) 

       ! Take the lowest level of the two cases
       
       notused = .TRUE.
       DO k = 2,nz-1
          ! New cloud + no old cloud
          act_mask(k,:,:) = MERGE( (actmask_oldcloud(k,:,:) .OR. actmask_newcloud(k,:,:))  &
                                .AND. cldp1(k,:,:), .FALSE., notused(:,:))
          notused(:,:) = notused(:,:) .AND. .NOT. act_mask(k,:,:)
          notused(:,:) =  (rh(k,:,:) < 0.99 .OR. notused(:,:))
       END DO

    END SELECT
     
  END SUBROUTINE maskactiv



end module util
