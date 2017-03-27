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
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module radiation

  ! 20151022: Made some modifications for parameters going to *rad* to
  ! avoid errors with SALSA. 
  ! Juha Tonttila, FMI


  use defs, only       : cp, rcp, cpr, rowt, roice, p00, pi, nv1, nv, SolarConstant
  use fuliou, only     : rad, set_random_offset, minSolarZenithCosForVis
  implicit none

  character (len=50) :: background = 'datafiles/dsrt.lay'
  LOGICAL :: McICA = .TRUE.

  LOGICAL :: ConstPress = .FALSE.
  REAL, ALLOCATABLE, SAVE :: exner(:), pres(:)

  logical, save     :: first_time = .True.
  real, allocatable, save ::  pp(:), pt(:), ph(:), po(:), pre(:), pde(:), &
       plwc(:), piwc(:), prwc(:), pgwc(:), fds(:), fus(:), fdir(:), fuir(:)

  integer :: k,i,j, npts
  real    :: ee, u0, day, time, alat, zz

  contains

    subroutine d4stream(n1, n2, n3, alat, time, sknt, sfc_albedo, dn0, pi0, pi1, dzm, &
         pip, tk, rv, rc, nc, tt, rflx, sflx, albedo, rr, ice, nice, grp, radsounding, useMcICA, ConstPrs)
      use mpi_interface, only: myid, pecount
      integer, intent (in) :: n1, n2, n3
      real, intent (in)    :: alat, time, sknt, sfc_albedo
      real, dimension (n1), intent (in)                 :: dn0, pi0, pi1, dzm
      real, dimension (n1,n2,n3), intent (in)           :: pip, tk, rv, rc, nc
      real, optional, dimension (n1,n2,n3), intent (in) :: rr, ice, nice, grp
      real, dimension (n1,n2,n3), intent (inout)        :: tt, rflx, sflx
      CHARACTER(len=50), OPTIONAL, INTENT(in)           :: radsounding
      LOGICAL, OPTIONAL, INTENT(in)                     :: useMcICA, ConstPrs
      real, intent (out)                                :: albedo(n2,n3)

      integer :: kk
      real    :: xfact, prw, pri, p0(n1)

      IF (PRESENT(radsounding)) background = radsounding
      IF (PRESENT(useMcICA)) McICA = useMcICA
      IF (PRESENT(ConstPrs)) ConstPress = ConstPrs

      if (first_time) then
         ! Possible to use constant LES pressure levels (fixed during the first call)
         ALLOCATE(exner(n1), pres(n1))
         exner(1:n1) = (pi0(1:n1)+pi1(1:n1))/cp
         pres(1:n1) = p00*( exner(1:n1) )**cpr
         IF (.TRUE.) THEN
            ! Works well with the LES model
            p0(1:n1) = 0.01*pres(1:n1)
            call setup_les(background,n1,nv1,nv,p0,tk(n1,3,3),rv(n1,3,3))
         ELSE
            p0(n1) = (p00*(pi0(n1)/cp)**cpr) / 100.
            p0(n1-1) = (p00*(pi0(n1-1)/cp)**cpr) / 100.
            call setup(background,n1,npts,nv1,nv,p0)
         ENDIF
         first_time = .False.
         if (allocated(pre))   pre(:) = 0.
         if (allocated(pde))   pde(:) = 0.
         if (allocated(piwc)) piwc(:) = 0.
         if (allocated(prwc)) prwc(:) = 0.
         if (allocated(plwc)) plwc(:) = 0.
         if (allocated(pgwc)) pgwc(:) = 0.
      end if
      !
      ! initialize surface albedo, emissivity and skin temperature.
      !
      ee = 1.0
      !
      ! determine the solar geometery, as measured by u0, the cosine of the
      ! solar zenith angle
      !
      u0 = zenith(alat,time)
      !
      ! Avoid identical random numbers by adding an offset for each PU
      !     myid=0,1,2,...
      !     (n3-4)*(n2-4) is the number of random numbers for each PU
      !     IR and optionally also visible wavelengths
      ! First, call random numbers to add offset for PUs before myid
      IF (McICA .and. myid>0) THEN
        if (u0 > minSolarZenithCosForVis) THEN
            kk=myid*(n3-4)*(n2-4)*2
        ELSE
            kk=myid*(n3-4)*(n2-4) ! Without solar wavelengths
        ENDIF
        CALL set_random_offset( kk )
      ENDIF
      !
      ! call the radiation
      !
      prw = (4./3.)*pi*rowt
      pri = (3.*sqrt(3.)/8.)*roice
      do j=3,n3-2
         do i=3,n2-2
            IF (.NOT.ConstPress) THEN
                ! Grid cell pressures in the LES model (Pa)
                exner(1:n1) = (pi0(1:n1)+pi1(1:n1)+pip(1:n1,i,j))/cp
                pres(1:n1) = p00*( exner(1:n1) )**cpr
            ENDIF

            ! LES and background pressure levels must be separate (LES pressure levels will change)
            !      Tomi Raatikainen 17.6.2016 & 21.12.2016
            pp(nv-n1+2) = pres(n1)/100. - 0.5*(pres(n1-1)-pres(n1)) / 100.
            if ( pp(nv-n1+2) < pp(nv-n1+1)+1.0 ) THEN
                ! Simple solution to the problem: remove sounding level nv-n1+1, which means that LES data is written over that
                nv1=nv1-1       ! This should do it (start from the previos location)
                nv=nv-1
                pp(nv-n1+2) = pres(n1)/100. - 0.5*(pres(n1-1)-pres(n1)) / 100. ! Update

                WRITE(*,*) 'Warning (radiation/d4stream): ovelapping pressure levels observed (i,j)!',i,j
            endif

            do k=2,n1
               kk = nv-(k-2)
               pp(kk+1) = 0.5*(pres(k-1)+pres(k)) / 100.
               pt(kk) = tk(k,i,j)
               ph(kk) = rv(k,i,j)

               ! Cloud water
               if ((rc(k,i,j).gt.0.) .and. (nc(k,i,j).gt.0.)) THEN
                  plwc(kk) = 1000.*dn0(k)*rc(k,i,j)
                  pre(kk)  = 1.e6*(plwc(kk)/(1000.*prw*nc(k,i,j)*dn0(k)))**(1./3.)
                  pre(kk)=min(max(pre(kk),4.18),31.23)
               ELSE
                  pre(kk) = 0.
                  plwc(kk) = 0.
               end if

               ! Precipitation (not used at the moment)
               if (present(rr)) then
                  prwc(kk) = 1000.*dn0(k)*rr(k,i,j)
               end if

               ! Ice
               if (present(ice)) then
                  if ((ice(k,i,j).gt.0.).and.(nice(k,i,j).gt.0.)) then
                     piwc(kk) = 1000.*dn0(k)*ice(k,i,j)
                     pde(kk)  = 1.e6*(piwc(kk)/(1000.*pri*nice(k,i,j)*dn0(k)))**(1./3.)
                     pde(kk)=min(max(pde(kk),20.),180.)
                  else
                     piwc(kk) = 0.0
                     pde(kk)  = 0.0
                  endif
               end if

               ! Graupel
               if (present(grp)) then
                  pgwc(kk) = 1000.*dn0(k)*grp(k,i,j)
               end if

            end do

            if (present(ice).and.present(rr).and.present(grp)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, McICA, plwc=plwc, pre=pre, piwc=piwc, pde=pde, prwc=prwc, pgwc=pgwc)
            ELSEif (present(ice).and.present(grp)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, McICA, plwc=plwc, pre=pre, piwc=piwc, pde=pde, pgwc=pgwc)
            ELSEif (present(ice).and.present(rr)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, McICA, plwc=plwc, pre=pre, piwc=piwc, pde=pde, prwc=prwc)
            ELSEif (present(ice)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, McICA, plwc=plwc, pre=pre, piwc=piwc, pde=pde)
            ELSEif (present(rr)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, McICA, plwc=plwc, pre=pre, prwc=prwc)
            else
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, McICA, plwc=plwc, pre=pre)
            end if

            do k=1,n1
               kk = nv1 - (k-1)
               sflx(k,i,j) = fus(kk)  - fds(kk)
               rflx(k,i,j) = sflx(k,i,j) + fuir(kk) - fdir(kk)
            end do

            if (u0 > minSolarZenithCosForVis) then
               albedo(i,j) = fus(1)/fds(1)
            else
               albedo(i,j) = -999.
            end if

            do k=2,n1-3
               xfact  = dzm(k)/(cp*dn0(k)*exner(k))
               tt(k,i,j) = tt(k,i,j) - (rflx(k,i,j) - rflx(k-1,i,j))*xfact
            end do

         end do
      end do

      ! Call random numbers to add offset for PUs after myid
      IF (McICA .and. myid<pecount-1) THEN
        if (u0 > minSolarZenithCosForVis) THEN
            kk=(pecount-1-myid)*(n3-4)*(n2-4)*2
        ELSE
            kk=(pecount-1-myid)*(n3-4)*(n2-4) ! Without solar wavelengths
        ENDIF
        CALL set_random_offset( kk )
      ENDIF
    end subroutine d4stream

  ! ---------------------------------------------------------------------------
  ! sets up the input data to extend through an atmosphere of appreiciable
  ! depth using a background souding specified as a parameter, match this to
  ! the original sounding using p0 as this does not depend on time and thus
  ! allows us to recompute the same background matching after a history start
  !
  subroutine setup(background,n1,npts,nv1,nv,zp)

    character (len=19), intent (in) :: background
    integer, intent (in) :: n1
    integer, intent (out):: npts,nv1,nv
    real, intent (in)    :: zp(n1)

    real, allocatable  :: sp(:), st(:), sh(:), so(:), sl(:)

    integer :: k, ns, norig, index
    logical :: blend
    real    :: pa, pb, ptop, ptest, test, dp1, dp2, dp3, Tsurf

    open ( unit = 08, file = background, status = 'old' )
    print *, 'Reading Background Sounding: ',background
    read (08,*) Tsurf, ns
    allocate ( sp(ns), st(ns), sh(ns), so(ns), sl(ns))
    do k=1,ns
       read ( 08, *) sp(k), st(k), sh(k), so(k), sl(k)
    enddo
    close (08)

    !
    ! identify what part, if any, of background sounding to use
    !
    ptop = zp(n1)
    if (sp(2) < ptop) then
       pa = sp(1)
       pb = sp(2)
       k = 3
       do while (sp(k) < ptop .and. k<ns)
          pa = pb
          pb = sp(k)
          k  = k+1
       end do
       k=k-1           ! identify first level above top of input
       blend = .True.
    else
       blend = .False.
    end if
    !
    ! if blend is true then the free atmosphere above the sounding will be
    ! specified based on the specified background climatology, here the
    ! pressure levels for this part of the sounding are determined
    !
    if (blend) then
       dp1 = pb-pa
       dp2 = ptop - pb
       dp3 = zp(n1-1) - zp(n1)
       if (dp1 > 2.*dp2) k = k-1 ! first level is too close, blend from prev
       npts  = k
       norig = k
       ptest = sp(k)
       test = ptop-ptest
       do while (test > 2*dp3)
          ptest = (ptest+ptop)*0.5
          test  = ptop-ptest
          npts  = npts + 1
       end do
       nv1 = npts + n1
    else
       nv1 = n1
    end if
    nv = nv1-1
    !
    ! allocate the arrays for the sounding data to be used in the radiation 
    ! profile and then fill them first with the sounding data, by afill, then
    ! by interpolating the background profile at pressures less than the
    ! pressure at the top of the sounding
    !
    allocate (pp(nv1),fds(nv1),fus(nv1),fdir(nv1),fuir(nv1))
    allocate (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),prwc(nv),piwc(nv),pgwc(nv))

    if (blend) then
       pp(1:norig) = sp(1:norig)
       pt(1:norig) = st(1:norig)
       ph(1:norig) = sh(1:norig)
       po(1:norig) = so(1:norig)

      do k=norig+1,npts
          pp(k) = (ptop + pp(k-1))*0.5
          index = getindex(sp,ns,pp(k))
          pt(k) =  intrpl(sp(index),st(index),sp(index+1),st(index+1),pp(k))
          ph(k) =  intrpl(sp(index),sh(index),sp(index+1),sh(index+1),pp(k))
          po(k) =  intrpl(sp(index),so(index),sp(index+1),so(index+1),pp(k))
       end do
       !
       ! set the ozone constant below the reference profile
       !
       do k=npts+1,nv
          po(k) =  po(npts)
       end do
    end if

  end subroutine setup
  !
  ! ---------------------------------------------------------------------------
  ! This version of setup is modified for LES simulations
  !    - Always using background soundings and LES data
  !    - Separate LES and background profiles
  !    - Weighting of the backgound towards smooth transition at the interface
  !
  ! Tomi Raatikainen, 22.12.2016
  !
  subroutine setup_les(background,n1,nv1,nv,zp,ttop,qtop)
    use mpi_interface, only: myid
    implicit none

    character (len=19), intent (in) :: background
    integer, intent (in) :: n1
    integer, intent (out):: nv1,nv
    real, intent (in)    :: zp(n1), ttop, qtop

    real, allocatable  :: sp(:), st(:), sh(:), so(:), sl(:)

    integer :: k, ind, ns, nb, nt
    real    :: ptop, dp, ptmp, Tsurf

    ! Read backgroud profiles (p [hPa], T [K], q [kg/kg], O3 [-] and an unused concentration)
    open ( unit = 08, file = background, status = 'old' )
    read (08,*) Tsurf, ns
    allocate ( sp(ns), st(ns), sh(ns), so(ns), sl(ns))
    do k=1,ns
       read ( 08, *) sp(k), st(k), sh(k), so(k), sl(k)
    enddo
    close (08)

    ! Merge background and LES pressure levels (p_LES > p_background)
    !
    ! The highest LES model pressure - defined for cell interface
    ptop=zp(n1)-0.5*(zp(n1-1)-zp(n1))
    !
    ! The first sounding level must be somewhat higher (lower pressure) then the last LES pressure level
    !  - Must avoid overlapping pressure levels when LES pressures change with time
    !  - 10 hPa should be large enough for most cases (LES has typically high pressure resolution)
    dp=MAX(10.0,zp(n1-1)-zp(n1))
    !
    k=1 ! Note: increasing pressure and decreasing altitude
    DO WHILE (sp(k) < ptop-dp )
        k=k+1
        IF (k==ns+1) EXIT
    END DO
    nb=k-1 ! Background soundings from 1 to k (nb points)

    ! Add layers between LES and soundings (mind the gap!)
    nt=0 ! Transition regime
    IF (nb>0) nt=FLOOR( (ptop-sp(nb))/dp )-1

    ! Include complete LES grid (n1 points), soundings (nb points) and the transition regime (nt points)
    nv1 = n1+nb+nt  ! Total number of pressure points (cell interfaces)
    nv = nv1-1   ! Total number of scalars (cell centers)

    ! allocate the arrays for the sounding data to be used in the radiation
    ! profile and then fill them first with the sounding data, by afill, then
    ! by interpolating the background profile at pressures less than the
    ! pressure at the top of the sounding
    !
    ! Note: LES pressure levels are increasing, but radiation calculations
    ! expect decreasing pressure grid (from TOA to surface)
    !
    allocate (pp(nv1),fds(nv1),fus(nv1),fdir(nv1),fuir(nv1)) ! Cell interfaces
    allocate (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),prwc(nv),piwc(nv),pgwc(nv)) ! Cell centers

    po=0.
    IF (nb>0) THEN
        ! Levels above LES domain: copy background soundings
        pp(1:nb) = sp(1:nb)
        pt(1:nb) = st(1:nb)
        ph(1:nb) = sh(1:nb)
        po(1:nb) = so(1:nb)
    ENDIF

    IF (nt>0) THEN
        ! Interpolated levels between pp(nb) and ptop
        DO k=1,nt
            pp(nb+k)=pp(nb+k-1)-(pp(nb)-ptop)/REAL(nt+1)
            ! Interpolate between LES model top and the lowest sounding
            !   y=y0+(y1-y0)(p1-p0)*(p-p0)
            ind = getindex(sp,ns,pp(nb+k)) ! Typically ind=nb
            pt(nb+k)=ttop+(st(ind)-ttop)/(sp(ind)-ptop)*(pp(nb+k)-ptop)
            ph(nb+k)=qtop+(sh(ind)-qtop)/(sp(ind)-ptop)*(pp(nb+k)-ptop)
            ! Interpolate ozone using background soundings
            po(nb+k)=so(ind+1)+(so(ind)-so(ind+1))/(sp(ind)-sp(ind+1))*(pp(nb+k)-sp(ind+1))
        ENDDO
    ENDIF

    ! Within the LES domain
    !  a) pressure, temperature and humidity will be obtained from the LES
    !  b) set the ozone to a constant value
    IF (nb+nt>0) po(nb+nt+1:nv) = po(nb+nt)

    ! Print information about sounding
    IF (myid==0) THEN
        WRITE(*,'(/,/,20A)') '-------------------------------------------------'
        WRITE(*,'(/2x,20A,/)') 'Reading Background Sounding: '//trim(background)
        WRITE(*,'(/2x,A7,4A10)') 'Source','p (hPa)','T (K)','q (g/kg)','O3 (ppm)'

        ! Calculate LES pressure levels ([zp]=hPa) and T and q for cell centers
        pp(nv-n1+2) = zp(n1)/100. - 0.5*(zp(n1-1)-zp(n1))
        do k=2,n1
             pp(nv-(k-2)+1) = 0.5*(zp(k-1)+zp(k))
        ENDDO

        DO k=1,nv
            IF (k<=nb .AND. nb>0) THEN
                ! Data from backgound soundings
                WRITE(*,'(2x,A7,2F10.2,2F10.3)') 'Backgr',pp(k),pt(k),ph(k)*1.e3,po(k)*1.e6
            ELSEIF (k<=nb+nt .AND. nt>0) THEN
                ! Interpolation
                WRITE(*,'(2x,A7,2F10.2,2F10.3)') 'Interp',pp(k),pt(k),ph(k)*1.e3,po(k)*1.e6
            ELSEIF (k==nb+nt+1) THEN
                ! The highest LES layer - other layers not yet available
                WRITE(*,'(2x,A7,2F10.2,2F10.3)') 'LES',ptop,ttop,qtop*1.e3,po(k)*1.e6
                !WRITE(*,'(2x,A7,3A10,F10.3)') 'LES','*','*','*',po(k)*1.e6
                WRITE(*,'(2x,A7,4A10)') 'LES','...','...','...','...'
            ENDIF
        ENDDO
        WRITE(*,*) ' '
    ENDIF

  end subroutine setup_les
  ! ---------------------------------------------------------------------------
  !  locate the index closest to a value
  !
  integer function getindex(x,n,xval)

    integer, intent (in)  :: n
    real,    intent (in)  :: x(n),xval

    integer :: ia, ib

    ia=1
    ib=n
    if (xval < x(1)) then
       getindex = 1
    elseif (xval > x(n)) then
       getindex = n-1
    else
       getindex = (ia+ib)/2
       do while (getindex /= ia .or. ib /= getindex+1)
          getindex = (ia+ib)/2
          if ((xval-x(getindex)) >= 0.0) then
             ia = getindex
          else
             ib = getindex
          end if
       end do
    endif

  end function getindex

  ! ---------------------------------------------------------------------------
  ! linear interpolation between two points,
  !
  real function intrpl(x1,y1,x2,y2,x)

    real, intent (in)  :: x1,y1,x2,y2,x

    real :: slope

    slope  = (y2-y1)/(x2 - x1 + epsilon(1.))
    intrpl = y1+slope*(x-x1)

  end function intrpl

  ! ---------------------------------------------------------------------------
  ! Return the cosine of the solar zenith angle give the decimal day and
  ! the latitude
  !
  real function zenith(alat,time)

    real, intent (in)  :: alat, time

    real :: lamda, d, sig, del, h, day

    day    = floor(time)
    lamda  = alat*pi/180.
    d      = 2.*pi*int(time)/365.
    sig    = d + pi/180.*(279.9340 + 1.914827*sin(d) - 0.7952*cos(d) &
         &                      + 0.019938*sin(2.*d) - 0.00162*cos(2.*d))
    del    = asin(sin(23.4439*pi/180.)*sin(sig))
    h      = 2.*pi*((time-day)-0.5)
    zenith = sin(lamda)*sin(del) + cos(lamda)*cos(del)*cos(h)

  end function zenith

end module radiation
