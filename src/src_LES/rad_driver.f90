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


  use defs, only       : cp, rcp, cpr, rowt, roice, p00, pi, nv1, nv
  use fuliou, only     : rad, rad_init, minSolarZenithCosForVis
  implicit none

  real, parameter :: SolarConstant = 1.365e+3

  ! NAMELIST parameters
  character (len=50) :: radsounding = 'datafiles/dsrt.lay'
  LOGICAL :: useMcICA = .TRUE.
  LOGICAL :: RadNewSetup = .TRUE. ! use the new radiation setup method
  REAL :: RadConstSZA = -360. ! constant solar zenith angle (values between -180 and 180 degrees)

  logical, save     :: first_time = .True.
  real, allocatable, save ::  pp(:), pt(:), ph(:), po(:), pre(:), pde(:), &
       plwc(:), piwc(:), pgwc(:), fds(:), fus(:), fdir(:), fuir(:)

  contains

    subroutine rad_new_setup(n1,pi0,th0,rt0)
      ! The new radiation setup method: initialize sounding profiles using basic state
      ! properties including pressure, temperature and total water mixing ratio. This
      ! means that the setup will be the same for restarts. This function initializes
      ! also the radiative transfer solver.
      integer, intent(in) :: n1
      real, intent(in) :: pi0(n1), th0(n1), rt0(n1) ! Basic state p, theta and rt
      real :: ttop, p0(n1)
      !
      ! Basic state pressure (hPa)
      p0(:) = 0.01*(p00*(pi0(:)/cp)**cpr)
      ! Basic state domain top temperature (K)
      ttop = th0(n1)*pi0(n1)/cp
      IF ( radsounding == 'auto' ) THEN
        ! Automatic profile generation
        !   - Upper atmosphere (P > 179 hPa) from dsrt.lay
        !   - log-log or lin-log interpolation between LES and the upper atmosphere P, T and rt
        !   - Boundary layer ozone concentration fixed to 50 ppt
        call setup_auto(n1,nv1,nv,p0,ttop,MAX(4e-7,rt0(n1)),50e-9)
      ELSE
        ! Works well with the LES model
        call setup_les(radsounding,n1,nv1,nv,p0,ttop,MAX(4e-7,rt0(n1)))
      ENDIF
      !
      ! Hydrometeor properties set to zero (nothing above LES domain)
      pre(:) = 0.
      pde(:) = 0.
      piwc(:) = 0.
      plwc(:) = 0.
      pgwc(:) = 0.
      !
      ! Initialize radiative transfer solver
      CALL rad_init
      !
      ! Setup done
      first_time = .False.
    end subroutine rad_new_setup

    subroutine d4stream(n1, n2, n3, alat, time, sknt, sfc_albedo, dn, pi0, pi1, dzm, &
         pip, tk, rv, rc, nc, tt, rflx, sflx, afus, afds, afuir, afdir, albedo, ice, nice, grp)
      integer, intent (in) :: n1, n2, n3
      real, intent (in)    :: alat, time, sknt, sfc_albedo
      real, dimension (n1), intent (in)                 :: pi0, pi1, dzm
      real, dimension (n1,n2,n3), intent (in)           :: dn, pip, tk, rv, rc, nc
      real, optional, dimension (n1,n2,n3), intent (in) :: ice, nice, grp
      real, dimension (n1,n2,n3), intent (inout)        :: tt, rflx, sflx
      real, dimension (n1+1,n2,n3), intent (inout)      :: afus, afds, afuir, afdir
      real, intent (out)                                :: albedo(n2,n3)

      integer :: kk, k, i, j, npts
      real    :: ee, u0, xfact, prw, pri, p0(n1), exner(n1), pres(n1)

      if (first_time) then
         p0(n1) = (p00*(pi0(n1)/cp)**cpr) / 100.
         p0(n1-1) = (p00*(pi0(n1-1)/cp)**cpr) / 100.
         call setup(radsounding,n1,npts,nv1,nv,p0)
         first_time = .False.
         if (allocated(pre))   pre(:) = 0.
         if (allocated(pde))   pde(:) = 0.
         if (allocated(piwc)) piwc(:) = 0.
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
      IF (-180. .LE. RadConstSZA .AND. RadConstSZA .LE. 180.) THEN
         ! Fixed solar zenith angle
         u0 = cos(RadConstSZA*pi/180.)
      ELSE
         ! Solar zenith angle based on latitude and decimal day of year
         u0 = zenith(alat,time)
      ENDIF
      !
      ! call the radiation
      !
      prw = (4./3.)*pi*rowt
      pri = (3.*sqrt(3.)/8.)*roice
      do j=3,n3-2
         do i=3,n2-2
            ! Grid cell pressures in the LES model (Pa)
            exner(1:n1) = (pi0(1:n1)+pi1(1:n1)+pip(1:n1,i,j))/cp
            pres(1:n1) = p00*( exner(1:n1) )**cpr

            ! LES and background pressure levels must be separate (LES pressure levels will change)
            !      Tomi Raatikainen 17.6.2016 & 21.12.2016
            pp(nv-n1+2) = pres(n1)/100. - 0.5*(pres(n1-1)-pres(n1)) / 100.
            if ( pp(nv-n1+2) < pp(nv-n1+1)+1.0 ) THEN
                ! Simple solution to the problem: remove sounding level nv-n1+1, which means that LES data is written over that
                nv1=nv1-1       ! This should do it (start from the previos location)
                nv=nv-1
                pp(nv-n1+2) = pres(n1)/100. - 0.5*(pres(n1-1)-pres(n1)) / 100. ! Update

                WRITE(*,*) 'Warning (radiation/d4stream): overlapping pressure levels observed (i,j)!',i,j
            endif

            do k=2,n1
               kk = nv-(k-2)
               pp(kk+1) = 0.5*(pres(k-1)+pres(k)) / 100.
               pt(kk) = tk(k,i,j)
               ph(kk) = rv(k,i,j)

               ! Cloud water
               if ((rc(k,i,j).gt.0.) .and. (nc(k,i,j).gt.0.)) THEN
                  plwc(kk) = 1000.*dn(k,i,j)*rc(k,i,j)
                  pre(kk)  = 1.e6*(plwc(kk)/(1000.*prw*nc(k,i,j)*dn(k,i,j)))**(1./3.)
                  pre(kk)=min(max(pre(kk),4.18),31.23)
               ELSE
                  pre(kk) = 0.
                  plwc(kk) = 0.
               end if

               ! Ice
               if (present(ice)) then
                  if ((ice(k,i,j).gt.0.).and.(nice(k,i,j).gt.0.)) then
                     piwc(kk) = 1000.*dn(k,i,j)*ice(k,i,j)
                     pde(kk)  = 1.e6*(piwc(kk)/(1000.*pri*nice(k,i,j)*dn(k,i,j)))**(1./3.)
                     pde(kk)=min(max(pde(kk),20.),180.)
                  else
                     piwc(kk) = 0.0
                     pde(kk)  = 0.0
                  endif
               end if

               ! Graupel
               if (present(grp)) then
                  pgwc(kk) = 1000.*dn(k,i,j)*grp(k,i,j)
               end if

            end do

            if (present(ice).and.present(grp)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, useMcICA, plwc=plwc, pre=pre, piwc=piwc, pde=pde, pgwc=pgwc)
            ELSEif (present(ice)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, useMcICA, plwc=plwc, pre=pre, piwc=piwc, pde=pde)
            else
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, useMcICA, plwc=plwc, pre=pre)
            end if

            do k=1,n1
               kk = nv1 - (k-1)
               afus(k,i,j) = fus(kk)
               afds(k,i,j) = fds(kk)
               afuir(k,i,j) = fuir(kk)
               afdir(k,i,j) = fdir(kk)
               sflx(k,i,j) = fus(kk)  - fds(kk)
               rflx(k,i,j) = sflx(k,i,j) + fuir(kk) - fdir(kk)
            end do

            ! TOA fluxes (k=n1+1)
            afus(k,i,j) = fus(1)
            afds(k,i,j) = fds(1)
            afuir(k,i,j) = fuir(1)
            afdir(k,i,j) = fdir(1)

            if (u0 > minSolarZenithCosForVis) then
               albedo(i,j) = fus(1)/fds(1)
            else
               albedo(i,j) = -999.
            end if

            do k=2,n1-3
               xfact  = dzm(k)/(cp*dn(k,i,j)*exner(k))
               tt(k,i,j) = tt(k,i,j) - (rflx(k,i,j) - rflx(k-1,i,j))*xfact
            end do

         end do
      end do

    end subroutine d4stream

  ! ---------------------------------------------------------------------------
  ! sets up the input data to extend through an atmosphere of appreiciable
  ! depth using a background souding specified as a parameter, match this to
  ! the original sounding using p0 as this does not depend on time and thus
  ! allows us to recompute the same background matching after a history start
  !
  subroutine setup(background,n1,npts,nv1,nv,zp)

    character (len=*), intent (in) :: background
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
    allocate (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),piwc(nv),pgwc(nv))

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

    character (len=*), intent (in) :: background
    integer, intent (in) :: n1
    integer, intent (out):: nv1,nv
    real, intent (in)    :: zp(n1), ttop, qtop

    real, allocatable  :: sp(:), st(:), sh(:), so(:), sl(:)

    integer :: k, ind, ns, nb, nt
    real    :: ptop, dp, Tsurf

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
    allocate (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),piwc(nv),pgwc(nv)) ! Cell centers

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
  !
  ! ---------------------------------------------------------------------------
  ! This version of setup is modified for LES simulations
  ! - Upper atmosphere (P > 179 hPa) from dsrt.lay (included here)
  ! - log-log or lin-log interpolation between LES and the upper atmosphere
  !
  ! Tomi Raatikainen, 19.4.2018
  !
  subroutine setup_auto(n1,nv1,nv,zp,ttop,qtop,otop)
    use mpi_interface, only: myid
    implicit none

    integer, intent (in) :: n1  ! LES nzp
    real, intent (in) :: zp(n1) ! LES pressure grid
    real, intent (in) :: ttop, qtop, otop ! T (K), q (kg/kg) and [O3] (-) at the LES column top
    integer, intent (out):: nv1,nv

    INTEGER, PARAMETER :: nb =22
    real  :: sp(nb), st(nb), sh(nb), so(nb)

    integer :: k, nt
    real    :: ptop, dp

    ! The first 22 values from dsrt.lay (the same as in esrt, hsrt, kmls and zh2o)
    !   p [hPa], T [K], q [kg/kg], O3 [-]
    sp = (/0.0709261, 0.13172, 0.253308, 0.48635, 0.952436, 1.76302, &
        3.33353, 6.5252, 13.172, 27.6612, 32.2207, 37.5908, 43.6702, 50.9655, &
        59.4766, 69.5076, 81.261, 95.041, 111.05, 129.997, 152.998, 179.038/)
    st = (/218, 231, 245, 260, 276, 270, 258, 245, 234, 224, 223, 222, 220, 219, &
        218, 217, 216, 216, 216, 216, 216, 216/)
    sh = (/4.32191e-07, 4.19289e-07, 4.24899e-07, 4.12278e-07, &
        4.59606e-07, 4.38849e-07, 4.08538e-07, 3.84971e-07, &
        3.57941e-07, 3.69835e-06, 3.66955e-06, 3.78973e-06, &
        3.69566e-06, 3.7091e-06, 3.69444e-06, 3.70354e-06, &
        3.71384e-06, 3.76013e-06, 3.72031e-06, 3.58904e-06, &
        3.74612e-06, 5.42635e-06/)
    so = (/6.9344e-07, 9.23243e-07, 1.2789e-06, 1.7038e-06, 2.6004e-06, &
        4.05764e-06, 6.29702e-06, 6.78494e-06 ,6.76637e-06, 6.4506e-06, &
        5.84206e-06, 5.48301e-06, 4.83796e-06, 4.14682e-06, 3.32304e-06, &
        2.66922e-06, 1.9921e-06, 1.46976e-06, 1.0866e-06, 8.09878e-07, &
        6.7769e-07, 4.90333e-07/)

    ! The highest LES domain pressure
    ptop=zp(n1)

    ! Interpolation between the lowest sounding and LES model top
    !   T=a+b*ln(P)
    !   ln(r)=a+b*ln(P)
    !   ln(o)=a+b*ln(P)
    ! Resolution at least dln(P)=0.12 (from the radiation soundings)
    dp=max( 0.12, log(zp(n1-1)/zp(n1)) )

    ! New layers between current soundings and LES model top
    nt=FLOOR( log(ptop/sp(nb))/dp )

    IF (nt<1 .AND. myid==0) THEN
        WRITE(*,*) ptop, sp(nb), dp, nt
        STOP 'Error in generating radiation profiles!'
    ENDIF

    ! Include complete LES grid (n1 points), soundings (nb points) and the transition regime (nt points)
    nv1 = n1+nb+nt  ! Total number of pressure points (cell interfaces)
    nv = nv1-1   ! Total number of scalars (cell centers)

    ! Sounding and LES data
    allocate (pp(nv1),fds(nv1),fus(nv1),fdir(nv1),fuir(nv1)) ! Cell interfaces
    allocate (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),piwc(nv),pgwc(nv)) ! Cell centers

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
            ! Logarithmic pressure grid
            pp(nb+k)=pp(nb+k-1)*exp(dp) !-(pp(nb)-ptop)/REAL(nt+1)
            ! Interpolate between LES model top and the lowest sounding (nb)
            pt(nb+k)=ttop+(st(nb)-ttop)/log(sp(nb)/ptop)*log(pp(nb+k)/ptop)
            ph(nb+k)=exp(log(qtop)+log(sh(nb)/qtop)/log(sp(nb)/ptop)*log(pp(nb+k)/ptop))
            po(nb+k)=exp(log(otop)+log(so(nb)/otop)/log(sp(nb)/ptop)*log(pp(nb+k)/ptop))
        ENDDO
    ENDIF

    ! Within the LES domain
    !  a) pressure, temperature and humidity will be obtained from the LES
    !  b) set the ozone to a constant value
    IF (nb+nt>0) po(nb+nt+1:nv) = otop

    ! Print information about sounding
    IF (myid==0) THEN
        WRITE(*,'(/,/,20A)') '-------------------------------------------------'
        WRITE(*,'(/2x,20A,/)') 'Reading Background Sounding: '//trim(radsounding)
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
                WRITE(*,'(2x,A7,4A10)') 'LES','...','...','...','...'
            ENDIF
        ENDDO
        WRITE(*,*) ' '
    ENDIF

  end subroutine setup_auto
  !
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
