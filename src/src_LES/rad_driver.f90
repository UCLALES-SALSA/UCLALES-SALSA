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


  use defs, only       : cp, rcp, cpr, rowt, roice, p00, pi, nv1, nv, ep2, Rd, g
  use fuliou, only     : rad, rad_init, minSolarZenithCosForVis, rad_tau
  implicit none

  real, parameter :: SolarConstant = 1.365e+3

  ! NAMELIST parameters
  character (len=50) :: radsounding = 'datafiles/dsrt.lay'
  LOGICAL :: useMcICA = .TRUE.
  LOGICAL :: RadNewSetup = .TRUE. ! use the new radiation setup method
  REAL :: RadConstSZA = -360. ! constant solar zenith angle (values between -180 and 180 degrees)
  REAL :: rad_lwp=0.0, rad_reff=10.0 ! Optional clouds to be added above LES domain
  REAL :: rad_iwp=0.0, rad_ieff=100.0
  INTEGER :: rad_nlev=2, rad_ilev=2
  REAL :: radsounding_pmax = -1.0 ! Optional pressure limit (hPa) for soundings

  INTEGER :: nt=0, nb=0 ! For radsounding_pmax>0

  logical, save :: calc_od = .FALSE.  ! Calculate optical depths
  real, allocatable, save :: tau_gas(:,:), tau_cloud(:,:), tau_rain(:,:), &
       tau_ice(:,:), tau_grp(:,:)

  logical, save     :: first_time = .True.
  real, allocatable, save ::  pp(:), pt(:), ph(:), po(:), pre(:), pde(:), &
       plwc(:), piwc(:), prwc(:), pgwc(:), fds(:), fus(:), fdir(:), fuir(:)

  contains

    subroutine rad_new_setup(n1,pp0,th0,rt0,zm,zt)
      ! The new radiation setup method: initialize sounding profiles using basic state
      ! properties including pressure, temperature and total water mixing ratio. This
      ! means that the setup will be the same for restarts. This function initializes
      ! also the radiative transfer solver.
      integer, intent(in) :: n1
      real, intent(in) :: pp0(n1), th0(n1), rt0(n1) ! Basic state p, theta and rt
      real, intent(in) :: zm(n1), zt(n1) ! Grid
      real :: ttop, p0(n1), dz, t0(n1)
      integer :: i
      !
      ! Basic state pressure (hPa)
      p0(:) = 0.01*pp0(:)
      ! Basic state temperature (K)
      t0(:) = th0(:)*(pp0(:)/p00)**rcp
      ttop = t0(n1)
      IF ( radsounding == 'auto' ) THEN
        ! Automatic profile generation
        !   - Upper atmosphere (P > 179 hPa) from dsrt.lay
        !   - log-log or lin-log interpolation between LES and the upper atmosphere P, T and rt
        !   - Boundary layer ozone concentration fixed to 50 ppt
        call setup_auto(n1,nv1,nv,p0,ttop,MAX(4e-7,rt0(n1)),50e-9)
      ELSEIF ( radsounding_pmax>0.0 ) THEN
        ! Similar to the previous, but
        !  - Interpolates sounding scalars for mid points
        !  - Linear region can be set between LES domain and pressure radsounding_pmax (hPa)
        !  - Updates the linear region for each column
        call setup_les_regrid(radsounding,n1,nv1,nv,p0,t0,rt0,zm,zt)
      ELSE
        ! Works well with the LES model
        call setup_les(radsounding,n1,nv1,nv,p0,ttop,MAX(4e-7,rt0(n1)))
      ENDIF
      !
      ! Hydrometeor properties set to zero (nothing above LES domain)
      pre(:) = 0.
      pde(:) = 0.
      piwc(:) = 0.
      prwc(:) = 0.
      plwc(:) = 0.
      pgwc(:) = 0.
      !
      ! Optional background clouds above LES domain
      IF (rad_lwp>1e-10) THEN ! Liquid cloud
        ! LES fills values from level nv-n1+2 onwards
        i=(nv-n1+1)-rad_nlev ! Leave rad_nlev levels between LES and the cloud
        !   LWP[g/m2] = LWC[g/m3]*dz[m]
        dz = (Rd/g) * pt(i)*(1.0+ep2*ph(i)) * alog( pp(i+1) / pp(i) ) ! dz (m)
        plwc(i) = rad_lwp/dz ! g/m3
        pre(i) = min(max(rad_reff,4.18),31.23) ! um
      ENDIF
      IF (rad_iwp>1e-10) THEN ! Ice cloud
        i=(nv-n1+1)-rad_ilev ! Leave rad_ilev levels between LES and the cloud
        dz = (Rd/g) * pt(i)*(1.0+ep2*ph(i)) * alog( pp(i+1) / pp(i) ) ! dz (m)
        piwc(i) = rad_iwp/dz ! g/m3
        pde(i) = min(max(rad_ieff,20.),120.) ! um
      ENDIF
      !
      ! Initialize radiative transfer solver
      CALL rad_init
      !
      ! Setup done
      first_time = .False.
    end subroutine rad_new_setup

    subroutine d4stream(n1, n2, n3, alat, time, sknt, sfc_albedo, dn, pi0, pi1, zm, zt, &
         pip, tk, rv, rc, nc, tt, rflx, sflx, afus, afds, afuir, afdir, albedo, rr, ice, nice, grp)
      integer, intent (in) :: n1, n2, n3
      real, intent (in)    :: alat, time, sknt, sfc_albedo
      real, dimension (n1), intent (in)                 :: pi0, pi1, zm, zt
      real, dimension (n1,n2,n3), intent (in)           :: dn, pip, tk, rv, rc, nc
      real, optional, dimension (n1,n2,n3), intent (in) :: rr, ice, nice, grp
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
            pres(1:n1) = p00*( exner(1:n1) )**cpr*0.01 ! hPa

            ! Calculates pressure grid, pp, for the momentum grid zm. Also, depending on
            ! the setup, interpolates temperature and humidity for the linear transition
            ! region between soundings and the LES domain.
            call set_pres_grid(n1,pres,zm,zt,tk(n1,i,j),rv(n1,i,j))

            ! LES and background pressure levels must be separate (LES pressure levels will change)
            !      Tomi Raatikainen 17.6.2016 & 21.12.2016
            if ( pp(nv-n1+2) < pp(nv-n1+1)+1.0 ) THEN
                ! Simple solution to the problem: remove sounding level nv-n1+1, which means that LES data is written over that
                nv1=nv1-1       ! This should do it (start from the previos location)
                nv=nv-1
                pp(nv-n1+2) = pres(n1) - 0.5*(pres(n1-1)-pres(n1)) ! Update

                WRITE(*,*) 'Warning (radiation/d4stream): overlapping pressure levels observed (i,j)!',i,j
            endif

            kk = nv1 ! Surface pressure
            do k=2,n1
               ! Scalars at mid points
               kk = kk-1
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

               ! Rain
               if (present(rr)) then
                  prwc(kk) = 1000.*dn(k,i,j)*rr(k,i,j)
               end if

               ! Ice
               if (present(ice)) then
                  if ((ice(k,i,j).gt.0.).and.(nice(k,i,j).gt.0.)) then
                     piwc(kk) = 1000.*dn(k,i,j)*ice(k,i,j)
                     pde(kk)  = 1.e6*(piwc(kk)/(1000.*pri*nice(k,i,j)*dn(k,i,j)))**(1./3.)
                     pde(kk)=min(max(pde(kk),20.),120.)
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

            if (present(ice).and.present(rr).and.present(grp)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, useMcICA, plwc=plwc, pre=pre, piwc=piwc, pde=pde, prwc=prwc, pgwc=pgwc)
            ELSEif (present(ice).and.present(grp)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, useMcICA, plwc=plwc, pre=pre, piwc=piwc, pde=pde, pgwc=pgwc)
            ELSEif (present(ice).and.present(rr)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, useMcICA, plwc=plwc, pre=pre, piwc=piwc, pde=pde, prwc=prwc)
            ELSEif (present(ice)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, useMcICA, plwc=plwc, pre=pre, piwc=piwc, pde=pde)
            ELSEif (present(rr)) then
                call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                     fds, fus, fdir, fuir, useMcICA, plwc=plwc, pre=pre, prwc=prwc)
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

            do k=2,n1
               xfact  = 1./((zm(k)-zm(k-1))*cp*dn(k,i,j)*exner(k))
               tt(k,i,j) = tt(k,i,j) - (rflx(k,i,j) - rflx(k-1,i,j))*xfact
            end do

            IF (calc_od) THEN
                ! Calculate optical depths
                !
                ! Allocate outputs
                IF (.NOT.ALLOCATED(tau_gas)) THEN
                    ALLOCATE(tau_gas(n2,n3), tau_cloud(n2,n3), tau_rain(n2,n3), &
                        tau_ice(n2,n3), tau_grp(n2,n3))
                    tau_gas=0.; tau_cloud=0.; tau_rain=0.; tau_ice=0.; tau_grp=0.
                ENDIF
                !
                IF (.NOT. present(rr)) prwc = 0.0
                IF (.NOT. present(ice)) THEN
                    piwc = 0.0
                    pde = 0.0
                ENDIF
                IF (.NOT. present(grp)) pgwc = 0.0
                !
                ! Calculations
                CALL rad_tau (pp, pt, ph, plwc, pre, prwc, piwc, pde, pgwc, &
                    tau_gas(i,j), tau_cloud(i,j), tau_rain(i,j), tau_ice(i,j), tau_grp(i,j))
             ENDIF

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
    allocate (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),prwc(nv),piwc(nv),pgwc(nv)) ! Cell centers

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
  ! This version of setup is modified for LES simulations
  ! - Any profile can be provided
  ! - Using pressure as interface and interpolating scalars for mid points
  ! - Ozone is interpolated for all levels and kept constant
  ! - Linear transition region can be added between the LES domain and the
  !   background sounding below pressure radsounding_pmax
  ! - Transtion region pressure and scalars (except ozone) are updated for each column
  !
  subroutine setup_les_regrid(background,n1,nv1,nv,zp,tt,qt,zm,zt)
    use mpi_interface, only: myid
    implicit none

    character (len=*), intent (in) :: background
    integer, intent (in) :: n1
    integer, intent (out):: nv1,nv
    real, intent (in)    :: zp(n1),zm(n1),zt(n1),tt(n1),qt(n1)
    real, allocatable  :: sp(:), st(:), sh(:), so(:)
    integer :: k, kk, ns
    real    :: ptop, dp, Tsurf, pmid

    ! Read backgroud profiles (p [hPa], T [K], q [kg/kg], O3 [-] and an unused concentration)
    open ( unit = 08, file = background, status = 'old' )
    ! Some files have surface temperature and dimensions in the first line, but these may not be correct
    read (08,*) Tsurf, ns
    ! ... so count the number of lines
    ns=0; k=0
    do while(k>=0)
        read (8,*,iostat=k) ptop
        ns=ns+1
    enddo
    ! ns is the total number of lines including the first line
    if (Tsurf>200.) ns=ns-1 ! Total number of sounding values

    ! Allocate data
    allocate (sp(ns), st(ns), sh(ns), so(ns))
    ! Read again
    rewind (8)
    if (Tsurf>200.) read (08,*) Tsurf, k
    do k=1,ns
        read (8,*) sp(k), st(k), sh(k), so(k)
    enddo
    close (08)

    ! The highest LES model pressure - defined for cell interface zm(n1)
    ptop=zp(n1)*(zp(n1-1)/zp(n1))**( (zm(n1)-zt(n1))/(zt(n1-1)-zt(n1)) )
    !
    ! The first sounding level must be somewhat higher (lower pressure) than the last LES pressure level
    !  - Must avoid overlapping pressure levels when LES pressures change with time
    !  - 10 hPa should be large enough for most cases (LES has typically high pressure resolution)
    !  - radsounding_pmax>0 [hPa] can be set to activate thicker transition layer
    dp=MAX(10.0,zp(n1-1)-zp(n1))
    !
    k=1 ! Note: increasing pressure and decreasing altitude
    do while (sp(k) < ptop-dp .AND. (radsounding_pmax<0.0 .OR. sp(k) < radsounding_pmax))
        k=k+1
        if (k==ns+1) exit
    enddo
    nb=k-1 ! Background soundings from 1 to k (nb points)
    !
    ! Add layers between LES and soundings: the transition regime
    nt=max(0,FLOOR( (ptop-sp(nb))/dp )-1)

    ! Include complete LES grid (n1 points), soundings (nb points) and the transition regime (nt points)
    nv1 = n1+nb+nt  ! Total number of pressure points (cell interfaces)
    nv = nv1-1   ! Total number of scalars (cell centers)

    ! allocate the arrays for the sounding data to be used in the radiation profile
    allocate (pp(nv1),fds(nv1),fus(nv1),fdir(nv1),fuir(nv1)) ! Cell interfaces
    allocate (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),prwc(nv),piwc(nv),pgwc(nv)) ! Cell centers

    ! Levels above the transition region are fixed, so fill these with sounding data
    do k=1,nb
        ! Pressure as is (interface)
        pp(k) = sp(k)
        ! Scalars interpolated for mid points (only up to nb-1)
        ! Linear T vs log(P)
        pt(k) = 0.5*(st(k)+st(k+1))
        ! Linear log(c) vs log(P)
        ph(k) = sqrt(sh(k)*sh(k+1))
        po(k) = sqrt(so(k)*so(k+1))
    enddo

    ! Set the remaining pressure gid for the given initial pressure profile.
    ! Also sets water vapor and temperature.
    call set_pres_grid(n1,zp,zm,zt,tt(n1),qt(n1))

    ! Transition region: interpolate ozone for the mid point(s) based on the background sounding
    do k=nb,nb+nt
        pmid=sqrt(pp(k)*pp(k+1))
        po(k) = interp_press_log(sp,so,ns,pmid)
    enddo

    ! The LES domain
    k=nb+nt
    do kk=n1,2,-1 ! Start from the top
        ! Use the initial grid cell pressures
        k=k+1
        pmid=zp(kk)
        ! Ozone is based on the background sounding
        po(k) = interp_press_log(sp,so,ns,pmid)
        ! Temperature and moisture come from the LES - here for initial information
        pt(k) = tt(kk)
        ph(k) = qt(kk)
    enddo

    ! Clean
    deallocate(sp, st, sh, so)

    ! Print information about sounding
    IF (myid==0) THEN
        WRITE(*,'(/,/,20A)') '-------------------------------------------------'
        WRITE(*,'(/2x,20A,/)') 'Reading Background Sounding: '//trim(background)
        WRITE(*,'(/2x,A7,5A10)') 'Source','p1 (hPa)','p2 (hPa)','T (K)','q (g/kg)','O3 (ppm)'

        DO k=1,nv
            IF (k<nb) THEN
                ! Data from backgound soundings
                WRITE(*,'(2x,A7,3F10.2,F10.5,F10.3)') 'Backgr',pp(k),pp(k+1),pt(k),ph(k)*1.e3,po(k)*1.e6
            ELSEIF (k<=nb+nt) THEN
                ! Interpolation
                WRITE(*,'(2x,A7,3F10.2,F10.5,F10.3)') 'Interp',pp(k),pp(k+1),pt(k),ph(k)*1.e3,po(k)*1.e6
            ELSEIF (k<=nb+nt+5) THEN
                ! LES top + four levels
                WRITE(*,'(2x,A7,3F10.2,F10.5,F10.3)') 'LES',pp(k),pp(k+1),pt(k),ph(k)*1.e3,po(k)*1.e6
            ENDIF
        ENDDO
        WRITE(*,*) ' '
    ENDIF
    !
    contains
      real function interp_press_log(pp,dat,n,pval)
        integer, intent (in) :: n
        real, intent (in) :: pp(n),dat(n),pval
        integer :: i
        ! pp is decreasing
        i=COUNT(pval>pp)
        if (i==0) then
            interp_press_log = dat(1)
        elseif (i==n) then
            interp_press_log = dat(n)
        else
            ! Interpolate between pp(i) and pp(i+1)
            interp_press_log = dat(i)*exp(log(dat(i)/dat(i+1))/log(pp(i)/pp(i+1))*log(pval/pp(i)))
        endif
      end function interp_press_log
  end subroutine setup_les_regrid
  !
  ! Updates radiation pressure grid for the LES and transition regions, and optionally sets
  ! water vapor mixing ratio and temperature values for the transition region.
  subroutine set_pres_grid(n1,pres,zm,zt,ttop,qtop)
    implicit none
    integer, intent (in) :: n1
    real, intent (in) :: pres(n1),zm(n1),zt(n1),ttop,qtop ! LES inputs
    integer :: kk, k
    real :: dp, ptop, pmid, psnd

    ! 1) LES pressure grid
    kk = nv1 ! Surface
    do k=2,n1
        ! Interpolate pressure (hPa) for level zm(k-1), which is between zt(k-1) and zt(k)
        pp(kk) = pres(k)*(pres(k-1)/pres(k))**( (zm(k-1)-zt(k))/(zt(k-1)-zt(k)) )
        kk = kk-1
    enddo
    ! Extrapolate pressure (hPa) for level zm(n1) based on pressures at zt(n1-1) and zt(n1).
    ! For radiation, this means pressure at the top of the LES domain, pp(nv-n1+2).
    ! It is assumed that pressure is hydrostatic, i.e., exponentially decreasing with altitude.
    pp(kk) = pres(n1)*(pres(n1-1)/pres(n1))**( (zm(n1)-zt(n1))/(zt(n1-1)-zt(n1)) )

    ! 2) Transition region between sounding base pp(nb) and LES top pp(kk)
    ptop=pres(n1) ! Pressure for ttop and qtop
    psnd=sqrt(pp(nb)*pp(nb-1)) ! The lowest sounding midpoint for pt(nb-1) and ph(nb-1)
    if (nt>0) then
        ! Interpolated logarithmic pressure levels between pp(nb) and pp(kk)
        dp=log(pp(kk)/pp(nb))/REAL(nt+1)
        do k=nb+1,nb+nt
            pp(k)=pp(k-1)*exp(dp)
        enddo
        !
        ! Interpolate temperature and moisture for the mid points based on their values at the
        ! LES model top (ptop, ttop, qtop) and the lowest sounding level (psnd, pt(nb-1), ph(nb-1)).
        do k=nb,nb+nt ! include nb
            pmid=sqrt(pp(k)*pp(k+1))
            pt(k)=ttop+(pt(nb-1)-ttop)/log(psnd/ptop)*log(pmid/ptop)
            ph(k)=qtop*exp( log(ph(nb-1)/qtop)/log(psnd/ptop)*log(pmid/ptop) )
        enddo
    elseif (nb>0) then
        ! Single mid point between pp(nb) and pp(kk)=pp(nb+1)
        pmid=sqrt(pp(nb)*pp(nb+1))
        ! Interpolate
        pt(nb)=ttop+(pt(nb-1)-ttop)/log(psnd/ptop)*log(pmid/ptop)
        ph(nb)=qtop*exp( log(ph(nb-1)/qtop)/log(psnd/ptop)*log(pmid/ptop) )
    endif

  end subroutine set_pres_grid
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
