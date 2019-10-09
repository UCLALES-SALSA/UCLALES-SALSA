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
MODULE radiation

  ! 20151022: Made some modifications for parameters going to *rad* to
  ! avoid errors with SALSA. 
  ! Juha Tonttila, FMI


  USE defs, ONLY       : cp, rcp, cpr, rowt, roice, p00, pi, nv1, nv, SolarConstant, R
  USE fuliou, ONLY     : rad, set_random_offset, minSolarZenithCosForVis
  USE mo_submctl, ONLY : nbins,spec,pi6
  IMPLICIT NONE

  CHARACTER (len=50) :: background = 'datafiles/dsrt.lay'
  LOGICAL :: McICA = .TRUE.

  LOGICAL :: ConstPress = .FALSE.
  REAL, ALLOCATABLE, SAVE :: exner(:), pres(:)

  LOGICAL, SAVE     :: first_time = .TRUE.
  REAL, ALLOCATABLE, SAVE ::  pp(:), pt(:), ph(:), po(:), pre(:), pde(:), &
                              plwc(:), piwc(:), prwc(:), pgwc(:), fds(:), fus(:), fdir(:), fuir(:), &
                              maerobin(:,:), naerobin(:,:)
  INTEGER :: k,i,j, npts
  REAL    :: ee, u0, day, zz !time, alat, Juha: Time and alat given as argument already! Potential bug, hopefully harmless.

  CONTAINS

    SUBROUTINE d4stream(n1, n2, n3, nspec, alat, time, sknt, sfc_albedo, dn0, pi0, pi1, dzm, &
                        pip, tk, rv, rc, nc, tt, rflx, sflx, afus, afds, afuir, afdir, &
                        albedo, rr, ice, nice, grp, radsounding, useMcICA, ConstPrs, maerop, naerop)
      USE mpi_interface, ONLY : myid, pecount
      INTEGER, INTENT (in) :: n1, n2, n3, nspec
      REAL, INTENT (in)    :: alat, time, sknt, sfc_albedo  ! Time in decimal days
      REAL, DIMENSION (n1), INTENT (in)                 :: dn0, pi0, pi1, dzm
      REAL, DIMENSION (n1,n2,n3), INTENT (in)           :: pip, tk, rv, rc, nc
      REAL, OPTIONAL, DIMENSION (n1,n2,n3), INTENT (in) :: rr, ice, nice, grp
      REAL, OPTIONAL, INTENT(in)                        :: maerop(n1,n2,n3,nspec*nbins),   & 
                                                           naerop(n1,n2,n3,nbins)
      REAL, DIMENSION (n1,n2,n3), INTENT (inout)        :: tt, rflx, sflx, afus, afds, afuir, afdir
      CHARACTER(len=50), OPTIONAL, INTENT(in)           :: radsounding
      LOGICAL, OPTIONAL, INTENT(in)                     :: useMcICA, ConstPrs
      !! NEED TO FIND A BETTER WAY WITH THESE INDICES BECAUSE THEY'RE ALL OVER THE PLACE NOW: 

      REAL, INTENT (out)                                :: albedo(n2,n3)

      INTEGER :: kk
      REAL    :: xfact, prw, pri, p0(n1)
      REAL    :: eps                        !Added to prevent division by zero, not entirely sure what the best value is - Aleksi

      eps = 1.e-25

      IF (present(radsounding)) background = radsounding
      IF (present(useMcICA)) McICA = useMcICA
      IF (present(ConstPrs)) ConstPress = ConstPrs

      IF (first_time) THEN
         ! Possible to use constant LES pressure levels (fixed during the first call)
         ALLOCATE(exner(n1), pres(n1))
         exner(1:n1) = (pi0(1:n1)+pi1(1:n1))/cp
         pres(1:n1) = p00*( exner(1:n1) )**cpr
         ! Works well with the LES model
         p0(1:n1) = 0.01*pres(1:n1)
         CALL setup_les(background,n1,nv1,nv,p0,tk(n1,3,3),rv(n1,3,3))
         first_time = .FALSE.
         IF (allocated(pre))   pre(:) = 0.
         IF (allocated(pde))   pde(:) = 0.
         IF (allocated(piwc)) piwc(:) = 0.
         IF (allocated(prwc)) prwc(:) = 0.
         IF (allocated(plwc)) plwc(:) = 0.
         IF (allocated(pgwc)) pgwc(:) = 0.
      END IF

      IF (allocated(maerobin)) maerobin(:,:) = 0.
      IF (allocated(naerobin)) naerobin(:,:) = 0.

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
      IF (McICA .AND. myid > 0) THEN

        IF (u0 > minSolarZenithCosForVis) THEN
            kk = myid*(n3-4)*(n2-4)*2
        ELSE
            kk = myid*(n3-4)*(n2-4) ! Without solar wavelengths
        END IF

        CALL set_random_offset( kk )
      END IF
      !
      ! call the radiation
      !
      prw = (4./3.)*pi*rowt
      pri = (3.*sqrt(3.)/8.)*roice
      DO j = 3, n3-2
         DO i = 3, n2-2
            IF (.NOT. ConstPress) THEN
                ! Grid cell pressures in the LES model (Pa)
                exner(1:n1) = (pi0(1:n1)+pi1(1:n1)+pip(1:n1,i,j))/cp
                pres(1:n1) = p00*( exner(1:n1) )**cpr
            END IF

            ! LES and background pressure levels must be separate (LES pressure levels will change)
            !      Tomi Raatikainen 17.6.2016 & 21.12.2016
            pp(nv-n1+2) = pres(n1)/100. - 0.5*(pres(n1-1)-pres(n1)) / 100.
            IF ( pp(nv-n1+2) < pp(nv-n1+1)+1.0 ) THEN
                ! Simple solution to the problem: remove sounding level nv-n1+1, which means that LES data is written over that
                nv1 = nv1-1       ! This should do it (start from the previos location)
                nv = nv-1
                pp(nv-n1+2) = pres(n1)/100. - 0.5*(pres(n1-1)-pres(n1)) / 100. ! Update

                WRITE(*,*) 'Warning (radiation/d4stream): overlapping pressure levels observed (i,j)!',i,j
            END IF

            DO k = 2, n1
               kk = nv-(k-2)
               pp(kk+1) = 0.5*(pres(k-1)+pres(k)) / 100.
               pt(kk) = tk(k,i,j)
               ph(kk) = rv(k,i,j)

               ! Cloud water
               IF ((rc(k,i,j) > 0.) .AND. (nc(k,i,j) > eps)) THEN
                  plwc(kk) = 1000.*dn0(k)*rc(k,i,j)
                  pre(kk)  = 1.e6*(plwc(kk)/(1000.*prw*nc(k,i,j)*dn0(k)))**((1./3.))
                  pre(kk) = min(max(pre(kk),4.18),31.23)
               ELSE
                  pre(kk) = 0.
                  plwc(kk) = 0.
               END IF

               ! Precipitation (not used at the moment)
               IF (present(rr)) THEN
                  prwc(kk) = 1000.*dn0(k)*rr(k,i,j)
               END IF

               ! Ice
               IF (present(ice)) THEN
                  IF ((ice(k,i,j) > 0.) .AND. (nice(k,i,j) > eps)) THEN
                     piwc(kk) = 1000.*dn0(k)*ice(k,i,j)
                     pde(kk)  = 1.e6*(piwc(kk)/(1000.*pri*nice(k,i,j)*dn0(k)))**((1./3.))
                     pde(kk) = min(max(pde(kk),20.),180.)
                  ELSE
                     piwc(kk) = 0.0
                     pde(kk)  = 0.0
                  END IF
               END IF

               ! Graupel
               IF (present(grp)) THEN
                  pgwc(kk) = 1000.*dn0(k)*grp(k,i,j)
               END IF

               ! Aerosol ! convert from /kg to /m3
               IF ( PRESENT(maerop) .AND. PRESENT(naerop) ) THEN
                  maerobin(kk,:) = maerop(k,i,j,:)*pres(k)/(R*tk(k,i,j))
                  naerobin(kk,:) = naerop(k,i,j,:)*pres(k)/(R*tk(k,i,j))
               END IF

            END DO

            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! REMEMBER TO ADD AEROSOL TO OTHER OPTIONS AS WELL 
            ! MAYBE BETTER TO DO SOME REFACTORING FIRST...
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !IF (PRESENT(ice).AND.PRESENT(rr).AND.PRESENT(grp)) THEN
            !   CALL rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
            !             fds, fus, fdir, fuir, McICA, nspec, plwc=plwc, pre=pre, &
            !             piwc=piwc, pde=pde, prwc=prwc, pgwc=pgwc)
            !ELSE IF (PRESENT(ice).AND.PRESENT(grp)) THEN
            !   CALL rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
            !             fds, fus, fdir, fuir, McICA, nspec, plwc=plwc, pre=pre, &
            !             piwc=piwc, pde=pde, pgwc=pgwc)
            !ELSE IF (PRESENT(ice).AND.PRESENT(rr)) THEN
            !   CALL rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
            !             fds, fus, fdir, fuir, McICA, nspec, plwc=plwc, pre=pre, &
            !             piwc=piwc, pde=pde, prwc=prwc)
            IF (PRESENT(ice)) THEN
               CALL rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                         fds, fus, fdir, fuir, McICA, nspec, plwc=plwc, pre=pre, &
                         piwc=piwc, pde=pde)
            !ELSE IF (PRESENT(rr)) THEN
            !   CALL rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
            !             fds, fus, fdir, fuir, McICA, nspec, plwc=plwc, pre=pre, &
            !             prwc=prwc)
            ELSE IF (PRESENT(maerop) .AND. PRESENT(naerop)) THEN
               CALL rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                         fds, fus, fdir, fuir, McICA, nspec, plwc=plwc, pre=pre, &
                         maerobin=maerobin, naerobin=naerobin )
            ELSE
               CALL rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                         fds, fus, fdir, fuir, McICA, nspec, plwc=plwc, pre=pre)
            END IF

            DO k = 1, n1
               kk = nv1 - (k-1)
               afus(k,i,j)  = fus(kk)
               afds(k,i,j)  = fds(kk)
               afuir(k,i,j) = fuir(kk)
               afdir(k,i,j) = fdir(kk)
               sflx(k,i,j)  = fus(kk)  - fds(kk)
               rflx(k,i,j)  = sflx(k,i,j) + fuir(kk) - fdir(kk)
            END DO

            IF (u0 > minSolarZenithCosForVis) THEN
               albedo(i,j) = fus(1)/fds(1)
            ELSE
               albedo(i,j) = -999.
            END IF

            DO k = 2, n1-3
               xfact = dzm(k)/(cp*dn0(k)*exner(k))
               tt(k,i,j) = tt(k,i,j) - (rflx(k,i,j) - rflx(k-1,i,j))*xfact
            END DO

         END DO
      END DO

      ! Call random numbers to add offset for PUs after myid
      IF (McICA .AND. myid < pecount-1) THEN
        IF (u0 > minSolarZenithCosForVis) THEN
            kk = (pecount-1-myid)*(n3-4)*(n2-4)*2
        ELSE
            kk = (pecount-1-myid)*(n3-4)*(n2-4) ! Without solar wavelengths
        END IF
        CALL set_random_offset( kk )
      END IF
    END SUBROUTINE d4stream

  ! ---------------------------------------------------------------------------
  ! sets up the input data to extend through an atmosphere of appreiciable
  ! depth using a background souding specified as a parameter, match this to
  ! the original sounding using p0 as this does not depend on time and thus
  ! allows us to recompute the same background matching after a history start
  !

  ! ---------------------------------------------------------------------------
  ! This version of setup is modified for LES simulations
  !    - Always using background soundings and LES data
  !    - Separate LES and background profiles
  !    - Weighting of the backgound towards smooth transition at the interface
  !
  ! Tomi Raatikainen, 22.12.2016
  !
  SUBROUTINE setup_les(background,n1,nv1,nv,zp,ttop,qtop)
    USE mpi_interface, ONLY: myid
    IMPLICIT NONE

    CHARACTER (len=19), INTENT (in) :: background
    INTEGER, INTENT (in)  :: n1
    INTEGER, INTENT (out) :: nv1,nv
    REAL, INTENT (in)     :: zp(n1), ttop, qtop

    REAL, ALLOCATABLE  :: sp(:), st(:), sh(:), so(:), sl(:)

    INTEGER :: k, ind, ns, nb, nt
    REAL    :: ptop, dp, Tsurf
    INTEGER :: nspec
    
    nspec = spec%getNSpec()
    nspec = MAX(1,nspec) ! If SALSA is not used, nspec would be 0 and that might cause an error in the allocation 
                         ! (which is done regardles if SALSA is used, to keep it simple)

    ! Read backgroud profiles (p [hPa], T [K], q [kg/kg], O3 [-] and an unused concentration)
    OPEN( unit = 08, file = background, status = 'old' )
    READ(08,*) Tsurf, ns
    ALLOCATE ( sp(ns), st(ns), sh(ns), so(ns), sl(ns))
    DO k = 1, ns
       READ( 08, *) sp(k), st(k), sh(k), so(k), sl(k)
    END DO
    CLOSE(08)

    ! Merge background and LES pressure levels (p_LES > p_background)
    !
    ! The highest LES model pressure - defined for cell interface
    ptop = zp(n1)-0.5*(zp(n1-1)-zp(n1))
    !
    ! The first sounding level must be somewhat higher (lower pressure) then the last LES pressure level
    !  - Must avoid overlapping pressure levels when LES pressures change with time
    !  - 10 hPa should be large enough for most cases (LES has typically high pressure resolution)
    dp = MAX(10.0,zp(n1-1)-zp(n1))
    !
    k = 1 ! Note: increasing pressure and decreasing altitude
    DO WHILE (sp(k) < ptop-dp )
        k = k+1
        IF (k == ns+1) EXIT
    END DO
    nb = k-1 ! Background soundings from 1 to k (nb points)

    ! Add layers between LES and soundings (mind the gap!)
    nt = 0 ! Transition regime
    IF (nb > 0) nt = FLOOR( (ptop-sp(nb))/dp )-1

    ! Include complete LES grid (n1 points), soundings (nb points) and the transition regime (nt points)
    nv1 = n1+nb+nt  ! Total number of pressure points (cell interfaces)
    nv = nv1-1   ! Total number of scalars (cell centers)

    ! Allocate the arrays for the sounding data to be used in the radiation
    ! profile and then fill them first with the sounding data, by afill, then
    ! by interpolating the background profile at pressures less than the
    ! pressure at the top of the sounding
    !
    ! Note: LES pressure levels are increasing, but radiation calculations
    ! expect decreasing pressure grid (from TOA to surface)
    !
    ALLOCATE (pp(nv1),fds(nv1),fus(nv1),fdir(nv1),fuir(nv1)) ! Cell interfaces
    ALLOCATE (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),prwc(nv),piwc(nv),pgwc(nv)) ! Cell centers
    ALLOCATE(maerobin(nv,nspec*nbins),naerobin(nv,nbins))
    
    po = 0.
    IF (nb > 0) THEN
        ! Levels above LES domain: copy background soundings
        pp(1:nb) = sp(1:nb)
        pt(1:nb) = st(1:nb)
        ph(1:nb) = sh(1:nb)
        po(1:nb) = so(1:nb)
    END IF

    IF (nt > 0) THEN
        ! Interpolated levels between pp(nb) and ptop
        DO k = 1, nt
            pp(nb+k) = pp(nb+k-1)-(pp(nb)-ptop)/REAL(nt+1)
            ! Interpolate between LES model top and the lowest sounding
            !   y=y0+(y1-y0)(p1-p0)*(p-p0)
            ind = getindex(sp,ns,pp(nb+k)) ! Typically ind=nb
            pt(nb+k) = ttop+(st(ind)-ttop)/(sp(ind)-ptop)*(pp(nb+k)-ptop)
            ph(nb+k) = qtop+(sh(ind)-qtop)/(sp(ind)-ptop)*(pp(nb+k)-ptop)
            ! Interpolate ozone using background soundings
            po(nb+k) = so(ind+1)+(so(ind)-so(ind+1))/(sp(ind)-sp(ind+1))*(pp(nb+k)-sp(ind+1))
        END DO
    END IF

    ! Within the LES domain
    !  a) pressure, temperature and humidity will be obtained from the LES
    !  b) set the ozone to a constant value
    IF (nb+nt > 0) po(nb+nt+1:nv) = po(nb+nt)

    ! Print information about sounding
    IF (myid == 0) THEN
        WRITE(*,'(/,/,20A)') '-------------------------------------------------'
        WRITE(*,'(/2x,20A,/)') 'Reading Background Sounding: '//trim(background)
        WRITE(*,'(/2x,A7,4A10)') 'Source','p (hPa)','T (K)','q (g/kg)','O3 (ppm)'

        ! Calculate LES pressure levels ([zp]=hPa) and T and q for cell centers
        pp(nv-n1+2) = zp(n1)/100. - 0.5*(zp(n1-1)-zp(n1))
        DO k = 2, n1
             pp(nv-(k-2)+1) = 0.5*(zp(k-1)+zp(k))
        END DO

        DO k = 1, nv
            IF (k <= nb .AND. nb > 0) THEN
                ! Data from backgound soundings
                WRITE(*,'(2x,A7,2F10.2,2F10.3)') 'Backgr',pp(k),pt(k),ph(k)*1.e3,po(k)*1.e6
            ELSE IF (k <= nb+nt .AND. nt > 0) THEN
                ! Interpolation
                WRITE(*,'(2x,A7,2F10.2,2F10.3)') 'Interp',pp(k),pt(k),ph(k)*1.e3,po(k)*1.e6
            ELSE IF (k == nb+nt+1) THEN
                ! The highest LES layer - other layers not yet available
                WRITE(*,'(2x,A7,2F10.2,2F10.3)') 'LES',ptop,ttop,qtop*1.e3,po(k)*1.e6
                !WRITE(*,'(2x,A7,3A10,F10.3)') 'LES','*','*','*',po(k)*1.e6
                WRITE(*,'(2x,A7,4A10)') 'LES','...','...','...','...'
            END IF
        END DO
        WRITE(*,*) ' '
    END IF

  END SUBROUTINE setup_les
  ! ---------------------------------------------------------------------------
  !  locate the index closest to a value
  !
  INTEGER FUNCTION getindex(x,n,xval)

    INTEGER, INTENT (in) :: n
    REAL,    INTENT (in) :: x(n),xval

    INTEGER :: ia, ib

    ia = 1
    ib = n
    IF (xval < x(1)) THEN
       getindex = 1
    ELSE IF (xval > x(n)) THEN
       getindex = n-1
    ELSE
       getindex = (ia+ib)/2
       DO WHILE (getindex /= ia .OR. ib /= getindex+1)
          getindex = (ia+ib)/2
          IF ((xval-x(getindex)) >= 0.0) THEN
             ia = getindex
          ELSE
             ib = getindex
          END IF
       END DO
    END IF

  END FUNCTION getindex

  ! ---------------------------------------------------------------------------
  ! linear interpolation between two points,
  !
  REAL FUNCTION intrpl(x1,y1,x2,y2,x)

    REAL, INTENT (in) :: x1,y1,x2,y2,x

    REAL :: slope

    slope  = (y2-y1)/(x2 - x1 + epsilon(1.))
    intrpl = y1+slope*(x-x1)

  END FUNCTION intrpl

  ! ---------------------------------------------------------------------------
  ! Return the cosine of the solar zenith angle give the decimal day and
  ! the latitude
  !
  REAL FUNCTION zenith(alat,time)

    REAL, INTENT (in) :: alat, time

    REAL :: lamda, d, sig, del, h, day

    day    = floor(time)
    lamda  = alat*pi/180.
    d      = 2.*pi*int(time)/365.
    sig    = d + pi/180.*(279.9340 + 1.914827*sin(d) - 0.7952*cos(d) &
         &                      + 0.019938*sin(2.*d) - 0.00162*cos(2.*d))
    del    = asin(sin(23.4439*pi/180.)*sin(sig))
    h      = 2.*pi*((time-day)-0.5)
    zenith = sin(lamda)*sin(del) + cos(lamda)*cos(del)*cos(h)

  END FUNCTION zenith

END MODULE radiation
