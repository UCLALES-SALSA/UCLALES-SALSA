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
! Copyright 1999, 2001, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE grid
  
  USE classFieldArray, ONLY : FieldArray
  USE mo_aux_state
  USE mo_submctl, ONLY : spec, nbins, ncld, nprc, nice
  
  IMPLICIT NONE

  SAVE
  
  CHARACTER(len=10), PARAMETER :: global_name = "grid"
  !
  INTEGER :: nxp = 132           ! number of x points
  INTEGER :: nyp = 132           ! number of y points
  INTEGER :: nzp = 105           ! number of z points
    
  LOGICAL :: nxpart = .TRUE.     ! number of processors in x
  
  REAL    :: deltax = 35.        ! dx for basic grid
  REAL    :: deltay = 35.        ! dy for basic grid
  REAL    :: deltaz = 17.5       ! dz for basic grid
  REAL    :: dzrat  = 1.02       ! grid stretching ratio
  REAL    :: dzmax  = 1200.      ! height to start grid-stretching
  REAL    :: dtlong = 10.0       ! long timestep
  REAL    :: th00   = 288.       ! basic state temperature
  
  REAL    :: CCN = 150.e6
  REAL    :: cntlat =  31.5      ! Latitude for radiation
  
  LOGICAL :: lbinanl = .FALSE.   ! Whether to write binned data to analysis files (takes a lot of space + mainly used for debugging)
  LOGICAL :: lbinprof = .TRUE.   ! The same for profile statistics
  LOGICAL :: lnudging = .FALSE.  ! Master switch for nudging scheme
  LOGICAL :: lemission = .FALSE. ! Master switch for aerosol emission
  
  INTEGER :: iradtyp
  INTEGER :: igrdtyp = 1         ! vertical grid type
  INTEGER :: isgstyp = 1         ! sgs model type
  INTEGER :: level   = 0         ! thermodynamic level
  INTEGER :: naddsc  = 0         ! number of additional scalars;
  INTEGER :: nsalsa  = 0         ! Number of tracers for SALSA
  INTEGER :: nfpt = 10           ! number of rayleigh friction points
  REAL    :: distim = 300.0      ! dissipation timescale
  
  REAL    :: sst = 283.   ! Surface temperature      added by Zubair Maalick
  REAL    :: W1  = 0.9   ! Water content
  REAL    :: W2  = 0.9
  REAL    :: W3  = 0.9
  
  LOGICAL            :: lsalsabbins = .FALSE.  ! This is brought here temporarily from stat.f90 
  CHARACTER (len=150) :: expnme = 'Default' ! Experiment name
  CHARACTER (len=150) :: filprf = 'x'       ! File Prefix
  CHARACTER (len=7)  :: runtype = 'INITIAL'! Run Type SELECTion

  ! Output file list; given here instead of mo_output.f90 to avoid cyclic dependencies
  CHARACTER(len=10) :: varlist_main(100),varlist_ps(100), varlist_ts(100)  ! 

  CHARACTER (len=7),  PRIVATE :: v_snm = 'sxx    '

  ! Grid definitions
  ! -----------------------------------------------------
  INTEGER           :: nz, nxyzp, nxyp
  REAL              :: dxi, dyi, dtl, dtlv, dtlt, umean, vmean, psrf

  REAL, ALLOCATABLE :: spng_wfct(:), spng_tfct(:)
  
  ! Some zero arrays ice with level < 5
  REAL, ALLOCATABLE, TARGET :: tmp_icep(:,:,:,:), tmp_icet(:,:,:,:)
  !
  ! wsave variables used in fft in x and y directons
  !
  REAL, ALLOCATABLE :: wsavex(:), wsavey(:)
  !
  ! For looping through scalar arrays
  REAL, POINTER :: a_sp(:,:,:),a_st(:,:,:) 
  
  ! Some stuff for cloud base activation (not recommended)
  REAL, ALLOCATABLE :: a_vactd(:,:,:,:), a_nactd(:,:,:,:)
  
  !---------------------------------------------------------------------------

  ! Main containers for the prognostic scalar values
  ! ----------------------------------------------------------------
  REAL, ALLOCATABLE, TARGET :: a_sclrp(:,:,:,:),a_sclrt(:,:,:,:)
  ! ----------------------------------------------------------------

  ! Other field arrays are initialized and stored in mo_field_state. The ones below are needed here
  ! but sine they don't contain any variable associated procedures, there is no risk for cyclic dependencies.
  TYPE(FieldArray) :: Axes              ! Contains the grid and size distribution axis vectors
  TYPE(FieldArray) :: outAxes           ! Subset from above assigned for output
  TYPE(FieldArray) :: outAxesPS         ! Subset from Axes assigned for ps-file output
  TYPE(FieldArray) :: BasicState        ! Basic state profiles (no particles)
  TYPE(FieldArray) :: outBasicState     ! Subset from above assigned for output
  
  ! Juha:
  ! Diagnostic variables needed to track mass conservation (of water).
  ! These are reset at every statistical output timestep (better use quite long statistical periods...)
  REAL :: mc_Mtot           ! Initial mass of water normalized by domain volume (== domain mean concentration)
  REAL :: mc_Matm           ! Atmospheric water content (instantaneous, domain mean concentration)
  REAL :: mc_Mevap          ! Evaporated water content, accumulated, Normalized by *domain surface area*/*domain volume*
  REAL :: mc_Mprec          ! Precipitated water content, - '' -
  REAL :: mc_Vdom           ! Domain volume
  REAL :: mc_Adom           ! Domain surface area
  REAL :: mc_ApVdom         ! Volume/Area
  
  !
  INTEGER :: nscl = 1
  INTEGER :: memsize
  !INTEGER, SAVE :: ncid0,ncid_s
  !
  
CONTAINS
   !
   !----------------------------------------------------------------------
   ! Subroutine define_vars
   !
   ! Modified for level 4
   ! Juha Tonttila, FMI, 2014.
   !
   SUBROUTINE define_vars

      USE mpi_interface, ONLY : myid

      CHARACTER(len=20), PARAMETER :: name = "define_vars"
      INTEGER :: zz
      INTEGER :: nc
      INTEGER :: st_salsa,en_salsa ! start and end indices for SALSA tracers
      
      ! Instanciate the field arrays
      BasicState = FieldArray()
      
      nc = 0
      ! Juha: Number of prognostic tracers for SALSA
      !       Aerosol bins + Cloud bins + gas compound tracers
      IF (level >= 4) THEN
         nc = spec%getNSpec(type="wet")
         nsalsa = (nc+1)*nbins + (nc+1)*ncld + (nc+1)*nprc + 5
         IF (level == 5) nsalsa = nsalsa + (nc+1+1)*nice ! (nc+1+1)*nice for RIMED ICE 
      END IF

      ! Initial condition vectors
      CALL setInitialProfiles(BasicState,nzp)
      CALL BasicState%getByOutputstatus(outBasicState)

      memsize = 2*nxyzp ! complexarray in pressure solver

!      ! Vector variables
!      ALLOCATE (a_up(nzp,nxp,nyp),a_vp(nzp,nxp,nyp),a_wp(nzp,nxp,nyp))
!      a_up(:,:,:) = 0.
!      a_vp(:,:,:) = 0.
!      a_wp(:,:,:) = 0.

!      ALLOCATE (a_uc(nzp,nxp,nyp),a_vc(nzp,nxp,nyp),a_wc(nzp,nxp,nyp))
!      a_uc(:,:,:) = 0.
!      a_vc(:,:,:) = 0.
!      a_wc(:,:,:) = 0.

!      ALLOCATE (a_ut(nzp,nxp,nyp),a_vt(nzp,nxp,nyp),a_wt(nzp,nxp,nyp))
!      a_ut(:,:,:) = 0.
!      a_vt(:,:,:) = 0.
!      a_wt(:,:,:) = 0.

      ! Juha: Allocate the main scalar arrays
      !-----------------------------------------------------
      IF (level < 4) THEN

         nscl = nscl+naddsc
         IF (level   > 0) nscl = nscl+1
         IF (level   > 2) nscl = nscl+2
         IF (isgstyp > 1) nscl = nscl+1

         ALLOCATE (a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
         a_sclrp(:,:,:,:) = 0.
         a_sclrt(:,:,:,:) = 0.

      ELSE IF (level >= 4) THEN
         
         ALLOCATE ( a_nactd(nzp,nxp,nyp,ncld), a_vactd(nzp,nxp,nyp,nc*ncld)  )
         a_nactd(:,:,:,:) = 0.
         a_vactd(:,:,:,:) = 0.

         ! Total number of prognostic scalars: temp + water vapor + tke(isgstyp>1) + SALSA
         nscl = 2 + nsalsa
         IF (isgstyp > 1) nscl = nscl+1

         ALLOCATE (a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
         a_sclrp(:,:,:,:) = 0.
         a_sclrt(:,:,:,:) = 0.

      END IF ! level
           
   END SUBROUTINE define_vars
   !
   !----------------------------------------------------------------------
   !
   SUBROUTINE define_grid

      USE mpi_interface, ONLY : xoffset, yoffset, wrxid, wryid, nxpg, nypg,   &
                                appl_abort, myid

      CHARACTER(len=20), PARAMETER :: name = "define_grid"
      INTEGER :: i,j,k,kmax,nchby
      REAL    :: dzrfm,dz,zb,dzmin
      REAL    :: zmnvc(-1:nzp+1)
      CHARACTER (len=51) :: &
         fm1 = '(//" ",49("-")/,"   grid DIMENSIONs:"/)            ',      &
         fm2 = '("   nxp-4 = ",i3,", dx, dx = ",f8.1,",",f8.1," m")',      &
         fm3 = '("   nyp-4 = ",i3,", dy, dy = ",f8.1,",",f8.1," m")',      &
         fm4 = '("   nzp   = ",i3,", dz, dz = ",f8.1,",",f8.1," m")',      &
         fm5 = '("   timestep: ",f7.3,"s ")                        ',      &
         fm6 = '("   thermo level: ",i3)                        '


      ! Initialize grid vectors
      Axes = FieldArray()
      CALL setGridSpacings(Axes,lbinanl,level,nzp,nxp,nyp)
      CALL Axes%getByOutputstatus(outAxes)
      CALL Axes%getByGroup("ps",outAxesPS)
      
      nxyzp = nxp*nyp*nzp
      nxyp  = nxp*nyp

      nz = nzp-1

      dxi = 1./deltax
      dyi = 1./deltay
      ALLOCATE(wsavex(4*nxpg+100),wsavey(4*nypg+100))
      wsavex = 0.0
      wsavey = 0.0

      ASSOCIATE(xm => xm%d(:), xt => xt%d(:), ym => ym%d(:), yt => yt%d(:),  &
                zm => zm%d(:), zt => zt%d(:), dzt => dzt%d(:), dzm => dzm%d(:)        )
      
        !
        ! define xm array for grid 1 from deltax
        !
        xm(1) = -float(max(nxpg-2,1))*.5*deltax+xoffset(wrxid)*deltax
        DO i = 2, nxp-1
           xm(i) = xm(i-1)+deltax
        END DO
        xm(nxp) = 2*xm(nxp-1)-xm(nxp-2)
        !
        ! define ym array for grid 1 from deltay
        !
        ym(1) = -float(max(nypg-2,1))*.5*deltay+yoffset(wryid)*deltay
        DO j = 2, nyp-1
           ym(j) = ym(j-1)+deltay
        END DO
        ym(nyp) = 2*ym(nyp-1)-ym(nyp-2)
        
        !
        !      define where the momentum points will lie in vertical
        !
        SELECT CASE (abs(igrdtyp))
           !
           ! Read in grid spacings from a file
           !
        CASE(3)
           OPEN(1,file='zm_grid_in',status='old',form='formatted')
           DO k = 1, nzp
              READ(1,*) zm(k)
           END DO
           CLOSE(1)
           IF (zm(1) /= 0.) THEN
              IF (myid == 0) PRINT *, 'ABORTING:  Error in input grid'
              CALL appl_abort(0)
           END IF
           !
           ! Tschebyschev Grid with vertical size given by dzmax
           !
        CASE(2)
           zm(1) = 0.
           nchby = nzp-3
           DO k = 1, nzp-1
              zm(k+1) = cos( ((2.*nchby - 1. - 2.*(k-1))*2.*asin(1.))/(2.*nchby))
              zm(k+1) = (zm(k+1)+1.)*dzmax/2.
           END DO
           zm(nzp-1) = dzmax
           zm(nzp)   = dzmax + zm(2)*zm(2)/(zm(3)-zm(2))
           !
           ! define zm array for grid 1 from deltaz and dzrat, if dzrat is
           ! negative compress grid so that dzmin is the grid spacing in a 100m
           ! interval below dzmax.  In both CASEs stretcvh grid uniformly by the
           ! ration |dzrat| above dzmax
           !
        CASE(1)
           dzmin=0.
           zm(1) = 0.
           zm(2) = deltaz
           zb = dzmax+100.
           IF (dzrat < 0.) THEN
              dzmin = -float(int(dzrat))
              dzrat =  dzrat+dzmin-1
              kmax = int(log(deltaz/dzmin)/log(abs(dzrat)))
              zb = dzmax-100.
              DO k = 1, kmax
                 zb = zb-dzmin*abs(dzrat)**k
              END DO
           END IF
           
           dz = deltaz
           DO k = 3, nzp
              IF(zm(k-1) > zb .AND. zm(k-1) < dzmax)then
                 dz = max(dzmin,dz/abs(dzrat))
              ELSE IF (zm(k-1) >= dzmax) THEN
                 dz = dz*abs(dzrat)
              END IF
              zm(k) = zm(k-1)+dz
           END DO
        CASE DEFAULT
           zm(1) = 0.
           DO k = 2, nzp ! Fixed: used to start from 1
              zm(k) = zm(k-1)+deltaz
           END DO
        END SELECT
        !
        ! Grid Points for Thermal Points (T-Grid):
        !
        DO i = 2, nxp
           xt(i) = .5*(xm(i)+xm(i-1))
        END DO
        xt(1) = 1.5*xm(1)-.5*xm(2)
        !
        DO j = 2, nyp
           yt(j) = .5*(ym(j)+ym(j-1))
        END DO
        yt(1) = 1.5*ym(1)-.5*ym(2)
        !
        IF (igrdtyp < 0) THEN
           !
           ! Read in grid spacings from a file
           !
           OPEN(2,file='zt_grid_in',status='old',form='formatted')
           DO k = 1, nzp
              READ(2,*) zt(k)
           END DO
           CLOSE(2)
        ELSE
           !
           ! calculate where the thermo points will lie based on geometric
           ! interpolation from the momentum points
           !
           DO k = 1, nzp
              zmnvc(k) = zm(k)
           END DO
           zmnvc(0) = -(zmnvc(2)-zmnvc(1))**2 /(zmnvc(3)-zmnvc(2))
           zmnvc(-1) = zmnvc(0)-(zmnvc(1)-zmnvc(0))**2 /(zmnvc(2)-zmnvc(1))
           zmnvc(nzp+1) = zmnvc(nzp)+(zmnvc(nzp)-zmnvc(nzp-1))**2              &
                /(zmnvc(nzp-1)-zmnvc(nzp-2))
           DO k = 1, nzp
              dzrfm = sqrt(sqrt((zmnvc(k+1)-zmnvc(k)) /(zmnvc(k-1)-zmnvc(k-2))))
              zt(k) = zmnvc(k-1)+(zmnvc(k)-zmnvc(k-1))/(1.+dzrfm)
           END DO
        END IF
        !
        ! compute other arrays based on the vertical grid.
        !   dzm: inverse of distance between thermal points k+1 and k
        !   dzt: inverse of distance between momentum points k and k-1
        !
        DO k = 1, nzp-1
           dzm(k) = 1./(zt(k+1)-zt(k))
        END DO
        dzm(nzp) = dzm(nzp-1)*dzm(nzp-1)/dzm(nzp-2)
        DO k = 2, nzp
           dzt(k) = 1./(zm(k)-zm(k-1))
        END DO
        dzt(1) = dzt(2)*dzt(2)/dzt(3)
        !
        ! set timesteps
        !
        dtl = dtlong
        dtlv = 2.*dtl
        dtlt = dtl
        !
        IF(myid == 0) THEN
           WRITE(6,fm1)
           WRITE(6,fm2) nxpg-4, deltax, 2.*xt(nxp-2)
           WRITE(6,fm3) nypg-4, deltay, 2.*yt(nyp-2)
           WRITE(6,fm4) nzp,zm(2)-zm(1),zm(nzp)
           WRITE(6,fm5) dtl
           WRITE(6,fm6) level
        END IF

      END ASSOCIATE
        
   END SUBROUTINE define_grid

 
   !
   ! ----------------------------------------------------------------------
   ! Subroutine newsclr:  This routine updates the scalar pointer to the
   ! value corresponding to the next scalar in the scalar table
   !
   SUBROUTINE newsclr(iscnum)

      INTEGER, INTENT(in) :: iscnum

      CHARACTER(len=20), PARAMETER :: name = "newsclr"

      a_sp => a_sclrp(:,:,:,iscnum)
      a_st => a_sclrt(:,:,:,iscnum)

      RETURN
   END SUBROUTINE newsclr



END MODULE grid

