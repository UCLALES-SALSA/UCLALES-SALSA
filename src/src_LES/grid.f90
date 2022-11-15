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
  USE mo_submctl, ONLY : spec, nbins, ncld, nprc, nice,     &
                         in1a, fn2a, in2b, fn2b, ica, icb, fca, fcb,  &
                         aerobins, cloudbins, precpbins, icebins,     &
                         ice_theta_dist,lssecice,lsicerimespln,lsicedropfrac
  
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
  
  LOGICAL :: lnudging = .FALSE.  ! Master switch for nudging scheme

  LOGICAL :: lemission = .FALSE. ! Master switch for aerosol emission
  
  LOGICAL :: lpback = .FALSE.       ! Master switch for running bulk microphysics
                                    ! in piggybacking mode, while level == 4
  INTEGER :: pbncsrc = 0            ! Source of CDNC for piggybacking microphysics.
                                    ! 0: use the constant CCN parameter
                                    ! 1: Use the CNDC taken from master microphysics (SALSA)
                                    !    gridpoint by gridpoint
                                    ! 2: Use the CDNC taken from master microphysics (SALSA) 
                                    !    as a domain mean

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
  
  CHARACTER (len=150) :: expnme = 'Default' ! Experiment name
  CHARACTER (len=150) :: filprf = 'x'       ! File Prefix
  CHARACTER (len=7)  :: runtype = 'INITIAL'! Run Type SELECTion

  ! Output file list; given here instead of mo_output.f90 to avoid cyclic dependencies
  CHARACTER(len=50) :: varlist_main(100),varlist_ps(100), varlist_ts(100)  ! 

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
  TYPE(FieldArray) :: AxesPS            ! All the axes that can be assigned to profile statistics files
  TYPE(FieldArray) :: AxesTS            ! Same for timeseries output
  TYPE(FieldArray) :: outAxes           ! Subset from above assigned for output
  TYPE(FieldArray) :: outAxesPS         ! Subset from Axes assigned for ps-file output
  TYPE(FieldArray) :: outAxesTS         ! Same for ts-files
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
     USE mo_aux_state, ONLY : setInitialProfiles

      CHARACTER(len=20), PARAMETER :: name = "define_vars"
      INTEGER :: nc
      
      ! Instanciate the field arrays
      BasicState = FieldArray()
      
      nc = 0
      ! Juha: Number of prognostic tracers for SALSA
      !       Aerosol bins + Cloud bins + gas compound tracers
      IF (level >= 4) THEN
         nc = spec%getNSpec(type="wet")
         nsalsa = (nc+1)*nbins + (nc+1)*ncld + (nc+1)*nprc + 5        ! (nc+1) for the mass tracers + number concentration
         IF (level == 5) nsalsa = nsalsa + (nc+1+1)*nice              ! (nc+1+1)*nice for RIMED ICE
         IF (level == 5 .AND. ice_theta_dist) nsalsa = nsalsa + nbins+ncld+nprc  ! If contact angle distributions for heterogeneous ice nucleation, 
                                                                                 ! add one more tracer for the "IN deficit fraction"
         IF (level == 5 .AND. lssecice%switch) nsalsa = nsalsa + 2.*nice

      END IF
         
      ! Initial condition vectors
      CALL setInitialProfiles(BasicState,nzp)
      CALL BasicState%getByOutputstatus(outBasicState)

      memsize = 2*nxyzp ! complexarray in pressure solver

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
         ! ... + bulk slave precip number and mass (lpback = .TRUE.)
         IF (lpback) nscl = nscl + 2
         
         
         ALLOCATE (a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
         a_sclrp(:,:,:,:) = 0.
         a_sclrt(:,:,:,:) = 0.

      END IF ! level
           
   END SUBROUTINE define_vars
   !
   !----------------------------------------------------------------------
   !
   SUBROUTINE define_grid
     USE mo_aux_state, ONLY : setGridSpacings,xt,xm,yt,ym,zt,zm,dzt,dzm,      &
                              aea, aeb, cla, clb, prc, ice, aetot, cltot
      USE mpi_interface, ONLY : xoffset, yoffset, wrxid, wryid, nxpg, nypg,   &
                                appl_abort, myid
      USE mo_structured_datatypes, ONLY : FloatArray1d
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
      CALL setGridSpacings(Axes,level,nzp,nxp,nyp)
      CALL Axes%getByOutputstatus(outAxes)
      CALL Axes%getByGroup("ps",AxesPS)
      CALL Axes%getByGroup("ts",AxesTS)
      CALL AxesPS%getByOutputstatus(outAxesPS)
      CALL AxesTS%getByOutputstatus(outAxesTS)

      nxyzp = nxp*nyp*nzp
      nxyp  = nxp*nyp

      nz = nzp-1

      dxi = 1./deltax
      dyi = 1./deltay
      ALLOCATE(wsavex(4*nxpg+100),wsavey(4*nypg+100))
      wsavex = 0.0
      wsavey = 0.0

      !
      ! define xm array for grid 1 from deltax
      !
      xm%d(1) = -float(max(nxpg-2,1))*.5*deltax+xoffset(wrxid)*deltax
      DO i = 2, nxp-1
         xm%d(i) = xm%d(i-1)+deltax
      END DO
      xm%d(nxp) = 2*xm%d(nxp-1)-xm%d(nxp-2)
      !
      ! define ym array for grid 1 from deltay
      !
      ym%d(1) = -float(max(nypg-2,1))*.5*deltay+yoffset(wryid)*deltay
      DO j = 2, nyp-1
         ym%d(j) = ym%d(j-1)+deltay
      END DO
      ym%d(nyp) = 2*ym%d(nyp-1)-ym%d(nyp-2)
      
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
            READ(1,*) zm%d(k)
         END DO
         CLOSE(1)
         IF (zm%d(1) /= 0.) THEN
            IF (myid == 0) PRINT *, 'ABORTING:  Error in input grid'
            CALL appl_abort(0)
         END IF
         !
         ! Tschebyschev Grid with vertical size given by dzmax
         !
      CASE(2)
         zm%d(1) = 0.
         nchby = nzp-3
         DO k = 1, nzp-1
            zm%d(k+1) = cos( ((2.*nchby - 1. - 2.*(k-1))*2.*asin(1.))/(2.*nchby))
            zm%d(k+1) = (zm%d(k+1)+1.)*dzmax/2.
         END DO
         zm%d(nzp-1) = dzmax
         zm%d(nzp)   = dzmax + zm%d(2)*zm%d(2)/(zm%d(3)-zm%d(2))
         !
         ! define zm array for grid 1 from deltaz and dzrat, if dzrat is
         ! negative compress grid so that dzmin is the grid spacing in a 100m
         ! interval below dzmax.  In both CASEs stretcvh grid uniformly by the
         ! ration |dzrat| above dzmax
         !
      CASE(1)
         dzmin=0.
         zm%d(:) = 0.
         zm%d(1) = 0.
         zm%d(2) = deltaz
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
            IF(zm%d(k-1) > zb .AND. zm%d(k-1) < dzmax)then
               dz = max(dzmin,dz/abs(dzrat))
            ELSE IF (zm%d(k-1) >= dzmax) THEN
               dz = dz*abs(dzrat)
            END IF
            zm%d(k) = zm%d(k-1)+dz
         END DO
      CASE DEFAULT
         zm%d(1) = 0.
         DO k = 2, nzp ! Fixed: used to start from 1
            zm%d(k) = zm%d(k-1)+deltaz
         END DO
      END SELECT

      !
      ! Grid Points for Thermal Points (T-Grid):
      !
      DO i = 2, nxp
         xt%d(i) = .5*(xm%d(i)+xm%d(i-1))
      END DO
      xt%d(1) = 1.5*xm%d(1)-.5*xm%d(2)
      !
      DO j = 2, nyp
         yt%d(j) = .5*(ym%d(j)+ym%d(j-1))
      END DO
      yt%d(1) = 1.5*ym%d(1)-.5*ym%d(2)
      !
      IF (igrdtyp < 0) THEN
         !
         ! Read in grid spacings from a file
         !
         OPEN(2,file='zt_grid_in',status='old',form='formatted')
         DO k = 1, nzp
            READ(2,*) zt%d(k)
         END DO
         CLOSE(2)
      ELSE
         !
         ! calculate where the thermo points will lie based on geometric
         ! interpolation from the momentum points
         !
         DO k = 1, nzp
            zmnvc(k) = zm%d(k)
         END DO
         zmnvc(0) = -(zmnvc(2)-zmnvc(1))**2 /(zmnvc(3)-zmnvc(2))
         zmnvc(-1) = zmnvc(0)-(zmnvc(1)-zmnvc(0))**2 /(zmnvc(2)-zmnvc(1))
         zmnvc(nzp+1) = zmnvc(nzp)+(zmnvc(nzp)-zmnvc(nzp-1))**2              &
              /(zmnvc(nzp-1)-zmnvc(nzp-2))
         DO k = 1, nzp
            dzrfm = sqrt(sqrt((zmnvc(k+1)-zmnvc(k)) /(zmnvc(k-1)-zmnvc(k-2))))
            zt%d(k) = zmnvc(k-1)+(zmnvc(k)-zmnvc(k-1))/(1.+dzrfm)
         END DO
      END IF
      !
      ! compute other arrays based on the vertical grid.
      !   dzm: inverse of distance between thermal points k+1 and k
      !   dzt: inverse of distance between momentum points k and k-1
      !
      DO k = 1, nzp-1
         dzm%d(k) = 1./(zt%d(k+1)-zt%d(k))
      END DO
      dzm%d(nzp) = dzm%d(nzp-1)*dzm%d(nzp-1)/dzm%d(nzp-2)
      DO k = 2, nzp
         dzt%d(k) = 1./(zm%d(k)-zm%d(k-1))
      END DO
      dzt%d(1) = dzt%d(2)*dzt%d(2)/dzt%d(3)
      !
      ! set timesteps
      !
      dtl = dtlong
      dtlv = 2.*dtl
      dtlt = dtl
      !

      ! Set bin diameter grids 
      IF (level >= 4) THEN         
         aea%d(:) = aerobins(in1a:fn2a)
         aeb%d(:) = aerobins(in2b:fn2b)
         aetot%d(:) = aerobins(in1a:fn2b)
         cla%d(:) = cloudbins(ica%cur:fca%cur)
         clb%d(:) = cloudbins(icb%cur:fcb%cur)
         cltot%d(:) = cloudbins(ica%cur:fcb%cur)
         prc%d(:) = precpbins(1:nprc)
      END IF
      IF (level > 4) THEN
         ice%d(:) = icebins(1:nice)
      END IF
      
      IF(myid == 0) THEN
         WRITE(6,fm1)
         WRITE(6,fm2) nxpg-4, deltax, 2.*xt%d(nxp-2)
         WRITE(6,fm3) nypg-4, deltay, 2.*yt%d(nyp-2)
         WRITE(6,fm4) nzp,zm%d(2)-zm%d(1),zm%d(nzp)
         WRITE(6,fm5) dtl
         WRITE(6,fm6) level
      END IF
        
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

