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
  
  USE ncio!, ONLY : open_nc, define_nc
  USE mo_structured_datatypes, ONLY : FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d
  USE classFieldArray, ONLY : FieldArray
  USE mo_diag_state
  USE mo_progn_state
  USE mo_aux_state
  USE mo_submctl, ONLY : spec, nbins

  
  IMPLICIT NONE
  
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
  
  LOGICAL            :: salsa_b_bins = .FALSE.  ! This is brought here temporarily from stat.f90 
  CHARACTER (len=7), ALLOCATABLE, SAVE :: sanal(:)
  CHARACTER (len=200) :: expnme = 'Default' ! Experiment name
  CHARACTER (len=200) :: filprf = 'x'       ! File Prefix
  CHARACTER (len=7)  :: runtype = 'INITIAL'! Run Type SELECTion
  
  CHARACTER (len=7),  PRIVATE :: v_snm = 'sxx    '
  CHARACTER (len=200), PRIVATE :: fname  
  INTEGER, PRIVATE, SAVE  ::  nrec0, nvar0, nbase=15

  ! Grid definitions
  ! -----------------------------------------------------
  INTEGER           :: nz, nxyzp, nxyp
  REAL              :: dxi, dyi, dtl, dtlv, dtlt, umean, vmean, psrf

  REAL, ALLOCATABLE :: spng_wfct(:), spng_tfct(:)
  
  ! Some zero arrays ice with level < 5
  REAL, ALLOCATABLE, TARGET :: tmp_icep(:,:,:,:), tmp_icet(:,:,:,:)
  !
  ! Prognostic vector variables (past, current and tendency)
  !
  REAL, ALLOCATABLE, TARGET :: a_up(:,:,:),a_uc(:,:,:),a_ut(:,:,:)
  REAL, ALLOCATABLE, TARGET :: a_vp(:,:,:),a_vc(:,:,:),a_vt(:,:,:)
  REAL, ALLOCATABLE, TARGET :: a_wp(:,:,:),a_wc(:,:,:),a_wt(:,:,:)
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
  
  ! Field arrays for organizing the prognostic and diagnostic variables and their attributes and output status
  ! ------------------------------------------------------------------------------------------------------------
  TYPE(FieldArray) :: Prog
  TYPE(FieldArray) :: Diag
  ! ------------------------------------------------------------------------------------------------------------

  ! Auxiliary FieldArray instances for pre-selected groups
  TYPE(FieldArray) :: SALSA_tracers_4d

  
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
  INTEGER, SAVE :: ncid0,ncid_s
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
      USE mo_submctl, ONLY : nbins,ncld,nprc,  & ! Number of aerosol and hydrometeor size bins for SALSA
                             nice          ! number of ice size bins for SALSA

      CHARACTER(len=20), PARAMETER :: name = "define_vars"
      INTEGER :: memsize
      INTEGER :: zz
      INTEGER :: nc
      INTEGER :: st_salsa,en_salsa ! start and end indices for SALSA tracers
      
      ! Instanciate the main field arrays
      Prog = FieldArray()
      Diag = FieldArray()

      nc = 0
      ! Juha: Number of prognostic tracers for SALSA
      !       Aerosol bins + Cloud bins + gas compound tracers
      IF (level >= 4) THEN
         nc = spec%getNSpec(type="wet")
         nsalsa = (nc+1)*nbins + (nc+1)*ncld + (nc+1)*nprc + 5
         IF (level == 5) nsalsa = nsalsa + (nc+1+1)*nice ! (nc+1+1)*nice for RIMED ICE 
      END IF

      ! Initial condition vectors
      CALL setInitialProfiles(nzp)

      memsize = 2*nxyzp ! complexarray in pressure solver

      ! Vector variables
      ALLOCATE (a_up(nzp,nxp,nyp),a_vp(nzp,nxp,nyp),a_wp(nzp,nxp,nyp))
      a_up(:,:,:) = 0.
      a_vp(:,:,:) = 0.
      a_wp(:,:,:) = 0.

      ALLOCATE (a_uc(nzp,nxp,nyp),a_vc(nzp,nxp,nyp),a_wc(nzp,nxp,nyp))
      a_uc(:,:,:) = 0.
      a_vc(:,:,:) = 0.
      a_wc(:,:,:) = 0.

      ALLOCATE (a_ut(nzp,nxp,nyp),a_vt(nzp,nxp,nyp),a_wt(nzp,nxp,nyp))
      a_ut(:,:,:) = 0.
      a_vt(:,:,:) = 0.
      a_wt(:,:,:) = 0.

      ! Diagnostic scalars
      CALL setDiagnosticVariables(Diag,memsize,level,iradtyp,nzp,nxp,nyp)

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

      CALL setPrognosticVariables(a_sclrp,a_sclrt,Prog,memsize,level,isgstyp,nzp,nxp,nyp,nscl)

      CALL Prog%getGroup("SALSA_4d",SALSA_tracers_4d)
      
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
      CALL setGridSpacings(level,nzp,nxp,nyp)
      
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
   ! Subroutine init_anal:  Defines the netcdf Analysis file
   !
   ! Modified for level 4.
   ! Juha Tonttila, FMI, 2014
   !
   !
   SUBROUTINE init_anal(time,salsa_b_bins)
      USE mpi_interface, ONLY : myid, ver, author, info
      REAL, INTENT(in) :: time
      LOGICAL, INTENT(in) :: salsa_b_bins
      
      CHARACTER(len=20), PARAMETER :: name = "init_anal"
     
   END SUBROUTINE init_anal
   !
   ! ----------------------------------------------------------------------
   ! Subroutine close_anal:  Closes netcdf anal file
   !
   INTEGER FUNCTION close_anal()
     USE netcdf
     
     CHARACTER(len=20), PARAMETER :: name = "close_anal"

     close_anal = nf90_close(ncid0)
      
   END FUNCTION close_anal

   !
   ! ----------------------------------------------------------------------
   ! Subroutine Write_anal:  Writes the netcdf Analysis file
   !
   ! Modified for levels 4 and 5
   ! Juha Tonttila, FMI, 2014
   !
   !
   SUBROUTINE write_anal(time)
      USE netcdf
      REAL, INTENT(in) :: time
      CHARACTER(len=20), PARAMETER :: name = "write_anal"

   END SUBROUTINE write_anal
   !
   ! ----------------------------------------------------------------------
   ! Subroutine write_hist:  This subroutine writes a binary history file
   !
   SUBROUTINE write_hist(htype, time)

     USE mpi_interface, ONLY : appl_abort, myid, wrxid, wryid
   
     !Ali
     !These are must be written and read from history file
     !for consistent nudging initialization
     USE nudg_defs, ONLY : theta_ref, rv_ref, u_ref, v_ref, aero_ref, &
                           ndg_theta, ndg_rv, ndg_u, ndg_v, ndg_aero
     
      INTEGER :: errcode = -17

      INTEGER, INTENT (in) :: htype
      REAL, INTENT (in)    :: time

      CHARACTER(len=20), PARAMETER :: name = "write_hist"

      CHARACTER (len=80) :: hname

      INTEGER :: n, iblank, ii, jj, kk,nn
      !
      ! create and open a new output file.
      !
      WRITE(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
      hname = trim(hname)//'.'//trim(filprf)

      SELECT CASE(htype)
         CASE DEFAULT
            hname = trim(hname)//'.iflg'
         CASE(0)
            hname = trim(hname)//'.R'
         CASE(1)
            hname = trim(hname)//'.rst'
         CASE(2)
            iblank=index(hname,' ')
            WRITE(hname(iblank:iblank+7),'(a1,i6.6,a1)') '.', int(time), 's'
      END SELECT
      !
      ! Write fields
      !
      IF (myid == 0) PRINT "(//' ',49('-')/,' ',/,'   History write to: ',A30)" &
                            ,hname
      OPEN(10,file=trim(hname), form='unformatted')

      WRITE(10) time,th00,umean,vmean,dtl,level,isgstyp,iradtyp,nzp,nxp,nyp,nscl
      WRITE(10) xt%d, xm%d, yt%d, ym%d, zt%d, zm%d, dn0%d, th0%d, u0%d, v0%d, pi0%d, &
                pi1%d, rt0%d, psrf,sst,W1,W2,W3 ! added by Zubair

      WRITE(10) a_ustar%d, a_tstar%d, a_rstar%d

      WRITE(10) a_pexnr%d
      WRITE(10) a_press%d
      WRITE(10) a_theta%d

      WRITE(10) a_up
      WRITE(10) a_vp
      WRITE(10) a_wp
      WRITE(10) a_uc
      WRITE(10) a_vc
      WRITE(10) a_wc

      DO n = 1, nscl
         CALL newsclr(n)  
         WRITE(10) a_sp
      END DO

      IF (ndg_theta%nudgetype > 0) THEN
        DO n = 1, nzp
          WRITE(10) theta_ref(n)
        END DO  
      END IF

      IF (ndg_rv%nudgetype > 0) THEN
        DO n = 1, nzp
          WRITE(10) rv_ref(n)
        END DO  
      END IF

      IF (ndg_u%nudgetype > 0) THEN
        DO n = 1, nzp
          WRITE(10) u_ref(n)
        END DO  
      END IF

      IF (ndg_v%nudgetype > 0) THEN
        DO n = 1, nzp
          WRITE(10) v_ref(n)
        END DO  
      END IF

      IF (level > 3 .AND. ndg_aero%nudgetype > 0) THEN
        WRITE(10) nbins
        DO n = 1, nzp
          DO nn = 1, nbins
            WRITE(10) aero_ref(n,nn)    
          END DO
        END DO
      END IF
    
      IF ( ASSOCIATED(a_rv%d)   ) WRITE(10) a_rv%d
      IF ( ASSOCIATED(a_rc%d)   ) WRITE(10) a_rc%d
      IF ( ASSOCIATED(a_rflx%d) ) WRITE(10) a_rflx%d
      CLOSE(10)

      IF (myid == 0 .AND. htype < 0) THEN
         PRINT *, 'CFL Violation'
         CALL appl_abort(errcode)
      END IF

      RETURN
   END SUBROUTINE write_hist
   !
   ! ----------------------------------------------------------------------
   ! Subroutine read_hist:  This subroutine reads a binary history file
   !
   !                        Modified for level 4
   !                Juha Tonttila, FMI, 20140828
   !

   SUBROUTINE read_hist(time, hfilin)

     USE mpi_interface, ONLY : appl_abort, myid, wrxid, wryid
     
      !Ali
      !These are must be written and read from history file
      !for consistent nudging initialization
      USE nudg_defs, ONLY : theta_ref, rv_ref, u_ref, v_ref, aero_ref, &
                           ndg_theta, ndg_rv, ndg_u, ndg_v, ndg_aero
      
      CHARACTER(len=80), INTENT(in) :: hfilin
      REAL, INTENT(out)             :: time

      CHARACTER(len=20), PARAMETER :: name = "read_hist"

      CHARACTER (len=80) :: hname
      INTEGER :: n, nxpx, nypx, nzpx, nsclx, iradx, isgsx, lvlx, ii, jj, kk
      LOGICAL :: exans
      REAL    :: umx, vmx, thx
      INTEGER :: nn, nnbins
      !
      ! open input file.
      !

      WRITE(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
      hname = trim(hname)//'.'//trim(hfilin)

      inquire(file=trim(hname),exist=exans)
      IF (.NOT. exans) THEN
         PRINT *,'ABORTING: History file', trim(hname),' not found'
         CALL appl_abort(0)
      ELSE
         OPEN(10,file=trim(hname),status='old',form='unformatted')
         READ(10) time,thx,umx,vmx,dtl,lvlx,isgsx,iradx,nzpx,nxpx,nypx,nsclx

         IF (nxpx /= nxp .OR. nypx /= nyp .OR. nzpx /= nzp)  THEN
            IF (myid == 0) PRINT *, nxp, nyp, nzp, nxpx, nypx, nzpx
            CALL appl_abort(-1)
         END IF

         READ(10) xt%d, xm%d, yt%d, ym%d, zt%d, zm%d, dn0%d, th0%d, u0%d, v0%d, pi0%d, pi1%d, rt0%d, psrf,sst,W1,W2,W3

         READ(10) a_ustar%d, a_tstar%d, a_rstar%d

         READ(10) a_pexnr%d
         READ(10) a_press%d
         READ(10) a_theta%d

         READ(10) a_up
         READ(10) a_vp
         READ(10) a_wp
         READ(10) a_uc
         READ(10) a_vc
         READ(10) a_wc

         DO n = 1, nscl
            CALL newsclr(n)
            IF (n <= nsclx) READ(10) a_sp
         END DO

         IF (ndg_theta%nudgetype > 0) THEN
           ALLOCATE(theta_ref(nzp))
           DO n = 1, nzp
              READ(10) theta_ref(n)
           END DO
         END IF

         IF (ndg_rv%nudgetype > 0) THEN
           ALLOCATE(rv_ref(nzp))
           DO n = 1, nzp
              READ(10) rv_ref(n)
           END DO
         END IF

         IF (ndg_u%nudgetype > 0) THEN
           ALLOCATE(u_ref(nzp))
           DO n = 1, nzp
              READ(10) u_ref(n)
           END DO
         END IF

         IF (ndg_v%nudgetype > 0) THEN
           ALLOCATE(v_ref(nzp))
           DO n = 1, nzp
              READ(10) v_ref(n)
           END DO
         END IF

         IF (level > 3 .AND. ndg_aero%nudgetype > 0) THEN
           READ(10) nnbins
           ALLOCATE(aero_ref(nzp,nnbins))
           DO n = 1, nzp
             DO nn = 1, nbins
               READ(10) aero_ref(n,nn)    
             END DO
           END DO
         END IF

         
         DO n = nscl+1, nsclx
            READ(10)
         END DO

         IF (lvlx > 0 .AND. lvlx < 4) THEN
            IF (level > 0 .AND. lvlx < 4) THEN
               READ(10) a_rv%d
            ELSE
               READ(10)
            END IF
         END IF
         IF (lvlx > 1) THEN
            IF (level > 1) THEN
               READ(10) a_rc%d
            ELSE
               READ(10)
            END IF
         END IF
         IF (iradx > 0) THEN
            IF (iradtyp > 0) THEN
               READ(10) a_rflx%d
            ELSE
               READ(10)
            END IF
         END IF

         CLOSE(10)
         !
         ! adjust namelist and basic state appropriately
         !
         IF (thx /= th00) THEN
            IF (myid == 0) PRINT "('  th00 changed  -  ',2f8.2)",th00,thx
            a_tp%d(:,:,:) = a_tp%d(:,:,:) + thx - th00
         END IF
         IF (umx /= umean) THEN
            IF (myid == 0) PRINT "('  umean changed  -  ',2f8.2)",umean,umx
            a_up = a_up + umx - umean
         END IF
         IF (vmx /= vmean) THEN
            IF (myid == 0) PRINT "('  vmean changed  -  ',2f8.2)",vmean,vmx
            a_vp = a_vp + vmx - vmean
         END IF
         dtlv = 2.*dtl
         dtlt = dtl

      END IF

   END SUBROUTINE read_hist
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
   !
   ! -----------------------------------
   ! Subroutine bulkMixrat: Find and calculate
   ! the total mixing ratio of a given compound
   ! in aerosol particles or hydrometeors
   !
   ! Juha Tonttila, FMI, 2015
   ! Jaakko Ahola, FMI, 2015
   SUBROUTINE bulkMixrat(icomp,ipart,itype,mixrat)

      USE mo_submctl, ONLY : ncld,nbins,nprc,      &
                             ica,fca,icb,fcb,      &
                             ira,fra,nice,iia,fia, &
                             in1a,in2b,fn2a,fn2b
      USE util, ONLY : getMassIndex

      CHARACTER(len=*), INTENT(in) :: icomp  ! This should be either:
                                             ! SO4,OC,NO,NH,BC,DU,SS,H2O.

      CHARACTER(len=*), INTENT(in) :: ipart  ! This should be either:
                                             ! aerosol,cloud,rain,ice,snow
      CHARACTER(len=*), INTENT(in) :: itype  ! Select bin regime: a or b

      REAL, INTENT(out) :: mixrat(nzp,nxp,nyp)

      CHARACTER(len=20), PARAMETER :: name = "bulkMixrat"

      INTEGER :: istr,iend,mm

      mixrat = 0.

      ! Determine multiplier
      mm = spec%getIndex(icomp)

      ! Given in kg/kg
      SELECT CASE(ipart)
         CASE('aerosol')
            IF (itype == 'ab') THEN
               istr = getMassIndex(nbins,in1a,mm)   
               iend = getMassIndex(nbins,fn2b,mm)    
            ELSE IF (itype == 'a') THEN
               istr = getMassIndex(nbins,in1a,mm)
               iend = getMassIndex(nbins,fn2a,mm)  
            ELSE IF (itype == 'b') THEN
               istr = getMassIndex(nbins,in2b,mm)
               iend = getMassIndex(nbins,fn2b,mm)
            ELSE
               STOP 'bulkMixrat: Invalid aerosol bin regime SELECTion'
            END IF
            mixrat(:,:,:) = SUM(a_maerop%d(:,:,:,istr:iend),DIM=4)
         CASE('cloud')
            IF (itype == 'ab') THEN
               istr = getMassIndex(ncld,ica%cur,mm)
               iend = getMassIndex(ncld,fcb%cur,mm)
            ELSE IF (itype == 'a') THEN
               istr = getMassIndex(ncld,ica%cur,mm)  
               iend = getMassIndex(ncld,fca%cur,mm)
            ELSE IF (itype == 'b') THEN
               istr = getMassIndex(ncld,icb%cur,mm)
               iend = getMassIndex(ncld,fcb%cur,mm)
            ELSE
               STOP 'bulkMixrat: Invalid cloud bin regime SELECTion'
            END IF
            mixrat(:,:,:) = SUM(a_mcloudp%d(:,:,:,istr:iend),DIM=4)
         CASE('precp')
            istr = getMassIndex(nprc,ira,mm)
            iend = getMassIndex(nprc,fra,mm)
            mixrat(:,:,:) = SUM(a_mprecpp%d(:,:,:,istr:iend),DIM=4)
         CASE('ice')
            istr = getMassIndex(nice,iia,mm)
            iend = getMassIndex(nice,fia,mm)
            mixrat(:,:,:) = SUM(a_micep%d(:,:,:,istr:iend),DIM=4)
      END SELECT

   END SUBROUTINE bulkMixrat
   !
   ! ----------------------------------------------
   ! Subroutine binSpecMixrat: Calculate the mixing
   ! ratio of selected aerosol species in individual
   ! bins.
   !
   ! Juha Tonttila, FMI, 2015
   SUBROUTINE binSpecMixrat(ipart,icomp,ibin,mixr)
      USE mo_submctl, ONLY : ncld, nbins, nprc, nice
      USE util, ONLY : getMassIndex

      CHARACTER(len=*), INTENT(in) :: icomp  ! This should be either:
                                             ! SO4,OC,NO,NH,BC,DU,SS,H2O.

      CHARACTER(len=*), INTENT(in) :: ipart  ! This should be either:
                                             ! aerosol,cloud,rain,ice
      INTEGER, INTENT(in) :: ibin

      REAL, INTENT(out)   :: mixr(nzp,nxp,nyp)

      CHARACTER(len=20), PARAMETER :: name = "binSpecMixrat"

      INTEGER :: mm

      ! Determine multipliers
      mm = spec%getIndex(icomp)

      SELECT CASE(ipart)
         CASE('aerosol')
            mixr(:,:,:) = a_maerop%d(:,:,:,getMassIndex(nbins,ibin,mm))
         CASE('cloud')
            mixr(:,:,:) = a_mcloudp%d(:,:,:,getMassIndex(ncld,ibin,mm))
         CASE('precp')
            mixr(:,:,:) = a_mprecpp%d(:,:,:,getMassIndex(nprc,ibin,mm))
         CASE('ice')
            mixr(:,:,:) = a_micep%d(:,:,:,getMassIndex(nice,ibin,mm))
      END SELECT

   END SUBROUTINE binSpecMixrat
   !
   ! ----------------------------------------------
   ! Subroutine binMixrat: Calculate the total dry or wet
   ! Mass concentration for individual bins
   !
   ! Juha Tonttila, FMI, 2015
   ! Tomi Raatikainen, FMI, 2016
   SUBROUTINE binMixrat(ipart,itype,ibin,ii,jj,kk,sumc)
      USE util, ONLY : getBinTotalMass
      USE mo_submctl, ONLY : ncld,nbins,nprc,nice
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(in) :: ipart
      CHARACTER(len=*), INTENT(in) :: itype
      INTEGER, INTENT(in) :: ibin,ii,jj,kk
      REAL, INTENT(out) :: sumc

      CHARACTER(len=20), PARAMETER :: name = "binMixrat"

      INTEGER :: iend

      REAL, POINTER :: tmp(:) => NULL()

      IF (itype == 'dry') THEN
         iend = spec%getNSpec(type="dry")   ! dry CASE
      ELSE IF (itype == 'wet') THEN
         iend = spec%getNSpec(type="wet")   ! wet CASE
         IF (ipart == 'ice') &
              iend = iend+1   ! For ice, take also rime
      ELSE
         STOP 'Error in binMixrat!'
      END IF

      SELECT CASE(ipart)
         CASE('aerosol')
            tmp => a_maerop%d(kk,ii,jj,1:iend*nbins)
            CALL getBinTotalMass(nbins,iend,ibin,tmp,sumc)
         CASE('cloud')
            tmp => a_mcloudp%d(kk,ii,jj,1:iend*ncld)
            CALL getBinTotalMass(ncld,iend,ibin,tmp,sumc)
         CASE('precp')
            tmp => a_mprecpp%d(kk,ii,jj,1:iend*nprc)
            CALL getBinTotalMass(nprc,iend,ibin,tmp,sumc)
         CASE('ice')
            tmp => a_micep%d(kk,ii,jj,1:iend*nice)
            CALL getBinTotalMass(nice,iend,ibin,tmp,sumc) 
         CASE DEFAULT
            STOP 'bin mixrat error'
      END SELECT

      tmp => NULL()
      
   END SUBROUTINE binMixrat

   !
   ! ----------------------------------------------
   ! Subroutine bulkNumc: Calculate the total number
   ! concentration of particles of given type
   !
   ! Juha Tonttila, FMI, 2015
   !
   SUBROUTINE bulkNumc(ipart,itype,numc)
      USE mo_submctl, ONLY : ica,fca,icb,fcb,ira,fra, &
                             iia,fia,in1a,in2b,fn2a,fn2b

      CHARACTER(len=*), INTENT(in) :: ipart
      CHARACTER(LEN=*), INTENT(in) :: itype
      REAL, INTENT(out) :: numc(nzp,nxp,nyp)
      
      CHARACTER(len=20), PARAMETER :: name = "bulkNumc"

      INTEGER :: istr,iend

      istr = 0
      iend = 0

      ! Outputs #/kg
      ! No concentration limits (nlim or prlim) for number

      SELECT CASE(ipart)
         CASE('aerosol')
            IF (itype == 'ab') THEN ! Note: 1a and 2a, 2b combined
               istr = in1a
               iend = fn2b
            ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
               istr = in1a
               iend = fn2a
            ELSE IF (itype == 'b') THEN ! 2b
               istr = in2b
               iend = fn2b
            END IF
            numc(:,:,:) = SUM(a_naerop%d(:,:,:,istr:iend),DIM=4)
         CASE('cloud')
            IF (itype == 'ab') THEN ! Note: 1a and 2a, 2b combined
               istr = ica%cur
               iend = fcb%cur
            ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
               istr = ica%cur
               iend = fca%cur
            ELSE IF (itype == 'b') THEN
               istr = icb%cur
               iend = fcb%cur
            END IF
            numc(:,:,:) = SUM(a_ncloudp%d(:,:,:,istr:iend),DIM=4)
         CASE('precp')
            istr = ira
            iend = fra
            numc(:,:,:) = SUM(a_nprecpp%d(:,:,:,istr:iend),DIM=4)
         CASE('ice')
            istr = iia
            iend = fia
            numc(:,:,:) = SUM(a_nicep%d(:,:,:,istr:iend),DIM=4)
      END SELECT

   END SUBROUTINE bulkNumc

   !
   ! -------------------------------------------------
   ! SUBROUTINE meanRadius
   ! Gets the mean wet (water=nspec+1) radius for particles - this function is for outputs only
   !
   SUBROUTINE meanRadius(ipart,itype,rad)
     USE mo_submctl, ONLY : spec,nbins,ncld,nprc,nice,     &
                            ica,fca,icb,fcb,ira,fra,       &
                            iia,fia,in1a,fn2a,in2b,fn2b,   &
                            nlim,prlim
     IMPLICIT NONE
     
     CHARACTER(len=*), INTENT(in) :: ipart
     CHARACTER(len=*), INTENT(in) :: itype
     REAL, INTENT(out) :: rad(nzp,nxp,nyp)
     CHARACTER(len=20), PARAMETER :: name = "meanRadius"
     
     INTEGER :: istr,iend
     INTEGER :: nspec
     
     nspec = spec%getNSpec(type="wet")
     
     rad = 0.
     
     SELECT CASE(ipart)
     CASE('aerosol')
        
        IF (itype == 'ab') THEN ! Note: 1a, 2a and 2b combined
           istr = in1a
           iend = fn2b
        ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
           istr = in1a
           iend = fn2a
        ELSE IF (itype == 'b') THEN
           istr = in2b
           iend = fn2b
        ELSE
           STOP 'meanRadius: Invalid bin regime selection (aerosol)'
        END IF
        
        CALL getRadius(istr,iend,nbins,nspec,a_naerop,a_maerop,nlim,rad,1)
        
     CASE('cloud')
        
        IF (itype == 'ab') THEN
           istr = ica%cur
           iend = fcb%cur
        ELSE IF (itype == 'a') THEN
           istr = ica%cur
           iend = fca%cur
        ELSE IF (itype == 'b') THEN
           istr = icb%cur
           iend = fcb%cur
        ELSE
           STOP 'meanRadius: Invalid bin regime selection (cloud)'
        END IF
        
        CALL getRadius(istr,iend,ncld,nspec,a_ncloudp,a_mcloudp,nlim,rad,2)
        
     CASE('precp')
        
        istr = ira
        iend = fra
        
        CALL getRadius(istr,iend,nprc,nspec,a_nprecpp,a_mprecpp,prlim,rad,3)
        
     CASE('ice')

        istr = iia
        iend = fia
         
        CALL getRadius(istr,iend,nice,nspec+1,a_nicep,a_micep,prlim,rad,4)
               
     END SELECT
     
   CONTAINS
     
     SUBROUTINE getRadius(zstr,zend,nb,ns,numc,mass,numlim,zrad,flag)
       USE util, ONLY : getBinMassArray
       USE mo_particle_external_properties, ONLY : calcDiamLES
       USE mo_submctl, ONLY : pi6
       IMPLICIT NONE
       
       INTEGER, INTENT(in) :: nb, ns ! Number of bins (nb) and compounds (ns)
       INTEGER, INTENT(in) :: zstr,zend  ! Start and end index for averaging
       TYPE(FloatArray4d), INTENT(in) :: numc
       TYPE(FloatArray4d), INTENT(in) :: mass
       REAL, INTENT(in) :: numlim
       INTEGER, INTENT(IN) :: flag
       REAL, INTENT(out) :: zrad(nzp,nxp,nyp)
       
       INTEGER :: k,i,j,bin
       REAL :: tot, rwet, tmp(ns)
       REAL :: zlm(nb*ns),zln(nb) ! Local grid point binned mass and number concentrations 
              
       zrad(:,:,:)=0.
       DO j = 3,nyp-2
          DO i = 3,nxp-2
             DO k = 1,nzp
                zlm(:) = mass%d(k,i,j,:)
                zln(:) = numc%d(k,i,j,:)
                tot=0.
                rwet=0.
                DO bin = zstr,zend                  
                   IF (zln(bin)>numlim) THEN
                      tot=tot+zln(bin)
                      tmp(:) = 0.
                      CALL getBinMassArray(nb,ns,bin,zlm,tmp)
                      rwet=rwet+0.5*calcDiamLES(ns,zln(bin),tmp,flag)*zln(bin)
                   ENDIF
                ENDDO
                IF (tot>numlim) THEN
                   zrad(k,i,j) = rwet/tot
                ENDIF
             END DO
          END DO
       END DO
       
     END SUBROUTINE getRadius
   END SUBROUTINE meanRadius
   
   !
   ! ---------------------------------------------------
   ! SUBROUTINE getBinRadius
   ! Calculates wet radius for each bin in the whole domain - this function is for outputs only
   SUBROUTINE getBinRadius(nb,ns,numc,mass,numlim,zrad,flag)
     USE util, ONLY : getBinMassArray
     USE mo_particle_external_properties, ONLY : calcDiamLES
     USE mo_submctl, ONLY : pi6
     IMPLICIT NONE
     
     INTEGER, INTENT(in) :: nb, ns ! Number of bins (nb) and aerosol species (ns)
     TYPE(FloatArray4d), INTENT(in) :: numc
     TYPE(FloatArray4d), INTENT(in) :: mass
     REAL, INTENT(in) :: numlim
     INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4)
     REAL, INTENT(out) :: zrad(nzp,nxp,nyp,nb)

     INTEGER :: k,i,j,bin
     REAL :: tmp(ns)
     REAL :: zlm(nb*ns), zln(nb)
     
     zrad(:,:,:,:)=0.
     DO j = 3,nyp-2
        DO i = 3,nxp-2
           DO k = 1,nzp
              zlm(:) = mass%d(k,i,j,:)
              zln(:) = numc%d(k,i,j,:)
              DO bin = 1,nb
                 IF (zln(bin)>numlim) THEN
                    tmp(:) = 0.
                    CALL getBinMassArray(nb,ns,bin,zlm,tmp)
                    zrad(k,i,j,bin)=0.5*calcDiamLES(ns,zln(bin),tmp,flag)
                 ENDIF
              END DO
           END DO
        END DO
     END DO

   END SUBROUTINE getBinRadius
   

END MODULE grid

