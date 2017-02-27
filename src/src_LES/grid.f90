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
module grid

  use ncio, only : open_nc, define_nc

  USE class_componentIndex, ONLY : componentIndex

  implicit none
  !
  integer           :: nxp = 132           ! number of x points
  integer           :: nyp = 132           ! number of y points
  integer           :: nzp = 105           ! number of z points


  logical           :: nxpart = .true.     ! number of processors in x

  real              :: deltax = 35.        ! dx for basic grid
  real              :: deltay = 35.        ! dy for basic grid
  real              :: deltaz = 17.5       ! dz for basic grid
  real              :: dzrat  = 1.02       ! grid stretching ratio
  real              :: dzmax  = 1200.      ! height to start grid-stretching
  real              :: dtlong = 10.0       ! long timestep
  real              :: th00   = 288.       ! basic state temperature

  real              :: CCN = 150.e6

  LOGICAL           :: lbinanl = .FALSE.   ! Whether to write binned data to analysis files (takes a lot of space + mainly used for debugging)
  integer           :: igrdtyp = 1         ! vertical grid type
  integer           :: isgstyp = 1         ! sgs model type
  integer           :: iradtyp = 0         ! radiation model type
  integer           :: level   = 0         ! thermodynamic level
  integer           :: naddsc  = 0         ! number of additional scalars;
  INTEGER           :: nsalsa  = 0         ! Number of tracers for SALSA
  integer           :: nfpt = 10           ! number of rayleigh friction points
  real              :: distim = 300.0      ! dissipation timescale

  real              :: sst=283.   ! Surface temperature      added by Zubair Maalick
  real            :: W1 = 0.9   !Water content
  real            :: W2 = 0.9
  real            :: W3 = 0.9



  character (len=7), allocatable, save :: sanal(:)
  character (len=80):: expnme = 'Default' ! Experiment name
  character (len=80):: filprf = 'x'       ! File Prefix
  character (len=7) :: runtype = 'INITIAL'! Run Type Selection

  REAL              :: Tspinup = 7200.    ! Spinup period in seconds (added by Juha)


  character (len=7),  private :: v_snm='sxx    '
  character (len=80), private :: fname

  integer, private, save  ::  nrec0, nvar0, nbase=15

  integer           :: nz, nxyzp, nxyp
  real              :: dxi, dyi, dtl, dtlv, dtlt, umean, vmean, psrf
  real, allocatable :: xt(:), xm(:), yt(:), ym(:), zt(:), zm(:), dzt(:), dzm(:)
  real, allocatable :: u0(:), v0(:), pi0(:), pi1(:), th0(:), dn0(:), rt0(:)
  real, allocatable :: spng_wfct(:), spng_tfct(:)
  !
  ! velocity variables (past, current and tendency)
  !
  real, allocatable, target :: a_up(:,:,:),a_uc(:,:,:),a_ut(:,:,:)
  real, allocatable, target :: a_vp(:,:,:),a_vc(:,:,:),a_vt(:,:,:)
  real, allocatable, target :: a_wp(:,:,:),a_wc(:,:,:),a_wt(:,:,:)
  !
  ! wsave variables used in fft in x and y directons
  !
  real, allocatable :: wsavex(:), wsavey(:)
  !
  ! prognostic scalar variables
  !
  real, pointer :: a_tp(:,:,:),a_tt(:,:,:)
  real, pointer :: a_rp(:,:,:),a_rt(:,:,:)  !Juha: In standard version this is the TOTAL water content.
                                            !      With SALSA this is taken as just the water VAPOUR content,
                                            !      in order not to over-specify the problem.
  real, pointer :: a_rpp(:,:,:),a_rpt(:,:,:)
  real, pointer :: a_npp(:,:,:),a_npt(:,:,:)
  real, pointer :: a_qp(:,:,:),a_qt(:,:,:)
  real, pointer :: a_sp(:,:,:),a_st(:,:,:)

  ! Juha: SALSA tracers
  !---------------------------------------------------------------------------
  ! -- Masses given in kg/kg, number concentrations in #/kg
  ! -- Each size bin/species will be treated as a separate tracer.
  ! -- NOTE: Volume concentration arrays are reduced to 4 dims.
  !          The 4th dim contains all the size bins sequentially for
  !          each aerosol species  + water
  !
  !          Gas tracers are contained sequentially in dimension
  !          4 as: 1. SO4, 2. HNO3, 3. NH3, 4. OCNV, 5. OCSV

  ! Prognostic tracers
  ! -- Number concentrations
  REAL, POINTER :: a_naerop(:,:,:,:), a_naerot(:,:,:,:),   &
                   a_ncloudp(:,:,:,:), a_ncloudt(:,:,:,:), &
                   a_nprecpp(:,:,:,:), a_nprecpt(:,:,:,:), &
                   a_nicep(:,:,:,:),   a_nicet(:,:,:,:), &
                   a_nsnowp(:,:,:,:),  a_nsnowt(:,:,:,:)
  ! -- Volume concentrations
  REAL, POINTER :: a_maerop(:,:,:,:), a_maerot(:,:,:,:),   &
                   a_mcloudp(:,:,:,:), a_mcloudt(:,:,:,:), &
                   a_mprecpp(:,:,:,:), a_mprecpt(:,:,:,:), &
                   a_micep(:,:,:,:), a_micet(:,:,:,:), &
                   a_msnowp(:,:,:,:), a_msnowt(:,:,:,:)
  ! -- Gas compound tracers
  REAL, POINTER :: a_gaerop(:,:,:,:), a_gaerot(:,:,:,:)

  ! Diagnostic tracers
  REAL, ALLOCATABLE :: a_Rawet(:,:,:,:), a_Radry(:,:,:,:),  &
                       a_Rcwet(:,:,:,:), a_Rcdry(:,:,:,:),  &
                       a_Rpwet(:,:,:,:), a_Rpdry(:,:,:,:),  &
                       a_Riwet(:,:,:,:), a_Ridry(:,:,:,:),  &
                       a_Rswet(:,:,:,:), a_Rsdry(:,:,:,:)

  ! Some stuff for tendency formulation
  REAL, ALLOCATABLE :: a_vactd(:,:,:,:), a_nactd(:,:,:,:)

  ! Particle component index tables
  TYPE(componentIndex) :: prtcl ! Contains "getIndex" which gives the index for a given
                               ! aerosol component name, i.e. 1:SO4, 2:OC, 3:BC, 4:DU,
                               ! 5:SS, 6:NO, 7:NH, 8:H2O

  !---------------------------------------------------------------------------

  real, allocatable, target :: a_sclrp(:,:,:,:),a_sclrt(:,:,:,:)
   !
  ! 3d diagnostic quantities
  !
  real, allocatable, target :: a_theta(:,:,:)  ! dry potential temp (k)
  real, allocatable :: a_pexnr(:,:,:)  ! perturbation exner func
  real, allocatable :: a_press(:,:,:)  ! pressure (hpa)
  real, allocatable :: a_rc(:,:,:)     ! Total cloud water
  real, allocatable :: a_ri(:,:,:)     ! Total ice cloud content
  real, allocatable :: a_rv(:,:,:)     ! water vapor (used only for levels < 4!)
  REAL, ALLOCATABLE :: a_srp(:,:,:)    ! Total rain water for use with LEVEL 4 (Diagnostic scalar!!)
  REAL, ALLOCATABLE :: a_snrp(:,:,:)   ! Total number of rain drops for use with LEVEL 4 (Diagnostic scalar!!)
  REAL, ALLOCATABLE :: a_srs(:,:,:)    ! Total snow for use with SALSA
  REAL, ALLOCATABLE :: a_snrs(:,:,:)   ! Total number of snow particles for use with LEVEL 5 (Diagnostic scalar!!)
  REAL, ALLOCATABLE :: a_rh(:,:,:)     ! Relative humidity
  REAL, ALLOCATABLE :: a_rsl(:,:,:)     ! water saturation vapor mixing ratio
  REAL, ALLOCATABLE :: a_rhi(:,:,:)     ! Relative humidity over ice
  REAL, ALLOCATABLE :: a_rsi(:,:,:)     ! ice saturation vapor mixing ratio
  REAL, ALLOCATABLE :: a_dn(:,:,:)     ! Air density (for normalizing concentrations according to mass, levels < 4!)
  !
  ! scratch arrays
  !
  real, allocatable, dimension (:,:,:) :: a_rflx, a_sflx, &
       a_temp, a_temp0 ! store temperatures of previous timestep
  !
  !
  real, allocatable :: a_ustar(:,:)
  real, allocatable :: a_tstar(:,:)
  real, allocatable :: a_rstar(:,:)
  real, allocatable :: uw_sfc(:,:)
  real, allocatable :: vw_sfc(:,:)
  real, allocatable :: ww_sfc(:,:)
  real, allocatable :: wt_sfc(:,:)
  real, allocatable :: wq_sfc(:,:)
  real, allocatable :: precip(:,:,:), snowin(:,:,:), albedo(:,:)

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
  integer :: nscl = 1
  integer, save :: ncid0,ncid_s
  !
contains
  !
  !----------------------------------------------------------------------
  ! SUBROUTINE define_vars
  !
  ! Modified for level 4
  ! Juha Tonttila, FMI, 2014.
  !
  subroutine define_vars

    use mpi_interface, only :myid
    USE mo_submctl, ONLY : nbins,ncld,nprc,  & ! Number of aerosol and hydrometeor size bins for SALSA
                               nice,nsnw,        & ! number of ice and snow size bins for SALSA
                               nspec, maxspec, listspec
    USE class_ComponentIndex, ONLY : ComponentIndexConstructor,  &
                                     GetNcomp, IsUsed


    integer :: memsize
    INTEGER :: zz
    INTEGER :: nc

    ! Juha: Number of prognostic tracers for SALSA
    !            Aerosol bins + Cloud bins + gas compound tracers

    IF (level >= 4) THEN
       ! Create index tables for different aerosol components (can be fetched by name using getIndex)
       CALL ComponentIndexConstructor(prtcl, nspec, maxspec, listspec)
       nc = GetNcomp(prtcl)

       nsalsa = (nc+2)*nbins + (nc+2)*ncld + (nc+2)*nprc + (nc+2)*nice + (nc+2)*nsnw + 5

    END IF

    ! Juha: Stuff that's allocated for all configurations
    !----------------------------------------------------------
    allocate (u0(nzp),v0(nzp),pi0(nzp),pi1(nzp),th0(nzp),dn0(nzp),rt0(nzp))

    memsize = 2*nxyzp ! complex array in pressure solver

    allocate (a_up(nzp,nxp,nyp),a_vp(nzp,nxp,nyp),a_wp(nzp,nxp,nyp))
    a_up(:,:,:) = 0.
    a_vp(:,:,:) = 0.
    a_wp(:,:,:) = 0.

    allocate (a_uc(nzp,nxp,nyp),a_vc(nzp,nxp,nyp),a_wc(nzp,nxp,nyp))
    a_uc(:,:,:) = 0.
    a_vc(:,:,:) = 0.
    a_wc(:,:,:) = 0.

    allocate (a_ut(nzp,nxp,nyp),a_vt(nzp,nxp,nyp),a_wt(nzp,nxp,nyp))
    a_ut(:,:,:) = 0.
    a_vt(:,:,:) = 0.
    a_wt(:,:,:) = 0.

    allocate (a_theta(nzp,nxp,nyp),a_pexnr(nzp,nxp,nyp),a_press(nzp,nxp,nyp))
    a_theta(:,:,:) = 0.
    a_pexnr(:,:,:) = 0.
    a_press(:,:,:) = 0.

    memsize = memsize + nxyzp*13 !

    if (iradtyp > 0 ) then
       allocate (a_rflx(nzp,nxp,nyp))
       a_rflx(:,:,:) = 0.
       memsize = memsize + nxyzp
    end if
    if (iradtyp >= 3) then
       allocate (a_sflx(nzp,nxp,nyp),albedo(nxp,nyp))
       a_sflx(:,:,:) = 0.
       albedo(:,:) = 0.
       memsize = memsize + nxyzp + nxyp
    end if

    allocate (a_temp(nzp,nxp,nyp),a_temp0(nzp,nxp,nyp),a_rsl(nzp,nxp,nyp))
    a_temp(:,:,:) = 0.
    a_temp0(:,:,:) = 0.
    a_rsl(:,:,:) = 0.
    memsize = memsize + nxyzp*3

    ! Juha: Stuff that's allocated if SALSA is NOT used
    !-----------------------------------------------------
    IF (level < 4) THEN

       if (level >= 0) then
          allocate (a_rv(nzp,nxp,nyp))
          a_rv(:,:,:) = 0.
          memsize = memsize + nxyzp
          if (level > 1) then
             allocate (a_rc(nzp,nxp,nyp))
             a_rc(:,:,:) = 0.
             memsize = memsize + nxyzp
          end if
       end if

       nscl = nscl+naddsc
       if (level   > 0) nscl = nscl+1
       if (level   > 2) nscl = nscl+2
       if (isgstyp > 1) nscl = nscl+1

       allocate (a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
       a_sclrp(:,:,:,:) = 0.
       a_sclrt(:,:,:,:) = 0.

       a_tp=>a_sclrp(:,:,:,1)
       a_tt=>a_sclrt(:,:,:,1)
       if (level >= 0) then
          a_rp=>a_sclrp(:,:,:,2)
          a_rt=>a_sclrt(:,:,:,2)
       end if
       if (level >= 3) then
          a_rpp=>a_sclrp(:,:,:,3)
          a_rpt=>a_sclrt(:,:,:,3)
          a_npp=>a_sclrp(:,:,:,4)
          a_npt=>a_sclrt(:,:,:,4)
       end if
       if (isgstyp > 1) then
          a_qp=>a_sclrp(:,:,:,nscl - naddsc)
          a_qt=>a_sclrt(:,:,:,nscl - naddsc)
       end if

    !Juha: Stuff that's allocated when SALSA is used
    !---------------------------------------------------
    ELSE IF (level >= 4) THEN

       nc = GetNcomp(prtcl) ! number of aerosol components used. For allocations + 1 for water.

       allocate (a_rc(nzp,nxp,nyp), a_srp(nzp,nxp,nyp), a_snrp(nzp,nxp,nyp),     &
                 a_Rawet(nzp,nxp,nyp,nbins),a_Radry(nzp,nxp,nyp,nbins),          &
                 a_Rcwet(nzp,nxp,nyp,ncld), a_Rcdry(nzp,nxp,nyp,ncld),           &
                 a_Rpwet(nzp,nxp,nyp,nprc), a_Rpdry(nzp,nxp,nyp,nprc),           &
                 a_rh(nzp,nxp,nyp),a_dn(nzp,nxp,nyp), &
                 a_nactd(nzp,nxp,nyp,ncld), a_vactd(nzp,nxp,nyp,(nc+1)*ncld)  )

       a_rc(:,:,:) = 0.
       a_srp(:,:,:) = 0.
       a_snrp(:,:,:) = 0.
       a_Rawet(:,:,:,:) = 1.e-10
       a_Radry(:,:,:,:) = 1.e-10
       a_Rcwet(:,:,:,:) = 1.e-10
       a_Rcdry(:,:,:,:) = 1.e-10
       a_Rpwet(:,:,:,:) = 1.e-10
       a_Rpdry(:,:,:,:) = 1.e-10
       a_rh(:,:,:) = 0.
       a_dn(:,:,:) = 0.
       a_nactd(:,:,:,:) = 0.
       a_vactd(:,:,:,:) = 0.
       memsize = memsize + 4*nxyzp + 3*nbins*nxyzp + 3*ncld*nxyzp + nxyzp*(nc+1)*ncld + 2*nprc*nxyzp

       allocate ( a_Riwet(nzp,nxp,nyp,nice), a_Ridry(nzp,nxp,nyp,nice),           &
                  a_Rswet(nzp,nxp,nyp,nsnw), a_Rsdry(nzp,nxp,nyp,nsnw),           &
                  a_ri(nzp,nxp,nyp), a_rsi(nzp,nxp,nyp), a_rhi(nzp,nxp,nyp),      &
                  a_srs(nzp,nxp,nyp), a_snrs(nzp,nxp,nyp)  )  ! ice'n'snow
       a_ri(:,:,:) = 0.
       a_rsi(:,:,:) = 0.
       a_rhi(:,:,:) = 0.
       a_srs(:,:,:) = 0.
       a_snrs(:,:,:) = 0.
       a_Riwet(:,:,:,:) = 1.e-10
       a_Ridry(:,:,:,:) = 1.e-10
       a_Rswet(:,:,:,:) = 1.e-10
       a_Rsdry(:,:,:,:) = 1.e-10
       memsize = memsize + 5*nxyzp + 2*nice*nxyzp + 2*nsnw*nxyzp

       ! Total number of prognostic scalars: temp + total water + SALSA + tke(?)
       nscl = 2 + nsalsa
       if (isgstyp > 1) nscl = nscl+1

       allocate (a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
       a_sclrp(:,:,:,:) = 0.
       a_sclrt(:,:,:,:) = 0.

       a_tp=>a_sclrp(:,:,:,1)
       a_tt=>a_sclrt(:,:,:,1)
       a_rp=>a_sclrp(:,:,:,2)
       a_rt=>a_sclrt(:,:,:,2)

       if (isgstyp > 1) then
          a_qp=>a_sclrp(:,:,:,nscl - nsalsa)
          a_qt=>a_sclrt(:,:,:,nscl - nsalsa)
       end if

       !JT: Set the pointers for prognostic SALSA variables (levels 4 & 5)
       zz = nscl-nsalsa
       a_naerop => a_sclrp(:,:,:,zz+1:zz+nbins)
       a_naerot => a_sclrt(:,:,:,zz+1:zz+nbins)

       zz = zz+nbins
       a_maerop => a_sclrp(:,:,:,zz+1:zz+(nc+1)*nbins)
       a_maerot => a_sclrt(:,:,:,zz+1:zz+(nc+1)*nbins)

       zz = zz+(nc+1)*nbins
       a_ncloudp => a_sclrp(:,:,:,zz+1:zz+ncld)
       a_ncloudt => a_sclrt(:,:,:,zz+1:zz+ncld)

       zz = zz+ncld
       a_mcloudp => a_sclrp(:,:,:,zz+1:zz+(nc+1)*ncld)
       a_mcloudt => a_sclrt(:,:,:,zz+1:zz+(nc+1)*ncld)

       zz = zz+(nc+1)*ncld
       a_nprecpp => a_sclrp(:,:,:,zz+1:zz+nprc)
       a_nprecpt => a_sclrt(:,:,:,zz+1:zz+nprc)

       zz = zz+nprc
       a_mprecpp => a_sclrp(:,:,:,zz+1:zz+(nc+1)*nprc)
       a_mprecpt => a_sclrt(:,:,:,zz+1:zz+(nc+1)*nprc)

       zz = zz+(nc+1)*nprc
       a_gaerop => a_sclrp(:,:,:,zz+1:zz+5)
       a_gaerot => a_sclrt(:,:,:,zz+1:zz+5)

       ! Level 5
       zz = zz+5
       a_nicep => a_sclrp(:,:,:,zz+1:zz+nice)
       a_nicet => a_sclrt(:,:,:,zz+1:zz+nice)

       zz = zz+nice

       a_micep => a_sclrp(:,:,:,zz+1:zz+(nc+1)*nice)
       a_micet => a_sclrt(:,:,:,zz+1:zz+(nc+1)*nice)

       zz = zz+(nc+1)*nice

       a_nsnowp => a_sclrp(:,:,:,zz+1:zz+nsnw)
       a_nsnowt => a_sclrt(:,:,:,zz+1:zz+nsnw)

       zz = zz+nsnw

       a_msnowp => a_sclrp(:,:,:,zz+1:zz+(nc+1)*nsnw)
       a_msnowt => a_sclrt(:,:,:,zz+1:zz+(nc+1)*nsnw)

       zz = zz+(nc+1)*nsnw

    END IF ! level

    !----------------------------------------------------

    allocate (a_ustar(nxp,nyp),a_tstar(nxp,nyp),a_rstar(nxp,nyp))
    allocate (uw_sfc(nxp,nyp),vw_sfc(nxp,nyp),ww_sfc(nxp,nyp))
    allocate (wt_sfc(nxp,nyp),wq_sfc(nxp,nyp))
    !allocate (ra(nxp,nyp))
    if (level >= 3) then
       allocate(precip(nzp,nxp,nyp))
       precip = 0.
       memsize = memsize + nxyzp
    end if

    allocate(snowin(nzp,nxp,nyp))
    memsize = memsize + nxyzp

    a_ustar(:,:) = 0.
    a_tstar(:,:) = 0.
    a_rstar(:,:) = 0.
    uw_sfc(:,:)  = 0.
    vw_sfc(:,:)  = 0.
    ww_sfc(:,:)  = 0.
    wt_sfc(:,:) = 0.
    wq_sfc(:,:) = 0.
    umean = 0.
    vmean = 0.

    memsize = memsize +  nxyzp*nscl*2 + 3*nxyp + nxyp*10

    if(myid == 0) then
       print "(//' ',49('-')/,' ',/3x,i3.3,' prognostic scalars')", nscl
       print "('   memory to be allocated  -  ',f8.3,' mbytes')", &
            memsize*1.e-6*kind(0.0)
    end if

  end subroutine define_vars
  !
  !----------------------------------------------------------------------
  !
  subroutine define_grid

    use mpi_interface, only: xoffset, yoffset, wrxid, wryid, nxpg, nypg,   &
         appl_abort, myid

    integer :: i,j,k,kmax,nchby
    real    :: dzrfm,dz,zb,dzmin
    real    :: zmnvc(-1:nzp+1)
    character (len=51) :: &
         fm1 = '(//" ",49("-")/,"   grid dimensions:"/)            ',      &
         fm2 = '("   nxp-4 = ",i3,", dx, dx = ",f8.1,",",f8.1," m")',      &
         fm3 = '("   nyp-4 = ",i3,", dy, dy = ",f8.1,",",f8.1," m")',      &
         fm4 = '("   nzp   = ",i3,", dz, dz = ",f8.1,",",f8.1," m")',      &
         fm5 = '("   timestep: ",f7.3,"s ")                        ',      &
         fm6 = '("   thermo level: ",i3)                        '

    nxyzp  = nxp*nyp*nzp
    nxyp   = nxp*nyp

    nz= nzp-1

    dxi=1./deltax
    dyi=1./deltay
    allocate(wsavex(4*nxpg+100),wsavey(4*nypg+100))
    wsavex=0.0
    wsavey=0.0

    !
    ! define xm array for grid 1 from deltax
    !
    allocate (xm(nxp))
    xm(1)=-float(max(nxpg-2,1))*.5*deltax+xoffset(wrxid)*deltax
    do i=2,nxp-1
       xm(i)=xm(i-1)+deltax
    end do
    xm(nxp)=2*xm(nxp-1)-xm(nxp-2)
    !
    ! define ym array for grid 1 from deltay
    !
    allocate (ym(nyp))
    ym(1)=-float(max(nypg-2,1))*.5*deltay+yoffset(wryid)*deltay
    do j=2,nyp-1
       ym(j)=ym(j-1)+deltay
    end do
    ym(nyp)=2*ym(nyp-1)-ym(nyp-2)

    !
    !      define where the momentum points will lie in vertical
    !
  allocate (zm(nzp))
  select case (abs(igrdtyp))
     !
     ! Read in grid spacings from a file
     !
  case(3)
     open (1,file='zm_grid_in',status='old',form='formatted')
     do k=1,nzp
        read (1,*) zm(k)
     end do
     close (1)
     if (zm(1) /= 0.) then
       if (myid == 0) print *, 'ABORTING:  Error in input grid'
       call appl_abort(0)
    end if
     !
     ! Tschebyschev Grid with vertical size given by dzmax
     !
  case(2)
     zm(1) = 0.
     nchby = nzp-3
     do k=1,nzp-1
        zm(k+1) = cos( ((2.*nchby - 1. - 2.*(k-1))*2.*asin(1.))/(2.*nchby))
        zm(k+1) = (zm(k+1)+1.)*dzmax/2.
     end do
     zm(nzp-1) = dzmax
     zm(nzp)   = dzmax + zm(2)*zm(2)/(zm(3)-zm(2))
     !
     ! define zm array for grid 1 from deltaz and dzrat, if dzrat is
     ! negative compress grid so that dzmin is the grid spacing in a 100m
     ! interval below dzmax.  In both cases stretcvh grid uniformly by the
     ! ration |dzrat| above dzmax
     !
  case(1)
     zm(1)=0.
     zm(2)=deltaz
     zb=dzmax+100.
     if (dzrat.lt.0.) then
        dzmin = -float(int(dzrat))
        dzrat =  dzrat+dzmin-1
        kmax = int(log(deltaz/dzmin)/log(abs(dzrat)))
        zb=dzmax-100.
        do k=1,kmax
           zb=zb-dzmin*abs(dzrat)**k
        end do
     end if

     dz=deltaz
     do k=3,nzp
        if(zm(k-1) > zb .and. zm(k-1) < dzmax)then
           dz=max(dzmin,dz/abs(dzrat))
        else if (zm(k-1) >= dzmax) then
           dz=dz*abs(dzrat)
        end if
        zm(k)=zm(k-1)+dz
     end do
  case default
     zm(1)=0.
     do k=2,nzp ! Fixed: used to start from 1
        zm(k)=zm(k-1)+deltaz
     end do
  end select
  !
  ! Grid Points for Thermal Points (T-Grid):
  !
  allocate (xt(nxp))
  do i=2,nxp
     xt(i)=.5*(xm(i)+xm(i-1))
  end do
  xt(1)=1.5*xm(1)-.5*xm(2)
  !
  allocate (yt(nyp))
  do j=2,nyp
     yt(j)=.5*(ym(j)+ym(j-1))
  end do
  yt(1)=1.5*ym(1)-.5*ym(2)
  !
  allocate (zt(nzp))
  if (igrdtyp .lt. 0) then
     !
     ! Read in grid spacings from a file
     !
     open (2,file='zt_grid_in',status='old',form='formatted')
     do k=1,nzp
        read (2,*) zt(k)
     end do
     close (2)
   else
     !
     ! calculate where the thermo points will lie based on geometric
     ! interpolation from the momentum points
     !
     do k=1,nzp
        zmnvc(k)=zm(k)
     end do
     zmnvc(0)=-(zmnvc(2)-zmnvc(1))**2 /(zmnvc(3)-zmnvc(2))
     zmnvc(-1)=zmnvc(0)-(zmnvc(1)-zmnvc(0))**2 /(zmnvc(2)-zmnvc(1))
     zmnvc(nzp+1)=zmnvc(nzp)+(zmnvc(nzp)-zmnvc(nzp-1))**2              &
                  /(zmnvc(nzp-1)-zmnvc(nzp-2))

     do k=1,nzp
       dzrfm=sqrt(sqrt((zmnvc(k+1)-zmnvc(k)) /(zmnvc(k-1)-zmnvc(k-2))))
       zt(k)=zmnvc(k-1)+(zmnvc(k)-zmnvc(k-1))/(1.+dzrfm)
     end do
  end if
  !
  ! compute other arrays based on the vertical grid.
  !   dzm: inverse of distance between thermal points k+1 and k
  !   dzt: inverse of distance between momentum points k and k-1
  !
  allocate (dzm(nzp))
  do k=1,nzp-1
     dzm(k)=1./(zt(k+1)-zt(k))
  end do
  dzm(nzp)=dzm(nzp-1)*dzm(nzp-1)/dzm(nzp-2)

  allocate (dzt(nzp))
  do k=2,nzp
     dzt(k)=1./(zm(k)-zm(k-1))
  end do
  dzt(1)=dzt(2)*dzt(2)/dzt(3)
  !
  ! set timesteps
  !
  dtl=dtlong
  dtlv=2.*dtl
  dtlt=dtl
  !
  if(myid == 0) then
     write(6,fm1)
     write(6,fm2) nxpg-4, deltax, 2.*xt(nxp-2)
     write(6,fm3) nypg-4, deltay, 2.*yt(nyp-2)

     write(6,fm4) nzp,zm(2)-zm(1),zm(nzp)
     write(6,fm5) dtl
     write(6,fm6) level
  endif

  end subroutine define_grid
  !
  ! ----------------------------------------------------------------------
  ! subroutine init_anal:  Defines the netcdf Analysis file
  !
  ! Modified for level 4.
  ! Juha Tonttila, FMI, 2014
  !
  !
  subroutine init_anal(time)

    use mpi_interface, only :myid, ver, author, info
    USE mo_submctl, ONLY : fn2a,fn2b,fca,fcb,fra, &
                               fia,fib,fsa
    USE class_ComponentIndex, ONLY : IsUsed
    integer, parameter :: nnames = 24
    integer, parameter :: salsa_nn = 104
    character (len=7), save :: sbase(nnames) =  (/ &
         'time   ','zt     ','zm     ','xt     ','xm     ','yt     '   ,& ! 1
         'ym     ','u0     ','v0     ','dn0    ','u      ','v      '   ,& ! 7
         'w      ','t      ','p      ','q      ','l      ','r      '   ,& ! 13
         'f      ','i      ','s      '                                 ,& ! 19 ice'n'snow
         'n      ','stke   ','rflx   '/)                                  ! 22 total 24
    ! Added for SALSA
    character(len=7), save :: salsa_sbase(salsa_nn) = (/ &
         'time   ','zt     ','zm     ','xt     ','xm     ','yt     ',  &  ! 1 
         'ym     ','aea    ','aeb    ','cla    ','clb    ','prc    ',  &  ! 7
         'ica    ','icb    ','snw    ','u0     ','v0     ','dn0    ',  &  ! 13
         'u      ','v      ','w      ','t      ','p      ','q      ',  &  ! 19
         'l      ','r      ','f      ','i      ','s      ',         &  ! 25
         'S_RH   ','S_RHI  ','S_Nact ','S_Na   ','S_Naba ','S_Rwaa ',  &  ! 30
         'S_Rwaba','S_Nb   ','S_Nabb ','S_Rwab ','S_Rwabb','S_Nc   ',  &  ! 36
         'S_Ncba ','S_Ncbb ','S_Rwca ','S_Rwcb ','S_Rwcba','S_Rwcbb',  &  ! 42
         'S_Np   ','S_Npba ','S_Rwpa ','S_Rwpba','S_Nic  ','S_Niba ',  &  ! 48
         'S_Nibb ','S_Rwia ','S_Rwib ','S_Rwiba','S_Rwibb','S_Ns   ',  &  ! 54
         'S_Nsba ','S_Rwsa ','S_Rwsba','S_aSO4a','S_aNHa ','S_aNOa ',  &  ! 60
         'S_aOCa ','S_aBCa ','S_aDUa ','S_aSSa ','S_aSO4b','S_aNHb ',  &  ! 66
         'S_aNOb ','S_aOCb ','S_aBCb ','S_aDUb ','S_aSSb ','S_cSO4a',  &  ! 72
         'S_cNHa ','S_cNOa ','S_cOCa ','S_cBCa ','S_cDUa ','S_cSSa ',  &  ! 78
         'S_cSO4b','S_cNHb ','S_cNOb ','S_cOCb ','S_cBCb ','S_cDUb ',  &  ! 84
         'S_cSSb ','S_iSO4a','S_iNHa ','S_iNOa ','S_iOCa ','S_iBCa ',  &  ! 90
         'S_iDUa ','S_iSSa ','S_iSO4b','S_iNHb ','S_iNOb ','S_iOCb ',  &  ! 96
         'S_iBCb ','S_iDUb ','S_iSSb ' /)
         ! total 104

    LOGICAL, SAVE :: salsabool(salsa_nn)

    real, intent (in) :: time
    integer           :: nbeg, nend

    IF (level < 4) THEN  ! Standard operation for levels 1-3

       nvar0 = nbase + naddsc
       if (level >= 1) nvar0 = nvar0+1
       if (level >= 2) nvar0 = nvar0+1
       if (level == 3) nvar0 = nvar0+2
       if (isgstyp > 1) nvar0 = nvar0+1
       if (iradtyp > 1) nvar0 = nvar0+1

       allocate (sanal(nvar0))
       sanal(1:nbase) = sbase(1:nbase)

       nvar0 = nbase
       !
       ! add liquid water, which is a diagnostic variable, first
       !
       if (level >= 2) then
          nvar0 = nvar0+1
          sanal(nvar0) = sbase(nbase+2)
       end if
       !
       ! add additional scalars, in the order in which they appear in scalar
       ! table
       !
       if (level >= 1) then
          nvar0 = nvar0+1
          sanal(nvar0) = sbase(nbase+1)
       end if

       if (level == 3) then
          nvar0 = nvar0+1
          sanal(nvar0) = sbase(nbase+3)
          nvar0 = nvar0+1
          sanal(nvar0) = sbase(nbase+4+3)
       end if

       if (isgstyp > 1) then
          nvar0 = nvar0+1
          sanal(nvar0) = sbase(nbase+5+3)
       end if

       if (iradtyp > 2) then
          nvar0 = nvar0+1
          sanal(nvar0) = sbase(nbase+6+3)
       end if

       nbeg = nvar0+1
       nend = nvar0+naddsc
       do nvar0 = nbeg, nend
          write(v_snm(2:3),'(i2.2)') nvar0-nbeg
          sanal(nvar0) = v_snm
       end do
       nvar0=nend

    ELSE IF (level >= 4) THEN ! Operation with SALSA

       ! Make a boolean array masking the output variables that are used.
       ! This mainly concerns unused aerosol species.
       salsabool(:) = .TRUE.
       IF (.NOT. lbinanl) THEN
          salsabool(8:15) = .FALSE.
          salsabool((/34,36,38,40,42,43,46,47,49,51,53,54,57,58,60,62/)) = .FALSE.
       END IF

       IF (level < 5 ) THEN
          salsabool(13:15) = .FALSE. ! ica, icb, snw
          salsabool(27:29) = .FALSE.
          salsabool(91:104) = .FALSE. ! aerosols in ice particles
          salsabool(31) = .FALSE.
          salsabool(52:62) = .FALSE.
       END IF

       IF (.NOT. IsUsed(prtcl,'SO4')) &
            salsabool((/ 63, 70, 77, 84, 91,  98 /)) = .FALSE.

       IF (.NOT. IsUsed(prtcl,'NH'))  &
            salsabool((/ 64, 71, 78, 85, 92,  99 /)) = .FALSE.

       IF (.NOT. IsUsed(prtcl,'NO'))  &
            salsabool((/ 65, 72, 79, 86, 93,  100 /)) = .FALSE.

       IF (.NOT. IsUsed(prtcl,'OC'))  &
            salsabool((/ 66, 73, 80, 87, 94,  101 /)) = .FALSE.

       IF (.NOT. IsUsed(prtcl,'BC'))  &
            salsabool((/ 67, 74, 81, 88, 95,  102 /)) = .FALSE.

       IF (.NOT. IsUsed(prtcl,'DU'))  &
            salsabool((/ 68, 75, 82, 89, 96,  103 /)) = .FALSE.

       IF (.NOT. IsUsed(prtcl,'SS'))  &
            salsabool((/ 69, 76, 83, 90, 97,  104 /)) = .FALSE.

       nvar0 = COUNT(salsabool) + naddsc
       ALLOCATE(sanal(nvar0))

       sanal = PACK(salsa_sbase,salsabool)

    END IF

    fname =  trim(filprf)
    if(myid == 0) print                                                  &
            "(//' ',49('-')/,' ',/,'   Initializing: ',A20)",trim(fname)
    call open_nc( fname, expnme, time, (nxp-4)*(nyp-4), ncid0, nrec0, ver, author, info)

    IF (level < 4 .OR. .NOT. lbinanl) THEN
       call define_nc( ncid0, nrec0, nvar0, sanal, n1=nzp, n2=nxp-4, n3=nyp-4)

    ELSE IF (level == 4 .AND. lbinanl) THEN
       call define_nc( ncid0, nrec0, nvar0, sanal, n1=nzp, n2=nxp-4, n3=nyp-4,  &
                       inae_a=fn2a, inae_b=fn2b-fn2a, incld_a=fca%cur,          &
                       incld_b=fcb%cur-fca%cur, inprc=fra )
    ELSE IF (level == 5 .AND. lbinanl) THEN
        call define_nc( ncid0, nrec0, nvar0, sanal, n1=nzp, n2=nxp-4, n3=nyp-4,  &
                        inae_a=fn2a,  inae_b =fn2b-fn2a, incld_a=fca%cur,        &
                        incld_b=fcb%cur-fca%cur, inprc=fra, inice_a=fia%cur,     &
                        inice_b=fib%cur-fia%cur, insnw=fsa )
    END IF
    if (myid == 0) print *,'   ...starting record: ', nrec0


  end subroutine init_anal
  !
  ! ----------------------------------------------------------------------
  ! subroutine close_anal:  Closes netcdf anal file
  !
  integer function close_anal()

    use netcdf

    close_anal = nf90_close(ncid0)

  end function close_anal
  !
  ! ----------------------------------------------------------------------
  ! Subroutine Write_anal:  Writes the netcdf Analysis file
  !
  ! Modified for levels 4 and 5
  ! Juha Tonttila, FMI, 2014
  !
  !
  subroutine write_anal(time)
    use netcdf
    use mpi_interface, only : myid, appl_abort
    USE class_ComponentIndex, ONLY : IsUsed
    USE mo_submctl, ONLY : in1a,fn2a,in2b,fn2b, &
                               ica,fca,icb,fcb,ira,fra,       &
                               iia,fia,iib,fib,isa,fsa,       &
                               aerobins, cloudbins, precpbins, &
                               icebins, snowbins

    real, intent (in) :: time

    integer :: iret, VarID, nn, n
    integer :: ibeg(4), icnt(4), i1, i2, j1, j2
    INTEGER :: ibegsd(5), icntaea(5), icntaeb(5), icntcla(5), icntclb(5), icntpra(5), & ! Juha: For sizedistribution variables
           icntica(5), icnticb(5), icntsna(5)
    REAL :: zsum(nzp,nxp,nyp) ! Juha: Helper for computing bulk output diagnostics
    REAL :: zvar(nzp,nxp,nyp)

    icnt = (/nzp, nxp-4, nyp-4, 1/)
    icntaea = (/nzp,nxp-4,nyp-4, fn2a, 1 /)
    icntaeb = (/nzp,nxp-4,nyp-4, fn2b-fn2a, 1/)
    icntcla = (/nzp,nxp-4,nyp-4, fca%cur, 1/)
    icntclb = (/nzp,nxp-4,nyp-4, fcb%cur-fca%cur, 1/)
    icntpra = (/nzp,nxp-4,nyp-4, fra, 1/)
    icntica = (/nzp,nxp-4,nyp-4, iia%cur, 1/)
    icnticb = (/nzp,nxp-4,nyp-4, iib%cur-iia%cur, 1/)
    icntsna = (/nzp,nxp-4,nyp-4, fsa, 1/)
    ibeg = (/1  ,1  ,1  ,nrec0/)
    ibegsd = (/1,1,1,1,nrec0/)

    i1 = 3
    i2 = nxp-2
    j1 = 3
    j2 = nyp-2

    iret = nf90_inq_varid(ncid0, sanal(1), VarID)
    iret = nf90_put_var(ncid0, VarID, time, start=(/nrec0/))

    if (nrec0 == 1) then
       iret = nf90_inq_varid(ncid0, sanal(2), VarID)
       iret = nf90_put_var(ncid0, VarID, zt, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(3), VarID)
       iret = nf90_put_var(ncid0, VarID, zm, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(4), VarID)
       iret = nf90_put_var(ncid0, VarID, xt(i1:i2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(5), VarID)
       iret = nf90_put_var(ncid0, VarID, xm(i1:i2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(6), VarID)
       iret = nf90_put_var(ncid0, VarID, yt(j1:j2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(7), VarID)
       iret = nf90_put_var(ncid0, VarID, ym(j1:j2), start = (/nrec0/))

       IF (level >= 4 .AND. lbinanl) THEN

          iret = nf90_inq_varid(ncid0,sanal(8), VarID)
          iret = nf90_put_var(ncid0, VarID, aerobins(in1a:fn2a), start = (/nrec0/))

          iret = nf90_inq_varid(ncid0,sanal(9), VarID)
          iret = nf90_put_var(ncid0,VarID, aerobins(in2b:fn2b), start = (/nrec0/))

          iret = nf90_inq_varid(ncid0,sanal(10), VarID)
          iret = nf90_put_var(ncid0,VarID, cloudbins(ica%cur:fca%cur), start = (/nrec0/))

          iret = nf90_inq_varid(ncid0,sanal(11), VarID)
          iret = nf90_put_var(ncid0,VarID, cloudbins(icb%cur:fcb%cur), start = (/nrec0/))

          iret = nf90_inq_varid(ncid0,sanal(12), VarID)
          iret = nf90_put_var(ncid0,VarID, precpbins(ira:fra), start = (/nrec0/))

          IF (level == 5) THEN
             iret = nf90_inq_varid(ncid0,sanal(13), VarID)
             iret = nf90_put_var(ncid0,VarID, icebins(iia%cur:fia%cur), start = (/nrec0/))
             
             iret = nf90_inq_varid(ncid0,sanal(14), VarID)
             iret = nf90_put_var(ncid0,VarID, icebins(iib%cur:fib%cur), start = (/nrec0/))

             iret = nf90_inq_varid(ncid0,sanal(15), VarID)
             iret = nf90_put_var(ncid0,VarID, snowbins(isa:fsa), start = (/nrec0/))
          END IF

       END IF
    end if

    IF (level < 4) THEN
       iret = nf90_inq_varid(ncid0, sanal(8), VarID)
       iret = nf90_put_var(ncid0, VarID, u0, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(9), VarID)
       iret = nf90_put_var(ncid0, VarID, v0, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(10), VarID)
       iret = nf90_put_var(ncid0, VarID, dn0, start = (/nrec0/))

       iret = nf90_inq_varid(ncid0, sanal(11), VarID)
       iret = nf90_put_var(ncid0, VarID, a_up(:,i1:i2,j1:j2), start=ibeg,    &
            count=icnt)
       iret = nf90_inq_varid(ncid0, sanal(12), VarID)
       iret = nf90_put_var(ncid0, VarID, a_vp(:,i1:i2,j1:j2), start=ibeg,    &
            count=icnt)
       iret = nf90_inq_varid(ncid0, sanal(13), VarID)
       iret = nf90_put_var(ncid0, VarID, a_wp(:,i1:i2,j1:j2), start=ibeg,    &
            count=icnt)
       iret = nf90_inq_varid(ncid0, sanal(14), VarID)
       iret = nf90_put_var(ncid0, VarID, a_theta(:,i1:i2,j1:j2), start=ibeg, &
            count=icnt)
       iret = nf90_inq_varid(ncid0, sanal(15), VarID)
       iret = nf90_put_var(ncid0, VarID, a_press(:,i1:i2,j1:j2), start=ibeg, &
            count=icnt)

    ELSE IF (level >= 4) THEN
       iret = nf90_inq_varid(ncid0, 'u0', VarID)
       iret = nf90_put_var(ncid0, VarID, u0, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, 'v0', VarID)
       iret = nf90_put_var(ncid0, VarID, v0, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, 'dn0', VarID)
       iret = nf90_put_var(ncid0, VarID, dn0, start = (/nrec0/))

       iret = nf90_inq_varid(ncid0, 'u', VarID)
       iret = nf90_put_var(ncid0, VarID, a_up(:,i1:i2,j1:j2), start=ibeg,    &
            count=icnt)
       iret = nf90_inq_varid(ncid0, 'v', VarID)
       iret = nf90_put_var(ncid0, VarID, a_vp(:,i1:i2,j1:j2), start=ibeg,    &
            count=icnt)
       iret = nf90_inq_varid(ncid0, 'w', VarID)
       iret = nf90_put_var(ncid0, VarID, a_wp(:,i1:i2,j1:j2), start=ibeg,    &
            count=icnt)
       iret = nf90_inq_varid(ncid0, 't', VarID)
       iret = nf90_put_var(ncid0, VarID, a_theta(:,i1:i2,j1:j2), start=ibeg, &
            count=icnt)
       iret = nf90_inq_varid(ncid0, 'p', VarID)
       iret = nf90_put_var(ncid0, VarID, a_press(:,i1:i2,j1:j2), start=ibeg, &
            count=icnt)

    END IF


    IF (level < 4) THEN ! Normal operation for levels 1-3

       nn = nbase
       if (level >= 2)  then
          nn = nn+1
          iret = nf90_inq_varid(ncid0, 'l', VarID)
          iret = nf90_put_var(ncid0, VarID, a_rc(:,i1:i2,j1:j2), start=ibeg, &
               count=icnt)
       end if

       do n = 2, nscl
          nn = nn+1
          call newsclr(n)
          iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
          iret = nf90_put_var(ncid0,VarID,a_sp(:,i1:i2,j1:j2), start=ibeg,   &
               count=icnt)
       end do

       if (isgstyp > 1)  then
          nn = nn+1
          iret = nf90_inq_varid(ncid0, 'stke', VarID)
          iret = nf90_put_var(ncid0, VarID, a_qp(:,i1:i2,j1:j2), start=ibeg, &
               count=icnt)
       end if

       if (iradtyp > 1)  then
          nn = nn+1
          iret = nf90_inq_varid(ncid0, 'rflx', VarID)
          iret = nf90_put_var(ncid0, VarID, a_rflx(:,i1:i2,j1:j2), start=ibeg, &
               count=icnt)
       end if

       if (nn /= nvar0) then
          if (myid == 0) print *, 'ABORTING:  Anal write error'
          call appl_abort(0)
       end if

    ELSE IF (level >= 4) THEN ! Operation with SALSA

       IF (.TRUE.) THEN

       ! Relative humidity
       iret = nf90_inq_varid(ncid0,'S_RH',VarID)
       iret = nf90_put_var(ncid0,VarID,a_rh(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)
       
       IF ( level == 5) THEN !
          ! Relative humidity
          iret = nf90_inq_varid(ncid0,'S_RHI',VarID)
          iret = nf90_put_var(ncid0,VarID,a_rhi(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
       END IF
       
       ! Total water mixing ratio
       zvar(:,:,:) = a_rp(:,:,:) + &   ! Water vapor
                         a_rc(:,:,:) + &   ! Liquid water
                         a_srp(:,:,:)        ! Rain
       
       iret = nf90_inq_varid(ncid0,'q',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)
       
       ! Liquid water mixing ratio
       zvar(:,:,:) = a_rc(:,:,:)
       iret = nf90_inq_varid(ncid0,'l',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)
       
       ! Rain water mixing ratio
       zvar(:,:,:) = a_srp(:,:,:)
       iret = nf90_inq_varid(ncid0,'r',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)
     
       IF ( level == 5) THEN !
          
          ! Total ice mixing ratio
          zvar(:,:,:) = a_ri(:,:,:) + & ! Ice cloud content
                           a_srs(:,:,:)      ! Snow
          iret = nf90_inq_varid(ncid0,'f',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          ! Ice mixing ratio
          zvar(:,:,:) = a_ri(:,:,:)
          iret = nf90_inq_varid(ncid0,'i',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          ! Snow ice mixing ratio
          zvar(:,:,:) = a_srs(:,:,:)
          iret = nf90_inq_varid(ncid0,'s',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
       END IF
       
       ! Total number of cloud droplets
       CALL bulkNumc('cloud','a',zvar(:,:,:))
       zsum = zvar
       CALL bulkNumc('cloud','b',zvar(:,:,:))
       zsum = zsum + zvar
       iret = nf90_inq_varid(ncid0,'S_Nc',VarID)
       iret = nf90_put_var(ncid0,VarID,zsum(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)
       
       IF (lbinanl) THEN
          ! Cloud droplet size distribution (regime A)
          iret = nf90_inq_varid(ncid0,'S_Ncba',VarID)
          iret = nf90_put_var(ncid0,VarID,a_ncloudp(:,i1:i2,j1:j2,ica%cur:fca%cur), &
               start=ibegsd,count=icntcla)
          
          ! Cloud droplet size distribution (regime B)
          iret = nf90_inq_varid(ncid0,'S_Ncbb',VarID)
          iret = nf90_put_var(ncid0,VarID,a_ncloudp(:,i1:i2,j1:j2,icb%cur:fcb%cur), &
               start=ibegsd,count=icntclb)
       END IF
       
       ! Mean cloud droplet wet radius (regime A)
       CALL meanRadius('cloud','a',zvar(:,:,:))
       iret = nf90_inq_varid(ncid0,'S_Rwca',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)
       
       ! Mean cloud droplet wet radius (regime B)
       CALL meanRadius('cloud','b',zvar(:,:,:))
       iret = nf90_inq_varid(ncid0,'S_Rwcb',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)
       
       IF (lbinanl) THEN
          ! Cloud droplet bin wet radius (regime A)
          iret = nf90_inq_varid(ncid0,'S_Rwcba',VarID)
          iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,ica%cur:fca%cur), &
               start=ibegsd,count=icntcla)
          
          ! Cloud droplet bin wet radius (regime B)
          iret = nf90_inq_varid(ncid0,'S_Rwcbb',VarID)
          iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,icb%cur:fcb%cur), &
               start=ibegsd,count=icntclb)
       END IF
       
       ! Number of rain droplets
       CALL bulkNumc('precp','a',zvar(:,:,:))
       iret = nf90_inq_varid(ncid0,'S_Np',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)
       
       IF (lbinanl) THEN
          ! Precipitation size distribution
          iret = nf90_inq_varid(ncid0,'S_Npba',VarID)
          iret = nf90_put_var(ncid0,VarID,a_nprecpp(:,i1:i2,j1:j2,ira:fra), &
               start=ibegsd,count=icntpra)
       END IF
       
       ! Mean rain drop wet radius
       CALL meanRadius('precp','a',zvar(:,:,:))
       iret = nf90_inq_varid(ncid0,'S_Rwpa',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,  &
            count=icnt)
       
       IF (lbinanl) THEN
          ! Rain drop bin wet radius
          iret = nf90_inq_varid(ncid0,'S_Rwpba',VarID)
          iret = nf90_put_var(ncid0,VarID,a_Rpwet(:,i1:i2,j1:j2,ira:fra),  &
               start=ibegsd,count=icntpra)
       END IF
       
       ! Number of soluble aerosols (regime A)
       CALL bulkNumc('aerosol','a',zvar(:,:,:))
       iret = nf90_inq_varid(ncid0,'S_Na',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)
       
       ! Number of insoluble aerosols (regime B)
       CALL bulkNumc('aerosol','b',zvar(:,:,:))
       iret = nf90_inq_varid(ncid0,'S_Nb',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)
       
       IF (lbinanl) THEN
          ! Aerosol size distribution (regime A)
          iret = nf90_inq_varid(ncid0,'S_Naba',VarID)
          iret = nf90_put_var(ncid0,VarId,a_naerop(:,i1:i2,j1:j2,in1a:fn2a), &
               start=ibegsd,count=icntaea)
          
          !Aerosol size distribution (regime B)
          iret = nf90_inq_varid(ncid0,'S_Nabb',VarID)
          iret = nf90_put_var(ncid0,VarID,a_naerop(:,i1:i2,j1:j2,in2b:fn2b), &
               start=ibegsd,count=icntaeb)
       END IF
       
       ! Mean aerosol wet radius (regime A)
       CALL meanRadius('aerosol','a',zvar(:,:,:))
       iret = nf90_inq_varid(ncid0,'S_Rwaa',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,  &
            count=icnt)
       
       ! Mean aerosol wet radius (regime B)
       CALL meanRadius('aerosol','b',zvar(:,:,:))
       iret = nf90_inq_varid(ncid0,'S_Rwab',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,  &
            count=icnt)
       
       IF (lbinanl) THEN
          ! Aerosol bin wet radius (regime A)
          iret = nf90_inq_varid(ncid0,'S_Rwaba',VarID)
          iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in1a:fn2a),  &
               start=ibegsd,count=icntaea)
          
          ! Aerosol bin wet radius (regime B)
          iret = nf90_inq_varid(ncid0,'S_Rwabb',VarID)
          iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in2b:fn2b),  &
               start=ibegsd,count=icntaeb)
       END IF
       
       IF (level == 5) THEN
          ! Number of ice particles
          CALL bulkNumc('ice','a',zvar(:,:,:))
          zsum = zvar
          CALL bulkNumc('ice','b',zvar(:,:,:))
          zsum = zsum + zvar
          iret = nf90_inq_varid(ncid0,'S_Nic',VarID)
          iret = nf90_put_var(ncid0,VarID,zsum(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          ! Number of snow droplets
          CALL bulkNumc('snow','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_Ns',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          ! Mean ice particle wet radius (a)
          CALL meanRadius('ice','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_Rwia',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          ! Mean ice particle wet radius (b)
          CALL meanRadius('ice','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_Rwib',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          ! Mean snow drop wet radius
          CALL meanRadius('snow','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_Rwsa',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,  &
               count=icnt)
          
          IF (lbinanl) THEN
             ! Ice particle size distribution reg. a
             iret = nf90_inq_varid(ncid0,'S_Niba',VarID)
             iret = nf90_put_var(ncid0,VarID,a_nicep(:,i1:i2,j1:j2,iia%cur:fia%cur), &
                  start=ibegsd,count=icntica)
             
             ! Ice particle size distribution reg. b
             iret = nf90_inq_varid(ncid0,'S_Nibb',VarID)
             iret = nf90_put_var(ncid0,VarID,a_nicep(:,i1:i2,j1:j2,iib%cur:fib%cur), &
                  start=ibegsd,count=icntica)
             
             ! Ice particle bin wet radius regime a
             iret = nf90_inq_varid(ncid0,'S_Rwiba',VarID)
             iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,iia%cur:fia%cur), &
                  start=ibegsd,count=icntica)
             
             ! Ice particle bin wet radius regime b
             iret = nf90_inq_varid(ncid0,'S_Rwibb',VarID)
             iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,iib%cur:fib%cur), &
                  start=ibegsd,count=icntica)
             
             ! Snow size distribution
             iret = nf90_inq_varid(ncid0,'S_Nsba',VarID)
             iret = nf90_put_var(ncid0,VarID,a_nsnowp(:,i1:i2,j1:j2,isa:fsa), &
                  start=ibegsd,count=icntsna)
             
             ! Snow drop bin wet radius
             iret = nf90_inq_varid(ncid0,'S_Rwsba',VarID)
             iret = nf90_put_var(ncid0,VarID,a_Rswet(:,i1:i2,j1:j2,isa:fsa),  &
                  start=ibegsd,count=icntsna)
             
          END IF !(lbinanl)

       END IF !level 5

       ! Number of newly activated
       zvar(:,:,:) = SUM(a_nactd(:,:,:,:), DIM=4)
       iret = nf90_inq_varid(ncid0,'S_Nact',VarID)
       iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
            count=icnt)

       ! Mass mixing ratios
       IF (IsUsed(prtcl,'SO4')) THEN

          ! --Sulphate (aerosol, regime A)
          CALL bulkMixrat('SO4','aerosol','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aSO4a',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)

          ! --Sulphate (aerosol, regime B)
          CALL bulkMixrat('SO4','aerosol','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aSO4b',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          ! --Sulphate (clouds, regime A)
          CALL bulkMixrat('SO4','cloud','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cSO4a',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          ! --Sulphate (clouds, regime B)
          CALL bulkMixrat('SO4','cloud','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cSO4b',VarID)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)

          IF (level == 5) THEN
             ! --Sulphate (ice, regime A)
             CALL bulkMixrat('SO4','ice','a',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iSO4a',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
             
             ! --Sulphate (ice, regime B)
             CALL bulkMixrat('SO4','ice','b',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iSO4b',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
          END IF ! level 5

       END IF

       IF (IsUSed(prtcl,'NH')) THEN

          !-- Ammonium (aerosol, regime A)
          CALL bulkMixrat('NH','aerosol','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aNH3a',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Ammonium (aerosol, regime B)
          CALL bulkMixrat('NH','aerosol','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aNH3b',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Ammonium (clouds, regime A)
          CALL bulkMixrat('NH','cloud','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cNH3a',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Ammonium (clouds, regime B)
          CALL bulkMixrat('NH','cloud','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cNH3b',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)

          IF (level == 5) THEN
             ! --Ammonium (ice, regime A)
             CALL bulkMixrat('NH','ice','a',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iNHa',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
             
             ! --Ammonium (ice, regime B)
             CALL bulkMixrat('NH','ice','b',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iNHb',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
          END IF ! level 5

       END IF

       IF (IsUsed(prtcl,'NO')) THEN

          !-- Nitrate (aerosol, regime A)
          CALL bulkMixrat('NO','aerosol','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aNO3a',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Nitrate (aerosol, regime B)
          CALL bulkMixrat('NO','aerosol','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aNO3b',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Nitrate (clouds, regime A)
          CALL bulkMixrat('NO','cloud','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cNO3a',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Nitrate (clouds, regime B)
          CALL bulkMixrat('NO','cloud','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cNO3b',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          IF (level == 5) THEN
             ! --Nitrate (ice, regime A)
             CALL bulkMixrat('NO','ice','a',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iNOa',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
             
             ! --Nitrate (ice, regime B)
             CALL bulkMixrat('NO','ice','b',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iNOb',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
          END IF ! level 5
          
       END IF

       IF (IsUsed(prtcl,'OC')) THEN
        
          !-- Organic Carbon (aerosol, regime A)
          CALL bulkMixrat('OC','aerosol','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aOCa',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Organic Carbon (aerosol, regime B)
          CALL bulkMixrat('OC','aerosol','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aOCb',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Organic Carbon (clouds, regime A)
          CALL bulkMixrat('OC','cloud','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cOCa',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Organic Carbon (clouds, regime B)
          CALL bulkMixrat('OC','cloud','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cOCb',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          IF (level == 5) THEN
             ! --Organic Carbon (ice, regime A)
             CALL bulkMixrat('OC','ice','a',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iOCa',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
             
             ! --Organic Carbon (ice, regime B)
             CALL bulkMixrat('OC','ice','b',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iOCb',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
          END IF ! level 5

       END IF
       
       IF (IsUsed(prtcl,'BC')) THEN

          !-- Black Carbon (aerosol, regime A)
          CALL bulkMixrat('BC','aerosol','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aBCa',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Black Carbon (aerosol, regime B)
          CALL bulkMixrat('BC','aerosol','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aBCb',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Black Carbon (clouds, regime A)
          CALL bulkMixrat('BC','cloud','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cBCa',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Black Carbon (clouds, regime B)
          CALL bulkMixrat('BC','cloud','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cBCb',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)

          IF (level == 5) THEN
             ! --Black Carbon (ice, regime A)
             CALL bulkMixrat('BC','ice','a',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iBCa',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
             
             ! --Black Carbon (ice, regime B)
             CALL bulkMixrat('BC','ice','b',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iBCb',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
          END IF ! level 5

       END IF

       IF (IsUsed(prtcl,'DU')) THEN

          !-- Dust (aerosol, regime A)
          CALL bulkMixrat('DU','aerosol','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aDUa',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Dust (aerosol, regime B)
          CALL bulkMixrat('DU','aerosol','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aDUb',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Dust (clouds, regime A)
          CALL bulkMixrat('DU','cloud','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cDUa',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Dust (clouds, regime B)
          CALL bulkMixrat('DU','cloud','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cDUb',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          IF (level == 5) THEN
             ! --Dust (ice, regime A)
             CALL bulkMixrat('DU','ice','a',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iDUa',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
             
             ! --Dust (ice, regime B)
             CALL bulkMixrat('DU','ice','b',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iDUb',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
          END IF ! level 5

       END IF
       
       IF (IsUsed(prtcl,'SS')) THEN

          !-- Sea Salt (aerosol, regime A)
          CALL bulkMixrat('SS','aerosol','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aSSa',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Sea Salt (aerosol, regime B)
          CALL bulkMixrat('SS','aerosol','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_aSSb',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Sea Salt (aerosol, regime A)
          CALL bulkMixrat('SS','cloud','a',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cSSa',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          !-- Sea Salt (aerosol, regime B)
          CALL bulkMixrat('SS','cloud','b',zvar(:,:,:))
          iret = nf90_inq_varid(ncid0,'S_cSSb',VarID)
          iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
               count=icnt)
          
          IF (level == 5) THEN
             ! -- Sea Salt (ice, regime A)
             CALL bulkMixrat('SS','ice','a',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iSSa',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
             
             ! -- Sea Salt (ice, regime B)
             CALL bulkMixrat('SS','ice','b',zvar(:,:,:))
             iret = nf90_inq_varid(ncid0,'S_iSSb',VarID)
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                  count=icnt)
          END IF ! level 5
          
          END IF

       END IF
       
    END IF

    if (myid==0) print "(//' ',12('-'),'   Record ',I3,' to: ',A60)",    &
         nrec0,fname

    iret  = nf90_sync(ncid0)
    nrec0 = nrec0+1

  end subroutine write_anal
  !
  ! ----------------------------------------------------------------------
  ! Subroutine write_hist:  This subroutine writes a binary history file
  !
  subroutine write_hist(htype, time)

    use mpi_interface, only : appl_abort, myid, wrxid, wryid
    integer :: errcode=-17

    integer, intent (in) :: htype
    real, intent (in)    :: time

    character (len=80) :: hname

    integer :: n, iblank
    !
    ! create and open a new output file.
    !
    write(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
    hname = trim(hname)//'.'//trim(filprf)

    select case(htype)
    case default
       hname = trim(hname)//'.iflg'
    case(0)
       hname = trim(hname)//'.R'
    case(1)
       hname = trim(hname)//'.rst'
    case(2)
       iblank=index(hname,' ')
       write (hname(iblank:iblank+7),'(a1,i6.6,a1)') '.', int(time), 's'
    end select
    !
    ! Write fields
    !
    if (myid == 0) print "(//' ',49('-')/,' ',/,'   History write to: ',A30)" &
         ,hname
    open(10,file=trim(hname), form='unformatted')

    write(10) time,th00,umean,vmean,dtl,level,isgstyp,iradtyp,nzp,nxp,nyp,nscl
    write(10) xt, xm, yt, ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf,sst,W1,W2,W3 ! added by Zubair

    write(10) a_ustar, a_tstar, a_rstar

    write(10) a_pexnr
    write(10) a_press
    write(10) a_theta

    write(10) a_up
    write(10) a_vp
    write(10) a_wp
    write(10) a_uc
    write(10) a_vc
    write(10) a_wc

    do n=1,nscl
       call newsclr(n)
       write(10) a_sp
    end do

    if ( allocated(a_rv)   ) write(10) a_rv
    if ( allocated(a_rc)   ) write(10) a_rc
    if ( allocated(a_rflx) ) write(10) a_rflx
    close(10)

    if (myid == 0 .and. htype < 0) then
       print *, 'CFL Violation'
       call appl_abort(errcode)
    end if

    return
  end subroutine write_hist
  !
  ! ----------------------------------------------------------------------
  ! Subroutine read_hist:  This subroutine reads a binary history file
  !
  !                        Modified for level 4
  !                Juha Tonttila, FMI, 20140828
  !

  subroutine read_hist(time, hfilin)

    use mpi_interface, only : appl_abort, myid, wrxid, wryid

    character(len=80), intent(in) :: hfilin
    real, intent(out)             :: time

    character (len=80) :: hname
    integer :: n, nxpx, nypx, nzpx, nsclx, iradx, isgsx, lvlx
    logical :: exans
    real :: umx, vmx, thx
    !
    ! open input file.
    !

    write(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
    hname = trim(hname)//'.'//trim(hfilin)

    inquire(file=trim(hname),exist=exans)
    if (.not.exans) then
       print *,'ABORTING: History file', trim(hname),' not found'
       call appl_abort(0)
    else
       open (10,file=trim(hname),status='old',form='unformatted')
       read (10) time,thx,umx,vmx,dtl,lvlx,isgsx,iradx,nzpx,nxpx,nypx,nsclx

       if (nxpx /= nxp .or. nypx /= nyp .or. nzpx /= nzp)  then
          if (myid == 0) print *, nxp, nyp, nzp, nxpx, nypx, nzpx
          call appl_abort(-1)
       end if

       read (10) xt, xm, yt, ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf,sst,W1,W2,W3

       read (10) a_ustar, a_tstar, a_rstar

       read (10) a_pexnr
       read (10) a_press
       read (10) a_theta

       read (10) a_up
       read (10) a_vp
       read (10) a_wp
       read (10) a_uc
       read (10) a_vc
       read (10) a_wc

       do n=1,nscl
          call newsclr(n)
          if (n <= nsclx) read (10) a_sp
       end do
       do n=nscl+1,nsclx
          read (10)
       end do

       if (lvlx > 0 .AND. lvlx < 4) then
          if (level > 0 .AND. lvlx < 4) then
             read (10) a_rv
          else
             read (10)
          end if
       end if
       if (lvlx > 1) then
          if (level > 1) then
             read (10) a_rc
          else
             read (10)
          end if
       end if
       if (iradx > 0) then
          if (iradtyp > 0) then
             read (10) a_rflx
          else
             read (10)
          end if
       end if

       close(10)
       !
       ! adjust namelist and basic state appropriately
       !
       if (thx /= th00) then
          if (myid == 0) print "('  th00 changed  -  ',2f8.2)",th00,thx
          a_tp(:,:,:) = a_tp(:,:,:) + thx - th00
       end if
       if (umx /= umean) then
          if (myid == 0) print "('  umean changed  -  ',2f8.2)",umean,umx
          a_up = a_up + umx - umean
       end if
       if (vmx /= vmean) then
          if (myid == 0) print "('  vmean changed  -  ',2f8.2)",vmean,vmx
          a_vp = a_vp + vmx - vmean
       end if
       dtlv=2.*dtl
       dtlt=dtl

    end if

  end subroutine read_hist
  !
  ! ----------------------------------------------------------------------
  ! Subroutine newsclr:  This routine updates the scalar pointer to the
  ! value corresponding to the next scalar in the scalar table
  !
  subroutine newsclr(iscnum)

    integer, intent(in) :: iscnum

    a_sp=>a_sclrp(:,:,:,iscnum)
    a_st=>a_sclrt(:,:,:,iscnum)

    return
  end subroutine newsclr
  !
  ! -----------------------------------
  ! Subroutine bulkMixrat: Find and calculate
  ! the total mixing ratio of a given compound
  ! in aerosol particles or hydrometeors
  !
  ! Juha Tonttila, FMI, 2015
  ! Jaakko Ahola, FMI, 2015
  SUBROUTINE bulkMixrat(icomp,ipart,itype,mixrat)
    USE mo_submctl, ONLY : ncld,nbins,nprc,   &
                               ica,fca,icb,fcb,   &
                               ira,fra,           &
                               nice,nsnw,         &
                               iia,fia,iib,fib,   &
                               isa,fsa,           &
                               in1a,in2b,         &
                               fn2a,fn2b

    USE class_ComponentIndex, ONLY : GetIndex, IsUsed

    CHARACTER(len=*), INTENT(in) :: icomp  ! This should be either:
                                           ! SO4,OC,NO,NH,BC,DU,SS,H2O.

    CHARACTER(len=*), INTENT(in) :: ipart  ! This should be either:
                                           ! aerosol,cloud,rain,ice,snow
    CHARACTER(len=*), INTENT(in) :: itype  ! Select bin regime: a or b

    REAL, INTENT(out) :: mixrat(nzp,nxp,nyp)

    INTEGER :: istr,iend, mm

    mixrat = 0.

    ! Determine multipliers
    mm = GetIndex(prtcl,icomp)

    ! Given in kg/kg
    SELECT CASE(ipart)
       CASE('aerosol')
          IF (itype == 'a') THEN
             istr = (mm-1)*nbins + in1a
             iend = (mm-1)*nbins + fn2a
          ELSE IF (itype == 'b') THEN
             istr = (mm-1)*nbins + in2b
             iend = (mm-1)*nbins + fn2b
          ELSE
             STOP 'bulkMixrat: Invalid aerosol bin regime selection'
          END IF
          mixrat(:,:,:) = SUM(a_maerop(:,:,:,istr:iend),DIM=4)
       CASE('cloud')
          IF (itype == 'a') THEN
             istr = (mm-1)*ncld + ica%cur
             iend = (mm-1)*ncld + fca%cur
          ELSE IF (itype == 'b') THEN
             istr = (mm-1)*ncld + icb%cur
             iend = (mm-1)*ncld + fcb%cur
          ELSE
             STOP 'bulkMixrat: Invalid cloud bin regime selection'
          END IF
          mixrat(:,:,:) = SUM(a_mcloudp(:,:,:,istr:iend),DIM=4)
       CASE('precp')
          istr = (mm-1)*nprc + ira
          iend = (mm-1)*nprc + fra
          mixrat(:,:,:) = SUM(a_mprecpp(:,:,:,istr:iend),DIM=4)
       CASE('ice')
          IF (itype == 'a') THEN
             istr = (mm-1)*nice + iia%cur
             iend = (mm-1)*nice + fia%cur
          ELSE IF (itype == 'b') THEN
             istr = (mm-1)*nice + iib%cur
             iend = (mm-1)*nice + fib%cur
          ELSE
             STOP 'bulkMixrat: Invalid ice bin regime selection'
          END IF
          mixrat(:,:,:) = SUM(a_micep(:,:,:,istr:iend),DIM=4)
       CASE('snow')
          istr = (mm-1)*nsnw + isa
          iend = (mm-1)*nsnw + fsa
          mixrat(:,:,:) = SUM(a_msnowp(:,:,:,istr:iend),DIM=4)
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
    USE mo_submctl, ONLY : ncld, nbins, nprc, nice, nsnw
    USE class_componentIndex, ONLY : GetIndex

    CHARACTER(len=*), INTENT(in) :: icomp  ! This should be either:
                                           ! SO4,OC,NO,NH,BC,DU,SS,H2O.

    CHARACTER(len=*), INTENT(in) :: ipart  ! This should be either:
                                           ! aerosol,cloud,rain,ice,snow
    INTEGER, INTENT(in) :: ibin

    REAL, INTENT(out) :: mixr(nzp,nxp,nyp)

    INTEGER :: mm

    ! Determine multipliers
    mm = GetIndex(prtcl,icomp)

    SELECT CASE(ipart)
       CASE('aerosol')
          mixr(:,:,:) = a_maerop(:,:,:,(mm-1)*nbins+ibin)
       CASE('cloud')
          mixr(:,:,:) = a_mcloudp(:,:,:,(mm-1)*ncld+ibin)
       CASE('precp')
          mixr(:,:,:) = a_mprecpp(:,:,:,(mm-1)*nprc+ibin)
       CASE('ice')
          mixr(:,:,:) = a_micep(:,:,:,(mm-1)*nice+ibin)
       CASE('snow')
          mixr(:,:,:) = a_msnowp(:,:,:,(mm-1)*nsnw+ibin)
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
    USE mo_submctl, ONLY : ncld,nbins,nprc,nice,nsnw
    USE class_ComponentIndex, ONLY : GetNcomp
    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: ipart
    CHARACTER(len=*), INTENT(in) :: itype
    INTEGER, INTENT(in) :: ibin,ii,jj,kk
    REAL, INTENT(out) :: sumc

    INTEGER :: iend

    ! Number of components-1
    IF (itype == 'dry') THEN
        iend=GetNcomp(prtcl)-1 ! dry case
    ELSEIF (itype == 'wet') THEN
        iend=GetNcomp(prtcl) ! wet case
    ELSE
        STOP 'Error in binMixrat!'
    ENDIF

    SELECT CASE(ipart)
       CASE('aerosol')
          sumc = SUM( a_maerop(kk,ii,jj,ibin:iend*nbins+ibin:nbins) )
       CASE('cloud')
          sumc = SUM( a_mcloudp(kk,ii,jj,ibin:iend*ncld+ibin:ncld) )
       CASE('precp')
          sumc = SUM( a_mprecpp(kk,ii,jj,ibin:iend*nprc+ibin:nprc) )
       CASE('ice')
          sumc = SUM( a_micep(kk,ii,jj,ibin:iend*nice+ibin:nice) )
       CASE('snow')
          sumc = SUM( a_msnowp(kk,ii,jj,ibin:iend*nsnw+ibin:nsnw) )
       CASE DEFAULT
          STOP 'bin mixrat error'
    END SELECT

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
                               iia,fia,iib,fib,isa,fsa, &
                               in1a,in2b,fn2a,fn2b

    CHARACTER(len=*), INTENT(in) :: ipart
    CHARACTER(LEN=*), INTENT(in) :: itype
    REAL, INTENT(out) :: numc(nzp,nxp,nyp)
    INTEGER :: istr,iend

    istr = 0
    iend = 0

    ! Outputs #/kg
    ! No concentration limits (nlim or prlim) for number

    SELECT CASE(ipart)
       CASE('aerosol')
          IF (itype == 'a') THEN ! Note: 1a and 2a combined
             istr = in1a
             iend = fn2a
          ELSE IF (itype == 'b') THEN ! 2b
             istr = in2b
             iend = fn2b
          END IF
          numc(:,:,:) = SUM(a_naerop(:,:,:,istr:iend),DIM=4)
       CASE('cloud')
          IF (itype == 'a') THEN
             istr = ica%cur
             iend = fca%cur
          ELSE IF (itype == 'b') THEN
             istr = icb%cur
             iend = fcb%cur
          END IF
          numc(:,:,:) = SUM(a_ncloudp(:,:,:,istr:iend),DIM=4)
       CASE('precp')
          istr = ira
          iend = fra
          numc(:,:,:) = SUM(a_nprecpp(:,:,:,istr:iend),DIM=4)
        CASE('ice')
          IF (itype == 'a') THEN
             istr = iia%cur
             iend = fia%cur
          ELSE IF (itype == 'b') THEN
             istr = iib%cur
             iend = fib%cur
          END IF
          numc(:,:,:) = SUM(a_nicep(:,:,:,istr:iend),DIM=4)
       CASE('snow')
          istr = isa
          iend = fsa
          numc(:,:,:) = SUM(a_nsnowp(:,:,:,istr:iend),DIM=4)
    END SELECT

  END SUBROUTINE bulkNumc


  !
  ! -------------------------------------------------
  ! SUBROUTINE meanRadius
  ! Gets the mean wet radius for particles.
  !
  SUBROUTINE meanRadius(ipart,itype,rad)
    USE mo_submctl, ONLY : nbins,ncld,nprc,               &
                               nice,nsnw,                     &
                               ica,fca,icb,fcb,ira,fra,       &
                               iia,fia,iib,fib,isa,fsa,       &
                               in1a,fn2a,in2b,fn2b, &
                               nlim,prlim
    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: ipart
    CHARACTER(len=*), INTENT(in) :: itype
    REAL, INTENT(out) :: rad(nzp,nxp,nyp)

    REAL :: zvar1(nzp,nxp,nyp)

    INTEGER :: istr,iend

    rad = 0.

    ! Get the total number concentration for selected particle type
    CALL bulkNumc(ipart,itype,zvar1)

    SELECT CASE(ipart)
    CASE('aerosol')

       IF (itype == 'a') THEN ! Note: 1a and 2a combined
          istr = in1a
          iend = fn2a
       ELSE IF (itype == 'b') THEN
          istr = in2b
          iend = fn2b
       ELSE
          STOP 'meanRadius: Invalid bin regime selection (aerosol)'
       END IF

       CALL getRadius(istr,iend,nbins,a_naerop,zvar1,nlim,a_Rawet,rad)

    CASE('cloud')

       IF (itype == 'a') THEN
          istr = ica%cur
          iend = fca%cur
       ELSE IF (itype == 'b') THEN
          istr = icb%cur
          iend = fcb%cur
       ELSE
          STOP 'meanRadius: Invalid bin regime selection (cloud)'
       END IF

       CALL getRadius(istr,iend,ncld,a_ncloudp,zvar1,nlim,a_Rcwet,rad)

    CASE('precp')

       istr = ira
       iend = fra

       CALL getRadius(istr,iend,nprc,a_nprecpp,zvar1,prlim,a_Rpwet,rad)

    CASE('ice')

       IF (itype == 'a') THEN
          istr = iia%cur
          iend = fia%cur
       ELSE IF (itype == 'b') THEN
          istr = iib%cur
          iend = fib%cur
       ELSE
          STOP 'meanRadius: Invalid bin regime selection (ice)'
       END IF

       CALL getRadius(istr,iend,nice,a_nicep,zvar1,prlim,a_Riwet,rad)

    CASE('snow')

       istr = isa
       iend = fsa

       CALL getRadius(istr,iend,nsnw,a_nsnowp,zvar1,prlim,a_Rswet,rad)

    END SELECT

  END SUBROUTINE meanRadius
  !
  ! ---------------------------------------------------
  ! SUBROUTINE getRadius
  !
  SUBROUTINE getRadius(zstr,zend,nb,numc,ntot,numlim,rpart,zrad)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: nb ! Number of bins for current particle distribution
    INTEGER, INTENT(in) :: zstr,zend  ! Start and end index for averaging
    REAL, INTENT(in) :: numc(nzp,nxp,nyp,nb)
    REAL, INTENT(in) :: ntot(nzp,nxp,nyp)
    REAL, INTENT(in) :: numlim
    REAL, INTENT(in) :: rpart(nzp,nxp,nyp,nb)

    REAL, INTENT(out) :: zrad(nzp,nxp,nyp)

    LOGICAL :: nlmask(nb)
    INTEGER :: k,i,j

    DO j = 3,nyp-2
       DO i = 3,nxp-2
          DO k = 3,nzp-2
             nlmask = .FALSE.
             nlmask(zstr:zend) = ( numc(k,i,j,zstr:zend) > numlim )

             IF (ntot(k,i,j) > numlim) THEN
                zrad(k,i,j) = SUM( rpart(k,i,j,zstr:zend) * &
                                   numc(k,i,j,zstr:zend),   &
                                   MASK=nlmask(zstr:zend)   ) / &
                                   ntot(k,i,j)
             ELSE
                zrad(k,i,j) = 0.
             END IF

          END DO
       END DO
    END DO

  END SUBROUTINE getRadius

end module grid

