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

  implicit none
  !
  integer           :: nxp = 132           ! number of x points
  integer           :: nyp = 132           ! number of y points
  integer           :: nzp = 105           ! number of z points


  logical           :: nxpart = .true.     ! number of processors in x

  real              :: deltax = 35.        ! dx for basic grid
  real              :: deltay = 35.        ! dy for basic grid
  real              :: deltaz = 17.5       ! dz for basic grid
  real              :: dzrat  = 1.0        ! grid stretching ratio
  real              :: dzmax  = 1200.      ! height to start grid-stretching
  real              :: dtlong = 10.0       ! long timestep
  real              :: th00   = 288.       ! basic state temperature

  real              :: CCN = 150.e6        ! Number of CCN per kg

  LOGICAL           :: lbinanl = .FALSE.   ! Whether to write binned data to analysis files (takes a lot of space + mainly used for debugging)
  LOGICAL           :: lbinprof = .TRUE.   ! The same for profile statistics
  LOGICAL           :: lmixranl = .FALSE.  ! Whether to write species mixing ratios to analysis files
  LOGICAL           :: lmixrprof = .TRUE.  ! The same for profile statistics
  LOGICAL           :: lremts = .TRUE.     ! Save removal rates for each species (time series)
  integer           :: igrdtyp = 1         ! vertical grid type
  integer           :: isgstyp = 1         ! sgs model type
  integer           :: iradtyp = 0         ! radiation model type
  integer           :: level   = 0         ! thermodynamic level
  integer           :: naddsc  = 0         ! number of additional scalars;
  integer           :: nfpt = 10           ! number of rayleigh friction points
  real              :: distim = 300.0      ! dissipation timescale

  real              :: sst=283.   ! Surface temperature      added by Zubair Maalick
  real            :: W1 = 0.9   !Water content
  real            :: W2 = 0.9
  real            :: W3 = 0.9

  LOGICAL :: sed_aero = .TRUE.,  &
             sed_cloud = .TRUE., &
             sed_precp = .TRUE., &
             sed_ice = .TRUE., &
             sed_snow = .TRUE.

  ! Save statistics about microphysical process rates
  LOGICAL :: stat_micro = .FALSE.    ! Analysis files (*.nc)
  LOGICAL :: stat_micro_ts = .FALSE. ! Profiles (*.ps.nc)
  LOGICAL :: stat_micro_ps = .FALSE. ! Time series (*.ts.nc)

  ! Nudging options
  !   1 = soft nudging with fixed nudging constant applied only for the specified time period
  !   2 = hard nudging with the same settings
  !   else = no nudging
  INTEGER :: nudge_theta=0, & ! (liquid water) potential temperature, depending on the microphysical level
    nudge_rv=0, & ! Water vapor mixing ratio (maintain total water)
    nudge_u=0, nudge_v=0, & ! Horizontal winds
    nudge_ccn=0 ! Sectional aerosol for levels 4 and 5 (maintain aerosol+cloud+ice)
  ! Parameters related to time, altitude range and and the nudging coefficient
  REAL :: nudge_theta_time=3600., nudge_theta_zmin=-1.e10, nudge_theta_zmax=1.e10, nudge_theta_tau=300.
  REAL :: nudge_rv_time=3600., nudge_rv_zmin=-1.e10, nudge_rv_zmax=1.e10, nudge_rv_tau=300.
  REAL :: nudge_u_time=3600., nudge_u_zmin=-1.e10, nudge_u_zmax=1.e10, nudge_u_tau=300.
  REAL :: nudge_v_time=3600., nudge_v_zmin=-1.e10, nudge_v_zmax=1.e10, nudge_v_tau=300.
  REAL :: nudge_ccn_time=3600., nudge_ccn_zmin=-1.e10, nudge_ccn_zmax=1.e10, nudge_ccn_tau=300.
  real, save, allocatable :: theta_ref(:), rv_ref(:), u_ref(:), v_ref(:), aero_ref(:,:)
  LOGICAL, SAVE :: nudge_init=.TRUE.

  character (len=80):: expnme = 'Default' ! Experiment name
  character (len=80):: filprf = 'x'       ! File Prefix
  character (len=7) :: runtype = 'INITIAL'! Run Type Selection

  REAL              :: Tspinup = 7200.    ! Spinup period in seconds (added by Juha)

  ! User control of analysis outputs
  INTEGER, PARAMETER :: maxn_list=100
  CHARACTER(len=7), dimension(maxn_list), SAVE :: anl_include='       ', anl_exclude='       '
  CHARACTER(len=7), dimension(maxn_list), SAVE :: out_an_list(maxn_list)='       '
  INTEGER, SAVE :: nv4_proc=0
  REAL, SAVE, ALLOCATABLE :: out_an_data(:,:,:,:)

  character (len=80), private :: fname

  integer, private, save  ::  nrec0

  integer           :: nz, nxyzp, nxyp
  real              :: dxi, dyi, dtl, umean, vmean, psrf
  real, allocatable :: xt(:), xm(:), yt(:), ym(:), zt(:), zm(:), dzt(:), dzm(:)
  real, allocatable :: u0(:), v0(:), pi0(:), pi1(:), th0(:), dn0(:), rt0(:)
  real, allocatable :: spng_wfct(:), spng_tfct(:)
  REAL, ALLOCATABLE, target :: tmp_prcp(:,:,:,:), tmp_prct(:,:,:,:)
  REAL, ALLOCATABLE, target :: tmp_icep(:,:,:,:), tmp_icet(:,:,:,:)
  REAL, ALLOCATABLE, target :: tmp_snwp(:,:,:,:), tmp_snwt(:,:,:,:)
  REAL, ALLOCATABLE, target :: tmp_gasp(:,:,:,:), tmp_gast(:,:,:,:)
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
  ! -- Local (LES) dimensions
  INTEGER :: nbins=0,ncld=0,nice=0,nprc=0,nsnw=0,nspec=0,ngases=0

  ! No prognostic b-bins for aerosol, cloud or ice
  LOGICAL :: no_b_bins = .FALSE.
  ! No prognostic rain
  LOGICAL :: no_prog_prc = .FALSE.
  ! No prognostic ice or snow (level=5)
  LOGICAL :: no_prog_ice = .FALSE.
  LOGICAL :: no_prog_snw = .FALSE.

  real, allocatable, target :: a_sclrp(:,:,:,:),a_sclrt(:,:,:,:)
   !
  ! 3d diagnostic quantities
  !
  real, allocatable, target :: a_theta(:,:,:)  ! dry potential temp (k)
  real, allocatable :: a_pexnr(:,:,:)  ! perturbation exner func
  real, allocatable :: a_press(:,:,:)  ! pressure (hpa)
  real, allocatable :: a_rc(:,:,:)     ! Cloud (level<=3) or aerosol+cloud (level>=4) water mixing ratio
  real, allocatable :: a_ri(:,:,:)     ! Total ice water mixing ratio
  real, allocatable :: a_rv(:,:,:)     ! Water vapor mixing ratio (levels < 4)
  REAL, ALLOCATABLE :: a_srp(:,:,:)    ! Total rain water mixing ratio (levels >= 4)
  REAL, ALLOCATABLE :: a_snrp(:,:,:)   ! Total rain drop number mixing ratio (levels >=4)
  REAL, ALLOCATABLE :: a_srs(:,:,:)    ! Total snow water mixing ratio (level 5)
  REAL, ALLOCATABLE :: a_snrs(:,:,:)   ! Total snow number mixing ratio (level 5)
  REAL, ALLOCATABLE :: a_rh(:,:,:)     ! Relative humidity
  REAL, ALLOCATABLE :: a_rsl(:,:,:)    ! Water vapor saturation mixing ratio
  REAL, ALLOCATABLE :: a_rhi(:,:,:)    ! Relative humidity over ice
  REAL, ALLOCATABLE :: a_rsi(:,:,:)    ! Water vapor saturation mixing ratio over ice
  REAL, ALLOCATABLE :: a_dn(:,:,:)     ! Air density
  !
  ! scratch arrays
  !
  real, allocatable, dimension (:,:,:) :: a_rflx, a_sflx, &
       a_fus, a_fds, a_fuir, a_fdir, &
       a_temp
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
  real, allocatable :: obl(:,:)
  real, allocatable :: aerin(:,:,:), cldin(:,:,:), precip(:,:,:), icein(:,:,:), snowin(:,:,:), albedo(:,:)

  ! Statistics for level 3, 4 and 5 microphysics
  REAL, dimension(:,:,:), allocatable :: &
            coag_ra, coag_na, coag_rc, coag_nc, coag_rr, coag_nr, &
            coag_ri, coag_ni, coag_rs, coag_ns, &
            cond_ra, cond_rc, cond_rr, cond_ri, cond_rs, &
            auto_rr, auto_nr, auto_rs, auto_ns, &
            cact_rc, cact_nc, nucl_ri, nucl_ni, &
            melt_ri, melt_ni, melt_rs, melt_ns, &
            sedi_ra, sedi_na, sedi_rc, sedi_nc, sedi_rr, sedi_nr, &
            sedi_ri, sedi_ni, sedi_rs, sedi_ns, &
            diag_ra, diag_na, diag_rc, diag_nc, diag_rr, diag_nr, &
            diag_ri, diag_ni, diag_rs, diag_ns, cond_nr

  !
  integer :: nscl = 1
  integer, save :: ncid0,ncid_s
  !
contains
  !
  ! Copy SALSA parameters to LES
  SUBROUTINE copy_salsa_pars(nspec_out,nprc_out,nsnw_out,ngases_out)
    USE mo_submctl, ONLY : nspec,nprc,nsnw,ngases
    INTEGER, INTENT(OUT) :: nspec_out, nprc_out, nsnw_out, ngases_out
    nspec_out=nspec
    nprc_out=nprc
    nsnw_out=nsnw
    ngases_out=ngases
  END SUBROUTINE copy_salsa_pars
  !
  !----------------------------------------------------------------------
  ! SUBROUTINE define_vars
  !
  ! Modified for level 4
  ! Juha Tonttila, FMI, 2014.
  !
  subroutine define_vars

    use mpi_interface, only :myid
    USE mo_submctl, ONLY : fn2a,fn2b,fnp2a,fnp2b ! Indexes for SALSA

    integer :: memsize
    INTEGER :: zz
    INTEGER :: nc
    INTEGER :: nsalsa

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
       allocate (a_fus(nzp,nxp,nyp),a_fds(nzp,nxp,nyp),a_fuir(nzp,nxp,nyp),a_fdir(nzp,nxp,nyp))
       a_fus(:,:,:) = 0.
       a_fds(:,:,:) = 0.
       a_fuir(:,:,:) = 0.
       a_fdir(:,:,:) = 0.
       memsize = memsize + 5*nxyzp + nxyp
    end if

    allocate (a_temp(nzp,nxp,nyp),a_rsl(nzp,nxp,nyp))
    a_temp(:,:,:) = 0.
    a_rsl(:,:,:) = 0.
    memsize = memsize + nxyzp*2

    ! Juha: Stuff that's allocated if SALSA is NOT used
    !-----------------------------------------------------
    IF (level < 4) THEN

       if (level < 1) then
            WRITE(*,*) 'Level < 1 not accepted!'
            STOP
       end if

       allocate (a_rv(nzp,nxp,nyp),a_rc(nzp,nxp,nyp),a_dn(nzp,nxp,nyp))
       a_rv(:,:,:) = 0.
       a_rc(:,:,:) = 0.
       a_dn(:,:,:) = 0.
       memsize = memsize + 3*nxyzp

       ! Prognostic scalars: temperature + total water + rain mass and number (level=3) + tke (isgstyp> 1) + additional scalars
       nscl = 2+naddsc
       if (level == 3) nscl = nscl+2 ! rain
       if (isgstyp > 1) nscl = nscl+1 ! tke

       allocate (a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
       a_sclrp(:,:,:,:) = 0.
       a_sclrt(:,:,:,:) = 0.
       memsize = memsize + 2*nscl*nxyzp

       a_tp=>a_sclrp(:,:,:,1)
       a_tt=>a_sclrt(:,:,:,1)
       a_rp=>a_sclrp(:,:,:,2)
       a_rt=>a_sclrt(:,:,:,2)
       if (level == 3) then
          a_rpp=>a_sclrp(:,:,:,3)
          a_rpt=>a_sclrt(:,:,:,3)
          a_npp=>a_sclrp(:,:,:,4)
          a_npt=>a_sclrt(:,:,:,4)
       end if
       if (isgstyp > 1) then
          a_qp=>a_sclrp(:,:,:,nscl - naddsc)
          a_qt=>a_sclrt(:,:,:,nscl - naddsc)
       end if

       IF (level == 3) THEN
          ! Allocate arrays for level 3 process rate statistics
          allocate ( coag_rr(nzp,nxp,nyp), coag_nr(nzp,nxp,nyp), &
             cond_rr(nzp,nxp,nyp), cond_nr(nzp,nxp,nyp), &
             auto_rr(nzp,nxp,nyp), auto_nr(nzp,nxp,nyp), &
             sedi_rc(nzp,nxp,nyp), sedi_rr(nzp,nxp,nyp), sedi_nr(nzp,nxp,nyp), &
             diag_rr(nzp,nxp,nyp), diag_nr(nzp,nxp,nyp) )
          coag_rr=0.; coag_nr=0.; cond_rr=0.; cond_nr=0.
          auto_rr=0.; auto_nr=0.; sedi_rc=0.; sedi_rr=0.; sedi_nr=0.
          diag_rr=0.; diag_nr=0.
          memsize = memsize + nxyzp*11
       END IF

    !Juha: Stuff that's allocated when SALSA is used
    !---------------------------------------------------
    ELSE IF (level >= 4) THEN

       allocate (a_rc(nzp,nxp,nyp), a_srp(nzp,nxp,nyp), a_snrp(nzp,nxp,nyp),     &
                 a_rh(nzp,nxp,nyp),a_dn(nzp,nxp,nyp)  )

       a_rc(:,:,:) = 0.
       a_srp(:,:,:) = 0.
       a_snrp(:,:,:) = 0.
       a_rh(:,:,:) = 0.
       a_dn(:,:,:) = 0.
       memsize = memsize + 5*nxyzp

       allocate (a_ri(nzp,nxp,nyp), a_rsi(nzp,nxp,nyp), a_rhi(nzp,nxp,nyp),      &
                  a_srs(nzp,nxp,nyp), a_snrs(nzp,nxp,nyp)  )
       a_ri(:,:,:) = 0.
       a_rsi(:,:,:) = 0.
       a_rhi(:,:,:) = 0.
       a_srs(:,:,:) = 0.
       a_snrs(:,:,:) = 0.
       memsize = memsize + 5*nxyzp

       ! Number of prognostic SALSA variables
       IF (no_b_bins) THEN
          ! No b-bins for LES, but full arrays are needed for SALSA
          nbins=fn2a
          ncld=fnp2a
          nice=fnp2a
       ELSE
          nbins=fn2b
          ncld=fnp2b
          nice=fnp2b
       ENDIF
       CALL copy_salsa_pars(nspec,nprc,nsnw,ngases)

       ! Total number of prognostic SALSA variables (number and mass for each aerosol component + gases)
       nc = nspec+2
       nsalsa = ngases + nc*nbins + nc*ncld
       IF (.NOT. no_prog_prc) nsalsa = nsalsa + nc*nprc
       IF (level>=5 .AND. .NOT. no_prog_ice) nsalsa = nsalsa + nc*nice
       IF (level>=5 .AND. .NOT. no_prog_snw) nsalsa = nsalsa + nc*nsnw

       ! Total number of prognostic scalars: temperature + water vapor + SALSA + tke (isgstyp > 1) [+ additional scalars]
       nscl = 2 + nsalsa + naddsc
       if (isgstyp > 1) nscl = nscl+1

       allocate (a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
       a_sclrp(:,:,:,:) = 0.
       a_sclrt(:,:,:,:) = 0.
       memsize = memsize + 2*nscl*nxyzp

       a_tp=>a_sclrp(:,:,:,1)
       a_tt=>a_sclrt(:,:,:,1)
       a_rp=>a_sclrp(:,:,:,2)
       a_rt=>a_sclrt(:,:,:,2)

       if (isgstyp > 1) then
          a_qp=>a_sclrp(:,:,:,nscl - naddsc)
          a_qt=>a_sclrt(:,:,:,nscl - naddsc)
       end if

       !JT: Set the pointers for prognostic SALSA variables (levels 4 & 5)
       nc = nspec+1
       zz = 2
       a_naerop => a_sclrp(:,:,:,zz+1:zz+nbins)
       a_naerot => a_sclrt(:,:,:,zz+1:zz+nbins)
       zz = zz+nbins
       a_maerop => a_sclrp(:,:,:,zz+1:zz+nc*nbins)
       a_maerot => a_sclrt(:,:,:,zz+1:zz+nc*nbins)
       zz = zz+nc*nbins

       a_ncloudp => a_sclrp(:,:,:,zz+1:zz+ncld)
       a_ncloudt => a_sclrt(:,:,:,zz+1:zz+ncld)
       zz = zz+ncld
       a_mcloudp => a_sclrp(:,:,:,zz+1:zz+nc*ncld)
       a_mcloudt => a_sclrt(:,:,:,zz+1:zz+nc*ncld)
       zz = zz+nc*ncld

       IF (.NOT. no_prog_prc) THEN
          ! Prognostic rain
          a_nprecpp => a_sclrp(:,:,:,zz+1:zz+nprc)
          a_nprecpt => a_sclrt(:,:,:,zz+1:zz+nprc)
          zz = zz+nprc
          a_mprecpp => a_sclrp(:,:,:,zz+1:zz+nc*nprc)
          a_mprecpt => a_sclrt(:,:,:,zz+1:zz+nc*nprc)
          zz = zz+nc*nprc
       ELSE
          ! Allocate zero arrays for pointers
          ALLOCATE (tmp_prcp(nzp,nxp,nyp,(nc+1)*nprc), &
                    tmp_prct(nzp,nxp,nyp,(nc+1)*nprc))
          tmp_prcp(:,:,:,:) = 0.
          tmp_prct(:,:,:,:) = 0.
          a_nprecpp => tmp_prcp(:,:,:,1:nprc)
          a_nprecpt => tmp_prct(:,:,:,1:nprc)
          a_mprecpp => tmp_prcp(:,:,:,nprc+1:(nc+1)*nprc)
          a_mprecpt => tmp_prct(:,:,:,nprc+1:(nc+1)*nprc)
       ENDIF

       IF (level>=5 .AND. .NOT.no_prog_ice) THEN
          ! Prognostic ice
          a_nicep => a_sclrp(:,:,:,zz+1:zz+nice)
          a_nicet => a_sclrt(:,:,:,zz+1:zz+nice)
          zz = zz+nice
          a_micep => a_sclrp(:,:,:,zz+1:zz+nc*nice)
          a_micet => a_sclrt(:,:,:,zz+1:zz+nc*nice)
          zz = zz+nc*nice
       ELSE
          ! Allocate zero arrays for pointers
          ALLOCATE (tmp_icep(nzp,nxp,nyp,(nc+1)*nice), &
                    tmp_icet(nzp,nxp,nyp,(nc+1)*nice))
          tmp_icep(:,:,:,:) = 0.
          tmp_icet(:,:,:,:) = 0.
          a_nicep => tmp_icep(:,:,:,1:nice)
          a_nicet => tmp_icet(:,:,:,1:nice)
          a_micep => tmp_icep(:,:,:,nice+1:(nc+1)*nice)
          a_micet => tmp_icet(:,:,:,nice+1:(nc+1)*nice)
       ENDIF

       IF (level>=5 .AND. .NOT.no_prog_snw) THEN
          ! Prognostic snow
          a_nsnowp => a_sclrp(:,:,:,zz+1:zz+nsnw)
          a_nsnowt => a_sclrt(:,:,:,zz+1:zz+nsnw)
          zz = zz+nsnw
          a_msnowp => a_sclrp(:,:,:,zz+1:zz+nc*nsnw)
          a_msnowt => a_sclrt(:,:,:,zz+1:zz+nc*nsnw)
          zz = zz+nc*nsnw
       ELSE
          ! Allocate zero arrays for pointers
          ALLOCATE (tmp_snwp(nzp,nxp,nyp,(nc+1)*nsnw), &
                    tmp_snwt(nzp,nxp,nyp,(nc+1)*nsnw))
          tmp_snwp(:,:,:,:) = 0.
          tmp_snwt(:,:,:,:) = 0.
          a_nsnowp => tmp_snwp(:,:,:,1:nsnw)
          a_nsnowt => tmp_snwt(:,:,:,1:nsnw)
          a_msnowp => tmp_snwp(:,:,:,nsnw+1:(nc+1)*nsnw)
          a_msnowt => tmp_snwt(:,:,:,nsnw+1:(nc+1)*nsnw)
       ENDIF

       IF (ngases>0) THEN
          a_gaerop => a_sclrp(:,:,:,zz+1:zz+ngases)
          a_gaerot => a_sclrt(:,:,:,zz+1:zz+ngases)
          zz = zz+ngases
       ELSE
          ! Gases not included so allocate zero arrays for pointers
          ALLOCATE (tmp_gasp(nzp,nxp,nyp,0),tmp_gast(nzp,nxp,nyp,0))
          tmp_gasp = 0.
          tmp_gast = 0.
          a_gaerop => tmp_gasp(:,:,:,:)
          a_gaerot => tmp_gast(:,:,:,:)
       ENDIF

        ! Allocate arrays for level 4 and 5 process rate statistics
        allocate ( coag_ra(nzp,nxp,nyp), coag_na(nzp,nxp,nyp), coag_rc(nzp,nxp,nyp), coag_nc(nzp,nxp,nyp), &
          coag_rr(nzp,nxp,nyp), coag_nr(nzp,nxp,nyp), &
          coag_ri(nzp,nxp,nyp), coag_ni(nzp,nxp,nyp), coag_rs(nzp,nxp,nyp), coag_ns(nzp,nxp,nyp), &
          cond_ra(nzp,nxp,nyp), cond_rc(nzp,nxp,nyp), cond_rr(nzp,nxp,nyp), cond_ri(nzp,nxp,nyp), cond_rs(nzp,nxp,nyp), &
          auto_rr(nzp,nxp,nyp), auto_nr(nzp,nxp,nyp), auto_rs(nzp,nxp,nyp), auto_ns(nzp,nxp,nyp), &
          cact_rc(nzp,nxp,nyp), cact_nc(nzp,nxp,nyp), nucl_ri(nzp,nxp,nyp), nucl_ni(nzp,nxp,nyp), &
          melt_ri(nzp,nxp,nyp), melt_ni(nzp,nxp,nyp), melt_rs(nzp,nxp,nyp), melt_ns(nzp,nxp,nyp), &
          sedi_ra(nzp,nxp,nyp), sedi_na(nzp,nxp,nyp), sedi_rc(nzp,nxp,nyp), sedi_nc(nzp,nxp,nyp), &
          sedi_rr(nzp,nxp,nyp), sedi_nr(nzp,nxp,nyp), &
          sedi_ri(nzp,nxp,nyp), sedi_ni(nzp,nxp,nyp), sedi_rs(nzp,nxp,nyp), sedi_ns(nzp,nxp,nyp), &
          diag_ra(nzp,nxp,nyp), diag_na(nzp,nxp,nyp), diag_rc(nzp,nxp,nyp), diag_nc(nzp,nxp,nyp), &
          diag_rr(nzp,nxp,nyp), diag_nr(nzp,nxp,nyp), &
          diag_ri(nzp,nxp,nyp), diag_ni(nzp,nxp,nyp), diag_rs(nzp,nxp,nyp), diag_ns(nzp,nxp,nyp) )
        coag_ra=0.; coag_na=0.; coag_rc=0.;coag_nc=0.; coag_rr=0.; coag_nr=0.
        coag_ri=0.; coag_ni=0.; coag_rs=0.; coag_ns=0.
        cond_ra=0.; cond_rc=0.; cond_rr=0.; cond_ri=0.; cond_rs=0.
        auto_rr=0.; auto_nr=0.; auto_rs=0.; auto_ns=0.
        cact_rc=0.; cact_nc=0.; nucl_ri=0.; nucl_ni=0.
        melt_ri=0.; melt_ni=0.; melt_rs=0.; melt_ns=0.
        sedi_ra=0.; sedi_na=0.; sedi_rc=0.; sedi_nc=0.; sedi_rr=0.; sedi_nr=0.
        sedi_ri=0.; sedi_ni=0.; sedi_rs=0.; sedi_ns=0.
        diag_ra=0.; diag_na=0.; diag_rc=0.; diag_nc=0.; diag_rr=0.; diag_nr=0.
        diag_ri=0.; diag_ni=0.; diag_rs=0.; diag_ns=0.
        memsize = memsize + nxyzp*47

    END IF ! level

    !----------------------------------------------------

    allocate (a_ustar(nxp,nyp),a_tstar(nxp,nyp),a_rstar(nxp,nyp))
    allocate (uw_sfc(nxp,nyp),vw_sfc(nxp,nyp),ww_sfc(nxp,nyp))
    allocate (wt_sfc(nxp,nyp),wq_sfc(nxp,nyp))
    allocate (obl(nxp,nyp))

    allocate(precip(nzp,nxp,nyp))
    precip = 0.
    memsize = memsize + nxyzp

    if (level >= 4) then
       allocate(aerin(nzp,nxp,nyp),cldin(nzp,nxp,nyp),icein(nzp,nxp,nyp),snowin(nzp,nxp,nyp))
       aerin = 0.
       cldin = 0.
       icein = 0.
       snowin = 0.
       memsize = memsize + nxyzp*4
    end if

    a_ustar(:,:) = 0.
    a_tstar(:,:) = 0.
    a_rstar(:,:) = 0.
    uw_sfc(:,:)  = 0.
    vw_sfc(:,:)  = 0.
    ww_sfc(:,:)  = 0.
    wt_sfc(:,:) = 0.
    wq_sfc(:,:) = 0.
    obl(:,:) = 0.
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
    dzmin = 0.
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
     do k=1,nzp-2
        zm(k+1) = cos( ((2.*nchby - 1. - 2.*(k-1))*2.*asin(1.))/(2.*nchby))
        zm(k+1) = (zm(k+1)+1.)*dzmax/2.
     end do
     zm(nzp-1) = dzmax
     zm(nzp)   = dzmax + (zm(nzp-1)-zm(nzp-2))
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
    USE mo_submctl, ONLY : fn1a,fn2a,fn2b, &
                nlcoag,nlcnd,nlauto,nlautosnow,nlactiv,nlicenucl,nlicmelt,&
                stat_b_bins,ice_target_opt,zspec,zgas
    IMPLICIT NONE
    real, intent (in) :: time
    ! Dimensions (time, x, y, x, and SALSA bins) and constants (u0, v0, dn0) are saved
    ! during initialization, and common variables (u, v, w, theta, p) are always saved.
    INTEGER, PARAMETER :: n_dims=14, n_base=13
    character(len=7) :: s_dims(n_dims) = (/ &
         'time   ','zt     ','zm     ','xt     ','xm     ','yt     ','ym     ', & ! 1-7
         'u0     ','v0     ','dn0    ','P_Rd12a','P_Rd2ab','P_Rwprc','P_Rwsnw'/)  ! 8-14
    character(len=7) :: s_base(n_base) = (/ &
         'u      ','v      ','w      ','theta  ','p      ','stke   ','rflx   ', & ! 1-7
         'q      ','l      ','r      ','n      ','i      ','s      '/)            ! 8-13
    LOGICAL, SAVE :: b_dims(n_dims)=.TRUE., b_base(n_base)=.TRUE.
    character (len=7), allocatable :: sanal(:)

    ! Process rates
    INTEGER, PARAMETER :: n_lvl3_rate = 11, n_salsa_rate = 47
    character(len=7) :: s_lvl3_rate(n_lvl3_rate) = (/ &
        'coag_rr','coag_nr','cond_rr','cond_nr','auto_rr','auto_nr', & ! 1-6
        'sedi_rc','sedi_rr','sedi_nr','diag_rr','diag_nr'/)   ! 7-11
    character(len=7) :: s_salsa_rate(n_salsa_rate) = (/ &
        'coag_ra','coag_na','coag_rc','coag_nc','coag_rr','coag_nr', & ! 1-6
        'cond_ra','cond_rc','cond_rr','auto_rr','auto_nr','act_rc ','act_nc ', & ! 7-13
        'sedi_ra','sedi_na','sedi_rc','sedi_nc','sedi_rr','sedi_nr', & ! 14-19
        'diag_ra','diag_na','diag_rc','diag_nc','diag_rr','diag_nr', & ! 20-25
        'coag_ri','coag_ni','coag_rs','coag_ns', & ! 26-29
        'cond_ri','cond_rs','auto_rs','auto_ns','nucl_ri','nucl_ni', & ! 30-35
        'sedi_ri','sedi_ni','sedi_rs','sedi_ns', & !36-39
        'diag_ri','diag_ni','diag_rs','diag_ns', & ! 40-43
        'melt_ri','melt_ni','melt_rs','melt_ns'/)  ! 44-47
    LOGICAL :: b_lvl3_rate(n_lvl3_rate), b_salsa_rate(n_salsa_rate)

    ! SALSA
    INTEGER, PARAMETER :: n_salsa = 16
    character(len=7) :: s_salsa(n_salsa) = (/ &
         'S_Naa  ','S_Rwaa ','S_Nab  ','S_Rwab ', & ! 1-4
         'S_Nca  ','S_Rwca ','S_Ncb  ','S_Rwcb ', & ! 5-8
         'S_Np   ','S_Rwp  ', & ! 9-10
         'S_Nia  ','S_Rwia ','S_Nib  ','S_Rwib ', & ! 11-14
         'S_Ns   ','S_Rws  '/) ! 15-16
    LOGICAL :: b_salsa(n_salsa) = .TRUE.

    ! SALSA species and bins
    CHARACTER (len=7), SAVE, ALLOCATABLE :: s_bin(:), s_mixr(:), stot(:)
    LOGICAL, SAVE, ALLOCATABLE :: b_bin(:), b_mixr(:), btot(:)

    ! Local variables
    INTEGER :: i, ii, e, ee, n, n_bin, n_mixr, nvar0
    character (len=7) :: v_snm='sxx    '
    LOGICAL :: found

    b_base(6) = (isgstyp > 1)
    b_base(7) = (iradtyp > 1)


    ! Allocate data for user selected process rate outputs (see init_stat in stat.f90)
    ALLOCATE ( out_an_data(nzp,nxp,nyp,nv4_proc) )
    out_an_data(:,:,:,:) = 0.

    IF (level < 4) THEN  ! Standard operation for levels 1-3
        b_dims(11:14) = .FALSE. ! SALSA bins
        b_base(10:11) = (level==3) ! Rain
        b_base(12:13) = .FALSE. ! Ice and snow
        b_lvl3_rate(:) = (level == 3 .AND. stat_micro) ! Process rate statistics

       ! Merge logical and name arrays
       i=n_dims+n_base+n_lvl3_rate+nv4_proc+naddsc
       ALLOCATE( btot(i), stot(i) )
       i=1; e=n_dims
       btot(i:e)=b_dims; stot(i:e)=s_dims
       i=e+1; e=e+n_base
       btot(i:e)=b_base; stot(i:e)=s_base
       i=e+1; e=e+n_lvl3_rate
       btot(i:e)=b_lvl3_rate; stot(i:e)=s_lvl3_rate
       IF (nv4_proc>0) THEN
          i=e+1; e=e+nv4_proc
          btot(i:e)=.TRUE.; stot(i:e)=out_an_list(1:nv4_proc)
       END IF
    ELSE IF (level >= 4) THEN ! Operation with SALSA
       ! Additional arrays for SALSA
       ! -- dimensions
       n_mixr = 5*(nspec+1)+ngases ! Total mixing ratios for aerosol, cloud, rain, ice and snow species, and gas concentrations
       n_bin  = 8*(nspec+1+2) ! Species mixing ratios for aerosol a/b, cloud a/b, rain, ice a/b and snow; also number and size
       ! -- allocate
       ALLOCATE ( s_mixr(n_mixr),s_bin(n_bin) )
       ALLOCATE ( b_mixr(n_mixr),b_bin(n_bin) )

       ! Create a boolean array for items that are actually used

       ! 1) Standard SALSA outputs
       b_dims(11:14) = lbinanl ! SALSA bins
       b_base(11)=.FALSE. ! Rain number
       b_base(12:13) = (level>4) ! Ice and snow
       b_salsa(11:16) = (level>4) ! Ice and snow

       ! b-bins are not always saved
       IF (.not. stat_b_bins) b_salsa((/3,4,7,8,13,14/)) = .FALSE.

       ! 2) Microphysicsal process rate statistics
       b_salsa_rate(:) = .FALSE.
       IF (stat_micro) THEN
          b_salsa_rate(:) = .FALSE.
          b_salsa_rate(1:6) = nlcoag    ! Coagulation
          b_salsa_rate(7:9) = nlcnd     ! Condensation
          b_salsa_rate(10:11) = nlauto  ! Autoconversion
          b_salsa_rate(12:13) = nlactiv ! Cloud activation
          b_salsa_rate(14:15) = sed_aero  ! Aerosol sedimentation
          b_salsa_rate(16:17) = sed_cloud ! Cloud sedimentation
          b_salsa_rate(18:19) = sed_precp ! Rain sedimentation
          b_salsa_rate(20:25) = .TRUE.  ! SALSA_diagnostics
          IF (level>=5) THEN
            b_salsa_rate(26:29) = nlcoag
            b_salsa_rate(30:31) = nlcnd
            b_salsa_rate(32:33) = nlautosnow .OR. (nlicenucl .AND. ice_target_opt>=0) ! Autoconversion or modelled ice nucleation
            b_salsa_rate(34:35) = nlicenucl
            b_salsa_rate(36:37) = sed_ice
            b_salsa_rate(38:39) = sed_snow
            b_salsa_rate(40:43) = .TRUE.   ! SALSA_diagnostics
            b_salsa_rate(44:47) = nlicmelt ! Melting
          ENDIF
       ENDIF

       ! 3) Species and bin dependent outputs - generate name and logical arrays
       b_mixr(:)=lmixranl
       b_bin(:)=lbinanl
       ! Bin dependent outputs
       s_bin(1:16)=(/'S_Naba ','S_Nabb ','S_Rwaba','S_Rwabb', & ! 1-4
            'S_Ncba ','S_Ncbb ','S_Rwcba','S_Rwcbb','S_Npb  ','S_Rwpb ', & ! 5-10
            'S_Niba ','S_Nibb ','S_Rwiba','S_Rwibb','S_Nsb  ','S_Rwsb '/) ! 11-16
       IF (.not.stat_b_bins) b_bin((/2,4,6,8,12,14/))=.FALSE.
       ! Species dependent outputs
       DO ee=1,nspec+1 ! Aerosol species and water
          ! Total mixing ratios
          ii = (ee-1)*5
          s_mixr(ii+1)='S_c'//TRIM(zspec(ee))//'a'
          s_mixr(ii+2)='S_c'//TRIM(zspec(ee))//'c'
          s_mixr(ii+3)='S_c'//TRIM(zspec(ee))//'p'
          s_mixr(ii+4)='S_c'//TRIM(zspec(ee))//'i'
          s_mixr(ii+5)='S_c'//TRIM(zspec(ee))//'s'
          b_mixr(ii+4:ii+5) = lmixranl .AND. (level>4)
          ! Mixing ratios for each bin
          ii=(ee-1)*8+16
          s_bin(ii+1)='S_'//TRIM(zspec(ee))//'aa'
          s_bin(ii+2)='S_'//TRIM(zspec(ee))//'ab'
          s_bin(ii+3)='S_'//TRIM(zspec(ee))//'ca'
          s_bin(ii+4)='S_'//TRIM(zspec(ee))//'cb'
          s_bin(ii+5)='S_'//TRIM(zspec(ee))//'pb'
          s_bin(ii+6)='S_'//TRIM(zspec(ee))//'ia'
          s_bin(ii+7)='S_'//TRIM(zspec(ee))//'ib'
          s_bin(ii+8)='S_'//TRIM(zspec(ee))//'sb'
          b_bin(ii+6:ii+8) = lbinanl .AND. (level>4)
          IF (.not.stat_b_bins) b_bin((/ii+2,ii+4,ii+7/)) = .FALSE.
       ENDDO

       ! 4) Gas mixing ratios
       DO ee=1,ngases
          ! Total mixing ratios
          ii = 5*(nspec+1)+ee
          s_mixr(ii)='S_c'//TRIM(zgas(ee))//'g'
          b_mixr(ii) = .TRUE.
       ENDDO

       ! Merge logical and name arrays
       i=n_dims+n_base+n_salsa_rate+n_salsa+n_mixr+n_bin+naddsc+nv4_proc
       ALLOCATE( btot(i), stot(i) )
       i=1; e=n_dims
       btot(i:e)=b_dims; stot(i:e)=s_dims
       i=e+1; e=e+n_base
       btot(i:e)=b_base; stot(i:e)=s_base
       i=e+1; e=e+n_salsa_rate
       btot(i:e)=b_salsa_rate; stot(i:e)=s_salsa_rate
       i=e+1; e=e+n_salsa
       btot(i:e)=b_salsa; stot(i:e)=s_salsa
       i=e+1; e=e+n_mixr
       btot(i:e)=b_mixr; stot(i:e)=s_mixr
       i=e+1; e=e+n_bin
       btot(i:e)=b_bin; stot(i:e)=s_bin
       IF (nv4_proc>0) THEN
          i=e+1; e=e+nv4_proc
          btot(i:e)=.TRUE.; stot(i:e)=out_an_list(1:nv4_proc)
       END IF
    END IF

    ! Addtional scalars
    IF (naddsc>0) THEN
        i=e+1; e=e+naddsc
        btot(i:e)=.TRUE.
        do i=1,naddsc
            write(v_snm(2:3),'(i2.2)') i
            stot(e-naddsc+i) = v_snm
        end do
    ENDIF

    ! Include and/or exclude outputs based on user inputs
    n=e ! Length of btot and stot
    DO i=1,maxn_list
        IF (LEN_TRIM(anl_include(i))>0) THEN
            found=.FALSE.
            DO e=1,n
                IF ( anl_include(i)==stot(e) ) THEN
                    btot(e)=.TRUE.
                    found=.TRUE.
                ENDIF
            ENDDO
            IF (.not.found) WRITE(*,*) 'Warning: can not save '//TRIM(anl_include(i))//'!',i
        ENDIF
        IF (LEN_TRIM(anl_exclude(i))>0) THEN
            found=.FALSE.
            DO e=1,n
                IF ( anl_exclude(i)==stot(e) ) THEN
                    btot(e)=.FALSE.
                    found=.TRUE.
                ENDIF
            ENDDO
            IF (.not.found) WRITE(*,*) 'Warning: can not exclude '//TRIM(anl_exclude(i))//'!',i
        ENDIF
    ENDDO


    nvar0 = COUNT(btot)
    ALLOCATE(sanal(nvar0))
    sanal = PACK(stot,btot)

    fname =  trim(filprf)
    if(myid == 0) print                                                  &
            "(//' ',49('-')/,' ',/,'   Initializing: ',A20,'  N=',I3)",trim(fname),COUNT(btot)
    call open_nc( fname, expnme, time, (nxp-4)*(nyp-4), ncid0, nrec0, ver, author, info)

    IF (level < 4 .OR. .NOT. lbinanl) THEN
       call define_nc( ncid0, nrec0, nvar0, sanal, n1=nzp, n2=nxp-4, n3=nyp-4)
    ELSE IF (lbinanl) THEN
       call define_nc( ncid0, nrec0, nvar0, sanal, n1=nzp, n2=nxp-4, n3=nyp-4,  &
                       n1a=fn1a,n2a=fn2a-fn1a, n2b=fn2b-fn2a, nprc=nprc, nsnw=nsnw )
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
    USE mo_submctl, ONLY : in1a,in2a,fn2a,in2b,fn2b, &
                               inp2a,fnp2a,inp2b,fnp2b, &
                               aerobins, precpbins, snowbins, &
                               nlim, prlim, &
                               zspec, zgas

    real, intent (in) :: time

    integer :: iret, VarID, bb, ee
    integer :: ibeg(4), icnt(4), i1, i2, j1, j2
    INTEGER :: ibegsd(5), icntaea(5), icntaeb(5), icntpr(5), icntsn(5)
    REAL :: zvar(nzp,nxp,nyp)
    REAL :: a_Rawet(nzp,nxp,nyp,nbins), a_Rcwet(nzp,nxp,nyp,ncld),a_Rpwet(nzp,nxp,nyp,nprc), &
          a_Riwet(nzp,nxp,nyp,nice),a_Rswet(nzp,nxp,nyp,nsnw)

    icnt = (/nzp, nxp-4, nyp-4, 1/)
    icntaea = (/nzp, nxp-4, nyp-4, fn2a, 1 /)
    icntaeb = (/nzp, nxp-4, nyp-4, fn2b-fn2a, 1/)
    icntpr = (/nzp, nxp-4, nyp-4, nprc, 1/)
    icntsn = (/nzp, nxp-4, nyp-4, nsnw, 1/)
    ibeg = (/1, 1, 1, nrec0/)
    ibegsd = (/1, 1, 1, 1, nrec0/)

    i1 = 3
    i2 = nxp-2
    j1 = 3
    j2 = nyp-2

    iret = nf90_inq_varid(ncid0, 'time', VarID)
    iret = nf90_put_var(ncid0, VarID, time, start=(/nrec0/))

    ! Dimensions and constants
    if (nrec0 == 1) then
       iret = nf90_inq_varid(ncid0, 'zt', VarID)
       iret = nf90_put_var(ncid0, VarID, zt, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, 'zm', VarID)
       iret = nf90_put_var(ncid0, VarID, zm, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, 'xt', VarID)
       iret = nf90_put_var(ncid0, VarID, xt(i1:i2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, 'xm', VarID)
       iret = nf90_put_var(ncid0, VarID, xm(i1:i2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, 'yt', VarID)
       iret = nf90_put_var(ncid0, VarID, yt(j1:j2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, 'ym', VarID)
       iret = nf90_put_var(ncid0, VarID, ym(j1:j2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, 'u0', VarID)
       iret = nf90_put_var(ncid0, VarID, u0, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, 'v0', VarID)
       iret = nf90_put_var(ncid0, VarID, v0, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, 'dn0', VarID)
       iret = nf90_put_var(ncid0, VarID, dn0, start = (/nrec0/))

       IF (level >= 4 .AND. lbinanl) THEN
          iret = nf90_inq_varid(ncid0,'S_Rd12a', VarID)
          iret = nf90_put_var(ncid0, VarID, aerobins(in1a:fn2a), start = (/nrec0/))

          iret = nf90_inq_varid(ncid0,'S_Rd2ab', VarID)
          iret = nf90_put_var(ncid0,VarID, aerobins(in2a:fn2a), start = (/nrec0/))

          iret = nf90_inq_varid(ncid0,'S_Rwprc', VarID)
          iret = nf90_put_var(ncid0,VarID, precpbins(1:nprc), start = (/nrec0/))

          IF (level == 5) THEN
             iret = nf90_inq_varid(ncid0,'S_Rwsnw', VarID)
             iret = nf90_put_var(ncid0,VarID, snowbins(1:nsnw), start = (/nrec0/))
          END IF
       END IF
    end if

    ! Always saved: u, v, w, theta, P
    iret = nf90_inq_varid(ncid0, 'u', VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0, VarID, a_up(:,i1:i2,j1:j2), start=ibeg, count=icnt)
    iret = nf90_inq_varid(ncid0, 'v', VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0, VarID, a_vp(:,i1:i2,j1:j2), start=ibeg, count=icnt)
    iret = nf90_inq_varid(ncid0, 'w', VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0, VarID, a_wp(:,i1:i2,j1:j2), start=ibeg, count=icnt)
    iret = nf90_inq_varid(ncid0, 'theta', VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0, VarID, a_theta(:,i1:i2,j1:j2), start=ibeg, count=icnt)
    iret = nf90_inq_varid(ncid0, 'p', VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0, VarID, a_press(:,i1:i2,j1:j2), start=ibeg, count=icnt)

    ! Common variables (optional)
    iret = nf90_inq_varid(ncid0,'l',VarID) ! Liquid water
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0, VarID, a_rc(:,i1:i2,j1:j2), start=ibeg, count=icnt)
    iret = nf90_inq_varid(ncid0, 'rflx', VarID) ! Total radiative flux
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0, VarID, a_rflx(:,i1:i2,j1:j2), start=ibeg, count=icnt)
    iret = nf90_inq_varid(ncid0, 'stke', VarID) ! Subgrid TKE
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0, VarID, a_qp(:,i1:i2,j1:j2), start=ibeg, count=icnt)

    ! User-selected process rate outputs
    IF (nv4_proc>0) THEN
        DO bb=1,nv4_proc
            iret = nf90_inq_varid(ncid0, out_an_list(bb), VarID)
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid0,VarID,out_an_data(:,i1:i2,j1:j2,bb),start=ibeg,count=icnt)
        ENDDO
        out_an_data(:,:,:,:) = 0.
    ENDIF

    IF (level < 4) THEN ! Normal operation for levels 1-3
       ! Total water
       iret = nf90_inq_varid(ncid0,'q',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_rp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       ! Rain water mixing ratio and number concentration
       iret = nf90_inq_varid(ncid0,'r',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_rpp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       iret = nf90_inq_varid(ncid0,'n',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_npp(:,i1:i2,j1:j2),start=ibeg,count=icnt)

       ! Process rate statistics
       IF (level==3 .AND. stat_micro) THEN
          ! Coagulation (2)
          iret = nf90_inq_varid(ncid0,'coag_rr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_rr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'coag_nr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_nr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ! Condensation (2)
          iret = nf90_inq_varid(ncid0,'cond_rr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,cond_rr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'cond_nr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,cond_nr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ! Sedimentation (3)
          iret = nf90_inq_varid(ncid0,'sedi_rc',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_rc(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_rr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_rr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_nr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_nr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ! Autoconversion (2)
          iret = nf90_inq_varid(ncid0,'auto_rr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,auto_rr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'auto_nr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,auto_nr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ! Diagnostics (2)
          iret = nf90_inq_varid(ncid0,'diag_rr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_rr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'diag_nr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_nr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

    ELSE IF (level >= 4) THEN ! Operation with SALSA

       ! Total water mixing ratio
       iret = nf90_inq_varid(ncid0,'q',VarID)
       IF (iret==NF90_NOERR) THEN
          zvar(:,:,:) = a_rp(:,:,:) + &   ! Water vapor
                            a_rc(:,:,:) + &   ! Aerosol + cloud water
                            a_srp(:,:,:) + &   ! Rain water
                            a_ri(:,:,:) + & ! Ice water (level 5)
                            a_srs(:,:,:)      ! Snow water (level 5)
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Rain water mixing ratio
       iret = nf90_inq_varid(ncid0,'r',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_srp(:,i1:i2,j1:j2),start=ibeg,count=icnt)

       ! Ice water mixing ratio
       iret = nf90_inq_varid(ncid0,'i',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_ri(:,i1:i2,j1:j2),start=ibeg,count=icnt)

       ! Snow water mixing ratio
       iret = nf90_inq_varid(ncid0,'s',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_srs(:,i1:i2,j1:j2),start=ibeg,count=icnt)

       ! Aerosol, cloud, rain, ice and snow

       ! Total number of aerosol (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Naa',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL bulkNumc('aerosol','a',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF
       iret = nf90_inq_varid(ncid0,'S_Nab',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL bulkNumc('aerosol','b',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Mean wet radius of aerosol (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Rwaa',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL meanRadius('aerosol','a',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF
       iret = nf90_inq_varid(ncid0,'S_Rwab',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL meanRadius('aerosol','b',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Aerosol bin number concentration (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Naba',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarId,a_naerop(:,i1:i2,j1:j2,in1a:fn2a),start=ibegsd,count=icntaea)
       iret = nf90_inq_varid(ncid0,'S_Nabb',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_naerop(:,i1:i2,j1:j2,in2b:fn2b),start=ibegsd,count=icntaeb)

       ! Aerosol bin wet radius (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Rwaba',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nbins,nspec+1,a_naerop,a_maerop,nlim,a_Rawet,1)
          iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in1a:fn2a),start=ibegsd,count=icntaea)
       END IF
       iret = nf90_inq_varid(ncid0,'S_Rwabb',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nbins,nspec+1,a_naerop,a_maerop,nlim,a_Rawet,1)
          iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in2b:fn2b),start=ibegsd,count=icntaeb)
       END IF

       ! Total number of cloud droplets (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Nca',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL bulkNumc('cloud','a',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF
       iret = nf90_inq_varid(ncid0,'S_Ncb',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL bulkNumc('cloud','b',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Mean wet radius of cloud droplets (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Rwca',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL meanRadius('cloud','a',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF
       iret = nf90_inq_varid(ncid0,'S_Rwcb',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL meanRadius('cloud','b',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Cloud droplet bin number concentration (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Ncba',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_ncloudp(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
       iret = nf90_inq_varid(ncid0,'S_Ncbb',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_ncloudp(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)

       ! Cloud droplet bin wet radius (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Rwcba',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(ncld,nspec+1,a_ncloudp,a_mcloudp,nlim,a_Rcwet,2)
          iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
       END IF
       iret = nf90_inq_varid(ncid0,'S_Rwcbb',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(ncld,nspec+1,a_ncloudp,a_mcloudp,nlim,a_Rcwet,2)
          iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)
       END IF

       ! Total number of rain drops
       iret = nf90_inq_varid(ncid0,'S_Np',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL bulkNumc('precp','a',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Mean wet radius of rain drops
       iret = nf90_inq_varid(ncid0,'S_Rwp',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL meanRadius('precp','a',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Rain drop bin number concentration
       iret = nf90_inq_varid(ncid0,'S_Npb',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_nprecpp(:,i1:i2,j1:j2,1:nprc),start=ibegsd,count=icntpr)

       ! Rain drop bin wet radius
       iret = nf90_inq_varid(ncid0,'S_Rwpb',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nprc,nspec+1,a_nprecpp,a_mprecpp,prlim,a_Rpwet,3)
          iret = nf90_put_var(ncid0,VarID,a_Rpwet(:,i1:i2,j1:j2,1:nprc),start=ibegsd,count=icntpr)
       END IF

       ! Total number of ice (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Nia',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL bulkNumc('ice','a',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF
       iret = nf90_inq_varid(ncid0,'S_Nib',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL bulkNumc('ice','b',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Mean wet radius of ice (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Rwia',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL meanRadius('ice','a',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF
       iret = nf90_inq_varid(ncid0,'S_Rwib',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL meanRadius('ice','b',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Ice bin number concentration (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Niba',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_nicep(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
       iret = nf90_inq_varid(ncid0,'S_Nibb',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_nicep(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)

       ! Ice bin wet radius (regimes A and B)
       iret = nf90_inq_varid(ncid0,'S_Rwiba',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nice,nspec+1,a_nicep,a_micep,prlim,a_Riwet,4)
          iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
       END IF
       iret = nf90_inq_varid(ncid0,'S_Rwibb',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nice,nspec+1,a_nicep,a_micep,prlim,a_Riwet,4)
          iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)
       END IF

       ! Total number of snow
       iret = nf90_inq_varid(ncid0,'S_Ns',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL bulkNumc('snow','a',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Mean wet radius of snow
       iret = nf90_inq_varid(ncid0,'S_Rws',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL meanRadius('snow','a',zvar(:,:,:))
          iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

       ! Snow bin number concentration
       iret = nf90_inq_varid(ncid0,'S_Nsb',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_nsnowp(:,i1:i2,j1:j2,1:nsnw),start=ibegsd,count=icntsn)

       ! Snow bin wet radius
       iret = nf90_inq_varid(ncid0,'S_Rws',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nsnw,nspec+1,a_nsnowp,a_msnowp,prlim,a_Rswet,5)
          iret = nf90_put_var(ncid0,VarID,a_Rswet(:,i1:i2,j1:j2,1:nsnw),start=ibegsd,count=icntsn)
       END IF


       DO ee=1,nspec+1 ! Aerosol species and water
          ! Total mixing ratios
          iret = nf90_inq_varid(ncid0,'S_c'//TRIM(zspec(ee))//'a',VarID)
          IF (iret==NF90_NOERR) THEN
             CALL bulkMixrat(ee,'aerosol','ab',zvar(:,:,:))
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ENDIF

          iret = nf90_inq_varid(ncid0,'S_c'//TRIM(zspec(ee))//'c',VarID)
          IF (iret==NF90_NOERR) THEN
             CALL bulkMixrat(ee,'cloud','ab',zvar(:,:,:))
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ENDIF

          iret = nf90_inq_varid(ncid0,'S_c'//TRIM(zspec(ee))//'p',VarID)
          IF (iret==NF90_NOERR) THEN
             CALL bulkMixrat(ee,'precip','ab',zvar(:,:,:))
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ENDIF

          iret = nf90_inq_varid(ncid0,'S_c'//TRIM(zspec(ee))//'i',VarID)
          IF (iret==NF90_NOERR) THEN
             CALL bulkMixrat(ee,'ice','ab',zvar(:,:,:))
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ENDIF

          iret = nf90_inq_varid(ncid0,'S_c'//TRIM(zspec(ee))//'s',VarID)
          IF (iret==NF90_NOERR) THEN
             CALL bulkMixrat(ee,'snow','ab',zvar(:,:,:))
             iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ENDIF

          ! Mixing ratios for each bin
          iret = nf90_inq_varid(ncid0,'S_'//TRIM(zspec(ee))//'aa',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = in1a,fn2a
                CALL binSpecMixrat('aerosol',ee,bb,a_Rawet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in1a:fn2a),start=ibegsd,count=icntaea)
          ENDIF
          iret = nf90_inq_varid(ncid0,'S_'//TRIM(zspec(ee))//'ab',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = in2b,fn2b
                CALL binSpecMixrat('aerosol',ee,bb,a_Rawet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in2b:fn2b),start=ibegsd,count=icntaeb)
          ENDIF

          iret = nf90_inq_varid(ncid0,'S_'//TRIM(zspec(ee))//'ca',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = inp2a,fnp2a
                CALL binSpecMixrat('cloud',ee,bb,a_Rcwet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
          ENDIF
          iret = nf90_inq_varid(ncid0,'S_'//TRIM(zspec(ee))//'cb',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = inp2b,fnp2b
                CALL binSpecMixrat('cloud',ee,bb,a_Rcwet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)
          ENDIF

          iret = nf90_inq_varid(ncid0,'S_'//TRIM(zspec(ee))//'pb',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = 1,nprc
                CALL binSpecMixrat('precip',ee,bb,a_Rpwet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rpwet(:,i1:i2,j1:j2,1:nprc),start=ibegsd,count=icntpr)
          ENDIF

          iret = nf90_inq_varid(ncid0,'S_'//TRIM(zspec(ee))//'ia',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = inp2a,fnp2a
                CALL binSpecMixrat('ice',ee,bb,a_Riwet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
          ENDIF
          iret = nf90_inq_varid(ncid0,'S_'//TRIM(zspec(ee))//'ib',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = inp2b,fnp2b
                CALL binSpecMixrat('ice',ee,bb,a_Riwet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)
          ENDIF

          iret = nf90_inq_varid(ncid0,'S_'//TRIM(zspec(ee))//'sb',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = 1,nsnw
                CALL binSpecMixrat('snow',ee,bb,a_Rswet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rswet(:,i1:i2,j1:j2,1:nsnw),start=ibegsd,count=icntsn)
          ENDIF
       ENDDO

       ! Gas mixing ratios
       DO ee=1,ngases
          iret = nf90_inq_varid(ncid0,'S_c'//TRIM(zgas(ee))//'g',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_gaerop(:,i1:i2,j1:j2,ee),start=ibeg,count=icnt)
       ENDDO

       ! Process rate statistics
       IF (stat_micro) THEN
          ! Coagulation (10)
          iret = nf90_inq_varid(ncid0,'coag_ra',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_ra(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'coag_na',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_na(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'coag_rc',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_rc(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'coag_nc',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_nc(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'coag_rr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_rr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'coag_nr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_nr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'coag_ri',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_ri(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'coag_ni',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_ni(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'coag_rs',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_rs(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'coag_ns',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,coag_ns(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ! Condensation (5)
          iret = nf90_inq_varid(ncid0,'cond_ra',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,cond_ra(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'cond_rc',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,cond_rc(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'cond_rr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,cond_rr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'cond_ri',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,cond_ri(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'cond_rs',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,cond_rs(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ! Sedimentation (10)
          iret = nf90_inq_varid(ncid0,'sedi_ra',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_ra(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_na',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_na(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_rc',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_rc(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_nc',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_nc(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_rr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_rr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_nr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_nr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_ri',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_ri(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_ni',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_ni(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_rs',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_rs(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'sedi_ns',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,sedi_ns(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ! Activation, ice nucleation and autoconversion (8)
          iret = nf90_inq_varid(ncid0,'act_rc',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,cact_rc(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'act_nc',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,cact_nc(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'auto_rr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,auto_rr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'auto_nr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,auto_nr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'nucl_ri',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,nucl_ri(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'nucl_ni',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,nucl_ni(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'auto_rs',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,auto_rs(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'auto_ns',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,auto_ns(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ! Diagnostics (10)
          iret = nf90_inq_varid(ncid0,'diag_ra',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_ra(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'diag_na',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_na(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'diag_rc',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_rc(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'diag_nc',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_nc(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'diag_rr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_rr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'diag_nr',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_nr(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'diag_ri',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_ri(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'diag_ni',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_ni(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'diag_rs',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_rs(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'diag_ns',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,diag_ns(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          ! Melting (4)
          iret = nf90_inq_varid(ncid0,'melt_ri',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,melt_ri(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'melt_ni',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,melt_ni(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'melt_rs',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,melt_rs(:,i1:i2,j1:j2),start=ibeg,count=icnt)
          iret = nf90_inq_varid(ncid0,'melt_ns',VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,melt_ns(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       END IF

    END IF

    if (myid==0) print "(/' ',12('-'),'   Record ',I3,' to: ',A60)",    &
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
    if (myid == 0) print "(/' ',12('-'),'   History write to: ',A30)" &
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

    IF (nudge_theta/=0) write(10) theta_ref
    IF (nudge_rv/=0) write(10) rv_ref
    IF (nudge_u/=0) write(10) u_ref
    IF (nudge_v/=0) write(10) v_ref
    IF (level>3 .AND. nudge_ccn/=0) write(10) aero_ref
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

       IF (nudge_theta/=0) THEN
          ALLOCATE(theta_ref(nzp))
          READ(10) theta_ref
       end if
       IF (nudge_rv/=0) THEN
          ALLOCATE(rv_ref(nzp))
          READ(10) rv_ref
       end if
       IF (nudge_u/=0) THEN
          ALLOCATE(u_ref(nzp))
          READ(10) u_ref
       end if
       IF (nudge_v/=0) THEN
          ALLOCATE(v_ref(nzp))
          READ(10) v_ref
       end if
       IF (level>3 .AND. nudge_ccn/=0) THEN
          ALLOCATE(aero_ref(nzp,nbins))
          READ(10) aero_ref
       end if
       nudge_init=.FALSE.

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
          u0 = u0 + umx - umean
       end if
       if (vmx /= vmean) then
          if (myid == 0) print "('  vmean changed  -  ',2f8.2)",vmean,vmx
          a_vp = a_vp + vmx - vmean
          v0 = v0 +vmx - vmean
       end if

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
  SUBROUTINE bulkMixrat(mm,ipart,itype,mixrat)
    USE mo_submctl, ONLY : fnp2a,inp2b,in2b,fn2a

    INTEGER, INTENT(in) :: mm ! Index to an active species (1, 2,... nspec+1)
    CHARACTER(len=*), INTENT(in) :: ipart  ! This should be aerosol, cloud, rain, ice or snow
    CHARACTER(len=*), INTENT(in) :: itype  ! Select bin regime: a or b
    REAL, INTENT(out) :: mixrat(nzp,nxp,nyp)

    INTEGER :: istr,iend

    mixrat = 0.

    ! Given in kg/kg
    SELECT CASE(ipart)
       CASE('aerosol')
          IF (itype == 'ab') THEN
             istr = (mm-1)*nbins + 1
             iend = mm*nbins
          ELSEIF (itype == 'a') THEN
             istr = (mm-1)*nbins + 1
             iend = (mm-1)*nbins + fn2a
          ELSE IF (itype == 'b') THEN
             istr = (mm-1)*nbins + in2b
             iend = mm*nbins
          ELSE
             STOP 'bulkMixrat: Invalid aerosol bin regime selection'
          END IF
          mixrat(:,:,:) = SUM(a_maerop(:,:,:,istr:iend),DIM=4)
       CASE('cloud')
          IF (itype == 'ab') THEN
             istr = (mm-1)*ncld + 1
             iend = mm*ncld
          ELSEIF (itype == 'a') THEN
             istr = (mm-1)*ncld + 1
             iend = (mm-1)*ncld + fnp2a
          ELSE IF (itype == 'b') THEN
             istr = (mm-1)*ncld + inp2b
             iend = mm*ncld
          ELSE
             STOP 'bulkMixrat: Invalid cloud bin regime selection'
          END IF
          mixrat(:,:,:) = SUM(a_mcloudp(:,:,:,istr:iend),DIM=4)
       CASE('precp')
          istr = (mm-1)*nprc + 1
          iend = (mm-1)*nprc + nprc
          mixrat(:,:,:) = SUM(a_mprecpp(:,:,:,istr:iend),DIM=4)
       CASE('ice')
          IF (itype == 'ab') THEN
             istr = (mm-1)*nice + 1
             iend = mm*nice
          ELSEIF (itype == 'a') THEN
             istr = (mm-1)*nice + 1
             iend = (mm-1)*nice + fnp2a
          ELSE IF (itype == 'b') THEN
             istr = (mm-1)*nice + inp2b
             iend = mm*nice
          ELSE
             STOP 'bulkMixrat: Invalid ice bin regime selection'
          END IF
          mixrat(:,:,:) = SUM(a_micep(:,:,:,istr:iend),DIM=4)
       CASE('snow')
          istr = (mm-1)*nsnw + 1
          iend = (mm-1)*nsnw + nsnw
          mixrat(:,:,:) = SUM(a_msnowp(:,:,:,istr:iend),DIM=4)
       CASE DEFAULT
          STOP 'bulkMixrat: Invalid particle type'
    END SELECT

  END SUBROUTINE bulkMixrat
  !
  ! ----------------------------------------------
  ! Subroutine binSpecMixrat: Calculate the mixing
  ! ratio of selected aerosol species in individual
  ! bins.
  !
  ! Juha Tonttila, FMI, 2015
  SUBROUTINE binSpecMixrat(ipart,mm,ibin,mixr)

    INTEGER, INTENT(in) :: mm ! Index to an active species (1, 2,... nspec+1)
    CHARACTER(len=*), INTENT(in) :: ipart  ! This should be aerosol, cloud, rain, ice or snow
    INTEGER, INTENT(in) :: ibin
    REAL, INTENT(out) :: mixr(nzp,nxp,nyp)

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
       CASE DEFAULT
          STOP 'binSpecMixrat: Invalid particle type'
    END SELECT

  END SUBROUTINE binSpecMixrat
  !
  ! ----------------------------------------------
  ! Subroutine bulkNumc: Calculate the total number
  ! concentration of particles of given type
  !
  ! Juha Tonttila, FMI, 2015
  !
  SUBROUTINE bulkNumc(ipart,itype,numc)
    USE mo_submctl, ONLY : fnp2a,inp2b,in2b,fn2a

    CHARACTER(len=*), INTENT(in) :: ipart
    CHARACTER(LEN=*), INTENT(in) :: itype
    REAL, INTENT(out) :: numc(nzp,nxp,nyp)
    INTEGER :: istr,iend

    ! Outputs #/kg
    ! No concentration limits (nlim or prlim) for number

    SELECT CASE(ipart)
       CASE('aerosol')
          IF (itype == 'ab') THEN ! Note: 1a, 2a and 2b combined
             istr = 1
             iend = nbins
          ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
             istr = 1
             iend = fn2a
          ELSE IF (itype == 'b') THEN ! 2b
             istr = in2b
             iend = nbins
          ELSE
             STOP 'bulkNumc: Invalid bin selection'
          END IF
          numc(:,:,:) = SUM(a_naerop(:,:,:,istr:iend),DIM=4)
       CASE('cloud')
          IF (itype == 'ab') THEN
             istr = 1
             iend = ncld
          ELSE IF (itype == 'a') THEN
             istr = 1
             iend = fnp2a
          ELSE IF (itype == 'b') THEN
             istr = inp2b
             iend = ncld
          ELSE
             STOP 'bulkNumc: Invalid bin selection'
          END IF
          numc(:,:,:) = SUM(a_ncloudp(:,:,:,istr:iend),DIM=4)
       CASE('precp')
          istr = 1
          iend = nprc
          numc(:,:,:) = SUM(a_nprecpp(:,:,:,istr:iend),DIM=4)
        CASE('ice')
          IF (itype == 'ab') THEN
             istr = 1
             iend = nice
          ELSE IF (itype == 'a') THEN
             istr = 1
             iend = fnp2a
          ELSE IF (itype == 'b') THEN
             istr = inp2b
             iend = nice
          ELSE
             STOP 'bulkNumc: Invalid bin selection'
          END IF
          numc(:,:,:) = SUM(a_nicep(:,:,:,istr:iend),DIM=4)
       CASE('snow')
          istr = 1
          iend = nsnw
          numc(:,:,:) = SUM(a_nsnowp(:,:,:,istr:iend),DIM=4)
       CASE DEFAULT
          STOP 'binSpecMixrat: Invalid particle type'
    END SELECT

  END SUBROUTINE bulkNumc


  !
  ! -------------------------------------------------
  ! SUBROUTINE meanRadius
  ! Gets the mean wet (total number of species = nspec+1) radius for particles.
  !
  SUBROUTINE meanRadius(ipart,itype,rad)
    USE mo_submctl, ONLY : fnp2a,inp2b,fn2a,in2b,nlim,prlim
    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: ipart
    CHARACTER(len=*), INTENT(in) :: itype
    REAL, INTENT(out) :: rad(nzp,nxp,nyp)

    INTEGER :: istr,iend

    rad = 0.

    SELECT CASE(ipart)
    CASE('aerosol')
       IF (itype == 'ab') THEN ! Note: 1a, 2a and 2b combined
          istr = 1
          iend = nbins
       ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
          istr = 1
          iend = fn2a
       ELSE IF (itype == 'b') THEN
          istr = in2b
          iend = nbins
       ELSE
          STOP 'meanRadius: Invalid bin regime selection (aerosol)'
       END IF
       CALL getRadius(istr,iend,nbins,nspec+1,a_naerop,a_maerop,nlim,rad,1)
    CASE('cloud')
       IF (itype == 'ab') THEN
          istr = 1
          iend = ncld
       ELSE IF (itype == 'a') THEN
          istr = 1
          iend = fnp2a
       ELSE IF (itype == 'b') THEN
          istr = inp2b
          iend = ncld
       ELSE
          STOP 'meanRadius: Invalid bin regime selection (cloud)'
       END IF
       CALL getRadius(istr,iend,ncld,nspec+1,a_ncloudp,a_mcloudp,nlim,rad,2)
    CASE('precp')
       istr = 1
       iend = nprc
       CALL getRadius(istr,iend,nprc,nspec+1,a_nprecpp,a_mprecpp,prlim,rad,3)
    CASE('ice')
       IF (itype == 'ab') THEN
          istr = 1
          iend = nice
       ELSE IF (itype == 'a') THEN
          istr = 1
          iend = fnp2a
       ELSE IF (itype == 'b') THEN
          istr = inp2b
          iend = nice
       ELSE
          STOP 'meanRadius: Invalid bin regime selection (ice)'
       END IF
       CALL getRadius(istr,iend,nice,nspec+1,a_nicep,a_micep,prlim,rad,4)
    CASE('snow')
       istr = 1
       iend = nsnw
       CALL getRadius(istr,iend,nsnw,nspec+1,a_nsnowp,a_msnowp,prlim,rad,5)
    CASE DEFAULT
          STOP 'meanRadius: Invalid particle type'
    END SELECT

  contains

   SUBROUTINE getRadius(zstr,zend,nn,n4,numc,mass,numlim,zrad,flag)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: nn, n4 ! Number of bins (nn) and aerosol species (n4)
    INTEGER, INTENT(in) :: zstr,zend  ! Start and end index for averaging
    REAL, INTENT(in) :: numc(nzp,nxp,nyp,nn)
    REAL, INTENT(in) :: mass(nzp,nxp,nyp,nn*n4)
    REAL, INTENT(in) :: numlim
    INTEGER, INTENT(IN) :: flag
    REAL, INTENT(out) :: zrad(nzp,nxp,nyp)

    INTEGER :: k,i,j,bin
    REAL :: tot, rwet, tmp(n4)

    zrad(:,:,:)=0.
    DO j = 3,nyp-2
      DO i = 3,nxp-2
        DO k = 1,nzp
          tot=0.
          rwet=0.
          DO bin = zstr,zend
            IF (numc(k,i,j,bin)>numlim) THEN
              tot=tot+numc(k,i,j,bin)
              tmp(:)=mass(k,i,j,bin:(n4-1)*nn+bin:nn)
              rwet=rwet+calc_eff_radius(n4,numc(k,i,j,bin),tmp,flag)*numc(k,i,j,bin)
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
  ! Calculates wet radius for each bin in the whole domain
  SUBROUTINE getBinRadius(nn,n4,numc,mass,numlim,zrad,flag)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: nn, n4 ! Number of bins (nn) and aerosol species (n4)
    REAL, INTENT(in) :: numc(nzp,nxp,nyp,nn)
    REAL, INTENT(in) :: mass(nzp,nxp,nyp,nn*n4)
    REAL, INTENT(in) :: numlim
    INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)
    REAL, INTENT(out) :: zrad(nzp,nxp,nyp,nn)

    INTEGER :: k,i,j,bin
    REAL :: tmp(n4)

    zrad(:,:,:,:)=0.
    DO j = 3,nyp-2
      DO i = 3,nxp-2
        DO k = 1,nzp
          DO bin = 1,nn
            IF (numc(k,i,j,bin)>numlim) THEN
              tmp(:)=mass(k,i,j,bin:(n4-1)*nn+bin:nn)
              zrad(k,i,j,bin)=calc_eff_radius(n4,numc(k,i,j,bin),tmp,flag)
            ENDIF
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE getBinRadius


  !********************************************************************
  !
  ! Function for calculating effective (wet) radius for any particle type
  ! - Aerosol, cloud and rain are spherical
  ! - Snow and ice can be irregular and their densities can be size-dependent
  !
  ! Edit this function when needed (also update CalcDimension in mo_submctl.f90)
  !
  ! Correct dimension is needed for irregular particles (e.g. ice and snow) for calculating fall speed (deposition and coagulation)
  ! and capacitance (condensation). Otherwise compact spherical structure can be expected,
  !
  REAL FUNCTION calc_eff_radius(n,numc,mass,flag)
    USE mo_submctl, ONLY : pi6, dens, dens_ice, dens_snow
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n ! Number of species
    INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)
    REAL, INTENT(IN) :: numc, mass(n)

    calc_eff_radius=0.

    ! Don't calculate if very low number concentration
    IF (numc<1e-15) RETURN

    IF (flag==4) THEN   ! Ice
        ! Spherical ice
        calc_eff_radius=0.5*( SUM(mass(:)/dens_ice(1:n))/numc/pi6)**(1./3.)
    ELSEIF (flag==5) THEN   ! Snow
        ! Spherical snow
        calc_eff_radius=0.5*( SUM(mass(:)/dens_snow(1:n))/numc/pi6)**(1./3.)
    ELSE
        ! Radius from total volume of a spherical particle or aqueous droplet
        calc_eff_radius=0.5*( SUM(mass(:)/dens(1:n))/numc/pi6)**(1./3.)
    ENDIF

  END FUNCTION calc_eff_radius
  !********************************************************************

end module grid

