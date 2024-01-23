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
  real              :: umean = 0.          ! Galilean transformation
  real              :: vmean = 0.          ! Galilean transformation
  integer           :: igrdtyp = 1         ! vertical grid type
  integer           :: isgstyp = 1         ! sgs model type
  integer           :: iradtyp = 0         ! radiation model type
  integer           :: level   = 0         ! thermodynamic level
  integer           :: lev_sb = 5          ! thermodynamic level for SB microphysics (level=0)
  integer           :: naddsc  = 0         ! number of additional scalars;
  integer           :: nfpt = 10           ! number of rayleigh friction points
  real              :: distim = 300.0      ! dissipation timescale

  real              :: sst=283.            ! Surface temperature

  LOGICAL :: sed_aero = .TRUE.,  &
             sed_cloud = .TRUE., &
             sed_precp = .TRUE., &
             sed_ice = .TRUE., &
             sed_snow = .TRUE.

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

  ! Marine emissions
  logical :: ifSeaSpray = .false. ! Aerosol
  logical :: ifSeaVOC = .false.   ! Isoprene and monoterpenes
  real :: sea_tspinup = 0.        ! Spin-up time (s) for marine emissions

  character (len=80):: expnme = 'Default' ! Experiment name
  character (len=80):: filprf = 'x'       ! File Prefix
  character (len=7) :: runtype = 'INITIAL'! Run Type Selection

  REAL              :: Tspinup = 7200.    ! Spinup period in seconds (added by Juha)

  ! User control of analysis outputs
  INTEGER, PARAMETER :: maxn_list=100
  CHARACTER(len=7), dimension(maxn_list), SAVE :: anl_include='       ', anl_exclude='       '
  CHARACTER(len=7), dimension(maxn_list), SAVE :: out_an_list='       ', user_an_list='       '
  INTEGER, SAVE :: nv4_proc=0, nv4_user=0
  REAL, SAVE, ALLOCATABLE :: out_an_data(:,:,:,:)

  integer           :: nz, nxyzp, nxyp
  real              :: dxi, dyi, dtl, psrf
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
  real, pointer :: a_tp(:,:,:),a_tt(:,:,:) ! (Ice-)liquid water potential temperature
  real, pointer :: a_rp(:,:,:),a_rt(:,:,:) ! Water vapour for SALSA; total water for SB
  ! Seifert & Beheng tracers: mass (kg/kg) and number (#/kg)
  real, pointer :: a_rpp(:,:,:),a_rpt(:,:,:),a_npp(:,:,:),a_npt(:,:,:) ! Rain
  real, pointer :: a_rip(:,:,:),a_rit(:,:,:),a_nip(:,:,:),a_nit(:,:,:) ! SB level 4 & 5 ice
  real, pointer :: a_rsp(:,:,:),a_rst(:,:,:),a_nsp(:,:,:),a_nst(:,:,:) ! Snow
  real, pointer :: a_rgp(:,:,:),a_rgt(:,:,:),a_ngp(:,:,:),a_ngt(:,:,:) ! Graupel
  real, pointer :: a_rhp(:,:,:),a_rht(:,:,:),a_nhp(:,:,:),a_nht(:,:,:) ! Hail
  ! Subgrid TKE and a scratch variable
  real, pointer :: a_qp(:,:,:),a_qt(:,:,:)
  real, pointer :: a_sp(:,:,:),a_st(:,:,:)
  ! SALSA tracers
  ! -- Number concentrations (#/kg)
  REAL, POINTER :: a_naerop(:,:,:,:), a_naerot(:,:,:,:),   &
                   a_ncloudp(:,:,:,:), a_ncloudt(:,:,:,:), &
                   a_nprecpp(:,:,:,:), a_nprecpt(:,:,:,:), &
                   a_nicep(:,:,:,:),   a_nicet(:,:,:,:), &
                   a_nsnowp(:,:,:,:),  a_nsnowt(:,:,:,:)
  ! -- Mass concentrations (kg/kg)
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
  REAL, ALLOCATABLE :: a_rsl(:,:,:)    ! Water vapor saturation mixing ratio
  REAL, ALLOCATABLE :: a_rsi(:,:,:)    ! Water vapor saturation mixing ratio over ice
  REAL, ALLOCATABLE :: a_dn(:,:,:)     ! Air density
  REAL, ALLOCATABLE :: a_edr(:,:,:)    ! Eddy dissipation rate
  REAL, ALLOCATABLE :: a_temp(:,:,:)   ! Air temperature
  !
  ! radiation
  real, allocatable, dimension (:,:,:) :: a_rflx, a_sflx, &
       a_fus, a_fds, a_fuir, a_fdir
  real, allocatable :: albedo(:,:)
  !
  ! surface
  real, allocatable :: a_ustar(:,:)
  real, allocatable :: a_tstar(:,:)
  real, allocatable :: a_rstar(:,:)
  real, allocatable :: uw_sfc(:,:)
  real, allocatable :: vw_sfc(:,:)
  real, allocatable :: ww_sfc(:,:)
  real, allocatable :: wt_sfc(:,:)
  real, allocatable :: wq_sfc(:,:)
  real, allocatable :: obl(:,:)
  !
  ! microphysics/precipitation
  real, allocatable, dimension(:,:,:) :: aerin, cldin, precip, icein, snowin, grin, hailin
  !
  integer :: nscl = 1
  integer, private, save :: ncid0
  integer, private, save :: nrec0
  character (len=80), private, save :: fname
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

    allocate (a_temp(nzp,nxp,nyp),a_dn(nzp,nxp,nyp),a_edr(nzp,nxp,nyp), &
        a_rc(nzp,nxp,nyp),a_rsl(nzp,nxp,nyp),a_ri(nzp,nxp,nyp),a_rsi(nzp,nxp,nyp))
    a_temp(:,:,:) = 0.
    a_dn(:,:,:) = 0.
    a_edr(:,:,:) = 0.
    a_rc(:,:,:) = 0.
    a_rsl(:,:,:) = 0.
    a_ri(:,:,:) = 0.
    a_rsi(:,:,:) = 0.
    memsize = memsize + nxyzp*7

    ! Juha: Stuff that's allocated if SALSA is NOT used
    !-----------------------------------------------------
    IF (level < 4) THEN
       if (level < 0) then
            WRITE(*,*) 'Level < 0 not accepted!'
            STOP
       end if

       allocate (a_rv(nzp,nxp,nyp))
       a_rv(:,:,:) = 0.
       memsize = memsize + nxyzp

       ! Prognostic scalars: temperature + total water + tke (isgstyp> 1) + additional scalars + ...
       !    rain mass and number (level=3) or
       !    rain mass and number (level=0 & lev_sb=3)
       !    rain and ice mass and number, and snow and graupel mass (level=0 & lev_sb=4)
       !    rain, ice, snow, graupel and hail mass and number (level=0 & lev_sb=5)
       nscl = 2+naddsc
       if (level == 3 .OR. level == 0) nscl = nscl+2 ! rain
       if (level == 0 .AND. lev_sb ==4) nscl = nscl+4 ! + ice and snow and graupel mass
       if (level == 0 .AND. lev_sb ==5) nscl = nscl+8 ! + ice, snow, graupel and hail
       if (isgstyp > 1) nscl = nscl+1 ! tke

       allocate (a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
       a_sclrp(:,:,:,:) = 0.
       a_sclrt(:,:,:,:) = 0.
       memsize = memsize + 2*nscl*nxyzp

       a_tp=>a_sclrp(:,:,:,1)
       a_tt=>a_sclrt(:,:,:,1)
       a_rp=>a_sclrp(:,:,:,2)
       a_rt=>a_sclrt(:,:,:,2)
       if (level == 3 .OR. level == 0) then
          a_rpp=>a_sclrp(:,:,:,3)
          a_rpt=>a_sclrt(:,:,:,3)
          a_npp=>a_sclrp(:,:,:,4)
          a_npt=>a_sclrt(:,:,:,4)
       else
          ALLOCATE (tmp_prcp(nzp,nxp,nyp,2),tmp_prct(nzp,nxp,nyp,2))
          tmp_prcp=0.; tmp_prct=0.
          a_rpp=>tmp_prcp(:,:,:,1)
          a_rpt=>tmp_prct(:,:,:,1)
          a_npp=>tmp_prcp(:,:,:,2)
          a_npt=>tmp_prct(:,:,:,2)
       end if
       if (level == 0) then
          if (lev_sb == 5) then
             a_rip => a_sclrp(:,:,:,5); a_rit => a_sclrt(:,:,:,5)
             a_nip => a_sclrp(:,:,:,6); a_nit => a_sclrt(:,:,:,6)
             a_rsp => a_sclrp(:,:,:,7); a_rst => a_sclrt(:,:,:,7)
             a_nsp => a_sclrp(:,:,:,8); a_nst => a_sclrt(:,:,:,8)
             a_rgp => a_sclrp(:,:,:,9); a_rgt => a_sclrt(:,:,:,9)
             a_ngp => a_sclrp(:,:,:,10); a_ngt => a_sclrt(:,:,:,10)
             a_rhp => a_sclrp(:,:,:,11); a_rht => a_sclrt(:,:,:,11)
             a_nhp => a_sclrp(:,:,:,12); a_nht => a_sclrt(:,:,:,12)
          elseif (lev_sb == 4) then
             a_nip => a_sclrp(:,:,:,5); a_nit => a_sclrt(:,:,:,5)
             a_rip => a_sclrp(:,:,:,6); a_rit => a_sclrt(:,:,:,6)
             a_rsp => a_sclrp(:,:,:,7); a_rst => a_sclrt(:,:,:,7)
             a_rgp => a_sclrp(:,:,:,8); a_rgt => a_sclrt(:,:,:,8)
             ALLOCATE (tmp_icep(nzp,nxp,nyp,4),tmp_icet(nzp,nxp,nyp,4))
             tmp_icep(:,:,:,:) = 0.; tmp_icet(:,:,:,:) = 0.
             a_nsp => tmp_icep(:,:,:,1); a_nst => tmp_icet(:,:,:,1)
             a_ngp => tmp_icep(:,:,:,2); a_ngt => tmp_icet(:,:,:,2)
             a_rhp => tmp_icep(:,:,:,3); a_rht => tmp_icet(:,:,:,3)
             a_nhp => tmp_icep(:,:,:,4); a_nht => tmp_icet(:,:,:,4)
          elseif (lev_sb == 3) then
             ALLOCATE (tmp_icep(nzp,nxp,nyp,8),tmp_icet(nzp,nxp,nyp,8))
             tmp_icep(:,:,:,:) = 0.; tmp_icet(:,:,:,:) = 0.
             a_rip => tmp_icep(:,:,:,1); a_rit => tmp_icet(:,:,:,1)
             a_nip => tmp_icep(:,:,:,2); a_nit => tmp_icet(:,:,:,2)
             a_rsp => tmp_icep(:,:,:,3); a_rst => tmp_icet(:,:,:,3)
             a_nsp => tmp_icep(:,:,:,4); a_nst => tmp_icet(:,:,:,4)
             a_rgp => tmp_icep(:,:,:,5); a_rgt => tmp_icet(:,:,:,5)
             a_ngp => tmp_icep(:,:,:,6); a_ngt => tmp_icet(:,:,:,6)
             a_rhp => tmp_icep(:,:,:,7); a_rht => tmp_icet(:,:,:,7)
             a_nhp => tmp_icep(:,:,:,8); a_nht => tmp_icet(:,:,:,8)
          else
             WRITE(*,*) 'Invalid SB level (lev_sb) for level=0:',lev_sb
             STOP
          end if
       end if
       if (isgstyp > 1) then
          a_qp=>a_sclrp(:,:,:,nscl - naddsc)
          a_qt=>a_sclrt(:,:,:,nscl - naddsc)
       end if

    !Juha: Stuff that's allocated when SALSA is used
    !---------------------------------------------------
    ELSE IF (level >= 4) THEN

       allocate ( a_srp(nzp,nxp,nyp), a_snrp(nzp,nxp,nyp), &
                 a_srs(nzp,nxp,nyp), a_snrs(nzp,nxp,nyp) )
       a_srp(:,:,:) = 0.
       a_snrp(:,:,:) = 0.
       a_srs(:,:,:) = 0.
       a_snrs(:,:,:) = 0.
       memsize = memsize + 4*nxyzp

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
    END IF ! level

    !----------------------------------------------------

    allocate (a_ustar(nxp,nyp),a_tstar(nxp,nyp),a_rstar(nxp,nyp))
    allocate (uw_sfc(nxp,nyp),vw_sfc(nxp,nyp),ww_sfc(nxp,nyp))
    allocate (wt_sfc(nxp,nyp),wq_sfc(nxp,nyp))
    allocate (obl(nxp,nyp))

    allocate(cldin(nzp,nxp,nyp),precip(nzp,nxp,nyp))
    precip = 0.
    cldin = 0.
    memsize = memsize + nxyzp*2

    if (level >= 4) then
       allocate(aerin(nzp,nxp,nyp),icein(nzp,nxp,nyp),snowin(nzp,nxp,nyp))
       aerin = 0.
       icein = 0.
       snowin = 0.
       memsize = memsize + nxyzp*3
    elseif (level == 0) then
       allocate(icein(nzp,nxp,nyp),snowin(nzp,nxp,nyp),grin(nzp,nxp,nyp),hailin(nzp,nxp,nyp))
       icein = 0.
       snowin = 0.
       grin = 0.
       hailin = 0.
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
    USE mo_submctl, ONLY : fn1a,fn2a,fn2b
    IMPLICIT NONE
    real, intent (in) :: time
    ! Dimensions (time, x, y, x, and SALSA bins) and constants (u0, v0, dn0) are saved
    ! during initialization, and common variables (u, v, w, theta, p) are always saved.
    INTEGER, PARAMETER :: n_dims=14, n_base=19
    character(len=7) :: s_dims(n_dims) = (/ &
         'time   ','zt     ','zm     ','xt     ','xm     ','yt     ','ym     ', & ! 1-7
         'u0     ','v0     ','dn0    ','B_Rd12a','B_Rd2ab','B_Rwprc','B_Rwsnw'/)  ! 8-14
    character(len=7) :: s_base(n_base) = (/ &
         'u      ','v      ','w      ','theta  ','p      ','stke   ','rflx   ', & ! 1-7
         'q      ','l      ','r      ','n      ','i      ','s      ','g      ', & ! 8-14
         'ni     ','h      ','ns     ','ng     ','nh     '/) ! 15-19
    LOGICAL, SAVE :: b_dims(n_dims)=.TRUE., b_base(n_base)=.TRUE.
    CHARACTER (len=7), ALLOCATABLE :: sanal(:), stot(:)
    LOGICAL, ALLOCATABLE :: btot(:)

    ! Local variables
    INTEGER :: i, e, n, nvar0
    character (len=7) :: v_snm='sxx    '
    LOGICAL :: found, lbinanl=.FALSE.

    b_base(6) = (isgstyp > 1)
    b_base(7) = (iradtyp > 1)


    ! Allocate data for user selected process rate outputs (see init_stat in stat.f90)
    ALLOCATE ( out_an_data(nzp,nxp,nyp,nv4_proc) )
    out_an_data(:,:,:,:) = 0.

    IF (level < 4) THEN  ! Standard operation for levels 1-3
        b_dims(11:14) = .FALSE. ! SALSA bins
        b_base(10:11) = (level==3 .OR. level==0) ! Rain
        b_base(12:15) = (level==0 .AND. lev_sb>=4) ! Ice, snow mass and graupel mass
        b_base(16:19) = (level==0 .AND. lev_sb==5) ! Ice, snow, graupel and hail

       ! Merge logical and name arrays
       i=n_dims+n_base+nv4_proc+nv4_user+naddsc
       ALLOCATE( btot(i), stot(i) )
       i=1; e=n_dims
       btot(i:e)=b_dims; stot(i:e)=s_dims
       i=e+1; e=e+n_base
       btot(i:e)=b_base; stot(i:e)=s_base
       IF (nv4_proc>0) THEN
          i=e+1; e=e+nv4_proc
          btot(i:e)=.TRUE.; stot(i:e)=out_an_list(1:nv4_proc)
       END IF
       IF (nv4_user>0) THEN
          i=e+1; e=e+nv4_user
          btot(i:e)=.TRUE.; stot(i:e)=user_an_list(1:nv4_user)
       END IF
    ELSE IF (level >= 4) THEN ! Operation with SALSA
       b_base(10)=.NOT.no_prog_prc ! Rain water
       b_base(11)=.FALSE. ! Rain number
       b_base(12) = (level>4 .AND. .NOT. no_prog_ice) ! Ice
       b_base(13) = (level>4 .AND. .NOT. no_prog_snw) ! Snow
       b_base(14:19) = .FALSE. ! SB outputs

       ! Dimensions for bin dependent outputs
       lbinanl = ANY(INDEX(user_an_list,'B_')>0)
       b_dims(11) = lbinanl
       b_dims(12) = lbinanl
       b_dims(13) = lbinanl .AND. (.NOT. no_prog_prc)
       b_dims(14) = lbinanl .AND. (.NOT. no_prog_snw) .AND. (level>4)

       ! Merge logical and name arrays
       i=n_dims+n_base+nv4_proc+nv4_user+naddsc
       ALLOCATE( btot(i), stot(i) )
       i=1; e=n_dims
       btot(i:e)=b_dims; stot(i:e)=s_dims
       i=e+1; e=e+n_base
       btot(i:e)=b_base; stot(i:e)=s_base
       IF (nv4_proc>0) THEN
          i=e+1; e=e+nv4_proc
          btot(i:e)=.TRUE.; stot(i:e)=out_an_list(1:nv4_proc)
       END IF
       IF (nv4_user>0) THEN
          i=e+1; e=e+nv4_user
          btot(i:e)=.TRUE.; stot(i:e)=user_an_list(1:nv4_user)
       END IF

       ! Remove bin-dependent outputs from the user_an_list, because
       ! these are calculated during the normal call to write_anal
       IF (nv4_user>0) THEN
           n=0
           DO i=1,nv4_user
               IF ('B_'/=user_an_list(i)(1:2)) THEN
                   n=n+1
                   IF (n<i) user_an_list(n)=user_an_list(i)
               ENDIF
           ENDDO
           nv4_user=n
       ENDIF
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
                               zspec

    real, intent (in) :: time

    integer :: iret, VarID, bb, ee
    integer :: ibeg(4), icnt(4), i1, i2, j1, j2
    INTEGER :: ibegsd(5), icntaea(5), icntaeb(5), icntpr(5), icntsn(5)
    REAL :: zvar(nzp,nxp,nyp)
    REAL :: a_Rawet(nzp,nxp,nyp,nbins), a_Rcwet(nzp,nxp,nyp,ncld),a_Rpwet(nzp,nxp,nyp,nprc), &
          a_Riwet(nzp,nxp,nyp,nice),a_Rswet(nzp,nxp,nyp,nsnw)
    CHARACTER(LEN=3) :: nam

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

       IF (level >= 4) THEN
          iret = nf90_inq_varid(ncid0,'B_Rd12a', VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0, VarID, aerobins(in1a:fn2a), start = (/nrec0/))

          iret = nf90_inq_varid(ncid0,'B_Rd2ab', VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID, aerobins(in2a:fn2a), start = (/nrec0/))

          iret = nf90_inq_varid(ncid0,'B_Rwprc', VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID, precpbins(1:nprc), start = (/nrec0/))

          iret = nf90_inq_varid(ncid0,'B_Rwsnw', VarID)
          IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID, snowbins(1:nsnw), start = (/nrec0/))
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

    ! Additional scalars
    IF (naddsc>0) THEN
        DO bb=1,naddsc
            write(nam,"('s',i2.2)") bb
            iret = nf90_inq_varid(ncid0, nam, VarID)
            ee=nscl-naddsc+bb
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_sclrp(:,i1:i2,j1:j2,ee), &
                start=ibeg,count=icnt)
        ENDDO
    ENDIF

    ! User-selected process rate outputs
    IF (nv4_proc>0) THEN
        DO bb=1,nv4_proc
            iret = nf90_inq_varid(ncid0, out_an_list(bb), VarID)
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid0,VarID,out_an_data(:,i1:i2,j1:j2,bb),start=ibeg,count=icnt)
        ENDDO
        out_an_data(:,:,:,:) = 0.
    ENDIF

    ! User-defined outputs
    IF (nv4_user>0)  call an_user_stats()

    IF (level < 4) THEN ! Normal operation for levels 1-3
       ! Total water
       iret = nf90_inq_varid(ncid0,'q',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_rp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       ! Rain water mixing ratio and number concentration
       iret = nf90_inq_varid(ncid0,'r',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_rpp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       iret = nf90_inq_varid(ncid0,'n',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_npp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       ! Ice, snow, graupel and hail
       iret = nf90_inq_varid(ncid0,'i',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_rip(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       iret = nf90_inq_varid(ncid0,'ni',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_nip(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       iret = nf90_inq_varid(ncid0,'s',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_rsp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       iret = nf90_inq_varid(ncid0,'ns',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_nsp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       iret = nf90_inq_varid(ncid0,'g',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_rgp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       iret = nf90_inq_varid(ncid0,'ng',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_ngp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       iret = nf90_inq_varid(ncid0,'h',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_rhp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
       iret = nf90_inq_varid(ncid0,'nh',VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid0,VarID,a_nhp(:,i1:i2,j1:j2),start=ibeg,count=icnt)
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


       !  All bin-dependent outputs are calculated here

       ! Aerosol bin number concentration (regimes A and B)
       iret = nf90_inq_varid(ncid0,'B_Naa',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarId,a_naerop(:,i1:i2,j1:j2,in1a:fn2a),start=ibegsd,count=icntaea)
       iret = nf90_inq_varid(ncid0,'B_Nab',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_naerop(:,i1:i2,j1:j2,in2b:fn2b),start=ibegsd,count=icntaeb)

       ! Aerosol bin wet radius (regimes A and B)
       iret = nf90_inq_varid(ncid0,'B_Rwaa',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nbins,nspec+1,a_naerop,a_maerop,nlim,a_Rawet,1)
          iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in1a:fn2a),start=ibegsd,count=icntaea)
       END IF
       iret = nf90_inq_varid(ncid0,'B_Rwab',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nbins,nspec+1,a_naerop,a_maerop,nlim,a_Rawet,1)
          iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in2b:fn2b),start=ibegsd,count=icntaeb)
       END IF

       ! Cloud droplet bin number concentration (regimes A and B)
       iret = nf90_inq_varid(ncid0,'B_Nca',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_ncloudp(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
       iret = nf90_inq_varid(ncid0,'B_Ncb',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_ncloudp(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)

       ! Cloud droplet bin wet radius (regimes A and B)
       iret = nf90_inq_varid(ncid0,'B_Rwca',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(ncld,nspec+1,a_ncloudp,a_mcloudp,nlim,a_Rcwet,2)
          iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
       END IF
       iret = nf90_inq_varid(ncid0,'B_Rwcb',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(ncld,nspec+1,a_ncloudp,a_mcloudp,nlim,a_Rcwet,2)
          iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)
       END IF

       ! Rain drop bin number concentration
       iret = nf90_inq_varid(ncid0,'B_Npt',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_nprecpp(:,i1:i2,j1:j2,1:nprc),start=ibegsd,count=icntpr)

       ! Rain drop bin wet radius
       iret = nf90_inq_varid(ncid0,'B_Rwpt',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nprc,nspec+1,a_nprecpp,a_mprecpp,prlim,a_Rpwet,3)
          iret = nf90_put_var(ncid0,VarID,a_Rpwet(:,i1:i2,j1:j2,1:nprc),start=ibegsd,count=icntpr)
       END IF

       ! Ice bin number concentration (regimes A and B)
       iret = nf90_inq_varid(ncid0,'B_Nia',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_nicep(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
       iret = nf90_inq_varid(ncid0,'B_Nib',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_nicep(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)

       ! Ice bin wet radius (regimes A and B)
       iret = nf90_inq_varid(ncid0,'B_Rwia',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nice,nspec+1,a_nicep,a_micep,prlim,a_Riwet,4)
          iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
       END IF
       iret = nf90_inq_varid(ncid0,'B_Rwib',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nice,nspec+1,a_nicep,a_micep,prlim,a_Riwet,4)
          iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)
       END IF

       ! Snow bin number concentration
       iret = nf90_inq_varid(ncid0,'B_Nst',VarID)
       IF (iret==NF90_NOERR) &
          iret = nf90_put_var(ncid0,VarID,a_nsnowp(:,i1:i2,j1:j2,1:nsnw),start=ibegsd,count=icntsn)

       ! Snow bin wet radius
       iret = nf90_inq_varid(ncid0,'B_Rwst',VarID)
       IF (iret==NF90_NOERR) THEN
          CALL getBinRadius(nsnw,nspec+1,a_nsnowp,a_msnowp,prlim,a_Rswet,5)
          iret = nf90_put_var(ncid0,VarID,a_Rswet(:,i1:i2,j1:j2,1:nsnw),start=ibegsd,count=icntsn)
       END IF


       DO ee=1,nspec+1 ! Aerosol species and water
          ! Mixing ratios for each bin
          iret = nf90_inq_varid(ncid0,'B_'//TRIM(zspec(ee))//'aa',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = in1a,fn2a
                CALL binSpecMixrat('aerosol',ee,bb,a_Rawet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in1a:fn2a),start=ibegsd,count=icntaea)
          ENDIF
          iret = nf90_inq_varid(ncid0,'B_'//TRIM(zspec(ee))//'ab',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = in2b,fn2b
                CALL binSpecMixrat('aerosol',ee,bb,a_Rawet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in2b:fn2b),start=ibegsd,count=icntaeb)
          ENDIF

          iret = nf90_inq_varid(ncid0,'B_'//TRIM(zspec(ee))//'ca',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = inp2a,fnp2a
                CALL binSpecMixrat('cloud',ee,bb,a_Rcwet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
          ENDIF
          iret = nf90_inq_varid(ncid0,'B_'//TRIM(zspec(ee))//'cb',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = inp2b,fnp2b
                CALL binSpecMixrat('cloud',ee,bb,a_Rcwet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)
          ENDIF

          iret = nf90_inq_varid(ncid0,'B_'//TRIM(zspec(ee))//'pt',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = 1,nprc
                CALL binSpecMixrat('precp',ee,bb,a_Rpwet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rpwet(:,i1:i2,j1:j2,1:nprc),start=ibegsd,count=icntpr)
          ENDIF

          iret = nf90_inq_varid(ncid0,'B_'//TRIM(zspec(ee))//'ia',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = inp2a,fnp2a
                CALL binSpecMixrat('ice',ee,bb,a_Riwet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,inp2a:fnp2a),start=ibegsd,count=icntaeb)
          ENDIF
          iret = nf90_inq_varid(ncid0,'B_'//TRIM(zspec(ee))//'ib',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = inp2b,fnp2b
                CALL binSpecMixrat('ice',ee,bb,a_Riwet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,inp2b:fnp2b),start=ibegsd,count=icntaeb)
          ENDIF

          iret = nf90_inq_varid(ncid0,'B_'//TRIM(zspec(ee))//'st',VarID)
          IF (iret==NF90_NOERR) THEN
             DO bb = 1,nsnw
                CALL binSpecMixrat('snow',ee,bb,a_Rswet(:,:,:,bb))
            END DO
            iret = nf90_put_var(ncid0,VarID,a_Rswet(:,i1:i2,j1:j2,1:nsnw),start=ibegsd,count=icntsn)
          ENDIF
       ENDDO

    END IF

    if (myid==0) print "(/' ',12('-'),'   Record ',I3,' to: ',A80)",    &
         nrec0,fname

    iret  = nf90_sync(ncid0)
    nrec0 = nrec0+1

  end subroutine write_anal
  !
  ! ----------------------------------------------------------------------
  ! User-defined outputs (given in NAMELIST/user_an_list)
  !
  subroutine an_user_stats()
    use netcdf
    INTEGER :: i, iret, VarID
    REAL :: output(nzp,nxp,nyp)
    LOGICAL :: fail
    !
    DO i=1,nv4_user
        ! Is this active output (should be)?
        iret = nf90_inq_varid(ncid0,user_an_list(i),VarID)
        IF (iret/=NF90_NOERR) CYCLE
        ! Yes, so do the calculations
        SELECT CASE (user_an_list(i))
        CASE ('CCN')
           ! Level 3 CCN as an example of output
            output(:,:,:)=CCN
        CASE DEFAULT
            ! Pre-defined 3D SALSA outputs
            fail = calc_user_data(user_an_list(i),output)
            IF (fail) THEN
                WRITE(*,*)" Error: failed to calculate '"//TRIM(user_an_list(i))//"' for analysis output!"
                STOP
            ENDIF
        END SELECT
        ! Save
        iret = nf90_put_var(ncid0, VarID, output(:,3:nxp-2,3:nyp-2), &
                            start=(/1,1,1,nrec0/), count=(/nzp,nxp-4,nyp-4,1/))
    ENDDO


    CONTAINS
      ! -------------------------------------------------------------------------
      ! Produce 3D grid cell mean outputs based on user provided variable name:
      !     name='C_'//<species>//<bin>
      !       species: aerosol (SO4, OC, BC, ..) or gas name for mass concentration,
      !                R for mean radius and N for number concentration
      !       bin: aa,ab,at,ca,cb,ct,rt,ia,ib,it,st,gt (aerosol/cloud/rain/ice/snow/gas a-bins/b-bins/total)
      LOGICAL FUNCTION calc_user_data(short_name,res)
        USE mo_submctl, ONLY : find_gas_id, find_spec_id
        CHARACTER(LEN=7), INTENT(IN) :: short_name ! Variable name
        REAL, INTENT(out) :: res(nzp,nxp,nyp)      ! Output data
        !
        CHARACTER(LEN=7) :: spec, bin, nam
        INTEGER :: i, k
        !
        ! Return .true. if failed
        calc_user_data=.TRUE.
        !
        ! String length (must be at least 5, e.g. C_Naa, and start with 'C_')
        i=LEN_TRIM(short_name)
        IF (i<5 .OR. short_name(1:2)/='C_') RETURN
        !
        ! Output species & bin: aa/ab/at/ca/cb/ct/rt/ia/ib/it/st
        SELECT CASE (short_name(i:i))
        CASE('t')
            bin='ab'
        CASE('a')
            bin='a'
        CASE('b')
            bin='b'
        case default
            RETURN
        END SELECT
        SELECT case (short_name(i-1:i-1))
        CASE('a')
            spec='aerosol'
        CASE('c')
            spec='cloud'
        CASE('r','p')
            spec='precp'
        CASE('i')
            spec='ice'
        CASE('s')
            spec='snow'
        CASE('g')
            spec='gas'
        case default
            RETURN
        END SELECT
        !
        ! Species name, N for number, or R for radius
        nam=short_name(3:i-2)
        !
        ! Data
        IF (spec=='gas') THEN
            ! Just mass mixing ratio, no radius or number
            k=find_gas_id(nam)
            IF (k<1) RETURN
            res(:,:,:)=a_gaerop(:,:,:,k)
        ELSEIF (nam=='N  ') THEN
            CALL bulkNumc(spec,bin,res)
        ELSEIF (nam=='Rw ') THEN
            CALL meanRadius(spec,bin,res)
        ELSE
            k=find_spec_id(nam)
            IF (k<1) RETURN
            CALL bulkMixrat(k,spec,bin,res)
        ENDIF
        !
        ! All done
        calc_user_data = .FALSE.
      END FUNCTION calc_user_data

  end subroutine an_user_stats
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
    integer, allocatable, dimension(:) :: seed
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

    call random_seed(size=n)
    allocate (seed(n))
    call random_seed(get=seed)

    !
    ! Write fields
    !
    if (myid == 0) print "(/' ',12('-'),'   History write to: ',A30)" &
         ,hname
    open(10,file=trim(hname), form='unformatted')

    write(10) time,th00,umean,vmean,dtl,level,isgstyp,iradtyp,nzp,nxp,nyp,nscl
    write(10) xt, xm, yt, ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf

    write(10) a_ustar, a_tstar, a_rstar

    write(10) a_pexnr
    write(10) a_press
    write(10) a_theta
    write(10) a_edr

    write(10) a_up
    write(10) a_vp
    write(10) a_wp
    write(10) a_uc
    write(10) a_vc
    write(10) a_wc

    write(10) n
    write(10) seed

    do n=1,nscl
       call newsclr(n)
       write(10) a_sp
    end do

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
    integer, dimension(:), allocatable :: seed
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

       if (nxpx/=nxp .or. nypx/=nyp .or. nzpx/=nzp .or. lvlx/=level .or. nsclx/=nscl) then
          if (myid == 0) print *, nxp, nyp, nzp, nxpx, nypx, nzpx, lvlx, level, nsclx, nscl
          call appl_abort(-1)
       end if

       read (10) xt, xm, yt, ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf

       read (10) a_ustar, a_tstar, a_rstar

       read (10) a_pexnr
       read (10) a_press
       read (10) a_theta
       read (10) a_edr

       read (10) a_up
       read (10) a_vp
       read (10) a_wp
       read (10) a_uc
       read (10) a_vc
       read (10) a_wc

       read (10) n
       allocate(seed(n))
       read(10) seed
       call random_seed(put=seed)

       do n=1,nscl
          call newsclr(n)
          read (10) a_sp
       end do

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
          STOP 'bulkNumc: Invalid particle type'
    END SELECT

  END SUBROUTINE bulkNumc


  !
  ! -------------------------------------------------
  ! SUBROUTINE meanRadius
  ! Gets the mean wet (total number of species = nspec+1) radius for particles.
  !
  SUBROUTINE meanRadius(ipart,itype,rad,ibin)
    USE mo_submctl, ONLY : fnp2a,inp2b,fn2a,in2b,nlim,prlim
    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: ipart
    CHARACTER(len=*), INTENT(in) :: itype
    INTEGER, OPTIONAL, INTENT(in) :: ibin
    REAL, INTENT(out) :: rad(nzp,nxp,nyp)

    INTEGER :: istr,iend

    rad = 0.

    SELECT CASE(ipart)
    CASE('aerosol')
       IF (present(ibin)) THEN ! bin ibin
          istr = ibin
          iend = ibin
       ELSE IF (itype == 'ab') THEN ! Note: 1a, 2a and 2b combined
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
       IF (present(ibin)) THEN
          istr = ibin
          iend = ibin
       ELSE IF (itype == 'ab') THEN
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
       IF (present(ibin)) THEN
          istr = ibin
          iend = ibin
       ELSE
          istr = 1
          iend = nprc
       END IF
       CALL getRadius(istr,iend,nprc,nspec+1,a_nprecpp,a_mprecpp,prlim,rad,3)
    CASE('ice')
       IF (present(ibin)) THEN
          istr = ibin
          iend = ibin
       ELSE IF (itype == 'ab') THEN
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
       IF (present(ibin)) THEN
          istr = ibin
          iend = ibin
       ELSE
          istr = 1
          iend = nsnw
       END IF
       CALL getRadius(istr,iend,nsnw,nspec+1,a_nsnowp,a_msnowp,prlim,rad,5)
    CASE DEFAULT
          STOP 'meanRadius: Invalid particle type'
    END SELECT

  contains

   SUBROUTINE getRadius(zstr,zend,nn,n4,numc,mass,numlim,zrad,flag)
    USE mo_submctl, ONLY : calc_eff_radius
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
              tmp(:)=mass(k,i,j,bin:(n4-1)*nn+bin:nn)/numc(k,i,j,bin)
              rwet=rwet+calc_eff_radius(n4,tmp,flag)*numc(k,i,j,bin)
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
    USE mo_submctl, ONLY : calc_eff_radius
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
              tmp(:)=mass(k,i,j,bin:(n4-1)*nn+bin:nn)/numc(k,i,j,bin)
              zrad(k,i,j,bin)=calc_eff_radius(n4,tmp,flag)
            ENDIF
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE getBinRadius

end module grid

