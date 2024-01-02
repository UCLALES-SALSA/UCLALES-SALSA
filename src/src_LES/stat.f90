!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributsed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module stat

  use ncio, only : open_nc, define_nc
  use grid, only : level, lev_sb, maxn_list
  use util, only : get_avg3, get_cor3, get_var3, get_avg_ts, get_avg2dh, get_3rd3, &
                   HistDistr, get_zi_dmax, get_max_val, get_pustat_scalar, get_pustat_vector
  use defs, only : rowt, pi

  implicit none
  private

  integer, parameter :: nvar1 = 37,               &
                        nv1_ice = 22,             &
                        nv1_lvl4 = 5,             &
                        nv1_lvl5 = 12,            &
                        nvar2 = 98,               &
                        nv2_ice = 21,             &
                        nv2_lvl4 = 0,             &
                        nv2_lvl5 = 11,            &
                        nv2_hist = 2,             &
                        nvar3 = 25,               &
                        nv3_ice = 20,             &
                        nv3_lvl4 = 0,             &
                        nv3_lvl5 = 10

  integer, save      :: nrec1, nrec2, nrec3, ncid1, ncid2, ncid3
  real, save         :: fsttm, lsttm, nsmp = 0

  logical            :: sflg = .false.
  LOGICAL            :: csflg = .FALSE.
  real               :: ssam_intvl = 30.   ! statistical sampling interval
  real               :: savg_intvl = 1800. ! statistical averaging interval

  ! User-provided include and exclude list for outputs (analysis lists in module grid)
  CHARACTER(len=7), SAVE :: cs_include(maxn_list)='       ', cs_exclude(maxn_list)='       ', &
    ps_include(maxn_list)='       ', ps_exclude(maxn_list)='       ', &
    ts_include(maxn_list)='       ', ts_exclude(maxn_list)='       '
  ! User-selected process rate output lists (analysis ouputs in module grid)
  CHARACTER(LEN=7), SAVE :: out_mcrp_list(maxn_list)='       ', &
    out_cs_list(maxn_list)='       ',out_ps_list(maxn_list)='       ',out_ts_list(maxn_list)='       '
  ! User-defined outputs. Define variables in ncio.f90 and calculations in subroutines *_user_stats().
  CHARACTER(len=7), SAVE :: user_cs_list(maxn_list)='       ', user_ps_list(maxn_list)='       ', &
    user_ts_list(maxn_list)='       '

  ! Default variables
  character (len=7), save :: s1(nvar1)=(/                           &
       'time   ','cfl    ','maxdiv ','zi1_bar','zi2_bar','zi3_bar', & ! 1
       'vtke   ','sfcbflx','wmax   ','tsrf   ','ustar  ','shf_bar', & ! 7
       'lhf_bar','zi_bar ','lwp_bar','lwp_var','zc     ','zb     ', & !13
       'cfrac  ','lmax   ','albedo ','nccnt  ','zcmn   ','zbmn   ', & !19
       'ncloud ','wvp_bar','rwp_bar','prcp   ','nrain  ','nrcnt  ', & !25
       'prcp_bc','tkeint ','thl_int','prcc_bc','wvar_bc','Rcloud ', & !31
       'Rrain  '/), & !37

       s1_ice(nv1_ice) = (/ &
       'iwp_bar','Rice   ','nice   ','nicnt  ','iprcp  ',  & ! 1-5
       'swp_bar','Rsnow  ','nsnow  ','nscnt  ','sprcp  ',  & ! 6-10
       'gwp_bar','Rgra   ','ngra   ','ngcnt  ','gprcp  ',  & ! 11-15
       'hwp_bar','Rhail  ','nhail  ','nhcnt  ','hprcp  ',  & ! 16-20
       'SSi_max','thi_int'/), & ! 21

       s1_lvl4(nv1_lvl4) = (/       &
       'SS_max ','flx_aer','flx_iso','flx_mt ','u10    '/), & ! 1-5

       s1_lvl5(nv1_lvl5) = (/  &
       'iwp_bar','Rice   ','nice   ','nicnt  ','iprcp  ',  & ! 1-5
       'swp_bar','Rsnow  ','nsnow  ','nscnt  ','sprcp  ',  & ! 6-10
       'SSi_max','thi_int'/), & ! 11-12

        s2(nvar2)=(/                                                 &
        'time   ','zt     ','zm     ','dn0    ','u0     ','v0     ', & ! 1
        'fsttm  ','lsttm  ','nsmp   ','u      ','v      ','thl    ', & ! 7
        'p      ','u_2    ','v_2    ','w_2    ','thl_2  ','w_3    ', & ! 13
        'thl_3  ','tot_tw ','sfs_tw ','tot_uw ','sfs_uw ','tot_vw ', & ! 19
        'sfs_vw ','tot_ww ','sfs_ww ','km     ','kh     ','lmbd   ', & ! 25
        'lmbde  ','sfs_tke','sfs_boy','sfs_shr','boy_prd','shr_prd', & ! 31
        'trans  ','diss   ','dff_u  ','dff_v  ','dff_w  ','adv_u  ', & ! 37
        'adv_v  ','adv_w  ','prs_u  ','prs_v  ','prs_w  ','prd_uw ', & ! 43
        'storage','q      ','q_2    ','q_3    ','tot_qw ','sfs_qw ', & ! 49
        'rflx   ','rflx2  ','sflx   ','sflx2  ','l      ','l_2    ', & ! 55
        'l_3    ','tot_lw ','sed_lw ','cs1    ','cnt_cs1','w_cs1  ', & ! 61
        'tl_cs1 ','tv_cs1 ','rt_cs1 ','rc_cs1 ','wtl_cs1','wtv_cs1', & ! 67
        'wrt_cs1','cs2    ','cnt_cs2','w_cs2  ','tl_cs2 ','tv_cs2 ', & ! 73
        'rt_cs2 ','rc_cs2 ','wtl_cs2','wtv_cs2','wrt_cs2','rv     ', & ! 79
        'crate  ','frac_ic','Nc_ic  ','evap   ','rr     ','rrate  ', & ! 85
        'frac_ir','Nr_ir  ','sw_up  ','sw_down','lw_up  ','lw_down', & ! 91
        'Rc_ic  ','Rr_ir  '/), & ! 97, total 98

        s2_ice(nv2_ice)=(/ &
        'ri     ','Ni_ii  ','Ri_ii  ','frac_ii','irate  ', & ! 1
        'rs     ','Ns_is  ','Rs_is  ','frac_is','srate  ', & ! 6
        'rg     ','Ng_ig  ','Rg_ig  ','frac_ig','grate  ', & ! 11
        'rh     ','Nh_ih  ','Rh_ih  ','frac_ih','hrate  ','thi    '/), & ! 16

        s2_lvl4(nv2_lvl4), & ! Not used

        s2_lvl5(nv2_lvl5) = (/ &
        'ri     ','Ni_ii  ','Ri_ii  ','frac_ii','irate  ', & ! 1-5
        'rs     ','Ns_is  ','Rs_is  ','frac_is','srate  ', & ! 6-10
        'thi    '/), & ! 11

        ! **** Cloud droplet and ice histograms
        s2_CldHist(nv2_hist) = (/'P_hNca ','P_hNcb '/), &
        s2_IceHist(nv2_hist) = (/'P_hNia ','P_hNib '/), &

        ! Column statistics
        s3(nvar3)=(/ & ! Warm clouds
        'time   ','xt     ','xm     ','yt     ','ym     ', & ! 1-5
        'lwp    ','Nc_ic  ','Rc_ic  ','cprcp  ','nccnt  ', & ! 6-10
        'rwp    ','Nr_ir  ','Rr_ir  ','rprcp  ','nrcnt  ', & ! 11-15
        'zb     ','zc     ','zi1    ','lmax   ','albedo ', & ! 16-20
        'wmin   ','wmax   ','wbase  ','prcp_bc','prcc_bc'/), & ! 21-25

        s3_ice(nv3_ice)=(/ & ! Level 0 ice
        'iwp    ','Ni_ii  ','Ri_ii  ','iprcp  ','nicnt  ', & ! 1-5
        'swp    ','Ns_is  ','Rs_is  ','sprcp  ','nscnt  ', & ! 6-10
        'gwp    ','Ng_ig  ','Rg_ig  ','gprcp  ','ngcnt  ', & ! 11-15
        'hwp    ','Nh_ih  ','Rh_ih  ','hprcp  ','nhcnt  '/), & ! 16-20

        s3_lvl4(nv3_lvl4), & ! Not used

        s3_lvl5(nv3_lvl5) = (/ & ! Level 5 ice
        'iwp    ','Ni_ii  ','Ri_ii  ','iprcp  ','nicnt  ', & ! 1-5
        'swp    ','Ns_is  ','Rs_is  ','sprcp  ','nscnt  '/) ! 6-10

  real, save, allocatable   :: tke_sgs(:), tke_res(:), tke0(:), wtv_sgs(:),  &
       wtv_res(:), wrl_sgs(:), thvar(:), svctr(:,:), ssclr(:),               &
       ssclr_ice(:), svctr_ice(:,:),                                         &
       ! Additional ssclr and svctr for BULK SALSA output
       svctr_lvl4(:,:), ssclr_lvl4(:),                                       &
       svctr_lvl5(:,:), ssclr_lvl5(:),                                       &
       ssclr_rem(:),                                                         &
       ! Bin dependent SALSA outputs
       svctr_bin(:,:),                                                       &
       ! Cloud and ice histograms
       svctr_ch(:,:,:), svctr_ih(:,:,:),                                     &
       ! User-selected process rate outputs
       out_cs_data(:,:,:), out_ps_data(:,:), out_ts_data(:),                 &
       out_mcrp_data(:,:,:,:),                                               &
       ! User-defined outputs
       user_ps_data(:,:), user_ts_data(:),                                   &
       ! Removal rates for columns
       scs_rm(:,:,:)

  ! Case-dependent dimensions: mixing ratios and removal rates and binned mixing ratio profiles for active species
  INTEGER, SAVE :: nv1_rem=0, nv2_bin=0, nv2_bin_len(maxn_list)=0, &
    nv1_proc=0, nv2_proc=0, nv3_proc=0, out_mcrp_nout=0, &
    nv1_user=0, nv2_user=0, nv3_user=0
  CHARACTER (len=7), SAVE, ALLOCATABLE :: s1_rem(:)
  CHARACTER (len=7), SAVE :: s2_bin(maxn_list)='       '

  ! SALSA cloud, precipitation, ice and snow masks
  LOGICAL, allocatable :: cloudmask(:,:,:), rainmask(:,:,:), icemask(:,:,:), snowmask(:,:,:)
  REAL :: ris_min=1e-8 ! Ice and snow mass limit (prlim is for number)

  ! Calculating particle sizes (SB microphysics)
  REAL, DIMENSION(-1:4) :: x_min=1e-20, x_max=1., a_geo=(6./(pi*rowt))**(1./3.), b_geo=1./3., eps=1e-20
  REAL :: ri_min=1e-10, ni_min=1e-10 ! Mass and number limits for ice statistics

  public :: sflg, ssam_intvl, savg_intvl, statistics, init_stat, write_ps,   &
       updtst, close_stat, fill_scalar, &
       tke_sgs, sgsflxs, comp_tke, acc_removal, cs_rem_set, csflg, &
       les_rate_stats, mcrp_var_save, out_cs_list, out_ps_list, out_ts_list, &
       cs_include, cs_exclude, ps_include, ps_exclude, ts_include, ts_exclude, &
       out_mcrp_data, out_mcrp_list, out_mcrp_nout, user_cs_list, user_ps_list, &
       user_ts_list, setSBradius

contains
  !
  ! ---------------------------------------------------------------------
  ! INIT_STAT:  This routine initializes the statistical arrays which
  ! are user/problem defined.
  !
  subroutine init_stat(time, filprf, expnme, nzp)

    use grid, only : nxp, nyp, nprc, nsnw, nspec, iradtyp, &
        no_b_bins, no_prog_prc, no_prog_ice, no_prog_snw, &
        sed_aero, sed_cloud, sed_precp, sed_ice, sed_snow, out_an_list, nv4_proc, &
        user_an_list, nv4_user, ifSeaSpray, ifSeaVOC
    use mpi_interface, only : myid, ver, author, info
    use mo_submctl, only : fn1a,fn2a,fn2b,fnp2a,nout_cld,nout_ice,zspec

    character (len=80), intent (in) :: filprf, expnme
    integer, intent (in)            :: nzp
    real, intent (in)               :: time

    INTEGER :: i,ii,e,ee
    character (len=80) :: fname

    ! Local boolean arrays
    LOGICAL :: s1_bool(nvar1), s1_lvl4_bool(nv1_lvl4), s1_lvl5_bool(nv1_lvl5)
    LOGICAL :: s2_bool(nvar2), s2_lvl4_bool(nv2_lvl4), s2_lvl5_bool(nv2_lvl5)
    LOGICAL :: s3_bool(nvar3), s3_lvl4_bool(nv3_lvl4), s3_lvl5_bool(nv3_lvl5)
    LOGICAL :: s1_rem_bool(5*(nspec+1))
    LOGICAL :: s2_CldHist_bool(nv2_hist), s2_IceHist_bool(nv2_hist)
    LOGICAL :: tmp_bool(100)

    ! SALSA dimensions
    INTEGER, PARAMETER :: nv2_ndims = 6
    LOGICAL :: s2_dims_bool(nv2_ndims)=.FALSE.
    CHARACTER (len=7) :: s2_dims(nv2_ndims)=(/'B_Rd12a','B_Rd2ab','B_Rwprc','B_Rwsnw','P_hRc  ','P_hRi  '/)

    ! Local combined arrays
    LOGICAL, ALLOCATABLE :: s1bool(:), s2bool(:)
    CHARACTER (len=7), ALLOCATABLE :: s1total(:), s2total(:)

    allocate (wtv_sgs(nzp),wtv_res(nzp),wrl_sgs(nzp))
    allocate (tke_res(nzp),tke_sgs(nzp),tke0(nzp),thvar(nzp))

    wtv_sgs(:) = 0.
    wtv_res(:) = 0.
    wrl_sgs(:) = 0.
    tke_res(:) = 0.
    tke_sgs(:) = 0.
    tke0(:)    = 0.
    thvar(:)   = 0.

    ! Default outputs
    ALLOCATE ( ssclr(nvar1), svctr(nzp,nvar2) )
    svctr(:,:) = 0.
    ssclr(:)   = 0.
    s1_bool(:) = .TRUE.
    s2_bool(:) = .TRUE.
    s2_bool(22:27) = .FALSE. ! Wind fluxes
    s2_bool(30:31) = .FALSE. ! Length scales
    s2_bool(39:48) = .FALSE. ! Winds from diffusion and advection
    s2_bool(64:83) = .FALSE. ! Conditional sampling

    ! User selected process rate outputs: check, count and order inputs
    CALL test_user_vars(out_ts_list,maxn_list,nv1_proc)
    CALL test_user_vars(out_ps_list,maxn_list,nv2_proc)
    CALL test_user_vars(out_an_list,maxn_list,nv4_proc)
    IF (csflg) CALL test_user_vars(out_cs_list,maxn_list,nv3_proc)
    ! Allocate data arrays (analysis data located in module grid)
    ALLOCATE ( out_cs_data(nxp,nyp,nv3_proc), out_ps_data(nzp,nv2_proc), out_ts_data(nv1_proc) )
    out_cs_data(:,:,:) = 0.; out_ps_data(:,:) = 0.; out_ts_data(:) = 0.
    ! Find those process rate outputs that require raw 3D data from microphysics (SALSA or S&B)
    out_mcrp_nout=0
    CALL add_mcrp_list(out_mcrp_nout,out_mcrp_list,out_ts_list,maxn_list)
    CALL add_mcrp_list(out_mcrp_nout,out_mcrp_list,out_ps_list,maxn_list)
    CALL add_mcrp_list(out_mcrp_nout,out_mcrp_list,out_an_list,maxn_list)
    IF (csflg) CALL add_mcrp_list(out_mcrp_nout,out_mcrp_list,out_cs_list,maxn_list)
    ALLOCATE(out_mcrp_data(nzp,nxp,nyp,out_mcrp_nout))

    ! User-defined outputs: check, count and order inputs
    CALL test_user_vars(user_ts_list,maxn_list,nv1_user)
    CALL test_user_vars(user_ps_list,maxn_list,nv2_user)
    CALL test_user_vars(user_an_list,maxn_list,nv4_user)
    IF (csflg) CALL test_user_vars(user_cs_list,maxn_list,nv3_user)
    ! Allocate data arrays
    ALLOCATE ( user_ps_data(nzp,nv2_user), user_ts_data(nv1_user) )
    user_ps_data(:,:) = 0.; user_ts_data(:) = 0.

    ! Cloud and rain masks for all levels
    ALLOCATE ( cloudmask(nzp,nxp,nyp), rainmask(nzp,nxp,nyp) )

    IF ( level < 4 ) THEN
       ! Merge logical and name arrays
       ! a) Time series
       i=nvar1+nv1_proc+nv1_user+nv1_ice
       ALLOCATE( s1bool(i), s1total(i) )
       i=1; e=nvar1
       s1bool(i:e)=s1_bool; s1total(i:e)=s1
       IF (nv1_proc>0) THEN
          i=e+1; e=e+nv1_proc
          s1bool(i:e)=.TRUE.; s1total(i:e)=out_ts_list(1:nv1_proc)
       ENDIF
       IF (nv1_user>0) THEN
          i=e+1; e=e+nv1_user
          s1bool(i:e)=.TRUE.; s1total(i:e)=user_ts_list(1:nv1_user)
       ENDIF
       IF (level==0 .AND. lev_sb==5) THEN
          tmp_bool=.TRUE.
       ELSEIF (level==0 .AND. lev_sb==4) THEN
          tmp_bool=.TRUE.
          tmp_bool((/8,13,16,17,18,19,20/))=.FALSE.
       ELSE
          tmp_bool=.FALSE.
       ENDIF
       i=e+1; e=e+nv1_ice
       s1bool(i:e)=tmp_bool(1:nv1_ice); s1total(i:e)=s1_ice(1:nv1_ice)
       ! b) Profiles
       i=nvar2+nv2_proc+nv2_user+nv2_ice
       ALLOCATE( s2bool(i), s2total(i) )
       i=1; e=nvar2
       s2bool(i:e)=s2_bool; s2total(i:e)=s2
       IF (nv2_proc>0) THEN
          i=e+1; e=e+nv2_proc
          s2bool(i:e)=.TRUE.; s2total(i:e)=out_ps_list(1:nv2_proc)
       ENDIF
       IF (nv2_user>0) THEN
          i=e+1; e=e+nv2_user
          s2bool(i:e)=.TRUE.; s2total(i:e)=user_ps_list(1:nv2_user)
       ENDIF
       IF (level==0 .AND. lev_sb==5) THEN
          tmp_bool=.TRUE.
       ELSEIF (level==0 .AND. lev_sb==4) THEN
          tmp_bool=.TRUE.
          tmp_bool((/7,8,12,13,16,17,18,19,20/))=.FALSE.
       ELSE
          tmp_bool=.FALSE.
       ENDIF
       i=e+1; e=e+nv2_ice
       s2bool(i:e)=tmp_bool(1:nv2_ice); s2total(i:e)=s2_ice(1:nv2_ice)
       !
       IF (level==0) THEN
           ALLOCATE ( ssclr_ice(nv1_ice), svctr_ice(nzp,nv2_ice) )
           svctr_ice(:,:)=0.
           ssclr_ice(:)=0.
       ENDIF
    ELSE IF ( level >= 4 ) THEN
       ! Additional arrays for SALSA
       ! -- dimensions
       nv1_rem = 5*(nspec+1)  ! Removal with aerosol, cloud, rain, ice and snow
       ! -- allocate
       ALLOCATE ( ssclr_lvl4(nv1_lvl4), svctr_lvl4(nzp,nv2_lvl4) )
       ALLOCATE ( ssclr_lvl5(nv1_lvl5), svctr_lvl5(nzp,nv2_lvl5) )
       ALLOCATE ( ssclr_rem(nv1_rem) )
       ALLOCATE ( s1_rem(nv1_rem) )
       ALLOCATE ( icemask(nzp,nxp,nyp), snowmask(nzp,nxp,nyp) )
       ! -- reset data arrays
       ssclr_lvl4(:) = 0.; svctr_lvl4(:,:) = 0.
       ssclr_lvl5(:) = 0.; svctr_lvl5(:,:) = 0.
       ssclr_rem(:) = 0.

       ! Histograms
       IF (nout_cld>0) THEN
          ALLOCATE ( svctr_ch(nzp,nout_cld,2) )
          svctr_ch(:,:,:) = 0.
       ENDIF
       IF (level >4 .AND. nout_ice>0) THEN
            ALLOCATE ( svctr_ih(nzp,nout_ice,2) )
            svctr_ih(:,:,:) = 0.
       ENDIF

       ! Create a boolean array for items that are actually used

       ! 1) Standard outputs
       s1_lvl4_bool(:) = .TRUE.
       s1_lvl5_bool(:) = (level>4)
       IF (no_prog_ice) s1_lvl5_bool(1:5)=.FALSE.
       IF (no_prog_snw) s1_lvl5_bool(6:10)=.FALSE.
       s2_lvl4_bool(:) = .TRUE.
       s2_lvl5_bool(:) = (level>4)
       IF (no_prog_ice) s2_lvl5_bool(1:5)=.FALSE.
       IF (no_prog_snw) s2_lvl5_bool(6:10)=.FALSE.

       ! Surface emissions
       s1_lvl4_bool(2) = ifSeaSpray ! Aerosol
       s1_lvl4_bool(3) = ifSeaVOC ! Isoprene
       s1_lvl4_bool(4) = ifSeaVOC ! Monoterpene
       s1_lvl4_bool(5) = ifSeaSpray .OR. ifSeaVOC ! u10
       ssclr_lvl4(5) = -999. ! Spin-up means zero emissions but u10 is not zero (undefined)

       ! 2) Microphysical process rate statistics (both ts and ps)

       ! 3) Removal; dry, cloud, precipitation, ice and snow
       DO ee=1,nspec+1 ! Aerosol species and water
          ii = (ee-1)*5
          s1_rem(ii+1)='rm'//TRIM(zspec(ee))//'dr'
          s1_rem(ii+2)='rm'//TRIM(zspec(ee))//'cl'
          s1_rem(ii+3)='rm'//TRIM(zspec(ee))//'pr'
          s1_rem(ii+4)='rm'//TRIM(zspec(ee))//'ic'
          s1_rem(ii+5)='rm'//TRIM(zspec(ee))//'sn'
          s1_rem_bool(ii+1) = sed_aero
          s1_rem_bool(ii+2) = sed_cloud
          s1_rem_bool(ii+3) = sed_precp .AND. (.NOT. no_prog_prc)
          s1_rem_bool(ii+4) = sed_ice .AND. (level>4) .AND. (.NOT. no_prog_ice)
          s1_rem_bool(ii+5) = sed_snow .AND. (level>4) .AND. (.NOT. no_prog_snw)
       ENDDO

       ! SALSA dimensions
       !    s2_dims(nv2_ndims)=(/'B_Rd12a','B_Rd2ab','B_Rwprc','B_Rwsnw','P_hRc  ','P_hRi  '/)
       ! Bin dependent outputs
       IF ( ANY(INDEX(user_ps_list,'B_')>0) ) THEN
            s2_dims_bool(1:3)=.TRUE.
            s2_dims_bool(4)=(level>4)
       ENDIF

       ! Cloud and ice histrograms
       s2_CldHist_bool=.FALSE.
       IF (nout_cld>0) THEN
          s2_CldHist_bool(1) = .TRUE. ! A-bins
          s2_CldHist_bool(2) = .NOT. no_b_bins ! B-bins
          ! Add dimension
          s2_dims_bool(5)=.TRUE.
       END IF
       s2_IceHist_bool=.FALSE.
       IF (nout_ice>0 .AND. level>=5 .AND. .NOT. no_prog_ice) THEN
          s2_IceHist_bool(1) = .TRUE.
          s2_IceHist_bool(2) = .NOT. no_b_bins
          ! Add dimension
          s2_dims_bool(6)=.TRUE.
       END IF

        ! Merge logical and name arrays
        ! a) Time series
        i=nvar1+nv1_lvl4+nv1_lvl5+nv1_rem+nv1_proc+nv1_user
        ALLOCATE( s1bool(i), s1total(i) )
        i=1; e=nvar1
        s1bool(i:e)=s1_bool; s1total(i:e)=s1
        i=e+1; e=e+nv1_lvl4
        s1bool(i:e)=s1_lvl4_bool; s1total(i:e)=s1_lvl4
        i=e+1; e=e+nv1_lvl5
        s1bool(i:e)=s1_lvl5_bool; s1total(i:e)=s1_lvl5
        i=e+1; e=e+nv1_rem
        s1bool(i:e)=s1_rem_bool; s1total(i:e)=s1_rem
        IF (nv1_proc>0) THEN
            i=e+1; e=e+nv1_proc
            s1bool(i:e)=.TRUE.; s1total(i:e)=out_ts_list(1:nv1_proc)
        ENDIF
        IF (nv1_user>0) THEN
            i=e+1; e=e+nv1_user
            s1bool(i:e)=.TRUE.; s1total(i:e)=user_ts_list(1:nv1_user)
        ENDIF
        ! b) Profiles
        i=nvar2+nv2_lvl4+nv2_lvl5+2*nv2_hist+nv2_ndims+nv2_proc+nv2_user
        ALLOCATE( s2bool(i), s2total(i) )
        i=1; e=nvar2
        s2bool(i:e)=s2_bool; s2total(i:e)=s2
        i=e+1; e=e+nv2_lvl4
        s2bool(i:e)=s2_lvl4_bool; s2total(i:e)=s2_lvl4
        i=e+1; e=e+nv2_lvl5
        s2bool(i:e)=s2_lvl5_bool; s2total(i:e)=s2_lvl5
        i=e+1; e=e+nv2_hist
        s2bool(i:e)=s2_CldHist_bool; s2total(i:e)=s2_CldHist
        i=e+1; e=e+nv2_hist
        s2bool(i:e)=s2_IceHist_bool; s2total(i:e)=s2_IceHist
        i=e+1; e=e+nv2_ndims
        s2bool(i:e)=s2_dims_bool; s2total(i:e)=s2_dims
        IF (nv2_proc>0) THEN
            i=e+1; e=e+nv2_proc
            s2bool(i:e)=.TRUE.; s2total(i:e)=out_ps_list(1:nv2_proc)
        ENDIF
        IF (nv2_user>0) THEN
            i=e+1; e=e+nv2_user
            s2bool(i:e)=.TRUE.; s2total(i:e)=user_ps_list(1:nv2_user)
        ENDIF
        !
        ! Move bin dependent profile outputs from user_ps_list to s2_bin.
        ! Also, determine the output size i.e. the number of bins.
        nv2_bin=0
        IF (nv2_user>0) THEN
            DO i=1,nv2_user
                IF ('B_'==user_ps_list(i)(1:2)) THEN
                    nv2_bin=nv2_bin+1
                    s2_bin(nv2_bin)=user_ps_list(i)
                    ! Dimensions
                    IF (INDEX(user_ps_list(i),'aa')>0) THEN
                        ! Aerosol a-bin (1a+2a)
                        nv2_bin_len(nv2_bin)=fn2a
                    ELSEIF (INDEX(user_ps_list(i),'rt')>0 .OR. INDEX(user_ps_list(i),'pt')>0) THEN
                        ! Rain
                        nv2_bin_len(nv2_bin)=nprc
                    ELSEIF (INDEX(user_ps_list(i),'st')>0) THEN
                        ! Snow
                        nv2_bin_len(nv2_bin)=nsnw
                    ELSE
                        ! Any 2a or 2b bin
                        nv2_bin_len(nv2_bin)=fnp2a
                    ENDIF
                ELSEIF (nv2_bin>0) THEN
                    user_ps_list(i-nv2_bin)=user_ps_list(i)
                ENDIF
            ENDDO
            nv2_user=nv2_user-nv2_bin
            ! Allocate output data
            IF (nv2_bin>0) THEN
                i=SUM(nv2_bin_len(1:nv2_bin))
                ALLOCATE(svctr_bin(nzp,i))
            ENDIF
        ENDIF

    END IF ! If level >=4

    ! User options for including and excluding outputs
    CALL apply_user_outputs(.FALSE.,ps_exclude,maxn_list,s2total,s2bool,SIZE(s2bool))
    CALL apply_user_outputs(.TRUE.,ps_include,maxn_list,s2total,s2bool,SIZE(s2bool))
    CALL apply_user_outputs(.FALSE.,ts_exclude,maxn_list,s1total,s1bool,SIZE(s1bool))
    CALL apply_user_outputs(.TRUE.,ts_include,maxn_list,s1total,s1bool,SIZE(s1bool))

    if (myid == 0) THEN
        fname =  trim(filprf)//'.ts'
        print "(//' ',49('-')/,' ',/,'  Initializing: ',A20,'  N=',I3)",trim(fname),COUNT(s1bool)
        call open_nc( fname, expnme, time, (nxp-4)*(nyp-4), ncid1, nrec1, ver, author, info, par=.FALSE.)
        ! Juha: Modified for SALSA output
        call define_nc( ncid1, nrec1, COUNT(s1bool), PACK(s1Total,s1bool))
        print *, '   ...starting record: ', nrec1
    ENDIF

    if (myid == 0) then
        fname =  trim(filprf)//'.ps'
        print "(//' ',49('-')/,' ',/,'  Initializing: ',A20,'  N=',I3)",trim(fname),COUNT(s2bool)
        call open_nc( fname, expnme, time,(nxp-4)*(nyp-4), ncid2, nrec2, ver, author, info, par=.FALSE.)
        ! Juha: Modified due to SALSA output
        IF (level<4) THEN
            call define_nc( ncid2, nrec2, COUNT(s2bool), PACK(s2Total,s2bool), n1=nzp)
        ELSEIF (level<5) THEN
            call define_nc( ncid2, nrec2, COUNT(s2bool), PACK(s2Total,s2bool), n1=nzp,  &
                n1a=fn1a, n2a=fn2a-fn1a, n2b=fn2b-fn2a, nprc=nprc, nchist=nout_cld)
        ELSE
            call define_nc( ncid2, nrec2, COUNT(s2bool), PACK(s2Total,s2bool), n1=nzp,  &
                n1a=fn1a, n2a=fn2a-fn1a, n2b=fn2b-fn2a, nprc=nprc, nchist=nout_cld, &
                nsnw=nsnw, nihist=nout_ice)
        ENDIF
        print *, '   ...starting record: ', nrec2
    ENDIF


    ! Optional column statistics
    IF (csflg) THEN
        ! Allocate data arrays
        IF (nv1_rem>0) THEN
            ! Removal with aerosol, cloud, rain, ice and snow (the same as for ts)
            ALLOCATE( scs_rm(nxp,nyp,nv1_rem) )
            scs_rm=0. ! Set to zero (during spinup)
        END IF

        ! Logical arrays
        s3_bool=.TRUE.
        s3_bool(20) = (iradtyp==3) ! Albedo
        ! Level 4
        s3_lvl4_bool=(level>=4)
        ! Level 5
        s3_lvl5_bool=(level>=5)
        IF (no_prog_ice) s3_lvl5_bool(1:5)=.FALSE.
        IF (no_prog_snw) s3_lvl5_bool(6:10)=.FALSE.
        ! Level 0 ice
        tmp_bool=(level==0 .AND. lev_sb>=4)
        tmp_bool(16:20)=(level==0 .AND. lev_sb>=5) ! hail is for sb_lev=5

        ! Merge logical and name arrays (reuse s1 arrays)
        DEALLOCATE( s1bool, s1total )
        ! Default (mostly levels 3 & 4) parameters, specific level 4 and 5 parameters,
        ! removal statistics, process rate statistics and user outputs.
        ! Note: removal outputs are the same as those in the time series outputs!
        i = nvar3 + nv3_lvl4 + nv3_lvl5 + nv3_ice + nv1_rem+nv3_proc+nv3_user
        ALLOCATE( s1bool(i), s1total(i) )
        i=1; e=nvar3
        s1total(i:e)=s3; s1bool(i:e)=s3_bool
        IF (nv3_lvl4>0) THEN
            i=e+1; e=e+nv3_lvl4
            s1bool(i:e)=s3_lvl4_bool; s1total(i:e)=s3_lvl4
        ENDIF
        IF (nv3_lvl5>0) THEN
            i=e+1; e=e+nv3_lvl5
            s1bool(i:e)=s3_lvl5_bool; s1total(i:e)=s3_lvl5
        ENDIF
        IF (nv3_ice>0) THEN
            i=e+1; e=e+nv3_ice
            s1bool(i:e)=tmp_bool(1:nv3_ice); s1total(i:e)=s3_ice
        ENDIF
        IF (nv1_rem>0) THEN
            i=e+1; e=e+nv1_rem
            ! Logical and name arrays are the same as those for ts
            s1bool(i:e)=s1_rem_bool; s1total(i:e)=s1_rem
        ENDIF
        IF (nv3_proc>0) THEN
            i=e+1; e=e+nv3_proc
            s1bool(i:e)=.TRUE.; s1total(i:e)=out_cs_list(1:nv3_proc)
        ENDIF
        IF (nv3_user>0) THEN
            i=e+1; e=e+nv3_user
            s1bool(i:e)=.TRUE.; s1total(i:e)=user_cs_list(1:nv3_user)
        ENDIF

        ! User options for including and excluding outputs
        CALL apply_user_outputs(.FALSE.,cs_exclude,maxn_list,s1total,s1bool,SIZE(s1bool))
        CALL apply_user_outputs(.TRUE.,cs_include,maxn_list,s1total,s1bool,SIZE(s1bool))

        fname =  trim(filprf)//'.cs'
        if(myid == 0) print                                                  &
            "(//' ',49('-')/,' ',/,'  Initializing: ',A20,'  N=',I3)",trim(fname),COUNT(s1bool)
        call open_nc( fname, expnme, time,(nxp-4)*(nyp-4), ncid3, nrec3, ver, author, info)
        call define_nc( ncid3, nrec3, COUNT(s1bool), PACK(s1Total,s1bool), n2=nxp-4, n3=nyp-4)
        if (myid == 0) print *, '   ...starting record: ', nrec3
    ENDIF

  end subroutine init_stat

  ! Check user-selected process rate outputs
  subroutine test_user_vars(list,n,nout)
    USE ncio, ONLY : ncinfo
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    CHARACTER(LEN=7), INTENT(INOUT) :: list(n)
    INTEGER, INTENT(OUT) :: nout
    ! Local
    INTEGER :: i, j, k
    LOGICAL :: fail
    CHARACTER(LEN=80) :: info
    !
    fail=.FALSE.
    j = 0
    DO i=1,n
        IF (LEN(TRIM(list(i)))>0) THEN
            ! Not an empty list item
            j=j+1
            ! Remove possible empty items so that the real data will be continuous
            IF (j<i) THEN
                list(j)=list(i)
                list(i)='       '
            ENDIF
            ! Examine that this is a valid output (if not then ncinfo will fail)
            info=ncinfo(0, list(j))
            ! Check that there are no duplicates
            DO k=1,j-1
                IF (list(k)==list(j)) THEN
                    WRITE(*,*)" Error: duplicate item #",i," '"//TRIM(list(k))//"' in the output list!"
                    fail=.TRUE.
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    ! Stop if failures
    IF (fail) stop
    !
    ! Return the true number of list items
    nout=j
    !
  end subroutine test_user_vars

  ! Add microphysics-related outputs to list "base" from the common list "add"
  SUBROUTINE add_mcrp_list(j,base,add,n)
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(INOUT) :: j
    INTEGER, INTENT(IN) :: n
    CHARACTER(LEN=7), INTENT(INOUT) :: base(n)
    CHARACTER(LEN=7), INTENT(IN) :: add(n)
    ! Local
    INTEGER :: i, k
    LOGICAL :: dupl
    CHARACTER(LEN=4) :: short
    !
    DO i=1,n
        IF (LEN(TRIM(add(i)))>0) THEN
            ! Do not include LES items (any prefix used in the calls from t_step/step.f90)
            short=add(i)(1:4)
            IF (short=='srfc' .OR. short=='diff' .OR. short=='forc' .OR. short=='mcrp' .OR. &
                short=='advf' .OR. short=='nudg' .OR. (level>3 .AND. (short=='diag' .OR. short=='sedi')) ) CYCLE
            ! Check that there are no duplicates
            dupl=.FALSE.
            DO k=1,j
                IF (base(k)==add(i)) dupl=.true.
            ENDDO
            IF (.not.dupl) THEN
                j=j+1
                base(j)=add(i)
            ENDIF
        ENDIF
    ENDDO
  END SUBROUTINE add_mcrp_list

  ! Set ouputs on or off based on list from user
  SUBROUTINE apply_user_outputs(set,list,n,out_list,out_flags,m)
    IMPLICIT NONE
    ! Inputs
    LOGICAL, INTENT(IN) :: set ! Set true or false
    INTEGER, INTENT(IN) :: n, m ! Maximum number of variables in the lists
    CHARACTER(LEN=7), INTENT(IN) :: list(n) ! List of variables to be included or excluded
    CHARACTER(LEN=7), INTENT(IN) :: out_list(m) ! Current output list
    LOGICAL, INTENT(INOUT) :: out_flags(m) ! Current output flags
    ! Local
    INTEGER :: i, j
    LOGICAL :: found
    !
    DO i=1,n ! Include/exclude list
        IF (LEN_TRIM(list(i))>0) THEN
            found=.FALSE.
            DO j=1,m ! Current output list
                IF ( list(i)==out_list(j) ) THEN
                    out_flags(j)=set
                    found=.TRUE.
                ENDIF
            ENDDO
            IF (.NOT.found) WRITE(*,*) 'Warning: '//TRIM(list(i))//' not found!',i
        ENDIF
    ENDDO
  END SUBROUTINE apply_user_outputs

  !
  ! ---------------------------------------------------------------------
  ! Subroutine Statistics:  This subroutine is the statistics driver
  ! it calls various other subroutines to compute and accumulate
  ! statistical quantities.  These are stored in two arrays:  SVCTR,
  ! and SSCLR (which accumulate scalar and vector statistics respectively
  !
  ! Modified for level 4
  ! Juha Tonttila, FMI, 2014
  !
  ! Modified for level 5
  ! Jaakko Ahola, FMI, 2016
  subroutine statistics(time)
    use mo_submctl, only : nlim, prlim
    use grid, only : a_up, a_vp, a_wp, a_rc, a_theta, a_temp, a_rv           &
         , a_rp, a_tp, a_press, nxp, nyp, nzp, dzm, dzt, zm, zt, th00, umean            &
         , vmean, dn0, a_dn, cldin, precip, a_rpp, a_npp, CCN, iradtyp, a_rflx, a_sflx  &
         , a_fus, a_fds, a_fuir, a_fdir, albedo, a_srp, a_snrp, a_ncloudp, a_ri, a_srs
    USE defs, ONLY : cp, alvi

    real, intent (in) :: time

    real :: rxt(nzp,nxp,nyp), rxl(nzp,nxp,nyp), rxv(nzp,nxp,nyp), rnt(nzp,nxp,nyp)
    REAL :: xrpp(nzp,nxp,nyp), xnpp(nzp,nxp,nyp), thl(nzp,nxp,nyp)

    SELECT CASE(level)
       CASE (0,3)
          rxt = a_rp ! Total water (vapor + condensed water and ice) = q
          rxl = a_rc ! Cloud water (+aerosol), but no precipitation or ice
          rxv = a_rv ! Water vapor
          xrpp = a_rpp ! Rain water
          xnpp = a_npp ! Rain number
          thl = a_tp ! Liquid water potential temperature
          IF (level==0) WHERE(a_temp>0.) thl = a_tp + (a_theta/a_temp)*alvi/cp*a_ri
       CASE(4)
          rxt = a_rp + a_rc + a_srp
          rxl = a_rc
          rxv = a_rp
          xrpp = a_srp
          xnpp = a_snrp
          thl = a_tp
       CASE(5)
          rxt = a_rp + a_rc + a_srp + a_ri + a_srs
          rxl = a_rc
          rxv = a_rp
          xrpp = a_srp
          xnpp = a_snrp
          thl = a_tp
          WHERE(a_temp>0) thl = a_tp + (a_theta/a_temp)*alvi/cp*(a_ri + a_srs)
       CASE DEFAULT
          rxt = a_rp
          rxl = a_rc
          rxv = a_rv
          xrpp = 0.
          xnpp = 0.
          thl = a_tp
    END SELECT

    ! CDNC
    rnt = CCN
    IF (level>3) rnt = SUM(a_ncloudp,DIM=4)

    ! Cloud and rain masks
    IF (level<4) THEN
        cloudmask(:,:,:) = ( rxl(:,:,:) > 1.e-5 )
        rainmask(:,:,:) = ( xrpp(:,:,:) > 1.e-8 )
    ELSE
        cloudmask(:,:,:) = ( rnt(:,:,:) > nlim .AND. rxl(:,:,:) > 1.e-5 )
        rainmask(:,:,:) = ( xnpp(:,:,:) > prlim .AND. xrpp(:,:,:) > 1.e-8 )
    ENDIF

    if (nsmp == 0.) fsttm = time
    nsmp=nsmp+1.
    ssclr(14:nvar1) = -999.
    !
    ! profile statistics
    !
    call accum_stat(nzp, nxp, nyp, a_up, a_vp, a_wp, thl, a_press, umean &
         ,vmean, th00)
    if (iradtyp == 3) then
       call accum_rad(nzp, nxp, nyp, a_rflx, sflx=a_sflx, sup=a_fus, sdwn=a_fds, &
         irup=a_fuir, irdwn=a_fdir, alb=albedo)
    elseif (iradtyp > 0) then
       call accum_rad(nzp, nxp, nyp, a_rflx)
    end if
    call accum_cld(nzp, nxp, nyp, th00, a_wp, a_theta, thl, rxl, rxt, rxv, cldin, rnt, xrpp, xnpp, precip)
    if (level ==0 .AND. (lev_sb==4 .OR. lev_sb==5)) call accum_ice(nzp, nxp, nyp)
    if (level >=4)  call accum_lvl4(nzp, nxp, nyp)
    if (level >=5)  call accum_lvl5(nzp, nxp, nyp)
    IF (nv2_user>0) call ps_user_stats()
    IF (nv2_bin>0)  call ps_user_bin_stats()
    !
    ! scalar statistics
    !
    call set_ts(nzp, nxp, nyp, a_wp, a_theta, thl, dn0, zt,zm,dzt,th00,time)
    CALL ts_cld(nzp, nxp, nyp, a_wp, a_dn, zm, zt, dzt, rxl, rxt, rxv, rnt, xrpp, xnpp, precip,cldin)
    IF ( level ==0 .AND. (lev_sb==4 .OR. lev_sb==5) ) CALL ts_ice(nzp, nxp, nyp)
    IF ( level >=4 ) CALL ts_lvl4(nzp, nxp, nyp)
    IF ( level >=5 ) CALL ts_lvl5(nzp, nxp, nyp)
    IF ( nv1_user>0 ) CALL ts_user_stats()

    call write_ts

    !
    ! Column statistics
    !
    IF (csflg) THEN
        ! Warm cloud statistics (rnt=CDNC, xrpp and xnpp from above)
        CALL set_cs_warm(nzp,nxp,nyp,a_rc,rnt,xrpp,xnpp,a_theta,a_wp,cldin,precip,a_dn,zt,dzm)

        ! Ice cloud statistics
        IF (level==5) THEN
            CALL set_cs_lvl5(nzp,nxp,nyp)
        ELSEIF (level==0 .AND. (lev_sb==4 .OR. lev_sb==5)) THEN
            CALL set_cs_ice(nzp,nxp,nyp)
        ENDIF

        ! User-defined outputs
        IF (nv3_user>0) CALL cs_user_stats()

        ! All other outputs
        CALL set_cs_other(time)
    ENDIF

  end subroutine statistics
  !
  ! -----------------------------------------------------------------------
  ! subroutines set_cs_warm, cs_rem_set, cs_rem_save and set_cs_other:
  ! write (and compute) column average statistics
  !
  ! Removal statistics (level>3): calculate values for further use
  SUBROUTINE cs_rem_set(n2,n3,n4,raer,rcld,rprc,rice,rsnw)

    USE grid, ONLY : nbins, ncld, nprc, nice, nsnw, nspec
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n2,n3,n4   ! Grid dimensions
    REAL, INTENT(in) :: raer(n2,n3,n4*nbins), & ! Removal arrays
                               rcld(n2,n3,n4*ncld), &
                               rprc(n2,n3,n4*nprc), &
                               rice(n2,n3,n4*nice), &
                               rsnw(n2,n3,n4*nsnw)

    INTEGER :: si, tt

    IF (.NOT.csflg) RETURN

    ! Calculate all removal fluxes and save those to scs_rm for later use
    DO si = 1,nspec+1 ! Aerosol species and water
        tt=(si-1)*5
        ! Removal by sedimentation of aerosol
        scs_rm(:,:,tt+1) = SUM(raer(:,:,(si-1)*nbins+1:si*nbins),DIM=3)
        ! Removal by sedimentation of cloud droplets
        scs_rm(:,:,tt+2) = SUM(rcld(:,:,(si-1)*ncld+1:si*ncld),DIM=3)
        ! Removal by precipitation
        scs_rm(:,:,tt+3) = SUM(rprc(:,:,(si-1)*nprc+1:si*nprc),DIM=3)
        IF (level>4) THEN
            ! Removal by sedimentation of ice particles
            scs_rm(:,:,tt+4) = SUM(rice(:,:,(si-1)*nice+1:si*nice),DIM=3)
            ! Removal by snow
            scs_rm(:,:,tt+5) = SUM(rsnw(:,:,(si-1)*nsnw+1:si*nsnw),DIM=3)
        ENDIF
    ENDDO

  END SUBROUTINE cs_rem_set
  !
  ! Calculate warm cloud statistics
  subroutine set_cs_warm(n1,n2,n3,rc,nc,rp,np,th,w,crate,rrate,dn,zt,dzm)

    use netcdf
    use grid, only : meanRadius

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: rc(n1,n2,n3),nc(n1,n2,n3),rp(n1,n2,n3),np(n1,n2,n3),th(n1,n2,n3)
    real, intent(in)    :: w(n1,n2,n3),crate(n1,n2,n3),rrate(n1,n2,n3),dn(n1,n2,n3),zt(n1),dzm(n1)
    REAL :: zb(n2,n3), zc(n2,n3), th1(n2,n3), prcp_bc(n2,n3), prcc_bc(n2,n3), wbase(n2,n3)
    REAL :: a(n2,n3), a1(n1,n2,n3)
    integer :: i, j, k, iret, VarID, ibase, itop
    real    :: sval, dmy

    ! Calculate stats
    zb=-999.    ! Cloud base height (m)
    zc=-999.    ! Cloud top height (m)
    th1=0.      ! Height of the maximum theta gradient (m)
    prcp_bc=-999. ! Rain precipitation rate at the cloud base (W/m^2)
    prcc_bc=-999. ! Cloud precipitation rate at the cloud base (W/m^2)
    wbase=-999.   ! Vertical velocity -||- (m/s)
    do j=3,n3-2
       do i=3,n2-2
          ibase=-1
          itop=-1
          sval = 0.
          do k=2,n1
             ! The first and last cloudy grid cell
             IF (cloudmask(k,i,j)) THEN
                IF (ibase<0) ibase=k
                itop=k
             END IF
             !
             ! Height of the maximum theta gradient
             if (k<=n1-5) then
                dmy = (th(k+1,i,j)-th(k,i,j))*dzm(k)
                if (dmy > sval ) then
                   sval = dmy
                   th1(i,j) = zt(k)
                end if
            ENDIF
          enddo
          IF (ibase>0) THEN
            ! Cloud base and top
            zb(i,j)=zt(ibase)
            zc(i,j)=zt(itop)
            ! Cloud base vertical velocity and precipitation (taken from the level just below cloud base)
            wbase(i,j) = w(max(2,ibase-1),i,j)
            prcp_bc(i,j) = rrate(max(2,ibase-1),i,j)
            prcc_bc(i,j) = crate(max(2,ibase-1),i,j)
          ENDIF
       end do
    end do

    ! Save the data

    iret = nf90_inq_varid(ncid3,'zb',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, zb(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'zc',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, zc(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'zi1',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, th1(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'wbase',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, wbase(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'prcp_bc',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, prcp_bc(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'prcc_bc',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, prcc_bc(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'lmax',VarID)
    IF (iret==NF90_NOERR) THEN
        a=MAXVAL(rc,DIM=1)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF

    iret = nf90_inq_varid(ncid3,'wmax',VarID)
    IF (iret==NF90_NOERR) THEN
        a=MAXVAL(w,DIM=1)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF

    iret = nf90_inq_varid(ncid3,'wmin',VarID)
    IF (iret==NF90_NOERR) THEN
        a=MINVAL(w,DIM=1)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF

    ! Liquid water path
    iret = nf90_inq_varid(ncid3,'lwp',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,rc,a,dens=dn)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Radius, conditional average
    iret = nf90_inq_varid(ncid3,'Rc_ic',VarID)
    IF (iret==NF90_NOERR) THEN
        ! Calculate wet radius
        IF (level>3) THEN
            CALL meanRadius('cloud','ab',a1)
        ELSE
            CALL getSBradius(n1,n2,n3,nc,rc,-1,a1)
        ENDIF
        CALL get_avg_cs(n1,n2,n3,a1,a,cond=cloudmask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Number concentration, conditional average
    iret = nf90_inq_varid(ncid3,'Nc_ic',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,nc,a,cond=cloudmask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Total number of cloudy grid cells
    iret = nf90_inq_varid(ncid3,'nccnt',VarID)
    IF (iret==NF90_NOERR) THEN
        a(:,:)=COUNT( cloudmask, DIM=1 )
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Surface precipitation rate
    iret = nf90_inq_varid(ncid3,'cprcp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, crate(2,3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    ! Rain water path
    iret = nf90_inq_varid(ncid3,'rwp',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,rp,a,dens=dn)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Radius, conditional average
    iret = nf90_inq_varid(ncid3,'Rr_ir',VarID)
    IF (iret==NF90_NOERR) THEN
        ! Calculate wet radius
        IF (level>3) THEN
            CALL meanRadius('precp','ab',a1)
        ELSE
            CALL getSBradius(n1,n2,n3,np,rp,0,a1)
        ENDIF
        CALL get_avg_cs(n1,n2,n3,a1,a,cond=rainmask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Number concentration, conditional average
    iret = nf90_inq_varid(ncid3,'Nr_ir',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,np,a,cond=rainmask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Total number of rainy grid cells
    iret = nf90_inq_varid(ncid3,'nrcnt',VarID)
    IF (iret==NF90_NOERR) THEN
        a(:,:)=COUNT( rainmask, DIM=1 )
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Surface precipitation rate
    iret = nf90_inq_varid(ncid3,'rprcp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, rrate(2,3:n2-2,3:n3-2), start=(/1,1,nrec3/))

  end subroutine set_cs_warm

  ! Calculate cold cloud statistics
  ! a) SALSA level 5
  subroutine set_cs_lvl5(n1,n2,n3)
    USE grid, ONLY : meanRadius, bulkNumc, a_dn, a_ri, a_srs, icein, snowin
    use netcdf
    integer, intent (in) :: n1,n2,n3
    real    :: a(n2,n3), a1(n1,n2,n3)
    INTEGER :: iret,VarID

    ! Ice water path
    iret = nf90_inq_varid(ncid3,'iwp',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,a_ri,a,dens=a_dn)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Radius, conditional average
    iret = nf90_inq_varid(ncid3,'Ri_ii',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL meanRadius('ice','ab',a1)
        a=-999.
        CALL get_avg_cs(n1,n2,n3,a1,a,cond=icemask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Number concentration, conditional average
    iret = nf90_inq_varid(ncid3,'Ni_ii',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL bulkNumc('ice','ab',a1)
        a=-999.
        CALL get_avg_cs(n1,n2,n3,a1,a,cond=icemask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Total number of icy grid cells
    iret = nf90_inq_varid(ncid3,'nicnt',VarID)
    IF (iret==NF90_NOERR) THEN
        a(:,:)=COUNT( icemask, DIM=1 )
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Surface precipitation rate
    iret = nf90_inq_varid(ncid3,'iprcp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, icein(2,3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    ! The same for snow
    iret = nf90_inq_varid(ncid3,'swp',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,a_srs,a,dens=a_dn)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Radius
    iret = nf90_inq_varid(ncid3,'Rs_is',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL meanRadius('snow','ab',a1)
        a=-999.
        CALL get_avg_cs(n1,n2,n3,a1,a,cond=snowmask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Number concentration
    iret = nf90_inq_varid(ncid3,'Ns_is',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL bulkNumc('snow','ab',a1)
        a=-999.
        CALL get_avg_cs(n1,n2,n3,a1,a,cond=snowmask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Total number
    iret = nf90_inq_varid(ncid3,'nscnt',VarID)
    IF (iret==NF90_NOERR) THEN
        a(:,:)=COUNT( snowmask, DIM=1 )
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Surface precipitation rate
    iret = nf90_inq_varid(ncid3,'sprcp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, snowin(2,3:n2-2,3:n3-2), start=(/1,1,nrec3/))

  end subroutine set_cs_lvl5

  ! b) Level 0 (SB 4 and 5) ice microphysics
  SUBROUTINE set_cs_ice(n1,n2,n3)
    USE grid, ONLY : a_rip, a_nip, a_rsp, a_nsp, a_rgp, a_ngp, a_rhp, a_nhp, &
        a_dn, icein, snowin, grin, hailin
    use netcdf
    integer, intent (in) :: n1,n2,n3
    real    :: a(n2,n3), a1(n1,n2,n3)
    INTEGER :: iret,VarID
    LOGICAL :: mask(n1,n2,n3)

    ! Ice water path
    iret = nf90_inq_varid(ncid3,'iwp',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,a_rip,a,dens=a_dn)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Mask for number concentration and radius
    mask = (a_nip > ni_min .OR. a_rip > ri_min)
    ! Radius, conditional average
    iret = nf90_inq_varid(ncid3,'Ri_ii',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL getSBradius(n1,n2,n3,a_nip,a_rip,1,a1)
        CALL get_avg_cs(n1,n2,n3,a1,a,cond=mask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Number concentration, conditional average
    iret = nf90_inq_varid(ncid3,'Ni_ii',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,a_nip,a,cond=mask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Total number of icy grid cells
    iret = nf90_inq_varid(ncid3,'nicnt',VarID)
    IF (iret==NF90_NOERR) THEN
        a(:,:)=COUNT( mask, DIM=1 )
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    ! Surface precipitation rate
    iret = nf90_inq_varid(ncid3,'iprcp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, icein(2,3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    ! The same for snow
    iret = nf90_inq_varid(ncid3,'swp',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,a_rsp,a,dens=a_dn)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    mask = (a_nsp > ni_min .OR. a_rsp > ri_min)
    iret = nf90_inq_varid(ncid3,'Rs_is',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL getSBradius(n1,n2,n3,a_nsp,a_rsp,2,a1)
        CALL get_avg_cs(n1,n2,n3,a1,a,cond=mask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    iret = nf90_inq_varid(ncid3,'Ns_is',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,a_nsp,a,cond=mask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    iret = nf90_inq_varid(ncid3,'nscnt',VarID)
    IF (iret==NF90_NOERR) THEN
        a(:,:)=COUNT( mask, DIM=1 )
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    iret = nf90_inq_varid(ncid3,'sprcp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, snowin(2,3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    ! The same for graupel
    iret = nf90_inq_varid(ncid3,'gwp',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,a_rgp,a,dens=a_dn)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    mask = (a_ngp > ni_min .OR. a_rgp > ri_min)
    iret = nf90_inq_varid(ncid3,'Rg_ig',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL getSBradius(n1,n2,n3,a_ngp,a_rgp,3,a1)
        CALL get_avg_cs(n1,n2,n3,a1,a,cond=mask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    iret = nf90_inq_varid(ncid3,'Ng_ig',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,a_ngp,a,cond=mask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    iret = nf90_inq_varid(ncid3,'ngcnt',VarID)
    IF (iret==NF90_NOERR) THEN
        a(:,:)=COUNT( mask, DIM=1 )
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    iret = nf90_inq_varid(ncid3,'gprcp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, grin(2,3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    ! The same for hail
    iret = nf90_inq_varid(ncid3,'hwp',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,a_rhp,a,dens=a_dn)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    mask = (a_nhp > ni_min .OR. a_rhp > ri_min)
    iret = nf90_inq_varid(ncid3,'Rh_ih',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL getSBradius(n1,n2,n3,a_nhp,a_rhp,4,a1)
        CALL get_avg_cs(n1,n2,n3,a1,a,cond=mask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    iret = nf90_inq_varid(ncid3,'Nh_ih',VarID)
    IF (iret==NF90_NOERR) THEN
        CALL get_avg_cs(n1,n2,n3,a_nhp,a,cond=mask)
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    iret = nf90_inq_varid(ncid3,'nhcnt',VarID)
    IF (iret==NF90_NOERR) THEN
        a(:,:)=COUNT( mask, DIM=1 )
        iret = nf90_put_var(ncid3, VarID, a(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
    ENDIF
    iret = nf90_inq_varid(ncid3,'hprcp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, hailin(2,3:n2-2,3:n3-2), start=(/1,1,nrec3/))

  END SUBROUTINE set_cs_ice

  ! User-defined outputs (given in NAMELIST/user_cs_list)
  subroutine cs_user_stats()
    use netcdf
    use grid, ONLY : CCN, nzp, nxp, nyp, a_dn, a_rv, a_rp, a_rsl, a_rsi, a_temp
    INTEGER :: ii, iret, VarID
    REAL :: output(nxp,nyp), a(nzp,nxp,nyp)
    LOGICAL :: fail, mask(nzp,nxp,nyp), mass
    !
    DO ii=1,nv3_user
        ! Is this active output (should be)?
        iret = nf90_inq_varid(ncid3,user_cs_list(ii),VarID)
        IF (iret/=NF90_NOERR) CYCLE
        ! Yes, so do the calculations
        SELECT CASE (user_cs_list(ii))
        CASE ('CCN')
            ! Level 3 CCN as an example of output
            output(:,:)=CCN
        CASE ('SS_max')
            ! Maximum supersaturation (%) over water
            a=0.
            IF (level<4) THEN
                WHERE(a_rsl>1e-10) a=a_rv/a_rsl
            ELSE
                WHERE(a_rsl>1e-10) a=a_rp/a_rsl
            ENDIF
            output(:,:)=(MAXVAL(a,DIM=1)-1.0)*100.
        CASE ('SSi_max')
            ! Maximum supersaturation (%) over ice
            a=0.
            IF (level<4) THEN
                WHERE(a_rsi>1e-10) a=a_rv/a_rsi
            ELSE
                WHERE(a_rsi>1e-10) a=a_rp/a_rsi
            ENDIF
            output(:,:)=(MAXVAL(a,DIM=1)-1.0)*100.
        CASE ('T_min')
            ! Minimum absolute temperature (K)
            output(:,:)=MINVAL(a_temp,DIM=1)
        CASE DEFAULT
            ! Pre-defined SALSA outputs
            fail = calc_user_data(user_cs_list(ii),a,mask,is_mass=mass)
            IF (fail) THEN
                WRITE(*,*)" Error: failed to calculate '"//TRIM(user_cs_list(ii))//"' for cs output!"
                STOP
            ENDIF
            ! Calculate vertical integral when mass concentration, otherwise mean
            IF (mass) THEN
                CALL get_avg_cs(nzp,nxp,nyp,a,output,cond=mask,dens=a_dn)
            ELSE
                CALL get_avg_cs(nzp,nxp,nyp,a,output,cond=mask)
            ENDIF
        END SELECT
        ! Save
        iret = nf90_put_var(ncid3, VarID, output(3:nxp-2,3:nyp-2), start=(/1,1,nrec3/))
    ENDDO
  end subroutine cs_user_stats

  ! Other column statistics
  subroutine set_cs_other(time)
    use netcdf
    USE grid, ONLY : nxp, nyp, xt, xm, yt, ym, albedo
    ! Inputs
    REAL, INTENT(IN) :: time
    ! Local variables
    integer :: i, iret, VarID

    ! Coordinates
    if (nrec3 == 1) then
       iret = nf90_inq_varid(ncid3, 'xt', VarID)
       iret = nf90_put_var(ncid3, VarID, xt(3:nxp-2), start = (/nrec3/))
       iret = nf90_inq_varid(ncid3, 'yt', VarID)
       iret = nf90_put_var(ncid3, VarID, yt(3:nyp-2), start = (/nrec3/))
       iret = nf90_inq_varid(ncid3, 'xm', VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, xm(3:nxp-2), start = (/nrec3/))
       iret = nf90_inq_varid(ncid3, 'ym', VarID)
       IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, ym(3:nyp-2), start = (/nrec3/))
    END IF

    ! Time
    iret = nf90_inq_varid(ncid3,'time',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, time, start=(/nrec3/))

    ! Albedo
    iret = nf90_inq_varid(ncid3,'albedo',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, albedo(3:nxp-2,3:nyp-2), start=(/1,1,nrec3/))

    ! Level 4 and 5 removal/precipitation statistics (variable names are the same as those for ts)
    DO i = 1,nv1_rem
        iret = nf90_inq_varid(ncid3,s1_rem(i),VarID)
        IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, scs_rm(3:nxp-2,3:nyp-2,i), start=(/1,1,nrec3/))
    ENDDO

    ! User-selected process rate outputs
    DO i = 1,nv3_proc
        iret = nf90_inq_varid(ncid3,out_cs_list(i),VarID)
        IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, out_cs_data(3:nxp-2,3:nyp-2,i), start=(/1,1,nrec3/))
    END DO

    iret = nf90_sync(ncid3)

    ! Update counter
    nrec3 = nrec3 + 1
  END subroutine set_cs_other
  !
  !---------------------------------------------------------------------
  ! GET_AVG_CS: gets average field over vertical dimension.
  !
  SUBROUTINE get_avg_cs(n1,n2,n3,a,res,cond,dens)
    USE grid, ONLY : zm ! Momentum points
    integer, intent (in) :: n1, n2, n3
    real, intent (in)   :: a(n1,n2,n3)
    real, intent (OUT)  :: res(n2,n3)
    logical, optional, intent(in) :: cond(n1,n2,n3) ! A flag for grid cells
    real, optional, intent(in) :: dens(n1,n2,n3) ! Air density
    integer :: i,j,k
    real :: ztot

    res(:,:) = 0.
    IF (PRESENT(cond) .AND. PRESENT(dens)) THEN
        ! Conditional vertical integral with density weights
        DO j=3,n3-2
            DO i=3,n2-2
                DO k = 2,n1
                    IF (cond(k,i,j)) THEN
                        res(i,j) = res(i,j) + a(k,i,j)*dens(k,i,j)*(zm(k)-zm(k-1))
                    END IF
                END DO
           END DO
       END DO
    ELSEIF (PRESENT(cond)) THEN
        ! Conditional average
        DO j=3,n3-2
            DO i=3,n2-2
                ztot = 0.
                DO k = 2,n1
                    IF (cond(k,i,j)) THEN
                        res(i,j) = res(i,j) + a(k,i,j)*(zm(k)-zm(k-1))
                        ztot = ztot + (zm(k)-zm(k-1))
                    END IF
                END DO
                ! Grid weighted vertical average for columns with at least one available value
                if (ztot>1e-10) THEN
                    res(i,j) = res(i,j)/ztot
                ELSE
                    res(i,j) = -999.
                END IF
          END DO
       END DO
    ELSEIF (PRESENT(dens)) THEN
       ! Vertical integral with density weights
       DO j=3,n3-2
          DO i=3,n2-2
             DO k = 2,n1
                res(i,j) = res(i,j) + a(k,i,j)*dens(k,i,j)*(zm(k)-zm(k-1))
             END DO
          END DO
       END DO
    ELSE
       ! Average
       DO j=3,n3-2
          DO i=3,n2-2
             ztot = 0.
             DO k = 2,n1
                res(i,j) = res(i,j) + a(k,i,j)*(zm(k)-zm(k-1))
             END DO
             ! Grid weighted vertical average
             res(i,j) = res(i,j)/(zm(n1)-zm(1))
          END DO
       END DO
    END IF

  end SUBROUTINE get_avg_cs
  !
  ! -----------------------------------------------------------------------
  ! subroutine set_ts: computes and writes time sequence stats
  !
  subroutine set_ts(n1,n2,n3,w,th,t,dn0,zt,zm,dzt,th00,time)

    USE grid, ONLY : th0

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: w(n1,n2,n3),th(n1,n2,n3),t(n1,n2,n3)
    real, intent(in)    :: dn0(n1),zt(n1),zm(n1),dzt(n1),th00,time

    integer :: i, j, k
    real    :: bf(n1), scr(n2,n3)

    ssclr(1) = time
    ssclr(4) = get_zi_dmax(n1, n2, n3, th, zt) ! height of the maximum theta gradient
    ssclr(5) = zt( MAXLOC(thvar,DIM=1) )       ! height of the maximum theta variance
    !
    ! buoyancy flux statistics
    !
    ssclr(7) = 0.
    ssclr(32) = 0.
    do k = 2,n1-2
       ssclr(7) = ssclr(7) + (tke_res(k)+tke_sgs(k))*dn0(k)/dzt(k)
       ssclr(32) = ssclr(32) + (tke_res(k)+tke_sgs(k))/dzt(k)
       svctr(k,33) = svctr(k,33) + wtv_sgs(k)*9.8/th00
    end do
    bf(:) = wtv_res(:) + wtv_sgs(:)
    ssclr(6) = zm( MINLOC(bf,DIM=1) ) ! height of the minimum buoyancy flux

    ssclr(8) = bf(2)
    ssclr(9) = get_max_val(n1, n2, n3, w)

    scr(:,:) = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1
             scr(i,j)=scr(i,j)+(t(k,i,j)+th00-th0(k))*(zm(k)-zm(k-1))
          enddo
       end do
    end do
    ssclr(33) = get_avg2dh(n2,n3,scr)

  end subroutine set_ts
  !
  ! -----------------------------------------------------------------------
  ! subroutine ts_cld: computes and writes time sequence stats
  !
  subroutine ts_cld(n1,n2,n3,w,dn,zm,zt,dzt,rc,rt,rv,nc,rr,nr,rrate,crate)
    use grid, only : meanRadius
    integer, intent(in) :: n1,n2,n3
    real, intent(in), DIMENSION(n1) :: zm, zt, dzt
    real, intent(in), DIMENSION(n1,n2,n3) :: w, dn, rc, rt, rv, nc, rr, nr, rrate, crate

    integer :: k,i,j,n,m,l
    real    :: scr(n2,n3), scr1(n2,n3), scr2(n2,n3), ct_sum, cb_sum, ct_max, cb_min, rrcb, rccb, wvar, a(n1,n2,n3)
    INTEGER :: ct_tmp, cb_tmp

    ssclr(14) = get_zi_dmax(n1, n2, n3, rt, zt) ! height of the maximum total water gradient

    scr(:,:) = 0.   ! LWP
    scr1(:,:) = 0.  ! RWP
    scr2(:,:) = 0.  ! WVP (water vapor path)
    n = 0           ! Number of cloudy colums
    m = 0           ! Number of cloudy grid cells
    l = 0           ! Number of rainly grid cells
    ct_sum = 0.     ! Sum of cloud top and base height for mean heigths
    cb_sum = 0.
    ct_max=  zt(1)  ! Maximum cloud top and minimum cloud base heights
    cb_min = zt(n1)
    rrcb = 0.       ! Rain and cloud water precipitation rates at the cloud base
    rccb = 0.
    wvar = 0.       ! Cloud base vertical velocity variance
    do j=3,n3-2
       do i=3,n2-2
          ct_tmp = 1
          cb_tmp = n1
          do k=2,n1
             ! Cloudy grid cell
             if (cloudmask(k,i,j)) then
                m = m+1
                ct_tmp = max(ct_tmp,k)
                cb_tmp = min(cb_tmp,k)
             end if
             ! Rainy grid cell
             if (rainmask(k,i,j)) then
                l = l +1
             end if
             ! LWP, RWP and WVP
             scr(i,j)=scr(i,j)+rc(k,i,j)*dn(k,i,j)*(zm(k)-zm(k-1))
             scr1(i,j)=scr1(i,j)+rr(k,i,j)*dn(k,i,j)*(zm(k)-zm(k-1))
             scr2(i,j)=scr2(i,j)+rv(k,i,j)*dn(k,i,j)*(zm(k)-zm(k-1))
          end do
          IF (ct_tmp>1) THEN
             ! Cloudy column
             n = n+1
             ct_sum = ct_sum + zt(ct_tmp)
             cb_sum = cb_sum + zt(cb_tmp)
             ct_max = max(ct_max,zt(ct_tmp))
             cb_min = min(cb_min,zt(cb_tmp))
             ! Precipitation from the lowest cloudy grid cell (no cloud = zero precipitation)
             rrcb = rrcb + rrate(cb_tmp,i,j)
             rccb = rccb + crate(cb_tmp,i,j)
             ! Cloud base vertical velocity variance (zero mean assumed)
             wvar = wvar + w(cb_tmp,i,j)**2
          ENDIF
       end do
    end do

    ! Fraction of cloudy columns
    ssclr(19) = get_pustat_scalar('avg',REAL(n)/REAL((n3-4)*(n2-4)))

    ! Maximum cloud water mixing ratio
    ssclr(20)  = get_max_val(n1,n2,n3,rc)

    ! Total number of cloudy grid cells
    ssclr(22) = get_pustat_scalar('sum',REAL(m))

    IF (ssclr(22)>0.) THEN
        ssclr(17) = get_pustat_scalar('max',ct_max) ! Maximum cloud top height
        ssclr(18) = get_pustat_scalar('min',cb_min) ! Minimum cloud base height
        ssclr(23) = get_pustat_scalar('avg',ct_sum,REAL(n)) ! Mean cloud top height
        ssclr(24) = get_pustat_scalar('avg',cb_sum,REAL(n)) ! Mean cloud base height
        ssclr(35) = get_pustat_scalar('avg',wvar,REAL(n))
    ELSE
        ssclr((/17,18,23,24,25,35/)) = -999.
    ENDIF

    ! liquid water path (without precipitation)
    ssclr(15) = get_avg2dh(n2,n3,scr)
    scr(:,:)=(scr(:,:)-ssclr(15))**2 ! For LWP variance
    ssclr(16) = get_avg2dh(n2,n3,scr)

    ! water vapor path
    ssclr(26) = get_avg2dh(n2,n3,scr2)

    ! rain water path
    ssclr(27) = get_avg2dh(n2,n3,scr1)

    ! surface precipitation
    scr(:,:) = rrate(2,:,:) + crate(2,:,:)
    ssclr(28) = get_avg2dh(n2,n3,scr)

    ! total number of rain grid cells
    ssclr(30) = get_pustat_scalar('sum',REAL(l))
    ! average cloud base precipitation rate
    ssclr(31) = get_pustat_scalar('avg',rrcb/REAL((n3-4)*(n2-4)))
    ssclr(34) = get_pustat_scalar('avg',rccb/REAL((n3-4)*(n2-4)))

    ! average cloud and rain drop number concentrations
    ssclr(25) = get_avg_ts(n1,n2,n3,nc,dzt,cloudmask)
    ssclr(29) = get_avg_ts(n1,n2,n3,nr,dzt,rainmask)
    ! average radius
    IF (level<4) THEN
        CALL getSBradius(n1,n2,n3,nc,rc,-1,a)
        ssclr(36) = get_avg_ts(n1,n2,n3,a,dzt,cloudmask)
        CALL getSBradius(n1,n2,n3,nr,rr,0,a)
        ssclr(37) = get_avg_ts(n1,n2,n3,a,dzt,rainmask)
    ELSE
        CALL meanRadius('cloud','ab',a)
        ssclr(36) = get_avg_ts(n1,n2,n3,a,dzt,cloudmask)
        CALL meanRadius('precp','ab',a)
        ssclr(37) = get_avg_ts(n1,n2,n3,a,dzt,rainmask)
    ENDIF

  end subroutine ts_cld
  !
  ! -----------------------------------------------------------------------
  ! subroutine ts_ice: computes and writes time sequence stats
  !
  subroutine ts_ice(n1,n2,n3)
    USE grid, ONLY : a_rip, a_nip, a_rsp, a_nsp, a_rgp, a_ngp, a_rhp, a_nhp, &
        a_rsi, a_rv, dzt, a_dn, zm, icein, snowin, grin, hailin, a_tp, th0, th00
    integer, intent (in) :: n1,n2,n3

    integer :: i, j, k
    real    :: scr(n2,n3), a1(n1,n2,n3)
    LOGICAL :: mask(n1,n2,n3)

    ! Ice water path
    ssclr_ice(1) = get_avg_ts(n1,n2,n3,a_rip,dzt,dens=a_dn)
    ! Mask for number concentration and radius
    mask = (a_nip > ni_min .OR. a_rip > ri_min)
    ! Radius
    CALL getSBradius(n1,n2,n3,a_nip,a_rip,1,a1)
    ssclr_ice(2) = get_avg_ts(n1,n2,n3,a1,dzt,mask)
    ! Number concentration
    ssclr_ice(3) = get_avg_ts(n1,n2,n3,a_nip,dzt,mask)
    ! Total number of icy grid cells
    k=COUNT( mask(2:n1,3:n2-2,3:n3-2) )
    ssclr_ice(4) = get_pustat_scalar('sum',REAL(k))
    ! Surface precipitation rate
    scr(:,:) = icein(2,:,:)
    ssclr_ice(5) = get_avg2dh(n2,n3,scr)

    ! The same for snow
    ssclr_ice(6) = get_avg_ts(n1,n2,n3,a_rsp,dzt,dens=a_dn)
    mask = (a_nsp > ni_min .OR. a_rsp > ri_min)
    CALL getSBradius(n1,n2,n3,a_nsp,a_rsp,2,a1)
    ssclr_ice(7) = get_avg_ts(n1,n2,n3,a1,dzt,mask)
    ssclr_ice(8) = get_avg_ts(n1,n2,n3,a_nsp,dzt,mask)
    k=COUNT( mask(2:n1,3:n2-2,3:n3-2) )
    ssclr_ice(9) = get_pustat_scalar('sum',REAL(k))
    scr(:,:) = snowin(2,:,:)
    ssclr_ice(10) = get_avg2dh(n2,n3,scr)

    ! The same for graupel
    ssclr_ice(11) = get_avg_ts(n1,n2,n3,a_rgp,dzt,dens=a_dn)
    mask = (a_ngp > ni_min .OR. a_rgp > ri_min)
    CALL getSBradius(n1,n2,n3,a_ngp,a_rgp,3,a1)
    ssclr_ice(12) = get_avg_ts(n1,n2,n3,a1,dzt,mask)
    ssclr_ice(13) = get_avg_ts(n1,n2,n3,a_ngp,dzt,mask)
    k=COUNT( mask(2:n1,3:n2-2,3:n3-2) )
    ssclr_ice(14) = get_pustat_scalar('sum',REAL(k))
    scr(:,:) = grin(2,:,:)
    ssclr_ice(15) = get_avg2dh(n2,n3,scr)

    ! The same for hail
    ssclr_ice(16) = get_avg_ts(n1,n2,n3,a_rhp,dzt,dens=a_dn)
    mask = (a_nhp > ni_min .OR. a_rhp > ri_min)
    CALL getSBradius(n1,n2,n3,a_nsp,a_rsp,4,a1)
    ssclr_ice(17) = get_avg_ts(n1,n2,n3,a1,dzt,mask)
    ssclr_ice(18) = get_avg_ts(n1,n2,n3,a_nhp,dzt,mask)
    k=COUNT( mask(2:n1,3:n2-2,3:n3-2) )
    ssclr_ice(19) = get_pustat_scalar('sum',REAL(k))
    scr(:,:) = hailin(2,:,:)
    ssclr_ice(20) = get_avg2dh(n2,n3,scr)

    ! Maximum supersaturation over ice
    a1=0.
    WHERE(a_rsi>1e-10) a1=a_rv/a_rsi
    ssclr_ice(21) = (get_max_val(n1,n2,n3,a1)-1.0)*100.

    ! Integrated ice-liquid water potential temperature - change from th0
    scr = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1
             scr(i,j)=scr(i,j)+(a_tp(k,i,j)+th00-th0(k))*(zm(k)-zm(k-1))
          end do
       end do
    end do
    ssclr_ice(22) = get_avg2dh(n2,n3,scr)

  end subroutine ts_ice
  !
  ! -----------------------------------------------------------------------
  ! subroutine ts_lvl4: computes and writes time sequence stats of Salsa variables --
  !  Implemented by Zubair Maalick 20/07/2015
  !  Some rewriting and adjusting by Juha Tonttila
  !
  SUBROUTINE ts_lvl4(n1,n2,n3)
    USE grid, ONLY : a_rp, a_rsl
    IMPLICIT NONE
    integer, intent(in) :: n1,n2,n3
    REAL :: a0(n1,n2,n3)

    ! Maximum supersaturation
    a0=0.
    WHERE(a_rsl>1e-10) a0=a_rp/a_rsl
    ssclr_lvl4(1) = (get_max_val(n1,n2,n3,a0)-1.0)*100.

  END SUBROUTINE ts_lvl4
  !
  ! -----------------------------------------------------------------------
  ! subroutine ts_lvl5: computes and writes time sequence stats of Salsa variables --
  !  Implemented by Jaakko Ahola 15/12/2016
  !
  SUBROUTINE ts_lvl5(n1,n2,n3)
    USE grid, ONLY : bulkNumc, meanRadius, dzt, a_dn, zm, a_ri, a_srs, icein, snowin, a_rp, a_rsi, a_tp, th0, th00
    IMPLICIT NONE
    integer, intent(in) :: n1,n2,n3
    REAL :: a0(n1,n2,n3), scr(n2,n3)
    integer :: i, j, k

    ! IWP and SWP
    ssclr_lvl5(1) = get_avg_ts(n1,n2,n3,a_ri,dzt,dens=a_dn)
    ssclr_lvl5(6) = get_avg_ts(n1,n2,n3,a_srs,dzt,dens=a_dn)

    ! Ice and snow mean radius
    CALL meanRadius('ice','ab',a0)
    ssclr_lvl5(2) = get_avg_ts(n1,n2,n3,a0,dzt,icemask)
    CALL meanRadius('snow','ab',a0)
    ssclr_lvl5(7) = get_avg_ts(n1,n2,n3,a0,dzt,snowmask)

    ! Ice and snow number concentrations
    CALL bulkNumc('ice','ab',a0)
    ssclr_lvl5(3) = get_avg_ts(n1,n2,n3,a0,dzt,icemask)
    CALL bulkNumc('snow','a',a0)
    ssclr_lvl5(8) = get_avg_ts(n1,n2,n3,a0,dzt,snowmask)

    ! Total number of icy and snowy grid cells
    k=COUNT( icemask(2:n1,3:n2-2,3:n3-2) )
    ssclr_lvl5(4) = get_pustat_scalar('sum',REAL(k))
    k = COUNT( snowmask(2:n1,3:n2-2,3:n3-2) )
    ssclr_lvl5(9) = get_pustat_scalar('sum',REAL(k))

    ! Surface precipitation rates
    scr(:,:) = icein(2,:,:)
    ssclr_lvl5(5) = get_avg2dh(n2,n3,scr)
    scr(:,:) = snowin(2,:,:)
    ssclr_lvl5(10) = get_avg2dh(n2,n3,scr)

    ! Maximum supersaturation over ice
    a0=0.
    WHERE(a_rsi>1e-10) a0=a_rp/a_rsi
    ssclr_lvl5(11) = (get_max_val(n1,n2,n3,a0)-1.0)*100.

    ! Integrated ice-liquid water potential temperature - change from th0
    scr = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1
             scr(i,j)=scr(i,j)+(a_tp(k,i,j)+th00-th0(k))*(zm(k)-zm(k-1))
          end do
       end do
    end do
    ssclr_lvl5(12) = get_avg2dh(n2,n3,scr)

    ! Removal statistics elsewhere ..

  END SUBROUTINE ts_lvl5

  !---------------------------------------------------------------------
  ! SUBROUTINE ts_user_stats: computes user-defined time series outputs.
  ! Variable names are given in NAMELIST/user_ts_list and these must be defined in ncio.f90.
  ! Outputs are calculated here to array user_ts_data(nv1_user).
  subroutine ts_user_stats()
    use grid, ONLY : CCN, nzp, nxp, nyp, dzt, a_dn, a_temp
    INTEGER :: i
    REAL :: a(nzp,nxp,nyp), a1
    LOGICAL :: fail, mask(nzp,nxp,nyp), mass
    !
    DO i=1,nv1_user
        SELECT CASE (user_ts_list(i))
        CASE ('CCN')
            ! Level 3 CCN as an example of output
            IF (level<4) THEN
                user_ts_data(i) = CCN
            ELSE
                user_ts_data(i) = -999.
            ENDIF
        CASE ('T_min')
            ! Minimum absolute temperature (K)
            a1=MINVAL(a_temp(:,3:nxp-2,3:nyp-2))
            user_ts_data(i) = get_pustat_scalar('min',a1)
         CASE DEFAULT
            ! Pre-defined SALSA outputs
            fail = calc_user_data(user_ts_list(i),a,mask,is_mass=mass)
            IF (fail) THEN
                WRITE(*,*)" Error: failed to calculate '"//TRIM(user_ts_list(i))//"' for ts output!"
                STOP
            ENDIF
            ! Calculate vertical integral when mass concentration, otherwise mean
            IF (mass) THEN
                user_ts_data(i) = get_avg_ts(nzp,nxp,nyp,a,dzt,cond=mask,dens=a_dn)
            ELSE
                user_ts_data(i) = get_avg_ts(nzp,nxp,nyp,a,dzt,cond=mask)
            ENDIF
        END SELECT
    ENDDO
    !
  end subroutine ts_user_stats

  !
  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_STAT: Accumulates various statistics over an
  ! averaging period for base (level 0) version of model
  !
  subroutine accum_stat(n1,n2,n3,u,v,w,t,p,um,vm,th00)

    integer, intent (in) :: n1,n2,n3
    real, dimension (n1,n2,n3), intent (in)    :: u, v, w, t, p
    real, intent (in)           :: um, vm, th00

    integer :: k
    real    :: a1(n1), b1(n1), c1(n1), d1(n1), a3(n1), b3(n1), tmp(n1)

    call get_avg3(n1,n2,n3, u,a1)
    call get_avg3(n1,n2,n3, v,b1)
    call get_avg3(n1,n2,n3, t,c1)
    call get_avg3(n1,n2,n3, p,d1)
    call get_var3(n1,n2,n3, t, c1, thvar)
    call get_3rd3(n1,n2,n3, t, c1, b3) ! Used to be (t-a1)**3
    tmp(:)=0.
    call get_3rd3(n1,n2,n3, w, tmp, a3) ! Now just w**3

    do k=1,n1
       svctr(k,10)=svctr(k,10) + a1(k) + um
       svctr(k,11)=svctr(k,11) + b1(k) + vm
       svctr(k,12)=svctr(k,12) + c1(k) + th00
       svctr(k,13)=svctr(k,13) + d1(k)
       svctr(k,17)=svctr(k,17) + thvar(k)
       svctr(k,18)=svctr(k,18) + a3(k)
       svctr(k,19)=svctr(k,19) + b3(k)
    end do

  end subroutine accum_stat
  !
  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_RAD: Accumulates various statistics over an
  ! averaging period for radiation variables
  !
  subroutine accum_rad(n1,n2,n3,rflx,sflx,sup,sdwn,irup,irdwn,alb)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: rflx(n1,n2,n3)
    real, optional, intent (in) :: sflx(n1,n2,n3), alb(n2,n3)
    real, optional, intent (in) :: sup(n1,n2,n3), sdwn(n1,n2,n3), irup(n1,n2,n3), irdwn(n1,n2,n3)

    integer :: k
    real    :: a1(n1),a2(n1)

    call get_avg3(n1,n2,n3,rflx,a1)
    call get_var3(n1,n2,n3,rflx,a1,a2)
    do k=1,n1
       svctr(k,55)=svctr(k,55) + a1(k)
       svctr(k,56)=svctr(k,56) + a2(k)
    end do

    if (present(sflx)) then
       call get_avg3(n1,n2,n3,sflx,a1)
       call get_var3(n1,n2,n3,sflx,a1,a2)
       do k=1,n1
          svctr(k,57)=svctr(k,57) + a1(k)
          svctr(k,58)=svctr(k,58) + a2(k)
       end do
       ssclr(21) = get_avg2dh(n2,n3,alb)
    end if

    if (present(sup)) then
        call get_avg3(n1,n2,n3,sup,a1)
        svctr(:,93)=svctr(:,93) + a1(:)
        call get_avg3(n1,n2,n3,sdwn,a1)
        svctr(:,94)=svctr(:,94) + a1(:)
        call get_avg3(n1,n2,n3,irup,a1)
        svctr(:,95)=svctr(:,95) + a1(:)
        call get_avg3(n1,n2,n3,irdwn,a1)
        svctr(:,96)=svctr(:,96) + a1(:)
    end if

  end subroutine accum_rad
  !
  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_CLD: Accumulates cloud statistics.
  !
  subroutine accum_cld(n1, n2, n3, th00, w, th, t, rl, rt, rv, crate, nc, rr, nr, rrate)
    use grid, only : meanRadius
    use defs, only : ep2

    integer, intent (in) :: n1,n2,n3
    real, intent (in)                      :: th00
    real, intent (in), dimension(n1,n2,n3) :: w, th, t, rl, rt, rv, crate, nc, rr, nr, rrate

    real, dimension(n1,n2,n3) :: tv    ! Local variable
    integer                   :: k, i, j, kp1
    real, dimension(n1)       :: a1, a2, a3, tvbar
    real, dimension(n1,n2,n3) :: xy1, xy2, tw, tvw, rtw
    LOGICAL :: cond(n1,n2,n3)

    ! total water statistics
    call get_avg3(n1,n2,n3,rt,a1)
    call get_var3(n1,n2,n3,rt,a1,a2)
    CALL get_3rd3(n1,n2,n3,rt,a1,a3)
    svctr(:,50)=svctr(:,50) + a1(:)
    svctr(:,51)=svctr(:,51) + a2(:)
    svctr(:,52)=svctr(:,52) + a3(:)

    ! liquid water statistics
    call get_avg3(n1,n2,n3,rl,a1)
    call get_var3(n1,n2,n3,rl,a1,a2)
    call get_3rd3(n1,n2,n3,rl,a1,a3)
    svctr(:,59)=svctr(:,59) + a1(:)
    svctr(:,60)=svctr(:,60) + a2(:)
    svctr(:,61)=svctr(:,61) + a3(:)

    ! water vapor mixing ratio
    call get_avg3(n1,n2,n3,rv,a1)
    svctr(:,84)=svctr(:,84) + a1(:)

    ! cloud water deposition flux
    call get_avg3(n1,n2,n3,crate,a1)
    svctr(:,85)=svctr(:,85) + a1(:)

    ! rain water mixing ratio
    call get_avg3(n1,n2,n3,rr,a1)
    svctr(:,89)=svctr(:,89) + a1(:) ! rr (kg/kg)

    ! rain water precipitation flux
    call get_avg3(n1,n2,n3,rrate,a1)
    svctr(:,90)=svctr(:,90)+a1(:)

    ! cloud and rain statistics
    ! Droplet number concentration
    CALL get_avg3(n1,n2,n3,nc,a1,cond=cloudmask)
    svctr(:,87) = svctr(:,87) + a1(:)
    CALL get_avg3(n1,n2,n3,nr,a1,cond=rainmask)
    svctr(:,92) = svctr(:,92) + a1(:)
    ! Average radius
    IF (level<4) THEN
        CALL getSBradius(n1,n2,n3,nc,rl,-1,xy1)
        CALL getSBradius(n1,n2,n3,nr,rr,0,xy2)
    ELSE
        CALL meanRadius('cloud','ab',xy1)
        CALL meanRadius('precp','ab',xy2)
    ENDIF
    CALL get_avg3(n1,n2,n3,xy1,a1,cond=cloudmask)
    svctr(:,97) = svctr(:,97) + a1(:)
    CALL get_avg3(n1,n2,n3,xy2,a1,cond=rainmask)
    svctr(:,98) = svctr(:,98) + a1(:)
    ! Fraction of cloudy/rainy grid cells
    WHERE(cloudmask)
        xy1=1.
    ELSEWHERE
        xy1=0.
    ENDWHERE
    CALL get_avg3(n1,n2,n3,xy1,a1)
    svctr(:,86) = svctr(:,86) + a1(:)
    WHERE(rainmask)
        xy1=1.
    ELSEWHERE
        xy1=0.
    ENDWHERE
    CALL get_avg3(n1,n2,n3,xy1,a1)
    svctr(:,91) = svctr(:,91) + a1(:)
    !
    ! do some conditional sampling statistics: cloud, cloud-core
    !
    tv(:,:,:) = th(:,:,:)*(1.+ep2*rt(:,:,:) - rl(:,:,:))
    call get_avg3(n1,n2,n3,tv,tvbar)
    !
    xy1=0.
    xy2=0.
    tvw=0.
    do k=1,n1
       kp1=k+1
       if (k==n1) kp1=k
       do j=3,n3-2
          do i=3,n2-2
             if (rl(k,i,j) > 1.e-5) then
                xy1(k,i,j)=1.
                if (tv(k,i,j) > tvbar(k)) THEN
                   xy2(k,i,j)=1.
                end if
                !
                tw(k,i,j)=(.5*(t(k,i,j)+t(kp1,i,j))+th00)*w(k,i,j)
                tvw(k,i,j)=(.5*(tv(k,i,j)+tv(kp1,i,j)))*w(k,i,j)
                rtw(k,i,j)=(.5*(rt(k,i,j)+rt(kp1,i,j)))*w(k,i,j)
             end if
          end do
       end do
    end do
    CALL get_avg3(n1,n2,n3,xy1,a1)
    svctr(:,64)=svctr(:,64)+a1(:)
    CALL get_avg3(n1,n2,n3,xy1,a1,normalize=.FALSE.)
    svctr(:,65)=svctr(:,65)+a1(:)   ! Counts
    cond(:,:,:)=xy1(:,:,:)>0.5
    CALL get_avg3(n1,n2,n3,w,a1,cond=cond)
    svctr(:,66)=svctr(:,66)+a1(:)
    CALL get_avg3(n1,n2,n3,t+th00,a1,cond=cond)
    svctr(:,67)=svctr(:,67)+a1(:)
    CALL get_avg3(n1,n2,n3,tv,a1,cond=cond)
    svctr(:,68)=svctr(:,68)+a1(:)
    CALL get_avg3(n1,n2,n3,rt,a1,cond=cond)
    svctr(:,69)=svctr(:,69)+a1(:)
    CALL get_avg3(n1,n2,n3,rl,a1,cond=cond)
    svctr(:,70)=svctr(:,70)+a1(:)
    CALL get_avg3(n1,n2,n3,tw,a1,cond=cond)
    svctr(:,71)=svctr(:,71)+a1(:)
    CALL get_avg3(n1,n2,n3,tvw,a1,cond=cond)
    svctr(:,72)=svctr(:,72)+a1(:)
    CALL get_avg3(n1,n2,n3,rtw,a1,cond=cond)
    svctr(:,73)=svctr(:,73)+a1(:)

    CALL get_avg3(n1,n2,n3,xy2,a1)
    svctr(:,74)=svctr(:,74)+a1(:)
    CALL get_avg3(n1,n2,n3,xy2,a1,normalize=.FALSE.)
    svctr(:,75)=svctr(:,75)+a1(:)   ! Counts
    cond(:,:,:)=xy2(:,:,:)>0.5
    CALL get_avg3(n1,n2,n3,w,a1,cond=cond)
    svctr(:,76)=svctr(:,76)+a1(:)
    CALL get_avg3(n1,n2,n3,t+th00,a1,cond=cond)
    svctr(:,77)=svctr(:,77)+a1(:)
    CALL get_avg3(n1,n2,n3,tv,a1,cond=cond)
    svctr(:,78)=svctr(:,78)+a1(:)
    CALL get_avg3(n1,n2,n3,rt,a1,cond=cond)
    svctr(:,79)=svctr(:,79)+a1(:)
    CALL get_avg3(n1,n2,n3,rl,a1,cond=cond)
    svctr(:,80)=svctr(:,80)+a1(:)
    CALL get_avg3(n1,n2,n3,tw,a1,cond=cond)
    svctr(:,81)=svctr(:,81)+a1(:)
    CALL get_avg3(n1,n2,n3,tvw,a1,cond=cond)
    svctr(:,82)=svctr(:,82)+a1(:)
    CALL get_avg3(n1,n2,n3,rtw,a1,cond=cond)
    svctr(:,83)=svctr(:,83)+a1(:)

  end subroutine accum_cld
  !
  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_ice: Accumulates ice statistics.
  !
  subroutine accum_ice(n1,n2,n3)
    use grid, ONLY : a_nip, a_rip, a_rsp, a_nsp, a_rgp, a_ngp, a_rhp, a_nhp, &
        a_tp, icein, snowin, grin, hailin, th00
    IMPLICIT NONE

    integer, intent (in) :: n1,n2,n3
    real :: a1(n1,n2,n3), col(n1)
    LOGICAL :: mask(n1,n2,n3)

    ! Ice mass
    CALL get_avg3(n1,n2,n3,a_rip,col)
    svctr_ice(:,1) = svctr_ice(:,1) + col(:)
    ! Ice number - requires mask
    mask = (a_nip > ni_min .OR. a_rip > ri_min)
    CALL get_avg3(n1,n2,n3,a_nip,col,cond=mask)
    svctr_ice(:,2) = svctr_ice(:,2) + col(:)
    ! Ice radius
    CALL getSBradius(n1,n2,n3,a_nip,a_rip,1,a1)
    CALL get_avg3(n1,n2,n3,a1,col,cond=mask)
    svctr_ice(:,3) = svctr_ice(:,3) + col(:)
    ! Fraction of icy grid cells
    a1 = merge(1., 0., mask)
    CALL get_avg3(n1,n2,n3,a1,col)
    svctr_ice(:,4) = svctr_ice(:,4) + col(:)
    ! Ice deposition flux
    call get_avg3(n1,n2,n3,icein,col)
    svctr_ice(:,5) = svctr_ice(:,5) + col(:)

    ! The same for snow
    CALL get_avg3(n1,n2,n3,a_rsp,col)
    svctr_ice(:,6) = svctr_ice(:,6) + col(:)
    mask = (a_nsp > ni_min .OR. a_rsp > ri_min)
    CALL get_avg3(n1,n2,n3,a_nsp,col,cond=mask)
    svctr_ice(:,7) = svctr_ice(:,7) + col(:)
    CALL getSBradius(n1,n2,n3,a_nsp,a_rsp,2,a1)
    CALL get_avg3(n1,n2,n3,a1,col,cond=mask)
    svctr_ice(:,8) = svctr_ice(:,8) + col(:)
    a1 = merge(1., 0., mask)
    CALL get_avg3(n1,n2,n3,a1,col)
    svctr_ice(:,9) = svctr_ice(:,9) + col(:)
    call get_avg3(n1,n2,n3,snowin,col)
    svctr_ice(:,10) = svctr_ice(:,10) + col(:)

    ! The same for graupel
    CALL get_avg3(n1,n2,n3,a_rgp,col)
    svctr_ice(:,11) = svctr_ice(:,11) + col(:)
    mask = (a_ngp > ni_min .OR. a_rgp > ri_min)
    CALL get_avg3(n1,n2,n3,a_nhp,col,cond=mask)
    svctr_ice(:,12) = svctr_ice(:,12) + col(:)
    CALL getSBradius(n1,n2,n3,a_ngp,a_rgp,3,a1)
    CALL get_avg3(n1,n2,n3,a1,col,cond=mask)
    svctr_ice(:,13) = svctr_ice(:,13) + col(:)
    a1 = merge(1., 0., mask)
    CALL get_avg3(n1,n2,n3,a1,col)
    svctr_ice(:,14) = svctr_ice(:,14) + col(:)
    call get_avg3(n1,n2,n3,grin,col)
    svctr_ice(:,15) = svctr_ice(:,15) + col(:)

    ! The same for hail
    CALL get_avg3(n1,n2,n3,a_rhp,col)
    svctr_ice(:,16) = svctr_ice(:,16) + col(:)
    mask = (a_nhp > ni_min .OR. a_rhp > ri_min)
    CALL get_avg3(n1,n2,n3,a_nhp,col,cond=mask)
    svctr_ice(:,17) = svctr_ice(:,17) + col(:)
    CALL getSBradius(n1,n2,n3,a_nhp,a_rhp,4,a1)
    CALL get_avg3(n1,n2,n3,a1,col,cond=mask)
    svctr_ice(:,18) = svctr_ice(:,18) + col(:)
    a1 = merge(1., 0., mask)
    CALL get_avg3(n1,n2,n3,a1,col)
    svctr_ice(:,19) = svctr_ice(:,19) + col(:)
    call get_avg3(n1,n2,n3,hailin,col)
    svctr_ice(:,20) = svctr_ice(:,20) + col(:)

    ! Ice-liquid water potential temperature
    call get_avg3(n1,n2,n3,a_tp,col)
    svctr_ice(:,21) = svctr_ice(:,21) + (col(:) + th00)

  end subroutine accum_ice

  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_LVL4: Accumulates specialized statistics that depend
  ! on level 4 variables.
  !
  subroutine accum_lvl4(n1,n2,n3)
    use mo_submctl, only : inp2a,fnp2a,inp2b,fnp2b,cldbinlim,nout_cld,nlim
    use grid, ONLY : getBinRadius, nspec, ncld, a_mcloudp, a_ncloudp

    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3

    REAL, DIMENSION(n1,n2,n3,fnp2b)     :: a_Rwet
    REAL, ALLOCATABLE                   :: hist(:,:)

    ! Cloud droplet histograms
    ! ------------------------
    IF (nout_cld>0) THEN
        ALLOCATE(hist(n1,nout_cld))
        ! Cloud droplet bin wet radius
        CALL getBinRadius(ncld,nspec+1,a_ncloudp,a_mcloudp,nlim,a_Rwet,2)
        ! Histograms (regime A)
        CALL HistDistr(n1,n2,n3,fnp2a,a_Rwet(:,:,:,inp2a:fnp2a),a_ncloudp(:,:,:,inp2a:fnp2a),cldbinlim,nout_cld,hist)
        svctr_ch(:,:,1) = svctr_ch(:,:,1) + hist(:,:)
        ! Histograms (regime B)
        IF (ncld==fnp2b) THEN
            CALL HistDistr(n1,n2,n3,fnp2b-fnp2a,a_Rwet(:,:,:,inp2b:fnp2b),a_ncloudp(:,:,:,inp2b:fnp2b),cldbinlim,nout_cld,hist)
            svctr_ch(:,:,2) = svctr_ch(:,:,2) + hist(:,:)
        ENDIF
        DEALLOCATE(hist)
    ENDIF
  end subroutine accum_lvl4

  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_LVL5: Accumulates specialized statistics that depend
  ! on level 5 variables.
  !
  subroutine accum_lvl5(n1,n2,n3)
    use mo_submctl, only : inp2a,fnp2a,inp2b,fnp2b,icebinlim,nout_ice, &
                     prlim ! Note: prlim in #/m^3, but close enough to #/kg for statistics
    use grid, ONLY : bulkNumc, bulkMixrat, meanRadius, binSpecMixrat, getBinRadius, &
                     a_ri, a_srs, nspec, nice, a_micep, a_nicep, icein, snowin, a_tp, th00

    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3

    REAL, DIMENSION(n1,n2,n3)           :: a1
    REAL, DIMENSION(n1)                 :: col
    REAL, DIMENSION(n1,n2,n3,fnp2b)     :: a_Rwet
    REAL, ALLOCATABLE                   :: hist(:,:)

    ! Ice mass
    CALL get_avg3(n1,n2,n3,a_ri,col)
    svctr_lvl5(:,1) = svctr_lvl5(:,1) + col(:)

    ! Generate SALSA ice mask
    CALL bulkNumc('ice','ab',a1)
    icemask(:,:,:) = ( a1(:,:,:) > prlim .AND. a_ri(:,:,:) > ris_min)

    ! Ice number concentration
    CALL get_avg3(n1,n2,n3,a1,col,cond=icemask)
    svctr_lvl5(:,2) = svctr_lvl5(:,2) + col(:)
    ! Mean radius
    CALL meanRadius('ice','ab',a1)
    CALL get_avg3(n1,n2,n3,a1,col,cond=icemask)
    svctr_lvl5(:,3) = svctr_lvl5(:,3) + col(:)
    ! Fraction of icy grid cells
    WHERE(icemask)
        a1=1.
    ELSEWHERE
        a1=0.
    ENDWHERE
    CALL get_avg3(n1,n2,n3,a1,col)
    svctr_lvl5(:,4) = svctr_lvl5(:,4) + col(:)

    ! Ice deposition flux
    call get_avg3(n1,n2,n3,icein,col)
    svctr_lvl5(:,5) = svctr_lvl5(:,5) + col(:)

    ! The same for snow (mass, mask, number concentration and fraction of snowy grid cells)
    CALL get_avg3(n1,n2,n3,a_srs,col)
    svctr_lvl5(:,6) = svctr_lvl5(:,6) + col(:)

    CALL bulkNumc('snow','a',a1)
    snowmask(:,:,:) = ( a1(:,:,:) > prlim .AND. a_srs(:,:,:) > ris_min)

    CALL get_avg3(n1,n2,n3,a1,col,cond=snowmask)
    svctr_lvl5(:,7) = svctr_lvl5(:,7) + col(:)

    CALL meanRadius('snow','ab',a1)
    CALL get_avg3(n1,n2,n3,a1,col,cond=snowmask)
    svctr_lvl5(:,8) = svctr_lvl5(:,8) + col(:)

    WHERE(snowmask)
        a1=1.
    ELSEWHERE
        a1=0.
    ENDWHERE
    CALL get_avg3(n1,n2,n3,a1,col)
    svctr_lvl5(:,9) = svctr_lvl5(:,9) + col(:)

    call get_avg3(n1,n2,n3,snowin,col)
    svctr_lvl5(:,10) = svctr_lvl5(:,10) + col(:)

    ! Ice-liquid water potential temperature
    call get_avg3(n1,n2,n3,a_tp,col)
    svctr_lvl5(:,11) = svctr_lvl5(:,11) + (col(:) + th00)

    ! Ice histograms
    ! --------------
    IF (nout_ice>0) THEN
        ALLOCATE(hist(n1,nout_ice))
        ! Cloud droplet bin wet radius
        CALL getBinRadius(nice,nspec+1,a_nicep,a_micep,prlim,a_Rwet,4)
        ! Histograms (regime A)
        CALL HistDistr(n1,n2,n3,fnp2a,a_Rwet(:,:,:,inp2a:fnp2a),a_nicep(:,:,:,inp2a:fnp2a),icebinlim,nout_ice,hist)
        svctr_ih(:,:,1) = svctr_ih(:,:,1) + hist(:,:)
        ! Histograms (regime B)
        IF (nice==fnp2b) THEN
            CALL HistDistr(n1,n2,n3,fnp2a,a_Rwet(:,:,:,inp2b:fnp2b),a_nicep(:,:,:,inp2b:fnp2b),icebinlim,nout_ice,hist)
            svctr_ih(:,:,2) = svctr_ih(:,:,2) + hist(:,:)
        ENDIF
        DEALLOCATE(hist)
    ENDIF
  end subroutine accum_lvl5

  !---------------------------------------------------------------------
  ! SUBROUTINE ps_user_stats: computes user-defined profile outputs.
  ! Variable names are given in NAMELIST/user_ps_list and these must be defined in ncio.f90.
  ! Outputs are calculated here to array user_ps_data(nzp,nv2_user).
  subroutine ps_user_stats()
    USE grid, ONLY : nzp, nxp, nyp, a_dn, a_rv, a_rp, a_rsl, a_rsi, a_temp
    INTEGER :: i
    LOGICAL :: fail, mask(nzp,nxp,nyp)
    REAL :: a(nzp,nxp,nyp), a1(nzp)
    !
    DO i=1,nv2_user
        SELECT CASE (user_ps_list(i))
        CASE('rho_air')
            ! Air density
            call get_avg3(nzp,nxp,nyp,a_dn,a1)
            user_ps_data(:,i)=user_ps_data(:,i)+a1(:)
        CASE ('SS_max')
            ! Maximum supersaturation (%) over water
            a=0.
            IF (level<4) THEN
                WHERE(a_rsl>1e-10) a=a_rv/a_rsl
            ELSE
                WHERE(a_rsl>1e-10) a=a_rp/a_rsl
            ENDIF
            a1(:)=MAXVAL(MAXVAL(a,DIM=3),DIM=2)
            CALL get_pustat_vector('max', nzp, a)
            user_ps_data(:,i) = user_ps_data(:,i) + (a1(:)-1.0)*100.
        CASE ('SSi_max')
            ! Maximum supersaturation (%) over ice
            a=0.
            IF (level<4) THEN
                WHERE(a_rsi>1e-10) a=a_rv/a_rsi
            ELSE
                WHERE(a_rsi>1e-10) a=a_rp/a_rsi
            ENDIF
            a1(:)=MAXVAL(MAXVAL(a,DIM=3),DIM=2)
            CALL get_pustat_vector('max', nzp, a)
            user_ps_data(:,i) = user_ps_data(:,i) + (a1(:)-1.0)*100.
        CASE ('SS_min')
            ! Minimum supersaturation (%) over water
            a=999.
            IF (level<4) THEN
                WHERE(a_rsl>1e-10) a=a_rv/a_rsl
            ELSE
                WHERE(a_rsl>1e-10) a=a_rp/a_rsl
            ENDIF
            a1(:)=MINVAL(MINVAL(a,DIM=3),DIM=2)
            CALL get_pustat_vector('min', nzp, a)
            user_ps_data(:,i) = user_ps_data(:,i) + (a1(:)-1.0)*100.
        CASE ('SSi_min')
            ! Minimum supersaturation (%) over ice
            a=999.
            IF (level<4) THEN
                WHERE(a_rsi>1e-10) a=a_rv/a_rsi
            ELSE
                WHERE(a_rsi>1e-10) a=a_rp/a_rsi
            ENDIF
            a1(:)=MINVAL(MINVAL(a,DIM=3),DIM=2)
            CALL get_pustat_vector('min', nzp, a)
            user_ps_data(:,i) = user_ps_data(:,i) + (a1(:)-1.0)*100.
        CASE ('T_max')
            ! Maximum absolute temperature (K)
            a1(:)=MAXVAL(MAXVAL(a_temp(:,3:nxp-2,3:nyp-2),DIM=3),DIM=2)
            CALL get_pustat_vector('max', nzp, a)
            user_ps_data(:,i) = user_ps_data(:,i) + a1(:)
        CASE ('T_min')
            ! Minimum absolute temperature (K)
            a1(:)=MINVAL(MINVAL(a_temp(:,3:nxp-2,3:nyp-2),DIM=3),DIM=2)
            CALL get_pustat_vector('min', nzp, a)
            user_ps_data(:,i) = user_ps_data(:,i) + a1(:)
        CASE DEFAULT
            ! Pre-defined SALSA outputs
            fail = calc_user_data(user_ps_list(i),a,mask)
            IF (fail) THEN
                WRITE(*,*)" Error: failed to calculate '"//TRIM(user_ps_list(i))//"' for ps output!"
                STOP
            ENDIF
            ! Averaging
            CALL get_avg3(nzp,nxp,nyp,a,a1,cond=mask)
            user_ps_data(:,i) = user_ps_data(:,i) + a1(:)
        END SELECT
    ENDDO
    !
  end subroutine ps_user_stats

  ! This is the same for bin dependent profile outputs
  subroutine ps_user_bin_stats()
    USE mo_submctl, ONLY : find_spec_id,fn2a,in2b,fn2b,fnp2a,inp2b,fnp2b,nprc,nsnw
    use grid, ONLY : nzp,nxp,nyp,binSpecMixrat,meanRadius,a_naerop,a_ncloudp,a_nprecpp,a_nicep,a_nsnowp
    INTEGER :: i, j, k, ind, bb
    CHARACTER(LEN=7) :: short_name, nam
    REAL :: a1(nzp,nxp,nyp), w(nzp,nxp,nyp), col(nzp)
    !
    ind=1
    DO j=1,nv2_bin
        ! Variable name
        short_name=s2_bin(j)
        i=LEN_TRIM(short_name)
        !
        w(:,:,:)=1.
        !
        ! Species name or N for number (radius is not valid)
        nam=short_name(3:i-2)
        IF (nam=='N  ') THEN
            ! Bin number concentration
            k=0
        ELSEIF (nam=='R  ' .OR. nam=='Rw ') THEN
            ! Bin radius
            k=-1
        ELSE
            ! Species name
            k=find_spec_id(nam)
            IF (k<1) THEN
                WRITE(*,*) 'Examining output '//TRIM(short_name)
                STOP 'Error in subroutine ps_user_bin_stats!'
            ENDIF
        ENDIF
        !
        ! Output species & bin: aa/ab/ca/cb/rt/ia/ib/st
        SELECT case (short_name(i-1:i))
        CASE('aa')
            ! Aerosol a-bins
            DO bb = 1,fn2a
                IF (k>0) THEN
                    ! Bin mixing ratios for all species
                    CALL binSpecMixrat('aerosol',k,bb,a1)
                ELSEIF (k<0) THEN
                    ! Bin radius (number mean)
                    CALL meanRadius('aerosol','x',a1,ibin=bb)
                    w(:,:,:)=a_naerop(:,:,:,bb)
                ELSE
                    ! Bin number concentration
                    a1(:,:,:)=a_naerop(:,:,:,bb)
                ENDIF
                CALL get_avg3(nzp,nxp,nyp,a1,col,weight=w)
                svctr_bin(:,ind) = svctr_bin(:,ind) + col(:)
                ind=ind+1
            END DO
        CASE('ab')
            ! Aerosol b-bins
            DO bb = in2b,fn2b
                IF (k>0) THEN
                    ! Bin mixing ratios for all species
                    CALL binSpecMixrat('aerosol',k,bb,a1)
                ELSEIF (k<0) THEN
                    ! Bin radius (number mean)
                    CALL meanRadius('aerosol','x',a1,ibin=bb)
                    w(:,:,:)=a_naerop(:,:,:,bb)
                ELSE
                    ! Bin number concentration
                    a1(:,:,:)=a_naerop(:,:,:,bb)
                ENDIF
                CALL get_avg3(nzp,nxp,nyp,a1,col,weight=w)
                svctr_bin(:,ind) = svctr_bin(:,ind) + col(:)
                ind=ind+1
            END DO
        CASE('ca')
            ! Cloud a-bins
            DO bb = 1,fnp2a
                IF (k>0) THEN
                    CALL binSpecMixrat('cloud',k,bb,a1)
                ELSEIF (k<0) THEN
                    CALL meanRadius('cloud','x',a1,ibin=bb)
                    w(:,:,:)=a_ncloudp(:,:,:,bb)
                ELSE
                    a1(:,:,:)=a_ncloudp(:,:,:,bb)
                ENDIF
                CALL get_avg3(nzp,nxp,nyp,a1,col,weight=w)
                svctr_bin(:,ind) = svctr_bin(:,ind) + col(:)
                ind=ind+1
            END DO
        CASE('cb')
            ! Cloud b-bins
            DO bb = inp2b,fnp2b
                IF (k>0) THEN
                    CALL binSpecMixrat('cloud',k,bb,a1)
                ELSEIF (k<0) THEN
                    CALL meanRadius('cloud','x',a1,ibin=bb)
                    w(:,:,:)=a_ncloudp(:,:,:,bb)
                ELSE
                    a1(:,:,:)=a_ncloudp(:,:,:,bb)
                ENDIF
                CALL get_avg3(nzp,nxp,nyp,a1,col,weight=w)
                svctr_bin(:,ind) = svctr_bin(:,ind) + col(:)
                ind=ind+1
            END DO
        CASE('rt','pt')
            ! Rain
            DO bb = 1,nprc
                IF (k>0) THEN
                    CALL binSpecMixrat('precp',k,bb,a1)
                ELSEIF (k<0) THEN
                    CALL meanRadius('precp','x',a1,ibin=bb)
                    w(:,:,:)=a_nprecpp(:,:,:,bb)
                ELSE
                    a1(:,:,:)=a_nprecpp(:,:,:,bb)
                ENDIF
                CALL get_avg3(nzp,nxp,nyp,a1,col,weight=w)
                svctr_bin(:,ind) = svctr_bin(:,ind) + col(:)
                ind=ind+1
            END DO
        CASE('ia')
            ! Ice a-bins
            DO bb = 1,fnp2a
                IF (k>0) THEN
                    CALL binSpecMixrat('ice',k,bb,a1)
                ELSEIF (k<0) THEN
                    CALL meanRadius('ice','x',a1,ibin=bb)
                    w(:,:,:)=a_nicep(:,:,:,bb)
                ELSE
                    a1(:,:,:)=a_nicep(:,:,:,bb)
                ENDIF
                CALL get_avg3(nzp,nxp,nyp,a1,col,weight=w)
                svctr_bin(:,ind) = svctr_bin(:,ind) + col(:)
                ind=ind+1
            END DO
        CASE('ib')
            ! Ice b-bins
            DO bb = inp2b,fnp2b
                IF (k>0) THEN
                    CALL binSpecMixrat('ice',k,bb,a1)
                ELSEIF (k<0) THEN
                    CALL meanRadius('ice','x',a1,ibin=bb)
                    w(:,:,:)=a_nicep(:,:,:,bb)
                ELSE
                    a1(:,:,:)=a_nicep(:,:,:,bb)
                ENDIF
                CALL get_avg3(nzp,nxp,nyp,a1,col,weight=w)
                svctr_bin(:,ind) = svctr_bin(:,ind) + col(:)
                ind=ind+1
            END DO
        CASE('st')
            ! Snow
            DO bb = 1,nsnw
                IF (k>0) THEN
                    CALL binSpecMixrat('snow',k,bb,a1)
                ELSEIF (k<0) THEN
                    CALL meanRadius('snow','x',a1,ibin=bb)
                    w(:,:,:)=a_nsnowp(:,:,:,bb)
                ELSE
                    a1(:,:,:)=a_nsnowp(:,:,:,bb)
                ENDIF
                CALL get_avg3(nzp,nxp,nyp,a1,col,weight=w)
                svctr_bin(:,ind) = svctr_bin(:,ind) + col(:)
                ind=ind+1
            END DO
        case default
            WRITE(*,*) 'Examining output '//TRIM(short_name)
            STOP 'Error in subroutine ps_user_bin_stats!'
        END SELECT
    ENDDO
    !
  END subroutine ps_user_bin_stats
  !
  !
  !
  subroutine comp_tke(n1,n2,n3,dzm,th00,u,v,w,s)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dzm(n1),th00,u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
    real, intent (inout) :: s(n1,n2,n3)

    integer :: k,kp1
    real    :: x1(n1), x2(n1)

    !
    ! ------
    ! Calculates buoyancy forcing
    !
    call get_buoyancy(n1,n2,n3,s,w,th00)
    !
    ! ------
    ! Estimates shear component of TKE budget
    !
    call get_shear(n1,n2,n3,u,v,w,dzm)
    !
    ! ------
    ! Calculates horizontal variances and resolved TKE
    !
    call get_avg3(n1,n2,n3,u**2,x1)
    call get_avg3(n1,n2,n3,u,x2)
    svctr(:,14) = svctr(:,14) + (x1(:)-x2(:)**2) ! Variance
    tke_res(:)=(x1(:)-x2(:)**2)

    call get_avg3(n1,n2,n3,v**2,x1)
    call get_avg3(n1,n2,n3,v,x2)
    svctr(:,15) = svctr(:,15) + (x1(:)-x2(:)**2)
    tke_res(:)  = tke_res(:) + (x1(:)-x2(:)**2)

    call get_avg3(n1,n2,n3,w**2,x1)
    svctr(:,16) = svctr(:,16)+x1(:) ! Raw moment
    do k=1,n1
       kp1 = min(k+1,n1)
       tke_res(k)  = 0.5*(0.5*(tke_res(k)+tke_res(kp1)) + x1(k))
    end do
    if (nsmp == 0) tke0(:) = tke_res(:)

  end subroutine comp_tke
  !
  ! ---------------------------------------------------------------------
  ! get_buoyancy:  estimates buoyancy production term in tke budget
  !
  subroutine get_buoyancy(n1,n2,n3,b,w,th00)

    use defs, only : g

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: w(n1,n2,n3),th00
    real, intent(inout) :: b(n1,n2,n3)

    integer :: i,j,k,kp1

    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             kp1 = min(k+1,n1)
             b(k,i,j) = (b(k,i,j) + b(kp1,i,j))
          end do
       end do
    end do
    call get_cor3(n1,n2,n3,b,w,wtv_res)
    do k=1,n1
       svctr(k,35) = svctr(k,35) + wtv_res(k)
       wtv_res(k) = wtv_res(k) * th00/g
    end do

  end subroutine get_buoyancy
  !
  ! ---------------------------------------------------------------------
  ! get_shear:  estimates shear production term in tke budget
  !
  subroutine get_shear(n1,n2,n3,u,v,w,dzm)

    integer, intent(in) :: n3,n2,n1
    real, intent(in)    :: w(n1,n2,n3),dzm(n1),u(n1,n2,n3),v(n1,n2,n3)

    real :: ub(n1), vb(n1)
    integer i,j,k
    real fact, uw_shear(n1), vw_shear(n1)

    fact = 0.25/float((n2-4)*(n3-4))

    call get_avg3(n1,n2,n3,u,ub)
    call get_avg3(n1,n2,n3,v,vb)

    uw_shear(:)=0.
    vw_shear(:)=0.
    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             uw_shear(k) = uw_shear(k) -(u(k,i,j)-ub(k))*fact*(               &
                  (w(k,i,j)  +w(k,i+1,j)  )*(ub(k+1)-ub(k)  )*dzm(k) +        &
                  (w(k-1,i,j)+w(k-1,i+1,j))*(ub(k)  -ub(k-1))*dzm(k-1))
             vw_shear(k) = vw_shear(k) -(v(k,i,j)-vb(k))*fact*(               &
                  (w(k,i,j)  +w(k,i,j+1)  )*(vb(k+1)-vb(k)  )*dzm(k) +        &
                  (w(k-1,i,j)+w(k-1,i,j+1))*(vb(k)  -vb(k-1))*dzm(k-1))
          end do
       end do
    end do

    ! Global
    CALL get_pustat_vector('avg',n1,uw_shear)
    CALL get_pustat_vector('avg',n1,vw_shear)

    svctr(:,48) = svctr(:,48)+uw_shear(:)
    svctr(:,36) = svctr(:,36)+uw_shear(:)+vw_shear(:)

  end subroutine get_shear
  !
  ! ----------------------------------------------------------------------
  ! Subroutine write_ts: writes the statistics file
  !
  subroutine write_ts

    use netcdf
    use mpi_interface, only : myid

    integer :: iret, n, VarID

    do n=1,nvar1
       iret = nf90_inq_varid(ncid1, s1(n), VarID)
       IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, ssclr(n), start=(/nrec1/))
    end do
    ssclr(:) = 0.

    IF (level == 0) THEN
       DO n = 1,nv1_ice
          iret = nf90_inq_varid(ncid1, s1_ice(n), VarID)
          IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, ssclr_ice(n), start=(/nrec1/))
       END DO
       ssclr_ice(:) = 0.
    END IF

    IF (level >= 4) THEN
       DO n = 1,nv1_lvl4
          iret = nf90_inq_varid(ncid1, s1_lvl4(n), VarID)
          IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, ssclr_lvl4(n), start=(/nrec1/))
       END DO
       ssclr_lvl4(:) = 0.
    END IF

    IF (level >= 5) THEN
       DO n = 1,nv1_lvl5
          iret = nf90_inq_varid(ncid1, s1_lvl5(n), VarID)
          IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, ssclr_lvl5(n), start=(/nrec1/))
       END DO
       ssclr_lvl5(:) = 0.
    END IF

    IF (nv1_rem>0) THEN
       DO n = 1,nv1_rem
          iret = nf90_inq_varid(ncid1, s1_rem(n), VarID)
          IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, ssclr_rem(n), start=(/nrec1/))
       END DO
       ssclr_rem(:) = 0.
    END IF

    ! User-selected process rate outputs
    IF (nv1_proc>0) THEN
       DO n=1,nv1_proc
          iret = nf90_inq_varid(ncid1, out_ts_list(n), VarID)
          IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, out_ts_data(n), start=(/nrec1/))
        ENDDO
        out_ts_data(:) = 0.
    ENDIF

    ! Other user-selected outputs
    IF (nv1_user>0) THEN
       DO n=1,nv1_user
          iret = nf90_inq_varid(ncid1, user_ts_list(n), VarID)
          IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, user_ts_data(n), start=(/nrec1/))
        ENDDO
        user_ts_data(:) = 0.
    ENDIF

    if (myid==0) print "(/' ',12('-'),'   Record ',I4,' to time series')",nrec1

    iret = nf90_sync(ncid1)
    nrec1 = nrec1 + 1

  end subroutine write_ts
  !
  ! ----------------------------------------------------------------------
  ! Subroutine write_ps: writes the time averaged elements of the
  ! statistics file
  !
  subroutine  write_ps(n1,dn0,u0,v0,zm,zt,time)

    use netcdf
    use defs, only : alvl, cp
    USE mo_submctl, ONLY : in1a,in2a,fn2a, &
                               aerobins,precpbins,snowbins, &
                               nout_cld, cldbinlim, nout_ice, icebinlim
    USE grid, ONLY : nprc,nsnw
    use mpi_interface, only : myid

    integer, intent (in) :: n1
    real, intent (in)    :: time
    real, intent (in)    :: dn0(n1), u0(n1), v0(n1), zm(n1), zt(n1)

    integer :: iret, VarID, k, n, kp1
    REAL, ALLOCATABLE :: rmid(:)

    lsttm = time
    do k=1,n1
       kp1 = min(n1,k+1)
       svctr(k,20) = (svctr(k,20)+svctr(k,21))*cp
       svctr(k,22) = svctr(k,22)+svctr(k,23)
       svctr(k,24) = svctr(k,24)+svctr(k,25)
       svctr(k,26) = svctr(k,26)+svctr(k,27)
       svctr(k,53) = (svctr(k,53)+svctr(k,54))*alvl
       svctr(k,21) = svctr(k,21)*cp
       svctr(k,54) = svctr(k,54)*alvl
       svctr(k,62) = svctr(k,62)*alvl
       svctr(k,37) = svctr(k,44) + svctr(k,47) +(                         &
            +svctr(k,45) + svctr(kp1,45) + svctr(k,46) + svctr(kp1,46)    &
            +svctr(k,42) + svctr(kp1,42) + svctr(k,43) + svctr(kp1,43)    &
            -svctr(k,36) - svctr(kp1,36)   )*0.5
       if (lsttm>fsttm) then
          svctr(k,49) = (tke_res(k) - tke0(k))/(lsttm-fsttm)
       else
          svctr(k,49) = 0.
       end if

       svctr(k,10:nvar2) = svctr(k,10:nvar2)/nsmp
    end do

    ! Time
    iret = nf90_inq_VarID(ncid2, s2(1), VarID)
    iret = nf90_put_var(ncid2, VarID, time, start=(/nrec2/))
    ! Constants (grid dimensions, bin limits,...)
    if (nrec2 == 1) then
       iret = nf90_inq_varid(ncid2, s2(2), VarID)
       iret = nf90_put_var(ncid2, VarID, zt, start = (/nrec2/))
       iret = nf90_inq_varid(ncid2, s2(3), VarID)
       iret = nf90_put_var(ncid2, VarID, zm, start = (/nrec2/))
       iret = nf90_inq_varid(ncid2, s2(4), VarID)
       iret = nf90_put_var(ncid2, VarID, dn0, start = (/nrec2/))
       iret = nf90_inq_varid(ncid2, s2(5), VarID)
       iret = nf90_put_var(ncid2, VarID, u0, start = (/nrec2/))
       iret = nf90_inq_varid(ncid2, s2(6), VarID)
       iret = nf90_put_var(ncid2, VarID, v0, start = (/nrec2/))
       ! Juha: For SALSA
       IF (level >= 4) THEN
          iret = nf90_inq_varid(ncid2,'B_Rd12a',VarID) ! 1a+2a
          IF (iret == NF90_NOERR) &
                iret = nf90_put_var(ncid2,VarID,aerobins(in1a:fn2a),start=(/nrec2/))
          iret = nf90_inq_varid(ncid2,'B_Rd2ab',VarID) ! 2a = 2b
          IF (iret == NF90_NOERR) &
                iret = nf90_put_var(ncid2,VarID,aerobins(in2a:fn2a),start=(/nrec2/))
          iret = nf90_inq_varid(ncid2,'B_Rwprc',VarID)
          IF (iret == NF90_NOERR) &
                iret = nf90_put_var(ncid2,VarID,precpbins(1:nprc),start=(/nrec2/))
          iret = nf90_inq_varid(ncid2,'B_Rwsnw',VarID)
          IF (iret == NF90_NOERR) &
                iret = nf90_put_var(ncid2,VarID,snowbins(1:nsnw),start=(/nrec2/))
          ! Histograms
          iret = nf90_inq_varid(ncid2,'P_hRc',VarID)
          IF (iret == NF90_NOERR) THEN
              ALLOCATE(rmid(nout_cld))
              rmid(:)=0.5*(cldbinlim(1:nout_cld)+cldbinlim(2:nout_cld+1))
              iret = nf90_put_var(ncid2,VarID,rmid,start=(/nrec2/))
              DEALLOCATE(rmid)
          END IF
          iret = nf90_inq_varid(ncid2,'P_hRi',VarID)
          IF (iret == NF90_NOERR) THEN
              ALLOCATE(rmid(nout_ice))
              rmid(:)=0.5*(icebinlim(1:nout_ice)+icebinlim(2:nout_ice+1))
              iret = nf90_put_var(ncid2,VarID,rmid,start=(/nrec2/))
              DEALLOCATE(rmid)
          END IF
       END IF
    end if

    iret = nf90_inq_VarID(ncid2, s2(7), VarID)
    iret = nf90_put_var(ncid2, VarID, fsttm, start=(/nrec2/))
    iret = nf90_inq_VarID(ncid2, s2(8), VarID)
    iret = nf90_put_var(ncid2, VarID, lsttm, start=(/nrec2/))
    iret = nf90_inq_VarID(ncid2, s2(9), VarID)
    iret = nf90_put_var(ncid2, VarID, nsmp,  start=(/nrec2/))

    do n=10,nvar2
       iret = nf90_inq_varid(ncid2, s2(n), VarID)
       IF (iret == NF90_NOERR) &
          iret = nf90_put_var(ncid2,VarID,svctr(:,n), start=(/1,nrec2/), count=(/n1,1/))
    end do

    IF (level==0) THEN
       svctr_ice(:,:) = svctr_ice(:,:)/nsmp
       do n=1,nv2_ice
          iret = nf90_inq_varid(ncid2, s2_ice(n), VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_ice(:,n), start=(/1,nrec2/), count=(/n1,1/))
       end do
       svctr_ice(:,:) = 0.
    ENDIF

    IF (level >= 4 .AND. nv2_lvl4>0) THEN
       ! SALSA level 4
       svctr_lvl4(:,:) = svctr_lvl4(:,:)/nsmp
       DO n = 1,nv2_lvl4
          iret = nf90_inq_varid(ncid2,s2_lvl4(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_lvl4(:,n), start=(/1,nrec2/), count=(/n1,1/))
       END DO
       svctr_lvl4(:,:) = 0.
    END IF

    IF (level >= 5 .AND. nv2_lvl5>0) THEN
       ! SALSA level 5
       svctr_lvl5(:,:) = svctr_lvl5(:,:)/nsmp
       DO n = 1,nv2_lvl5
          iret = nf90_inq_varid(ncid2,s2_lvl5(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_lvl5(:,n), start=(/1,nrec2/), count=(/n1,1/))
       END DO
       svctr_lvl5(:,:) = 0.
    END IF

    IF (nv2_bin>0) THEN
       ! SALSA binned data
       svctr_bin(:,:) = svctr_bin(:,:)/nsmp
       k=0
       DO n = 1,nv2_bin
          kp1=nv2_bin_len(n)
          iret = nf90_inq_varid(ncid2,s2_bin(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_bin(:,k+1:k+kp1), start=(/1,1,nrec2/), count=(/n1,kp1,1/))
          k=k+kp1
       END DO
       svctr_bin(:,:) = 0.
    ENDIF

    IF (level >= 4 .AND. nout_cld>0) THEN
        iret = nf90_inq_varid(ncid2, s2_CldHist(1), VarID)
        IF (iret == NF90_NOERR) &
            iret = nf90_put_var(ncid2,VarID,svctr_ch(:,:,1), start=(/1,1,nrec2/), count=(/n1,nout_cld,1/))
        iret = nf90_inq_varid(ncid2, s2_CldHist(2), VarID)
        IF (iret == NF90_NOERR) &
            iret = nf90_put_var(ncid2,VarID,svctr_ch(:,:,2), start=(/1,1,nrec2/), count=(/n1,nout_cld,1/))
        svctr_ch(:,:,:) = 0.
    END IF

    IF (level >= 5 .AND. nout_ice>0) THEN
        iret = nf90_inq_varid(ncid2, s2_IceHist(1), VarID)
        IF (iret == NF90_NOERR) &
            iret = nf90_put_var(ncid2,VarID,svctr_ih(:,:,1), start=(/1,1,nrec2/), count=(/n1,nout_ice,1/))
        iret = nf90_inq_varid(ncid2, s2_IceHist(2), VarID)
        IF (iret == NF90_NOERR) &
            iret = nf90_put_var(ncid2,VarID,svctr_ih(:,:,2), start=(/1,1,nrec2/), count=(/n1,nout_ice,1/))
        svctr_ih(:,:,:) = 0.
    END IF

    IF (nv2_proc>0) THEN
        DO n=1,nv2_proc
            iret = nf90_inq_varid(ncid2, out_ps_list(n), VarID)
            IF (iret == NF90_NOERR) THEN
                ! Instantaneous, so no need to divide by nsmp: out_ps_data(:,n)=out_ps_data(:,n)/nsmp
                iret = nf90_put_var(ncid2,VarID,out_ps_data(:,n), start=(/1,nrec2/),count=(/n1,1/))
            ENDIF
        ENDDO
        out_ps_data(:,:) = 0.
    ENDIF

    IF (nv2_user>0) THEN
        DO n=1,nv2_user
            iret = nf90_inq_varid(ncid2, user_ps_list(n), VarID)
            IF (iret == NF90_NOERR) THEN
                user_ps_data(:,n)=user_ps_data(:,n)/nsmp
                iret = nf90_put_var(ncid2,VarID,user_ps_data(:,n),start=(/1,nrec2/),count=(/n1,1/))
            ENDIF
        ENDDO
        user_ps_data(:,:) = 0.
    ENDIF

    if (myid==0) print "(/' ',12('-'),'   Record ',I4,' to profiles')",nrec2

    iret  = nf90_sync(ncid2)
    nrec2 = nrec2+1
    nsmp  = 0.

    svctr(:,:) = 0.

  end subroutine write_ps
  !
  ! ----------------------------------------------------------------------
  !
  ! subroutine: fills scalar array based on variable name.
  !
  subroutine fill_scalar(xval,vname,op,wg)
    ! Inputs
    real, intent(in) :: xval
    character(len=7), intent (in) :: vname
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: op ! Optional operation; default='avg'
    REAL, INTENT(IN), OPTIONAL :: wg ! Optional weight; default=none
    ! Local
    INTEGER :: i
    CHARACTER(LEN=3) :: myop
    !
    myop='avg'
    IF (PRESENT(op)) myop=op
    !
    ! s1(1:nvar1)
    DO i=1,nvar1
        IF ( vname == s1(i) ) THEN
            ! Output array found, calculate result over all PUs
            IF (present(wg)) THEN
                ssclr(i) = get_pustat_scalar(myop,xval,wg)
            ELSE
                ssclr(i) = get_pustat_scalar(myop,xval)
            ENDIF
            RETURN
        ENDIF
    ENDDO
    !
    ! s1_lvl4(1:nv1_lvl4)
    !DO i=1,nv1_lvl4
    !    IF ( vname == s1_lvl4(i) ) THEN
    !        ! Output array found, calculate result over all PUs
    !        IF (present(wg)) THEN
    !            ssclr(i) = get_pustat_scalar(myop,xval,wg)
    !        ELSE
    !            ssclr(i) = get_pustat_scalar(myop,xval)
    !        ENDIF
    !        RETURN
    !    ENDIF
    !ENDDO
    !
    ! Could add s1_ice(1:nv1_ice), s1_lvl5(1:nv1_lvl5),
    ! user_ts_list(1:nv1_user) and out_ts_list(1:nv1_proc)
  end subroutine fill_scalar
  !
  ! --------------------------------------------------------------------------
  ! SGSFLXS: estimates the sgs rl and tv flux from the sgs theta_l and sgs r_t
  ! fluxes
  !
  subroutine sgsflxs(n1,n2,n3,level,rl,rv,th,flx,type)

    use defs, only : alvl, cp, rm, ep2

    integer, intent(in) :: n1,n2,n3,level
    real, intent(in)    :: rl(n1,n2,n3),rv(n1,n2,n3)
    real, intent(in)    :: th(n1,n2,n3),flx(n1,n2,n3)
    character (len=2)   :: type

    integer :: k,i,j
    real    :: rnpts      ! reciprical of number of points and
    real    :: fctl, fctt ! factors for liquid (l) and tv (t) fluxes

    if (type == 'tl') then
       wrl_sgs(:) = 0.
       wtv_sgs(:) = 0.
    end if
    rnpts = 1./real((n2-4)*(n3-4))
    !
    ! calculate fluxes assuming liquid water.
    !
       do j = 3,n3-2
          do i = 3,n2-2
             do k = 1,n1-1
                if (rl(k+1,i,j) > 0.) then
                   fctt = rnpts*(1. + rv(k,i,j)*(1.+ep2 +ep2*rv(k,i,j)*alvl   &
                        /(rm*th(k,i,j))))                                     &
                        /(1.+(rv(k,i,j)*(alvl/th(k,i,j))**2)/(rm*cp))
                   select case (type)
                   case ('tl')
                      fctl =-rnpts/(rm*th(k,i,j)**2/(rv(k,i,j)*alvl)+alvl/cp)
                   case ('rt')
                      fctl =rnpts/(1.+(rv(k,i,j)*alvl**2)/(cp*rm*th(k,i,j)**2))
                      fctt = (alvl*fctt/cp - th(k,i,j)*rnpts)
                   end select
                   wrl_sgs(k) = wrl_sgs(k) + fctl*flx(k,i,j)
                   wtv_sgs(k) = wtv_sgs(k) + fctt*flx(k,i,j)
                else
                   select case (type)
                   case ('tl')
                      fctt = rnpts*(1. + ep2*rv(k,i,j))
                   case ('rt')
                      fctt = rnpts*(ep2*th(k,i,j))
                   end select
                   wtv_sgs(k) = wtv_sgs(k) + fctt*flx(k,i,j)
                end if
             end do
          end do
       end do

    ! Global
    if (type == 'rt') then
       CALL get_pustat_vector('avg', n1, wtv_sgs)
       CALL get_pustat_vector('avg', n1, wrl_sgs)
    ENDIF

  end subroutine sgsflxs
  !
  !---------------------------------------------------------------------
  ! subroutine updtst: updates appropriate statistical arrays
  !
  subroutine updtst(n1,values,ic,nam)

    integer, intent(in)            :: n1,ic
    real, intent (in)              :: values(n1)
    character (len=7), intent (in) :: nam

    integer :: nn
    REAL :: tmp(n1)

    DO nn=1,nvar2
        IF (nam==s2(nn)) THEN
            if (ic == 0) svctr(:,nn)=0.
            ! Global
            tmp(:)=values(:)
            CALL get_pustat_vector('avg',n1,tmp)
            svctr(:,nn)=svctr(:,nn)+tmp(:)
            RETURN
        ENDIF
    ENDDO

  end subroutine updtst

  ! ----------------------------------------------------------------------
  ! Subroutines for calculating radius for Seifert & Beheng hydrometeors
  ! Local parameters:
  !     REAL, DIMENSION(0:4) :: x_min=1e-20, x_max=1., a_geo=(6./(pi*rowt))**(1./3.), b_geo=1./3., eps=1e-20
  ! a) Set parameters
  SUBROUTINE setSBradius(idno,my_x_min,my_x_max,my_a_geo,my_b_geo,my_eps)
    INTEGER, INTENT(IN) :: idno
    REAL, INTENT(IN) :: my_x_min,my_x_max,my_a_geo,my_b_geo,my_eps
    x_min(idno)=my_x_min ! x_min < mass < x_max
    x_max(idno)=my_x_max
    a_geo(idno)=my_a_geo ! D=a_geo*(mass/num)**b_geo
    b_geo(idno)=my_b_geo
    eps(idno)=my_eps
  END SUBROUTINE setSBradius
  ! b) Get radius
  SUBROUTINE getSBradius(n1,n2,n3,numc,mass,idno,zrad)
    INTEGER, INTENT(IN) :: n1,n2,n3,idno
    REAL, DIMENSION(n1,n2,n3), INTENT(IN) :: numc,mass
    REAL, INTENT(OUT) :: zrad(n1,n2,n3)
    ! Local parameters
    INTEGER :: i, j, k
    REAL :: x_c
    ! Calculate
    DO j = 3,n3-2
        DO i = 3,n2-2
            DO k = 1,n1
                ! Mean mass
                x_c = MIN(MAX(mass(k,i,j)/(numc(k,i,j)+eps(idno)),x_min(idno)),x_max(idno))
                ! Radius
                zrad(k,i,j) = 0.5*a_geo(idno)*x_c**b_geo(idno)
            END DO
        END DO
    END DO
  END SUBROUTINE getSBradius

  ! -------------------------------------------------------------------------
  ! Produce 3D outputs for further averaging based on user provided variable name:
  !     name=<condition>//<species>//<bin>
  !       condition: tc=total, ic=in cloud, ir=in rain, ii=in ice, is=in snow
  !       species: aerosol (SO4, OC, BC, ..) or gas name for mass concentration,
  !                R for mean radius and N for number concentration
  !       bin: aa,ab,at,ca,cb,ct,rt,ia,ib,it,st,gt (aerosol/cloud/rain/ice/snow/gas a-bins/b-bins/total)
  LOGICAL FUNCTION calc_user_data(short_name,res,mask,is_mass)
    use grid, ONLY : nzp,nxp,nyp,bulkNumc,bulkMixrat,meanRadius,a_gaerop
    USE mo_submctl, ONLY : find_gas_id, find_spec_id
    CHARACTER(LEN=7), INTENT(IN) :: short_name ! Variable name
    REAL, INTENT(out) :: res(nzp,nxp,nyp)      ! Output data
    LOGICAL, INTENT(out) :: mask(nzp,nxp,nyp)  ! ... and mask
    LOGICAL, INTENT(out), OPTIONAL :: is_mass  ! True if the ouput is mass concentration
    !
    CHARACTER(LEN=7) :: spec, bin, nam
    LOGICAL :: mass
    INTEGER :: i, k
    !
    ! Return .true. if failed
    calc_user_data=.TRUE.
    !
    ! String length (must be at least 5, e.g. tcNaa)
    i=LEN(TRIM(short_name))
    IF (i<5) RETURN
    !
    ! Condition (2 chars): tc=total, ic=in-cloud, ir=in rain, ii=in ice, is=in snow
    SELECT CASE (short_name(1:2))
    CASE('tc','ia') ! Total or aerosol
        mask=.TRUE.
    CASE('ic')
        mask=cloudmask
    CASE('ir')
        mask=rainmask
    CASE('ii')
        mask=icemask
    CASE('is')
        mask=snowmask
    case default
        RETURN
    END SELECT
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
    mass=.TRUE. ! Output is mass concentration?
    IF (spec=='gas') THEN
        ! Just mass mixing ratio, no radius or number
        k=find_gas_id(nam)
        IF (k<1) RETURN
        res(:,:,:)=a_gaerop(:,:,:,k)
    ELSEIF (nam=='N  ') THEN
        CALL bulkNumc(spec,bin,res)
        mass=.FALSE.
    ELSEIF (nam=='R  ' .OR. nam=='Rw ') THEN
        CALL meanRadius(spec,bin,res)
        ! Mean radius not defined (arbitarily set to zero) when there are no particles, so these must be ignored
        WHERE(res<1e-20) mask=.FALSE.
        mass=.FALSE.
    ELSE
        k=find_spec_id(nam)
        IF (k<1) RETURN
        CALL bulkMixrat(k,spec,bin,res)
    ENDIF
    !
    ! Optional output
    IF (present(is_mass)) is_mass=mass
    !
    ! All done
    calc_user_data = .FALSE.
  END FUNCTION calc_user_data

  ! -------------------------------------------------------------------------
  ! Produce user-requested (out_ts_list, out_ps_list, out_cs_list and out_an_list) outputs that give
  ! information about the effects of LES and microphysics function calls on prognostic variables.
  ! These calculations are instantaneous and the outputs describe the rate of change (tendency).
  !
  ! a) Data from SALSA or S&B microphysics is ready for calculations
  SUBROUTINE mcrp_var_save()
    USE grid, ONLY : nxp, nyp, nzp
    IMPLICIT NONE
    INTEGER :: i
    !
    ! 4D array out_mcrp_data contains 3D data arrays whose names are
    ! specified in the out_mcrp_list containing out_mcrp_nout items
    DO i=1,out_mcrp_nout
        ! Calculate different outputs
        CALL scalar_rate_stat(out_mcrp_list(i),nzp,nxp,nyp,out_mcrp_data(:,:,:,i))
    ENDDO
    !
  END SUBROUTINE mcrp_var_save
  !
  ! b) LES is has data (tendencies) for outputs (several function calls from t_step in step.f90)
  SUBROUTINE les_rate_stats(prefix)
    USE grid, ONLY : out_an_list, nxp, nyp, nzp, level, &
                     a_tt, a_rt,a_rpt, a_npt
    use defs, only : cp
    IMPLICIT NONE
    ! Input
    character (len=4), intent (in) :: prefix ! 'srfc', 'diff', 'forc', 'mcrp', 'sedi', 'advf',...
    ! Local
    INTEGER :: i
    !
    ! a) LES variables
    ! [Ice-liquid water] potential temperature tendency:
    ! convert to heating rate [W/kg] by multiplying with heat capacity
    CALL scalar_rate_stat(prefix//'_tt',nzp,nxp,nyp,a_tt,factor=cp)
    !
    ! Concentrations
    IF (level==3 .OR. level==0) THEN
        ! Level 3: total water (a_rt) and rain water mass and droplet number
        CALL scalar_rate_stat(prefix//'_rt',nzp,nxp,nyp,a_rt)
        CALL scalar_rate_stat(prefix//'_nr',nzp,nxp,nyp,a_npt)
        CALL scalar_rate_stat(prefix//'_rr',nzp,nxp,nyp,a_rpt)
    ELSEIF (level>=4) THEN
        ! Level 4 and 5: water vapor (a_rt)
        CALL scalar_rate_stat(prefix//'_rg',nzp,nxp,nyp,a_rt)
    ENDIF
    !
    ! b) LES variables related to SALSA
    IF (level<4) RETURN
    ! There are so many possible outputs, that it is faster to examine the requested outputs
    DO i=1,maxn_list
        ! The first four characters contain process name, which should match with the given prefix.
        ! The 5th character is just '_', but the 6th character indicates the output type (n=number,
        ! r=water mixing ratio or an integer i referring to i:th component in condensed or gas phase).
        ! The last (7th) character is phase (a, c, r, i, s or g).
        IF (prefix == out_an_list(i)(1:4)) THEN
            ! Analysis data
            CALL calc_salsa_rate(i,4,out_an_list(i)(6:6),out_an_list(i)(7:7))
        ENDIF
        IF (csflg .AND. prefix == out_cs_list(i)(1:4)) THEN
            ! Column statistics
            CALL calc_salsa_rate(i,3,out_cs_list(i)(6:6),out_cs_list(i)(7:7))
        ENDIF
        IF (prefix == out_ps_list(i)(1:4)) THEN
            ! Profiles
            CALL calc_salsa_rate(i,2,out_ps_list(i)(6:6),out_ps_list(i)(7:7))
        ENDIF
        IF (prefix == out_ts_list(i)(1:4)) THEN
            ! Time series
            CALL calc_salsa_rate(i,1,out_ts_list(i)(6:6),out_ts_list(i)(7:7))
        ENDIF
    ENDDO
    !
  END SUBROUTINE les_rate_stats

  ! This is the actual function for calculating level 4 and 5 outputs
  SUBROUTINE calc_salsa_rate(out_ind,out_id,tchar,pchar)
    USE grid, ONLY : dzt, a_dn, out_an_data, &
                     nxp, nyp, nzp, nbins, ncld, nprc, nice, nsnw, &
                     a_naerot, a_maerot, a_ncloudt, a_mcloudt, a_nprecpt, a_mprecpt, &
                     a_nicet,  a_micet, a_nsnowt, a_msnowt, a_gaerot
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: out_ind, out_id ! Index and type of output
    CHARACTER, INTENT(IN) :: tchar, pchar ! Type and phase characters
    ! Local
    INTEGER :: j, k
    REAL ::tend(nzp,nxp,nyp), tmp(nzp,nxp,nyp), col(nzp), area(nxp,nyp)
    !
    SELECT CASE (pchar)
        CASE ('a')
            ! Aerosol
            SELECT CASE (tchar)
                CASE ('n')
                   ! Number concentration
                   tend=SUM(a_naerot,DIM=4)
                CASE('r','1')
                   ! Water (component 1) mass mixing ratio
                   tend=SUM(a_maerot(:,:,:,1:nbins),DIM=4)
                CASE DEFAULT
                   ! Mass mixing ratio of component j
                   READ(UNIT=tchar,FMT='(I1)') j
                   tend=SUM(a_maerot(:,:,:,(j-1)*nbins+1:j*nbins),DIM=4)
                END SELECT
        CASE ('c')
            ! Cloud
            SELECT CASE (tchar)
                CASE ('n')
                   tend=SUM(a_ncloudt,DIM=4)
                CASE('r','1')
                   tend=SUM(a_mcloudt(:,:,:,1:ncld),DIM=4)
                CASE DEFAULT
                   READ(UNIT=tchar,FMT='(I1)') j
                   tend=SUM(a_mcloudt(:,:,:,(j-1)*ncld+1:j*ncld),DIM=4)
                END SELECT
        CASE ('r')
            ! Rain
            SELECT CASE (tchar)
                CASE ('n')
                   ! Number concentration
                   tend=SUM(a_nprecpt,DIM=4)
                CASE('r','1')
                   tend=SUM(a_mprecpt(:,:,:,1:nprc),DIM=4)
                CASE DEFAULT
                   READ(UNIT=tchar,FMT='(I1)') j
                   tend=SUM(a_mprecpt(:,:,:,(j-1)*nprc+1:j*nprc),DIM=4)
                END SELECT
        CASE ('i')
            ! Ice
            SELECT CASE (tchar)
                CASE ('n')
                   tend=SUM(a_nicet,DIM=4)
                CASE('r','1')
                   tend=SUM(a_micet(:,:,:,1:nice),DIM=4)
                CASE DEFAULT
                   READ(UNIT=tchar,FMT='(I1)') j
                   tend=SUM(a_micet(:,:,:,(j-1)*nice+1:j*nice),DIM=4)
                END SELECT
        CASE ('s')
            ! Snow
            SELECT CASE (tchar)
                CASE ('n')
                   ! Number concentration
                   tend=SUM(a_nsnowt,DIM=4)
                CASE('r','1')
                   ! Water (component 1) mass mixing ratio
                   tend=SUM(a_msnowt(:,:,:,1:nsnw),DIM=4)
                CASE DEFAULT
                   ! Mass mixing ratio of component j
                   READ(UNIT=tchar,FMT='(I1)') j
                   tend=SUM(a_msnowt(:,:,:,(j-1)*nsnw+1:j*nsnw),DIM=4)
                END SELECT
        CASE ('t')
            ! Total = aerosol+cloud+rain+ice+snow
            SELECT CASE (tchar)
                CASE ('n')
                   ! Number concentration
                   tend=SUM(a_naerot,DIM=4)+SUM(a_ncloudt,DIM=4)+SUM(a_nprecpt,DIM=4)
                   IF (level>4) tend=tend+SUM(a_nicet,DIM=4)+SUM(a_nsnowt,DIM=4)
                CASE('r','1')
                   ! Water (component 1) mass mixing ratio
                   tend=SUM(a_maerot(:,:,:,1:nbins),DIM=4)+SUM(a_mcloudt(:,:,:,1:ncld),DIM=4)+ &
                        SUM(a_mprecpt(:,:,:,1:nprc),DIM=4)
                   IF (level>4) tend=tend+SUM(a_micet(:,:,:,1:nice),DIM=4)+SUM(a_msnowt(:,:,:,1:nsnw),DIM=4)
                CASE DEFAULT
                   ! Mass mixing ratio of component j
                   READ(UNIT=tchar,FMT='(I1)') j
                   tend=SUM(a_maerot(:,:,:,(j-1)*nbins+1:j*nbins),DIM=4)+ &
                        SUM(a_mcloudt(:,:,:,(j-1)*ncld+1:j*ncld),DIM=4)+ &
                        SUM(a_mprecpt(:,:,:,(j-1)*nprc+1:j*nprc),DIM=4)
                   IF (level>4) tend=tend+SUM(a_micet(:,:,:,(j-1)*nice+1:j*nice),DIM=4)+ &
                                          SUM(a_msnowt(:,:,:,(j-1)*nsnw+1:j*nsnw),DIM=4)
                END SELECT
        CASE ('g')
            ! Gas: scalar so just mass mixing ratio of component j
            READ(UNIT=tchar,FMT='(I1)') j
            IF (j==0) j=10 ! Allow tenth gas
            tend=a_gaerot(:,:,:,j)
    END SELECT
    !
    ! Save the data; calculate averages when needed
    IF (out_id==4) THEN
        ! Analysis data as is (per mass of air)
        out_an_data(:,:,:,out_ind) = tend(:,:,:)
    ELSEIF (out_id==3) THEN
        ! Column outputs are integrals over vertical dimension (per volume of air)
        area(:,:) = 0.
        do k=2,nzp
            area(:,:)=area(:,:)+tend(k,:,:)*a_dn(k,:,:)/dzt(k)
        ENDDO
        out_cs_data(:,:,out_ind) = area(:,:)
    ELSEIF (out_id==2) THEN
        ! Profiles: averaged tendencies are weighted by air density
        tmp(:,:,:) = tend(:,:,:)*a_dn(:,:,:)
        ! Average over horizontal dimensions
        CALL get_avg3(nzp,nxp,nyp,tmp,col)
        ! Profile
        out_ps_data(:,out_ind) = col(:)
    ELSE
        ! Time series: averaged tendencies are weighted by air density
        tmp(:,:,:) = tend(:,:,:)*a_dn(:,:,:)
        ! Average over horizontal dimensions
        CALL get_avg3(nzp,nxp,nyp,tmp,col)
        ! Integrate over vertical dimension to get the domain mean
        out_ts_data(out_ind) = SUM( col(2:nzp)/dzt(2:nzp) )
    ENDIF
    !
  END SUBROUTINE calc_salsa_rate
  !
  ! Scalar 3D arrays
  SUBROUTINE scalar_rate_stat(vname,n1,n2,n3,tend,factor)
    USE grid, ONLY : dzt, a_dn, out_an_list, out_an_data
    IMPLICIT NONE
    ! Inputs
    character (len=7), intent (in) :: vname ! Variable name such as 'diag_rr'
    INTEGER, INTENT(IN) :: n1,n2,n3 ! Dimensions
    REAL, INTENT(in) :: tend(n1,n2,n3) ! Tendency array
    REAL, OPTIONAL, INTENT(in) :: factor ! Scaling factor for e.g. unit conversions
    ! Local
    INTEGER :: i, k
    REAL :: tmp(n1,n2,n3), col(n1), area(n2,n3), fact
    !
    fact=1.0
    IF (PRESENT(factor)) fact=factor
    !
    ! Is this output selected
    DO i=1,maxn_list
        IF ( vname == out_an_list(i) ) THEN
            ! Analysis data as is (per mass of air)
            out_an_data(:,:,:,i) = tend(:,:,:)*fact
        ENDIF
        IF ( vname == out_ts_list(i) .OR. vname == out_ps_list(i) ) THEN
            ! Averaged tendencies are weighted by air density
            tmp(:,:,:) = tend(:,:,:)*a_dn(:,:,:)
            ! Average over horizontal dimensions
            CALL get_avg3(n1,n2,n3,tmp,col)
            ! Profile
            IF ( vname == out_ps_list(i) ) out_ps_data(:,i) = col(:)*fact
            ! Integrate over vertical dimension to get the domain mean
            IF ( vname == out_ts_list(i) ) out_ts_data(i) = SUM( col(2:n1)/dzt(2:n1) )*fact
        ENDIF
        IF ( csflg .AND. vname == out_cs_list(i) ) THEN
            ! Column outputs are integrals over vertical dimension (per volume of air)
            area(:,:) = 0.
            do k=2,n1
                area(:,:)=area(:,:)+tend(k,:,:)*a_dn(k,:,:)/dzt(k)
            ENDDO
            out_cs_data(:,:,i) = area(:,:)*fact
        ENDIF
    ENDDO
  END SUBROUTINE scalar_rate_stat
  !
  !
  ! -------------------------------------------------------------------------
  ! Similar to updtst but intended for making temporal statistics of the
  ! aerosol removal processes.
  ! Juha Tonttila, FMI, 2015

  ! Jaakko Ahola, FMI, 2016
  ! Modified for ice'n'snow
  !
  ! -------------------------------------------------------------------------
  !
  SUBROUTINE acc_removal(n2,n3,n4,raer,rcld,rprc,rice,rsnw)
    USE grid, ONLY : nbins, ncld, nprc, nice, nsnw, nspec
    IMPLICIT NONE

    INTEGER, INTENT(in)           :: n2,n3,n4                     ! Grid dimensions
    REAL, INTENT(in)              :: raer(n2,n3,n4*nbins), &      ! Arrays containing the binned 2d-fields
                                     rcld(n2,n3,n4*ncld), &
                                     rprc(n2,n3,n4*nprc), &
                                     rice(n2,n3,n4*nice), &
                                     rsnw(n2,n3,n4*nsnw)

    INTEGER :: si, tt, end, str

    DO si = 1,nspec+1 ! Aerosol species and water
        tt=(si-1)*5+1

        ! Removal by sedimentation of aerosol
        str = (si-1)*nbins+1
        end = si*nbins
        ssclr_rem(tt) = get_avg2dh( n2,n3,SUM(raer(:,:,str:end),DIM=3) )
        tt=tt+1

        ! Removal by sedimentation of cloud droplets
        str = (si-1)*ncld+1
        end = si*ncld
        ssclr_rem(tt) = get_avg2dh( n2,n3,SUM(rcld(:,:,str:end),DIM=3) )
        tt=tt+1

        ! Removal by precipitation
        str = (si-1)*nprc+1
        end = si*nprc
        ssclr_rem(tt) = get_avg2dh( n2,n3,SUM(rprc(:,:,str:end),DIM=3) )
        tt=tt+1

        IF (level<5) CYCLE

        ! Removal by sedimentation of ice particles
        str = (si-1)*nice+1
        end = si*nice
        ssclr_rem(tt) = get_avg2dh( n2,n3,SUM(rice(:,:,str:end),DIM=3) )
        tt=tt+1

        ! Removal by snow
        str = (si-1)*nsnw+1
        end = si*nsnw
        ssclr_rem(tt) = get_avg2dh( n2,n3,SUM(rsnw(:,:,str:end),DIM=3) )
        tt=tt+1

    END DO

  END SUBROUTINE acc_removal
  !
  ! -------------------------------------------------------------------------
  !
  integer function close_stat()

    use netcdf

    close_stat = nf90_close(ncid1) + nf90_close(ncid2)
    IF (csflg) close_stat = close_stat + nf90_close(ncid3)

  end function close_stat


end module stat

