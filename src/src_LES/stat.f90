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
  use grid, only : level, lbinprof, lmixrprof, lremts, stat_micro, stat_micro_ts, stat_micro_ps, maxn_list
  use util, only : get_avg3, get_cor3, get_var3, get_avg_ts, &
                   get_avg2dh, get_3rd3, HistDistr

  implicit none
  private

  integer, parameter :: nvar1 = 35,               &
                        nv1_lvl4 = 13,            &
                        nv1_lvl5 = 18,            &
                        nvar2 = 96,               &
                        nv2_lvl4 = 12,            &
                        nv2_lvl5 = 9,             &
                        nv2_hist = 2,             &
                        nvar3 = 20,               &
                        nv3_lvl4 = 0,             &
                        nv3_lvl5 = 6

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

  ! New output variables:
  ! 1. Add them to one of the arrays below
  ! 2. Update the dimensions of the output arrays accordingly
  ! 3. Make a subroutine for accumulating the new data and to plane them
  !    in the output arrays
  ! 4. Add the new variables in ncio.f90 list of variables
  ! 5. Make sure stat_init is up to date (boolean arrays etc).


  character (len=7), save :: s1(nvar1)=(/                           &
       'time   ','cfl    ','maxdiv ','zi1_bar','zi2_bar','zi3_bar', & ! 1
       'vtke   ','sfcbflx','wmax   ','tsrf   ','ustar  ','shf_bar', & ! 7
       'lhf_bar','zi_bar ','lwp_bar','lwp_var','zc     ','zb     ', & !13
       'cfrac  ','lmax   ','albedo ','rwp_bar','prcp   ','pfrac  ', & !19
       'CCN    ','nrain  ','nrcnt  ','nccnt  ','prcp_bc','zcmn   ', & !25
       'zbmn   ','tkeint ','thl_int','wvp_bar','wvp_var'/),         & !31

       ! **** Bulk temporal statistics for SALSA ****
       s1_lvl4(nv1_lvl4) = (/       &
       'Nc_ic  ','Rc_ic  ','Nca_ica','Rca_ica','Ncb_icb','Rcb_icb', & !1
       'Na_int ','Ra_int ','Naa_int','Raa_int','Nab_int','Rab_int', & !7
       'SS_max '/), & !13

       s1_lvl5(nv1_lvl5) = (/  &
       'Ni_ii  ','Ri_ii  ','Nia_iia','Ria_iia','Nib_iib','Rib_iib', & ! 1-6
       'Ns_is  ','Rs_is  ',  & ! 7-8
       'iwp_bar','imax   ','nicnt  ',  & ! 9-11
       'swp_bar','smax   ','nscnt  ',  & ! 12-14
       'sfrac  ','sprcp  ','SSi_max',  & ! 15-17
       'thi_int'/), & ! 18

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
        'rt_cs2 ','rc_cs2 ','wtl_cs2','wtv_cs2','wrt_cs2','Nc     ', & ! 79
        'Nr     ','rr     ','rrate  ','evap   ','frc_prc','prc_prc', & ! 85
        'frc_ran','hst_srf','sw_up  ','sw_down','lw_up  ','lw_down'/), & ! 91, total 96

        ! **** BULK PROFILE OUTPUT FOR SALSA ****
        s2_lvl4(nv2_lvl4) = (/                                       &
        'P_Naa  ','P_Nab  ','P_Nca  ','P_Ncb  ','P_Np   ','P_Rwaa ', & ! 1-6
        'P_Rwab ','P_Rwca ','P_Rwcb ','P_Rwp  ','P_rv   ','crate  '/), & ! 7-12

        s2_lvl5(nv2_lvl5) = (/ &
        'P_Nia  ','P_Nib  ','P_Ns   ','P_Rwia ','P_Rwib ','P_Rws  ', & ! 1-6
        'irate  ','srate  ','thi    '/), & ! 7-9

        ! **** Cloud droplet and ice histograms
        s2_CldHist(nv2_hist) = (/'P_hNca ','P_hNcb '/), &
        s2_IceHist(nv2_hist) = (/'P_hNia ','P_hNib '/), &

        ! Column statistics
       s3(nvar3)=(/ &
       'time   ','xt     ','xm     ','yt     ','ym     ','lwp    ','rwp    ', & ! 1-7
       'Nc     ','Nr     ','nccnt  ','nrcnt  ','zb     ','zc     ','zi1    ', & ! 8-14
       'lmax   ','wmax   ','wbase  ','prcp   ','prcp_bc','albedo '/), & ! 15-20
       s3_lvl4(nv3_lvl4), & ! Not used
       s3_lvl5(nv3_lvl5) = (/ &
       'iwp    ','swp    ','Ni     ','Ns     ','nicnt  ','nscnt  '/) ! 1-6

  real, save, allocatable   :: tke_sgs(:), tke_res(:), tke0(:), wtv_sgs(:),  &
       wtv_res(:), wrl_sgs(:), thvar(:), svctr(:,:), ssclr(:),               &
       ! Additional ssclr and svctr for BULK SALSA output
       svctr_lvl4(:,:), ssclr_lvl4(:),                                       &
       svctr_lvl5(:,:), ssclr_lvl5(:),                                       &
       ssclr_mixr(:), ssclr_rem(:), svctr_mixr(:,:),                         &
       ! Additional ssclr and svctr for BINNED SALSA output.
       svctr_aa(:,:,:), svctr_ca(:,:,:), svctr_p(:,:,:),                     &
       svctr_ab(:,:,:), svctr_cb(:,:,:),                                     &
       svctr_ia(:,:,:), svctr_ib(:,:,:), svctr_s(:,:,:),                     &
       ! Cloud and ice histograms
       svctr_ch(:,:,:), svctr_ih(:,:,:),                                     &
       ! Gas phase concentrations
       svctr_gas(:,:), ssclr_gas(:),                                         &
       ! User-selected process rate outputs
       out_cs_data(:,:,:), out_ps_data(:,:), out_ts_data(:),                 &
       out_mcrp_data(:,:,:,:),                                               &
       ! Removal rates for columns
       scs_rm(:,:,:)

  ! Case-dependent dimensions: mixing ratios and removal rates and binned mixing ratio profiles for active species
  INTEGER, SAVE :: nv1_mixr=0, nv1_rem=0, nv2_mixr=0, nv2_bin=0, &
    nv1_proc=0, nv2_proc=0, nv3_proc=0, out_mcrp_nout=0
  CHARACTER (len=7), SAVE, ALLOCATABLE :: s1_mixr(:), s1_rem(:), s2_mixr(:), &
    s2_aa(:), s2_ab(:), s2_ca(:), s2_cb(:), s2_p(:), s2_ia(:), s2_ib(:), s2_s(:), &
    s1_gas(:), s2_gas(:)

  ! SALSA cloud, precipitation, ice and snow masks
  LOGICAL, allocatable :: cloudmask(:,:,:), drizzmask(:,:,:), icemask(:,:,:), snowmask(:,:,:)

  public :: sflg, ssam_intvl, savg_intvl, statistics, init_stat, write_ps,   &
       acc_tend, updtst, sfc_stat, close_stat, fill_scalar, tke_sgs, sgsflxs,&
       sgs_vel, comp_tke, get_zi, acc_removal, cs_rem_set, csflg, &
       les_rate_stats, mcrp_var_save, out_cs_list, out_ps_list, out_ts_list, &
       cs_include, cs_exclude, ps_include, ps_exclude, ts_include, ts_exclude, &
       out_mcrp_data, out_mcrp_list, out_mcrp_nout

contains
  !
  ! ---------------------------------------------------------------------
  ! INIT_STAT:  This routine initializes the statistical arrays which
  ! are user/problem defined.  Note that svctr is given 100 elements, and
  ! elements 90 and above are used for computing the TKE budget. Hence
  ! if (nvar2 >= 90 the program stops
  !
  subroutine init_stat(time, filprf, expnme, nzp)

    use grid, only : nxp, nyp, nprc, nsnw, nspec, ngases, iradtyp, &
        sed_aero, sed_cloud, sed_precp, sed_ice, sed_snow, out_an_list, nv4_proc
    use mpi_interface, only : myid, ver, author, info
    use mo_submctl, only : fn1a,fn2a,fn2b,fnp2a,fnp2b,stat_b_bins, &
        nlcoag,nlcnd,nlauto,nlautosnow,nlactiv,nlicenucl,nlicmelt,ice_target_opt,nout_cld,nout_ice, &
        zspec, zgas

    character (len=80), intent (in) :: filprf, expnme
    integer, intent (in)            :: nzp
    real, intent (in)               :: time

    INTEGER :: i,ii,e,ee
    character (len=80) :: fname
    CHARACTER (len=7) :: tmp_list(maxn_list)='       '

    ! Local boolean arrays
    LOGICAL :: s1_bool(nvar1), s1_lvl4_bool(nv1_lvl4), s1_lvl5_bool(nv1_lvl5), s1_gas_bool(ngases)
    LOGICAL :: s2_bool(nvar2), s2_lvl4_bool(nv2_lvl4), s2_lvl5_bool(nv2_lvl5), s2_gas_bool(ngases)
    LOGICAL :: s3_bool(nvar3), s3_lvl4_bool(nv3_lvl4), s3_lvl5_bool(nv3_lvl5)
    LOGICAL :: s1_mixr_bool(4*(nspec+1)), s1_rem_bool(5*(nspec+1)), s2_mixr_bool(5*(nspec+1))
    LOGICAL :: s2_aa_bool(2+nspec), s2_ab_bool(2+nspec), s2_ca_bool(2+nspec), s2_cb_bool(2+nspec), &
        s2_p_bool(2+nspec), s2_ia_bool(2+nspec), s2_ib_bool(2+nspec), s2_s_bool(2+nspec)
    LOGICAL :: s2_CldHist_bool(nv2_hist), s2_IceHist_bool(nv2_hist)

    ! SALSA dimensions
    INTEGER, PARAMETER :: nv2_ndims = 6
    LOGICAL :: s2_dims_bool(nv2_ndims)=.FALSE.
    CHARACTER (len=7) :: s2_dims(nv2_ndims)=(/'P_Rd12a','P_Rd2ab','P_Rwprc','P_Rwsnw','P_hRc  ','P_hRi  '/)

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

    ! Default outputs
    ALLOCATE ( ssclr(nvar1), svctr(nzp,nvar2) )
    svctr(:,:) = 0.
    ssclr(:)   = 0.
    s1_bool(:) = .TRUE.
    s2_bool(:) = .TRUE.

    ! User selected process rate outputs: check, count and order inputs
    CALL test_user_vars(out_ts_list,maxn_list,nv1_proc)
    CALL test_user_vars(out_ps_list,maxn_list,nv2_proc)
    CALL test_user_vars(out_an_list,maxn_list,nv4_proc)
    IF (csflg) CALL test_user_vars(out_cs_list,maxn_list,nv3_proc)
    ! The old rate outputs as defaults
    IF (stat_micro .OR. stat_micro_ts .OR. stat_micro_ps) THEN
        ! The old defaults
        IF (level==3) THEN
            tmp_list(1:6)=(/'coag_rr','coag_nr','cond_rr','cond_nr','auto_rr','auto_nr'/)
            IF (sed_cloud) tmp_list(7)='sedi_rc'
            IF (sed_precp) tmp_list(8:9)=(/'sedi_rr','sedi_nr'/)
            tmp_list(10:11)=(/'diag_rr','diag_nr'/) ! 7-11
        ELSEIF (level>3) THEN
            IF (nlcoag) tmp_list(1:6)=(/'coag_ra','coag_na','coag_rc','coag_nc','coag_rr','coag_nr'/)
            IF (nlcnd) tmp_list(7:9)=(/'cond_ra','cond_rc','cond_rr'/)
            IF (nlauto) tmp_list(10:11)=(/'auto_rr','auto_nr'/)
            IF (nlactiv) tmp_list(12:13)=(/'cact_rc','cact_nc'/)
            IF (sed_aero) tmp_list(14:15)=(/'sedi_ra','sedi_na'/)
            IF (sed_cloud) tmp_list(16:17)=(/'sedi_rc','sedi_nc'/)
            IF (sed_precp) tmp_list(18:19)=(/'sedi_rr','sedi_nr'/)
            tmp_list(20:25)=(/'diag_ra','diag_na','diag_rc','diag_nc','diag_rr','diag_nr'/)
            IF (level>=5) THEN
                IF (nlcoag) tmp_list(26:29)=(/'coag_ri','coag_ni','coag_rs','coag_ns'/)
                IF (nlcnd) tmp_list(30:31)=(/'cond_ri','cond_rs'/)
                IF (nlautosnow .OR. (nlicenucl .AND. ice_target_opt>=0)) tmp_list(32:33)=(/'auto_rs','auto_ns'/)
                IF (nlicenucl) tmp_list(34:35)=(/'nucl_ri','nucl_ni'/)
                IF (sed_ice) tmp_list(36:37)=(/'sedi_ri','sedi_ni'/)
                IF (sed_snow) tmp_list(38:39)=(/'sedi_rs','sedi_ns'/)
                tmp_list(40:43)=(/'diag_ri','diag_ni','diag_rs','diag_ns'/)
                IF (nlicmelt) tmp_list(44:47)=(/'melt_ri','melt_ni','melt_rs','melt_ns'/)
            ENDIF
        ENDIF
        ! Set lists, but only when an empty list was given
        IF (stat_micro_ts .AND. nv1_proc==0) THEN
            out_ts_list = tmp_list
            CALL test_user_vars(out_ts_list,maxn_list,nv1_proc)
        ENDIF
        IF (stat_micro_ps .AND. nv2_proc==0) THEN
            out_ps_list= tmp_list
            CALL test_user_vars(out_ps_list,maxn_list,nv2_proc)
        ENDIF
        IF (stat_micro .AND. nv4_proc==0) THEN
            out_an_list= tmp_list
            CALL test_user_vars(out_an_list,maxn_list,nv4_proc)
        ENDIF
    ENDIF

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

    IF ( level < 4 ) THEN
       ! Merge logical and name arrays
       ! a) Time series
       i=nvar1+nv1_proc
       ALLOCATE( s1bool(i), s1total(i) )
       i=1; e=nvar1
       s1bool(i:e)=s1_bool; s1total(i:e)=s1
       IF (nv1_proc>0) THEN
          i=e+1; e=e+nv1_proc
          s1bool(i:e)=.TRUE.; s1total(i:e)=out_ts_list(1:nv1_proc)
       ENDIF
       ! b) Profiles
       i=nvar2+nv2_proc
       ALLOCATE( s2bool(i), s2total(i) )
       i=1; e=nvar2
       s2bool(i:e)=s2_bool; s2total(i:e)=s2
       IF (nv2_proc>0) THEN
          i=e+1; e=e+nv2_proc
          s2bool(i:e)=.TRUE.; s2total(i:e)=out_ps_list(1:nv2_proc)
       ENDIF
    ELSE IF ( level >= 4 ) THEN
       ! Additional arrays for SALSA
       ! -- dimensions
       nv1_mixr = 4*(nspec+1) ! Mixing ratios; in-cloud, intersitial, in-ice and in-snow
       nv1_rem = 5*(nspec+1)  ! Removal with aerosol, cloud, rain, ice and snow
       nv2_mixr = 5*(nspec+1) ! Mixing ratio profiles for aerosol, cloud, rain, ice and snow
       nv2_bin  = 2+nspec     ! Bin number concentrations and species mixing ratios
       ! -- allocate
       ALLOCATE ( ssclr_lvl4(nv1_lvl4), svctr_lvl4(nzp,nv2_lvl4) )
       ALLOCATE ( ssclr_lvl5(nv1_lvl5), svctr_lvl5(nzp,nv2_lvl5) )
       ALLOCATE ( ssclr_mixr(nv1_mixr), ssclr_rem(nv1_rem), svctr_mixr(nzp,nv2_mixr) )
       ALLOCATE ( s1_mixr(nv1_mixr),s1_rem(nv1_rem), s2_mixr(nv2_mixr) )
       ALLOCATE ( svctr_aa(nzp,fn2a,nv2_bin), svctr_ab(nzp,fn2b-fn2a,nv2_bin), &
                  svctr_ca(nzp,fnp2a,nv2_bin), svctr_cb(nzp,fnp2b-fnp2a,nv2_bin), &
                  svctr_p(nzp,nprc,nv2_bin), svctr_ia(nzp,fnp2a,nv2_bin), &
                  svctr_ib(nzp,fnp2b-fnp2a,nv2_bin), svctr_s(nzp,nsnw,nv2_bin) )
       ALLOCATE ( s2_aa(nv2_bin), s2_ab(nv2_bin), s2_ca(nv2_bin), s2_cb(nv2_bin), &
                  s2_p(nv2_bin), s2_ia(nv2_bin), s2_ib(nv2_bin), s2_s(nv2_bin) )
       ALLOCATE ( cloudmask(nzp,nxp,nyp), drizzmask(nzp,nxp,nyp), &
                  icemask(nzp,nxp,nyp), snowmask(nzp,nxp,nyp) )
       ALLOCATE ( ssclr_gas(ngases), svctr_gas(nzp,ngases) )
       ALLOCATE ( s1_gas(ngases), s2_gas(ngases) )
       ! -- reset data arrays
       ssclr_lvl4(:) = 0.; svctr_lvl4(:,:) = 0.
       ssclr_lvl5(:) = 0.; svctr_lvl5(:,:) = 0.
       ssclr_mixr(:) = 0.; svctr_mixr(:,:) = 0.
       ssclr_rem(:) = 0.
       svctr_aa(:,:,:) = 0.; svctr_ab(:,:,:) = 0.
       svctr_ca(:,:,:) = 0.; svctr_cb(:,:,:) = 0.
       svctr_p(:,:,:) = 0.
       svctr_ia(:,:,:) = 0.; svctr_ib(:,:,:) = 0.
       svctr_s(:,:,:) = 0.
       ssclr_gas(:)=0.; svctr_gas(:,:) = 0.

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
       s2_lvl4_bool(:) = .TRUE.
       s2_lvl5_bool(:) = (level>4)

       ! b-bins are not always saved
       IF (.not. stat_b_bins) THEN
          ! Concentrations and sizes
          s1_lvl4_bool(3:6) = .FALSE.
          s1_lvl4_bool(9:12) = .FALSE.
          s1_lvl5_bool(3:6) = .FALSE.
          s2_lvl4_bool((/2,4,7,9/)) = .FALSE.
          s2_lvl5_bool((/2,5/)) = .FALSE.
       ENDIF

       ! 2) Microphysical process rate statistics (both ts and ps)

       ! 3) Species or bin dependent outputs
       s1_mixr_bool(:)=.TRUE.
       s1_rem_bool(:)=lremts
       s2_mixr_bool(:)=lmixrprof

       s2_aa_bool(:)=lbinprof
       s2_ab_bool(:)=lbinprof .AND. stat_b_bins
       s2_ca_bool(:)=lbinprof
       s2_cb_bool(:)=lbinprof .AND. stat_b_bins
       s2_p_bool(:)=lbinprof
       s2_ia_bool(:)=lbinprof .AND. (level>4)
       s2_ib_bool(:)=lbinprof .AND. (level>4) .AND. stat_b_bins
       s2_s_bool(:)=lbinprof .AND. (level>4)

       DO ee=1,nspec+1 ! Aerosol species and water
          ! Mixing ratios
          ii = (ee-1)*4
          s1_mixr(ii+1)=TRIM(zspec(ee))//'_ic'
          s1_mixr(ii+2)=TRIM(zspec(ee))//'_int'
          s1_mixr(ii+3)=TRIM(zspec(ee))//'_ii'
          s1_mixr(ii+4)=TRIM(zspec(ee))//'_is'
          s1_mixr_bool(ii+3:ii+4) = .TRUE. .AND. (level>4)

          ! Removal; dry, cloud, precipitation and snow
          ii = (ee-1)*5
          s1_rem(ii+1)='rm'//TRIM(zspec(ee))//'dr'
          s1_rem(ii+2)='rm'//TRIM(zspec(ee))//'cl'
          s1_rem(ii+3)='rm'//TRIM(zspec(ee))//'pr'
          s1_rem(ii+4)='rm'//TRIM(zspec(ee))//'ic'
          s1_rem(ii+5)='rm'//TRIM(zspec(ee))//'sn'
          s1_rem_bool(ii+1) = lremts .AND. sed_aero
          s1_rem_bool(ii+2) = lremts .AND. sed_cloud
          s1_rem_bool(ii+3) = lremts .AND. sed_precp
          s1_rem_bool(ii+4) = lremts .AND. sed_ice .AND. (level>4)
          s1_rem_bool(ii+5) = lremts .AND. sed_snow .AND. (level>4)

          ! Total mixing ratios, profiles
          ii = (ee-1)*5
          s2_mixr(ii+1)='P_c'//TRIM(zspec(ee))//'a'
          s2_mixr(ii+2)='P_c'//TRIM(zspec(ee))//'c'
          s2_mixr(ii+3)='P_c'//TRIM(zspec(ee))//'p'
          s2_mixr(ii+4)='P_c'//TRIM(zspec(ee))//'i'
          s2_mixr(ii+5)='P_c'//TRIM(zspec(ee))//'s'
          s2_mixr_bool(ii+4:ii+5) = lmixrprof .AND. (level>4)

          ! Mixing ratios for each bin (profiles)
          IF (ee==1) THEN
             ! Number concentrations first
             ii=1
             s2_aa(ii)='P_Naba'
             s2_ab(ii)='P_Nabb'
             s2_ca(ii)='P_Ncba'
             s2_cb(ii)='P_Ncbb'
             s2_p(ii)='P_Npb'
             s2_ia(ii)='P_Niba'
             s2_ib(ii)='P_Nibb'
             s2_s(ii)='P_Nsb'
          ENDIF
          ! Add species
          ii=1+ee
          s2_aa(ii)='P_'//TRIM(zspec(ee))//'aa'
          s2_ab(ii)='P_'//TRIM(zspec(ee))//'ab'
          s2_ca(ii)='P_'//TRIM(zspec(ee))//'ca'
          s2_cb(ii)='P_'//TRIM(zspec(ee))//'cb'
          s2_p(ii)='P_'//TRIM(zspec(ee))//'pb'
          s2_ia(ii)='P_'//TRIM(zspec(ee))//'ia'
          s2_ib(ii)='P_'//TRIM(zspec(ee))//'ib'
          s2_s(ii)='P_'//TRIM(zspec(ee))//'sb'
       ENDDO

       ! SALSA dimensions
       !    s2_dims(nv2_ndims)=(/'P_Rd12a','P_Rd2ab','P_Rwprc','P_Rwsnw','P_hRc  ','P_hRi  '/)
       IF (lbinprof) THEN
            s2_dims_bool(1:3)=.TRUE.
            s2_dims_bool(4)=(level>4)
       ENDIF

       ! Cloud and ice histrograms
       s2_CldHist_bool=.FALSE.
       IF (nout_cld>0) THEN
          s2_CldHist_bool(1) = .TRUE. ! A-bins
          s2_CldHist_bool(2) = stat_b_bins ! B-bins
          ! Add dimension
          s2_dims_bool(5)=.TRUE.
       END IF
       s2_IceHist_bool=.FALSE.
       IF (nout_ice>0 .AND. level>=5) THEN
          s2_IceHist_bool(1) = .TRUE.
          s2_IceHist_bool(2) = stat_b_bins
          ! Add dimension
          s2_dims_bool(6)=.TRUE.
       END IF

       ! Gas phase concentrations
       DO i=1,ngases
          s1_gas_bool(i) = .FALSE.
          s1_gas(i) = 'int'//TRIM(zgas(i))//'g'
          s2_gas_bool(i) = .TRUE.
          s2_gas(i) = 'P_c'//TRIM(zgas(i))//'g'
       ENDDO

        ! Merge logical and name arrays
        ! a) Time series
        i=nvar1+nv1_lvl4+nv1_lvl5+nv1_mixr+nv1_rem+ngases+nv1_proc
        ALLOCATE( s1bool(i), s1total(i) )
        i=1; e=nvar1
        s1bool(i:e)=s1_bool; s1total(i:e)=s1
        i=e+1; e=e+nv1_lvl4
        s1bool(i:e)=s1_lvl4_bool; s1total(i:e)=s1_lvl4
        i=e+1; e=e+nv1_lvl5
        s1bool(i:e)=s1_lvl5_bool; s1total(i:e)=s1_lvl5
        i=e+1; e=e+nv1_mixr
        s1bool(i:e)=s1_mixr_bool; s1total(i:e)=s1_mixr
        i=e+1; e=e+nv1_rem
        s1bool(i:e)=s1_rem_bool; s1total(i:e)=s1_rem
        IF (ngases>0) THEN
            i=e+1; e=e+ngases
            s1bool(i:e)=s1_gas_bool; s1total(i:e)=s1_gas
        ENDIF
        IF (nv1_proc>0) THEN
            i=e+1; e=e+nv1_proc
            s1bool(i:e)=.TRUE.; s1total(i:e)=out_ts_list(1:nv1_proc)
        ENDIF
        ! b) Profiles
        i=nvar2+nv2_lvl4+nv2_lvl5+nv2_mixr+8*nv2_bin+2*nv2_hist+nv2_ndims+ngases+nv2_proc
        ALLOCATE( s2bool(i), s2total(i) )
        i=1; e=nvar2
        s2bool(i:e)=s2_bool; s2total(i:e)=s2
        i=e+1; e=e+nv2_lvl4
        s2bool(i:e)=s2_lvl4_bool; s2total(i:e)=s2_lvl4
        i=e+1; e=e+nv2_lvl5
        s2bool(i:e)=s2_lvl5_bool; s2total(i:e)=s2_lvl5
        i=e+1; e=e+nv2_mixr
        s2bool(i:e)=s2_mixr_bool; s2total(i:e)=s2_mixr
        i=e+1; e=e+nv2_bin
        s2bool(i:e)=s2_aa_bool; s2total(i:e)=s2_aa
        i=e+1; e=e+nv2_bin
        s2bool(i:e)=s2_ab_bool; s2total(i:e)=s2_ab
        i=e+1; e=e+nv2_bin
        s2bool(i:e)=s2_ca_bool; s2total(i:e)=s2_ca
        i=e+1; e=e+nv2_bin
        s2bool(i:e)=s2_cb_bool; s2total(i:e)=s2_cb
        i=e+1; e=e+nv2_bin
        s2bool(i:e)=s2_p_bool; s2total(i:e)=s2_p
        i=e+1; e=e+nv2_bin
        s2bool(i:e)=s2_ia_bool; s2total(i:e)=s2_ia
        i=e+1; e=e+nv2_bin
        s2bool(i:e)=s2_ib_bool; s2total(i:e)=s2_ib
        i=e+1; e=e+nv2_bin
        s2bool(i:e)=s2_s_bool; s2total(i:e)=s2_s
        i=e+1; e=e+nv2_hist
        s2bool(i:e)=s2_CldHist_bool; s2total(i:e)=s2_CldHist
        i=e+1; e=e+nv2_hist
        s2bool(i:e)=s2_IceHist_bool; s2total(i:e)=s2_IceHist
        i=e+1; e=e+nv2_ndims
        s2bool(i:e)=s2_dims_bool; s2total(i:e)=s2_dims
        IF (ngases>0) THEN
            i=e+1; e=e+ngases
            s2bool(i:e)=s2_gas_bool; s2total(i:e)=s2_gas
        ENDIF
        IF (nv2_proc>0) THEN
            i=e+1; e=e+nv2_proc
            s2bool(i:e)=.TRUE.; s2total(i:e)=out_ps_list(1:nv2_proc)
        ENDIF
    END IF ! If level >=4

    ! User options for including and excluding outputs
    CALL apply_user_outputs(.FALSE.,ps_exclude,maxn_list,s2total,s2bool,SIZE(s2bool))
    CALL apply_user_outputs(.TRUE.,ps_include,maxn_list,s2total,s2bool,SIZE(s2bool))
    CALL apply_user_outputs(.FALSE.,ts_exclude,maxn_list,s1total,s1bool,SIZE(s1bool))
    CALL apply_user_outputs(.TRUE.,ts_include,maxn_list,s1total,s1bool,SIZE(s1bool))

    fname =  trim(filprf)//'.ts'
    if(myid == 0) print                                                  &
         "(//' ',49('-')/,' ',/,'  Initializing: ',A20,'  N=',I3)",trim(fname),COUNT(s1bool)
    call open_nc( fname, expnme, time, (nxp-4)*(nyp-4), ncid1, nrec1, ver, author, info)
    ! Juha: Modified for SALSA output
    call define_nc( ncid1, nrec1, COUNT(s1bool), PACK(s1Total,s1bool))
    if (myid == 0) print *, '   ...starting record: ', nrec1


    fname =  trim(filprf)//'.ps'
    if(myid == 0) print                                                  &
         "(//' ',49('-')/,' ',/,'  Initializing: ',A20,'  N=',I3)",trim(fname),COUNT(s2bool)
    call open_nc( fname, expnme, time,(nxp-4)*(nyp-4), ncid2, nrec2, ver, author, info)
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
    if (myid == 0) print *, '   ...starting record: ', nrec2


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

        ! Merge logical and name arrays (reuse s1 arrays)
        DEALLOCATE( s1bool, s1total )
        ! Default (mostly levels 3 & 4) parameters, specific level 4 and 5 parameters,
        ! removal statistics, gases, and process rate statistics.
        ! Note: removal and gas outputs are the same as those in the time series outputs!
        i = nvar3 + nv3_lvl4 + nv3_lvl5 + nv1_rem+ngases+nv3_proc
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
        IF (nv1_rem>0) THEN
            i=e+1; e=e+nv1_rem
            ! Logical and name arrays are the same as those for ts
            s1bool(i:e)=s1_rem_bool; s1total(i:e)=s1_rem
        ENDIF
        IF (ngases>0) THEN
            i=e+1; e=e+ngases
            ! Name array is the same as that for ts, but logical is not (set to true)
            s1bool(i:e)=.TRUE.; s1total(i:e)=s1_gas
        ENDIF
        IF (nv3_proc>0) THEN
            i=e+1; e=e+nv3_proc
            s1bool(i:e)=.TRUE.; s1total(i:e)=out_cs_list(1:nv3_proc)
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
    CHARACTER(LEN=7), DIMENSION(*), INTENT(IN) :: out_list ! Current output list
    LOGICAL, DIMENSION(*) ,INTENT(INOUT) :: out_flags ! Current output flags
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
  ! Modified for level 4: rxt contains a_rp if level < 4 and a_rp+a_rc
  ! if level == 4
  ! Juha Tonttila, FMI, 2014
  !
  ! Modified for level 5
  ! Jaakko Ahola, FMI, 2016
  subroutine statistics(time)

    use grid, only : a_up, a_vp, a_wp, a_rc, a_theta, a_temp, a_rv           &
         , a_rp, a_tp, a_press, nxp, nyp, nzp, dzm, dzt, zm, zt, th00, umean            &
         , vmean, dn0, precip, a_rpp, a_npp, CCN, iradtyp, a_rflx, a_sflx               &
         , a_fus, a_fds, a_fuir, a_fdir, albedo, a_srp, a_snrp, a_ncloudp, a_ri, a_nicep, a_srs, a_snrs
    USE defs, ONLY : cp, alvi

    real, intent (in) :: time

    real :: rxt(nzp,nxp,nyp), rxl(nzp,nxp,nyp), rxv(nzp,nxp,nyp), rnt(nzp,nxp,nyp)
    REAL :: xrpp(nzp,nxp,nyp), xnpp(nzp,nxp,nyp), thl(nzp,nxp,nyp)

    SELECT CASE(level)
       CASE (3)
          rxt = a_rp ! Total water (vapor + condensed water and ice) = q
          rxl = a_rc ! Cloud water (+aerosol), but no precipitation or ice
          rxv = a_rv ! Water vapor
          xrpp = a_rpp ! Rain water
          xnpp = a_npp ! Rain number
          thl = a_tp ! Liquid water potential temperature
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
          WHERE(a_temp>0) thl = a_tp + (a_theta/a_temp)*alvi/cp*(a_ri + a_srs)
       CASE DEFAULT
          rxt = a_rp
          rxl = a_rc
          rxv = a_rv
          xrpp = 0.
          xnpp = 0.
          thl = a_tp
    END SELECT

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
    if (level >=1) call accum_lvl1(nzp, nxp, nyp, rxt)
    if (level >=2) call accum_lvl2(nzp, nxp, nyp, th00, a_wp,        &
                                   a_theta, thl, rxv, rxl, rxt   )
    if (level >=3) call accum_lvl3(nzp, nxp, nyp, dn0, zm, rxl, xrpp,  &
                                   xnpp, precip, CCN                    )
    if (level >=4)  call accum_lvl4(nzp, nxp, nyp)
    if (level >=5)  call accum_lvl5(nzp, nxp, nyp)
    !
    ! scalar statistics
    !
    call set_ts(nzp, nxp, nyp, a_wp, a_theta, thl, dn0, zt,zm,dzt,dzm,th00,time)
    IF ( level >=1 ) CALL ts_lvl1(nzp, nxp, nyp, dn0, zt, dzm, rxt)
    IF ( level >=2 ) CALL ts_lvl2(nzp, nxp, nyp, dn0, zm, zt, rxl, rxv)
    IF ( level >=4 ) CALL ts_lvl4(nzp, nxp, nyp)
    IF ( level >=5 ) CALL ts_lvl5(nzp, nxp, nyp)

    call write_ts

    !
    ! Column statistics
    !
    IF (csflg) THEN
        ! Warm cloud statistics (xrpp and xnpp from above; CDNC is calculated below)
        rnt = CCN
        IF (level>3) rnt = SUM(a_ncloudp,DIM=4)
        CALL set_cs_warm(nzp,nxp,nyp,a_rc,rnt,xrpp,xnpp,a_theta,a_wp,precip,dn0,zm,zt,dzm)

        ! Ice cloud statistics
        IF (level==5) THEN
            rnt = SUM(a_nicep,DIM=4)
            CALL set_cs_cold(nzp,nxp,nyp,a_ri,rnt,a_srs,a_snrs,dn0,zm)
        ENDIF

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
  subroutine set_cs_warm(n1,n2,n3,rc,nc,rp,np,th,w,rrate,dn0,zm,zt,dzm)

    use netcdf

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: rc(n1,n2,n3),nc(n1,n2,n3),rp(n1,n2,n3),np(n1,n2,n3),th(n1,n2,n3)
    real, intent(in)    :: w(n1,n2,n3),rrate(n1,n2,n3),dn0(n1),zm(n1),zt(n1),dzm(n1)
    REAL :: lwp(n2,n3), ncld(n2,n3), rwp(n2,n3), nrain(n2,n3), zb(n2,n3), zc(n2,n3), &
                th1(n2,n3), lmax(n2,n3), prcp_bc(n2,n3), wbase(n2,n3), wmax(n2,n3)
    INTEGER :: ncloudy(n2,n3), nrainy(n2,n3)
    integer :: i, j, k, iret, VarID
    real    :: cld, rn, sval, dmy

    ! No outputs for level 1
    IF (level<2) RETURN

    ! Calculate stats
    lwp=0.      ! LWP (kg/m^2)
    ncld=0.     ! Average CDNC (#/kg)
    ncloudy=0 ! Number of cloudy grid cells
    rwp=0.      ! RWP (kg/m^2)
    nrain=0.    ! Average RDNC (#/kg)
    nrainy=0    ! Number of cloudy grid cells
    zb=zm(n1)+100.  ! Cloud base (m)
    zc=0.           ! Cloud top (m)
    lmax=0.     ! Liquid water mixing ratio (kg/kg)
    th1=0.      ! Height of the maximum theta gradient
    prcp_bc = -999. ! Precipitation at the cloud base (W/m^2)
    wbase = -999.   ! Vertical velocity -||- (m/s)
    wmax = 0.   ! Maximum absolute vertical velocity (m/s)
    do j=3,n3-2
       do i=3,n2-2
          cld=0.
          rn=0.
          sval = 0.
          do k=2,n1
             lwp(i,j)=lwp(i,j)+rc(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
             IF (rc(k,i,j)>0.01e-3) THEN
                ! Cloudy grid
                ! Volume weighted average of the CDNC
                ncld(i,j)=ncld(i,j)+nc(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                cld=cld+dn0(k)*(zm(k)-zm(k-1))
                ! Number of cloudy pixels
                ncloudy(i,j)=ncloudy(i,j)+1
                ! Cloud base and top
                zb(i,j)=min(zt(k),zb(i,j))
                zc(i,j)=max(zt(k),zc(i,j))
                ! Cloud base vertical velocity and precipitation
                IF (ncloudy(i,j)==1) THEN
                    ! Take those from level k-1 (>=2), which is just below cloud base
                    wbase(i,j) = w(max(2,k-1),i,j)
                    prcp_bc(i,j) = rrate(max(2,k-1),i,j)
                ENDIF
             END IF
             !
             rwp(i,j)=rwp(i,j)+rp(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
             if (rp(k,i,j) > 0.001e-3) then
                ! Rainy grid cell
                ! Volume weighted average of the RDNC
                nrain(i,j)=nrain(i,j)+np(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                rn=rn+dn0(k)*(zm(k)-zm(k-1))
                ! Number of rainy pixels
                nrainy(i,j)=nrainy(i,j)+1
             end if
             ! Maximum liquid water mixing ratio
             lmax(i,j) = max(lmax(i,j),rc(k,i,j))
             ! Maximum absolute vertical velocity
             IF (abs(wmax(i,j))<ABS(w(k,i,j))) wmax(i,j) = w(k,i,j)
             ! Height of the maximum theta gradient
             if (k<=n1-5) then
                dmy = (th(k+1,i,j)-th(k,i,j))*dzm(k)
                if (dmy > sval ) then
                   sval = dmy
                   th1(i,j) = zt(k)
                end if
            ENDIF
          enddo
          IF (cld>0.) THEN
            ncld(i,j)=ncld(i,j)/cld
          ELSE
            zb(i,j)=-999.
            zc(i,j)=-999.
          END IF
          IF (rn>0.) THEN
            nrain(i,j)=nrain(i,j)/rn
          END IF
       end do
    end do

    ! Save the data

    iret = nf90_inq_varid(ncid3,'lwp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, lwp(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'rwp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, rwp(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'Nc',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, ncld(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'Nr',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nrain(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'nccnt',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, ncloudy(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'nrcnt',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nrainy(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'zb',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, zb(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'zc',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, zc(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'zi1',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, th1(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'lmax',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, lmax(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'wmax',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, wmax(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'wbase',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, wbase(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'prcp_bc',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, prcp_bc(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

  end subroutine set_cs_warm

  ! Calculate cold cloud statistics
  subroutine set_cs_cold(n1,n2,n3,ri,ni,rs,ns,dn0,zm)

    use netcdf

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: ri(n1,n2,n3),ni(n1,n2,n3),rs(n1,n2,n3),ns(n1,n2,n3)
    real, intent(in)    :: dn0(n1),zm(n1)
    REAL :: iwp(n2,n3), nice(n2,n3), swp(n2,n3), nsnow(n2,n3), imax(n2,n3)
    INTEGER :: nicy(n2,n3), nsnowy(n2,n3)
    integer :: i, j, k, iret, VarID
    real    :: ice, sn

    ! No outputs for levels less than 5
    IF (level<5) RETURN

    ! Calculate stats
    iwp=0.  ! IWP (kg/m^2)
    nice=0. ! Average ice number concentration (#/kg)
    nicy=0  ! Number of icy grid cells
    swp=0.  ! SWP (kg/m^2)
    nsnow=0.    ! Average snow number concentration (#/kg)
    nsnowy=0    ! Number of snowy grid cells
    imax=0. ! Maximum ice water mixing ratio (kg/kg)
    do j=3,n3-2
       do i=3,n2-2
          ice=0.
          sn=0.
          do k=2,n1
             iwp(i,j)=iwp(i,j)+ri(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
             IF (icemask(k,i,j)) THEN
                ! Icy grid cell
                ! Volume weighted average of the ice number concentration
                nice(i,j)=nice(i,j)+ni(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                ice=ice+dn0(k)*(zm(k)-zm(k-1))
                ! Number of icy pixels
                nicy(i,j)=nicy(i,j)+1
             END IF
             swp(i,j)=swp(i,j)+rs(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
             if (snowmask(k,i,j)) then
                ! Snowy grid cell
                ! Volume weighted average of the snow number concentration
                nsnow(i,j)=nsnow(i,j)+ns(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                sn=sn+dn0(k)*(zm(k)-zm(k-1))
                ! Number of snowy pixels
                nsnowy(i,j)=nsnowy(i,j)+1
             end if
             ! Maximum ice water mixing ratio
             imax(i,j) = max(imax(i,j),ri(k,i,j))
          enddo
          IF (ice>0.) THEN
            nice(i,j)=nice(i,j)/ice
          END IF
          IF (sn>0.) THEN
            nsnow(i,j)=nsnow(i,j)/sn
          END IF
       end do
    end do

    iret = nf90_inq_varid(ncid3,'iwp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, iwp(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'swp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, swp(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'Ni',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nice(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'Ns',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nsnow(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'nicnt',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nicy(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'nscnt',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nsnowy(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'imax',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, imax(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

  end subroutine set_cs_cold

  ! Other column statistics
  subroutine set_cs_other(time)
    use netcdf
    USE grid, ONLY : nzp, nxp, nyp, xt, xm, yt, ym, a_dn, zm, albedo, ngases, a_gaerop, precip
    ! Inputs
    REAL, INTENT(IN) :: time
    ! Local variables
    REAL :: gwp(nxp,nyp)
    integer :: i, k, iret, VarID

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

    ! Surface precipitation
    iret = nf90_inq_varid(ncid3,'prcp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, precip(2,3:nxp-2,3:nyp-2), start=(/1,1,nrec3/))

    ! Level 4 and 5 removal/precipitation statistics (variable names are the same as those for ts)
    DO i = 1,nv1_rem
        iret = nf90_inq_varid(ncid3,s1_rem(i),VarID)
        IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, scs_rm(3:nxp-2,3:nyp-2,i), start=(/1,1,nrec3/))
    ENDDO

    ! Gases (names are the same as those for ts)
    DO i=1,ngases
        iret = nf90_inq_varid(ncid3,s1_gas(i),VarID)
        IF (iret==NF90_NOERR) THEN
            ! Calculate gas path (kg/m2)
            gwp(:,:) = 0.
            do k=2,nzp
                gwp(:,:) = gwp(:,:) + a_gaerop(k,:,:,i)*a_dn(k,:,:)*(zm(k)-zm(k-1))
            ENDDO
            ! Save the data
            iret = nf90_put_var(ncid3, VarID, gwp(3:nxp-2,3:nyp-2), start=(/1,1,nrec3/))
        ENDIF
    END DO

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
  ! -----------------------------------------------------------------------
  ! subroutine set_ts: computes and writes time sequence stats
  !
  subroutine set_ts(n1,n2,n3,w,th,t,dn0,zt,zm,dzt,dzm,th00,time)

    use defs, only : cp
    USE grid, ONLY : th0

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: w(n1,n2,n3),th(n1,n2,n3),t(n1,n2,n3)
    real, intent(in)    :: dn0(n1),zt(n1),zm(n1),dzt(n1),dzm(n1),th00,time

    integer :: i, j, k
    real    :: bf(n1), scr(n2,n3)

    ssclr(1) = time
    ssclr(4) = get_zi(n1, n2, n3, 2, th, dzm, zt, 1.)   ! maximum gradient
    ssclr(5) = get_zi(n1, n2, n3, 3, th, thvar, zt, 1.) ! maximum variance
    !
    ! buoyancy flux statistics
    !
    ssclr(7) = 0.
    ssclr(32) = 0.
    do k = 2,n1-2
       bf(k) = wtv_res(k) + wtv_sgs(k)
       ssclr(7) = ssclr(7) + (tke_res(k)+tke_sgs(k))*dn0(k)/dzt(k)
       ssclr(32) = ssclr(32) + (tke_res(k)+tke_sgs(k))/dzt(k)
       svctr(k,33) = svctr(k,33) + wtv_sgs(k)*9.8/th00
    end do
    ssclr(6) = get_zi(n1, n2, n3, 4, th, bf, zm, 1.) ! minimum buoyancy flux

    ssclr(8) = bf(2)
    ssclr(9) = maxval(w)

    ssclr(12) = ssclr(12)*cp*(dn0(1)+dn0(2))*0.5

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
  ! subroutine ts_lvl1: computes and writes time sequence stats; for the
  ! zi calculation setting itype=1 selects a concentration threshold
  !
  subroutine ts_lvl1(n1,n2,n3,dn0,zt,dzm,q)

    use defs, only : alvl

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: q(n1,n2,n3)
    real, intent(in)    :: dn0(n1),zt(n1),dzm(n1)

    ssclr(13) = ssclr(13)*alvl*(dn0(1)+dn0(2))*0.5
    ssclr(14) = get_zi(n1, n2, n3, 2, q, dzm, zt, 0.5e-3)

  end subroutine ts_lvl1
  !
  ! -----------------------------------------------------------------------
  ! subroutine ts_lvl2: computes and writes time sequence stats
  !
  subroutine ts_lvl2(n1,n2,n3,dn0,zm,zt,rc,rv)

    integer, intent(in) :: n1,n2,n3
    real, intent(in) :: dn0(n1), zm(n1), zt(n1), rc(n1,n2,n3), rv(n1,n2,n3)

    integer :: k,i,j
    real    :: cpnt, unit, scr(n2,n3), scr2(n2,n3), ct_tmp, cb_tmp

    ssclr(18)  = zt(n1)
    ssclr(19)  = 0.
    ssclr(20)  = 0.
    ssclr(28)  = 0.
    ssclr(30)  = 0.
    ssclr(31)  = 0.

    unit = 1./real((n2-4)*(n3-4))
    scr(:,:) = 0.   ! LWP
    scr2(:,:) = 0.  ! WVP (water vapor path)
    do j=3,n3-2
       do i=3,n2-2
          cpnt  = 0.
          ct_tmp =0.
          cb_tmp =zt(n1)
          do k=2,n1
             if (rc(k,i,j) > 1.e-5) then
                ssclr(17) = max(ssclr(17),zt(k))
                ssclr(18) = min(ssclr(18),zt(k))
                ct_tmp = max(ct_tmp,zt(k))
                cb_tmp = min(cb_tmp,zt(k))
                cpnt = unit
                ssclr(20) = max(ssclr(20), rc(k,i,j))
                ssclr(28) = ssclr(28) + 1.
             end if
             scr(i,j)=scr(i,j)+rc(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
             scr2(i,j)=scr2(i,j)+rv(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
          end do
          ssclr(19) = ssclr(19) + cpnt
          if (cpnt>0.) then
             ssclr(30) = ssclr(30)+ct_tmp
             ssclr(31) = ssclr(31)+cb_tmp
          end if
       end do
    end do
    if (ssclr(19)>0.) then
       ssclr(30) =ssclr(30)/ssclr(19)*unit
       ssclr(31) =ssclr(31)/ssclr(19)*unit
    else
       ssclr(30) =-999.
       ssclr(31) =-999.
    end if

    ! liquid water path (without precipitation)
    ssclr(15) = get_avg2dh(n2,n3,scr)
    scr(:,:)=(scr(:,:)-ssclr(15))**2 ! For LWP variance
    ssclr(16) = get_avg2dh(n2,n3,scr)

    ! water vapor path
    ssclr(34) = get_avg2dh(n2,n3,scr2)
    scr2(:,:)=(scr2(:,:)-ssclr(34))**2 ! For WVP variance
    ssclr(35) = get_avg2dh(n2,n3,scr2)

    if (ssclr(18) == zt(n1)) ssclr(18) = -999.

  end subroutine ts_lvl2
  !
  ! -----------------------------------------------------------------------
  ! subroutine ts_lvl4: computes and writes time sequence stats of Salsa variables --
  !  Implemented by Zubair Maalick 20/07/2015
  !  Some rewriting and adjusting by Juha Tonttila
  !
  SUBROUTINE ts_lvl4(n1,n2,n3)
    use mo_submctl, only : nlim ! Note: #/m^3, but close enough to #/kg for statistics
    USE grid, ONLY : bulkNumc, bulkMixrat, meanradius, dzt, a_rh, &
        a_gaerop, nspec, ngases

    IMPLICIT NONE

    integer, intent(in) :: n1,n2,n3

    REAL :: a0(n1,n2,n3), a1(n1,n2,n3), a2(n1,n2,n3)
    integer :: ii, i
    LOGICAL :: mask(n1,n2,n3)

    ! Clouds; combined and a and b separately for a and b regions, respectively
    CALL bulkNumc('cloud','a',a0)
    CALL bulkNumc('cloud','b',a1)
    ssclr_lvl4(1) = get_avg_ts(n1,n2,n3,a0+a1,dzt,cloudmask)
    CALL meanRadius('cloud','ab',a2)
    ssclr_lvl4(2) = get_avg_ts(n1,n2,n3,a2,dzt,cloudmask)

    mask(:,:,:) = ( a0(:,:,:) > nlim .AND. cloudmask(:,:,:) )
    ssclr_lvl4(3) = get_avg_ts(n1,n2,n3,a0,dzt,mask)
    CALL meanRadius('cloud','a',a2)
    ssclr_lvl4(4) = get_avg_ts(n1,n2,n3,a2,dzt,mask)

    mask(:,:,:) = ( a1(:,:,:) > nlim .AND. cloudmask(:,:,:) )
    ssclr_lvl4(5) = get_avg_ts(n1,n2,n3,a1,dzt,mask)
    CALL meanRadius('cloud','b',a2)
    ssclr_lvl4(6) = get_avg_ts(n1,n2,n3,a2,dzt,mask)

    ! Interstitial aerosol; combined and a and b separately
    CALL bulkNumc('aerosol','a',a0)
    CALL bulkNumc('aerosol','b',a1)
    ssclr_lvl4(7) = get_avg_ts(n1,n2,n3,a0+a1,dzt,cloudmask)
    CALL meanRadius('aerosol','ab',a2)
    ssclr_lvl4(8) = get_avg_ts(n1,n2,n3,a2,dzt,cloudmask)

    ssclr_lvl4(9) = get_avg_ts(n1,n2,n3,a0,dzt,cloudmask)
    CALL meanRadius('aerosol','a',a2)
    ssclr_lvl4(10) = get_avg_ts(n1,n2,n3,a2,dzt,cloudmask)

    ssclr_lvl4(11) = get_avg_ts(n1,n2,n3,a1,dzt,cloudmask)
    CALL meanRadius('aerosol','b',a2)
    ssclr_lvl4(12) = get_avg_ts(n1,n2,n3,a2,dzt,cloudmask)

    ! Maximum supersaturation
    ssclr_lvl4(13) = (MAXVAL(a_rh(2:n1,3:n2-2,3:n3-2))-1.0)*100.

    ! Species mixing ratios in cloud droplets and interstitial aerosol
    DO i=1,nspec+1 ! Aerosol species and water
        ii = (i-1)*4 +1
        CALL bulkMixrat(i,'cloud','ab',a0)
        ssclr_mixr(ii) = get_avg_ts(n1,n2,n3,a0,dzt,cloudmask)
        CALL bulkMixrat(i,'aerosol','ab',a0)
        ssclr_mixr(ii+1) = get_avg_ts(n1,n2,n3,a0,dzt,cloudmask)
    END DO

    ! Gases
    DO i=1,ngases
        a0(:,:,:) = a_gaerop(:,:,:,i)
        ssclr_gas(i) = get_avg_ts(n1,n2,n3,a0,dzt)
    END DO

  END SUBROUTINE ts_lvl4
  !
  ! -----------------------------------------------------------------------
  ! subroutine ts_lvl5: computes and writes time sequence stats of Salsa variables --
  !  Implemented by Jaakko Ahola 15/12/2016
  !
  SUBROUTINE ts_lvl5(n1,n2,n3)
    USE mo_submctl, only : prlim ! Note: #/m^3, but close enough to #/kg for statistics
    USE grid, ONLY : bulkNumc, bulkMixrat,meanRadius, dzt, &
        dn0, zm, a_ri, a_srs, snowin, a_rhi, a_tp, th0, th00, nspec

    IMPLICIT NONE

    integer, intent(in) :: n1,n2,n3

    REAL :: a0(n1,n2,n3), a1(n1,n2,n3), a2(n1,n2,n3)
    REAL :: scr(n2,n3), scr2(n2,n3), scr3(n2,n3)
    integer :: i, j, k, ii, sscnt
    LOGICAL :: mask(n1,n2,n3)

    ! Ice; combined and a and b separately for a and b regions, respectively
    CALL bulkNumc('ice','a',a0)
    CALL bulkNumc('ice','b',a1)
    ssclr_lvl5(1) = get_avg_ts(n1,n2,n3,a0+a1,dzt,icemask)
    CALL meanRadius('ice','ab',a2)
    ssclr_lvl5(2) = get_avg_ts(n1,n2,n3,a2,dzt,icemask)

    mask(:,:,:) = ( a0(:,:,:) > prlim .AND. icemask(:,:,:) )
    ssclr_lvl5(3) = get_avg_ts(n1,n2,n3,a0,dzt,mask)
    CALL meanRadius('ice','a',a2)
    ssclr_lvl5(4) = get_avg_ts(n1,n2,n3,a2,dzt,mask)

    mask(:,:,:) = ( a1(:,:,:) > prlim .AND. icemask(:,:,:) )
    ssclr_lvl5(5) = get_avg_ts(n1,n2,n3,a1,dzt,mask)
    CALL meanRadius('ice','b',a2)
    ssclr_lvl5(6) = get_avg_ts(n1,n2,n3,a2,dzt,mask)

    ! Snow
    CALL bulkNumc('snow','a',a0)
    ssclr_lvl5(7) = get_avg_ts(n1,n2,n3,a0,dzt,snowmask)
    CALL meanRadius('snow','ab',a0)
    ssclr_lvl5(8) = get_avg_ts(n1,n2,n3,a0,dzt,snowmask)


    ! IWP, SWP, max(ice) and max(snow), and the number of ice/snow cells
    ssclr_lvl5(10:11) = 0.
    ssclr_lvl5(13:14) = 0.
    scr = 0.   ! IWP
    scr2 = 0.   ! SWP
    scr3 = 0.   ! Integrated ice-liquid water potential temperature
    sscnt = 0
    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1
             if (icemask(k,i,j)) then
                ssclr_lvl5(10) = max(ssclr_lvl5(10), a_ri(k,i,j))
                ssclr_lvl5(11) = ssclr_lvl5(11) + 1.
             end if
             if (snowmask(k,i,j)) then
                ssclr_lvl5(13) = max(ssclr_lvl5(13), a_srs(k,i,j))
                ssclr_lvl5(14) = ssclr_lvl5(14) + 1.
             end if
             !
             ! Ice-water path (without snow)
             scr(i,j)=scr(i,j)+a_ri(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
             !
             ! Snow-water path
             scr2(i,j)=scr2(i,j)+a_srs(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
             !
             ! Integrated ice-liquid water potential temperature - change from th0
             scr3(i,j)=scr3(i,j)+(a_tp(k,i,j)+th00-th0(k))*(zm(k)-zm(k-1))
          end do
          !
          ! Surface snow rate
          if (snowmask(2,i,j)) sscnt=sscnt+1
          !
       end do
    end do
    ssclr_lvl5(9) = get_avg2dh(n2,n3,scr(:,:)) ! IWP
    ssclr_lvl5(12) = get_avg2dh(n2,n3,scr2(:,:)) ! SWP
    !
    ssclr_lvl5(15) = REAL(sscnt)/REAL( (n3-4)*(n2-4) )
    scr2(:,:) = snowin(2,:,:)
    ssclr_lvl5(16) = get_avg2dh(n2,n3,scr2)
    !
    ssclr_lvl5(18) = get_avg2dh(n2,n3,scr3)
    !
    ! Maximum supersaturation over ice
    ssclr_lvl5(17) = (MAXVAL(a_rhi(2:n1,3:n2-2,3:n3-2))-1.0)*100.

    ! Species mixing ratios in ice and snow
    DO i=1,nspec+1 ! Aerosol species and water
        ii = (i-1)*4 + 3
        CALL bulkMixrat(i,'ice','ab',a0)
        ssclr_mixr(ii) = get_avg_ts(n1,n2,n3,a0,dzt,icemask)
        CALL bulkMixrat(i,'snow','ab',a0)
        ssclr_mixr(ii+1) = get_avg_ts(n1,n2,n3,a0,dzt,snowmask)
    END DO

    ! Removal statistics elsewhere ..

  END SUBROUTINE ts_lvl5
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
  ! SUBROUTINE ACCUM_STAT: Accumulates various statistics over an
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
  ! SUBROUTINE ACCUM_LVL1: Accumulates various statistics over an
  ! averaging period for moisture variable (smoke or total water)
  !
  subroutine accum_lvl1(n1,n2,n3,q)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)  :: q(n1,n2,n3)

    integer :: k
    real    :: a1(n1),a2(n1),a3(n1)

    call get_avg3(n1,n2,n3,q,a1)
    call get_var3(n1,n2,n3,q,a1,a2)
    CALL get_3rd3(n1,n2,n3,q,a1,a3)

    do k=1,n1
       svctr(k,50)=svctr(k,50) + a1(k)
       svctr(k,51)=svctr(k,51) + a2(k)
       svctr(k,52)=svctr(k,52) + a3(k)
    end do

  end subroutine accum_lvl1
  !
  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_LVL2: Accumulates specialized statistics that depend
  ! on level 2 variables.
  !
  subroutine accum_lvl2(n1, n2, n3, th00, w, th, t, &
       rv, rl, rt)

    use defs, only : ep2

    integer, intent (in) :: n1,n2,n3
    real, intent (in)                       :: th00
    real, intent (in), dimension(n1,n2,n3)  :: w, th, t, rv, rl, rt

    real, dimension(n1,n2,n3) :: tv    ! Local variable
    integer                   :: k, i, j, kp1
    real, dimension(n1)       :: a1, a2, a3, tvbar
    real, dimension(n1,n2,n3) :: xy1, xy2, tw, tvw, rtw
    LOGICAL :: cond(n1,n2,n3)

    !
    ! liquid water statistics
    !
    call get_avg3(n1,n2,n3,rl,a1)
    call get_var3(n1,n2,n3,rl,a1,a2)
    call get_3rd3(n1,n2,n3,rl,a1,a3)
    svctr(:,59)=svctr(:,59) + a1(:)
    svctr(:,60)=svctr(:,60) + a2(:)
    svctr(:,61)=svctr(:,61) + a3(:)

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

  end subroutine accum_lvl2
  !
  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_LVL3: Accumulates specialized statistics that depend
  ! on level 3 variables.
  !
  subroutine accum_lvl3(n1, n2, n3, dn0, zm, rc, rr, nr, rrate, CCN)

    IMPLICIT NONE

    integer, intent (in) :: n1,n2,n3
    real, intent (in)                      :: CCN
    real, intent (in), dimension(n1)       :: zm, dn0
    real, intent (in), dimension(n1,n2,n3) :: rc, rr, nr, rrate

    integer                :: k, i, j
    real                   :: nrsum, nrcnt, rrcnt, rrcb
    real                   :: rmax, rmin
    real, dimension(n1)    :: a1
    real, dimension(n2,n3) :: scr2
    REAL :: mask(n1,n2,n3), tmp(n1,n2,n3)
    LOGICAL :: below

    !
    ! Average rain water mixing ratio
    !
    call get_avg3(n1,n2,n3,rr,a1)
    svctr(:,86)=svctr(:,86) + a1(:) ! rr (kg/kg)
    svctr(:,84)=svctr(:,84) + CCN   ! nc (#/kg)

    !
    ! conditionally average rain droplet concentrations
    !
    WHERE (rr > 0.001e-3)
       mask = 1.
    ELSEWHERE
       mask = 0.
    END WHERE

    call get_avg3(n1,n2,n3,nr,a1,cond=(mask>0.5))
    svctr(:,85)=svctr(:,85)+a1(:)
    call get_avg3(n1,n2,n3,mask,a1)
    svctr(:,91)=svctr(:,91)+a1(:)

    !
    ! precipitation flux
    !
    call get_avg3(n1,n2,n3,rrate,a1)
    svctr(:,87)=svctr(:,87)+a1(:)

    !
    ! conditionally average precip fluxes
    !
    WHERE (rrate > 3.65e-5)
       mask = 1.
    ELSEWHERE
       mask = 0.
    END WHERE

    tmp(:,:,:)=mask(:,:,:)*rrate(:,:,:)
    call get_avg3(n1,n2,n3,tmp,a1)
    svctr(:,90)=svctr(:,90)+a1(:)
    call get_avg3(n1,n2,n3,mask,a1)
    svctr(:,89)=svctr(:,89)+a1(:)

    !
    ! Histogram of surface rain rates
    !
    mask = 0.
    do k=1,n1
       rmin = max(6.2e-8,(k-1)*3.421e-5)
       rmax =  k * 3.421e-5
       do j=3,n3-2
          do i=3,n2-2
             if (rrate(2,i,j) > rmin .and. rrate(2,i,j) <= rmax) mask(k,i,j)=1.
          end do
       end do
    end do
    call get_avg3(n1,n2,n3,mask,a1,normalize=.FALSE.)
    svctr(:,92)=svctr(:,92)+a1(:)

    !
    ! Temporal statistics
    !
    scr2(:,:) = 0.
    nrsum = 0.
    nrcnt = 0.
    rrcnt = 0.
    rrcb = 0.
    do j=3,n3-2
       do i=3,n2-2
          below=.TRUE.
          do k=2,n1

             ! RWP
             scr2(i,j)=scr2(i,j)+rr(k,i,j)*dn0(k)*(zm(k)-zm(k-1))

             ! Rainy grid cell
             if (rr(k,i,j) > 0.001e-3) then
                nrsum = nrsum + nr(k,i,j)
                nrcnt = nrcnt + 1.
             end if

             ! Surface precipitation for this column
             if (k==2 .AND. rrate(k,i,j) > 3.65e-5) rrcnt = rrcnt + 1.

             ! Precpitation at cloud base (no cloud = no precip.)
             if (rc(k,i,j) > 1.e-5 .AND. below) then
                ! Take precpitation from level k-1 (>=2), which is just below cloud base
                rrcb = rrcb + rrate(max(2,k-1),i,j)
                below=.FALSE.
             ENDIF
          end do
       end do
    end do
    ssclr(24) = rrcnt/REAL( (n3-4)*(n2-4) )
    ssclr(22) = get_avg2dh(n2,n3,scr2)
    scr2(:,:) = rrate(2,:,:)
    ssclr(23) = get_avg2dh(n2,n3,scr2)
    ssclr(25) = CCN
    IF (nrcnt>0.) ssclr(26) = nrsum/nrcnt
    ssclr(27) = nrcnt
    ssclr(29) = rrcb/REAL( (n3-4)*(n2-4) )

  end subroutine accum_lvl3

  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_LVL4: Accumulates specialized statistics that depend
  ! on level 4 variables.
  !
  subroutine accum_lvl4(n1,n2,n3)
    use mo_submctl, only : fn2a,fn2b,inp2a,fnp2a,inp2b,fnp2b,cldbinlim,nout_cld, &
                     nlim,prlim ! Note: nlim and prlim in #/m^3, but close enough to #/kg for statistics
    use grid, ONLY : bulkNumc, bulkMixrat, meanRadius, binSpecMixrat, getBinRadius, &
                     a_rc, a_srp, a_rp, cldin, nspec, nbins, ncld, nprc, ngases, &
                     a_naerop, a_mcloudp, a_ncloudp, a_nprecpp, a_gaerop

    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3
    INTEGER :: ii,bb,i

    REAL, DIMENSION(n1,n2,n3)           :: a1,a12
    REAL, DIMENSION(n1,5)               :: a2
    REAL, DIMENSION(n1)                 :: col(n1)
    REAL, DIMENSION(n1,fn2a)            :: a3_a
    REAL, DIMENSION(n1,fn2b-fn2a)       :: a3_b
    REAL, DIMENSION(n1,fnp2a)           :: a4_a
    REAL, DIMENSION(n1,fnp2b-fnp2a)     :: a4_b
    REAL, DIMENSION(n1,nprc)            :: a5
    REAL, DIMENSION(n1,n2,n3,fnp2b)     :: a_Rwet
    REAL, ALLOCATABLE                   :: hist(:,:)


    ! *************************
    ! Bulk output for SALSA
    ! *************************
    CALL bulkNumc('aerosol','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,1))

    CALL bulkNumc('aerosol','b',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,2))

    CALL bulkNumc('cloud','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,3))

    CALL bulkNumc('cloud','b',a12)
    CALL get_avg3(n1,n2,n3,a12,a2(:,4))

    ! Generate cloud mask (grid cells with cloud)
    cloudmask(:,:,:) = ( a1(:,:,:)+a12(:,:,:) > nlim .AND. a_rc(:,:,:) > 1.e-5 )

    CALL bulkNumc('precp','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,5))

    ! Generate drizzle mask (grid cells with precipitation)
    drizzmask(:,:,:) = ( a1(:,:,:) > prlim .AND. a_srp(:,:,:) > 1.e-8 )

    svctr_lvl4(:,1:5) = svctr_lvl4(:,1:5) + a2(:,1:5)

    ! Particle radius
    CALL meanRadius('aerosol','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,1))

    CALL meanRadius('aerosol','b',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,2))

    ! In-cloud
    CALL meanRadius('cloud','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,3),cond=cloudmask)

    CALL meanRadius('cloud','b',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,4),cond=cloudmask)

    ! In-drizzle
    CALL meanRadius('precp','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,5),cond=drizzmask)

    svctr_lvl4(:,6:10) = svctr_lvl4(:,6:10) + a2(:,1:5)

    ! Water vapor mixing ratio
    CALL get_avg3(n1,n2,n3,a_rp,a2(:,1))

    ! Cloud water deposition flux
    call get_avg3(n1,n2,n3,cldin,a2(:,2))

    svctr_lvl4(:,11:12) = svctr_lvl4(:,11:12) + a2(:,1:2)

    ! Bin number concentrations
    ! -------------------------
    IF (lbinprof) THEN
        DO bb = 1,nbins
           IF (bb<=fn2a) THEN
              CALL get_avg3(n1,n2,n3,a_naerop(:,:,:,bb),a3_a(:,bb))
           ELSE
              CALL get_avg3(n1,n2,n3,a_naerop(:,:,:,bb),a3_b(:,bb-fn2a))
           ENDIF
        END DO
        DO bb = 1,ncld
           IF (bb<=fnp2a) THEN
              CALL get_avg3(n1,n2,n3,a_ncloudp(:,:,:,bb),a4_a(:,bb))
           ELSE
              CALL get_avg3(n1,n2,n3,a_ncloudp(:,:,:,bb),a4_b(:,bb-fnp2a))
           ENDIF
        END DO
        DO bb = 1,nprc
           CALL get_avg3(n1,n2,n3,a_nprecpp(:,:,:,bb),a5(:,bb))
        END DO

        svctr_aa(:,:,1) = svctr_aa(:,:,1) + a3_a(:,:)
        svctr_ab(:,:,1) = svctr_ab(:,:,1) + a3_b(:,:)
        svctr_ca(:,:,1) = svctr_ca(:,:,1) + a4_a(:,:)
        svctr_cb(:,:,1) = svctr_cb(:,:,1) + a4_b(:,:)
        svctr_p(:,:,1) = svctr_p(:,:,1) + a5(:,:)
    END IF

    ! Species mixing ratios
    ! ---------------------
    DO i=1,nspec+1 ! Aerosol species and water
        IF(lmixrprof) THEN
            ! Total mass mixing ratios
            CALL bulkMixrat(i,'aerosol','ab',a1)
            CALL get_avg3(n1,n2,n3,a1,a2(:,1))

            ! In-cloud
            CALL bulkMixrat(i,'cloud','ab',a1)
            CALL get_avg3(n1,n2,n3,a1,a2(:,2),cond=cloudmask)

            ! In-drizzle
            CALL bulkMixrat(i,'precp','a',a1)
            CALL get_avg3(n1,n2,n3,a1,a2(:,3),cond=drizzmask)

            ii = (i-1)*5+1 ! Order: aer, cld, prc, ice, snw
            svctr_mixr(:,ii:ii+2) = svctr_mixr(:,ii:ii+2) + a2(:,1:3)
        END IF

        ! Bin mixing ratios for all species
        IF (lbinprof) THEN
            DO bb = 1,nbins
                CALL binSpecMixrat('aerosol',i,bb,a1)
                IF (bb<=fn2a) THEN
                   CALL get_avg3(n1,n2,n3,a1,a3_a(:,bb))
                ELSE
                   CALL get_avg3(n1,n2,n3,a1,a3_b(:,bb-fn2a))
                ENDIF
            END DO
            DO bb = 1,ncld
                CALL binSpecMixrat('cloud',i,bb,a1)
                IF (bb<=fnp2a) THEN
                   CALL get_avg3(n1,n2,n3,a1,a4_a(:,bb),cond=cloudmask)
                ELSE
                   CALL get_avg3(n1,n2,n3,a1,a4_b(:,bb-fnp2a),cond=cloudmask)
                ENDIF
            END DO
            DO bb = 1,nprc
                CALL binSpecMixrat('precp',i,bb,a1)
                CALL get_avg3(n1,n2,n3,a1,a5(:,bb),cond=drizzmask)
            END DO

            svctr_aa(:,:,i+1) = svctr_aa(:,:,i+1) + a3_a(:,:)
            svctr_ab(:,:,i+1) = svctr_ab(:,:,i+1) + a3_b(:,:)
            svctr_ca(:,:,i+1) = svctr_ca(:,:,i+1) + a4_a(:,:)
            svctr_cb(:,:,i+1) = svctr_cb(:,:,i+1) + a4_b(:,:)
            svctr_p(:,:,i+1) = svctr_p(:,:,i+1) + a5(:,:)
        END IF
    END DO

    ! Gases
    DO i=1,ngases
        CALL get_avg3(n1,n2,n3,a_gaerop(:,:,:,i),col)
        svctr_gas(:,i) = svctr_gas(:,i) + col(:)
    END DO

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
                     a_ri, a_srs, nspec, nice, nsnw, a_micep, a_nicep, a_nsnowp, icein, snowin, a_tp, th00

    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3
    INTEGER :: ii,bb,i

    REAL, DIMENSION(n1,n2,n3)           :: a1,a12
    REAL, DIMENSION(n1,3)               :: a2
    REAL, DIMENSION(n1,fnp2a)           :: a4_a
    REAL, DIMENSION(n1,fnp2b-fnp2a)     :: a4_b
    REAL, DIMENSION(n1,nsnw)            :: a5
    REAL, DIMENSION(n1,n2,n3,fnp2b)     :: a_Rwet
    REAL, ALLOCATABLE                   :: hist(:,:)

    ! *************************
    ! Bulk output for SALSA
    ! *************************
    CALL bulkNumc('ice','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,1))

    CALL bulkNumc('ice','b',a12)
    CALL get_avg3(n1,n2,n3,a12,a2(:,2))

    ! Generate ice mask (grid cells with ice)
    icemask(:,:,:) = ( a1(:,:,:)+a12(:,:,:) > prlim .AND. a_ri(:,:,:) > 1.e-8)

    CALL bulkNumc('snow','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,3))

    ! Generate snow mask (grid cells with snow)
    snowmask(:,:,:) = ( a1(:,:,:) > prlim .AND. a_srs(:,:,:) > 1.e-11 )

    svctr_lvl5(:,1:3) = svctr_lvl5(:,1:3) + a2(:,1:3)

    ! Ice radius
    CALL meanRadius('ice','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,1),cond=icemask)

    CALL meanRadius('ice','b',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,2),cond=icemask)

    ! Snow radius
    CALL meanRadius('snow','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,3),cond=snowmask)

    svctr_lvl5(:,4:6) = svctr_lvl5(:,4:6) + a2(:,1:3)

    ! Ice deposition flux
    call get_avg3(n1,n2,n3,icein,a2(:,1))

    ! Snow deposition flux
    call get_avg3(n1,n2,n3,snowin,a2(:,2))

    ! Ice-liquid water potential temperature
    call get_avg3(n1,n2,n3,a_tp,a2(:,3))
    a2(:,3) = a2(:,3) + th00

    svctr_lvl5(:,7:9)=svctr_lvl5(:,7:9)+a2(:,1:3)

    ! Bin number concentrations
    ! -------------------------
    IF (lbinprof) THEN
        DO bb = 1,nice
            IF (bb<=fnp2a) THEN
                CALL get_avg3(n1,n2,n3,a_nicep(:,:,:,bb),a4_a(:,bb))
            ELSE
                CALL get_avg3(n1,n2,n3,a_nicep(:,:,:,bb),a4_b(:,bb-fnp2a))
            ENDIF
        END DO
        DO bb = 1,nsnw
           CALL get_avg3(n1,n2,n3,a_nsnowp(:,:,:,bb),a5(:,bb))
        END DO

        svctr_ia(:,:,1) = svctr_ia(:,:,1) + a4_a(:,:)
        svctr_ib(:,:,1) = svctr_ib(:,:,1) + a4_b(:,:)
        svctr_s(:,:,1) = svctr_s(:,:,1) + a5(:,:)
    ENDIF

    ! Species mixing ratios
    ! ---------------------
    DO i=1,nspec+1 ! Aerosol species and water
        IF(lmixrprof) THEN
            ! In-ice
            CALL bulkMixrat(i,'ice','ab',a1)
            CALL get_avg3(n1,n2,n3,a1,a2(:,1),cond=icemask)

            ! In-snow
            CALL bulkMixrat(i,'snow','a',a1)
            CALL get_avg3(n1,n2,n3,a1,a2(:,2),cond=snowmask)

            ii = (i-1)*5+4 ! Order: aer, cld, prc, ice, snw
            svctr_mixr(:,ii:ii+1) = svctr_mixr(:,ii:ii+1) + a2(:,1:2)
        ENDIF

        ! Binned mixing ratios for all species
        IF (lbinprof) THEN
            DO bb = 1,nice
                CALL binSpecMixrat('ice',i,bb,a1)
                IF (bb<=fnp2a) THEN
                    CaLL get_avg3(n1,n2,n3,a1,a4_a(:,bb),cond=icemask)
                ELSE
                    CALL get_avg3(n1,n2,n3,a1,a4_b(:,bb-fnp2a),cond=icemask)
                ENDIF
            END DO
            DO bb = 1,nsnw
                CALL binSpecMixrat('snow',i,bb,a1)
                CALL get_avg3(n1,n2,n3,a1,a5(:,bb),cond=snowmask)
            END DO

            svctr_ia(:,:,i+1) = svctr_ia(:,:,i+1) + a4_a(:,:)
            svctr_ib(:,:,i+1) = svctr_ib(:,:,i+1) + a4_b(:,:)
            svctr_s(:,:,i+1) = svctr_s(:,:,i+1) + a5(:,:)
       END IF

    END DO

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
    real fact, uw_shear, vw_shear

    fact = 0.25/float((n2-4)*(n3-4))

    call get_avg3(n1,n2,n3,u,ub)
    call get_avg3(n1,n2,n3,v,vb)

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             uw_shear = -(u(k,i,j)-ub(k))*fact*(                              &
                  (w(k,i,j)  +w(k,i+1,j)  )*(ub(k+1)-ub(k)  )*dzm(k) +        &
                  (w(k-1,i,j)+w(k-1,i+1,j))*(ub(k)  -ub(k-1))*dzm(k-1))
             if (j > 1) vw_shear = -(v(k,i,j)-vb(k))*fact*(                   &
                  (w(k,i,j)  +w(k,i,j+1)  )*(vb(k+1)-vb(k)  )*dzm(k) +        &
                  (w(k-1,i,j)+w(k-1,i,j+1))*(vb(k)  -vb(k-1))*dzm(k-1))

             svctr(k,48) = svctr(k,48)+uw_shear
             svctr(k,36) = svctr(k,36)+uw_shear+vw_shear
          end do
       end do
    end do

  end subroutine get_shear
  !
  ! ----------------------------------------------------------------------
  ! Subroutine write_ts: writes the statistics file
  !
  subroutine write_ts

    use netcdf
    use mpi_interface, only : myid
    USE grid, ONLY : ngases

    integer :: iret, n, VarID

    do n=1,nvar1
       iret = nf90_inq_varid(ncid1, s1(n), VarID)
       IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, ssclr(n), start=(/nrec1/))
    end do
    ssclr(:) = 0.

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

    IF (nv1_mixr>0) THEN
       DO n = 1,nv1_mixr
          iret = nf90_inq_varid(ncid1, s1_mixr(n), VarID)
          IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, ssclr_mixr(n), start=(/nrec1/))
       END DO
       ssclr_mixr(:) = 0.
    END IF

    IF (nv1_rem>0) THEN
       DO n = 1,nv1_rem
          iret = nf90_inq_varid(ncid1, s1_rem(n), VarID)
          IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, ssclr_rem(n), start=(/nrec1/))
       END DO
       ssclr_rem(:) = 0.
    END IF

    ! Gas concentrations
    DO n=1,ngases
       iret = nf90_inq_varid(ncid1, s1_gas(n), VarID)
       IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, ssclr_gas(n), start=(/nrec1/))
       ssclr_gas(n) = 0.
    END DO

    ! User-selected process rate outputs
    IF (nv1_proc>0) THEN
       DO n=1,nv1_proc
          iret = nf90_inq_varid(ncid1, out_ts_list(n), VarID)
          IF (iret == NF90_NOERR) iret = nf90_put_var(ncid1, VarID, out_ts_data(n), start=(/nrec1/))
        ENDDO
        out_ts_data(:) = 0.
    ENDIF

    if (myid==0) print "(/' ',12('-'),'   Record ',I3,' to time series')",nrec1

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
    USE mo_submctl, ONLY : in1a,in2a,fn2a,fn2b,fnp2a,fnp2b, &
                               aerobins,precpbins,snowbins, &
                               nout_cld, cldbinlim, nout_ice, icebinlim
    USE grid, ONLY : nprc,nsnw,ngases
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
       IF (level >=4 ) THEN
          svctr_lvl4(k,:) = svctr_lvl4(k,:)/nsmp
          svctr_mixr(k,:) = svctr_mixr(k,:)/nsmp
          svctr_aa(k,:,:) = svctr_aa(k,:,:)/nsmp
          svctr_ab(k,:,:) = svctr_ab(k,:,:)/nsmp
          svctr_ca(k,:,:) = svctr_ca(k,:,:)/nsmp
          svctr_cb(k,:,:) = svctr_cb(k,:,:)/nsmp
          svctr_p(k,:,:) = svctr_p(k,:,:)/nsmp

          ! Replace level 3 CDNC=CCN with that from SALSA (a + b bins)
          svctr(k,84)=svctr_lvl4(k,3)+svctr_lvl4(k,4)
       END IF
       IF (level >= 5) THEN
          svctr_lvl5(k,:) = svctr_lvl5(k,:)/nsmp
          svctr_ia(k,:,:) = svctr_ia(k,:,:)/nsmp
          svctr_ib(k,:,:) = svctr_ib(k,:,:)/nsmp
          svctr_s(k,:,:) = svctr_s(k,:,:)/nsmp
       ENDIF

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
          iret = nf90_inq_varid(ncid2,'P_Rd12a',VarID) ! 1a+2a
          IF (iret == NF90_NOERR) &
                iret = nf90_put_var(ncid2,VarID,aerobins(in1a:fn2a),start=(/nrec2/))
          iret = nf90_inq_varid(ncid2,'P_Rd2ab',VarID) ! 2a = 2b
          IF (iret == NF90_NOERR) &
                iret = nf90_put_var(ncid2,VarID,aerobins(in2a:fn2a),start=(/nrec2/))
          iret = nf90_inq_varid(ncid2,'P_Rwprc',VarID)
          IF (iret == NF90_NOERR) &
                iret = nf90_put_var(ncid2,VarID,precpbins(1:nprc),start=(/nrec2/))
          iret = nf90_inq_varid(ncid2,'P_Rwsnw',VarID)
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
       iret = nf90_put_var(ncid2,VarID,svctr(:,n), start=(/1,nrec2/),    &
            count=(/n1,1/))
    end do

    IF (level >= 4) THEN
       ! SALSA level 4
       DO n = 1,nv2_lvl4
          iret = nf90_inq_varid(ncid2,s2_lvl4(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_lvl4(:,n), start=(/1,nrec2/), count=(/n1,1/))
       END DO

       ! Mixing ratios
       DO n = 1,nv2_mixr
          iret = nf90_inq_varid(ncid2,s2_mixr(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_mixr(:,n), start=(/1,nrec2/), count=(/n1,1/))
       END DO

       ! Binned aerosols, regime a
       DO n = 1,nv2_bin
          iret = nf90_inq_varid(ncid2,s2_aa(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_aa(:,:,n), start=(/1,1,nrec2/), count=(/n1,fn2a,1/))
       END DO

       ! Aerosols, regime b
       DO n = 1,nv2_bin
          iret = nf90_inq_varid(ncid2,s2_ab(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_ab(:,:,n), start=(/1,1,nrec2/), count=(/n1,fn2b-fn2a,1/))
       END DO

       ! Cloud droplets, regime a
       DO n = 1,nv2_bin
          iret = nf90_inq_varid(ncid2,s2_ca(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_ca(:,:,n), start=(/1,1,nrec2/), count=(/n1,fnp2a,1/))
       END DO

       ! Cloud droplets, regime b
       DO n = 1,nv2_bin
          iret = nf90_inq_varid(ncid2,s2_cb(n),VarID)
          IF (iret == NF90_NOERR) &
            iret = nf90_put_var(ncid2,VarID,svctr_cb(:,:,n), start=(/1,1,nrec2/), count=(/n1,fnp2b-fnp2a,1/))
       END DO

       ! Precipitation
       DO n = 1,nv2_bin
          iret = nf90_inq_varid(ncid2,s2_p(n),VarID)
          IF (iret == NF90_NOERR) &
            iret = nf90_put_var(ncid2,VarID,svctr_p(:,:,n), start=(/1,1,nrec2/), count=(/n1,nprc,1/))
       END DO

    END IF

    IF (level >= 5) THEN
       ! SALSA level 5
       DO n = 1,nv2_lvl5
          iret = nf90_inq_varid(ncid2,s2_lvl5(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_lvl5(:,n), start=(/1,nrec2/), count=(/n1,1/))
       END DO

       ! Ice, regime a
       DO n = 1,nv2_bin
          iret = nf90_inq_varid(ncid2,s2_ia(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_ia(:,:,n), start=(/1,1,nrec2/), count=(/n1,fnp2a,1/))
       END DO

       ! Ice, regime b
       DO n = 1,nv2_bin
          iret = nf90_inq_varid(ncid2,s2_ib(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_ib(:,:,n), start=(/1,1,nrec2/), count=(/n1,fnp2b-fnp2a,1/))
       END DO

       ! Snow
       DO n = 1,nv2_bin
          iret = nf90_inq_varid(ncid2,s2_s(n),VarID)
          IF (iret == NF90_NOERR) &
             iret = nf90_put_var(ncid2,VarID,svctr_s(:,:,n), start=(/1,1,nrec2/), count=(/n1,nsnw,1/))
       END DO

    END IF

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

    IF (level >= 4 .AND. ngases>0) THEN
        DO n=1,ngases
            iret = nf90_inq_varid(ncid2, s2_gas(n), VarID)
            IF (iret == NF90_NOERR) THEN
                svctr_gas(:,n)=svctr_gas(:,n)/nsmp
                iret = nf90_put_var(ncid2,VarID,svctr_gas(:,n), start=(/1,nrec2/),count=(/n1,1/))
            ENDIF
       END DO
       svctr_gas(:,:)=0.
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

    if (myid==0) print "(/' ',12('-'),'   Record ',I3,' to profiles')",nrec2

    iret  = nf90_sync(ncid2)
    nrec2 = nrec2+1
    nsmp  = 0.

    svctr(:,:) = 0.
    IF (level >= 4) THEN
        svctr_lvl4(:,:) = 0.
        svctr_mixr(:,:) = 0.
        svctr_aa(:,:,:) = 0.
        svctr_ab(:,:,:) = 0.
        svctr_ca(:,:,:) = 0.
        svctr_cb(:,:,:) = 0.
        svctr_p(:,:,:) = 0.
    END IF
    IF (level >=5 ) THEN
        svctr_lvl5(:,:) = 0.
        svctr_ia(:,:,:) = 0.
        svctr_ib(:,:,:) = 0.
        svctr_s(:,:,:) = 0.
    ENDIF

  end subroutine write_ps
  !
  ! ----------------------------------------------------------------------
  ! subroutine: sfc_stat:  Updates statistical arrays with surface flux
  ! variables
  !
  subroutine sfc_stat(n2,n3,tflx,qflx,ustar,sst)

    integer, intent(in) :: n2,n3
    real, intent(in), dimension(n2,n3) :: tflx, qflx, ustar
    real, intent(in)    :: sst

    ssclr(10) = sst
    ssclr(11) = get_avg2dh(n2,n3,ustar)

    ssclr(12) = get_avg2dh(n2,n3,tflx)
    if (level >= 1) ssclr(13) = get_avg2dh(n2,n3,qflx)

  end subroutine sfc_stat
  !
  ! ----------------------------------------------------------------------
  ! subroutine: fills scalar array based on index
  ! 1: cfl; 2 max divergence
  !
  subroutine fill_scalar(index,xval)

    integer, intent(in) :: index
    real, intent (in)   :: xval

    select case(index)
    case(1)
       ssclr(2) = xval
    case(2)
       ssclr(3) = xval
    end select

  end subroutine fill_scalar
  !
  ! ----------------------------------------------------------------------
  ! subroutine: calculates the dissipation for output diagnostics, if
  ! isgstyp equals 2 then le is passed in via diss
  !
  subroutine sgs_vel(n1,n2,n3,v1,v2,v3)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: v1(n1),v2(n1),v3(n1)

    svctr(:,23)=svctr(:,23)+v1(:)/float((n2-4)*(n3-4))
    svctr(:,25)=svctr(:,25)+v2(:)/float((n2-4)*(n3-4))
    svctr(:,27)=svctr(:,27)+v3(:)/float((n2-4)*(n3-4))

  end subroutine sgs_vel
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
    ! calculate fluxes assuming the possibility of liquid water.  if liquid
    ! water does not exist sgs_rl = 0.
    !
    if ( level >= 2 ) then
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
       !
       ! calculate fluxes for dry thermodynamics, i.e., wrl_sgs is by def
       ! zero
       !
    else
       do k = 1,n1
          wrl_sgs(k) = 0.
       end do
       do j = 3,n3-2
          do i = 3,n2-2
             do k = 1,n1-1
                if ( level >= 1) then
                   select case (type)
                   case ('tl')
                      fctt = rnpts * (1. + ep2*rv(k,i,j))
                   case ('rt')
                      fctt = rnpts * ep2*th(k,i,j)
                   end select
                else
                   fctt = rnpts
                end if
                wtv_sgs(k) = wtv_sgs(k) + fctt*flx(k,i,j)
             end do
          end do
       end do
    end if

  end subroutine sgsflxs
  !
  ! ----------------------------------------------------------------------
  ! subroutine fill_tend: fills arrays with current value of tendencies
  !
  subroutine acc_tend(n1,n2,n3,f1,f2,f3,t1,t2,t3,v1,v2,v3,ic,routine)

    integer, intent(in) :: n1,n2,n3,ic
    real, intent(in)    :: f1(n1,n2,n3),f2(n1,n2,n3),f3(n1,n2,n3)
    real, intent(in)    :: t1(n1,n2,n3),t2(n1,n2,n3),t3(n1,n2,n3)
    real, intent(inout) :: v1(n1),v2(n1),v3(n1)
    character (len=3)   :: routine

    integer :: k,ii
    real    :: x1(n1),x2(n1),x3(n1)

    call get_cor3(n1,n2,n3,f1,t1,x1)
    call get_cor3(n1,n2,n3,f2,t2,x2)
    call get_cor3(n1,n2,n3,f3,t3,x3)

    select case (routine)
    case ('sgs')
       ii = 39
    case ('adv')
       ii = 42
    end select

    select case (ic)
    case (1)
       do k=1,n1
          v1(k) = x1(k)
          v2(k) = x2(k)
          v3(k) = x3(k)
       end do
    case (2)
       do k=1,n1
          svctr(k,ii)   = svctr(k,ii)   + (x1(k)-v1(k))
          svctr(k,ii+1) = svctr(k,ii+1) + (x2(k)-v2(k))
          svctr(k,ii+2) = svctr(k,ii+2) + (x3(k)-v3(k))
       end do
    end select

  end subroutine acc_tend
  !
  !---------------------------------------------------------------------
  ! subroutine updtst: updates appropriate statistical arrays
  !
  subroutine updtst(n1,routine,nfld,values,ic)

    integer, intent(in)            :: n1,nfld,ic
    real, intent (in)              :: values(n1)
    character (len=3), intent (in) :: routine

    integer :: nn,k

    select case (routine)
    case("sgs")
       select case (nfld)
       case (-6)
          nn = 31 ! dissipation length-scale
       case (-5)
          nn = 30 ! mixing length
       case (-4)
          nn = 29 ! eddy diffusivity
       case (-3)
          nn = 28 ! eddy viscosity
       case (-2)
          nn = 38 ! dissipation
       case (-1)
          nn = 32 ! estimated sgs energy
       case (1)
          nn = 21 ! sgs tl flux
       case (2)
          nn = 54 ! sgs rt flux
       case default
          nn = 0
       end select
    case("adv")
       select case (nfld)
       case (-3)
          nn = 26 ! adv w flux
       case (-2)
          nn = 24 ! adv v flux
       case (-1)
          nn = 22 ! adv u flux
       case (0)
          nn = 62 ! adv rl flux
       case (1)
          nn = 20 ! adv tl flux
       case (2)
          nn = 53 ! adv rt flux
       case default
          nn = 0
       end select
    case("prs")
       select case (nfld)
       case (1)
          nn = 45 ! dpdx u corr
       case (2)
          nn = 46 ! dpdy v corr
       case (3)
          nn = 47 ! dpdz w corr
       case default
          nn = 0
       end select
    case("prc")
       select case (nfld)
       case (2)
          nn = 88
       case (3)
          nn = 63
       case default
          nn = 0
       end select
    case default
       nn = 0
    end select

    if (nn > 0) then
       if (ic == 0) svctr(:,nn)=0.
       do k=1,n1
          svctr(k,nn)=svctr(k,nn)+values(k)
       enddo
    end if

  end subroutine updtst


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
    IF (level==3) THEN
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
        CASE ('g')
            ! Gas: scalar so just mass mixing ratio of component j
            READ(UNIT=tchar,FMT='(I1)') j
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
  !
  ! -------------------------------------------------------------------------
  !
  real function get_zi (n1, n2, n3, itype, sx, xx, z, threshold)

    integer, intent (in) :: n1, n2, n3, itype
    real, intent (in)    :: xx(n1), z(n1), sx(n1,n2,n3), threshold

    integer :: i, j, k, kk
    real    :: zibar, sval, dmy, scr(n2,n3)

    get_zi = -999.
    select case(itype)
    case (1)
       !
       ! find level at which sx=threshold (xx is one over grid spacing)
       !
       zibar = 0.
       do j=3,n3-2
          do i=3,n2-2
             k = 2
             do while (k < n1-2 .and. sx(k,i,j) > threshold)
                k = k+1
             end do
             if (k == n1-2) zibar = -999.
             if (zibar /= -999.) zibar = zibar + z(k-1) +  &
                  (threshold - sx(k-1,i,j))/xx(k-1)     /  &
                  (sx(k,i,j) - sx(k-1,i,j) + epsilon(1.))
          end do
       end do
       if (zibar /= -999.) get_zi = zibar/real((n3-4)*(n2-4))

    case(2)
       !
       ! find level of maximum gradient (xx is one over grid spacing)
       !
       scr=0.
       do j=3,n3-2
          do i=3,n2-2
             sval = 0.
             do k=2,n1-5
                dmy = (sx(k+1,i,j)-sx(k,i,j))*xx(k)
                if (dmy > sval) then
                   sval = dmy
                   scr(i,j) = z(k)
                end if
             end do
          end do
       end do
       get_zi = get_avg2dh(n2,n3,scr)

    case(3)
       !
       ! find level where xx is a maximum
       !
       sval = -huge(1.)
       kk = 1
       do k=2,n1
          if (xx(k) > sval) then
             kk = k
             sval = xx(k)
          end if
       end do
       get_zi = z(kk)

    case(4)
       !
       ! find level where xx is a minimum
       !
       sval = huge(1.)
       kk = 1
       do k=2,n1-2
          if (xx(k) < sval) then
             kk = k
             sval = xx(k)
          end if
       end do
       get_zi = z(kk)
    end select

  end function get_zi

end module stat

