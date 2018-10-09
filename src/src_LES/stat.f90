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
MODULE stat

   USE mo_submctl, ONLY : spec
   USE ncio, ONLY : open_nc, define_nc, define_nc_cs
   USE grid, ONLY : level, lbinprof   ! IT WOULD BE BEST TO IMPORT EVERYTHING NEEDED HERE ONCE AND USE THAT FOR THE WHOLE MODULE
   USE util, ONLY : get_avg3, get_cor3, get_var3, get_avg_ts, &
                    get_avg2dh, get_3rd3, getMassIndex


   IMPLICIT NONE
   PRIVATE

  integer, parameter :: nvar1 = 29,               &
                        nv1_lvl4 = 53,            &
                        nv1MB = 4,                &
                        nv1_lvl5 = 49,            &
                        nvar2 = 96,               &
                        nv2_lvl4 = 43,            &
                        nv2_lvl5 = 29, &
                        nv2saa = 8, nv2sab = 8,   &
                        nv2sca = 8, nv2scb = 8,   &
                        nv2sp = 8, &
                        nv2sia = 8, nv2sib = 8,   &
                        nv2ss = 8

  ! DEPRECATED
  ! All SALSA species
  CHARACTER(len=3), PARAMETER :: zspec(8) = (/'SO4','OC ','BC ','DU ','SS ','NO ','NH ','H2O'/)
  ! Active SALSA species
  character (len=3), save :: actspec(8)
  integer, save      :: nspec_depr=0

   INTEGER, SAVE      :: nrec1, nrec2, nrec3, ncid1, ncid2, ncid3, nv1=nvar1, nv2 = nvar2
   REAL, SAVE         :: fsttm, lsttm, nsmp = 0
   REAL, SAVE         :: avgtime = 0

   LOGICAL            :: sflg = .FALSE.
   LOGICAL            :: mcflg = .FALSE.
   LOGICAL            :: csflg = .FALSE.
   LOGICAL            :: salsa_b_bins = .FALSE.
   LOGICAL            :: cloudy_col_stats = .FALSE.
   REAL               :: ssam_intvl = 30.   ! statistical sampling interval
   REAL               :: savg_intvl = 1800. ! statistical averaging interval

   ! New output variables:
   ! 1. Add them to one of the arrays below
   ! 2. Update the dimensions of the output arrays accordingly
   ! 3. Make a subroutine for accumulating the new data and to plane them
   !    in the output arrays
   ! 4. Add the new variables in ncio.f90 list of variables
   ! 5. Make sure stat_init is up to date (boolean arrays etc).


   CHARACTER (len=7), SAVE :: s1(nvar1)=(/                         &
      'time   ','cfl    ','maxdiv ','zi1_bar','zi2_bar','zi3_bar', & ! 1
      'vtke   ','sfcbflx','wmax   ','tsrf   ','ustar  ','shf_bar', & ! 7
      'lhf_bar','zi_bar ','lwp_bar','lwp_var','zc     ','zb     ', & !13
      'cfrac  ','lmax   ','albedo ','rwp_bar','prcp   ','pfrac  ', & !19
      'CCN    ','nrain  ','nrcnt  ','nccnt  ','prcp_bc'/),         & !25

       ! **** Bulk temporal statistics for SALSA ****
       s1_lvl4(nv1_lvl4) = (/       &
       'Nc_ic  ','Rc_ic  ','Nca_ica','Rca_ica','Ncb_icb','Rcb_icb', & !1
       'Na_int ','Ra_int ','Naa_int','Raa_int','Nab_int','Rab_int', & !7
       'SO4_ic ','SO4_int','OC_ic  ','OC_int ',   & !13
       'BC_ic  ','BC_int ','DU_ic  ','DU_int ',   & !17
       'SS_ic  ','SS_int ','NO_ic  ','NO_int ',   & !21
       'NH_ic  ','NH_int ','H2O_ic ','H2O_int',   & !25
       'rmSO4dr','rmSO4cl','rmSO4pr',   & !29
       'rmOCdr ','rmOCcl ','rmOCpr ',   & !32
       'rmBCdr ','rmBCcl ','rmBCpr ',   & !35
       'rmDUdr ','rmDUcl ','rmDUpr ',   & !38
       'rmSSdr ','rmSScl ','rmSSpr ',   & !41
       'rmNOdr ','rmNOcl ','rmNOpr ',   & !44
       'rmNHdr ','rmNHcl ','rmNHpr ',   & !47
       'rmH2Oae','rmH2Ocl','rmH2Opr',   & !50
       'SS_max '/), & !53

       s1_lvl5(nv1_lvl5) = (/  &
       'Ni_ii  ','Ri_ii  ','Nia_iia','Ria_iia','Nib_iib','Rib_iib', & ! 1-6
       'Ns_is  ','Rs_is  ',  & ! 7-8
       'iwp_bar','imax   ','nicnt  ',  & ! 9-11
       'swp_bar','smax   ','nscnt  ',  & ! 12-14
       'sfrac  ','sprcp  ','SSi_max',  & ! 15-17
       'SO4_ii ','SO4_is ','OC_ii  ','OC_is  ',   & ! 18-21
       'BC_ii  ','BC_is  ','DU_ii  ','DU_is  ',   & ! 22-25
       'SS_ii  ','SS_is  ','NO_ii  ','NO_is  ',   & ! 26-29
       'NH_ii  ','NH_is  ','H2O_ii ','H2O_is ',   & ! 30-33
       'rmSO4ic','rmSO4sn','rmOCic ','rmOCsn ',   & ! 34-37
       'rmBCic ','rmBCsn ','rmDUic ','rmDUsn ',   & ! 38-41
       'rmSSic ','rmSSsn ','rmNOic ','rmNOsn ',   & ! 42-45
       'rmNHic ','rmNHsn ','rmH2Oic','rmH2Osn'/), & ! 46-49

        s2(nvar2)=(/                                                 &
        'time   ','zt     ','zm     ','dn0    ','u0     ','v0     ', & ! 1
        'fsttm  ','lsttm  ','nsmp   ','u      ','v      ','theta  ', & ! 7
        'p      ','u_2    ','v_2    ','w_2    ','theta_2','w_3    ', & ! 13
        'theta_3','tot_tw ','sfs_tw ','tot_uw ','sfs_uw ','tot_vw ', & ! 19
        'sfs_vw ','tot_ww ','sfs_ww ','km     ','kh     ','lmbd   ', & ! 25
        'lmbde  ','sfs_tke','sfs_boy','sfs_shr','boy_prd','shr_prd', & ! 31
        'trans  ','diss   ','dff_u  ','dff_v  ','dff_w  ','adv_u  ', & ! 37
        'adv_v  ','adv_w  ','prs_u  ','prs_v  ','prs_w  ','prd_uw ', & ! 43
        'storage','q      ','q_2    ','q_3    ','tot_qw ','sfs_qw ', & ! 49
        'rflx   ','rflx2  ','sflx   ','sflx2  ','l      ','l_2    ', & ! 55
        'l_3    ','tot_lw ','sed_lw ','cs1    ','cnt_cs1','w_cs1  ', & ! 61
        't_cs1  ','tv_cs1 ','rt_cs1 ','rc_cs1 ','wt_cs1 ','wtv_cs1', & ! 67
        'wrt_cs1','cs2    ','cnt_cs2','w_cs2  ','t_cs2  ','tv_cs2 ', & ! 73
        'rt_cs2 ','rc_cs2 ','wt_cs2 ','wtv_cs2','wrt_cs2','Nc     ', & ! 79
        'Nr     ','rr     ','rrate  ','evap   ','frc_prc','prc_prc', & ! 85
        'frc_ran','hst_srf','sw_up  ','sw_down','lw_up  ','lw_down'/), & ! 91, total 96

        ! **** BULK PROFILE OUTPUT FOR SALSA ****
        s2_lvl4(nv2_lvl4) = (/                                       &
        'aea    ','aeb    ','cla    ','clb    ','prc    ',           & !1
        'P_Naa  ','P_Nab  ','P_Nca  ','P_Ncb  ','P_Np   ',           & !6
        'P_Rwaa ','P_Rwab ','P_Rwca ','P_Rwcb ','P_Rwp  ',           & !11
        'P_cSO4a','P_cSO4c','P_cSO4p',                               & !16
        'P_cOCa ','P_cOCc ','P_cOCp ',                               & !19
        'P_cBCa ','P_cBCc ','P_cBCp ',                               & !22
        'P_cDUa ','P_cDUc ','P_cDUp ',                               & !25
        'P_cSSa ','P_cSSc ','P_cSSp ',                               & !28
        'P_cNOa ','P_cNOc ','P_cNOp ',                               & !31
        'P_cNHa ','P_cNHc ','P_cNHp ',                               & !34
        'P_cH2Oa','P_cH2Oc','P_cH2Op',                               & !37
        'P_rl   ','P_rr   ','P_rv   ','P_RH   '/),                   & !40-43

        s2_lvl5(nv2_lvl5) = (/ &
        'ica    ','icb    ','snw    ',           & ! 1-3
        'P_Nia  ','P_Nib  ','P_Ns   ','P_Rwia ','P_Rwib ','P_Rws  ', & ! 4-9
        'P_cSO4i','P_cSO4s','P_cOCi ','P_cOCs ', & ! 10-13
        'P_cBCi ','P_cBCs ','P_cDUi ','P_cDUs ', & ! 14-17
        'P_cSSi ','P_cSSs ','P_cNOi ','P_cNOs ', & ! 18-21
        'P_cNHi ','P_cNHs ','P_cH2Oi','P_cH2Os', & ! 22-25
        'P_ri   ','P_rs   ', 'P_RHi  ','srate  '/),    & ! 26-29

        ! **** BINNED PROFILE OUTPUT FOR SALSA ****
        ! **** Aerosols
        s2Aeroa(nv2saa) = (/                                         &
        'P_Naba ','P_SO4aa','P_OCaa ','P_BCaa ',                     &
        'P_DUaa ','P_SSaa ','P_NOaa ','P_NHaa '/),                   &
        s2Aerob(nv2sab) = (/                                         &
        'P_Nabb ','P_SO4ab','P_OCab ','P_BCab ',                     &
        'P_DUab ','P_SSab ','P_NOab ','P_NHab '/),                   &

        ! **** Clouds
        s2Clouda(nv2sca) = (/                                        &
        'P_Ncba ','P_SO4ca','P_OCca ','P_BCca ',                     &
        'P_DUca ','P_SSca ','P_NOca ','P_NHca '/),                   &
        s2Cloudb(nv2scb) = (/                                        &
        'P_Ncbb ','P_SO4cb','P_OCcb ','P_BCcb ',                     &
        'P_DUcb ','P_SScb ','P_NOcb ','P_NHcb '/),                   &

        ! **** Precip
        s2Precp(nv2sp) = (/                                         &
        'P_Npb  ','P_SO4pb','P_OCpb ','P_BCpb ',                     &
        'P_DUpb ','P_SSpb ','P_NOpb ','P_NHpb '/),                   &

        ! **** Ice
        s2Icea(nv2sia) = (/                                        &
        'P_Niba ','P_SO4ia','P_OCia ','P_BCia ',                     &
        'P_DUia ','P_SSia ','P_NOia ','P_NHia '/),                   &
        s2Iceb(nv2sib) = (/                                        &
        'P_Nibb ','P_SO4ib','P_OCib ','P_BCib ',                     &
        'P_DUib ','P_SSib ','P_NOib ','P_NHib '/),                   &

        ! **** Snow
        s2Snow(nv2ss) = (/                                         &
        'P_Nsb  ','P_SO4sb','P_OCsb ','P_BCsb ',                     &
        'P_DUsb ','P_SSsb ','P_NOsb ','P_NHsb '/),                   &

        s1Total(nvar1+nv1_lvl4+nv1_lvl5),                            &
        s2Total(nvar2+nv2_lvl4+nv2_lvl5+nv2saa+nv2sab+nv2sca+nv2scb+nv2sp+nv2sia+nv2sib+nv2ss)


  LOGICAL, save :: s2bool(nvar2+nv2_lvl4+nv2_lvl5+nv2saa+nv2sab+nv2sca+nv2scb+nv2sp+nv2sia+nv2sib+nv2ss)
  LOGICAL, save :: s1bool(nvar1+nv1_lvl4+nv1_lvl5)

  real, save, allocatable   :: tke_sgs(:), tke_res(:), tke0(:), wtv_sgs(:),  &
       wtv_res(:), wrl_sgs(:), thvar(:), svctr(:,:), ssclr(:),               &
       ! Additional ssclr and svctr for BULK SALSA output
       svctr_b(:,:), ssclr_b(:),                                             &
       ! Additional ssclr and svctr for BINNED SALSA output.
       svctr_aa(:,:,:), svctr_ca(:,:,:), svctr_p(:,:,:),                     &
       svctr_ab(:,:,:), svctr_cb(:,:,:),                                     &
       ! The same for SALSA level 5
       svctr_lvl5(:,:), ssclr_lvl5(:), &
       svctr_ia(:,:,:), svctr_ib(:,:,:), svctr_s(:,:,:), &
       ! Mass budget arrays
       massbdg(:), scs_rm(:,:,:)

  ! SALSA cloud, precipitation, ice and snow masks
  LOGICAL, allocatable :: cloudmask(:,:,:), drizzmask(:,:,:), icemask(:,:,:), snowmask(:,:,:)

  public :: sflg, ssam_intvl, savg_intvl, statistics, init_stat, write_ps,   &
       acc_tend, updtst, sfc_stat, close_stat, fill_scalar, tke_sgs, sgsflxs,&
       sgs_vel, comp_tke, get_zi, acc_removal, cs_rem_set, acc_massbudged, write_massbudged, mcflg, csflg, &
       salsa_b_bins, cloudy_col_stats

CONTAINS
   !
   ! ---------------------------------------------------------------------
   ! INIT_STAT:  This routine initializes the statistical arrays which
   ! are user/problem defined.  Note that svctr is given 100 elements, and
   ! elements 90 and above are used for computing the TKE budget. Hence
   ! if (nvar2 >= 90 the program stops
   !
   SUBROUTINE init_stat(time, filprf, expnme, nzp)

      USE grid, ONLY          : nxp, nyp, iradtyp
      USE mpi_interface, ONLY : myid, ver, author, info
      USE mo_submctl, ONLY    : nprc, nsnw, fn2a,fn2b,fca,fcb,fra,fia,fib,fsa,spec

      CHARACTER (len=200), INTENT (in) :: filprf, expnme
      INTEGER, INTENT (in)            :: nzp
      REAL, INTENT (in)               :: time

      INTEGER :: i,e
      CHARACTER (len=200) :: fname
      INTEGER :: nspec,nspec_depr ! This whole thing needs to be rewritten...

      nspec = spec%getNSpec()
      ALLOCATE (wtv_sgs(nzp),wtv_res(nzp),wrl_sgs(nzp))
      ALLOCATE (tke_res(nzp),tke_sgs(nzp),tke0(nzp),thvar(nzp))

    ! Combine name arrays
    i = 1; e = nvar1
    s1Total(i:e) = s1
    i = e + 1; e = e + nv1_lvl4
    s1Total(i:e) = s1_lvl4
    i = e + 1; e = e + nv1_lvl5
    s1Total(i:e) = s1_lvl5
    ! ---
    i = 1; e = nvar2
    s2Total(i:e) = s2
    i = e + 1; e = e + nv2_lvl4
    s2Total(i:e) = s2_lvl4
    i = e + 1; e = e + nv2_lvl5
    s2Total(i:e) = s2_lvl5
    i = e + 1; e = e + nv2saa
    s2Total(i:e) = s2Aeroa
    i = e + 1; e = e + nv2sab
    s2Total(i:e) = s2Aerob
    i = e + 1; e = e + nv2sca
    s2Total(i:e) = s2Clouda
    i = e + 1; e = e + nv2scb
    s2Total(i:e) = s2Cloudb
    i = e + 1; e = e + nv2sp
    s2Total(i:e) = s2Precp
    i = e + 1; e = e + nv2sia
    s2Total(i:e) = s2Icea
    i = e + 1; e = e + nv2sib
    s2Total(i:e) = s2Iceb
    i = e + 1; e = e + nv2ss
    s2Total(i:e) = s2Snow

      wtv_sgs(:) = 0.
      wtv_res(:) = 0.
      wrl_sgs(:) = 0.
      tke_res(:) = 0.
      tke_sgs(:) = 0.
      tke0(:)    = 0.

      select case(level)
      case (0)
         nv1 = 13
         nv2 = 58
      case (1)
         nv1 = 14
         nv2 = 58
      case (2)
         nv1 = 20
         nv2 = 83
         if (iradtyp == 3) nv1=21
      ! For SALSA
      case (4,5)
         nv1 = nvar1
         nv2 = nvar2
      case default
         nv1 = nvar1
         nv2 = nvar2
      end select

      IF ( level < 4 ) THEN
         ALLOCATE ( ssclr(nvar1), svctr(nzp,nvar2) )
         IF (mcflg) &
            ALLOCATE ( massbdg(nv1MB) ) ! Mass budged array; ALL CALCULATIONS ARE NOT IMPLEMENTED FOR LEVEL < 4
       svctr(:,:) = 0.
       ssclr(:)   = 0.
       IF (mcflg) massbdg(:) = 0.
       s1bool(1:nvar1) = .TRUE.
       s2bool(1:nvar2) = .TRUE.
    ELSE IF ( level >= 4 ) THEN
       ! Additional arrays for SALSA
       ALLOCATE ( ssclr(nvar1), svctr(nzp,nvar2) )
       ALLOCATE ( ssclr_b(nv1_lvl4), svctr_b(nzp,nv2_lvl4))
       ALLOCATE ( svctr_aa(nzp,fn2a,nv2saa), svctr_ab(nzp,fn2b-fn2a,nv2sab),          &
                  svctr_ca(nzp,fca%cur,nv2sca), svctr_cb(nzp,fcb%cur-fca%cur,nv2scb), &
                  svctr_p(nzp,nprc,nv2sp)   )
       ALLOCATE ( cloudmask(nzp,nxp,nyp), drizzmask(nzp,nxp,nyp), &
                  icemask(nzp,nxp,nyp), snowmask(nzp,nxp,nyp) )
       svctr(:,:) = 0.
       ssclr(:) = 0.
       svctr_b(:,:) = 0.
       ssclr_b(:) = 0.
       svctr_aa(:,:,:) = 0.; svctr_ab(:,:,:) = 0.
       svctr_ca(:,:,:) = 0.; svctr_cb(:,:,:) = 0.
       svctr_p(:,:,:) = 0.
       IF (mcflg) THEN
           ALLOCATE ( massbdg(nv1MB) ) ! Mass budged array
           massbdg(:) = 0.
       END IF
       IF (level >=5 ) THEN
           ALLOCATE( ssclr_lvl5(nv1_lvl5), svctr_lvl5(nzp,nv2_lvl5) )
           ALLOCATE( svctr_ia(nzp,fia%cur,nv2sia), svctr_ib(nzp,fib%cur-fia%cur,nv2sib), &
                  svctr_s(nzp,nsnw,nv2ss) )
           ssclr_lvl5(:) = 0.
           svctr_lvl5(:,:) = 0.
           svctr_ia(:,:,:) = 0.; svctr_ib(:,:,:) = 0.
           svctr_s(:,:,:) = 0.
       END IF

         ! Create a boolean array for items that are actually used
         s2bool(:) = .FALSE.
         s1bool(:) = .FALSE.

         s1bool(1:nvar1) = .TRUE.
         s2bool(1:nvar2) = .TRUE.     ! Original LES vars (assume always used...)

       s1bool(nvar1+1:nvar1+12) = .TRUE. ! Aerosol and cloud number concentrations and sizes (total, a and b)
       s2bool(nvar2+1:nvar2+15) = .TRUE. ! Bin dimensions, number concentrations and radius
       IF (level>=5) THEN
          s1bool(nvar1+nv1_lvl4+1:nvar1+nv1_lvl4+17) = .TRUE.
          s2bool(nvar2+nv2_lvl4+1:nvar2+nv2_lvl4+9) = .TRUE.
          s2bool(nvar2+nv2_lvl4+26:nvar2+nv2_lvl4+29) = .TRUE.
       ENDIF
       ! Bin number concentrations
       i = nvar2+nv2_lvl4+nv2_lvl5+1  ! binned
       s2bool(i) = lbinprof
       i = i+nv2saa
       s2bool(i) = lbinprof .AND. salsa_b_bins
       i = i+nv2sab
       s2bool(i) = lbinprof
       i = i+nv2sca
       s2bool(i) = lbinprof .AND. salsa_b_bins
       i =i+nv2scb
       s2bool(i) = lbinprof
       IF (level>=5) THEN
            i = i+nv2sp
            s2bool(i) = lbinprof
            i = i + nv2sia
            s2bool(i) = lbinprof .AND. salsa_b_bins
            i = i + nv2sib
            s2bool(i) = lbinprof
         END IF

       nspec_depr = 0
       DO e=1,8 ! With water
          IF ( .NOT. spec%isUsed(zspec(e)) ) CYCLE

          ! List of active species (including water, which is the last species)
          nspec_depr = nspec_depr+1
          actspec(nspec_depr)=zspec(e)

          ! Mixing rations in cloud droplets and interstitial aerosol (level 4)
          i=nvar1+13+(e-1)*2
          s1bool(i:i+1)=.TRUE.
          ! Removal with aerosol, cloud and precipitation (level 4)
          i=nvar1+29+(e-1)*3
          s1bool(i:i+2)=.TRUE.

          IF (level>=5) THEN
             ! Mixing ratios in ice and snow
             i=nvar1+nv1_lvl4+18+(e-1)*2
             s1bool(i:i+1)=.TRUE.

             ! Removal with ice and snow
             i=nvar1+nv1_lvl4+34+(e-1)*2
             s1bool(i:i+1)=.TRUE.
          ENDIF


          ! Bulk mixing ratios
          i = nvar2+16+(e-1)*3
          s2bool(i:i+2) = .TRUE.
          IF (level>=5) THEN
             i = nvar2+nv2_lvl4+10+(e-1)*2
             s2bool(i:i+1) = .TRUE.
          ENDIF
          ! Bin number concentrations
          IF (e==8) CYCLE ! Not for water
          i = nvar2+nv2_lvl4+nv2_lvl5+1+e
          s2bool(i) = lbinprof
          i = i+nv2saa
          s2bool(i) = lbinprof .AND. salsa_b_bins
          i = i+nv2sab
          s2bool(i) = lbinprof
          i = i+nv2sca
          s2bool(i) = lbinprof .AND. salsa_b_bins
          i = i+nv2scb
          s2bool(i) = lbinprof
          IF (level>=5) THEN
             i = i+nv2sp
             s2bool(i) = lbinprof
             i = i + nv2sia
             s2bool(i) = lbinprof .AND. salsa_b_bins
             i = i + nv2sib
             s2bool(i) = lbinprof
          ENDIF
       ENDDO

       s1bool(nvar1+53)=.TRUE.  ! Maximum supersaturation

       s2bool(nvar2+40:nvar2+43) = .TRUE.     ! Water mixing ratios and RH

       ! b-bins are not always saved
       IF (.not. salsa_b_bins) THEN
          ! Concentrations and sizes
          s1bool(nvar1+3:nvar1+6) = .FALSE.
          s1bool(nvar1+9:nvar1+12) = .FALSE.
          s1bool(nvar1+nv1_lvl4+3:nvar1+nv1_lvl4+6) = .FALSE.
          s2bool(nvar2+7) = .FALSE.
          s2bool(nvar2+9) = .FALSE.
          s2bool(nvar2+12) = .FALSE.
          s2bool(nvar2+14) = .FALSE.
          s2bool(nvar2+nv2_lvl4+5) = .FALSE.
          s2bool(nvar2+nv2_lvl4+8) = .FALSE.
       ENDIF

       IF (csflg .and. level>3) THEN
           ! Allocate array for level 4 removal rate column statistics
           ! Total number of outputs is 3 for warm (aerosol, cloud and precipitation)
           ! and 5 (add ice and snow) for each species including water
           IF (level==4) THEN
              ALLOCATE( scs_rm(3*nspec,nxp,nyp) )
           ELSE
              ALLOCATE( scs_rm(5*nspec,nxp,nyp) )
          ENDIF
          scs_rm=0. ! Set to zero (during spinup)
       END IF
    END IF ! If level >=4

    fname =  trim(filprf)//'.ts'
    if(myid == 0) print                                                  &
         "(//' ',49('-')/,' ',/,'  Initializing: ',A20)",trim(fname)
      CALL open_nc( fname, expnme, time, (nxp-4)*(nyp-4), ncid1, nrec1, ver, author, info)
      ! Juha: Modified for SALSA output
      CALL define_nc( ncid1, nrec1, COUNT(s1bool), PACK(s1Total,s1bool))
      IF (myid == 0) PRINT *, '   ...starting record: ', nrec1

      fname = trim(filprf)//'.ps'
      IF(myid == 0) PRINT                                                  &
         "(//' ',49('-')/,' ',/,'  Initializing: ',A20)",trim(fname)
      CALL open_nc( fname, expnme, time,(nxp-4)*(nyp-4), ncid2, nrec2, ver, author, info)
      ! Juha: Modified due to SALSA output
      CALL define_nc( ncid2, nrec2, COUNT(s2bool), PACK(s2Total,s2bool), n1=nzp, inae_a=fn2a, inae_b=fn2b-fn2a, &
                      incld_a=fca%cur, incld_b=fcb%cur-fca%cur, inice_a=fia%cur, inice_b=fib%cur-fia%cur, inprc=fra, insnw=fsa)
      IF (myid == 0) PRINT *, '   ...starting record: ', nrec2

      ! Optional column statistics
      IF (csflg) THEN
         fname =  TRIM(filprf)//'.cs'
         IF(myid == 0) PRINT "(//' ',49('-')/,' ',/,'  Initializing: ',A20)",trim(fname)
         CALL open_nc( fname, expnme, time,(nxp-4)*(nyp-4), ncid3, nrec3, ver, author, info)
         IF (ncid3>=0) CALL define_nc_cs(ncid3, nrec3, nxp-4, nyp-4, level, iradtyp, spec%names, nspec)
         IF (myid == 0) PRINT *, '   ...starting record: ', nrec3
      END IF

   END SUBROUTINE init_stat
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
   SUBROUTINE statistics(time)

      USE grid, ONLY : a_up, a_vp, a_wp, a_rc, a_theta, a_rv,                                     &
                       a_rp, a_tp, a_press, nxp, nyp, nzp, dzm, dzt, zm, zt, th00, umean,    &
                       vmean, dn0, precip, a_rpp, a_npp, CCN, iradtyp, a_rflx, a_sflx,       &
                       a_fus, a_fds, a_fuir, a_fdir, albedo, a_srp, a_snrp, a_ncloudp, xt, yt, a_ri, a_nicep, a_srs, a_snrs, snowin

      REAL, INTENT (in) :: time

      REAL :: rxt(nzp,nxp,nyp), rxl(nzp,nxp,nyp), rxv(nzp,nxp,nyp), rnt(nzp,nxp,nyp)
      REAL :: xrpp(nzp,nxp,nyp), xnpp(nzp,nxp,nyp)

    SELECT CASE(level)
       CASE(1,2,3)
          rxt = a_rp ! Total water (vapor + condensed water and ice) = q
          rxl = a_rc-a_rpp ! Cloud water (+aerosol), but no precipitation or ice
          rxv = a_rv ! Water vapor
          xrpp = a_rpp ! Rain water
          xnpp = a_npp ! Rain number
       CASE(4,5)
          rxt = a_rp + a_rc + a_srp + a_ri + a_srs
          rxl = a_rc
          rxv = a_rp
          xrpp = a_srp
          xnpp = a_snrp
    END SELECT

      IF (nsmp == 0.) fsttm = time
      nsmp = nsmp+1.
      ssclr(14:nvar1) = -999.
      !
      ! profile statistics
      !
      CALL accum_stat(nzp, nxp, nyp, a_up, a_vp, a_wp, a_theta, a_press, umean, &
                      vmean)
      IF (iradtyp == 3) THEN
         CALL accum_rad(nzp, nxp, nyp, a_rflx, sflx=a_sflx, sup=a_fus, sdwn=a_fds, &
                        irup=a_fuir, irdwn=a_fdir, alb=albedo)
      ELSE IF (iradtyp > 0) THEN
         CALL accum_rad(nzp, nxp, nyp, a_rflx)
      END IF
      IF (level >= 1) CALL accum_lvl1(nzp, nxp, nyp, rxt)
      IF (level >= 2) CALL accum_lvl2(nzp, nxp, nyp, th00, a_wp,        &
                                     a_theta, a_tp, rxv, rxl, rxt   )
      IF (level >= 3) CALL accum_lvl3(nzp, nxp, nyp, dn0, zm, rxl, xrpp,  &
                                     xnpp, precip, CCN                    )
      IF (level >= 4) CALL accum_lvl4(nzp, nxp, nyp)
      IF (level >= 5) CALL accum_lvl5(nzp,nxp,nyp)
      !
      ! scalar statistics
      !
      CALL set_ts(nzp, nxp, nyp, a_wp, a_theta, dn0, zt,zm,dzt,dzm,th00,time)
      IF ( level >= 1 ) CALL ts_lvl1(nzp, nxp, nyp, dn0, zt, dzm, rxt)
      IF ( level >= 2 ) CALL ts_lvl2(nzp, nxp, nyp, dn0, zm, zt, rxl)
      IF ( level >= 4 ) CALL ts_lvl4(nzp, nxp, nyp)
      IF ( level >= 5 ) CALL ts_lvl5(nzp, nxp, nyp)
      
      CALL write_ts
      
      !
      ! Column statistics
      !
      IF (csflg) THEN
         ! Radiation
         IF (iradtyp == 3) CALL set_cs_any(nxp,nyp,albedo,'albedo')
         
         ! Deposition statistics
         IF (level == 3) THEN
            CALL set_cs_any(nxp,nyp,precip(2,:,:),'prcp')
         ELSE IF (level > 3) THEN
            ! Surface deposition fluxes
            CALL cs_rem_save(nxp,nyp)
         END IF

         ! Warm cloud statistics
         rxt = 0.    ! Condensate
         rnt = 0.
         xrpp = 0.   ! Precipitate
         xnpp = 0.
         IF (level == 1) THEN
             ! No clouds or precipitation
         ELSE IF (level == 2) THEN
            ! Clouds available
            rxt = a_rc
            rnt = CCN
         ELSE IF (level == 3) THEN
            ! Clouds and precipitation available
            rxt = a_rc
            rnt = CCN
            xrpp = a_rpp
            xnpp = a_npp
         ELSE IF (level >= 4) THEN
            ! Levels 4 warm clouds
            rxt = a_rc
            rnt = SUM(a_ncloudp,DIM=4)
            xrpp = a_srp
            xnpp = a_snrp
         END IF
         CALL set_cs_warm(nzp,nxp,nyp,rxt,rnt,xrpp,xnpp,a_theta,dn0,zm,zt,dzm,xt,yt,time)

         ! Ice cloud statistics
         IF (level==5) THEN
            rxt = a_ri
            rnt = SUM(a_nicep,DIM=4)
            xrpp = a_srs
            xnpp = a_snrs
            CALL set_cs_cold(nzp,nxp,nyp,rxt,rnt,xrpp,xnpp,dn0,zm,xt,yt,time)
         END IF

         nrec3 = nrec3 + 1
      END IF

   END SUBROUTINE statistics
   !
   ! -----------------------------------------------------------------------
   ! Subroutines set_cs_warm, cs_rem_set, cs_rem_save and set_cs_any:
   ! write (and compute) column average statistics
   !
   ! Save named data (already available)
   SUBROUTINE set_cs_any(n2,n3,r,nam)
      USE netcdf
      INTEGER, INTENT(in) :: n2,n3
      REAL, INTENT(IN)    :: r(n2,n3)
      CHARACTER (LEN=*) :: nam
      INTEGER :: iret, VarID

      IF (csflg) THEN
         iret = nf90_inq_varid(ncid3,trim(nam),VarID)
         IF (iret == NF90_NOERR) THEN
            iret = nf90_put_var(ncid3, VarID, r(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
            iret = nf90_sync(ncid3)
         END IF
      END IF

    END SUBROUTINE set_cs_any
    !
    ! Removal statistics (level>3): calculate values for further use
    SUBROUTINE cs_rem_set(n2,n3,n4,raer,rcld,rprc,rice,rsnw)
      
      USE mo_submctl, ONLY : nbins, ncld, nprc, nice,  nsnw
      IMPLICIT NONE
      
      INTEGER, INTENT(in) :: n2,n3,n4   ! Grid dimensions
      REAL, INTENT(in) :: raer(n2,n3,n4*nbins), & ! Removal arrays
                          rcld(n2,n3,n4*ncld),  &
                          rprc(n2,n3,n4*nprc),  &
                          rice(n2,n3,n4*nice),  &
                          rsnw(n2,n3,n4*nsnw)

      INTEGER :: si, sii, i, end,str
      
      IF (.NOT. csflg) RETURN
      
      ! Calculate all removal fluxes and save those to scs_rm for later use
      i = 1
      DO si = 1,nspec_depr
         ! Removal by sedimentation of aerosol
         ! Get the real mass index
         sii = spec%getIndex(actspec(si),notFoundValue=0)
         IF (sii == 0) STOP "stat cs_rem_set"
         
         str = getMassIndex(nbins,1,sii)
         end = getMassIndex(nbins,nbins,sii)
         scs_rm(i,:,:) = SUM(raer(:,:,str:end),DIM=3)
         i = i+1

         ! Removal by sedimentation of cloud droplets
         str = getMassIndex(ncld,1,sii) 
         end = getMassIndex(ncld,ncld,sii)
         scs_rm(i,:,:) = SUM(rcld(:,:,str:end),DIM=3)
         i = i+1

         ! Removal by precipitation
         str = getMassIndex(nprc,1,sii)
         end = getMassIndex(nprc,nprc,sii)
         scs_rm(i,:,:) = SUM(rprc(:,:,str:end),DIM=3)
         i = i+1

         IF (level > 4) THEN
            ! Removal by sedimentation of ice particles
            str = getMassIndex(nice,1,sii)
            end = getMassIndex(nice,nice,sii)
            scs_rm(i,:,:) = SUM(rice(:,:,str:end),DIM=3)
            i = i+1

            ! Removal by snow
            str = getMassIndex(nsnw,1,sii)
            end = getMassIndex(nsnw,nsnw,sii)
            scs_rm(i,:,:) = SUM(rsnw(:,:,str:end),DIM=3)
            i = i+1
         END IF
      END DO

   END SUBROUTINE cs_rem_set
   !
   ! Removal statistics (level>3): save values
   SUBROUTINE cs_rem_save(n2,n3)

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n2,n3

      INTEGER :: si, i
      CHARACTER(LEN=3) :: nam

      IF (.NOT. csflg) RETURN

      ! Save all previously calculated removal fluxes
      !   Note: fluxes not calculated during spinup, so saving zeros
      i=1
      DO si = 1,nspec_depr
         nam=spec%names(si)
         
         ! Removal by sedimentation of aerosol
         CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//trim(nam)//'dr') ! 'dr' should be for aerosol and 'ae' for water
         i = i+1
         
         ! Removal by sedimentation of cloud droplets
         CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//trim(nam)//'cl')
         i = i+1
         
         ! Removal by precipitation
         CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//trim(nam)//'pr')
         i = i+1

         IF (level > 4) THEN
            ! Removal by sedimentation of ice particles
            CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//TRIM(nam)//'ic')
            i=i+1

            ! Removal by snow
            CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//TRIM(nam)//'sn')
            i=i+1
        ENDIF
    ENDDO

  END SUBROUTINE cs_rem_save
  !
  ! Calculate warm cloud statistics
  subroutine set_cs_warm(n1,n2,n3,rc,nc,rp,np,th,dn0,zm,zt,dzm,xt,yt,time)

    use netcdf

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: rc(n1,n2,n3),nc(n1,n2,n3),rp(n1,n2,n3),np(n1,n2,n3),th(n1,n2,n3)
    real, intent(in)    :: dn0(n1),zm(n1),zt(n1),dzm(n1),xt(n2),yt(n3),time
    REAL :: lwp(n2,n3), ncld(n2,n3), rwp(n2,n3), nrain(n2,n3), zb(n2,n3), zc(n2,n3), &
                th1(n2,n3), lmax(n2,n3)
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
    do j=3,n3-2
       do i=3,n2-2
          cld=0.
          rn=0.
          sval = 0.
          do k=2,n1
             IF (rc(k,i,j)>0.01e-3) THEN
                ! Cloudy grid
                lwp(i,j)=lwp(i,j)+rc(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                ! Volume weighted average of the CDNC
                ncld(i,j)=ncld(i,j)+nc(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                cld=cld+dn0(k)*(zm(k)-zm(k-1))
                ! Number of cloudy pixels
                ncloudy(i,j)=ncloudy(i,j)+1
                ! Cloud base and top
                zb(i,j)=min(zt(k),zb(i,j))
                zc(i,j)=max(zt(k),zc(i,j))
             END IF
             if (rp(k,i,j) > 0.001e-3) then
                ! Rainy grid cell
                rwp(i,j)=rwp(i,j)+rp(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                ! Volume weighted average of the RDNC
                nrain(i,j)=nrain(i,j)+np(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                rn=rn+dn0(k)*(zm(k)-zm(k-1))
                ! Number of rainy pixels
                nrainy(i,j)=nrainy(i,j)+1
             end if
             ! Maximum liquid water mixing ratio
             lmax(i,j) = max(lmax(i,j),rc(k,i,j))
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
    iret = nf90_inq_varid(ncid3,'time',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, time, start=(/nrec3/))

    if (nrec3 == 1) then
       iret = nf90_inq_varid(ncid3, 'xt', VarID)
       iret = nf90_put_var(ncid3, VarID, xt(3:n2-2), start = (/nrec3/))
       iret = nf90_inq_varid(ncid3, 'yt', VarID)
       iret = nf90_put_var(ncid3, VarID, yt(3:n3-2), start = (/nrec3/))
    END IF

      iret = nf90_inq_varid(ncid3,'lwp',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, lwp(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'rwp',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, rwp(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'Nc',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, ncld(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'Nr',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nrain(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'nccnt',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, ncloudy(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'nrcnt',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nrainy(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'zb',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, zb(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'zc',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, zc(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'zi1',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, th1(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'lmax',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, lmax(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_sync(ncid3)

   END SUBROUTINE set_cs_warm

   !
   ! Calculate cold cloud statistics
   SUBROUTINE set_cs_cold(n1,n2,n3,ri,ni,rs,ns,dn0,zm,xt,yt,time)

     USE netcdf
     USE mo_submctl, only : prlim

     INTEGER, INTENT(in) :: n1,n2,n3
     REAL, INTENT(in)    :: ri(n1,n2,n3),ni(n1,n2,n3),rs(n1,n2,n3),ns(n1,n2,n3)
     REAL, INTENT(in)    :: dn0(n1),zm(n1),xt(n2),yt(n3),time
     REAL :: iwp(n2,n3), nice(n2,n3), swp(n2,n3), nsnow(n2,n3), imax(n2,n3)
     INTEGER :: nicy(n2,n3), nsnowy(n2,n3)
     INTEGER :: i, j, k, iret, VarID
     REAL    :: ice, sn

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
             IF (icemask(k,i,j)) THEN
                ! Icy grid cell
                iwp(i,j)=iwp(i,j)+ri(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                ! Volume weighted average of the ice number concentration
                nice(i,j)=nice(i,j)+ni(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                ice=ice+dn0(k)*(zm(k)-zm(k-1))
                ! Number of icy pixels
                nicy(i,j)=nicy(i,j)+1
             END IF
             if (snowmask(k,i,j)) then
                ! Snowy grid cell
                swp(i,j)=swp(i,j)+rs(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
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

     ! Save the data
     iret = nf90_inq_varid(ncid3,'time',VarID)
     IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, time, start=(/nrec3/)) 

     IF (nrec3 == 1) THEN
        iret = nf90_inq_varid(ncid3, 'xt', VarID)
        iret = nf90_put_var(ncid3, VarID, xt(3:n2-2), start = (/nrec3/))
        iret = nf90_inq_varid(ncid3, 'yt', VarID)
        iret = nf90_put_var(ncid3, VarID, yt(3:n3-2), start = (/nrec3/))
     END IF

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

     iret = nf90_sync(ncid3)

   END SUBROUTINE set_cs_cold
   !
   ! -----------------------------------------------------------------------
   ! Subroutine set_ts: computes and writes time sequence stats
   !
   SUBROUTINE set_ts(n1,n2,n3,w,th,dn0,zt,zm,dzt,dzm,th00,time)

      USE defs, ONLY : cp

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: w(n1,n2,n3),th(n1,n2,n3)
      REAL, INTENT(in)    :: dn0(n1),zt(n1),zm(n1),dzt(n1),dzm(n1),th00,time

      INTEGER :: k
      REAL    :: bf(n1)

      ssclr(1) = time
      ssclr(4) = get_zi(n1, n2, n3, 2, th, dzm, zt, 1.)   ! maximum gradient
      ssclr(5) = get_zi(n1, n2, n3, 3, th, thvar, zt, 1.) ! maximum variance
      !
      ! buoyancy flux statistics
      !
      ssclr(7) = 0.
      DO k = 2, n1-2
         bf(k) = wtv_res(k) + wtv_sgs(k)
         ssclr(7) = ssclr(7) + (tke_res(k)+tke_sgs(k))*dn0(k)/dzt(k)
         svctr(k,33) = svctr(k,33) + wtv_sgs(k)*9.8/th00
      END DO
      ssclr(6) = get_zi(n1, n2, n3, 4, th, bf, zm, 1.) ! minimum buoyancy flux

      ssclr(8) = bf(2)
      ssclr(9) = maxval(w)

      ssclr(12) = ssclr(12)*cp*(dn0(1)+dn0(2))*0.5

   END SUBROUTINE set_ts
   !
   ! -----------------------------------------------------------------------
   ! Subroutine ts_lvl1: computes and writes time sequence stats; for the
   ! zi calculation setting itype=1 selects a concentration threshold
   !
   SUBROUTINE ts_lvl1(n1,n2,n3,dn0,zt,dzm,q)

      USE defs, ONLY : alvl

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: q(n1,n2,n3)
      REAL, INTENT(in)    :: dn0(n1),zt(n1),dzm(n1)

      ssclr(13) = ssclr(13)*alvl*(dn0(1)+dn0(2))*0.5
      ssclr(14) = get_zi(n1, n2, n3, 2, q, dzm, zt, 0.5e-3)

   END SUBROUTINE ts_lvl1
   !
   ! -----------------------------------------------------------------------
   ! Subroutine ts_lvl2: computes and writes time sequence stats
   !
   SUBROUTINE ts_lvl2(n1,n2,n3,dn0,zm,zt,rc)

    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(in) :: dn0(n1), zm(n1), zt(n1), rc(n1,n2,n3)

      INTEGER :: k,i,j
      REAL    :: cpnt, unit, scr(n2,n3)

      ssclr(18) = zt(n1)
      ssclr(19) = 0.
      ssclr(20) = 0.
      ssclr(28) = 0.

      unit = 1./REAL((n2-4)*(n3-4))
      scr(:,:) = 0. ! LWP
      DO j = 3, n3-2
         DO i = 3, n2-2
            cpnt = 0.
            DO k = 2, n1-2
               IF (rc(k,i,j) > 1.e-5) THEN
                  ssclr(17) = max(ssclr(17),zt(k))
                  ssclr(18) = min(ssclr(18),zt(k))
                  cpnt = unit
                  ssclr(20) = max(ssclr(20), rc(k,i,j))
                  ssclr(28) = ssclr(28) + 1.
               END IF
               scr(i,j)=scr(i,j)+rc(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
            END DO
            ssclr(19) = ssclr(19) + cpnt
         END DO
      END DO

    ! liquid water path (without precipitation)
    ssclr(15) = get_avg2dh(n2,n3,scr)
    scr(:,:)=(scr(:,:)-ssclr(15))**2 ! For LWP variance
    ssclr(16) = get_avg2dh(n2,n3,scr)

    if (ssclr(18) == zt(n1)) ssclr(18) = -999.

   END SUBROUTINE ts_lvl2
   !
   ! -----------------------------------------------------------------------
   ! Subroutine ts_lvl4: computes and writes time sequence stats of Salsa variables --
   !  Implemented by Zubair Maalick 20/07/2015
   !  Some rewriting and adjusting by Juha Tonttila
   !
   SUBROUTINE ts_lvl4(n1,n2,n3)
      USE mo_submctl, ONLY : nlim,spec
      USE grid, ONLY       : bulkNumc, bulkMixrat, meanradius,dzt,a_rh

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1,n2,n3

      REAL    :: a0(n1,n2,n3), a1(n1,n2,n3), a2(n1,n2,n3)
      INTEGER :: ii,ss
      LOGICAL :: mask(n1,n2,n3)

    ! Clouds; combined and a and b separately for a and b regions, respectively
    CALL bulkNumc('cloud','a',a0)
    CALL bulkNumc('cloud','b',a1)
    ssclr_b(1) = get_avg_ts(n1,n2,n3,a0+a1,dzt,cloudmask)
    CALL meanRadius('cloud','ab',a2)
    ssclr_b(2) = get_avg_ts(n1,n2,n3,a2,dzt,cloudmask)

    mask(:,:,:) = ( a0(:,:,:) > nlim .AND. cloudmask(:,:,:) )
    ssclr_b(3) = get_avg_ts(n1,n2,n3,a0,dzt,mask)
    CALL meanRadius('cloud','a',a2)
    ssclr_b(4) = get_avg_ts(n1,n2,n3,a2,dzt,mask)

    mask(:,:,:) = ( a1(:,:,:) > nlim .AND. cloudmask(:,:,:) )
    ssclr_b(5) = get_avg_ts(n1,n2,n3,a1,dzt,mask)
    CALL meanRadius('cloud','b',a2)
    ssclr_b(6) = get_avg_ts(n1,n2,n3,a2,dzt,mask)

    ! Interstitial aerosol; combined and a and b separately
    CALL bulkNumc('aerosol','a',a0)
    CALL bulkNumc('aerosol','b',a1)
    ssclr_b(7) = get_avg_ts(n1,n2,n3,a0+a1,dzt,cloudmask)
    CALL meanRadius('aerosol','ab',a2)
    ssclr_b(8) = get_avg_ts(n1,n2,n3,a2,dzt,cloudmask)

    ssclr_b(9) = get_avg_ts(n1,n2,n3,a0,dzt,cloudmask)
    CALL meanRadius('aerosol','a',a2)
    ssclr_b(10) = get_avg_ts(n1,n2,n3,a2,dzt,cloudmask)

    ssclr_b(11) = get_avg_ts(n1,n2,n3,a1,dzt,cloudmask)
    CALL meanRadius('aerosol','b',a2)
    ssclr_b(12) = get_avg_ts(n1,n2,n3,a2,dzt,cloudmask)

    ! Species mixing ratios in cloud droplets and interstitial aerosol
    ii = 13
    DO ss = 1,8
       IF (spec%isUsed(zspec(ss))) THEN
          CALL bulkMixrat(zspec(ss),'cloud','ab',a0)
          ssclr_b(ii) = get_avg_ts(n1,n2,n3,a0,dzt,cloudmask)
          ii = ii + 1
          
          CALL bulkMixrat(zspec(ss),'aerosol','ab',a0)
          ssclr_b(ii) = get_avg_ts(n1,n2,n3,a0,dzt,cloudmask)
          ii = ii + 1
       ELSE
          ii = ii + 2
       END IF
    END DO
    
    ! Removal statistics elsewhere ..

    ! Maximum supersaturation
    ssclr_b(53) = (MAXVAL(a_rh(2:n1,3:n2-2,3:n3-2))-1.0)*100.

  END SUBROUTINE ts_lvl4
  !
  ! -----------------------------------------------------------------------
  ! subroutine ts_lvl5: computes and writes time sequence stats of Salsa variables --
  !  Implemented by Jaakko Ahola 15/12/2016
  !
  SUBROUTINE ts_lvl5(n1,n2,n3)
    USE mo_submctl, only : nlim,prlim
    USE grid, ONLY : bulkNumc, bulkMixrat,meanRadius,dzt, &
        dn0, zm, a_ri, a_srs, snowin, a_rhi

    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3

    REAL :: a0(n1,n2,n3), a1(n1,n2,n3), a2(n1,n2,n3)
    REAL :: scr(n2,n3), scr2(n2,n3)
    REAL :: sscnt
    INTEGER :: i, j, k, ii,ss
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
    ! Maximum supersaturation over ice
    ssclr_lvl5(17) = (MAXVAL(a_rhi(2:n1,3:n2-2,3:n3-2))-1.0)*100.

    ! Species mixing ratios in ice and snow
    ii = 18
    DO ss = 1,8
       IF ( spec%isUsed(zspec(ss)) ) THEN
          CALL bulkMixrat(zspec(ss),'ice','ab',a0)
          ssclr_lvl5(ii) = get_avg_ts(n1,n2,n3,a0,dzt,icemask)
          ii = ii + 1

          CALL bulkMixrat(zspec(ss),'snow','ab',a0)
          ssclr_lvl5(ii) = get_avg_ts(n1,n2,n3,a0,dzt,snowmask)
          ii = ii + 1
       ELSE
          ii = ii + 2
       END IF
    END DO

    ! Removal statistics elsewhere ..

  END SUBROUTINE ts_lvl5
   !
   !---------------------------------------------------------------------
   ! Subroutine ACCUM_STAT: Accumulates various statistics over an
   ! averaging period for base (level 0) version of model
   !
   SUBROUTINE accum_stat(n1,n2,n3,u,v,w,t,p,um,vm)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, DIMENSION (n1,n2,n3), INTENT (in) :: u, v, w, t, p
      REAL, INTENT (in) :: um, vm

      INTEGER :: k
      REAL    :: a1(n1), b1(n1), c1(n1), d1(n1), a3(n1), b3(n1), tmp(n1)

      CALL get_avg3(n1,n2,n3, u,a1)
      CALL get_avg3(n1,n2,n3, v,b1)
      CALL get_avg3(n1,n2,n3, t,c1)
      CALL get_avg3(n1,n2,n3, p,d1)
      CALL get_var3(n1,n2,n3, t, c1, thvar)
      CALL get_3rd3(n1,n2,n3, t, c1, b3) ! Used to be (t-a1)**3
      tmp(:) = 0.
      CALL get_3rd3(n1,n2,n3, w, tmp, a3) ! Now just w**3

      DO k = 1, n1
         svctr(k,10) = svctr(k,10) + a1(k) + um
         svctr(k,11) = svctr(k,11) + b1(k) + vm
         svctr(k,12) = svctr(k,12) + c1(k)
         svctr(k,13) = svctr(k,13) + d1(k)
         svctr(k,17) = svctr(k,17) + thvar(k)
         svctr(k,18) = svctr(k,18) + a3(k)
         svctr(k,19) = svctr(k,19) + b3(k)
      END DO

   END SUBROUTINE accum_stat
   !
   !---------------------------------------------------------------------
   ! Subroutine ACCUM_RAD: Accumulates various statistics over an
   ! averaging period for radiation variables
   !
   SUBROUTINE accum_rad(n1,n2,n3,rflx,sflx,sup,sdwn,irup,irdwn,alb)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: rflx(n1,n2,n3)
      REAL, OPTIONAL, INTENT (in) :: sflx(n1,n2,n3), alb(n2,n3)
      REAL, OPTIONAL, INTENT (in) :: sup(n1,n2,n3), sdwn(n1,n2,n3), irup(n1,n2,n3), irdwn(n1,n2,n3)

      INTEGER :: k
      REAL    :: a1(n1),a2(n1)

      CALL get_avg3(n1,n2,n3,rflx,a1)
      CALL get_var3(n1,n2,n3,rflx,a1,a2)
      DO k = 1, n1
         svctr(k,55) = svctr(k,55) + a1(k)
         svctr(k,56) = svctr(k,56) + a2(k)
      END DO

      IF (present(sflx)) THEN
         CALL get_avg3(n1,n2,n3,sflx,a1)
         CALL get_var3(n1,n2,n3,sflx,a1,a2)
         DO k = 1, n1
            svctr(k,57) = svctr(k,57) + a1(k)
            svctr(k,58) = svctr(k,58) + a2(k)
         END DO
         ssclr(21) = get_avg2dh(n2,n3,alb)
      END IF

      IF (present(sup)) THEN
        CALL get_avg3(n1,n2,n3,sup,a1)
        svctr(:,93)=svctr(:,93) + a1(:)
        CALL get_avg3(n1,n2,n3,sdwn,a1)
        svctr(:,94)=svctr(:,94) + a1(:)
        CALL get_avg3(n1,n2,n3,irup,a1)
        svctr(:,95)=svctr(:,95) + a1(:)
        CALL get_avg3(n1,n2,n3,irdwn,a1)
        svctr(:,96)=svctr(:,96) + a1(:)
      END IF

   END SUBROUTINE accum_rad
   !
   !---------------------------------------------------------------------
   ! Subroutine ACCUM_LVL1: Accumulates various statistics over an
   ! averaging period for moisture variable (smoke or total water)
   !
   SUBROUTINE accum_lvl1(n1,n2,n3,q)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: q(n1,n2,n3)

      INTEGER :: k
      REAL    :: a1(n1),a2(n1),a3(n1)

      CALL get_avg3(n1,n2,n3,q,a1)
      CALL get_var3(n1,n2,n3,q,a1,a2)
      CALL get_3rd3(n1,n2,n3,q,a1,a3)

      DO k = 1, n1
         svctr(k,50) = svctr(k,50) + a1(k)
         svctr(k,51) = svctr(k,51) + a2(k)
         svctr(k,52) = svctr(k,52) + a3(k)
      END DO

   END SUBROUTINE accum_lvl1
   !
   !---------------------------------------------------------------------
   ! Subroutine ACCUM_LVL2: Accumulates specialized statistics that depend
   ! on level 2 variables.
   !
   SUBROUTINE accum_lvl2(n1, n2, n3, th00, w, th, t, &
                         rv, rl, rt)

      USE defs, ONLY : ep2

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)                      :: th00
      REAL, INTENT (in), DIMENSION(n1,n2,n3) :: w, th, t, rv, rl, rt

      REAL, DIMENSION(n1,n2,n3) :: tv    ! Local variable
      INTEGER                   :: k, i, j, kp1
      REAL, DIMENSION(n1)       :: a1, a2, a3, tvbar
      REAL, DIMENSION(n1,n2,n3) :: scr, xy1, xy2, tw, tvw, rtw
      REAL, DIMENSION(n2,n3) :: scr2d
      LOGICAL :: cond(n1,n2,n3)

      !
      ! liquid water statistics
      !
      CALL get_avg3(n1,n2,n3,rl,a1)
      CALL get_var3(n1,n2,n3,rl,a1,a2)
      CALL get_3rd3(n1,n2,n3,rl,a1,a3)
      svctr(:,59) = svctr(:,59) + a1(:)
      svctr(:,60) = svctr(:,60) + a2(:)
      svctr(:,61) = svctr(:,61) + a3(:)

      !
      ! Do some conditional sampling statistics: cloud, cloud-core
      !
      tv(:,:,:) = th(:,:,:)*(1. + ep2*rv(:,:,:) - rl(:,:,:))
      CALL get_avg3(n1,n2,n3,tv,tvbar)
      !
      xy1 = 0.
      xy2 = 0.
      tvw = 0.
      DO k = 1, n1
         kp1 = k+1
         IF (k == n1) kp1 = k
         DO j = 3, n3-2
            DO i = 3, n2-2
               IF (rl(k,i,j) > 1.e-5) THEN
                  xy1(k,i,j) = 1.
                  IF (tv(k,i,j) > tvbar(k)) THEN
                     xy2(k,i,j) = 1.
                  END IF
                  !
                  tw(k,i,j) = (.5*(t(k,i,j)+t(kp1,i,j))+th00)*w(k,i,j)
                  tvw(k,i,j) = (.5*(tv(k,i,j)+tv(kp1,i,j)))*w(k,i,j)
                  rtw(k,i,j) = (.5*(rt(k,i,j)+rt(kp1,i,j)))*w(k,i,j)
               END IF
            END DO
         END DO
      END DO
      CALL get_avg3(n1,n2,n3,xy1,a1)
      svctr(:,64) = svctr(:,64)+a1(:)
      CALL get_avg3(n1,n2,n3,xy1,a1,normalize=.FALSE.)
      svctr(:,65) = svctr(:,65)+a1(:)   ! Counts
      cond(:,:,:) = xy1(:,:,:)>0.5
      CALL get_avg3(n1,n2,n3,w,a1,cond=cond)
      svctr(:,66) = svctr(:,66)+a1(:)
      CALL get_avg3(n1,n2,n3,t+th00,a1,cond=cond)
      svctr(:,67) = svctr(:,67)+a1(:)
      CALL get_avg3(n1,n2,n3,tv,a1,cond=cond)
      svctr(:,68) = svctr(:,68)+a1(:)
      CALL get_avg3(n1,n2,n3,rt,a1,cond=cond)
      svctr(:,69) = svctr(:,69)+a1(:)
      CALL get_avg3(n1,n2,n3,rl,a1,cond=cond)
      svctr(:,70) = svctr(:,70)+a1(:)
      CALL get_avg3(n1,n2,n3,tw,a1,cond=cond)
      svctr(:,71) = svctr(:,71)+a1(:)
      CALL get_avg3(n1,n2,n3,tvw,a1,cond=cond)
      svctr(:,72) = svctr(:,72)+a1(:)
      CALL get_avg3(n1,n2,n3,rtw,a1,cond=cond)
      svctr(:,73) = svctr(:,73)+a1(:)
      CALL get_avg3(n1,n2,n3,xy2,a1)
      svctr(:,74) = svctr(:,74)+a1(:)
      CALL get_avg3(n1,n2,n3,xy2,a1,normalize=.FALSE.)
      svctr(:,75) = svctr(:,75)+a1(:)   ! Counts
      cond(:,:,:) = xy2(:,:,:)>0.5
      CALL get_avg3(n1,n2,n3,w,a1,cond=cond)
      svctr(:,76) = svctr(:,76)+a1(:)
      CALL get_avg3(n1,n2,n3,t+th00,a1,cond=cond)
      svctr(:,77) = svctr(:,77)+a1(:)
      CALL get_avg3(n1,n2,n3,tv,a1,cond=cond)
      svctr(:,78) = svctr(:,78)+a1(:)
      CALL get_avg3(n1,n2,n3,rt,a1,cond=cond)
      svctr(:,79) = svctr(:,79)+a1(:)
      CALL get_avg3(n1,n2,n3,rl,a1,cond=cond)
      svctr(:,80) = svctr(:,80)+a1(:)
      CALL get_avg3(n1,n2,n3,tw,a1,cond=cond)
      svctr(:,81) = svctr(:,81)+a1(:)
      CALL get_avg3(n1,n2,n3,tvw,a1,cond=cond)
      svctr(:,82) = svctr(:,82)+a1(:)
      CALL get_avg3(n1,n2,n3,rtw,a1,cond=cond)
      svctr(:,83) = svctr(:,83)+a1(:)

   END SUBROUTINE accum_lvl2
   !
   !---------------------------------------------------------------------
   ! Subroutine ACCUM_LVL3: Accumulates specialized statistics that depend
   ! on level 3 variables.
   !
   SUBROUTINE accum_lvl3(n1, n2, n3, dn0, zm, rc, rr, nr, rrate, CCN)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)                      :: CCN
      REAL, INTENT (in), DIMENSION(n1)       :: zm, dn0
      REAL, INTENT (in), DIMENSION(n1,n2,n3) :: rc, rr, nr, rrate

      INTEGER                :: k, i, j
      REAL                   :: nrsum, nrcnt, rrcnt, rrcb
      REAL                   :: rmax, rmin
      REAL, DIMENSION(n1)    :: a1
      REAL, DIMENSION(n2,n3) :: scr2
      REAL                   :: mask(n1,n2,n3), tmp(n1,n2,n3)
      LOGICAL :: below

      !
      ! Average rain water mixing ratio
      !
      CALL get_avg3(n1,n2,n3,rr,a1)
      svctr(:,86) = svctr(:,86) + a1(:) ! rr (kg/kg)
      svctr(:,84) = svctr(:,84) + CCN   ! nc (#/kg)

      !
      ! conditionally average rain droplet concentrations
      !
      WHERE (rr > 0.001e-3)
         mask = 1.
      ELSE WHERE
         mask = 0.
      END WHERE

      CALL get_avg3(n1,n2,n3,nr,a1,cond=(mask>0.5))
      svctr(:,85) = svctr(:,85)+a1(:)
      CALL get_avg3(n1,n2,n3,mask,a1)
      svctr(:,91) = svctr(:,91)+a1(:)

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
      ELSE WHERE
         mask = 0.
      END WHERE

      tmp(:,:,:) = mask(:,:,:)*rrate(:,:,:)
      CALL get_avg3(n1,n2,n3,tmp,a1)
      svctr(:,90) = svctr(:,90)+a1(:)
      CALL get_avg3(n1,n2,n3,mask,a1)
      svctr(:,89) = svctr(:,89)+a1(:)

      !
      ! Histogram of surface rain rates
      !
      mask = 0.
      DO k = 1, n1
         rmin = max(6.2e-8,(k-1)*3.421e-5)
         rmax =  k * 3.421e-5
         DO j = 3, n3-2
            DO i = 3, n2-2
               IF (rrate(2,i,j) > rmin .AND. rrate(2,i,j) <= rmax) mask(k,i,j) = 1.
            END DO
         END DO
      END DO
      CALL get_avg3(n1,n2,n3,mask,a1,normalize = .FALSE.)
      svctr(:,92) = svctr(:,92)+a1(:)

      !
      ! Temporal statistics
      !
      scr2(:,:) = 0.
      nrsum = 0.
      nrcnt = 0.
      rrcnt = 0.
      rrcb = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            below = .TRUE.
            DO k=2, n1

               ! RWP
               scr2(i,j) = scr2(i,j)+rr(k,i,j)*dn0(k)*(zm(k)-zm(k-1))

               ! Rainy grid cell
               IF (rr(k,i,j) > 0.001e-3) THEN
                  nrsum = nrsum + nr(k,i,j)
                  nrcnt = nrcnt + 1.
               END IF

              ! Surface precipitation for this column
               IF (k == 2 .AND. rrate(k,i,j) > 3.65e-5) rrcnt = rrcnt + 1.

              ! Precpitation at cloud base (no cloud = no precip.)
               IF (rc(k,i,j) > 1.e-5 .AND. below) THEN
                 ! Take precpitation from level k-1 (>=2), which is just below cloud base
                 rrcb = rrcb + rrate(max(2,k-1),i,j)
                 below = .FALSE.
               END IF
            END DO
         END DO
      END DO
      ssclr(24) = rrcnt/REAL( (n3-4)*(n2-4) )
      ssclr(22) = get_avg2dh(n2,n3,scr2)
      scr2(:,:) = rrate(2,:,:)
      ssclr(23) = get_avg2dh(n2,n3,scr2)
      ssclr(25) = CCN
      IF (nrcnt > 0.) ssclr(26) = nrsum/nrcnt
      ssclr(27) = nrcnt
      ssclr(29) = rrcb/REAL( (n3-4)*(n2-4) )

   END SUBROUTINE accum_lvl3

   !---------------------------------------------------------------------
   ! Subroutine ACCUM_LVL4: Accumulates specialized statistics that depend
   ! on level 4 variables.
   !
   SUBROUTINE accum_lvl4(n1,n2,n3)
      USE mo_submctl, ONLY : in1a,in2b,fn2a,fn2b,     &
                             ica,fca,icb,fcb,ira,fra, &
                             nprc,nlim,prlim,spec
      USE grid, ONLY : bulkNumc, bulkMixrat, meanRadius, binSpecMixrat, &
                       a_rc, a_srp, a_rp, a_rh,    &
                       a_naerop, a_ncloudp, a_nprecpp

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1,n2,n3
      INTEGER :: ii,ss,bb

      REAL, DIMENSION(n1,n2,n3)           :: a1,a12
      REAL, DIMENSION(n1,5)               :: a2
      REAL, DIMENSION(n1,fn2a)            :: a3_a
      REAL, DIMENSION(n1,fn2b-fn2a)       :: a3_b
      REAL, DIMENSION(n1,fca%cur)         :: a4_a
      REAL, DIMENSION(n1,fcb%cur-fca%cur) :: a4_b
      REAL, DIMENSION(n1,nprc)            :: a5

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
    drizzmask(:,:,:) = ( a1(:,:,:) > prlim .AND. a_srp(:,:,:) > 1.e-6 )

      svctr_b(:,6:10) = svctr_b(:,6:10) + a2(:,1:5)

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

      svctr_b(:,11:15) = svctr_b(:,11:15) + a2(:,1:5)

      ! Bin number concentrations
      ! -------------------------------------------
      IF (lbinprof) THEN
         DO bb = in1a, fn2a
            CALL get_avg3(n1,n2,n3,a_naerop(:,:,:,bb),a3_a(:,bb))
         END DO
         DO bb = in2b, fn2b
            CALL get_avg3(n1,n2,n3,a_naerop(:,:,:,bb),a3_b(:,bb-fn2a))
         END DO
         DO bb = ica%cur, fca%cur
            CALL get_avg3(n1,n2,n3,a_ncloudp(:,:,:,bb),a4_a(:,bb))
         END DO
         DO bb = icb%cur, fcb%cur
            CALL get_avg3(n1,n2,n3,a_ncloudp(:,:,:,bb),a4_b(:,bb-fca%cur))
         END DO
         DO bb = ira, fra
            CALL get_avg3(n1,n2,n3,a_nprecpp(:,:,:,bb),a5(:,bb))
         END DO

         svctr_aa(:,:,1) = svctr_aa(:,:,1) + a3_a(:,:)
         svctr_ab(:,:,1) = svctr_ab(:,:,1) + a3_b(:,:)
         svctr_ca(:,:,1) = svctr_ca(:,:,1) + a4_a(:,:)
         svctr_cb(:,:,1) = svctr_cb(:,:,1) + a4_b(:,:)
         svctr_p(:,:,1) = svctr_p(:,:,1) + a5(:,:)
      END IF

      ! Species mixing ratios
      ! -------------------------------------------
      ii = 16
      DO ss = 1, 8  
         IF (spec%isUsed(zspec(ss))) THEN
            ! Total mass mixing ratios
            CALL bulkMixrat(zspec(ss),'aerosol','ab',a1)
            CALL get_avg3(n1,n2,n3,a1,a2(:,1))

            ! In-cloud
            CALL bulkMixrat(zspec(ss),'cloud','ab',a1)
            CALL get_avg3(n1,n2,n3,a1,a2(:,2),cond=cloudmask)

            ! In-drizzle
            CALL bulkMixrat(zspec(ss),'precp','a',a1)
            CALL get_avg3(n1,n2,n3,a1,a2(:,3),cond=drizzmask)

            svctr_b(:,ii:ii+2) = svctr_b(:,ii:ii+2) + a2(:,1:3)

            ! Binned mixing ratios
            IF (lbinprof .AND. ss < 8) THEN
               DO bb = in1a, fn2a
                  CALL binSpecMixrat('aerosol',zspec(ss),bb,a1) ! z,x,y-field for bin bb
                  CALL get_avg3(n1,n2,n3,a1,a3_a(:,bb))         ! average profile for bin bb for species ss
               END DO
               DO bb = in2b, fn2b
                  CALL binSpecMixrat('aerosol',zspec(ss),bb,a1) ! z,x,y-field for bin bb
                  CALL get_avg3(n1,n2,n3,a1,a3_b(:,bb-fn2a))    ! average profile for bin bb for species ss
               END DO

               DO bb = ica%cur, fca%cur
                  CALL binSpecMixrat('cloud',zspec(ss),bb,a1)  ! z,x,y-field for bin bb
                  CaLL get_avg3(n1,n2,n3,a1,a4_a(:,bb),cond=cloudmask)        ! average profile for bin bb for species ss
               END DO
               DO bb = icb%cur, fcb%cur
                  CALL binSpecMixrat('cloud',zspec(ss),bb,a1)  ! z,x,y-field for bin bb
                  CALL get_avg3(n1,n2,n3,a1,a4_b(:,bb-fca%cur),cond=cloudmask)! average profile for bin bb for species ss
               END DO

               DO bb = 1, nprc
                  CALL binSpecMixrat('precp',zspec(ss),bb,a1)  ! z,x,y-field for bin bb
                  CALL get_avg3(n1,n2,n3,a1,a5(:,bb),cond=drizzmask)          ! average profile for bin bb for species ss
               END DO

               svctr_aa(:,:,ss+1) = svctr_aa(:,:,ss+1) + a3_a(:,:)
               svctr_ab(:,:,ss+1) = svctr_ab(:,:,ss+1) + a3_b(:,:)
               svctr_ca(:,:,ss+1) = svctr_ca(:,:,ss+1) + a4_a(:,:)
               svctr_cb(:,:,ss+1) = svctr_cb(:,:,ss+1) + a4_b(:,:)
               svctr_p(:,:,ss+1)  = svctr_p(:,:,ss+1)  + a5(:,:)
            END IF

         END IF ! IsUsed

         ii = ii + 3

      END DO ! ss

      ! Liquid water mixing ratio
      CALL get_avg3(n1,n2,n3,a_rc,a2(:,1))

      ! Precipitation mixing ratio
      CALL get_avg3(n1,n2,n3,a_srp,a2(:,2))

      ! Water vapor mixing ratio
      CALL get_avg3(n1,n2,n3,a_rp,a2(:,3))

      ! Relative humidity
      CALL get_avg3(n1,n2,n3,a_rh,a2(:,4))
      a2(:,4) = a2(:,4)*100.0 ! RH in %

      svctr_b(:,40:43) = svctr_b(:,40:43) + a2(:,1:4)

   END SUBROUTINE accum_lvl4

   !---------------------------------------------------------------------
   ! SUBROUTINE ACCUM_LVL5: Accumulates specialized statistics that depend
   ! on level 5 variables.
   !
   subroutine accum_lvl5(n1,n2,n3)
     USE mo_submctl, only : iia,fia,iib,fib,isa,fsa,nsnw,prlim, spec
     USE grid, ONLY : bulkNumc, bulkMixrat, meanRadius, binSpecMixrat, &
                      a_ri, a_srs, a_rhi, a_nicep, a_nsnowp, snowin

     IMPLICIT NONE

     INTEGER, INTENT(in) :: n1,n2,n3
     INTEGER :: ii,ss,bb

     REAL, DIMENSION(n1,n2,n3) :: a1,a12
     REAL, DIMENSION(n1,3) :: a2
     REAL, DIMENSION(n1,fia%cur) :: a4_a
     REAL, DIMENSION(n1,fib%cur-fia%cur) :: a4_b
     REAL, DIMENSION(n1,nsnw) :: a5

     ! *************************
     ! Bulk output for SALSA
     ! *************************
     CALL bulkNumc('ice','a',a1)
     CALL get_avg3(n1,n2,n3,a1,a2(:,1))

     CALL bulkNumc('ice','b',a12)
     CALL get_avg3(n1,n2,n3,a12,a2(:,2))

    ! Generate ice mask (grid cells with ice)
    !   - One 60 um ice particle would be about 1e-10 kg and concentrations can be as low as 1e-6 #/kg (=prlim)
    icemask(:,:,:) = ( a1(:,:,:)+a12(:,:,:) > prlim .AND. a_ri(:,:,:) > 1.e-15)

     CALL bulkNumc('snow','a',a1)
     CALL get_avg3(n1,n2,n3,a1,a2(:,3))

    ! Generate snow mask (grid cells with snow)
    ! - One small 250 um snow particle would be about 1e-8 kg, but concentrations are often lower than that of ice
    snowmask(:,:,:) = ( a1(:,:,:) > prlim .AND. a_srs(:,:,:) > 1.e-20 )

     svctr_lvl5(:,4:6) = svctr_lvl5(:,4:6) + a2(:,1:3)

     ! Particle radius

     ! Ice
     CALL meanRadius('ice','a',a1)
     CALL get_avg3(n1,n2,n3,a1,a2(:,1),cond=icemask)

     CALL meanRadius('ice','b',a1)
     CALL get_avg3(n1,n2,n3,a1,a2(:,2),cond=icemask)

     ! Snow
     CALL meanRadius('snow','a',a1)
     CALL get_avg3(n1,n2,n3,a1,a2(:,3),cond=snowmask)

     svctr_lvl5(:,7:9) = svctr_lvl5(:,7:9) + a2(:,1:3)

     ! Bin number concentrations
     ! -------------------------------------------
     IF (lbinprof) THEN
        DO bb = iia%cur,fia%cur
           CALL get_avg3(n1,n2,n3,a_nicep(:,:,:,bb),a4_a(:,bb))
        END DO
        DO bb = iib%cur,fib%cur
           CALL get_avg3(n1,n2,n3,a_nicep(:,:,:,bb),a4_b(:,bb-fia%cur))
        END DO
        DO bb = isa,fsa           
           CALL get_avg3(n1,n2,n3,a_nsnowp(:,:,:,bb),a5(:,bb))
        END DO
        
        svctr_ia(:,:,1) = svctr_ia(:,:,1) + a4_a(:,:)
        svctr_ib(:,:,1) = svctr_ib(:,:,1) + a4_b(:,:)
        svctr_s(:,:,1) = svctr_s(:,:,1) + a5(:,:)
     END IF

     ! Species mixing ratios
     ! -------------------------------------------
     ii=10 ! 'P_cSO4i'
     DO ss = 1,8  ! Including water
        IF (spec%isUsed(zspec(ss)) ) THEN
           ! Total mass mixing ratios

           ! In-ice
           CALL bulkMixrat(zspec(ss),'ice','ab',a1)
           CALL get_avg3(n1,n2,n3,a1,a2(:,1),cond=icemask)

           ! In-snow
           CALL bulkMixrat(zspec(ss),'snow','a',a1)
           CALL get_avg3(n1,n2,n3,a1,a2(:,2),cond=snowmask)

           svctr_lvl5(:,ii:ii+1) = svctr_lvl5(:,ii:ii+1) + a2(:,1:2)

           ! Binned mixing ratios
           IF (lbinprof .AND. ss<8) THEN
              DO bb = iia%cur,fia%cur
                 CALL binSpecMixrat('ice',zspec(ss),bb,a1)
                 CaLL get_avg3(n1,n2,n3,a1,a4_a(:,bb),cond=icemask)
              END DO
              DO bb = iib%cur,fib%cur
                 CALL binSpecMixrat('ice',zspec(ss),bb,a1)
                 CALL get_avg3(n1,n2,n3,a1,a4_b(:,bb-fia%cur),cond=icemask)
              END DO

              DO bb = 1,nsnw
                 CALL binSpecMixrat('snow',zspec(ss),bb,a1)
                 CALL get_avg3(n1,n2,n3,a1,a5(:,bb),cond=snowmask)
              END DO

              svctr_ia(:,:,ss+1) = svctr_ia(:,:,ss+1) + a4_a(:,:)
              svctr_ib(:,:,ss+1) = svctr_ib(:,:,ss+1) + a4_b(:,:)
              svctr_s(:,:,ss+1) = svctr_s(:,:,ss+1) + a5(:,:)
           END IF

        END IF ! IsUsed

        ii = ii + 2

     END DO ! ss

     ! Ice water mixing ratio
     CALL get_avg3(n1,n2,n3,a_ri,a2(:,1))

     ! Snow water mixing ratio
     CALL get_avg3(n1,n2,n3,a_srs,a2(:,2))

     ! Relative humidity ove ice
     CALL get_avg3(n1,n2,n3,a_rhi,a2(:,3))
     a2(:,3)=a2(:,3)*100.0 ! RH in %

     svctr_lvl5(:,26:28) = svctr_lvl5(:,26:28) + a2(:,1:3)

     ! Snow deposition flux
     call get_avg3(n1,n2,n3,snowin,a2(:,1))
     svctr_lvl5(:,29)=svctr(:,29)+a2(:,1)

   end subroutine accum_lvl5

   !
   !
   !
   SUBROUTINE comp_tke(n1,n2,n3,dzm,th00,u,v,w,s)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: dzm(n1),th00,u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
      REAL, INTENT (inout) :: s(n1,n2,n3)

      INTEGER :: k,kp1
      REAL    :: x1(n1), x2(n1)

      !
      ! ------
      ! Calculates buoyancy forcing
      !
      CALL get_buoyancy(n1,n2,n3,s,w,th00)
      !
      ! ------
      ! Estimates shear component of TKE budget
      !
      CALL get_shear(n1,n2,n3,u,v,w,dzm)
      !
      ! ------
      ! Calculates horizontal variances and resolved TKE
      !
      CALL get_avg3(n1,n2,n3,u**2,x1)
      CALL get_avg3(n1,n2,n3,u,x2)
      svctr(:,14) = svctr(:,14) + (x1(:)-x2(:)**2) ! Variance
      tke_res(:) = (x1(:)-x2(:)**2)

      CALL get_avg3(n1,n2,n3,v**2,x1)
      CALL get_avg3(n1,n2,n3,v,x2)
      svctr(:,15) = svctr(:,15) + (x1(:)-x2(:)**2)
      tke_res(:)  = tke_res(:) + (x1(:)-x2(:)**2)

      CALL get_avg3(n1,n2,n3,w**2,x1)
      svctr(:,16) = svctr(:,16)+x1(:) ! Raw moment
      DO k = 1, n1
         kp1 = min(k+1,n1)
         tke_res(k)  = 0.5*(0.5*(tke_res(k)+tke_res(kp1)) + x1(k))
      END DO
      IF (nsmp == 0) tke0(:) = tke_res(:)

   END SUBROUTINE comp_tke
   !
   ! ---------------------------------------------------------------------
   ! get_buoyancy:  estimates buoyancy production term in tke budget
   !
   SUBROUTINE get_buoyancy(n1,n2,n3,b,w,th00)

      USE defs, ONLY : g

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: w(n1,n2,n3),th00
      REAL, INTENT(inout) :: b(n1,n2,n3)

      INTEGER :: i,j,k,kp1

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               kp1 = min(k+1,n1)
               b(k,i,j) = (b(k,i,j) + b(kp1,i,j))
            END DO
         END DO
      END DO
      CALL get_cor3(n1,n2,n3,b,w,wtv_res)
      DO k = 1, n1
         svctr(k,35) = svctr(k,35) + wtv_res(k)
         wtv_res(k) = wtv_res(k) * th00/g
      END DO

   END SUBROUTINE get_buoyancy
   !
   ! ---------------------------------------------------------------------
   ! get_shear:  estimates shear production term in tke budget
   !
   SUBROUTINE get_shear(n1,n2,n3,u,v,w,dzm)

      INTEGER, INTENT(in) :: n3,n2,n1
      REAL, INTENT(in)    :: w(n1,n2,n3),dzm(n1),u(n1,n2,n3),v(n1,n2,n3)

      REAL    :: ub(n1), vb(n1)
      INTEGER :: i,j,k
      REAL    :: fact, uw_shear, vw_shear

      fact = 0.25/float((n2-4)*(n3-4))

      CALL get_avg3(n1,n2,n3,u,ub)
      CALL get_avg3(n1,n2,n3,v,vb)

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 2, n1-1
               uw_shear = -(u(k,i,j)-ub(k))*fact*(                              &
                           (w(k,i,j)  +w(k,i+1,j)  )*(ub(k+1)-ub(k)  )*dzm(k) +        &
                           (w(k-1,i,j)+w(k-1,i+1,j))*(ub(k)  -ub(k-1))*dzm(k-1))
               IF (j > 1) vw_shear = -(v(k,i,j)-vb(k))*fact*(                   &
                                      (w(k,i,j)  +w(k,i,j+1)  )*(vb(k+1)-vb(k)  )*dzm(k) +  &
                                      (w(k-1,i,j)+w(k-1,i,j+1))*(vb(k)  -vb(k-1))*dzm(k-1))

               svctr(k,48) = svctr(k,48)+uw_shear
               svctr(k,36) = svctr(k,36)+uw_shear+vw_shear
            END DO
         END DO
      END DO

   END SUBROUTINE get_shear
   !
   ! ----------------------------------------------------------------------
   ! Subroutine write_ts: writes the statistics file
   !
   SUBROUTINE write_ts

      USE netcdf
      USE grid, ONLY : level

      INTEGER :: iret, n, VarID

      DO n = 1, nv1
         iret = nf90_inq_varid(ncid1, s1(n), VarID)
         iret = nf90_put_var(ncid1, VarID, ssclr(n), start=(/nrec1/))
         ssclr(n) = 0.
      END DO

      IF (level >= 4) THEN
         DO n = 1,nv1_lvl4
            iret = nf90_inq_varid(ncid1, s1_lvl4(n), VarID)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid1, VarID, ssclr_b(n), start=(/nrec1/))
            ssclr_b(n) = 0.
         END DO
      END IF

      IF (level >= 5) THEN
         DO n = 1,nv1_lvl5
            iret = nf90_inq_varid(ncid1, s1_lvl5(n), VarID)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid1, VarID, ssclr_lvl5(n), start=(/nrec1/))
            ssclr_lvl5(n) = 0.
         END DO
      END IF

      iret = nf90_sync(ncid1)
      nrec1 = nrec1 + 1

   END SUBROUTINE write_ts
   !
   ! ----------------------------------------------------------------------
   ! Subroutine write_ps: writes the time averaged elements of the
   ! statistics file
   !
   SUBROUTINE  write_ps(n1,dn0,u0,v0,zm,zt,time)

      USE netcdf
      USE defs, ONLY : alvl, cp
      USE mo_submctl, ONLY : in1a,in2b,fn2a,fn2b,fca,ica,fcb,icb,fra,ira, &
                             iia, fia, iib, fib, isa, fsa, &
                             aerobins,cloudbins,precpbins,icebins,snowbins

      INTEGER, INTENT (in) :: n1
      REAL, INTENT (in)    :: time
      REAL, INTENT (in)    :: dn0(n1), u0(n1), v0(n1), zm(n1), zt(n1)

      INTEGER :: iret, VarID, k, n, kp1

      lsttm = time
      DO k = 1, n1
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
         IF (lsttm > fsttm) THEN
            svctr(k,49) = (tke_res(k) - tke0(k))/(lsttm-fsttm)
         ELSE
            svctr(k,49) = 0.
         END IF

         svctr(k,10:nv2) = svctr(k,10:nv2)/nsmp
         IF (level >= 4 ) THEN
            svctr_b(k,:)    = svctr_b(k,:)/nsmp
            svctr_aa(k,:,:) = svctr_aa(k,:,:)/nsmp
            svctr_ab(k,:,:) = svctr_ab(k,:,:)/nsmp
            svctr_ca(k,:,:) = svctr_ca(k,:,:)/nsmp
            svctr_cb(k,:,:) = svctr_cb(k,:,:)/nsmp
            svctr_p(k,:,:)  = svctr_p(k,:,:)/nsmp

            ! Replace level 3 CDNC=CCN with that from SALSA (a + b bins)
            svctr(k,84) = svctr_b(k,8)+svctr_b(k,9)
         END IF
         IF (level >= 5) THEN
            svctr_lvl5(k,:) = svctr_lvl5(k,:)/nsmp
            svctr_ia(k,:,:) = svctr_ia(k,:,:)/nsmp
            svctr_ib(k,:,:) = svctr_ib(k,:,:)/nsmp
            svctr_s(k,:,:) = svctr_s(k,:,:)/nsmp
         ENDIF

      END DO

      iret = nf90_inq_VarID(ncid2, s2(1), VarID)
      iret = nf90_put_var(ncid2, VarID, time, start=(/nrec2/))
      IF (nrec2 == 1) THEN
         iret = nf90_inq_varid(ncid2, s2(2), VarID)
         iret = nf90_put_var(ncid2, VarID, zt, start = (/nrec2/))
         iret = nf90_inq_varid(ncid2, s2(3), VarID)
         iret = nf90_put_var(ncid2, VarID, zm, start = (/nrec2/))
         iret = nf90_inq_varid(ncid2, s2(4), VarID)
         iret = nf90_put_var(ncid2, VarID, dn0, start = (/nrec2/))
         ! Juha: For SALSA
         IF (level >= 4) THEN
            iret = nf90_inq_varid(ncid2,s2_lvl4(1),VarID)
            iret = nf90_put_var(ncid2,VarID,aerobins(in1a:fn2a),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,s2_lvl4(2),VarID)
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid2,VarID,aerobins(in2b:fn2b),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,s2_lvl4(3),VarID)
            iret = nf90_put_var(ncid2,VarID,cloudbins(ica%cur:fca%cur),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,s2_lvl4(4),VarID)
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid2,VarID,cloudbins(icb%cur:fcb%cur),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,s2_lvl4(5),VarID)
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid2,VarID,precpbins(ira:fra),start=(/nrec2/))
         END IF
         IF (level >= 5) THEN
            iret = nf90_inq_varid(ncid2,s2_lvl5(1),VarID)
            iret = nf90_put_var(ncid2,VarID,cloudbins(iia%cur:fia%cur),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,s2_lvl5(2),VarID)
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid2,VarID,icebins(iib%cur:fib%cur),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,s2_lvl5(3),VarID)
            iret = nf90_put_var(ncid2,VarID,snowbins(isa:fsa),start=(/nrec2/))
         END IF
         ! \\ SALSA
         iret = nf90_inq_varid(ncid2, s2(5), VarID)
         iret = nf90_put_var(ncid2, VarID, u0, start = (/nrec2/))
         iret = nf90_inq_varid(ncid2, s2(6), VarID)
         iret = nf90_put_var(ncid2, VarID, v0, start = (/nrec2/))
      END IF

      iret = nf90_inq_VarID(ncid2, s2(7), VarID)
      iret = nf90_put_var(ncid2, VarID, fsttm, start=(/nrec2/))
      iret = nf90_inq_VarID(ncid2, s2(8), VarID)
      iret = nf90_put_var(ncid2, VarID, lsttm, start=(/nrec2/))
      iret = nf90_inq_VarID(ncid2, s2(9), VarID)
      iret = nf90_put_var(ncid2, VarID, nsmp,  start=(/nrec2/))

      DO n = 10, nvar2
         iret = nf90_inq_varid(ncid2, s2(n), VarID)
         iret = nf90_put_var(ncid2,VarID,svctr(:,n), start=(/1,nrec2/),    &
                             count=(/n1,1/))
      END DO

      IF (level >= 4) THEN
         ! Bulk SALSA
         DO n = 6, nv2_lvl4
            iret = nf90_inq_varid(ncid2,s2_lvl4(n),VarID)
            IF (iret /= NF90_NOERR) CYCLE ! This is intended to mainly keep track of which variables are used
            iret = nf90_put_var(ncid2,VarID,svctr_b(:,n), start=(/1,nrec2/),  &
                                count=(/n1,1/))
         END DO

         ! Binned aerosols, regime a
         DO n = 1, nv2saa
            iret = nf90_inq_varid(ncid2,s2Aeroa(n),VarID)
            IF (iret /= NF90_NOERR) CYCLE ! keep track of used vars
            iret = nf90_put_var(ncid2,VarID,svctr_aa(:,:,n), start=(/1,1,nrec2/),  &
                                count=(/n1,fn2a,1/))
         END DO

         ! Aerosols, regime b
         DO n = 1, nv2sab
            iret = nf90_inq_varid(ncid2,s2Aerob(n),VarID)
            IF (iret /= NF90_NOERR) CYCLE ! Keep track of used variables
            iret = nf90_put_var(ncid2,VarID,svctr_ab(:,:,n), start=(/1,1,nrec2/),  &
                                count=(/n1,fn2b-fn2a,1/))
         END DO

         ! Cloud droplets, regime a
         DO n = 1, nv2sca
            iret = nf90_inq_varid(ncid2,s2Clouda(n),VarID)
            IF (iret /= NF90_NOERR) CYCLE ! Keep track of used variables
            iret = nf90_put_var(ncid2,VarID,svctr_ca(:,:,n), start=(/1,1,nrec2/),  &
                                count=(/n1,fca%cur,1/))
         END DO

         ! Cloud droplets, regime b
         DO n = 1, nv2scb
            iret = nf90_inq_varid(ncid2,s2Cloudb(n),VarID)
            IF (iret /= NF90_NOERR) CYCLE ! Keep track of used variables
            iret = nf90_put_var(ncid2,VarID,svctr_cb(:,:,n), start=(/1,1,nrec2/),  &
                                count=(/n1,fcb%cur-fca%cur,1/))
         END DO

         ! Precipitation
         DO n = 1, nv2sp
            iret = nf90_inq_varid(ncid2,s2Precp(n),VarID)
            IF (iret /= NF90_NOERR) CYCLE ! Keep track of used variables
            iret = nf90_put_var(ncid2,VarID,svctr_p(:,:,n), start=(/1,1,nrec2/),   &
                                count=(/n1,fra,1/))
         END DO

      END IF

      IF (level >= 5) THEN
         ! SALSA level 5
         DO n = 4,nv2_lvl5
            iret = nf90_inq_varid(ncid2,s2_lvl5(n),VarID)
            IF (iret /= NF90_NOERR) CYCLE ! This is intended to mainly keep track of which variables are used
            iret = nf90_put_var(ncid2,VarID,svctr_lvl5(:,n), start=(/1,nrec2/),  &
                                count=(/n1,1/))
         END DO

         ! Ice, regime a
         DO n = 1,nv2sia
            iret = nf90_inq_varid(ncid2,s2Icea(n),VarID)
            IF (iret /= NF90_NOERR) CYCLE
            iret = nf90_put_var(ncid2,VarID,svctr_ia(:,:,n), start=(/1,1,nrec2/),  &
                                count=(/n1,fia%cur,1/))
         END DO

         ! Ice, regime b
         DO n = 1,nv2sib
            iret = nf90_inq_varid(ncid2,s2Iceb(n),VarID)
            IF (iret /= NF90_NOERR) CYCLE
            iret = nf90_put_var(ncid2,VarID,svctr_ib(:,:,n), start=(/1,1,nrec2/),  &
                                count=(/n1,fib%cur-fia%cur,1/))
         END DO

         ! Snow
         DO n = 1,nv2ss
            iret = nf90_inq_varid(ncid2,s2Snow(n),VarID)
            IF (iret /= NF90_NOERR) CYCLE
            iret = nf90_put_var(ncid2,VarID,svctr_s(:,:,n), start=(/1,1,nrec2/),   &
                                count=(/n1,fsa,1/))
         END DO

      END IF

      iret  = nf90_sync(ncid2)
      nrec2 = nrec2+1
      nsmp  = 0.

      svctr(:,:) = 0.
      IF (level >= 4) THEN 
         svctr_b(:,:) = 0.
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

   END SUBROUTINE write_ps
   !
   ! ----------------------------------------------------------------------
   ! Subroutine: sfc_stat:  Updates statistical arrays with surface flux
   ! variables
   !
   SUBROUTINE sfc_stat(n2,n3,tflx,qflx,ustar,sst)

      INTEGER, INTENT(in) :: n2,n3
      REAL, INTENT(in), DIMENSION(n2,n3) :: tflx, qflx, ustar
      REAL, INTENT(in)    :: sst

      ssclr(10) = sst
      ssclr(11) = get_avg2dh(n2,n3,ustar)

      ssclr(12) = get_avg2dh(n2,n3,tflx)
      IF (level >= 1) ssclr(13) = get_avg2dh(n2,n3,qflx)

   END SUBROUTINE sfc_stat
   !
   ! ----------------------------------------------------------------------
   ! Subroutine: fills scalar array based on index
   ! 1: cfl; 2 max divergence
   !
   SUBROUTINE fill_scalar(index,xval)

      INTEGER, INTENT(in) :: index
      REAL, INTENT (in)   :: xval

      SELECT CASE(index)
         CASE(1)
            ssclr(2) = xval
         CASE(2)
            ssclr(3) = xval
      END SELECT

   END SUBROUTINE fill_scalar
   !
   ! ----------------------------------------------------------------------
   ! Subroutine: calculates the dissipation for output diagnostics, if
   ! isgstyp equals 2 then le is passed in via diss
   !
   SUBROUTINE sgs_vel(n1,n2,n3,v1,v2,v3)

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: v1(n1),v2(n1),v3(n1)

      svctr(:,23) = svctr(:,23)+v1(:)/float((n2-2)*(n3-2))
      svctr(:,25) = svctr(:,25)+v2(:)/float((n2-2)*(n3-2))
      svctr(:,27) = svctr(:,27)+v3(:)/float((n2-2)*(n3-2))

   END SUBROUTINE sgs_vel
   !
   ! --------------------------------------------------------------------------
   ! SGSFLXS: estimates the sgs rl and tv flux from the sgs theta_l and sgs r_t
   ! fluxes
   !
   SUBROUTINE sgsflxs(n1,n2,n3,level,rl,rv,th,flx,type)

      USE defs, ONLY : alvl, cp, rm, ep2

      INTEGER, INTENT(in) :: n1,n2,n3,level
      REAL, INTENT(in)    :: rl(n1,n2,n3),rv(n1,n2,n3)
      REAL, INTENT(in)    :: th(n1,n2,n3),flx(n1,n2,n3)
      CHARACTER (len=2)   :: type

      INTEGER :: k,i,j
      REAL    :: rnpts      ! reciprical of number of points and
      REAL    :: fctl, fctt ! factors for liquid (l) and tv (t) fluxes

      IF (type == 'tl') THEN
         wrl_sgs(:) = 0.
         wtv_sgs(:) = 0.
      END IF
      rnpts = 1./REAL((n2-4)*(n3-4))
      !
      ! calculate fluxes assuming the possibility of liquid water.  if liquid
      ! water does not exist sgs_rl = 0.
      !
      IF ( level >= 2 ) THEN
         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1-1
                  IF (rl(k+1,i,j) > 0.) THEN
                     fctt = rnpts*(1. + rv(k,i,j)*(1.+ep2 +ep2*rv(k,i,j)*alvl     &
                                   /(rm*th(k,i,j))))                              &
                            /(1.+(rv(k,i,j)*(alvl/th(k,i,j))**2)/(rm*cp))
                     SELECT CASE (type)
                        CASE ('tl')
                           fctl = -rnpts/(rm*th(k,i,j)**2/(rv(k,i,j)*alvl)+alvl/cp)
                        CASE ('rt')
                           fctl = rnpts/(1.+(rv(k,i,j)*alvl**2)/(cp*rm*th(k,i,j)**2))
                           fctt = (alvl*fctt/cp - th(k,i,j)*rnpts)
                     END SELECT
                     wrl_sgs(k) = wrl_sgs(k) + fctl*flx(k,i,j)
                     wtv_sgs(k) = wtv_sgs(k) + fctt*flx(k,i,j)
                  ELSE
                     SELECT CASE (type)
                        CASE ('tl')
                           fctt = rnpts*(1. + ep2*rv(k,i,j))
                        CASE ('rt')
                           fctt = rnpts*(ep2*th(k,i,j))
                     END SELECT
                     wtv_sgs(k) = wtv_sgs(k) + fctt*flx(k,i,j)
                  END IF
               END DO
            END DO
         END DO
         !
         ! calculate fluxes for dry thermodynamics, i.e., wrl_sgs is by def
         ! zero
         !
      ELSE
         DO k = 1, n1
            wrl_sgs(k) = 0.
         END DO
         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1-1
                  IF ( level >= 1) THEN
                     SELECT CASE (type)
                        CASE ('tl')
                           fctt = rnpts * (1. + ep2*rv(k,i,j))
                        CASE ('rt')
                           fctt = rnpts * ep2*th(k,i,j)
                     END SELECT
                  ELSE
                     fctt = rnpts
                  END IF
                  wtv_sgs(k) = wtv_sgs(k) + fctt*flx(k,i,j)
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE sgsflxs
   !
   ! ----------------------------------------------------------------------
   ! Subroutine fill_tend: fills arrays with current value of tendencies
   !
   SUBROUTINE acc_tend(n1,n2,n3,f1,f2,f3,t1,t2,t3,v1,v2,v3,ic,routine)

      INTEGER, INTENT(in) :: n1,n2,n3,ic
      REAL, INTENT(in)    :: f1(n1,n2,n3),f2(n1,n2,n3),f3(n1,n2,n3)
      REAL, INTENT(in)    :: t1(n1,n2,n3),t2(n1,n2,n3),t3(n1,n2,n3)
      REAL, INTENT(inout) :: v1(n1),v2(n1),v3(n1)
      CHARACTER (len=3)   :: routine

      INTEGER :: k,ii
      REAL    :: x1(n1),x2(n1),x3(n1)

      CALL get_cor3(n1,n2,n3,f1,t1,x1)
      CALL get_cor3(n1,n2,n3,f2,t2,x2)
      CALL get_cor3(n1,n2,n3,f3,t3,x3)

      SELECT CASE (routine)
         CASE ('sgs')
            ii = 39
         CASE ('adv')
            ii = 42
      END SELECT

      SELECT CASE (ic)
         CASE (1)
            DO k = 1, n1
               v1(k) = x1(k)
               v2(k) = x2(k)
               v3(k) = x3(k)
            END DO
         CASE (2)
            DO k = 1, n1
               svctr(k,ii)   = svctr(k,ii)   + (x1(k)-v1(k))
               svctr(k,ii+1) = svctr(k,ii+1) + (x2(k)-v2(k))
               svctr(k,ii+2) = svctr(k,ii+2) + (x3(k)-v3(k))
            END DO
      END SELECT

   END SUBROUTINE acc_tend
   !
   !---------------------------------------------------------------------
   ! Subroutine updtst: updates appropriate statistical arrays
   !
   SUBROUTINE updtst(n1,routine,nfld,values,ic)

      INTEGER, INTENT(in)            :: n1,nfld,ic
      REAL, INTENT (in)              :: values(n1)
      CHARACTER (len=3), INTENT (in) :: routine

      INTEGER :: nn,k

      SELECT CASE (routine)
         CASE("sgs")
            SELECT CASE (nfld)
               CASE (-6)
                  nn = 31 ! dissipation length-scale
               CASE (-5)
                  nn = 30 ! mixing length
               CASE (-4)
                  nn = 29 ! eddy diffusivity
               CASE (-3)
                  nn = 28 ! eddy viscosity
               CASE (-2)
                  nn = 38 ! dissipation
               CASE (-1)
                  nn = 32 ! estimated sgs energy
               CASE (1)
                  nn = 21 ! sgs tl flux
               CASE (2)
                  nn = 54 ! sgs rt flux
               CASE DEFAULT
                  nn = 0
            END SELECT
         CASE("adv")
            SELECT CASE (nfld)
               CASE (-3)
                  nn = 26 ! adv w flux
               CASE (-2)
                  nn = 24 ! adv v flux
               CASE (-1)
                  nn = 22 ! adv u flux
               CASE (0)
                  nn = 62 ! adv rl flux
               CASE (1)
                  nn = 20 ! adv tl flux
               CASE (2)
                  nn = 53 ! adv rt flux
               CASE DEFAULT
                  nn = 0
            END SELECT
         CASE("prs")
            SELECT CASE (nfld)
               CASE (1)
                  nn = 45 ! dpdx u corr
               CASE (2)
                  nn = 46 ! dpdy v corr
               CASE (3)
                  nn = 47 ! dpdz w corr
               CASE DEFAULT
                  nn = 0
            END SELECT
         CASE("prc")
            SELECT CASE (nfld)
               CASE (2)
                  nn = 88
               CASE (3)
                  nn = 63
               CASE DEFAULT
                  nn = 0
            END SELECT
         CASE DEFAULT
            nn = 0
      END SELECT

      IF (nn > 0) THEN
         IF (ic == 0) svctr(:,nn) = 0.
         DO k = 1, n1
            svctr(k,nn) = svctr(k,nn)+values(k)
         END DO
      END IF

   END SUBROUTINE updtst

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
      USE mo_submctl, ONLY : nbins, ncld, nprc, nice,  nsnw, spec
      IMPLICIT NONE

      INTEGER, INTENT(in)        :: n2,n3,n4                     ! Grid dimensions
      REAL, INTENT(in)           :: raer(n2,n3,n4*nbins)        ! Array containing the binned 2d-field
      REAL,           INTENT(in) :: rcld(n2,n3,n4*ncld), &
                                    rprc(n2,n3,n4*nprc), &    
                                    rice(n2,n3,n4*nice), &
                                    rsnw(n2,n3,n4*nsnw)

      INTEGER :: ss, si
      INTEGER :: tt
      INTEGER :: end,str
      REAL :: zavg


      DO ss = 1,8
         IF ( .NOT. spec%isUsed(zspec(ss)) ) CYCLE

         si = spec%getIndex(zspec(ss))

         ! Index to ssclr_b and s1SalsaBulk
         tt = 29 +(ss-1)*3

         ! Removal by sedimentation of aerosol
         str = getMassIndex(nbins,1,si)
         end = getMassIndex(nbins,nbins,si)
         zavg = get_avg2dh( n2,n3,SUM(raer(:,:,str:end),DIM=3) )
         ssclr_b(tt) = ssclr_b(tt) + zavg
         tt=tt+1

         ! Removal by sedimentation of cloud droplets
         str = getMassIndex(ncld,1,si)
         end = getMassIndex(ncld,ncld,si)
         zavg = get_avg2dh( n2,n3,SUM(rcld(:,:,str:end),DIM=3) )
         ssclr_b(tt) = ssclr_b(tt) + zavg
         tt=tt+1

         ! Removal by precipitation
         str = getMassIndex(nprc,1,si)
         end = getMassIndex(nprc,nprc,si)
         zavg = get_avg2dh( n2,n3,SUM(rprc(:,:,str:end),DIM=3) )
         ssclr_b(tt) = ssclr_b(tt) + zavg
         tt=tt+1

         IF (level<5) CYCLE

         ! Index to ssclr_lvl5 and s1_lvl5
         tt = 34 +(ss-1)*2

         ! Removal by sedimentation of ice particles
         str = getMassIndex(nice,1,si)
         end = getMassIndex(nice,nice,si)
         zavg = get_avg2dh( n2,n3,SUM(rice(:,:,str:end),DIM=3) )
         ssclr_lvl5(tt) = ssclr_lvl5(tt) + zavg
         tt=tt+1

         ! Removal by snow
         str = getMassIndex(nsnw,1,si)
         end = getMassIndex(nsnw,nsnw,si)
         zavg = get_avg2dh( n2,n3,SUM(rsnw(:,:,str:end),DIM=3) )
         ssclr_lvl5(tt) = ssclr_lvl5(tt) + zavg

      END DO

   END SUBROUTINE acc_removal
   !
   !--------------------------------------------------------------------------
   !
   ! Accumulate mass budget terms. Must be done for every timestep. Is called from step,srfc,mcrp
   ! Juha Tonttila, FMI, 2016
   !
   SUBROUTINE acc_massbudged(n1,n2,n3,type,tstep,dz,dn,    &
                             rv,rc,prc,revap,rdep,ApVdom   )
      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1,n2,n3  !
      INTEGER, INTENT(in) :: type      ! 1: Atmospheric water, 2: Water evaporated to the atm,
                                       ! 3: Water removed by precip/deposition
      REAL, INTENT(in) :: tstep        ! Current timestep length

      REAL, INTENT(in) :: dz(n1)
      REAL, INTENT(in) :: dn(n1,n2,n3) ! Air density

      REAL, INTENT(in), OPTIONAL :: rv(n1,n2,n3)      ! Water vapor mixing ratio
      REAL, INTENT(in), OPTIONAL :: rc(n1,n2,n3), &   ! Cloud water mixing ratio
                                    prc(n1,n2,n3)     ! Precipitation mixing ratio

      REAL, INTENT(in), OPTIONAL :: revap(n2,n3), &   ! Gain of water through evaporation at the surface
                                    rdep(n2,n3)       ! Loss of water through deposition

      REAL, INTENT(in), OPTIONAL :: ApVdom            ! Domain surface area / Domain atm volume

      REAL :: a1,a2,a3

      SELECT CASE(type)
         CASE(0)
            IF ( .NOT. present(rv) .OR. .NOT. present(rc) .OR. .NOT. present(prc) ) &
               STOP 'acc_massbudget (stat): ERROR - for atm water q,rc and prc must be present'

            a1 = 0.; a2 = 0.; a3 = 0.
            a1 = get_avg_ts(n1,n2,n3,rv*dn,dz)
            a2 = get_avg_ts(n1,n2,n3,rc*dn,dz)
            a3 = get_avg_ts(n1,n2,n3,prc*dn,dz)
            massbdg(1) = (a1 + a2 + a3)

         CASE(1)
            IF ( .NOT. present(rv) .OR. .NOT. present(rc) .OR. .NOT. present(prc) ) &
               STOP 'acc_massbudget (stat): ERROR - for atm water q,rc and prc must be present'

            a1 = 0.; a2 = 0.; a3 = 0.
            a1 = get_avg_ts(n1,n2,n3,rv*dn,dz)
            a2 = get_avg_ts(n1,n2,n3,rc*dn,dz)
            a3 = get_avg_ts(n1,n2,n3,prc*dn,dz)
            massbdg(2) = (a1 + a2 + a3) ! Not accumulated

         CASE(2)
            IF ( .NOT. present(revap) .OR. .NOT. present(ApVdom) ) &
               STOP 'acc_massbudget (stat): ERROR - for evaporation stats revap must be present'

            a1 = 0.
            a1 = get_avg2dh(n2,n3,revap)
            massbdg(3) = massbdg(3) + a1*tstep*ApVdom

         CASE(3)
            IF ( .NOT. present(rdep) .OR. .NOT. present(ApVdom) ) &
               STOP 'acc_massbudget (stat): ERROR - for deposition stats rdep must be present'

            a1 = 0.
            a1 = get_avg2dh(n2,n3,rdep)
            massbdg(4) = massbdg(4) + a1*tstep*ApVdom

      END SELECT

      avgtime = avgtime + tstep ! Accumulate time for normalization

   END SUBROUTINE acc_massbudged
   !
   ! -------------------------------------------------------------------------
   !
   SUBROUTINE write_massbudged()
      IMPLICIT NONE

      OPEN(88,FILE='MASSBUDGED.TXT')
      WRITE(88,*) 'Initial mass (atm), Final mass (atm), Final-Initial, Total evpaoration, Total deposition, Evap-Dep  '
      WRITE(88,*) massbdg(1), massbdg(2), massbdg(2)-massbdg(1),  &
                  massbdg(3), massbdg(4), massbdg(3)-massbdg(4)
      CLOSE(88)

   END SUBROUTINE write_massbudged
   !
   ! -------------------------------------------------------------------------
   !
   INTEGER FUNCTION close_stat()

      USE netcdf

      close_stat = nf90_close(ncid1) + nf90_close(ncid2)
      IF (csflg) close_stat = close_stat + nf90_close(ncid3)

   END FUNCTION close_stat
   !
   ! -------------------------------------------------------------------------
   !
   REAL FUNCTION get_zi (n1, n2, n3, itype, sx, xx, z, threshold)

      INTEGER, INTENT (in) :: n1, n2, n3, itype
      REAL, INTENT (in)    :: xx(n1), z(n1), sx(n1,n2,n3), threshold

      INTEGER :: i, j, k, kk
      REAL    :: zibar, sval, dmy, scr(n2,n3)

      get_zi = -999.
      SELECT CASE(itype)
         CASE (1)
            !
            ! find level at which sx=threshold (xx is one over grid spacing)
            !
            zibar = 0.
            DO j = 3, n3-2
               DO i = 3, n2-2
                  k = 2
                  DO WHILE (k < n1-2 .AND. sx(k,i,j) > threshold)
                     k = k+1
                  END DO
                  IF (k == n1-2) zibar = -999.
                  IF (zibar /= -999.) zibar = zibar + z(k-1) +  &
                                             (threshold - sx(k-1,i,j))/xx(k-1)     /  &
                                             (sx(k,i,j) - sx(k-1,i,j) + epsilon(1.))
               END DO
            END DO
            IF (zibar /= -999.) get_zi = zibar/REAL((n3-4)*(n2-4))

         CASE(2)
            !
            ! find level of maximum gradient (xx is one over grid spacing)
            !
            scr = 0.
            DO j = 3, n3-2
               DO i = 3, n2-2
                  sval = 0.
                  DO k = 2, n1-5
                     dmy = (sx(k+1,i,j)-sx(k,i,j))*xx(k)
                     IF (dmy > sval) THEN
                        sval = dmy
                        scr(i,j) = z(k)
                     END IF
                  END DO
               END DO
            END DO
            get_zi = get_avg2dh(n2,n3,scr)

         CASE(3)
            !
            ! find level where xx is a maximum
            !
            sval = -huge(1.)
            kk = 1
            DO k = 2, n1
               IF (xx(k) > sval) THEN
                  kk = k
                  sval = xx(k)
               END IF
            END DO
            get_zi = z(kk)

         CASE(4)
            !
            ! find level where xx is a minimum
            !
            sval = huge(1.)
            kk = 1
            DO k = 2, n1-2
               IF (xx(k) < sval) THEN
                  kk = k
                  sval = xx(k)
               END IF
            END DO
            get_zi = z(kk)
      END SELECT

   END FUNCTION get_zi

END MODULE stat

