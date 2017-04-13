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

  use ncio, only : open_nc, define_nc, define_nc_cs
  use grid, only : level
  use util, only : get_avg3, get_cor3, get_var3, get_avg_ts, &
                   get_avg2dh, get_3rd3

  implicit none
  private

  integer, parameter :: nvar1 = 28,               &
                        nv1sbulk = 62,            &
                        nv1MB = 4,                &
                        nvar2 = 92,               &
                        nv2sbulk = 40,            &
                        nv2saa = 8, nv2sab = 8,   &
                        nv2sca = 8, nv2scb = 8,   &
                        nv2sp = 8

  integer, save      :: nvar_spec3=0

  integer, save      :: nrec1, nrec2, nrec3, ncid1, ncid2, ncid3, nv1=nvar1, nv2 = nvar2
  real, save         :: fsttm, lsttm, nsmp = 0
  REAL, SAVE         :: avgtime = 0

  logical            :: sflg = .false.
  LOGICAL            :: mcflg = .FALSE.
  LOGICAL            :: csflg = .FALSE.
  real               :: ssam_intvl = 30.   ! statistical sampling interval
  real               :: savg_intvl = 1800. ! statistical averaging interval

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
       'CCN    ','nrain  ','nrcnt  ','nccnt  '/),                   & !25

       ! **** Bulk temporal statistics for SALSA ****
       s1SalsaBulk(nv1sbulk) = (/                                    &
       'Nc_ic  ','Na_int ','Na_oc  ',                                & !1
       'SO4_ic ','SO4_int','SO4_oc ',                                & !4
       'OC_ic  ','OC_int ','OC_oc  ',                                & !7
       'BC_ic  ','BC_int ','BC_oc  ',                                & !10
       'DU_ic  ','DU_int ','DU_oc  ',                                & !13
       'SS_ic  ','SS_int ','SS_oc  ',                                & !16
       'NH3_ic ','NH3_int','NH3_oc ',                                & !19
       'NO3_ic ','NO3_int','NO3_oc ',                                & !22
       'rmH2Oae','rmH2Ocl','rmH2Opr',                                & !25
       'rmSO4dr','rmSO4cl','rmSO4pr','rmSO4wt','rmSO4tt',            & !28
       'rmOCdr ','rmOCcl ','rmOCpr ','rmOCwt ','rmOCtt ',            & !33
       'rmBCdr ','rmBCcl ','rmBCpr ','rmBCwt ','rmBCtt ',            & !38
       'rmDUdr ','rmDUcl ','rmDUpr ','rmDUwt ','rmDUtt ',            & !43
       'rmSSdr ','rmSScl ','rmSSpr ','rmSSwt ','rmSStt ',            & !48
       'rmNH3dr','rmNH3cl','rmNH3pr','rmNH3wt','rmNH3tt',            & !53
       'rmNO3dr','rmNO3cl','rmNO3pr','rmNO3wt','rmNO3tt'             & !58, total 62
       /),                                                           &

        s2(nvar2)=(/                                                 &
        'time   ','zt     ','zm     ','dn0    ','u0     ','v0     ', & ! 1
        'fsttm  ','lsttm  ','nsmp   ','u      ','v      ','t      ', & ! 7
        'p      ','u_2    ','v_2    ','w_2    ','t_2    ','w_3    ', & ! 13
        't_3    ','tot_tw ','sfs_tw ','tot_uw ','sfs_uw ','tot_vw ', & ! 19
        'sfs_vw ','tot_ww ','sfs_ww ','km     ','kh     ','lmbd   ', & ! 25
        'lmbde  ','sfs_tke','sfs_boy','sfs_shr','boy_prd','shr_prd', & ! 31
        'trans  ','diss   ','dff_u  ','dff_v  ','dff_w  ','adv_u  ', & ! 37
        'adv_v  ','adv_w  ','prs_u  ','prs_v  ','prs_w  ','prd_uw ', & ! 43
        'storage','q      ','q_2    ','q_3    ','tot_qw ','sfs_qw ', & ! 49
        'rflx   ','rflx2  ','sflx   ','sflx2  ','l      ','l_2    ', & ! 55
        'l_3    ','tot_lw ','sed_lw ','cs1    ','cnt_cs1','w_cs1  ', & ! 61
        'tl_cs1 ','tv_cs1 ','rt_cs1 ','rl_cs1 ','wt_cs1 ','wv_cs1 ', & ! 67
        'wr_cs1 ','cs2    ','cnt_cs2','w_cs2  ','tl_cs2 ','tv_cs2 ', & ! 73
        'rt_cs2 ','rl_cs2 ','wt_cs2 ','wv_cs2 ','wr_cs2 ','Nc     ', & ! 79
        'Nr     ','rr     ','precip ','evap   ','frc_prc','prc_prc', & ! 85
        'frc_ran','hst_srf'/),                                       & ! 91, total 92

        ! **** BULK PROFILE OUTPUT FOR SALSA ****
        s2SalsaBulk(nv2sbulk) = (/                                   &
        'aea    ','aeb    ','cla    ','clb    ','prc    ',           & !1
        'P_Naa  ','P_Nab  ','P_Nca  ','P_Ncb  ','P_Np   ',           & !6
        'P_Rwaa ','P_Rwab ','P_Rwca ','P_Rwcb ','P_Rwp  ',           & !11
        'P_cSO4a','P_cSO4c','P_cSO4p',                               & !16
        'P_cOCa ','P_cOCc ','P_cOCp ',                               & !19
        'P_cBCa ','P_cBCc ','P_cBCp ',                               & !22
        'P_cDUa ','P_cDUc ','P_cDUp ',                               & !25
        'P_cSSa ','P_cSSc ','P_cSSp ',                               & !28
        'P_cNH3a','P_cNH3c','P_cNH3p',                               & !31
        'P_cNO3a','P_cNO3c','P_cNO3p',                               & !34
        'P_rl   ','P_rr   ','P_rv   ','P_RH   '/),                   & !37, total 40

        ! **** BINNED PROFILE OUTPUT FOR SALSA ****
        ! **** Aerosols
        s2Aeroa(nv2saa) = (/                                         &
        'P_Naba ','P_SO4aa','P_OCaa ','P_BCaa ',                     &
        'P_DUaa ','P_SSaa ','P_NH3aa','P_NO3aa'/),                   &
        s2Aerob(nv2sab) = (/                                         &
        'P_Nabb ','P_SO4ab','P_OCab ','P_BCab ',                     &
        'P_DUab ','P_SSab ','P_NH3ab','P_NO3ab'/),                   &

        ! **** Clouds
        s2Clouda(nv2sca) = (/                                        &
        'P_Ncba ','P_SO4ca','P_OCca ','P_BCca ',                     &
        'P_DUca ','P_SSca ','P_NH3ca','P_NO3ca'/),                   &
        s2Cloudb(nv2scb) = (/                                        &
        'P_Ncbb ','P_SO4cb','P_OCcb ','P_BCcb ',                     &
        'P_DUcb ','P_SScb ','P_NH3cb','P_NO3cb'/),                   &

        ! **** Precip
        s2Precp(nv2sp) = (/                                         &
        'P_Npb  ','P_SO4pb','P_OCpb ','P_BCpb ',                     &
        'P_DUpb ','P_SSpb ','P_NH3pb','P_NO3pb'/),                   &

        s1Total(nvar1+nv1sbulk),                                     &
        s2Total(nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+nv2scb+nv2sp)

  character (len=7), save :: spec3(7)


  LOGICAL, save :: s2bool(nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+nv2scb+nv2sp)
  LOGICAL, save :: s1bool(nvar1+nv1sbulk)

  real, save, allocatable   :: tke_sgs(:), tke_res(:), tke0(:), wtv_sgs(:),  &
       wtv_res(:), wrl_sgs(:), thvar(:), svctr(:,:), ssclr(:),               &
       ! Additional ssclr and svctr for BULK SALSA output
       svctr_b(:,:), ssclr_b(:),                                             &
       ! Additional ssclr and svctr for BINNED SALSA output.
       svctr_aa(:,:,:), svctr_ca(:,:,:), svctr_p(:,:,:),                     &
       svctr_ab(:,:,:), svctr_cb(:,:,:),                                     &
       ! Mass budget arrays
       massbdg(:), scs_rm(:,:,:)

  public :: sflg, ssam_intvl, savg_intvl, statistics, init_stat, write_ps,   &
       acc_tend, updtst, sfc_stat, close_stat, fill_scalar, tke_sgs, sgsflxs,&
       sgs_vel, comp_tke, get_zi, acc_removal, cs_rem_set, acc_massbudged, write_massbudged, mcflg, csflg

contains
  !
  ! ---------------------------------------------------------------------
  ! INIT_STAT:  This routine initializes the statistical arrays which
  ! are user/problem defined.  Note that svctr is given 100 elements, and
  ! elements 90 and above are used for computing the TKE budget. Hence
  ! if (nvar2 >= 90 the program stops
  !
  subroutine init_stat(time, filprf, expnme, nzp)

    use grid, only : nxp, nyp, iradtyp, prtcl
    use mpi_interface, only : myid, ver, author, info
    use mo_submctl, only : nprc, fn2a,fn2b,fca,fcb,fra
    USE class_ComponentIndex, ONLY : IsUsed

    character (len=80), intent (in) :: filprf, expnme
    integer, intent (in)            :: nzp
    real, intent (in)               :: time

    INTEGER :: i,e
    character (len=80) :: fname

    allocate (wtv_sgs(nzp),wtv_res(nzp),wrl_sgs(nzp))
    allocate (tke_res(nzp),tke_sgs(nzp),tke0(nzp),thvar(nzp))

    ! Combine name arrays
    i = 1; e = nvar1
    s1Total(i:e) = s1
    i = e + 1; e = e + nv1sbulk
    s1Total(i:e) = s1SalsaBulk
    ! ---
    i = 1; e = nvar2
    s2Total(i:e) = s2
    i = e + 1; e = e + nv2sbulk
    s2Total(i:e) = s2SalsaBulk
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
    case (4)
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
       ALLOCATE ( ssclr_b(nv1sbulk), svctr_b(nzp,nv2sbulk))
       ALLOCATE ( svctr_aa(nzp,fn2a,nv2saa), svctr_ab(nzp,fn2b-fn2a,nv2sab),          &
                  svctr_ca(nzp,fca%cur,nv2sca), svctr_cb(nzp,fcb%cur-fca%cur,nv2scb), &
                  svctr_p(nzp,nprc,nv2sp)   )
       IF (mcflg) &
            ALLOCATE ( massbdg(nv1MB) ) ! Mass budged array
       svctr(:,:) = 0.
       ssclr(:) = 0.
       svctr_b(:,:) = 0.
       ssclr_b(:) = 0.
       svctr_aa(:,:,:) = 0.; svctr_ab(:,:,:) = 0.
       svctr_ca(:,:,:) = 0.; svctr_cb(:,:,:) = 0.
       svctr_p(:,:,:) = 0.
       IF (mcflg) massbdg(:) = 0.

       ! Create a boolean array for items that are actually used
       s2bool(:) = .FALSE.
       s1bool(:) = .FALSE.

       s1bool(1:nvar1) = .TRUE.
       s2bool(1:nvar2) = .TRUE.     ! Original LES vars (assume always used...)

       s1bool(nvar1+1:nvar1+3) = .TRUE.  ! Number concentrations
       s2bool(nvar2+1:nvar2+15) = .TRUE. ! Bin dimensions, number concentrations and radius

       ! Bin number concentrations
       i = nvar2+nv2sbulk+1  ! binned
       s2bool(i) = .TRUE.
       i = nvar2+nv2sbulk+nv2saa+1
       s2bool(i) = .TRUE.
       i = nvar2+nv2sbulk+nv2saa+nv2sab+1
       s2bool(i) = .TRUE.
       i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+1
       s2bool(i) = .TRUE.
       i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+nv2scb+1
       s2bool(i) = .TRUE.

       ! Water removal temporal statistics
       i = nvar1+25; e = nvar1+27
       s1bool(i:e) = .TRUE.

       nvar_spec3 = 0

       IF (IsUsed(prtcl,'SO4')) THEN

          i = nvar1+4; e = nvar1+6
          s1bool(i:e) = .TRUE.
          i = nvar1+28; e = nvar1+32
          s1bool(i:e) = .TRUE.

          i = nvar2+16; e = nvar2+18
          s2bool(i:e) = .TRUE.  ! Bulk sulphate
          i = nvar2+nv2sbulk+2  ! binned
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+2
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+2
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+2
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+nv2scb+2
          s2bool(i) = .TRUE.

         nvar_spec3 = nvar_spec3 +1
         spec3(nvar_spec3)='SO4'
       END IF
       IF (IsUsed(prtcl,'OC')) THEN

          i = nvar1+7; e = nvar1+9
          s1bool(i:e) = .TRUE.
          i = nvar1+33; e = nvar1+37
          s1bool(i:e) = .TRUE.

          i = nvar2+19; e = nvar2+21
          s2bool(i:e) = .TRUE.  ! Bulk OC
          i = nvar2+nv2sbulk+3  ! binned
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+3
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+3
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+3
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+nv2scb+3
          s2bool(i) = .TRUE.

         nvar_spec3 = nvar_spec3 +1
         spec3(nvar_spec3)='OC'
       END IF
       IF (IsUsed(prtcl,'BC')) THEN

          i = nvar1+10; e = nvar1+12
          s1bool(i:e) = .TRUE.
          i = nvar1+38; e = nvar1+42
          s1bool(i:e) = .TRUE.

          i = nvar2+22; e = nvar2+24
          s2bool(i:e) = .TRUE.  ! Bulk BC
          i = nvar2+nv2sbulk+4  ! binned
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+4
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+4
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+4
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+nv2scb+4
          s2bool(i) = .TRUE.

         nvar_spec3 = nvar_spec3 +1
         spec3(nvar_spec3)='BC'
       END IF
       IF (IsUsed(prtcl,'DU')) THEN

          i = nvar1+13; e = nvar1 + 15
          s1bool(i:e) = .TRUE.
          i = nvar1+43; e = nvar1+47
          s1bool(i:e) = .TRUE.

          i = nvar2+25; e = nvar2+27
          s2bool(i:e) = .TRUE.  ! Bulk DU
          i = nvar2+nv2sbulk+5  ! binned
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+5
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+5
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+5
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+nv2scb+5
          s2bool(i) = .TRUE.

         nvar_spec3 = nvar_spec3 +1
         spec3(nvar_spec3)='DU'
       END IF
       IF (IsUsed(prtcl,'SS')) THEN

          i = nvar1+16; e = nvar1+18
          s1bool(i:e) = .TRUE.
          i = nvar1+48; e = nvar1+52
          s1bool(i:e) = .TRUE.

          i = nvar2+28; e = nvar2+30
          s2bool(i:e) = .TRUE.  ! Bulk SS
          i = nvar2+nv2sbulk+6  ! binned
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+6
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+6
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+6
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+nv2scb+6
          s2bool(i) = .TRUE.

         nvar_spec3 = nvar_spec3 +1
         spec3(nvar_spec3)='SS'
       END IF
       IF (IsUsed(prtcl,'NH')) THEN

          i = nvar1+19; e = nvar1+21
          s1bool(i:e) = .TRUE.
          i = nvar1+53; e = nvar1+57
          s1bool(i:e) = .TRUE.

          i = nvar2+31; e = nvar2+33
          s2bool(i:e) = .TRUE.  ! Bulk NH
          i = nvar2+nv2sbulk+7  ! binned
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+7
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+7
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+7
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+nv2scb+7
          s2bool(i) = .TRUE.

         nvar_spec3 = nvar_spec3 +1
         spec3(nvar_spec3)='NH'
       END IF
       IF (IsUsed(prtcl,'NO')) THEN

          i = nvar1+22; e = nvar1+24
          s1bool(i:e) = .TRUE.
          i = nvar1+58; e = nvar1+62
          s1bool(i:e) = .TRUE.

          i = nvar2+34; e = nvar2+36
          s2bool(i:e) = .TRUE.  ! Bulk NO
          i = nvar2+nv2sbulk+8  ! binned
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+8
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+8
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+8
          s2bool(i) = .TRUE.
          i = nvar2+nv2sbulk+nv2saa+nv2sab+nv2sca+nv2scb+8
          s2bool(i) = .TRUE.

         nvar_spec3 = nvar_spec3 +1
         spec3(nvar_spec3)='NO'
       END IF

       s2bool(nvar2+37:nvar2+40) = .TRUE.     ! Water mixing ratios

       IF (csflg .and. level>3) THEN
           ! Allocate array for level 4 removal rate column statistics
           ! Total number of ouputs is 3 for warm (aerosol, cloud and precipitation)
           ! and 5 (add ice and snow) for each species including water
           IF (level==4) THEN
              ALLOCATE( scs_rm(3*(nvar_spec3+1),nxp,nyp) )
           ELSE
              ALLOCATE( scs_rm(5*(nvar_spec3+1),nxp,nyp) )
          ENDIF
          scs_rm=0. ! Set to zero (during spinup)
       END IF
    END IF ! If level >=4

    fname =  trim(filprf)//'.ts'
    if(myid == 0) print                                                  &
         "(//' ',49('-')/,' ',/,'  Initializing: ',A20)",trim(fname)
    call open_nc( fname, expnme, time, (nxp-4)*(nyp-4), ncid1, nrec1, ver, author, info)
    ! Juha: Modified for SALSA output
    call define_nc( ncid1, nrec1, COUNT(s1bool), PACK(s1Total,s1bool))
    if (myid == 0) print *, '   ...starting record: ', nrec1



    fname =  trim(filprf)//'.ps'
    if(myid == 0) print                                                  &
         "(//' ',49('-')/,' ',/,'  Initializing: ',A20)",trim(fname)
    call open_nc( fname, expnme, time,(nxp-4)*(nyp-4), ncid2, nrec2, ver, author, info)
    ! Juha: Modified due to SALSA output
    call define_nc( ncid2, nrec2, COUNT(s2bool), PACK(s2Total,s2bool), n1=nzp, inae_a=fn2a, inae_b=fn2b-fn2a, &
                    incld_a=fca%cur, incld_b=fcb%cur-fca%cur, inprc=fra)
    if (myid == 0) print *, '   ...starting record: ', nrec2



    ! Optional column statistics
    IF (csflg) THEN
        fname =  trim(filprf)//'.cs'
        if(myid == 0) print "(//' ',49('-')/,' ',/,'  Initializing: ',A20)",trim(fname)
        call open_nc( fname, expnme, time,(nxp-4)*(nyp-4), ncid3, nrec3, ver, author, info)
        IF (ncid3>=0) CALL define_nc_cs(ncid3, nrec3, nxp-4, nyp-4, level, iradtyp, spec3, nvar_spec3)
        if (myid == 0) print *, '   ...starting record: ', nrec3
    ENDIF

  end subroutine init_stat
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

    use grid, only : a_up, a_vp, a_wp, a_rc, a_theta, a_rsl           &
         , a_rp, a_tp, a_press, nxp, nyp, nzp, dzm, dzt, zm, zt, th00, umean            &
         , vmean, dn0, precip, a_rpp, a_npp, CCN, iradtyp, a_rflx               &
         , a_sflx, albedo, a_srp, a_snrp, a_ncloudp, a_nprecpp, xt, yt

    real, intent (in) :: time

    real :: rxt(nzp,nxp,nyp), rnt(nzp,nxp,nyp)
    REAL :: xrpp(nzp,nxp,nyp), xnpp(nzp,nxp,nyp)

    SELECT CASE(level)
       CASE(1,2,3)
          rxt = a_rp
          xrpp = a_rpp
          xnpp = a_npp
       CASE(4,5)
          rxt = a_rp + a_rc
          xrpp = a_srp
          xnpp = a_snrp
    END SELECT

    if (nsmp == 0.) fsttm = time
    nsmp=nsmp+1.
    ssclr(14:nvar1) = -999.
    !
    ! profile statistics
    !
    call accum_stat(nzp, nxp, nyp, a_up, a_vp, a_wp, a_tp, a_press, umean &
         ,vmean,th00)
    if (iradtyp == 3) then
       call accum_rad(nzp, nxp, nyp, a_rflx, sflx=a_sflx, alb=albedo)
    elseif (iradtyp > 0) then
       call accum_rad(nzp, nxp, nyp, a_rflx)
    end if
    if (level >=1) call accum_lvl1(nzp, nxp, nyp, rxt)
    if (level >=2) call accum_lvl2(nzp, nxp, nyp, th00, dn0, zm, a_wp,        &
                                   a_theta, a_tp, a_rc, a_rsl, rxt   )
    if (level >=3) call accum_lvl3(nzp, nxp, nyp, dn0, zm, a_rc, xrpp,  &
                                   xnpp, precip, CCN                    )
    if (level >=4)  call accum_lvl4(nzp, nxp, nyp)
     !for Salsa output in ps files .. by Zubair Maalick

    !
    ! scalar statistics
    !
    call set_ts(nzp, nxp, nyp, a_wp, a_theta, dn0, zt,zm,dzt,dzm,th00,time)
    IF ( level >=1 ) CALL ts_lvl1(nzp, nxp, nyp, dn0, zt, dzm, rxt)
    IF ( level >=2 ) CALL ts_lvl2(nzp, nxp, nyp, rxt, a_rsl, zt)
    IF ( level >=4 ) CALL ts_lvl4(nzp, nxp, nyp, a_rc)

    call write_ts

    !
    ! Column statistics
    !
    IF (csflg) THEN
        ! Radiation
        if (iradtyp == 3) CALL set_cs_any(nxp,nyp,albedo,'albedo')

        ! Deposition statistics
        IF (level == 3) THEN
            CALL set_cs_any(nxp,nyp,precip(2,:,:),'prcp')
        ELSEIF (level>3) THEN
            ! Surface deposition fluxes
            CALL cs_rem_save(nxp,nyp)
        ENDIF

        ! Warm cloud statistics
        rxt = 0.    ! Condensate
        rnt = 0.
        xrpp = 0.   ! Precipitate
        xnpp = 0.
        IF (level==1) THEN
            ! No clouds or precipitation
        ELSEIF (level==2) THEN
            ! Clouds available
            rxt = a_rc
            rnt = CCN
        ELSEIF (level==3) THEN
            ! Clouds and precipitation available
            rxt = a_rc
            rnt = CCN
            xrpp = a_rpp
            xnpp = a_npp
        ELSE
            ! Levels 4 and 5
            rxt = a_rc
            rnt = SUM(a_ncloudp,DIM=4)
            xrpp = a_srp
            xnpp = a_snrp
        ENDIF
        CALL set_cs_warm(nzp,nxp,nyp,rxt,rnt,xrpp,xnpp,a_theta,dn0,zm,zt,dzm,xt,yt,time)
    ENDIF

  end subroutine statistics
  !
  ! -----------------------------------------------------------------------
  ! subroutines set_cs_warm, cs_rem_set, cs_rem_save and set_cs_any:
  ! write (and compute) column average statistics
  !
  ! Save named data (already available)
  subroutine set_cs_any(n2,n3,r,nam)
    use netcdf
    integer, intent(in) :: n2,n3
    REAL, INTENT(IN) :: r(n2,n3)
    CHARACTER (LEN=*) :: nam
    INTEGER :: iret, VarID

    IF (csflg) THEN
        iret = nf90_inq_varid(ncid3,trim(nam),VarID)
        IF (iret==NF90_NOERR) THEN
            iret = nf90_put_var(ncid3, VarID, r(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
            iret = nf90_sync(ncid3)
        ENDIF
    ENDIF

  END subroutine set_cs_any
  !
  ! Removal statistics (level>3): calculate values for further use
  SUBROUTINE cs_rem_set(n2,n3,n4,raer,rcld,rprc,rice,rsnw)

    USE mo_submctl, ONLY : nbins, ncld, nprc, nice,  nsnw
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n2,n3,n4   ! Grid dimensions
    REAL, INTENT(in) :: raer(n2,n3,n4*nbins), & ! Removal arrays
                               rcld(n2,n3,n4*ncld), &
                               rprc(n2,n3,n4*nprc), &
                               rice(n2,n3,n4*nice), &
                               rsnw(n2,n3,n4*nsnw)

    INTEGER :: si, i, end,str

    IF (.NOT.csflg) RETURN

    ! Calculate all removal fluxes and save those to scs_rm for later use
    i=1
    DO si = 1,nvar_spec3+1 ! +1 for water
        ! Removal by sedimentation of aerosol
        str = (si-1)*nbins+1
        end = si*nbins
        scs_rm(i,:,:) = SUM(raer(:,:,str:end),DIM=3)
        i=i+1

        ! Removal by sedimentation of cloud droplets
        str = (si-1)*ncld+1
        end = si*ncld
        scs_rm(i,:,:) = SUM(rcld(:,:,str:end),DIM=3)
        i=i+1

        ! Removal by precipitation
        str = (si-1)*nprc+1
        end = si*nprc
        scs_rm(i,:,:) = SUM(rprc(:,:,str:end),DIM=3)
        i=i+1

        IF (level>4) THEN
            ! Removal by sedimentation of ice particles
            str = (si-1)*nice+1
            end = si*nice
            scs_rm(i,:,:) = SUM(rice(:,:,str:end),DIM=3)
            i=i+1

            ! Removal by snow
            str = (si-1)*nsnw+1
            end = si*nsnw
            scs_rm(i,:,:) = SUM(rsnw(:,:,str:end),DIM=3)
            i=i+1
        ENDIF
    ENDDO

  END SUBROUTINE cs_rem_set
  !
  ! Removal statistics (level>3): save values
  SUBROUTINE cs_rem_save(n2,n3)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: n2,n3

    INTEGER :: si, i, end,str
    CHARACTER(LEN=3) :: nam

    IF (.NOT.csflg) RETURN

    ! Save all previously calculated removal fluxes
    !   Note: fluxes not calculated during spinup, so saving zeros
    i=1
    DO si = 1,nvar_spec3+1 ! +1 for water
        IF (si==nvar_spec3+1) THEN
            nam='H2O'
        else
            nam=spec3(si)
        ENDIF

        ! Removal by sedimentation of aerosol
        CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//nam//'dr') ! 'dr' should be for aerosol and 'ae' for water
        i=i+1

        ! Removal by sedimentation of cloud droplets
        CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//nam//'cl')
        i=i+1

        ! Removal by precipitation
        CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//nam//'pr')
        i=i+1

        IF (level>4) THEN
            ! Removal by sedimentation of ice particles
            !CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//name//'ic')
            i=i+1

            ! Removal by snow
            !CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//nam//'sn')
            i=i+1
        ENDIF
    ENDDO

  END SUBROUTINE cs_rem_save
  !
  ! Calculate warm cloud statistics
  subroutine set_cs_warm(n1,n2,n3,rc,nc,rp,np,th,dn0,zm,zt,dzm,xt,yt,time)

    use netcdf
    integer :: iret, n, VarID

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: rc(n1,n2,n3),nc(n1,n2,n3),rp(n1,n2,n3),np(n1,n2,n3),th(n1,n2,n3)
    real, intent(in)    :: dn0(n1),zm(n1),zt(n1),dzm(n1),xt(n2),yt(n3),time
    REAL :: lwp(n2,n3), ncld(n2,n3), rwp(n2,n3), nrain(n2,n3), zb(n2,n3), zc(n2,n3), &
                th1(n2,n3), lmax(n2,n3)
    INTEGER :: ncloudy(n2,n3), nrainy(n2,n3)
    integer :: i, j, k
    real    :: bf(n1), cld, rn, sval, dmy

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
                nrain(i,j)=ncld(i,j)+nc(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
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

    iret = nf90_sync(ncid3)
    nrec3 = nrec3 + 1

  end subroutine set_cs_warm
  !
  ! -----------------------------------------------------------------------
  ! subroutine set_ts: computes and writes time sequence stats
  !
  subroutine set_ts(n1,n2,n3,w,th,dn0,zt,zm,dzt,dzm,th00,time)

    use defs, only : cp

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: w(n1,n2,n3),th(n1,n2,n3)
    real, intent(in)    :: dn0(n1),zt(n1),zm(n1),dzt(n1),dzm(n1),th00,time

    integer :: k
    real    :: bf(n1)

    ssclr(1) = time
    ssclr(4) = get_zi(n1, n2, n3, 2, th, dzm, zt, 1.)   ! maximum gradient
    ssclr(5) = get_zi(n1, n2, n3, 3, th, thvar, zt, 1.) ! maximum variance
    !
    ! buoyancy flux statistics
    !
    ssclr(7) = 0.
    do k = 2,n1-2
       bf(k) = wtv_res(k) + wtv_sgs(k)
       ssclr(7) = ssclr(7) + (tke_res(k)+tke_sgs(k))*dn0(k)/dzt(k)
       svctr(k,33) = svctr(k,33) + wtv_sgs(k)*9.8/th00
    end do
    ssclr(6) = get_zi(n1, n2, n3, 4, th, bf, zm, 1.) ! minimum buoyancy flux

    ssclr(8) = bf(2)
    ssclr(9) = maxval(w)

    ssclr(12) = ssclr(12)*cp*(dn0(1)+dn0(2))*0.5

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
    ssclr(14) = get_zi(n1, n2, n3, 1, q, dzm, zt, 0.5e-3)

  end subroutine ts_lvl1
  !
  ! -----------------------------------------------------------------------
  ! subroutine ts_lvl2: computes and writes time sequence stats
  !
  subroutine ts_lvl2(n1,n2,n3,rt,rs,zt)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: rt(n1,n2,n3),rs(n1,n2,n3), zt(n1)

    integer :: k,i,j
    real    :: cpnt, unit, xaqua

    ssclr(18)  = zt(n1)
    ssclr(19)  = 0.
    ssclr(20)  = 0.
    ssclr(28)  = 0.

    unit = 1./real((n2-4)*(n3-4))
    do j=3,n3-2
       do i=3,n2-2
          cpnt  = 0.
          do k=2,n1-2
             xaqua = rt(k,i,j) - rs(k,i,j)
             if (xaqua > 1.e-5) then
                ssclr(17) = max(ssclr(17),zt(k))
                ssclr(18) = min(ssclr(18),zt(k))
                cpnt = unit
                ssclr(20) = max(ssclr(20), xaqua)
                ssclr(28) = ssclr(28) + 1.
             end if
          end do
          ssclr(19) = ssclr(19) + cpnt
       end do
    end do

    if (ssclr(18) == zt(n1)) ssclr(18) = -999.

  end subroutine ts_lvl2
 !
  ! -----------------------------------------------------------------------
  ! subroutine ts_lvl4: computes and writes time sequence stats of Salsa variables --
  !  Implemented by Zubair Maalick 20/07/2015
  !  Some rewriting and adjusting by Juha Tonttila
  !
  SUBROUTINE ts_lvl4(n1,n2,n3,rc)
    use mo_submctl, only : nlim
    USE grid, ONLY : prtcl, bulkNumc, bulkMixrat,dzt
    USE class_componentIndex, ONLY : IsUsed

    IMPLICIT NONE

    integer, intent(in) :: n1,n2,n3
    REAL, INTENT(in) :: rc(n1,n2,n3)

    REAL :: a0(n1,n2,n3), a1(n1,n2,n3)
    integer :: ii,ss
    LOGICAL :: cond_ic(n1,n2,n3), cond_oc(n1,n2,n3)
    CHARACTER(len=3), PARAMETER :: zspec(7) = (/'SO4','OC ','BC ','DU ','SS ','NH ','NO '/)

    CALL bulkNumc('cloud','a',a0)
    CALL bulkNumc('cloud','b',a1)
    cond_ic(:,:,:) = ( a0(:,:,:) + a1(:,:,:) > nlim .AND. rc(:,:,:) > 1.e-5 )
    cond_oc = .NOT. cond_ic

    ssclr_b(1) = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_ic)
    CALL bulkNumc('aerosol','a',a0)
    CALL bulkNumc('aerosol','b',a1)
    ssclr_b(2) = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_ic)
    ssclr_b(3) = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_oc)

    ii = 4
    DO ss = 1,7

       IF (IsUsed(prtcl,zspec(ss))) THEN
          CALL bulkMixrat(zspec(ss),'cloud','a',a0)
          CALL bulkMixrat(zspec(ss),'cloud','b',a1)
          ssclr_b(ii) = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_ic)
          ii = ii + 1

          CALL bulkMixrat(zspec(ss),'aerosol','a',a0)
          CALL bulkMixrat(zspec(ss),'aerosol','b',a1)

          ssclr_b(ii) = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_ic)
          ii = ii + 1

          ssclr_b(ii) = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_oc)
          ii = ii + 1
       ELSE
          ii = ii + 3
       END IF
    END DO

  END SUBROUTINE ts_lvl4
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
  subroutine accum_rad(n1,n2,n3,rflx,sflx,alb)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: rflx(n1,n2,n3)
    real, optional, intent (in) :: sflx(n1,n2,n3), alb(n2,n3)

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

  end subroutine accum_rad
  !
  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_LVL1: Accumulates various statistics over an
  ! averaging period for moisture variable (smoke or total water)
  !
  subroutine accum_lvl1(n1,n2,n3,rt)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)  :: rt(n1,n2,n3)

    integer :: k
    real    :: a1(n1),a2(n1),a3(n1)

    call get_avg3(n1,n2,n3,rt,a1)
    call get_var3(n1,n2,n3,rt,a1,a2)
    CALL get_3rd3(n1,n2,n3,rt,a1,a3)

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
  subroutine accum_lvl2(n1, n2, n3, th00, dn0, zm, w, th, tl, &
       rl, rs, rt)

    use defs, only : ep2

    integer, intent (in) :: n1,n2,n3
    real, intent (in)                       :: th00
    real, intent (in), dimension(n1)        :: zm, dn0
    real, intent (in), dimension(n1,n2,n3)  :: w, th, tl, rl, rs, rt

    real, dimension(n1,n2,n3) :: tv    ! Local variable
    integer                   :: k, i, j, kp1
    real, dimension(n1)       :: a1, a2, a3, tvbar
    real, dimension(n1,n2,n3) :: scr, xy1, xy2, tlw, tvw, rtw
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
             if (rt(k,i,j) > rs(k,i,j) + 0.01e-3) then
                xy1(k,i,j)=1.
                if (tv(k,i,j) > tvbar(k)) THEN
                   xy2(k,i,j)=1.
                end if
                !
                tlw(k,i,j)=(.5*(tl(k,i,j)+tl(kp1,i,j))+th00)*w(k,i,j)
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
    CALL get_avg3(n1,n2,n3,tl+th00,a1,cond=cond)
    svctr(:,67)=svctr(:,67)+a1(:)
    CALL get_avg3(n1,n2,n3,tv,a1,cond=cond)
    svctr(:,68)=svctr(:,68)+a1(:)
    CALL get_avg3(n1,n2,n3,rt,a1,cond=cond)
    svctr(:,69)=svctr(:,69)+a1(:)
    CALL get_avg3(n1,n2,n3,rl,a1,cond=cond)
    svctr(:,70)=svctr(:,70)+a1(:)
    CALL get_avg3(n1,n2,n3,tlw,a1,cond=cond)
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
    CALL get_avg3(n1,n2,n3,tl+th00,a1,cond=cond)
    svctr(:,77)=svctr(:,77)+a1(:)
    CALL get_avg3(n1,n2,n3,tv,a1,cond=cond)
    svctr(:,78)=svctr(:,78)+a1(:)
    CALL get_avg3(n1,n2,n3,rt,a1,cond=cond)
    svctr(:,79)=svctr(:,79)+a1(:)
    CALL get_avg3(n1,n2,n3,rl,a1,cond=cond)
    svctr(:,80)=svctr(:,80)+a1(:)
    CALL get_avg3(n1,n2,n3,tlw,a1,cond=cond)
    svctr(:,81)=svctr(:,81)+a1(:)
    CALL get_avg3(n1,n2,n3,tvw,a1,cond=cond)
    svctr(:,82)=svctr(:,82)+a1(:)
    CALL get_avg3(n1,n2,n3,rtw,a1,cond=cond)
    svctr(:,83)=svctr(:,83)+a1(:)

    !
    ! liquid water path
    !
    do j=3,n3-2
       do i=3,n2-2
          scr(1,i,j) = 0.
          do k=2,n1
             scr(1,i,j)=scr(1,i,j)+rl(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
          enddo
       end do
    end do
    ssclr(15) = get_avg2dh(n2,n3,scr(1,:,:))
    scr(1,:,:)=(scr(1,:,:)-ssclr(15))**2 ! For LWP variance
    ssclr(16) = get_avg2dh(n2,n3,scr(1,:,:))
  end subroutine accum_lvl2
  !
  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_LVL3: Accumulates specialized statistics that depend
  ! on level 3 variables.
  !
  subroutine accum_lvl3(n1, n2, n3, dn0, zm, rc, rr, nr, rrate, CCN)

    use defs, only : alvl

    integer, intent (in) :: n1,n2,n3
    real, intent (in)                      :: CCN
    real, intent (in), dimension(n1)       :: zm, dn0
    real, intent (in), dimension(n1,n2,n3) :: rc, rr, nr, rrate

    integer                :: k, i, j
    real                   :: nrsum, nrcnt, rrsum, rrcnt
    real                   :: rmax, rmin
    real, dimension(n1)    :: a1
    real, dimension(n2,n3) :: scr2
    REAL :: mask(n1,n2,n3), tmp(n1,n2,n3)

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
    do k=2,n1
       rrsum = 0.
       rrcnt = 0.
       do j=3,n3-2
          do i=3,n2-2
             ! RWP
             scr2(i,j)=scr2(i,j)+rr(k,i,j)*dn0(k)*(zm(k)-zm(k-1))

             ! Rainy grid cell
             if (rr(k,i,j) > 0.001e-3) then
                nrsum = nrsum + nr(k,i,j)
                nrcnt = nrcnt + 1.
             end if

             ! Precipitating grid cell
             if (rrate(k,i,j) > 3.65e-5) then
                rrsum = rrsum + rrate(k,i,j)
                rrcnt = rrcnt + 1.
             end if
          end do
       end do
       if (k == 2 ) ssclr(24) = rrcnt/REAL( (n3-4)*(n2-4) )
    end do
    ssclr(22) = get_avg2dh(n2,n3,scr2)
    scr2(:,:) = rrate(2,:,:)
    ssclr(23) = get_avg2dh(n2,n3,scr2)
    ssclr(25) = CCN
    IF (nrcnt>0.) ssclr(26) = nrsum/nrcnt
    ssclr(27) = nrcnt

  end subroutine accum_lvl3

  !---------------------------------------------------------------------
  ! SUBROUTINE ACCUM_LVL4: Accumulates specialized statistics that depend
  ! on level 4 variables.
  !
  subroutine accum_lvl4(n1,n2,n3)
    use mo_submctl, only : in1a,in2b,fn2a,fn2b, &
                               ica,fca,icb,fcb,ira,fra, &
                               nprc,nlim,prlim
    use grid, ONLY : bulkNumc, bulkMixrat, meanRadius, binSpecMixrat, &
                     a_rc, a_srp, a_rp, a_rh, prtcl,    &
                     a_naerop, a_ncloudp, a_nprecpp
    USE class_ComponentIndex, ONLY : IsUsed

    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3
    INTEGER :: ii,ss,k,bb

    LOGICAL :: cloudmask(n1,n2,n3)
    LOGICAL :: drizzmask(n1,n2,n3)

    REAL :: a0
    REAL, DIMENSION(n1,n2,n3)           :: a1,a12
    REAL, DIMENSION(n1,5)               :: a2
    REAL, DIMENSION(n1,fn2a)            :: a3_a
    REAL, DIMENSION(n1,fn2b-fn2a)       :: a3_b
    REAL, DIMENSION(n1,fca%cur)         :: a4_a
    REAL, DIMENSION(n1,fcb%cur-fca%cur) :: a4_b
    REAL, DIMENSION(n1,nprc)            :: a5

    CHARACTER(len=3), PARAMETER :: zspec(7) = (/'SO4','OC ','BC ','DU ','SS ','NH ','NO '/)


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

    ! In cloud mask
    WHERE (a1+a12 > nlim)
       cloudmask = .TRUE.
    ELSEWHERE
       cloudmask = .FALSE.
    END WHERE

    CALL bulkNumc('precp','a',a1)
    CALL get_avg3(n1,n2,n3,a1,a2(:,5))

    ! In drizzle mask
    WHERE (a1 > prlim)
       drizzmask = .TRUE.
    ELSEWHERE
       drizzmask = .FALSE.
    END WHERE

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
    DO bb = in1a,fn2a
       CALL get_avg3(n1,n2,n3,a_naerop(:,:,:,bb),a3_a(:,bb))
    END DO
    DO bb = in2b,fn2b
       CALL get_avg3(n1,n2,n3,a_naerop(:,:,:,bb),a3_b(:,bb-fn2a))
    END DO
    DO bb = ica%cur,fca%cur
       CALL get_avg3(n1,n2,n3,a_ncloudp(:,:,:,bb),a4_a(:,bb))
    END DO
    DO bb = icb%cur,fcb%cur
       CALL get_avg3(n1,n2,n3,a_ncloudp(:,:,:,bb),a4_b(:,bb-fca%cur))
    END DO
    DO bb = ira,fra
       CALL get_avg3(n1,n2,n3,a_nprecpp(:,:,:,bb),a5(:,bb))
    END DO

    svctr_aa(:,:,1) = svctr_aa(:,:,1) + a3_a(:,:)
    svctr_ab(:,:,1) = svctr_ab(:,:,1) + a3_b(:,:)
    svctr_ca(:,:,1) = svctr_ca(:,:,1) + a4_a(:,:)
    svctr_cb(:,:,1) = svctr_cb(:,:,1) + a4_b(:,:)
    svctr_p(:,:,1) = svctr_p(:,:,1) + a5(:,:)


    IF (.TRUE.) THEN

    ! Species mixing ratios
    ! -------------------------------------------
    ii = 16
    DO ss = 1,7
       IF (IsUsed(prtcl,zspec(ss))) THEN
          ! Total mass mixing ratios
          CALL bulkMixrat(zspec(ss),'aerosol','a',a1)
          CALL bulkMixrat(zspec(ss),'aerosol','b',a12)
          CALL get_avg3(n1,n2,n3,a1+a12,a2(:,1))

          ! Binned mixing ratios
          DO bb = in1a,fn2a
             CALL binSpecMixrat('aerosol',zspec(ss),bb,a1) ! z,x,y-field for bin bb
             CALL get_avg3(n1,n2,n3,a1,a3_a(:,bb))         ! average profile for bin bb for species ss
          END DO
          DO bb = in2b,fn2b
             CALL binSpecMixrat('aerosol',zspec(ss),bb,a1) ! z,x,y-field for bin bb
             CALL get_avg3(n1,n2,n3,a1,a3_b(:,bb-fn2a))    ! average profile for bin bb for species ss
          END DO

          ! In-cloud
          CALL bulkMixrat(zspec(ss),'cloud','a',a1)
          CALL bulkMixrat(zspec(ss),'cloud','b',a12)
          CALL get_avg3(n1,n2,n3,a1+a12,a2(:,2),cond=cloudmask)

          ! Binned mixing ratios
          DO bb = ica%cur,fca%cur
             CALL binSpecMixrat('cloud',zspec(ss),bb,a1)  ! z,x,y-field for bin bb
             CaLL get_avg3(n1,n2,n3,a1,a4_a(:,bb),cond=cloudmask)        ! average profile for bin bb for species ss
          END DO
          DO bb = icb%cur,fcb%cur
             CALL binSpecMixrat('cloud',zspec(ss),bb,a1)  ! z,x,y-field for bin bb
             CALL get_avg3(n1,n2,n3,a1,a4_b(:,bb-fca%cur),cond=cloudmask)! average profile for bin bb for species ss
          END DO

          ! In-drizzle
          CALL bulkMixrat(zspec(ss),'precp','a',a1)
          CALL get_avg3(n1,n2,n3,a1,a2(:,3),cond=drizzmask)

          ! Binned mixing ratios
          DO bb = 1,nprc
             CALL binSpecMixrat('precp',zspec(ss),bb,a1)  ! z,x,y-field for bin bb
             CALL get_avg3(n1,n2,n3,a1,a5(:,bb),cond=drizzmask)          ! average profile for bin bb for species ss
          END DO

          svctr_b(:,ii:ii+2) = svctr_b(:,ii:ii+2) + a2(:,1:3)
          svctr_aa(:,:,ss+1) = svctr_aa(:,:,ss+1) + a3_a(:,:)
          svctr_ab(:,:,ss+1) = svctr_ab(:,:,ss+1) + a3_b(:,:)
          svctr_ca(:,:,ss+1) = svctr_ca(:,:,ss+1) + a4_a(:,:)
          svctr_cb(:,:,ss+1) = svctr_cb(:,:,ss+1) + a4_b(:,:)
          svctr_p(:,:,ss+1) = svctr_p(:,:,ss+1) + a5(:,:)

       END IF ! IsUsed

       ii = ii + 3

    END DO ! ss

    END IF

    ! Liquid water mixing ratio
    CALL get_avg3(n1,n2,n3,a_rc,a2(:,1))

    ! Precipitation mixing ratio
    CALL get_avg3(n1,n2,n3,a_srp,a2(:,2))

    ! Water vapor mixing ratio
    CALL get_avg3(n1,n2,n3,a_rp,a2(:,3))

    ! Relative humidity
    CALL get_avg3(n1,n2,n3,a_rh,a2(:,4))

    svctr_b(:,37:40) = svctr_b(:,37:40) + a2(:,1:4)

  end subroutine accum_lvl4

  !
  !
  !
  subroutine comp_tke(n1,n2,n3,dzm,th00,u,v,w,s)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dzm(n1),th00,u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
    real, intent (inout) :: s(n1,n2,n3)

    integer :: k,kp1,i,j
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
    USE grid, ONLY : level

    integer :: iret, n, VarID

    do n=1,nv1
       iret = nf90_inq_varid(ncid1, s1(n), VarID)
       iret = nf90_put_var(ncid1, VarID, ssclr(n), start=(/nrec1/))
       ssclr(n) = 0.
    end do

    IF (level >= 4) THEN
       DO n = 1,nv1sbulk
          iret = nf90_inq_varid(ncid1, s1SalsaBulk(n), VarID)
          IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
          iret = nf90_put_var(ncid1, VarID, ssclr_b(n), start=(/nrec1/))
          ssclr_b(n) = 0.
       END DO
    END IF

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
    USE mo_submctl, ONLY : in1a,in2b,fn2a,fn2b,fca,ica,fcb,icb,fra,ira, &
                               aerobins,cloudbins,precpbins

    integer, intent (in) :: n1
    real, intent (in)    :: time
    real, intent (in)    :: dn0(n1), u0(n1), v0(n1), zm(n1), zt(n1)

    integer :: iret, VarID, k, n, kp1

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

       svctr(k,10:nv2) = svctr(k,10:nv2)/nsmp
       IF (level >=4 ) THEN
          svctr_b(k,:) = svctr_b(k,:)/nsmp
          svctr_aa(k,:,:) = svctr_aa(k,:,:)/nsmp
          svctr_ab(k,:,:) = svctr_ab(k,:,:)/nsmp
          svctr_ca(k,:,:) = svctr_ca(k,:,:)/nsmp
          svctr_cb(k,:,:) = svctr_cb(k,:,:)/nsmp
          svctr_p(k,:,:) = svctr_p(k,:,:)/nsmp
       END IF

    end do

    iret = nf90_inq_VarID(ncid2, s2(1), VarID)
    iret = nf90_put_var(ncid2, VarID, time, start=(/nrec2/))
    if (nrec2 == 1) then
       iret = nf90_inq_varid(ncid2, s2(2), VarID)
       iret = nf90_put_var(ncid2, VarID, zt, start = (/nrec2/))
       iret = nf90_inq_varid(ncid2, s2(3), VarID)
       iret = nf90_put_var(ncid2, VarID, zm, start = (/nrec2/))
       iret = nf90_inq_varid(ncid2, s2(4), VarID)
       iret = nf90_put_var(ncid2, VarID, dn0, start = (/nrec2/))
       ! Juha: For SALSA
       IF (level >= 4) THEN
          iret = nf90_inq_varid(ncid2,s2SalsaBulk(1),VarID)
          iret = nf90_put_var(ncid2,VarID,aerobins(in1a:fn2a),start=(/nrec2/))
          iret = nf90_inq_varid(ncid2,s2SalsaBulk(2),VarID)
          iret = nf90_put_var(ncid2,VarID,aerobins(in2b:fn2b),start=(/nrec2/))
          iret = nf90_inq_varid(ncid2,s2SalsaBulk(3),VarID)
          iret = nf90_put_var(ncid2,VarID,cloudbins(ica%cur:fca%cur),start=(/nrec2/))
          iret = nf90_inq_varid(ncid2,s2SalsaBulk(4),VarID)
          iret = nf90_put_var(ncid2,VarID,cloudbins(icb%cur:fcb%cur),start=(/nrec2/))
          iret = nf90_inq_varid(ncid2,s2SalsaBulk(5),VarID)
          iret = nf90_put_var(ncid2,VarID,precpbins(ira:fra),start=(/nrec2/))
       END IF
       ! \\ SALSA
       iret = nf90_inq_varid(ncid2, s2(5), VarID)
       iret = nf90_put_var(ncid2, VarID, u0, start = (/nrec2/))
       iret = nf90_inq_varid(ncid2, s2(6), VarID)
       iret = nf90_put_var(ncid2, VarID, v0, start = (/nrec2/))
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
       ! Bulk SALSA
       DO n = 6,nv2sbulk
          iret = nf90_inq_varid(ncid2,s2SalsaBulk(n),VarID)
          IF (iret /= NF90_NOERR) CYCLE ! This is intended to mainly keep track of which variables are used
          iret = nf90_put_var(ncid2,VarID,svctr_b(:,n), start=(/1,nrec2/),  &
               count=(/n1,1/))
       END DO

       ! Binned aerosols, regime a
       DO n = 1,nv2saa
          iret = nf90_inq_varid(ncid2,s2Aeroa(n),VarID)
          IF (iret /= NF90_NOERR) CYCLE ! keep track of used vars
          iret = nf90_put_var(ncid2,VarID,svctr_aa(:,:,n), start=(/1,1,nrec2/),  &
               count=(/n1,fn2a,1/))
       END DO

       ! Aerosols, regime b
       DO n = 1,nv2sab
          iret = nf90_inq_varid(ncid2,s2Aerob(n),VarID)
          IF (iret /= NF90_NOERR) CYCLE ! Keep track of used variables
          iret = nf90_put_var(ncid2,VarID,svctr_ab(:,:,n), start=(/1,1,nrec2/),  &
               count=(/n1,fn2b-fn2a,1/))
       END DO

       ! Cloud droplets, regime a
       DO n = 1,nv2sca
          iret = nf90_inq_varid(ncid2,s2Clouda(n),VarID)
          IF (iret /= NF90_NOERR) CYCLE ! Keep track of used variables
          iret = nf90_put_var(ncid2,VarID,svctr_ca(:,:,n), start=(/1,1,nrec2/),  &
               count=(/n1,fca%cur,1/))
       END DO

       ! Cloud droplets, regime b
       DO n = 1,nv2scb
          iret = nf90_inq_varid(ncid2,s2Cloudb(n),VarID)
          IF (iret /= NF90_NOERR) CYCLE ! Keep track of used variables
          iret = nf90_put_var(ncid2,VarID,svctr_cb(:,:,n), start=(/1,1,nrec2/),  &
               count=(/n1,fcb%cur-fca%cur,1/))
       END DO

       ! Precipitation
       DO n = 1,nv2sp
          iret = nf90_inq_varid(ncid2,s2Precp(n),VarID)
          IF (iret /= NF90_NOERR) CYCLE ! Keep track of used variables
          iret = nf90_put_var(ncid2,VarID,svctr_p(:,:,n), start=(/1,1,nrec2/),   &
               count=(/n1,fra,1/))
       END DO

    END IF

    iret  = nf90_sync(ncid2)
    nrec2 = nrec2+1
    nsmp  = 0.

    do k=1,n1
       svctr(k,:) = 0.
       IF (level >= 4) THEN
          svctr_b(k,:) = 0.
          svctr_aa(k,:,:) = 0.
          svctr_ab(k,:,:) = 0.
          svctr_ca(k,:,:) = 0.
          svctr_cb(k,:,:) = 0.
          svctr_p(k,:,:) = 0.
       END IF
    end do

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

    svctr(:,23)=svctr(:,23)+v1(:)/float((n2-2)*(n3-2))
    svctr(:,25)=svctr(:,25)+v2(:)/float((n2-2)*(n3-2))
    svctr(:,27)=svctr(:,27)+v3(:)/float((n2-2)*(n3-2))

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
       case (1)
          nn = 87
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
    USE grid, ONLY : prtcl
    USE mo_submctl, ONLY : nbins, ncld, nprc, nice,  nsnw
    USE class_componentIndex, ONLY : IsUsed, GetIndex
    IMPLICIT NONE

    INTEGER, INTENT(in)           :: n2,n3,n4                     ! Grid dimensions
    REAL, INTENT(in)              :: raer(n2,n3,n4*nbins)        ! Array containing the binned 2d-field
    REAL, OPTIONAL, INTENT(in)    :: rcld(n2,n3,n4*ncld), &
                                     rprc(n2,n3,n4*nprc), &     ! 2 optional arrays for calculating total removals
                                     rice(n2,n3,n4*nice), &
                                     rsnw(n2,n3,n4*nsnw)

    CHARACTER(len=3), PARAMETER :: zspec(7) = (/'SO4','OC ','BC ','DU ','SS ','NH ','NO '/)
    INTEGER, PARAMETER          :: sstrt(7) = (/28,    33,   38,   43,   48,   53,   58/)

    REAL :: zavg

    INTEGER :: ss, si
    INTEGER :: tt
    INTEGER :: end,str

    ! Removal of water first
    si = GetIndex(prtcl,'H2O')
    ! Aerosols
    zavg = 0.
    tt = 25
    str = (si-1)*nbins+1
    end = si*nbins
    zavg = get_avg2dh( n2,n3,SUM(raer(:,:,str:end),DIM=3) )
    ssclr_b(tt) = ssclr_b(tt) + zavg

    ! Cloud
    zavg = 0.
    tt = 26
    str = (si-1)*ncld+1
    end = si*ncld
    zavg = get_avg2dh( n2,n3,SUM(rcld(:,:,str:end),DIM=3) )
    ssclr_b(tt) = ssclr_b(tt) + zavg

    ! Precipitation
    zavg = 0.
    tt = 27
    str = (si-1)*nprc+1
    end = si*nprc
    zavg = get_avg2dh( n2,n3,SUM(rprc(:,:,str:end),DIM=3) )
    ssclr_b(tt) = ssclr_b(tt) + zavg

    ! Ice
    zavg = 0.
    tt = 28
    str = (si-1)*nice+1
    end = si*nice
    zavg = get_avg2dh( n2,n3,SUM(rice(:,:,str:end),DIM=3) )
    ssclr_b(tt) = ssclr_b(tt) + zavg

    ! Precipitation
    zavg = 0.
    tt = 29
    str = (si-1)*nsnw+1
    end = si*nsnw
    zavg = get_avg2dh( n2,n3,SUM(rsnw(:,:,str:end),DIM=3) )
    ssclr_b(tt) = ssclr_b(tt) + zavg

    DO ss = 1,7

       IF ( IsUsed(prtcl,zspec(ss)) ) THEN

          si = GetIndex(prtcl,zspec(ss))

          ! Dry removal
          zavg = 0.
          tt = sstrt(ss)
          str = (si-1)*nbins+1
          end = si*nbins
          zavg = get_avg2dh( n2,n3,SUM(raer(:,:,str:end),DIM=3) )
          ssclr_b(tt) = ssclr_b(tt) + zavg

          ! Removal by sedimentation of cloud droplets
          zavg = 0.
          tt = sstrt(ss) + 1
          str = (si-1)*ncld+1
          end = si*ncld
          zavg = get_avg2dh( n2,n3,SUM(rcld(:,:,str:end),DIM=3) )
          ssclr_b(tt) = ssclr_b(tt) + zavg

          ! Removal by precipitation
          zavg = 0.
          tt = sstrt(ss) + 2
          str = (si-1)*nprc+1
          end = si*nprc
          zavg = get_avg2dh( n2,n3,SUM(rprc(:,:,str:end),DIM=3) )
          ssclr_b(tt) = ssclr_b(tt) + zavg

          ! Total wet removal
          zavg = 0.
          tt = sstrt(ss) + 3
          str = (si-1)*ncld+1
          end = si*ncld
          zavg = zavg + get_avg2dh( n2,n3,SUM(rcld(:,:,str:end),DIM=3) )
          str = (si-1)*nprc+1
          end = si*nprc
          zavg = zavg + get_avg2dh( n2,n3,SUM(rprc(:,:,str:end),DIM=3) )
          ssclr_b(tt) = ssclr_b(tt) + zavg

          ! Removal by sedimentation of ice particles
          zavg = 0.
          tt = sstrt(ss) + 4
          str = (si-1)*nice+1
          end = si*nice
          zavg = get_avg2dh( n2,n3,SUM(rice(:,:,str:end),DIM=3) )
          ssclr_b(tt) = ssclr_b(tt) + zavg

          ! Removal by snowing
          zavg = 0.
          tt = sstrt(ss) + 5
          str = (si-1)*nsnw+1
          end = si*nsnw
          zavg = get_avg2dh( n2,n3,SUM(rsnw(:,:,str:end),DIM=3) )
          ssclr_b(tt) = ssclr_b(tt) + zavg

          ! Total icy removal
          zavg = 0.
          tt = sstrt(ss) + 6
          str = (si-1)*nice+1
          end = si*nice
          zavg = zavg + get_avg2dh( n2,n3,SUM(rice(:,:,str:end),DIM=3) )
          str = (si-1)*nsnw+1
          end = si*nsnw
          zavg = zavg + get_avg2dh( n2,n3,SUM(rsnw(:,:,str:end),DIM=3) )
          ssclr_b(tt) = ssclr_b(tt) + zavg

          ! Total removal
          zavg = 0.
          tt = sstrt(ss) + 7
          str = (si-1)*nbins+1
          end = si*nbins
          zavg = zavg + get_avg2dh( n2,n3,SUM(raer(:,:,str:end),DIM=3) )
          str = (si-1)*ncld+1
          end = si*ncld
          zavg = zavg + get_avg2dh( n2,n3,SUM(rcld(:,:,str:end),DIM=3) )
          str = (si-1)*nprc+1
          end = si*nprc
          zavg = zavg + get_avg2dh( n2,n3,SUM(rprc(:,:,str:end),DIM=3) )
          str = (si-1)*nice+1
          end = si*nice
          zavg = zavg + get_avg2dh( n2,n3,SUM(rice(:,:,str:end),DIM=3) )
          str = (si-1)*nsnw+1
          end = si*nsnw
          zavg = zavg + get_avg2dh( n2,n3,SUM(rprc(:,:,str:end),DIM=3) )
          ssclr_b(tt) = ssclr_b(tt) + zavg

       END IF

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
    USE mo_submctl, ONLY : rhowa
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
          IF ( .NOT. PRESENT(rv) .OR. .NOT. PRESENT(rc) .OR. .NOT. PRESENT(prc) ) &
               STOP 'acc_massbudget (stat): ERROR - for atm water q,rc and prc must be present'

          a1 = 0.; a2 = 0.; a3 = 0.
          a1 = get_avg_ts(n1,n2,n3,rv*dn,dz)
          a2 = get_avg_ts(n1,n2,n3,rc*dn,dz)
          a3 = get_avg_ts(n1,n2,n3,prc*dn,dz)
          massbdg(1) = (a1 + a2 + a3)

       CASE(1)
          IF ( .NOT. PRESENT(rv) .OR. .NOT. PRESENT(rc) .OR. .NOT. PRESENT(prc) ) &
               STOP 'acc_massbudget (stat): ERROR - for atm water q,rc and prc must be present'

          a1 = 0.; a2 = 0.; a3 = 0.
          a1 = get_avg_ts(n1,n2,n3,rv*dn,dz)
          a2 = get_avg_ts(n1,n2,n3,rc*dn,dz)
          a3 = get_avg_ts(n1,n2,n3,prc*dn,dz)
          massbdg(2) = (a1 + a2 + a3) ! Not accumulated

       CASE(2)
          IF ( .NOT. PRESENT(revap) .OR. .NOT. PRESENT(ApVdom) ) &
               STOP 'acc_massbudget (stat): ERROR - for evaporation stats revap must be present'

          a1 = 0.
          a1 = get_avg2dh(n2,n3,revap)
          massbdg(3) = massbdg(3) + a1*tstep*ApVdom

       CASE(3)
          IF ( .NOT. PRESENT(rdep) .OR. .NOT. PRESENT(ApVdom) ) &
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

