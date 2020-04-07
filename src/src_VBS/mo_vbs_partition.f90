!****************************************************************
! Module mo_vbs_partition: subroutines for VOC oxidation and vapor-liquid partitioning
!
! To be used with UCLALES-SALSA only. VBS species are produced when gas phase
! volatile organic vapors are oxidized. Vapor-liquid partitioning is calculated
! for all semi-volatile species including VBS and aqSOA species and simple
! sulfate and organic vapors.
!
! Code developers
!  Thomas. Kuehn (UEF) - initial setup
!  Joonas Merikanto (FMI) - aqueous phase SOA code
!  Tomi Raatikainen (FMI) - modified for UCLALES-SALSA (2019)
!
!****************************************************************

MODULE mo_vbs_partition

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
       vbs_gas_phase_chem,& ! VOC oxidation
       vbs_condensation     ! Vapor-liquid partitioning

CONTAINS

  SUBROUTINE vbs_gas_phase_chem(kbdim, klev,&
       ptemp, tstep, etime, pc_gas, ngas)

    ! -----------------------------------------------------------------------
    !
    ! SUBROUTINE vbs_gas_phase_chem
    !
    ! Computes the production of VBS species due to oxidation of VOC
    ! precursors in the gas phase. Considered oxidants are OH, O3 and NO3
    !
    ! Authors:
    ! --------
    ! Declan O'Donnell, MPI-Met
    ! Kai Zhang, MPI-Met
    ! Thomas Kuehn, UEF, 6/2015 --
    ! Joonas Merikanto, FMI 7/2015
    !
    ! -----------------------------------------------------------------------

    ! use statements
    USE mo_submctl, ONLY : id_oh, id_o3, id_no3, &
        model_lat, start_doy, avog

    USE mo_vbsctl, ONLY: &
         t_vbs_group, t_voc_prec, t_aq_soa, & ! the VBS, VOC and aqSOA data structures
         vbs_set, aqsoa_set,vbs_voc_set,   &
         aqsoa_ngroup, vbs_ngroup, vbs_nvocs
		 
    USE mo_vbs, ONLY : rate_o3_o1d_ave, maxdayfac, zenith

    ! input/output parameters
    INTEGER,  INTENT(in) :: kbdim                    ! geographic block max number of locations
    INTEGER,  INTENT(in) :: klev                     ! numer of levels
    INTEGER,  INTENT(in) :: ngas                     ! number of gases
    REAL, INTENT(in) :: ptemp (kbdim,klev)           ! air temperature [K]
    REAL, INTENT(inOUT) :: pc_gas(kbdim,klev,ngas)   ! gas phase concentrations [mol/m3]
    REAL, INTENT(IN) :: tstep, etime                 ! time step and elapsed time [s]

    ! local variables
    TYPE(t_voc_prec), POINTER :: zprec  ! local copy for easier reading
    TYPE(t_vbs_group), POINTER :: zbase ! local copy for easier reading
    TYPE(t_aq_soa), POINTER :: zaqsoa   ! local copy for easier reading

    ! Gas phase reactant concentrations
    REAL ::     &
         zc_oh(kbdim, klev),  & ! OH concentration at t+dt
         zc_o3(kbdim, klev),  & ! O3 concentration at t+dt
         zc_no3(kbdim, klev), & ! NO3 concentration at t+dt
         zc_prec,             & ! VOC precursor concentation at t+dt
         zc_aqsoa_prec          ! aqsoa concentation at t+dt

    ! Gas phase reaction rates
    REAL ::     &
         zr_rate_oh_prec,   & ! voc precursor reaction rate with OH [m3/(mol*s)]
         zc_rate_oh_prec,   & ! voc precursor conversion rate due to OH [1/s]
         zr_rate_o3_prec,   & ! voc precursor reaction rate with O3 [m3/(mol*s)]
         zc_rate_o3_prec,   & ! voc precursor conversion rate due to O3 [1/s]
         zr_rate_no3_prec,  & ! voc precursor reaction rate with NO3 [m3/(mol*s)]
         zc_rate_no3_prec,  & ! voc precursor conversion rate due to NO3 [1/s]
         zc_rate_tot_prec,  & ! voc precursor total conversion rate
         zr_rate_oh_aqsoa,  & ! aqsoa reaction rate with OH [m3/(mol*s)]
         zc_rate_oh_aqsoa,  & ! aqsoa conversion rate due to OH [1/s]
         zr_rate_o3_aqsoa,  & ! aqsoa reaction rate with O3 [m3/(mol*s)]
         zc_rate_o3_aqsoa,  & ! aqsoa conversion rate due to O3 [1/s]
         zr_rate_no3_aqsoa, & ! aqsoa reaction rate with NO3 [m3/(mol*s)]
         zc_rate_no3_aqsoa    ! aqsoa conversion rate due to NO3 [1/s]

    ! zloss/zproduction rates
    REAL :: &
         zloss_prec,           & ! total zloss of precursor VOC
         zloss_aqsoa,          & ! total zloss of aqsoa
         zprodrate_aqsoa,      & ! zproduction rate for aqsoa
         zlossoh_aqsoa,        & ! OH zloss of aqsoa
         zlosso3_aqsoa,        & ! O3 zloss of aqsoa
         zlossno3_aqsoa,       & ! NO3 zloss of aqsoa
         zlossphoto_aqsoa        ! photodissociation zloss of aqsoa

    REAL :: zdayfac, u0, znightfac
    INTEGER :: jk, jl, jb, jv
    INTEGER :: idt_prec, idt_base, idt_aqsoa

    ! -----------------------------------------------------------------------
    ! executable procedure
    ! -----------------------------------------------------------------------

    ! For LES time scales it is either day or night (cos(zenith_angle) < 0)
    ! Cosine of the solar zenith angle
    u0=zenith(model_lat, start_doy+etime/86400.)
    !zdayfac=1.
    IF (u0<=0.)then
    	zdayfac=0.
    else
        zdayfac = 6.073e-5 * u0**1.743 * exp(-0.474 / u0)
        zdayfac = zdayfac / rate_o3_o1d_ave
		znightfac = (maxdayfac - zdayfac) / (maxdayfac - 1.)
    endif

    ! Obtain oxidant concentrations and scale for day/night time length. Also,
    ! convert mol/m^3 to #/cm^3, because the unit of prefactor k0 is cm^3/#/s.
    zc_oh(:,:)=0.
    zc_o3(:,:)=0.
    zc_no3(:,:)=0.
    IF (id_oh>0)  zc_oh(:,:)=pc_gas(:,:,id_oh)*zdayfac*avog*1e-6
    IF (id_o3>0)  zc_o3(:,:)=pc_gas(:,:,id_o3)*avog*1e-6   !*zdayfac
    IF (id_no3>0) zc_no3(:,:)=pc_gas(:,:,id_no3)*znightfac*avog*1e-6

    DO jv = 1,vbs_nvocs

       ! marking the precursor we currently process
       zprec => vbs_voc_set(jv)
       idt_prec = zprec%id_gas

       DO jk = 1,klev
          DO jl = 1,kbdim
!print *, jv, pc_gas(jl,jk,idt_prec), pc_gas(jl,jk,id_oh), pc_gas(jl,jk,id_o3), pc_gas(jl,jk,id_no3), zdayfac !,zprec%k_0_OH,zprec%k_0_O3,zprec%k_0_NO3,ptemp(jl,jk)
              IF (pc_gas(jl,jk,idt_prec)<1e-20) CYCLE

              ! the precursor concentration in this grid cell [mol/m3]
              zc_prec = pc_gas(jl,jk,idt_prec)

              !---reaction rates k=k_0*exp(E_p/T) [cm^3/s]
              zr_rate_oh_prec  = zprec%k_0_OH  * EXP(zprec%Eact_p_OH/ptemp(jl,jk))
              zr_rate_o3_prec  = zprec%k_0_O3  * EXP(zprec%Eact_p_O3/ptemp(jl,jk))
              zr_rate_no3_prec = zprec%k_0_NO3 * EXP(zprec%Eact_p_NO3/ptemp(jl,jk))

              !---conversion rates k_x*c [1/s]
              zc_rate_oh_prec  = zr_rate_oh_prec  * zc_oh(jl,jk)
              zc_rate_o3_prec  = zr_rate_o3_prec  * zc_o3(jl,jk)
              zc_rate_no3_prec = zr_rate_no3_prec * zc_no3(jl,jk)

              !---total conversion rate [1/s]
              zc_rate_tot_prec = zc_rate_oh_prec + zc_rate_o3_prec + zc_rate_no3_prec

              !---amount of reacted voc [mol/m3]
              zloss_prec = zc_prec * (1. - EXP(-zc_rate_tot_prec*tstep))
			  
              !---correcting concentration [mol/m3]
              pc_gas(jl,jk,idt_prec) = pc_gas(jl,jk,idt_prec) - zloss_prec

              !---calculating the production of products
              DO jb = 1,vbs_ngroup
                  ! marking the product we currently process
                 zbase => vbs_set(jb)
                 idt_base = vbs_set(jb)%id_gas

                 !---correcting concentration
                 pc_gas(jl,jk,idt_base) = pc_gas(jl,jk,idt_base) + zloss_prec * zprec%stoich_coeff(jb)
              END DO !vbs_ngroup

              DO jb = 1 , aqsoa_ngroup
                 zaqsoa => aqsoa_set(jb)
                 idt_aqsoa = aqsoa_set(jb)%id_gas

                 ! the zproduction rate for jb
                 zprodrate_aqsoa = zc_prec*(1.-EXP(-zc_rate_oh_prec*tstep))/tstep  &
                              * zprec%stoich_coeff(jb+vbs_ngroup)
                 if(jb.eq.2.and.jv.eq.2) zprodrate_aqsoa =  zprodrate_aqsoa +  & !Hard coded for ISOP+O3 -> 0.01 GLYX
                            zc_prec*(1.-EXP(-zc_rate_o3_prec*tstep))/tstep &
                             * zprec%stoich_coeff(jb+vbs_ngroup)*0.5
                 if(jb.eq.2) zprodrate_aqsoa = zprodrate_aqsoa + 0.24*zlossoh_aqsoa/tstep !hard coded for EIPOX+OH->0.24*Glyoxal

                 pc_gas(jl,jk,idt_aqsoa) = pc_gas(jl,jk,idt_aqsoa)+ zprodrate_aqsoa*tstep

                 ! the gas phase loss rates for aqSOA
                 !concentration
                 zc_aqsoa_prec=pc_gas(jl,jk,idt_aqsoa)
                 !---reaction rates k=k_0*exp(E_p/T); E_p=E/R
                 zr_rate_oh_aqsoa   = zaqsoa%k_0_OH  * EXP(zaqsoa%Eact_p_OH /ptemp(jl,jk))
                 zr_rate_o3_aqsoa   = zaqsoa%k_0_O3  * EXP(zaqsoa%Eact_p_O3 /ptemp(jl,jk))
                 zr_rate_no3_aqsoa  = zaqsoa%k_0_NO3 * EXP(zaqsoa%Eact_p_NO3/ptemp(jl,jk))

                 zc_rate_oh_aqsoa  = zr_rate_oh_aqsoa  * zc_oh(jl,jk)
                 zc_rate_o3_aqsoa  = zr_rate_o3_aqsoa  * zc_o3(jl,jk)
                 zc_rate_no3_aqsoa = zr_rate_no3_aqsoa * zc_no3(jl,jk)

                 zlossoh_aqsoa    = zc_aqsoa_prec * (1. - EXP(-zc_rate_oh_aqsoa*tstep))
                 zlosso3_aqsoa    = zc_aqsoa_prec * (1. - EXP(-zc_rate_o3_aqsoa*tstep))
                 zlossno3_aqsoa   = zc_aqsoa_prec * (1. - EXP(-zc_rate_no3_aqsoa*tstep))
                 zlossphoto_aqsoa = zc_aqsoa_prec * (1. - EXP(-zaqsoa%photodis *tstep))*(1.-zdayfac)
                 zloss_aqsoa = zlossoh_aqsoa + zlosso3_aqsoa + zlossno3_aqsoa + zlossphoto_aqsoa

                 pc_gas(jl,jk,idt_aqsoa) = pc_gas(jl,jk,idt_aqsoa) - zloss_aqsoa
              END DO
          END DO !kbdim
       END DO !klev
  END DO !vbs_nvocs

  END SUBROUTINE vbs_gas_phase_chem

  SUBROUTINE vbs_condensation(kbdim,  klev,       &
            paero,  pcloud, pprecp, pice,  psnow, &
            pc_gas, ngas,   ptemp,  ppres, ptstep )
    ! SALSA hydrometeors
    USE mo_submctl, ONLY : t_section, nbins, ncld, nprc, nice, nsnw
    IMPLICIT NONE
    !-- Input and output variables ----------
    INTEGER, INTENT(IN) :: kbdim, klev, ngas
    REAL, INTENT(IN) :: ptemp(kbdim,klev), ppres(kbdim,klev), ptstep
    REAL, INTENT(INOUT) :: pc_gas(kbdim,klev,ngas)      ! gas concentrations [mol/m3]
    TYPE(t_section), INTENT(INOUT) :: paero(kbdim,klev,nbins), pcloud(kbdim,klev,ncld), &
        pprecp(kbdim,klev,nprc), pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    !-- Local variables ----------------------
    INTEGER :: ii, jj

    DO jj = 1,klev
        DO ii = 1,kbdim
            ! This calculates equilibrium for VBS bins, aqSOA, a non-volatile organic vapor and sulfate
            call vbs_condensation_salsa( &
                ptemp(ii,jj),ppres(ii,jj),ptstep,pc_gas(ii,jj,:),ngas, &
                paero(ii,jj,:),pcloud(ii,jj,:),pprecp(ii,jj,:), &
                pice(ii,jj,:),psnow(ii,jj,:) )
        END DO ! kbdim
    END DO ! klev

  END SUBROUTINE vbs_condensation


  SUBROUTINE vbs_condensation_salsa(&
       ptemp,ppres,                         &
       pt_step,                             &
       pc_gas,ngas,                         &
       paero,pcloud,pprecp,pice,psnow       &
       )

    ! -----------------------------------------------------------------------
    !
    ! SUBROUTINE vbs_condensation_salsa
    !
    ! Computes the exchange of volatility basis set (VBS) and optionally
    ! also aqSOA, non-volatile organic vapor and sulfate groups between
    ! gas- and particle phases and their distribution among the aerosol
    ! size bins. Changes are only calculated for one grid cell, looping
    ! over grid and height happens in the calling routine ( *condensation*
    ! in *mo_salsa_dynamics* ).
    !
    ! Authors:
    ! --------
    ! Thomas Kühn, UEF    6/2015 --
    ! Joonas Merikanto, FMI    6/2015 --
    ! Tomi Raatikainen, FMI    11/2019 --
    !
    ! vbs_condensation_salsa is called
    ! from *condensation* in *mo_salsa_dynamics*
    !
    ! Related literature:
    ! -Jacobson (2005), Fundamentals of Atmospheric Modelling, 2nd Edition
    ! -Farina et. al (2010), Modeling global secondary organic aerosol
    !  formation and processing with the volatility basis set: Implications
    !  for anthropogenic secondary organic aerosol
    ! -Nguyen, T. B., et al. "Organic aerosol formation from the reactive uptake of isoprene
    !  epoxydiols (IEPOX) onto non-acidified inorganic seeds." Atmospheric Chemistry and Physics
    !  14.7 (2014): 3497-3510.
    !
    ! -----------------------------------------------------------------------

    USE mo_vbsctl, ONLY: &
         t_vbs_group, vbs_ngroup, vbs_set, & ! VBS data structure
         t_aq_soa, aqsoa_set, aqsoa_ngroup   ! aqSOA data structure

    USE mo_submctl, ONLY : &
         pstand, pi, rg, surfw0,               & ! Constants
         dens, mws, rhowa, rhoic, rhosn,       & ! Physics properties
         laqsoa,                               & ! aqSOA partitioning (optional)
         ioc, iocg, moc, part_ocnv,            & ! NVOA partitioning (optional)
         iso, isog, msu, part_h2so4,           & ! Sulfate partitioning (optional)
         nspec, nbins, ncld, nprc, nice, nsnw, & ! Number of species and bins
         t_section, nlim, prlim,               & ! Data type for the bin representation and concentration tresholds
         CalcDimension                           ! Function for updating wet sizes

    ! -----------------------------------------------------------------------

    ! input / output variables

    INTEGER, INTENT(IN) :: ngas ! number of gases

    REAL, INTENT(IN) :: &
         ptemp,                                & ! ambient temperature [K]
         ppres,                                & ! ambient pressure [K]
         pt_step                                 ! timestep [s]

    REAL, INTENT(INOUT) :: &
         pc_gas(ngas)     ! Gas phase tracer concentrations [mol/m3]

    TYPE(t_section), INTENT(INOUT) :: &
         paero(nbins),pcloud(ncld),pprecp(nprc), & ! Aerosol, cloud, rain, ice and snow properties: wet diameters [m]
         pice(nice),psnow(nsnw)                    ! and bin number [#/m3] and per species volume concentrations [m3/m3]

    ! -----------------------------------------------------------------------

    ! local variables

    ! indexing VBS and aqsoa groups
    TYPE(t_vbs_group), POINTER :: zgroup
    TYPE(t_aq_soa), POINTER :: zgroupa

    ! intermediates for calculation
    REAL :: &
        zeps = 1.e-30,                                     & ! small number
        zc_part_all(nbins+ncld+nprc+nice+nsnw,nspec+1),    & ! particle/droplet phase concentrations per species [mol/m3]
        zk_diam(nbins+ncld+nprc+nice+nsnw),                & ! particle/droplet diameter [m]
        zk_numc(nbins+ncld+nprc+nice+nsnw),                & ! particle/droplet number concentration [#/m3]
        zk_mass(nbins+ncld+nprc+nice+nsnw),                & ! mass transfer coefficient [1/s]
        zk_mass_tot,                                       & ! total mass transfer coefficient
        zs_prime(vbs_ngroup,nbins+ncld+nprc+nice+nsnw),    & ! Kelvin effect
        zknud(nbins+ncld+nprc+nice+nsnw),                  & ! Knudsen number
        zbeta(nbins+ncld+nprc+nice+nsnw),                  & ! transitional correction factor
        zc_sat(vbs_ngroup),                                & ! uncorrected saturation vapor concentration [mol/m3]
        zc_sat_part(nbins+ncld+nprc+nice+nsnw),            & ! mole-fraction corrected sat. vap. conc. per class
        zc_gas(vbs_ngroup),                                & ! current gas phase concentration [mol/m3]
        zc_part(vbs_ngroup,nbins+ncld+nprc+nice+nsnw),     & ! currect particle phase concentration [mol/m3]
        zc_tot(vbs_ngroup)                                   ! total concentration [mol/m3]

    ! we divide the time step into increasing
    ! sub-intervals, with
    ! sum(dt_n) = pt_step and
    ! dt_(n+1)=dt_n * ztime_fac
    INTEGER, PARAMETER :: zn_steps = 10          ! amount of sub-time steps
    REAL,PARAMETER :: ztime_fac = 1.5     ! sub-time step growth factor
    REAL ::           zh                     ! current sub-time step length

    ! loop indices
    INTEGER :: &
         jt,                                   & ! time sub-step
         jc,                                   & ! class
         jvc,                                  & ! vbs class
         jg                                      ! vbs group

    ! other
    INTEGER :: nbin_tot, nbin_aq, nbin_aer, iaqv, nn
    REAL :: ztemp_sum, zdfvap, zmfp, p_lwc, p_awc, p_cwc, Heff, zdvap, zcvap_new

    ! -----------------------------------------------------------------------
    ! executable procedure
    ! -----------------------------------------------------------------------

    ! Total number of species in volc
    nn = nspec+1


    ! ********* Calculating particle phase concentrations [mol/m3] for all species including water *********
    ! Find all bins that have particles, but keep record of the bin types (aerosol,
    ! aerosol + cloud + rain [i.e. typically aqueous], and all bins)
    jvc = 0
    zc_part_all(:,:)=0.
    CALL CalcDimension(nbins,paero,nlim,1) ! Update wet diameters (%dwet)
    DO jc = 1,nbins ! Aerosol
       IF (paero(jc)%numc>nlim) THEN
          ! increment size bin counter
          jvc = jvc + 1
          ! concentrations for all particle phase species [mol/m3]
          zc_part_all(jvc,:)=paero(jc)%volc(1:nn)*dens(1:nn)/mws(1:nn)
          ! aerosol bin number concentration [#/m3]
          zk_numc(jvc) = paero(jc)%numc
          ! aerosol bin wet diameter [m]
          zk_diam(jvc) = paero(jc)%dwet
       END IF
    END DO
    nbin_aer = jvc ! Number of aerosol bins
    CALL CalcDimension(ncld,pcloud,nlim,2)
    DO jc = 1,ncld ! Cloud droplets
       IF (pcloud(jc)%numc>nlim) THEN
          jvc = jvc + 1
          zc_part_all(jvc,:)=pcloud(jc)%volc(1:nn)*dens(1:nn)/mws(1:nn)
          zk_numc(jvc) = pcloud(jc)%numc
          zk_diam(jvc) = pcloud(jc)%dwet
       END IF
    END DO
    CALL CalcDimension(nprc,pprecp,prlim,3)
    DO jc = 1,nprc ! Rain drops
       IF (pprecp(jc)%numc>prlim) THEN
          jvc = jvc + 1
          zc_part_all(jvc,:)=pprecp(jc)%volc(1:nn)*dens(1:nn)/mws(1:nn)
          zk_numc(jvc) = pprecp(jc)%numc
          zk_diam(jvc) = pprecp(jc)%dwet
       END IF
    END DO
    nbin_aq = jvc ! Total number of the commonly water soluble bins
    CALL CalcDimension(nice,pice,prlim,4)
    DO jc = 1,nice ! Ice
       IF (pice(jc)%numc>prlim) THEN
          jvc = jvc + 1
          zc_part_all(jvc,:)=pice(jc)%volc(1:nn)*dens(1:nn)/mws(1:nn)
          zc_part_all(jvc,1)=zc_part_all(jvc,1)*rhoic/rhowa ! Correction for ice density
          zk_numc(jvc) = pice(jc)%numc
          zk_diam(jvc) = pice(jc)%dwet
       END IF
    END DO
    CALL CalcDimension(nsnw,psnow,prlim,5)
    DO jc = 1,nsnw ! Snow
       IF (psnow(jc)%numc>prlim) THEN
          jvc = jvc + 1
          zc_part_all(jvc,:)=psnow(jc)%volc(1:nn)*dens(1:nn)/mws(1:nn)
          zc_part_all(jvc,1)=zc_part_all(jvc,1)*rhosn/rhowa ! Correction for snow density
          zk_numc(jvc) = psnow(jc)%numc
          zk_diam(jvc) = psnow(jc)%dwet
       END IF
    END DO
    nbin_tot=jvc ! Total number of bins


    ! ********* VBS partitioning *********
    ! Note 1: VBS partitioning in the default SALSA is calculated only for the soluble a-bins, but here
    !       solubility is not predefined. Therefore, include all bins and for the known solids (ice and snow)
    !       just ignore Raoult and Kelvin effects.
    ! Note 2: aqSOA partitioning is calculated for all bins that can contain liquid water (aerosol, cloud and rain)

    ! Calculate the mass-transfer coefficient (Eq. 6.64) for a model substance (alpha-pinene) representing all
    ! semi-volatile organics (VBS, aqSOA and non-volatile organics).
    ! Reference:
    !   Mark Z. Jacobson: Fundamentals of Atmospheric Modeling
    !   Second Edition, Cambridge University Press, 2005

    ! Diffusion coefficient [m2/s] for an organic vapor in air: Fuller's method (Eq. 11-4.4 in Reid et al.) gives
    ! D=0.0610519 cm2/s for alpha-pinene (C10H16; M=136 g/mol) at P=1 bar and T=298 K.
    !   Niinemets and Reichstein (Global Biogeochem. Cycles, 16(4), 1110, doi:10.1029/2002GB001927, 2002)
    !   report that alpha-pinene D=5.812e-6 m2/s (T=25 C), which is quite close.
    zdfvap = 0.0610519e-4*(ptemp/298.)**1.75*(1e5/ppres)

    ! Mean free path for alpha-pinene vapor in air (Eq. 16.23) [m]
    zmfp = 64.*zdfvap/(5.*pi)*29./(29.+136.)*sqrt(pi*136e-3/(8.*rg*ptemp))

    ! Knudsen number (Eq. 16.20)
    zknud(1:jvc) = 2.*zmfp/zk_diam(1:jvc)
    ! transitional correction factor (Eq. 16.19), without mass accomodation coefficient
    zbeta(1:jvc) = 1./(1.+zknud(1:jvc)*(1.33+0.71/zknud(1:jvc))/(1.+1./zknud(1:jvc)) )

    ! mass transfer coefficient (Eqs 16.56 and 16.64), without ventilation correction [1/s]
    zk_mass(1:jvc) = zk_numc(1:jvc)*2.*pi*zk_diam(1:jvc)*zdfvap*zbeta(1:jvc)

    ! total mass transfer coefficient [1/s]
    zk_mass_tot = sum(zk_mass(1:jvc))

    ! mapping to VBS
    zs_prime(:,:) = 1. ! Kelvin effect (one for ice and snow)
    DO jg = 1,vbs_ngroup
       ! easier access
       zgroup => vbs_set(jg)

       ! Kelvin effect for the typically aqueous droplets (aerosol, cloud and rain)
       zs_prime(jg,1:nbin_aq) = exp(4.0*surfw0*mws(zgroup%id_vols)/dens(zgroup%id_vols)/(rg*ptemp*zk_diam(1:nbin_aq)))

       ! computing the uncorrected saturation vapor concentration [mol/m3]
       ! after Farina et al (2010), Equation (4)
       zc_sat(jg) = zgroup%C0*(zgroup%T0/ptemp)*exp(zgroup%Hvap_eff*(1./zgroup%T0-1./ptemp))

       ! copying the particle volume concentrations [mol/m3]
       zc_part(jg,:) = zc_part_all(:,zgroup%id_vols)

       ! gas phase concentration [mol/m3]
       zc_gas(jg) = pc_gas(zgroup%id_gas)

       ! total concentration [mol/m3]
       zc_tot(jg) = zc_gas(jg)+sum(zc_part(jg,:))
    END DO ! jg

    ! Starting the sub-time step iteration
    zh=pt_step*(ztime_fac-1.)/(ztime_fac**real(zn_steps)-1.)

    DO jt = 1, zn_steps ! loop over sub-time steps
       ! the particle phase concentrations [mol/m3] and
       ! the gas phase concentrations [mol/m3] before condensation:
       DO jg = 1,vbs_ngroup

          ! if concentrations are too small, we skip
          IF (zc_tot(jg) < zeps) CYCLE

          ! easier access
          zgroup => vbs_set(jg)

          ! computing the new mole-fraction weighted, per-bin equilibrium vapor concentration
          zc_sat_part(:) = zc_sat(jg) ! x=1 for ice and snow
          DO jvc = 1,nbin_aq
             ! Note: assuming an organic phase, so ignore all other than VBS species (POA could be included).
             IF (sum(zc_part(:,jvc)) > 0.0) &
                zc_sat_part(jvc) = zc_sat(jg)*zc_part(jg,jvc)/sum(zc_part(:,jvc))
          END DO ! jvc

          ! computing the new gas phase concentration (16.71)
          zc_gas(jg) = min(&
               (zc_gas(jg)+zh*sum(zk_mass(1:nbin_tot)*zs_prime(jg,1:nbin_tot)*zc_sat_part(1:nbin_tot)))/(1.0+zh*zk_mass_tot),&
               zc_tot(jg)&
               )

          ! computing the new particle phase concentrations and correcting for negative
          ! concentrations (16.69) + (16.72a)
          DO jvc = 1, nbin_tot
             zc_part(jg,jvc)=max(&
                  zc_part(jg,jvc)+zh*zk_mass(jvc)*(zc_gas(jg)-zs_prime(jg,jvc)*zc_sat_part(jvc)),&
                  0.&
                  )
          END DO ! jvc

          ! correcting for too high particle concentrations
          ztemp_sum = sum(zc_part(jg,:))
          IF (ztemp_sum < zeps) THEN
             zc_part(jg,:) = 0.0
             zc_gas(jg)    = zc_tot(jg)
          ELSE
             zc_part(jg,:) = zc_part(jg,:)*(zc_tot(jg)-zc_gas(jg))/ztemp_sum
          END IF
       END DO ! jg

       ! increasing the time sub-step size
       zh = zh*ztime_fac
    END DO ! jt

    ! mapping back to gas and particles
    DO jg = 1,vbs_ngroup
       ! easier access
       zgroup => vbs_set(jg)

       ! gas phase concentration [mol/m3]
       pc_gas(zgroup%id_gas) = zc_gas(jg)

       ! particle phase concentration [mol/m3]
       zc_part_all(:,zgroup%id_vols) = zc_part(jg,:)
    END DO ! jg


    ! ********* aqSOA partitioning *********
    ! NOTE: Currently very simple scheme based on the use of effective Henry's coefficient
    !       Can be updated later on

    ! Total relevant liquid water content [mol/m3]
    p_lwc=sum(zc_part_all(1:nbin_aq,1))

    IF (laqsoa .AND. nbin_aq > 0 .AND. p_lwc>1e-30) THEN
       ! Aerosol and cloud/rain liquid water contents [mol/m3]
       p_awc=sum(zc_part_all(1:nbin_aer,1))
       p_cwc=max(0.,p_lwc-p_awc)

       DO jc = 1,aqsoa_ngroup
          ! easier access
          zgroupa => aqsoa_set(jc)
          iaqv = zgroupa%id_vols ! Index to volc

          ! Total aqsoa concentration (gas+liquid phase) [mol/m3]
          ztemp_sum = pc_gas(zgroupa%id_gas) + sum(zc_part_all(1:nbin_aq,iaqv))

          ! First calculate new gas phase concentration according to total aerosol water content and
          ! "effective Henry's constant", ref. Nguyen et al. (2014)
          !   aerosol and cloud/rain water weighted Henry's constant
          Heff = (p_awc*zgroupa%Eff_henry_aerosol + p_cwc*zgroupa%Eff_henry_cloud)/p_lwc ! [M/atm]
          pc_gas(zgroupa%id_gas) = ztemp_sum/(1.0+Heff*0.08205736*ptemp*p_lwc*mws(1)/dens(1))

          ! Then distribute remaining aqsoa according to aerosol water content [mol/m3]
          zc_part_all(1:nbin_aq,iaqv) = (ztemp_sum-pc_gas(zgroupa%id_gas))*zc_part_all(1:nbin_aq,1)/p_lwc
       END DO
    END IF


    ! ********* non-volatile organic vapor and sulfate partitioning *********
    !   - simplification of the VBS partitioning
    !   - no nucleation
    !   - for details, see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling,
    !     Second Edition (2005), equation (16.73)
    IF (pc_gas(iocg) > 1.e-20 .and. zk_mass_tot>1e-30 .and. part_ocnv) THEN
        ! Assume the same chemical composition (alpha-pinene) as above in the VBS partitioning
        ! => no need to update mass trasfer coeffcients

        zcvap_new = pc_gas(iocg)/(1.+pt_step* zk_mass_tot) ! new gas phase concentration [mol/m3]
        zdvap = pc_gas(iocg) - zcvap_new                   ! change in gas concentration [mol/m3]
        pc_gas(iocg) = zcvap_new                           ! updating vapour concentration [mol/m3]

        ! New concentration (mol/m3) in aerosol, cloud droplets, rain drops, ice and snow
        zc_part_all(:,ioc)=zc_part_all(:,ioc) + zk_mass(:)/zk_mass_tot*zdvap
    END IF

    IF (pc_gas(isog) > 1.e-20 .and. zk_mass_tot > 1.e-30 .and. part_h2so4) THEN
        ! Update parameters for sulfate
        zdfvap = 5.1111e-10*ptemp**1.75*pstand/ppres  ! Diffusion coefficient [m2/s]
        zmfp   = 3.*zdfvap*sqrt(pi*msu/(8.*rg*ptemp)) ! Mean free path [m]

        jvc = nbin_tot
        ! Knudsen number (Eq. 16.20)
        zknud(1:jvc) = 2.*zmfp/zk_diam(1:jvc)
        ! transitional correction factor (Eq. 16.19), without mass accomodation coefficient
        zbeta(1:jvc) = 1./(1.+zknud(1:jvc)*(1.33+0.71/zknud(1:jvc))/(1.+1./zknud(1:jvc)) )
        ! mass transfer coefficient (Eqs 16.56 and 16.64), without ventilation correction [1/s]
        zk_mass(1:jvc) = zk_numc(1:jvc)*2.*pi*zk_diam(1:jvc)*zdfvap*zbeta(1:jvc)
        ! total mass transfer coefficient [1/s]
        zk_mass_tot = sum(zk_mass(1:jvc))

        zcvap_new = pc_gas(isog) /(1.+pt_step*zk_mass_tot) ! new gas phase concentration [mol/m3]
        zdvap = pc_gas(isog) - zcvap_new                   ! change in gas concentration [mol/m3]
        pc_gas(isog) = zcvap_new                           ! updating vapour concentration [mol/m3]

        ! New concentration (mol/m3) in aerosol, cloud droplets, rain drops, ice and snow
        zc_part_all(:,iso)=zc_part_all(:,iso) + zk_mass(:)/zk_mass_tot*zdvap
    END IF



    ! ********* update particle phase volume concentrations (m3/m3) *********
    jvc = 0
    DO jc = 1,nbins ! Aerosol
       IF (paero(jc)%numc>nlim) THEN
          ! increment size bin counter
          jvc = jvc + 1
          ! update VBS species
          paero(jc)%volc(1:nn) = zc_part_all(jvc,1:nn)*mws(1:nn)/dens(1:nn)
       END IF
    END DO
    DO jc = 1,ncld ! Cloud droplets
       IF (pcloud(jc)%numc>nlim) THEN
          jvc = jvc + 1
          pcloud(jc)%volc(1:nn) = zc_part_all(jvc,1:nn)*mws(1:nn)/dens(1:nn)
       END IF
    END DO
    DO jc = 1,nprc ! Rain drops
       IF (pprecp(jc)%numc>prlim) THEN
          jvc = jvc + 1
          pprecp(jc)%volc(1:nn) = zc_part_all(jvc,1:nn)*mws(1:nn)/dens(1:nn)
       END IF
    END DO
    DO jc = 1,nice ! Ice
       IF (pice(jc)%numc>prlim) THEN
          jvc = jvc + 1
          pice(jc)%volc(1:nn) = zc_part_all(jvc,1:nn)*mws(1:nn)/dens(1:nn)
       END IF
    END DO
    DO jc = 1,nsnw ! Snow
       IF (psnow(jc)%numc>prlim) THEN
          jvc = jvc + 1
          psnow(jc)%volc(1:nn) = zc_part_all(jvc,1:nn)*mws(1:nn)/dens(1:nn)
       END IF
    END DO

  END SUBROUTINE vbs_condensation_salsa

END MODULE mo_vbs_partition
