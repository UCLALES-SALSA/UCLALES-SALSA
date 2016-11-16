MODULE mo_salsa_cloud

  !*********************************************************
  !  MO_SALSA_CLOUD
  !*********************************************************
  !
  ! Purpose: Calculates the number of activated cloud
  ! droplets according to parameterizations by:
  !
  ! Abdul-Razzak et al: "A parameterization of aerosol activation -
  !                      3. Sectional representation"
  !                      J. Geophys. Res. 107, 10.1029/2001JD000483, 2002.
  !                      [Part 3]
  !
  ! Abdul Razzak et al: "A parameterization of aerosol activation -
  !                      1. Single aerosol type"
  !                      J. Geophys. Res. 103, 6123-6130, 1998.
  !                      [Part 1]
  !
  !
  ! Interface:
  ! ----------
  ! Called from main aerosol model
  !
  !
  ! Coded by:
  ! ---------
  ! T. Anttila (FMI)     2007
  ! H. Kokkola (FMI)     2007
  ! A.-I. Partanen (FMI) 2007
  !
  !*********************************************************

CONTAINS

  SUBROUTINE cloud_activation(kproma, kbdim, klev,   &
                              temp,   pres,  rv,     &
                              rs,     w,     paero,  &
                              pcloud, pactd          )

    USE mo_constants,   ONLY : g

    USE mo_submctl, ONLY :                      &
         lsactbase,                                 & ! Activation at cloud base
         lsactintst,                                & ! Activation of interstitial particles
         rg,                                        & ! molar gas constant
                                                      ! [J/(mol K)]
         slim,                                      &
         surfw0,                                    & ! surface tension
                                                      ! of water [J/m2]
         !nbin,                                      & ! number of size bins
                                                      ! in subranges
         nlim,                                      & ! lowest possible particle conc. in a bin [#/m3]
         rhosu, msu,                                & ! properties of compounds
         rhooc, moc,                                &
         rhono, mno,                                &
         rhonh, mnh,                                &
         rhobc, mbc,                                &
         rhoss, mss,                                &
         rhowa, mwa,                                &
         rhodu, mdu,                                &
         pi,                                        &
         pi6,                                       &
         cpa,                                       &
         mair,                                      &
         in1a,in2a,in2b,fn2a, fn2b,            & ! size regime bin indices
         t_section,                                 & ! Data type for cloud/rain drops
         t_parallelbin,                             & ! Data type for mapping indices between parallel size bins
         ncld                                         ! Total number of cloud bins

    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::              &
             kproma,                    & ! number of horiz. grid points
             kbdim,                     & ! dimension for arrays
             klev                       ! number of vertical levels

    REAL, INTENT(in) ::             &
             pres(kbdim,klev),          &
             temp(kbdim,klev),          &
             w(kbdim,klev)

    REAL, INTENT(inout) :: rv(kbdim,klev) ! Water vapor mixing ratio
    REAL, INTENT(in)    :: rs(kbdim,klev) ! Saturation vapor mixing ratio

    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld),  &
                                      paero(kbdim,klev,fn2b)

    ! Properties of newly activate particles
    TYPE(t_section), INTENT(out) :: pactd(kbdim,klev,ncld)

    !! Properties of newly activated interstitial particles (for diagnostics)
    !TYPE(t_section) :: pactd_intst(kbdim,klev,ncld)


    !-- local variables --------------
    INTEGER :: ii, jj, kk             ! loop indices

    INTEGER :: bcrita(kbdim,klev),bcritb(kbdim,klev) ! Index of the critical aerosol bin for regimes a and b

    REAL ::                         &
             sil,                       & ! critical supersaturation
                                          !     at the upper bound of the bin
             siu,                       & !  "  at the lower bound of the bin
             scrit(fn2b),               & !  "  at the center of the bin
             aa(kbdim,klev),                        & ! curvature (Kelvin) effect [m]
             bb,                        & ! solute (Raoult) effect [m3]
             ns(kbdim,klev,fn2b),                  & ! number of moles of solute
             nshi,                      & !             " at the upper bound of the bin
             nslo,                      & !             " at the lower bound of the bin
             ratio,                     & ! volume ratio
             vmiddle(kbdim,klev,fn2b),  & ! volume in the middle of the bin [m3]
             s_max,                     & ! maximum supersaturation
             s_eff,                     & ! effective supersaturation
             x, x1, x2, x3, a1, sum1,   & ! technical variables
             ka1,                       & ! thermal conductivity
             dv1,                       & ! diffusion coefficient
             Gc,                        & ! growth coefficient
             alpha,                     & ! see Abdul-Razzak and Ghan, part 3
             gamma,                     & ! see Abdul-Razzak and Ghan, part 3
             L,                         & ! latent heat of evaporation
             ps,                        & ! saturation vapor pressure of water [Pa]
             khi,                       & ! see Abdul-Razzak and Ghan, part 3
             theta,                     & ! see Abdul-Razzak and Ghan, part 3
             frac(kbdim,klev,fn2b),                & ! fraction of activated droplets in a bin
             ntot,                      & ! total number conc of particles [#/m3]
             dinsol(fn2b),              & ! diameter of the insoluble fraction [m]
             dinsolhi,                  & !    "   at the upper bound of a bin [m]
             dinsollo,                  & !    "   at the lower bound of a bin [m]
             zdcrit(kbdim,klev,fn2b),   & ! critical diameter [m]
             zdcstar(kbdim,klev),                   & ! Critical diameter corresponding to Smax
             !zdcint,                    &  ! Critical diameter for calculating activation of interstitial particles
             !zvcrhi,zvcrlo,             & ! critical volume boundaries for partially activated bins
             zdcrhi(kbdim,klev,fn2b),   & ! Critical diameter at the high end of a bin
             zdcrlo(kbdim,klev,fn2b),   & ! Critical diameter at the low end of a bin
             V,                         & ! updraft velocity [m/s]
             rref,                      & ! reference radius [m]
             A, dmx,           & !
             vlo, k,                    &
             !zvact,               & ! Total volume concentration of the newly activated droplets
             !zcoreact,            & ! Volume concentration of the newly activated dry CCN
             !zvtot,               &
             !zcore,               &
             zrho(8)


    ! ------------------------------------------------------------------
    ! Initialization
    ! ------------------------------------------------------------------
    zrho(1:8) = (/rhosu,rhooc,rhono,rhonh,rhobc,rhodu,rhoss,rhowa/)

    bb = 6.*mwa/(pi*rhowa)             ! Raoult effect [m3]
                                          ! NOTE!
                                          ! bb must be multiplied
                                          ! by the number of moles of
                                          ! solute
    zdcrit(:,:,:) = 0.
    zdcrlo(:,:,:) = 0.
    zdcrhi(:,:,:) = 0.
    frac(:,:,:) = 0.
    bcrita(:,:) = fn2a
    bcritb(:,:) = fn2b

    DO jj = 1,klev    ! vertical grid
       DO ii = 1,kproma ! horizontal grid

          vmiddle(ii,jj,in1a:fn2b) = pi6*paero(ii,jj,in1a:fn2b)%dmid**3
          DO kk = 1,ncld
             pactd(ii,jj,kk)%volc(:) = 0.
             pactd(ii,jj,kk)%numc = 0.
          END DO
          aa(ii,jj) = 4.*mwa*surfw0/(rg*rhowa*temp(ii,jj)) ! Kelvin effect [m]

       END DO
    END DO

    ! Get moles of solute at the middle of the bin
    CALL getSolute(kproma,kbdim,klev,paero,ns)

    ! ----------------------------------------------------------------


    ! -------------------------------------
    ! Interstitial activation
    ! -------------------------------------
    IF ( lsactintst ) THEN

       !write(*,*) 'ACT'
       CALL actInterst(kproma,kbdim,klev,paero,pcloud,rv,rs,aa)

       ! Update the moles of solute after interstitial activation
       CALL getSolute(kproma,kbdim,klev,paero,ns)
    END IF


    ! -----------------------------------
    ! Activation at cloud base
    ! -----------------------------------
    IF ( lsactbase ) THEN

       DO jj = 1,klev    ! vertical grid
          DO ii = 1,kproma ! horizontal grid

             A  = aa(ii,jj) * 1.e6                            !     "         [um]
             x  = 4.*aa(ii,jj)**3/(27.*bb)

             ! Get the critical supersaturation for aerosol bins for necessary summation terms (sum1 & ntot)
             ntot = 0.
             sum1 = 0.

             scrit(in1a:fn2b) = exp(sqrt(x/max(epsilon(1.0),ns(ii,jj,in1a:fn2b)))) - 1.

             !-- sums in equation (8), part 3
             ntot = ntot + SUM(paero(ii,jj,in1a:fn2b)%numc)
             sum1 = sum1 + SUM(paero(ii,jj,in1a:fn2b)%numc/scrit(in1a:fn2b)**(2./3.))

             IF(ntot < nlim) CYCLE
             IF (w(ii,jj) <= 0.) CYCLE
             V  = w(ii,jj)

             !-- latent heat of evaporation [J/kg]
             L     = 2.501e6-2370.*(temp(ii,jj)-273.15)

             !-- saturation vapor pressure of water [Pa]
             a1    = 1.-(373.15/temp(ii,jj))
             ps    = 101325.*                                                 &
                  exp(13.3185*a1-1.976*a1**2-0.6445*a1**3-0.1299*a1**4)

             !-- part 1, eq (11)
             alpha = g*mwa*L/(cpa*rg*temp(ii,jj)**2)-                            &
                  g*mair/(rg*temp(ii,jj))

             !-- part 1, eq (12)
             gamma = rg*temp(ii,jj)/(ps*mwa) &
                  + mwa*L**2/(cpa*pres(ii,jj)*mair*temp(ii,jj))

             !-- diffusivity [m2/s], Seinfeld and Pandis (15.65)
             x1 = pres(ii,jj) / 101325.
             dv1= 1.e-4 * (0.211/x1) * ((temp(ii,jj)/273.)**1.94)

             rref = 10.e-9
             !-- corrected diffusivity, part 1, eq (17)
             ! dv = dv1 / (rref/(rref + deltaV) + (dv1/(rref * alphac)) *        &
             !     SQRT(2.*pi*mwa/(rg*temp(ii,jj))))

             !-- thermal conductivity [J/(m s K)], Seinfeld and Pandis (15.75)
             ka1= 1.e-3 * (4.39 + 0.071 * temp(ii,jj))

             !-- growth coefficient, part 1, eq (16)
             !-- (note: here uncorrected diffusivities and conductivities are used
             !    based on personal communication with H. Abdul-Razzak, 2007)
             Gc = 1./(rhowa*rg*temp(ii,jj)/(ps*dv1*mwa) +                      &
                  L*rhowa/(ka1*temp(ii,jj)) * (L*mwa/(temp(ii,jj)*rg)-1.))

             !-- effective critical supersaturation: part 3, eq (8)
             s_eff = (ntot/sum1)**(3./2.)

             !-- part 3, equation (5)

             theta = ((alpha*V/Gc)**(3./2.))/(2.*pi*rhowa*gamma*ntot)

             !-- part 3, equation (6)
             khi = (2./3.)*aa(ii,jj)*SQRT(alpha*V/Gc)

             !-- maximum supersaturation of the air parcel: part 3, equation (9)
             s_max = s_eff / SQRT(0.5*(khi/theta)**(3./2.)              &
                  + ((s_eff**2)/(theta+3.*khi))**(3./4.))

             !-- Juha: Get the critical diameter corresponding to the maximum supersaturation
             ! ------ NÄISTÄ LASKUISTA PUUTTUU LIUKENEMATTOMAN COREN VAIKUTUS?
             zdcstar = 2.*aa(ii,jj)/(3.*s_max)

             !WRITE(*,*) s_max
             !WRITE(*,*) zdcstar

             DO kk = in1a, fn2b

                IF (paero(ii,jj,kk)%numc < nlim) CYCLE

                !-- moles of solute in particle at the upper bound of the bin
                nshi = ns(ii,jj,kk)*paero(ii,jj,kk)%vratiohi

                !-- critical supersaturation
                sil = exp(sqrt(x/nshi)) - 1.

                IF(s_max < sil) CYCLE

                !-- moles of solute at the lower bound of the bin:
                nslo = ns(ii,jj,kk)*paero(ii,jj,kk)%vratiolo

                !-- critical supersaturation
                siu = exp(sqrt(x/nslo)) - 1.

                !-- fraction of activated in a bin, eq (13), part 3
                frac(ii,jj,kk) = min(1.,log(s_max/sil)/log(siu/sil))

                !-- Critical diameters for each bin and bin edges
                zdcrlo(ii,jj,kk) = 2.*sqrt(3.*nslo*bb/aa(ii,jj))
                zdcrhi(ii,jj,kk) = 2.*sqrt(3.*nshi*bb/aa(ii,jj))
                zdcrit(ii,jj,kk) = 2.*sqrt(3.*ns(ii,jj,kk)*bb/aa(ii,jj))

             END DO ! kk

             ! Find critical bin
             DO kk = in1a,fn2a
                IF (frac(ii,jj,kk) < 1. .AND. frac(ii,jj,kk) > 0.) THEN
                   bcrita(ii,jj) = kk
                   EXIT
                END IF
             END DO
             DO kk = in2b,fn2b
                IF (frac(ii,jj,kk) < 1. .AND. frac(ii,jj,kk) > 0.) THEN
                   bcritb(ii,jj) = kk
                   EXIT
                END IF
             END DO

          END DO ! ii

       END DO ! jj

       CALL activate3(kproma,kbdim,klev,paero,bcrita,bcritb,  &
                      zdcrit, zdcrlo, zdcrhi, zdcstar, pactd  )

    END IF ! lsactbase

 
  END SUBROUTINE cloud_activation

! -----------------------------------------------------------------

  SUBROUTINE getSolute(kproma,kbdim,klev,paero,pns)
    
    USE mo_submctl, ONLY : t_section,nlim,       &
                               in1a,fn1a,            &
                               in2a,fn2a,            &
                               in2b,fn2b,            &
                               rhosu, rhooc, rhobc,  &
                               rhonh, rhono, rhodu,  &
                               rhoss,                &
                               msu, moc, mbc,        &
                               mnh, mno, mdu,        &
                               mss
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kproma,kbdim,klev
    TYPE(t_section), INTENT(IN) :: paero(kbdim,klev,fn2b)
    REAL, INTENT(OUT) :: pns(kbdim,klev,fn2b)

    INTEGER :: ii,jj,kk

    pns = 0.

    DO jj = 1,klev

       DO ii = 1,kproma

          !-- subrange 1a

          !-- calculation of critical superaturation in the middle of the bin

          !-- volume in the middle of the bin

          DO kk = in1a, fn1a
             IF (paero(ii,jj,kk)%numc > nlim) THEN

                !-- number of moles of solute in one particle [mol]
                pns(ii,jj,kk) = (3.*paero(ii,jj,kk)%volc(1)*rhosu/msu  +   &
                     paero(ii,jj,kk)%volc(2)*rhooc/moc)/                      &
                     paero(ii,jj,kk)%numc

             END IF
          END DO

          !-- subrange 2a

          DO kk = in2a, fn2a
             IF (paero(ii,jj,kk)%numc > nlim) THEN

                pns(ii,jj,kk) = (3.*paero(ii,jj,kk)%volc(1)*rhosu/msu  +   &
                     paero(ii,jj,kk)%volc(2) * rhooc/moc +                    &
                     paero(ii,jj,kk)%volc(6) * rhono/mno +                    &
                     paero(ii,jj,kk)%volc(7) * rhonh/mnh +                    &
                     2.*paero(ii,jj,kk)%volc(5) * rhoss/mss)/              &
                     paero(ii,jj,kk)%numc

             END IF
          END DO

          !-- subrange 2b

          DO kk = in2b, fn2b
             IF (paero(ii,jj,kk)%numc > nlim) THEN

                pns(ii,jj,kk) = (3.*paero(ii,jj,kk)%volc(1)*rhosu/msu  +   &
                        paero(ii,jj,kk)%volc(2) * rhooc/moc +                 &
                        paero(ii,jj,kk)%volc(6) * rhono/mno +                 &
                        paero(ii,jj,kk)%volc(7) * rhonh/mnh +                 &
                        2.*paero(ii,jj,kk)%volc(5) * rhoss/mss)/           &
                        paero(ii,jj,kk)%numc

             END IF
          END DO



       END DO

    END DO


  END SUBROUTINE getSolute

! -----------------------------------------------------------------

  SUBROUTINE actInterst(kproma,kbdim,klev,paero,pcloud,prv,prs,paa)
    !
    ! Activate interstitial aerosols if they've reached their critical size
    ! HUOM TEKEE VAAN A-BINIT NYT!!
    !
    !
    ! 1. Formulate the profiles of Dwet within bins
    !      - Get the slopes between adjacent bin mids (known)
    !      - Interpolate Dwet to bin edges
    ! 2. Based on Dcrit and Dwet slopes, estimate the Ddry where Dwet becomes > Dcrit
    ! 3. Formulate the slopes for number concentration
    ! 4. Use the Dry limits from (2) as the integration limits if they are defined

    
    USE mo_submctl, ONLY : t_section,nlim,pi6,ica,fca,icb,fcb, &
                               in1a,fn1a,in2a,fn2a,                &
                               nbins,ncld
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kproma,kbdim,klev
    TYPE(t_section), INTENT(INOUT) :: paero(kbdim,klev,nbins),  &
                                      pcloud(kbdim,klev,ncld)
    REAL, INTENT(IN) :: prv(kbdim,klev),prs(kbdim,klev),    &  ! Water vapour and saturation mixin ratios
                            paa(kbdim,klev)                        ! Coefficient for Kelvin effect

    REAL :: zdcstar,zvcstar   ! Critical diameter/volume corresponding to S_LES
    REAL :: zvcint            ! Integration limit volume based on above
    REAL :: zactvol           ! Total volume of the activated particles

    TYPE(t_section) :: zactd(ncld)   ! Temporary storage for activated particles
    REAL :: Nact, Vact(8)        ! Helper variables for transferring the activated particles

    REAL :: Nmid, Nim1, Nip1     ! Bin number concentrations in current and adjacent bins
    REAL :: dNmid, dNim1, dNip1  ! Density function value of the number distribution for current and adjacent bins

    REAL :: Vmid, Vim1, Vip1     ! Dry particle volume in the middle of the bin
    REAL :: Vlo, Vhi             ! Dry particle volume scaled to bin edges
    REAL :: Vlom1,Vhim1          ! - '' - For adjacent bins
    REAL :: Vlop1,Vhip1          !

    REAL :: Vwmid, Vwim1, Vwip1  ! Wet particle volume in the middle of the bin
    REAL :: Vwlo,Vwhi            ! Wet particle volume at bin edges
    REAL :: Vwlom1, Vwhim1       ! Wet particle volumes for adjacent bins
    REAL :: Vwlop1, Vwhip1       ! - '' -

    REAL :: zs1,zs2           ! Slopes for number concetration distributions within bins

    REAL :: N01,N02           ! Origin values for number distribution slopes
    REAL :: V01,V02           ! Origin values for wet particle volume slopes
    REAL :: Nnorm, Vnorm      ! Normalization factors for number and volume integrals

    REAL :: vcut,vint1,vint2  ! cut volume, integration limit volumes
    LOGICAL :: intrange(4)    ! Logical table for integration ranges depending on the shape of the wet size profile:
                                  ! [Vlo -- vint1][vint1 -- Vmid][Vmid -- vint2][vint1 -- Vhi]
    INTEGER :: cb,ab, ii,jj,ss


    DO jj = 1,klev

       DO ii = 1,kproma

          IF ( prv(ii,jj)/prs(ii,jj) > 1.000 ) THEN

             zactd(1:ncld)%numc = 0.
             DO ss = 1,8
                zactd(1:ncld)%volc(ss) = 0.
             END DO

             ! Determine Dstar == critical diameter corresponding to the host model S
             zdcstar = 2.*paa(ii,jj)/( 3.*( (prv(ii,jj)/prs(ii,jj))-1. ) )
             zvcstar = pi6*zdcstar**3

             ! Loop over cloud droplet (and aerosol) bins
             DO cb = ica%cur, fca%cur
                ab = ica%par + (cb-ica%cur)

                IF ( paero(ii,jj,ab)%numc < nlim) CYCLE

                intrange = .FALSE.

                ! Define some parameters
                Nmid = paero(ii,jj,ab)%numc     ! Number concentration at the current bin center
                Vwmid = SUM(paero(ii,jj,ab)%volc(:))/Nmid  ! Wet volume at the current bin center
                Vmid = SUM(paero(ii,jj,ab)%volc(1:7))/Nmid ! Dry volume at the current bin center
                Vlo = Vmid*paero(ii,jj,ab)%vratiolo        ! Dry vol at low limit
                Vhi = Vmid*paero(ii,jj,ab)%vratiohi        ! Dry vol at high limit

                ! Number concentrations and volumes at adjacent bins (check for sizedistribution boundaries)
                IF ( ab > in1a ) THEN
                   Nim1 = paero(ii,jj,ab-1)%numc
                   IF (Nim1 > nlim) THEN
                      Vim1 = SUM(paero(ii,jj,ab-1)%volc(1:7))/Nim1
                      Vwim1 = SUM(paero(ii,jj,ab-1)%volc(:))/Nim1
                   ELSE
                      Vim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                      Vwim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                   END IF
                   Vlom1 = Vim1*paero(ii,jj,ab-1)%vratiolo
                   Vhim1 = Vim1*paero(ii,jj,ab-1)%vratiohi
                ELSE ! ab == in1a
                   Nim1 = nlim
                   Vim1 = Vlo/2.
                   Vlom1 = 0.
                   Vhim1 = Vlo
                   Vwim1 = Vwmid/3. ! Tää ny o vähä tämmöne
                END IF
                IF ( ab < fn2a ) THEN
                   Nip1 = paero(ii,jj,ab+1)%numc
                   IF (Nip1 > nlim) THEN
                      Vip1 = SUM(paero(ii,jj,ab+1)%volc(1:7))/Nip1
                      Vwip1 = SUM(paero(ii,jj,ab+1)%volc(:))/Nip1
                   ELSE
                      Vip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                      Vwip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                   END IF
                   Vlop1 = Vip1*paero(ii,jj,ab+1)%vratiolo
                   Vhip1 = Vip1*paero(ii,jj,ab+1)%vratiohi
                ELSE ! ab == fn2a
                   Nip1 = nlim
                   Vip1 = Vhi + 0.5*(Vhi-Vlo)
                   Vlop1 = Vhi
                   Vhip1 = Vhi + (Vhi-Vlo)
                   Vwip1 = Vhip1  ! ....
                END IF

                ! Keeping thins smooth...
                Vip1 = MAX(Vhi,Vip1)
                Vim1 = MIN(Vlo,Vim1)

                ! First, make profiles of particle wet radius in
                ! order to determine the integration boundaries
                zs1 = (Vwmid - MAX(Vwim1,0.))/(Vmid - Vim1)
                zs2 = (MAX(Vwip1,0.) - Vwmid)/(Vip1 - Vmid)

                ! Get the origin values for slope equations
                V01 = Vwmid - zs1*Vmid
                V02 = Vwmid - zs2*Vmid

                ! Get the wet sizes at bins edges
                Vwlo = MAX(V01 + zs1*Vlo, 0.)
                Vwhi = MAX(V02 + zs2*Vhi, 0.)

                ! Find out dry vol integration boundaries based on *zvcstar*:
                IF ( zvcstar < Vwlo .AND. zvcstar < Vwmid .AND. zvcstar < Vwhi ) THEN
                   ! Whole bin activates
                   vint1 = Vlo
                   vint2 = Vhi

                   intrange(1:4) = .TRUE.

                ELSE IF ( zvcstar > Vwlo .AND. zvcstar > Vwmid .AND. zvcstar > Vwhi) THEN
                   ! None activates
                   vint1 = 999.
                   vint2 = 999.

                   intrange(1:4) = .FALSE.

                ELSE
                   ! Partial activation:
                   ! Slope1
                   vcut = (zvcstar - V01)/zs1  ! Where the wet size profile intersects the critical size (slope 1)
                   IF (vcut < Vlo .OR. vcut > Vmid) THEN
                      ! intersection volume outside the current size range -> set as the lower limit
                      vint1 = Vlo
                   ELSE
                      vint1 = vcut
                   END IF

                   ! Slope2
                   vcut = (zvcstar - V02)/zs2  ! Where the wet size profile intersects the critical size (slope 2)
                   IF (vcut < Vmid .OR. vcut > Vhi) THEN
                      ! Intersection volume outside the current size range -> set as the lower limit
                      vint2 = Vmid
                   ELSE
                      vint2 = vcut
                   END IF

                   ! Determine which size ranges have wet volume larger than the critical
                   intrange(1) = ( Vwlo > zvcstar )
                   intrange(2) = ( Vwmid > zvcstar )
                   intrange(3) = ( Vwmid > zvcstar )
                   intrange(4) = ( Vwhi > zvcstar )

                END IF
                ! DONE WITH INTEGRATION LIMITS
                ! -------------------------------------

                ! Number concentration profiles within bins and integration for number of activated:
                ! -----------------------------------------------------------------------------------
                ! get density distribution values for number concentration
                dNim1 = Nim1/(Vhim1-Vlom1)
                dNip1 = Nip1/(Vhip1-Vlop1)
                dNmid = Nmid/(Vhi-Vlo)

                ! Get slopes
                zs1 = ( dNmid - dNim1 )/( Vmid - Vim1 )
                zs2 = ( dNip1 - dNmid )/( Vip1 - Vmid )

                N01 = dNmid - zs1*Vmid  ! Origins
                N02 = dNmid - zs2*Vmid  !

                ! Define normalization factors
                Nnorm = intgN(zs1,N01,Vlo,Vmid) + intgN(zs2,N02,Vmid,Vhi)
                Vnorm = intgV(zs1,N01,Vlo,Vmid) + intgV(zs2,N02,Vmid,Vhi)

                ! Accumulated variables
                zactvol = 0.
                Nact = 0.
                Vact(:) = 0.

                ! Integration over each size range within a bin
                IF ( intrange(1) ) THEN
                   Nact = Nact + (Nmid/Nnorm)*intgN(zs1,N01,Vlo,vint1)
                   zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs1,N01,Vlo,vint1)
                END IF

                IF ( intrange(2) ) THEN
                   Nact = Nact + (Nmid/Nnorm)*intgN(zs1,N01,vint1,Vmid)
                   zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs1,N01,vint1,Vmid)
                END IF

                IF ( intrange(3) ) THEN
                   Nact = Nact + (Nmid/Nnorm)*intgN(zs2,N02,Vmid,vint2)
                   zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs2,N02,Vmid,vint2)
                END IF

                IF ( intrange(4) ) THEN
                   Nact = Nact + (Nmid/Nnorm)*intgN(zs2,N02,vint2,Vhi)
                   zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs2,N02,vint2,Vhi)
                END IF

                DO ss = 1,8
                   Vact(ss) = zactvol*( paero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
                END DO

                ! Store the number concentration and mass of activated particles for current bins
                zactd(cb)%numc = MIN(Nact,Nmid)
                zactd(cb)%volc(:) = MIN(Vact(:),paero(ii,jj,ab)%volc(:))
                

                IF (zactd(cb)%numc < 0.) THEN
                   WRITE(*,*) Nim1,Nmid,Nip1
                   WRITE(*,*) dNim1,dNmid,dNip1
                   WRITE(*,*) Vim1,Vlo,Vmid,Vhi,Vip1
                   WRITE(*,*) N01,N02,zs1,zs2
                   WRITE(*,*) vint1,vint2
                END IF

             END DO ! cb

             IF (ANY(zactd(:)%numc < 0.)) THEN
                WRITE(*,*) 'NEGATIVE ACTIVATED'
                WRITE(*,*) zactd(:)%numc
                WRITE(*,*) '---------'
                WRITE(*,*) 'Setting MAX(Nact,0)'
                zactd(:)%numc = MAX(zactd(:)%numc, 0.)
                DO ss = 1,8
                   WHERE(zactd(:)%numc == 0.) zactd(:)%volc(ss) = 0.
                END DO

             END IF
             
             ! Apply the number and mass activated to aerosol and cloud bins
             !WRITE(*,*) zactd(1:7)%numc/paero(ii,jj,4:10)%numc
             paero(ii,jj,ica%par:fca%par)%numc =   &
                  MAX(0., paero(ii,jj,ica%par:fca%par)%numc - zactd(ica%cur:fca%cur)%numc)
             pcloud(ii,jj,ica%cur:fca%cur)%numc = pcloud(ii,jj,ica%cur:fca%cur)%numc + zactd(ica%cur:fca%cur)%numc
             DO ss = 1,8
                paero(ii,jj,ica%par:fca%par)%volc(ss) =  &
                     MAX(0., paero(ii,jj,ica%par:fca%par)%volc(ss) - zactd(ica%cur:fca%cur)%volc(ss))
                pcloud(ii,jj,ica%cur:fca%cur)%volc(ss) = pcloud(ii,jj,ica%cur:fca%cur)%volc(ss) + zactd(ica%cur:fca%cur)%volc(ss)
             END DO
             
          END IF ! RH limit

       END DO ! ii

    END DO ! jj

  END SUBROUTINE actInterst


  ! ----------------------------------------------

  SUBROUTINE activate3(kproma,kbdim,klev,paero,pbcrita,pbcritb, &
                       pdcrit, pdcrlo, pdcrhi, pdcstar, pactd   )
    !
    ! Gets the number and mass activated in the critical aerosol size bin
    ! TEKEE NYT VAIN a-binit!!!!
    !
    
    USE mo_submctl, ONLY : t_parallelbin, t_section, pi6, nlim, fn2b, ncld,  &
                               in1a,in2a,fn1a,fn2a, ica,fca,icb,fcb
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kproma,kbdim,klev
    TYPE(t_section), INTENT(IN) :: paero(kbdim,klev,fn2b)
    INTEGER, INTENT(IN) :: pbcrita(kbdim,klev),         & ! Index of the critical aerosol bin in regime a
                           pbcritb(kbdim,klev)            ! Index of the critical aerosol bin in regime b
    REAL, INTENT(IN) :: pdcrit(kbdim,klev,fn2b),    & ! Bin middle critical diameter
                            pdcrlo(kbdim,klev,fn2b),    & ! Critical diameter at low limit
                            pdcrhi(kbdim,klev,fn2b)       ! Critical diameter at high limit
    REAL, INTENT(IN) :: pdcstar(kbdim,klev)           ! Critical diameter corresponding to Smax
    TYPE(t_section), INTENT(OUT) :: pactd(kbdim,klev,ncld) ! Properties of the maximum amount of newly activated droplets


    REAL :: zvcstar, zvcint
    REAL :: zs1,zs2             ! Slopes
    REAL :: Nmid, Nim1,Nip1
    REAL :: dNmid, dNim1, dNip1
    REAL :: Vmid,Vlo,Vhi
    REAL :: Vim1,Vlom1,Vhim1
    REAL :: Vip1,Vlop1,Vhip1
    REAL :: Nnorm,Vnorm,N01,N02
    REAL :: zactvol

    REAL :: zhlp
    REAL :: zratio

    INTEGER :: ii,jj,ss,cb,ab

    DO jj = 1,klev
       DO ii = 1,kproma

          ! This means in practice that vertical velocity is <= 0 or Ntot == 0
          IF ( ALL(pdcrit(ii,jj,:) < epsilon(1.0)) ) CYCLE
          
          zvcstar = 0.

          IF ( paero(ii,jj,pbcrita(ii,jj))%numc < nlim ) THEN
             Vmid = pi6*paero(ii,jj,pbcrita(ii,jj))%dmid**3
          ELSE
             Vmid = SUM( paero(ii,jj,pbcrita(ii,jj))%volc(1:7) )/MAX(nlim,paero(ii,jj,pbcrita(ii,jj))%numc)
          END IF
          Vhi = paero(ii,jj,pbcrita(ii,jj))%vhilim
          Vlo = paero(ii,jj,pbcrita(ii,jj))%vlolim


          IF ( pdcstar(ii,jj) >= pdcrit(ii,jj,pbcrita(ii,jj)) ) THEN
             zhlp = ( pi6*pdcstar(ii,jj)**3 - pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 ) / &
                    MAX( epsilon(1.0), pi6*pdcrhi(ii,jj,pbcrita(ii,jj))**3 - pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 )

             zvcstar = Vmid + zhlp*(Vhi-Vmid)
          ELSE IF (pdcstar(ii,jj) < pdcrit(ii,jj,pbcrita(ii,jj)) ) THEN

             zhlp = ( pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 - pi6*pdcstar(ii,jj)**3 ) / &
                    MAX( epsilon(1.0), pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 - pi6*pdcrlo(ii,jj,pbcrita(ii,jj))**3 )

             zvcstar = Vmid - zhlp*(Vmid-Vlo)
          END IF

          zvcstar = MAX( zvcstar, paero(ii,jj,pbcrita(ii,jj))%vlolim )
          zvcstar = MIN( zvcstar, paero(ii,jj,pbcrita(ii,jj))%vhilim ) 

          ! Loop over cloud droplet (and aerosol) bins
          DO cb = ica%cur,fca%cur
             ab = ica%par + (cb-ica%cur)

             IF ( paero(ii,jj,ab)%numc < nlim) CYCLE

             ! Formulate a slope for Wet particle size within bins and integrate over
             ! the particles larger than zvcstar

             Nmid = MAX(paero(ii,jj,ab)%numc, nlim)
             Vmid = SUM(paero(ii,jj,ab)%volc(1:7))/Nmid ! Dry bin mid volume
             Vlo = paero(ii,jj,ab)%vlolim      ! Mid dry volume scaled to bin low limit (this is mostly an educated guess... )
             Vhi = paero(ii,jj,ab)%vhilim      ! Same for high limit

             IF ( ab > in1a ) THEN
                Nim1 = MAX(paero(ii,jj,ab-1)%numc, nlim)
                IF (Nim1 > nlim) THEN
                   Vim1 = SUM(paero(ii,jj,ab-1)%volc(1:7))/Nim1
                ELSE
                   Vim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                END IF
                Vlom1 = paero(ii,jj,ab-1)%vlolim
                Vhim1 = paero(ii,jj,ab-1)%vhilim
             ELSE ! ab == in1a
                Nim1 = nlim
                Vim1 = Vlo/2.
                Vlom1 = 0.
                Vhim1 = Vlo
             END IF

             IF ( ab < fn2a ) THEN
                Nip1 = MAX(paero(ii,jj,ab+1)%numc, nlim)
                IF (Nip1 > nlim) THEN
                   Vip1 = SUM(paero(ii,jj,ab+1)%volc(1:7))/Nip1
                ELSE
                   Vip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                END IF
                Vlop1 = paero(ii,jj,ab+1)%vlolim
                Vhip1 = paero(ii,jj,ab+1)%vhilim
             ELSE ! ab == fn2a
                Nip1 = nlim
                Vip1 = Vhi + 0.5*(Vhi-Vlo)
                Vlop1 = Vhi
                Vhip1 = Vhi + (Vhi-Vlo)
             END IF

             Vip1 = MAX(Vlop1,MIN(Vip1,Vhip1))
             Vim1 = MAX(Vlom1,MIN(Vim1,Vhim1))
             

             ! get density distribution values for
             dNim1 = Nim1/(Vhim1-Vlom1)
             dNip1 = Nip1/(Vhip1-Vlop1)
             dNmid = Nmid/(Vhi-Vlo)

             ! Get slopes
             zs1 = ( dNmid - dNim1 )/( Vmid - Vim1 )
             zs2 = ( dNip1 - dNmid )/( Vip1 - Vmid )

             N01 = dNmid - zs1*Vmid  ! Origins
             N02 = dNip1 - zs2*Vip1  !

             ! Define normalization factors
             Nnorm = intgN(zs1,N01,Vlo,Vmid) + intgN(zs2,N02,Vmid,Vhi)
             Vnorm = intgV(zs1,N01,Vlo,Vmid) + intgV(zs2,N02,Vmid,Vhi)

             IF (zvcstar < Vmid) THEN

                ! Use actual critical volume only in the critical bin, otherwise current bin limits
                zvcint = MAX(zvcstar, Vlo)

                pactd(ii,jj,cb)%numc = (Nmid/Nnorm) * ( intgN(zs1,N01,zvcint,Vmid) + intgN(zs2,N02,Vmid,Vhi) )
                ! For different species, assume the mass distribution identical in particles within the bin
                zactvol = (Nmid*Vmid/Vnorm) * ( intgV(zs1,N01,zvcint,Vmid) + intgV(zs2,N02,Vmid,Vhi) )
                DO ss = 1,8
                   pactd(ii,jj,cb)%volc(ss) = zactvol*( paero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
                END DO

             ELSE IF (zvcstar >= Vmid) THEN

                ! Use actual critical volume only in the critical bin, otherwise current bin limits
                zvcint = MIN(zvcstar,Vhi)

                pactd(ii,jj,cb)%numc = (Nmid/Nnorm) * ( intgN(zs2,N02,zvcint,Vhi) )
                zactvol = (Nmid*Vmid/Vnorm) * ( intgV(zs2,N02,zvcint,Vhi) )
                DO ss = 1,8
                   pactd(ii,jj,cb)%volc(ss) = zactvol*( paero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
                END DO

             END IF

             !IF ( pactd(ii,jj,cb)%numc < 0. .OR. ANY(pactd(ii,jj,cb)%volc(:) < 0.) ) WRITE(*,*) 'HEPHEPHPE'
             IF ( pactd(ii,jj,cb)%numc < -1.e-10 ) THEN
                WRITE(*,*) 'activate3: negative numc, ',pactd(ii,jj,cb)%numc, ab
                WRITE(*,*) dNim1,dNip1,dNmid
                WRITE(*,*) Nnorm,Vnorm
                WRITE(*,*) Vip1, Vhip1,Vlop1
                WRITE(*,*) Vlo,Vhi,Vmid,zvcint
                WRITE(*,*) zs1,zs2,N01,N02
             END IF

             pactd(ii,jj,cb)%numc = MAX(0., pactd(ii,jj,cb)%numc)
             DO ss = 1,8
                pactd(ii,jj,cb)%volc(ss) = MAX(0., pactd(ii,jj,cb)%volc(ss))
             END DO

             ! "Artificially" adjust the wet size of newly activated a little bit to prevent them from being
             ! evaporated immediately
             pactd(ii,jj,cb)%volc(8) = pactd(ii,jj,cb)%numc*pi6*(pdcrit(ii,jj,ab)**3) *  &
                                       MIN(2.,(3.e-6/max(epsilon(1.0),pdcrit(ii,jj,ab)))**2)

          END DO ! cb

       END DO ! ii
    END DO ! jj

  END SUBROUTINE activate3
  ! ------------------------------------------------
  REAL FUNCTION intgN(ikk,icc,ilow,ihigh)
    ! Gets the integral over a (linear) number concentration distribution
    !
    
    IMPLICIT NONE
    REAL, INTENT(in) :: ikk,icc,ilow,ihigh
    intgN = 0.5*ikk*MAX(ihigh**2 - ilow**2,0.) + icc*MAX(ihigh - ilow,0.)
  END FUNCTION intgN
  ! ------------------------------------------------
  REAL FUNCTION intgV(ikk,icc,ilow,ihigh)
    ! Gets the integral over a volume volume distribution based on a linear
    ! number concentration distribution
    
    IMPLICIT NONE
    REAL, INTENT(in) :: ikk,icc,ilow,ihigh
    intgV = (1./3.)*ikk*MAX(ihigh**3 - ilow**3,0.) + 0.5*icc*MAX(ihigh**2 - ilow**2,0.)
  END FUNCTION intgV




  !-----------------------------------------
  SUBROUTINE autoconv2(kproma,kbdim,klev,   &
                      pcloud,pprecp         )
  !
  ! Uses a more straightforward method for converting cloud droplets to drizzle.
  ! Assume a lognormal cloud droplet distribution for each bin. Sigma_g is an adjustable
  ! parameter and is set to 1.2 by default
  !
    
    USE mo_submctl, ONLY : t_section,   &
                               ncld,        &
                               nprc,        &
                               rhowa,       &
                               pi6,         &
                               nlim
    USE mo_constants, ONLY : rd
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kproma,kbdim,klev
    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld)
    TYPE(t_section), INTENT(inout) :: pprecp(kbdim,klev,nprc)

    REAL :: Vrem, Nrem, Vtot, Ntot
    REAL :: dvg,dg

    REAL, PARAMETER :: zd0 = 50.e-6
    REAL, PARAMETER :: sigmag = 1.2

    INTEGER :: ii,jj,cc,ss

    ! Find the cloud bins where the mean droplet diameter is above 50 um
    ! Do some fitting...
    DO jj = 1,klev
       DO ii = 1,kproma
          DO cc = 1,ncld

             Vrem = 0.
             Nrem = 0.
             Ntot = 0.
             Vtot = 0.
             dvg = 0.
             dg = 0.

             Ntot = pcloud(ii,jj,cc)%numc
             Vtot = SUM(pcloud(ii,jj,cc)%volc(:))

             IF ( Ntot > 0. .AND. Vtot > 0. ) THEN

                ! Volume geometric mean diameter
                dvg = ((Vtot/Ntot/pi6)**(1./3.))*EXP( (3.*LOG(sigmag)**2)/2. )
                dg = dvg*EXP( -3.*LOG(sigmag)**2 )

                !testi = cumlognorm(2.19e-4,sigmag,zd0)
                !WRITE(*,*) 'TESTAAN: ',testi

                Vrem = Vtot*( 1. - cumlognorm(dvg,sigmag,zd0) )
                Nrem = Ntot*( 1. - cumlognorm(dg,sigmag,zd0) )



                IF ( Vrem < 0. ) WRITE(*,*) 'ERROR Vrem < 0', Vrem, Vtot
                IF ( Vrem > Vtot ) WRITE(*,*) 'ERROR Vrem > Vtot', cumlognorm(dvg,sigmag,zd0),dvg,dg,Vtot
                IF ( Nrem < 0. ) WRITE(*,*) 'ERROR Nrem < 0', Nrem
                IF ( Nrem > Ntot ) WRITE(*,*) 'ERROR Nrem > Ntot', cumlognorm(dg,sigmag,zd0),dvg,dg,Ntot

                IF ( Vrem > 0. .AND. Nrem > 0.) THEN

                   ! Put the mass and number to the first precipitation bin and remover from
                   ! cloud droplets
                   DO ss = 1,7
                      pprecp(ii,jj,1)%volc(ss) = pprecp(ii,jj,1)%volc(ss) + pcloud(ii,jj,cc)%volc(ss)*(Nrem/Ntot)
                      pcloud(ii,jj,cc)%volc(ss) = pcloud(ii,jj,cc)%volc(ss)*(1. - (Nrem/Ntot))
                   END DO
                   
                   pprecp(ii,jj,1)%volc(8) = pprecp(ii,jj,1)%volc(8) + pcloud(ii,jj,cc)%volc(8)*(Vrem/Vtot)
                   pcloud(ii,jj,cc)%volc(8) = pcloud(ii,jj,cc)%volc(8)*(1. - (Vrem/Vtot))

                   pprecp(ii,jj,1)%numc = pprecp(ii,jj,1)%numc + Nrem
                   pcloud(ii,jj,cc)%numc = pcloud(ii,jj,cc)%numc - Nrem

                END IF ! Nrem Vrem

             END IF ! Ntot Vtot

          END DO ! cc
       END DO ! ii
    END DO ! jj

  END SUBROUTINE autoconv2


  !***********************************************
  !
  ! heterogenous nucleation according to Morrison et al. 2005 (JAS 62:1665-1677) eq. (25)
  ! as of referenced as [Mor05]
  !
  !***********************************************

  SUBROUTINE ice_het_nucl(kproma,kbdim,klev,   &
                      pcloud,pice,paero,ppres, &
                      ptemp,prv,prs,ptstep ) !'debugkebab'



    
    USE mo_submctl, ONLY : t_section,   &
                               in2b,fn2b,   &
                               ica,fca,     &
                               icb,fcb,     &
                               ncld,        &
                               iia,fia,     &
                               iib,fib,     &
                               nice,        &
                               rhowa,       &
                               rhoic,       &
                               planck,      &
                               pi6,         &
                               pi,          &
                               nlim, prlim,eps, debug
    USE mo_constants, ONLY : rd, alf, avo

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,kproma,klev
    REAL, INTENT(in) :: ptstep
    REAL, INTENT(in) :: ppres(kbdim,klev),  &
                            ptemp(kbdim,klev),  &
                            prv(kbdim,klev),    &
                            prs(kbdim,klev)

    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), &
                                      pice(kbdim,klev,nice),  &
                                      paero(kbdim,klev,fn2b)

    REAL, PARAMETER :: zd00 = 30.e-6,  &  ! Drizzle onset mean diameter
                           zd0 = 50.e-6,   &  ! Initial drizzle droplet diameter
                           ze1 = 2.47,     &
                           ze2 = -1.79,    &
                           mult = 1.e-3
    REAL :: dmixr,   &  ! Change in cloud droplet mixing ratio
                cmixr,   &  ! Current mixing ratio
                dnumb,   &  ! Change in droplet number concentration
                frvol       ! Fractional change in volume concentration





    INTEGER :: ii,jj,kk,ss
    INTEGER :: hh
    REAL :: phf = 0., & ! probability of homogeneous freezing of a wet aerosol particle
                rn, & !radius of the insoluble portion of the aerosol
                rdry,qv,jcf
    LOGICAL :: freez
    REAL :: Vrem, Nrem, Vtot, Ntot, frac


    DO kk = in2b, fn2b ! insoluble materials !1,nice
       DO ii = 1,kproma

          DO jj = 1,klev


              rdry = (3.*sum(paero(ii,jj,kk)%volc) /paero(ii,jj,kk)%numc/4./pi)**(1./3.) !! dry radius of particle  !!huomhuom pcloud vai paero
              qv = (1.-sum( paero(ii,jj,kk)%volc(3:4) ))/(1. - sum(paero(ii,jj,kk)%volc(1:7))) !!! the volume soluble fraction of the aerosol huomhuom tarkista tämän laskeminen
              rn = rdry*(1-qv)**(1./3.)
              jcf = calc_JCF( rn,ptemp(ii,jj), ppres(ii,jj), prv(ii,jj), prs(ii,jj) )
              phf = 1 - exp( -jcf*ptstep )
!              if (phf > 0.)  write(*,*) 'phf ', phf, ' heterogenous debugkebab'
              Ntot = pcloud(ii,jj,kk)%numc
              Vtot = SUM(pcloud(ii,jj,kk)%volc(:))

!              if ( phf> 0. ) write (*,*) 'phf ',phf , ' debugkebab'

              frac = MAX(0.,MIN(1.,phf))

              DO ss = 1,8
                      pice(ii,jj,kk)%volc(ss) = max(0., pice(ii,jj,kk)%volc(ss) + &
                                                           pcloud(ii,jj,kk)%volc(ss)*frac)

                      pcloud(ii,jj,kk)%volc(ss) = max(0., pcloud(ii,jj,kk)%volc(ss)*(1. - frac))
               END DO

               pice(ii,jj,kk)%numc = max( 0., pice(ii,jj,kk)%numc + pcloud(ii,jj,kk)%numc*frac )
               pcloud(ii,jj,kk)%numc = max(0., pcloud(ii,jj,kk)%numc*(1. - frac) )

          END DO
       END DO
    END DO
IF (debug)             write(*,*)  'nyt on jäänukleoitu ', ' debugkebab'

  END SUBROUTINE ice_het_nucl
  !***********************************************
  !
  ! homogenous nucleation according to Morrison et al. 2005 (JAS 62:1665-1677) eq. (27)
  ! as of referenced as [Mor05]
  !
  !***********************************************
  SUBROUTINE ice_hom_nucl(kproma,kbdim,klev,   &
                      pcloud,pice,paero,ppres, &
                      ptemp,prv,prs,ptstep ) 

    
    USE mo_submctl, ONLY : t_section,   &
                               in2b,fn2b,   &
                               ica,fca,     &
                               icb,fcb,     &
                               ncld,        &
                               iia,fia,     &
                               iib,fib,     &
                               nice,        &
                               rhowa,       &
                               rhoic,       &
                               planck,      &
                               pi6,         &
                               pi,          &
                               nlim, prlim,eps
    USE mo_constants, ONLY : rd, alf, avo

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,kproma,klev
    REAL, INTENT(in) :: ptstep
    REAL, INTENT(in) :: ppres(kbdim,klev),  &
                            ptemp(kbdim,klev),  &
                            prv(kbdim,klev),    &
                            prs(kbdim,klev)

    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), &
                                      pice(kbdim,klev,nice),  &
                                      paero(kbdim,klev,fn2b)

    REAL, PARAMETER :: zd00 = 30.e-6,  &  ! Drizzle onset mean diameter
                           zd0 = 50.e-6,   &  ! Initial drizzle droplet diameter
                           ze1 = 2.47,     &
                           ze2 = -1.79,    &
                           mult = 1.e-3
    REAL :: dmixr,   &  ! Change in cloud droplet mixing ratio
                cmixr,   &  ! Current mixing ratio
                dnumb,   &  ! Change in droplet number concentration
                frvol       ! Fractional change in volume concentration


    INTEGER :: ii,jj,kk,ss
    INTEGER :: hh
    REAL :: phf = 0., & ! probability of homogeneous freezing of a wet aerosol particle
                rn, & !radius of the insoluble portion of the aerosol
                rdry, qv,rw,NL,jhf
    LOGICAL :: freez
    REAL :: Vrem, Nrem, Vtot, Ntot,frac


    DO kk = 1,nice
       DO ii = 1,kproma

          DO jj = 1,klev

              rdry = (3.*sum(paero(ii,jj,kk)%volc(1:7)) /paero(ii,jj,kk)%numc/4./pi)**(1./3.)
              qv = (1.-sum(paero(ii,jj,kk)%volc(3:4)))/(1. - sum(paero(ii,jj,kk)%volc(1:7))) !!! the volume soluble fraction of the aerosol huomhuom tarkista tämän laskeminen #arvo
              rn = rdry*(1-qv)**(1./3.) !! radius of the insoluble portion of particle
              rw = paero(ii,jj,kk)%dwet  ! #arvo
              NL = paero(ii,jj,kk)%numc !! #arvo
              if (ptemp(ii,jj) > 243 ) cycle
              jhf = calc_JHF( NL , ptemp(ii,jj))
              phf = 1 - exp( -jhf*pi6*( rw**3 - rn**3 )*ptstep)
!              if (phf > 0.)  write(*,*) 'phf ', phf, ' homogenous debugkebab'
              Ntot = paero(ii,jj,kk)%numc
              Vtot = SUM(paero(ii,jj,kk)%volc(:))

              frac = MAX(0.,MIN(1.,phf))

              DO ss = 1,8
                      pice(ii,jj,kk)%volc(ss) = max(0., pice(ii,jj,kk)%volc(ss) + &
                                                           pcloud(ii,jj,kk)%volc(ss)*frac)

                      pcloud(ii,jj,kk)%volc(ss) = max(0., pcloud(ii,jj,kk)%volc(ss)*(1. - frac))
               END DO

               pice(ii,jj,kk)%numc = max( 0., pice(ii,jj,kk)%numc + pcloud(ii,jj,kk)%numc*frac )
               pcloud(ii,jj,kk)%numc = max(0., pcloud(ii,jj,kk)%numc*(1. - frac) )


              DO ss = 1,8
                      pice(ii,jj,kk)%volc(ss) = max(0., pice(ii,jj,kk)%volc(ss) + &
                                                           paero(ii,jj,kk)%volc(ss)*frac)

                      paero(ii,jj,kk)%volc(ss) = max(0., paero(ii,jj,kk)%volc(ss)*(1. - frac))
               END DO

               pice(ii,jj,kk)%numc = max( 0., pice(ii,jj,kk)%numc + paero(ii,jj,kk)%numc*frac )
               paero(ii,jj,kk)%numc = max(0., paero(ii,jj,kk)%numc*(1. - frac) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          END DO
       END DO
    END DO


  END SUBROUTINE ice_hom_nucl

  !***********************************************
  !
  ! heterogenous immersion nucleation according to Dieh & Wurzler 2006 JAS 61:2063-2072
  ! as of referenced as [Mor05]
  !
  !***********************************************
  SUBROUTINE ice_immers_nucl(kproma,kbdim,klev,   &
                      pcloud,pice,ppres, &
                      ptemp,ptt,prv,prs,ptstep, time )



    
    USE mo_submctl, ONLY : t_section,   &
                               ica,fca,     &
                               icb,fcb,     &
                               ncld,        &
                               iia,fia,     &
                               iib,fib,     &
                               nice,        &
                               rhowa,       &
                               rhoic,       &
                               planck,      &
                               pi6,         &
                               pi,          &
                               nlim, prlim,eps, &
                               debug
    USE mo_constants, ONLY : rd, alf, avo

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,kproma,klev
    REAL, INTENT(in) :: ptstep, time
    REAL, INTENT(in) :: ppres(kbdim,klev),  &
                            ptemp(kbdim,klev),  &
                            ptt(kbdim,klev),    &
                            prv(kbdim,klev),    &
                            prs(kbdim,klev)

    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), &
                                      pice(kbdim,klev,nice)

    REAL, PARAMETER :: zd00 = 30.e-6,  &  ! Drizzle onset mean diameter
                           zd0 = 50.e-6,   &  ! Initial drizzle droplet diameter
                           ze1 = 2.47,     &
                           ze2 = -1.79,    &
                           mult = 1.e-3
    REAL :: dmixr,   &  ! Change in cloud droplet mixing ratio
                cmixr,   &  ! Current mixing ratio
                dnumb,   &  ! Change in droplet number concentration
                frvol       ! Fractional change in volume concentration


    INTEGER :: ii,jj,kk,ss
    INTEGER :: hh
    REAL :: frac_numc = 0., & ! probability of homogeneous freezing of a wet aerosol particle
                frac_volc = 0., &
                frac,              &
                Ts
    LOGICAL :: freez
    REAL :: Vrem, Nrem, Vtot, Ntot, &
                a_kiehl, B_kiehl, nucl_rate, &
                Nicetot, Vicetot, Vinsolub,Temp_tend

      B_kiehl=1.0e-6
      a_kiehl=1.0
open(23, file="immersoitu", position='append')
    DO kk =1,nice  ! insoluble materials !
       DO ii = 1,kproma

          DO jj = 1,klev

              Ntot = pcloud(ii,jj,kk)%numc
              Vtot = SUM(pcloud(ii,jj,kk)%volc(:))
              write(23,'(A12,A3,I3,A5,F5.1, 3(A5,ES9.1E3))') &
              'tulostusta','kk', kk, 'aika',time, 'Vtot',Vtot,'Ntot',Ntot, 'dwet', pcloud(ii,jj,kk)%dwet
              Vinsolub = SUM(pcloud(ii,jj,kk)%volc(3:4))

              Nicetot = pice(ii,jj,kk)%numc
              Vicetot = SUM(pice(ii,jj,kk)%volc(:))

              Ts = 273.15-ptemp(ii,jj)
              Temp_tend = ptt(ii,jj)
              if (Temp_tend > 0. .and. abs(Ntot)<eps .and.  Vinsolub < eps)  then
                    write(23,*) 'temp_0_N_0_insolub_0', time
              elseif(Temp_tend > 0. .and. abs(Ntot)<eps .and.  Vinsolub > eps) then
                    write(23,*) 'temp_0_N_0_insolub_1', time
              elseif(Temp_tend > 0. .and. abs(Ntot)>eps .and.  Vinsolub < eps) then
                    write(23,*) 'temp_0_N_1_insolub_0', time
              elseif(Temp_tend > 0. .and. abs(Ntot)>eps .and.  Vinsolub > eps) then
                    write(23,*) 'temp_0_N_1_insolub_1', time, Temp_tend
              elseif(Temp_tend < 0. .and. abs(Ntot)<eps .and.  Vinsolub < eps) then
                    write(23,*) 'temp_1_N_0_insolub_0', time, Temp_tend
              elseif(Temp_tend < 0. .and. abs(Ntot)<eps .and.  Vinsolub > eps) then
                    write(23,*) 'temp_1_N_0_insolub_1', time, Temp_tend
              elseif(Temp_tend < 0. .and. abs(Ntot)>eps .and.  Vinsolub < eps) then
                    write(23,*) 'temp_1_N_1_insolub_0', time, Temp_tend
              elseif(Temp_tend < 0. .and. abs(Ntot)>eps .and.  Vinsolub > eps) then
                    write(23,'(A25,A5,F7.1, A5, I3, 9(A15,ES10.1E3))') 'temp_1_N_1_insolub_1', &
                    'aika',time,'bini', kk,&
                    'tend', Temp_tend,              'Ntot', Ntot,           'Nicetot', Nicetot,&
                    'Vtot', Vtot,                   'Vinsolub', Vinsolub,   'Vtot/Ntot', Vtot/Ntot, &
                    'Vinsolub/Vtot', Vinsolub/Vtot, 'dwet cloud', pcloud(ii,jj,kk)%dwet,&
                    'dwet ice', pice(ii,jj,kk)%dwet
              end if


              if (Temp_tend > 0. .or. abs(Ntot)<eps .or.  Vinsolub < eps) cycle

              if (Temp_tend > 0. .or. abs(Ntot)<eps .or.  Vinsolub < eps) cycle
!              write(*,*) 'nyt vois tulla jääpilvee debugkebab'
              nucl_rate = Vtot*(-a_kiehl)*B_kiehl*exp(a_kiehl*Ts)*Temp_tend*ptstep

              frac = (nucl_rate+Nicetot)/max(eps,Ntot)

              frac = MAX(0.,MIN(1.,frac))

              DO ss = 1,8
                      pice(ii,jj,kk)%volc(ss) = max(0., pice(ii,jj,kk)%volc(ss) + &
                                                           pcloud(ii,jj,kk)%volc(ss)*frac)

                      pcloud(ii,jj,kk)%volc(ss) = max(0., pcloud(ii,jj,kk)%volc(ss)*(1. - frac))
               END DO

               pice(ii,jj,kk)%numc = max( 0., pice(ii,jj,kk)%numc + pcloud(ii,jj,kk)%numc*frac )
               pcloud(ii,jj,kk)%numc = max(0., pcloud(ii,jj,kk)%numc*(1. - frac) )

               if(pice(ii,jj,kk)%numc-Nicetot>0.0) write(23,*) 'jäätä', pice(ii,jj,kk)%numc-Nicetot, frac, Ts, Temp_tend






          END DO
       END DO
    END DO
close(23)
    IF (debug)              write(*,*)  'nyt on immersoitu ', ' debugkebab'

  END SUBROUTINE ice_immers_nucl

  ! ------------------------------------------------------------

  REAL FUNCTION calc_JCF(rn,temp,ppres,prv,prs) ! heterogenous (condensation) freezing  !!check  [Mor05] eq. (26)
                      !the rate of germ formation per volume of solution
        
        USE mo_submctl, ONLY : boltz, planck,pi
        REAL, INTENT(in) :: rn,  &
                              temp,ppres, prv,prs
        REAL :: c_1s, psi
        psi = 1.
        c_1s = 1.e19 !! 10**15 cm^-2 !! concentration of water molecules adsorbed on 1 cm^-2 of surface
        calc_JCF= boltz*temp/planck*psi*c_1s*4.*pi*rn**2.*&
                  exp((-calc_act_energy(temp,'het')-calc_crit_energy(rn,prv,prs,temp))/(boltz*temp))

  END FUNCTION calc_JCF

  ! ------------------------------------------------------------

  REAL FUNCTION calc_JHF(NL,temp) ! homogenous freezing !! Khovosrotyanov & Sassen 1998 [KS98] eq. (7)

    
    USE mo_submctl, ONLY : boltz, planck,surfi0,pi
    REAL, intent(in) :: NL, & !  number of water molecules per unit volume of the liquid
                            temp
    REAL :: r_cr

    r_cr = calc_r_cr(temp)

    calc_JHF = NL*boltz*temp/planck*exp((-4.*pi/3.*surfi0*r_cr**2-calc_act_energy(temp,'hom'))/(boltz*temp))

  END FUNCTION calc_JHF

  ! ------------------------------------------------------------

  REAL FUNCTION calc_act_energy(temp,nucltype) ! activation energy of solution ice interface  !!check

    
    REAL, INTENT(in) :: temp
    CHARACTER(len=*), INTENT(in) :: nucltype
    REAL :: Tc
    Tc = temp-273.15

    calc_act_energy = 0.

    select case(nucltype)

        case('hom')
            ![KC00] p. 4084 beginning of chapter 3.2.
            calc_act_energy = 0.694e-12 * (1.000+ 0.027*(Tc+30.000)*exp(0.010*(Tc+30.000)))
        case('het')
            ! Pruppacher & Klett 1997 [PK97] eq. (3-22)
            a0 = 5.550
            a1 = -8.423e-3
            a2 = 6.384e-4
            a3 = 7.891e-6
            calc_act_energy = 4178.800*a0*exp(a1*Tc + a2*Tc**2 + a3*Tc**3)
     end select
  END FUNCTION calc_act_energy

  ! ------------------------------------------------------------

  REAL FUNCTION calc_crit_energy(rn,prv,prs,temp) ! critical energy KC[00] (eq. 2.10) !!huomhuom #arvo
    
    USE mo_submctl, ONLY : surfi0
    REAL, INTENT(in) :: rn, prv,prs, temp
    REAL :: mis, r_g, x, sigma_is,sigma_ns,sigma_ni

    sigma_is = surfi0!! surface tension of ice-solution interface #arvo
    sigma_ns = 1 ! #arvo
    sigma_ni = 0.5 ! #arvo
    alpha = 0 ! #arvo

    mis = (sigma_ns - sigma_ni )/sigma_is
    r_g = calc_r_g(sigma_is,prv,prs,temp)
    x = rn/r_g
    calc_crit_energy = 4.*pi/3.*sigma_is*r_g**2*calc_shapefactor(mis,x)-alpha*(1-mis)*rn**2

  END FUNCTION calc_crit_energy

  ! ------------------------------------------------------------

  ! [KC00] eq. (2.9)
  REAL FUNCTION calc_shapefactor(m,x) !! according to Khvorostyanov & Curry, Geophysical Research letters 27(24):4081-4084, 
                                      !december 2000  !!check
                                 !! as of referenced as [KC00]
    
    REAL, INTENT(IN) :: m,x
    REAL :: psi,fii
    fii = (1.-2.*m*x+x**2)**(0.5)
    psi = (x-m)/fii

    calc_shapefactor = 1. + ( ( 1.-m*x)/fii)**3 + (2.-3.*psi-psi**3)*x**3 + 3.*m*(psi-1.)*x**2

    calc_shapefactor = 0.5*calc_shapefactor

  END FUNCTION calc_shapefactor

  ! ------------------------------------------------------------

  ! [KC00] eq. (2.6)
  REAL FUNCTION calc_r_g(sigma_is,prv,prs,temp) !! calculate ice germ radius [KC00] !!huomhuom !! sigma_is #arvo
    
    USE mo_submctl, ONLY : rhoic,rg,mwa
    REAL, intent(in) :: sigma_is, prv, prs, temp
    REAL :: Late, epsi,temp00,GG,C

    Late = calc_Lefm(temp)
    GG = rg*temp*Late/mwa
    C = 1.7e10 !! 1.7*10^10 Pa == 1.7*10^11 dyn cm^-2
    epsi = 0.025 ! 2.5%
    temp00 = calc_temp00(temp)
    calc_r_g = 2.*sigma_is/( rhoic*Late*log((temp00/temp)*(prv/prs)**GG) -C*epsi**2) !! huomhuom täydennä

  END FUNCTION calc_r_g

  ! ------------------------------------------------------------

  ! [KS98] eq. (8)
  REAL FUNCTION calc_r_cr(temp) !! calculate ice embryo radius [KC00] !!check
    
    USE mo_submctl, ONLY : rhoic,surfi0
    REAL, intent(in) :: temp
    REAL :: Late, temp00

    Late = calc_Lefm(temp)


    temp00 = calc_temp00(temp)
    calc_r_cr = 2*surfi0/( rhoic*Late*log(temp00/temp)) !! huomhuom täydennä

  END FUNCTION calc_r_cr

  ! ------------------------------------------------------------

  ! Harri Kokkola pilvikurssi eq. (2.43)
  REAL FUNCTION calc_Lefm(temp) !! Latent heat of fusion !!check
    
    REAL, intent(in) :: temp
    REAL :: Tc ! temperature in celsius degrees
    Tc = temp-273.15
    calc_Lefm = 2.83458e+6-Tc*(340.+10.46*Tc)

  END FUNCTION calc_Lefm

  ! ------------------------------------------------------------

  REAL function calc_temp00(temp) !!freezing point depression !! huomhuom korjaa parametrit ja täydennä #arvo
    
    REAL, intent(in) :: temp

    calc_temp00 = 273.15

  END FUNCTION calc_temp00

  ! ------------------------------------------------------------

  SUBROUTINE ice_melt(kproma,kbdim,klev,   &
                      pcloud,pice,pprecp,psnow,ppres, &
                      ptemp,prv,prs,ptstep )



    
    USE mo_submctl, ONLY : t_section,   &
                               ica,fca,     &
                               icb,fcb,     &
                               ncld,        &
                               iia,fia,     &
                               iib,fib,     &
                               nice,        &
                               nsnw,        &
                               nprc,        &
                               rhowa,       &
                               pi6,         &
                               nlim, prlim
    USE mo_constants, ONLY : rd

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,kproma,klev
    REAL, INTENT(in) :: ptstep
    REAL, INTENT(in) :: ppres(kbdim,klev),  &
                            ptemp(kbdim,klev),  &
                            prv(kbdim,klev),    &
                            prs(kbdim,klev)

    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), &
                                      pice(kbdim,klev,nice),   &
                                      psnow(kbdim,klev,nsnw),  &
                                      pprecp(kbdim,klev,nprc)


    INTEGER :: ii,jj,kk,ss
    INTEGER :: hh
    REAL :: zrh

    DO ii = 1,kproma
       DO jj = 1,klev
          if (ptemp(ii,jj) < 273.15 ) cycle !!huomhuom add effect of freezing point depression
          DO kk = 1,nice

             DO ss = 1,8
                pcloud(ii,jj,kk)%volc(ss) = pice(ii,jj,kk)%volc(ss) + pcloud(ii,jj,kk)%volc(ss)
                pice(ii,jj,kk)%volc(ss) = 0.

             END DO

                DO ss = 1,8
                        pcloud(ii,jj,kk)%volc(ss) = pice(ii,jj,kk)%volc(ss) + pcloud(ii,jj,kk)%volc(ss)
                        pice(ii,jj,kk)%volc(ss) = 0.

                END DO

               pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc + pice(ii,jj,kk)%numc
               pice(ii,jj,kk)%numc = 0.


            END DO

            DO kk =1,nsnw
                DO ss = 1,8
                        pprecp(ii,jj,kk)%volc(ss) = psnow(ii,jj,kk)%volc(ss) + pprecp(ii,jj,kk)%volc(ss)
                        psnow(ii,jj,kk)%volc(ss) = 0.

                END DO

               pprecp(ii,jj,kk)%numc = pprecp(ii,jj,kk)%numc + psnow(ii,jj,kk)%numc
               psnow(ii,jj,kk)%numc = 0.
            END DO
       END DO
    END DO

  END SUBROUTINE ice_melt


  SUBROUTINE autosnow(kproma,kbdim,klev,   &
                      pice,psnow         )
  !
  ! Uses a more straightforward method for converting cloud droplets to drizzle.
  ! Assume a lognormal cloud droplet distribution for each bin. Sigma_g is an adjustable
  ! parameter and is set to 1.2 by default
  !
    
    USE mo_submctl, ONLY : t_section,   &
                               iia,fia,     &
                               iib,fib,     &
                               isa,fsa,     &
                               nice,        &
                               nsnw,        &
                               rhowa,       &
                               pi6,         &
                               nlim, debug
    USE mo_constants, ONLY : rd
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kproma,kbdim,klev
    TYPE(t_section), INTENT(inout) :: pice(kbdim,klev,nice)
    TYPE(t_section), INTENT(inout) :: psnow(kbdim,klev,nsnw)

    REAL :: Vrem, Nrem, Vtot, Ntot
    REAL :: dvg,dg, frac

    REAL, PARAMETER :: zd0 = 250.e-6  !! huomhuom tarvii varmaan tuunata
    REAL, PARAMETER :: sigmag = 1.2   !! huomhuom tarvii varmaan tuunata

    INTEGER :: ii,jj,cc,ss

    ! Find the ice particle bins where the mean droplet diameter is above 250 um
    ! Do some fitting...
    DO jj = 1,klev
       DO ii = 1,kproma
          DO cc = 1,nice

             Vrem = 0.
             Nrem = 0.
             Ntot = 0.
             Vtot = 0.
             dvg = 0.
             dg = 0.

             Ntot = pice(ii,jj,cc)%numc
             Vtot = SUM(pice(ii,jj,cc)%volc(:))

             IF ( Ntot > 0. .AND. Vtot > 0. ) THEN
!               write(*,*) 'autosnow Ntot, Vtot > 0 debugkebab'
                ! Volume geometric mean diameter
                dvg = pice(ii,jj,cc)%dwet*EXP( (3.*LOG(sigmag)**2)/2. )
                dg = dvg*EXP( -3.*LOG(sigmag)**2 )

                !testi = cumlognorm(2.19e-4,sigmag,zd0)
                !WRITE(*,*) 'TESTAAN: ',testi

                Vrem = Max(0., Vtot*( 1. - cumlognorm(dvg,sigmag,zd0) ) )
                Nrem = Max(0., Ntot*( 1. - cumlognorm(dg,sigmag,zd0) )  )

                IF ( Vrem < 0. ) WRITE(*,*) 'ERROR Vrem < 0', Vrem, Vtot
                IF ( Vrem > Vtot ) WRITE(*,*) 'ERROR Vrem > Vtot', cumlognorm(dvg,sigmag,zd0),dvg,dg,Vtot
                IF ( Nrem < 0. ) WRITE(*,*) 'ERROR Nrem < 0', Nrem
                IF ( Nrem > Ntot ) WRITE(*,*) 'ERROR Nrem > Ntot', cumlognorm(dg,sigmag,zd0),dvg,dg,Ntot

                IF ( Vrem > 0. .AND. Nrem > 0.) THEN
!                   write(*,*) 'autosnow Vrem, Nrem > 0 debugkebab'
                   ! Put the mass and number to the first precipitation bin and remover from
                   ! cloud droplets

                   frac = MAX(0.,MIN(1.,min(Vrem/Vtot, Nrem/Ntot))) ! huomhuom


                   DO ss = 1,8
                      psnow(ii,jj,cc)%volc(ss) = max(0., psnow(ii,jj,cc)%volc(ss) + &
                                                           pice(ii,jj,cc)%volc(ss)*frac)
                      pice(ii,jj,cc)%volc(ss) = max(0., pice(ii,jj,cc)%volc(ss)*(1. - frac))
                    END DO

                    psnow(ii,jj,cc)%numc = max( 0., psnow(ii,jj,cc)%numc + pice(ii,jj,cc)%numc*frac )
                    pice(ii,jj,cc)%numc = max(0., pice(ii,jj,cc)%numc*(1. - frac) )


                END IF ! Nrem Vrem

             END IF ! Ntot Vtot

          END DO ! cc
       END DO ! ii
    END DO ! jj
    IF (debug)              write(*,*)  'nyt on autosnowattu ', ' debugkebab'

  END SUBROUTINE autosnow

  !
  ! -----------------------------------------------------------------
  !
  REAL FUNCTION cumlognorm(dg,sigmag,dpart)
    
    USE mo_submctl, ONLY : pi
    IMPLICIT NONE
    ! Cumulative lognormal function
    REAL, INTENT(in) :: dg
    REAL, INTENT(in) :: sigmag
    REAL, INTENT(in) :: dpart

    REAL :: hlp1,hlp2

    !cumlognorm = 0.

    hlp1 = ( LOG(dpart) - LOG(dg) )
    hlp2 = SQRT(2.)*LOG(sigmag)
    cumlognorm = 0.5 + 0.5*ERF( hlp1/hlp2 )

  END FUNCTION cumlognorm
  !
  ! ----------------------------------------------------------------
  !
  REAL FUNCTION errf(x)
    
    IMPLICIT NONE
    ! (Approximative) Error function.
    ! This is available as an intrinsic function as well but the implementation is somewhat compiler-specific
    ! so an approximation is given here explicitly...

    ! EI TOIMI OIKEIN ÄÄRIPÄISSÄ


    REAL, INTENT(in) :: x
    REAL :: hlp

    hlp = 1. + 0.0705230784*x + 0.0422820123*(x**2) +  &
          0.0092705272*(x**3) + 0.0001520143*(x**4) +       &
          0.0002765672*(x**5) + 0.0000430638*(x**6)

    errf = 1. - 1./(hlp**16)

    WRITE(*,*) 'TESTAAN 2: ', x, errf

  END FUNCTION errf




END MODULE mo_salsa_cloud
