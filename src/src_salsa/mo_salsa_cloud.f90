MODULE mo_salsa_cloud
    USE mo_salsa_types, ONLY : aero, cloud, precp
    USE mo_submctl, ONLY : nbins, ncld, spec,          &
                           in1a,in2a,in2b,fn2a,fn2b,   &
                           ica,icb,fca,fcb,            &
                           pi, pi6, grav, rg,          &
                           surfw0, cpa, mair
    USE util, ONLY : cumlognorm
  IMPLICIT NONE

  !*********************************************************
  !  MO_SALSA_CLOUD
  !*********************************************************
  !
  ! Purpose: Calculates the number of activated cloud droplets
  !
  !
  ! Interface:
  ! ----------
  ! Called from the main aerosol model
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
                              rs,     w,     pactd          )

    USE classSection, ONLY : Section
    USE mo_submctl, ONLY : lsactiv
    
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
    
    ! Properties of newly activate particles
    TYPE(Section), INTENT(out) :: pactd(kbdim,klev,ncld)
    
    INTEGER :: ii, jj, kk
    INTEGER :: nspec
    
    nspec = spec%getNSpec()
    
    ! This is needed for cloud base activation, but must be set to zero for interstitial activation
    DO jj = 1, klev    ! vertical grid
       DO ii = 1, kbdim ! horizontal grid
          ! Reset activated
          DO kk = 1, ncld
             pactd(ii,jj,kk)%volc(1:nspec) = 0.
             pactd(ii,jj,kk)%numc = 0.
          END DO
       END DO
    END DO
    
    ! -------------------------------------
    ! Interstitial activation
    ! -------------------------------------
    IF ( lsactiv%mode == 1 ) THEN
       
       CALL actInterst(kproma,kbdim,klev,rv,rs,temp)
       
    END IF
    
    ! -----------------------------------
    ! Activation at cloud base
    ! -----------------------------------
    IF ( lsactiv%mode == 2 ) THEN
       
       CALL ActCloudBase(kproma,kbdim,klev,pres,temp,w,pactd)
       
    END IF
    
  END SUBROUTINE cloud_activation
  
  
  ! -----------------------------------------------------------------
  ! Calculates the number of moles of dissolved solutes in one particle
  !
  SUBROUTINE getSolute(kproma,kbdim,klev,pns)

    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: kproma,kbdim,klev
    REAL, INTENT(OUT) :: pns(kbdim,klev,fn2b)
    
    REAL :: diss,rho,mm
    INTEGER :: ii,jj,kk,nn
    
    pns = 0.
    
    DO jj = 1, klev
       
       DO ii = 1, kbdim
          
          !-- subranges 1a, 2a and 2b
          
          DO kk = 1, fn2b
             IF (aero(ii,jj,kk)%numc > aero(ii,jj,kk)%nlim) THEN
                
                !-- number of moles of solute in one particle [mol]
                !   BC and dust are insoluble - ignored
                !   SS or NaCl produces 2 ions
                !   SO or H2SO4 produces 3 ions
                !  --- This information is contained in the class object "spec"
                !      so we can loop. The dissolution coefficient is 0 for non-soluble compounds
                !      so they'll get ignored.
                DO nn = 1,spec%getNSpec(type="dry")
                   ! Get parameters
                   diss = spec%diss(nn)
                   mm = spec%MM(nn)
                   rho = spec%rholiq(nn)
                   pns(ii,jj,kk) = pns(ii,jj,kk) + diss*aero(ii,jj,kk)%volc(nn)*rho/mm
                END DO
                pns(ii,jj,kk) = pns(ii,jj,kk)/aero(ii,jj,kk)%numc
                
             END IF
          END DO
          
       END DO
       
    END DO
    
  END SUBROUTINE getSolute
  
  ! -----------------------------------------------------------------

  SUBROUTINE ActCloudBase(kproma,kbdim,klev,pres,temp,w,pactd)
    ! Cloud base activation following:
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
    ! Note: updates pactd, but does not change pcloud?
    ! Note: insoluble species are not properly accounted for
    !
    USE classSection, ONLY : Section

    
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
    
    ! Properties of newly activate particles
    TYPE(Section), INTENT(out) :: pactd(kbdim,klev,ncld)
    
    !-- local variables --------------
    INTEGER :: ii, jj, kk             ! loop indices
    
    INTEGER :: bcrita(kbdim,klev),bcritb(kbdim,klev) ! Index of the critical aerosol bin for regimes a and b
    
    REAL ::                         &
         sil,                       & ! critical supersaturation
         !     at the upper bound of the bin
         siu,                       & !  "  at the lower bound of the bin
         scrit(fn2b),               & !  "  at the center of the bin
         aa,                        & ! curvature (Kelvin) effect [m]
         bb,                        & ! solute (Raoult) effect [m3]
         ns(kbdim,klev,fn2b),                  & ! number of moles of solute
         nshi,                      & !             " at the upper bound of the bin
         nslo,                      & !             " at the lower bound of the bin
         s_max,                     & ! maximum supersaturation
         s_eff,                     & ! effective supersaturation
         x, x1, a1, sum1,       & ! technical variables
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
         zdcrit(kbdim,klev,fn2b),   & ! critical diameter [m]
         zdcstar(kbdim,klev),                   & ! Critical diameter corresponding to Smax
         zdcrhi(kbdim,klev,fn2b),   & ! Critical diameter at the high end of a bin
         zdcrlo(kbdim,klev,fn2b),   & ! Critical diameter at the low end of a bin
         V,                         & ! updraft velocity [m/s]
         rref                       ! reference radius [m]

    ! ------------------------------------------------------------------
    ! Initialization
    ! ------------------------------------------------------------------
    
    bb = 6.*spec%mwa/(pi*spec%rhowa)   ! Raoult effect [m3]
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
    
    ! Get moles of solute at the middle of the bin
    CALL getSolute(kproma,kbdim,klev,ns)
    
    ! ----------------------------------------------------------------
    
    DO jj = 1, klev    ! vertical grid
       DO ii = 1, kbdim ! horizontal grid
          ! Reset activated - done already
          !   pactd(ii,jj,kk)%volc(:) = 0.
          !   pactd(ii,jj,kk)%numc = 0.
          
          ! Positive updraft velocity required
          IF (w(ii,jj) <= 0.) CYCLE
          
          aa = 4.*spec%mwa*surfw0/(rg*spec%rhowa*temp(ii,jj)) ! Kelvin effect [m]
          
          x  = 4.*aa**3/(27.*bb)
          
          ! Get the critical supersaturation for aerosol bins for necessary summation terms (sum1 & ntot)
          
          scrit(in1a:fn2b) = exp(SQRT(x/max(epsilon(1.0),ns(ii,jj,in1a:fn2b)))) - 1.
          
          !-- sums in equation (8), part 3
          ntot = SUM(aero(ii,jj,in1a:fn2b)%numc)
          sum1 = SUM(aero(ii,jj,in1a:fn2b)%numc/scrit(in1a:fn2b)**(2./3.))
          
          IF(ntot < aero(ii,jj,1)%nlim) CYCLE
          V = w(ii,jj)
          
          !-- latent heat of evaporation [J/kg]
          L = 2.501e6-2370.*(temp(ii,jj)-273.15)
          
          !-- saturation vapor pressure of water [Pa]
          a1 = 1.-(373.15/temp(ii,jj))
          ps = 101325.*                                                 &
               exp(13.3185*a1-1.976*a1**2-0.6445*a1**3-0.1299*a1**4)
          
          !-- part 1, eq (11)
          alpha = grav*spec%mwa*L/(cpa*rg*temp(ii,jj)**2)-                            &
               grav*mair/(rg*temp(ii,jj))
          
          !-- part 1, eq (12)
          gamma = rg*temp(ii,jj)/(ps*spec%mwa) &
               + spec%mwa*L**2/(cpa*pres(ii,jj)*mair*temp(ii,jj))
          
          !-- diffusivity [m2/s], Seinfeld and Pandis (15.65)
          x1  = pres(ii,jj) / 101325.
          dv1 = 1.e-4 * (0.211/x1) * ((temp(ii,jj)/273.)**1.94)
          
          rref = 10.e-9
          
          !-- thermal conductivity [J/(m s K)], Seinfeld and Pandis (15.75)
          ka1 = 1.e-3 * (4.39 + 0.071 * temp(ii,jj))
          
          !-- growth coefficient, part 1, eq (16)
          !-- (note: here uncorrected diffusivities and conductivities are used
          !    based on personal communication with H. Abdul-Razzak, 2007)
          Gc = 1./(spec%rhowa*rg*temp(ii,jj)/(ps*dv1*spec%mwa) +                      &
               L*spec%rhowa/(ka1*temp(ii,jj)) * (L*spec%mwa/(temp(ii,jj)*rg)-1.))
          
          !-- effective critical supersaturation: part 3, eq (8)
          s_eff = SQRT( (ntot/sum1) **3)
          
          !-- part 3, equation (5)
          
          theta = ( SQRT( (alpha*V/Gc)**3 ) ) /(2.*pi*spec%rhowa*gamma*ntot)
          
          !-- part 3, equation (6)
          khi = (2./3.)*aa*SQRT(alpha*V/Gc)
          
          !-- maximum supersaturation of the air parcel: part 3, equation (9)
          s_max = s_eff / SQRT(  0.5* SQRT( (khi/theta)**3 ) &
               +  ((s_eff**2)/(theta+3.*khi))**(0.75)  )
          
          !-- Juha: Get the critical diameter corresponding to the maximum supersaturation
          zdcstar = 2.*aa/(3.*s_max)
          
          DO kk = in1a, fn2b
             
             IF (aero(ii,jj,kk)%numc < aero(ii,jj,kk)%nlim) CYCLE
             
             !-- moles of solute in particle at the upper bound of the bin
             nshi = ns(ii,jj,kk)*aero(ii,jj,kk)%vratiohi
             
             !-- critical supersaturation
             sil = exp(SQRT(x/nshi)) - 1.
             
             IF(s_max < sil) CYCLE
             
             !-- moles of solute at the lower bound of the bin:
             nslo = ns(ii,jj,kk)*aero(ii,jj,kk)%vratiolo
             
             !-- critical supersaturation
             siu = exp(SQRT(x/nslo)) - 1.
             
             !-- fraction of activated in a bin, eq (13), part 3
             frac(ii,jj,kk) = min(1.,log(s_max/sil)/log(siu/sil))
             
             !-- Critical diameters for each bin and bin edges
             zdcrlo(ii,jj,kk) = 2.*SQRT(3.*nslo*bb/aa)
             zdcrhi(ii,jj,kk) = 2.*SQRT(3.*nshi*bb/aa)
             zdcrit(ii,jj,kk) = 2.*SQRT(3.*ns(ii,jj,kk)*bb/aa)
             
          END DO ! kk
          
          ! Find critical bin
          DO kk = in1a, fn2a
             IF (frac(ii,jj,kk) < 1. .AND. frac(ii,jj,kk) > 0.) THEN
                bcrita(ii,jj) = kk
                EXIT
             END IF
          END DO
          DO kk = in2b, fn2b
             IF (frac(ii,jj,kk) < 1. .AND. frac(ii,jj,kk) > 0.) THEN
                bcritb(ii,jj) = kk
                EXIT
             END IF
          END DO
          
       END DO ! ii
       
    END DO ! jj
    
    CALL activate3(kproma, kbdim, klev, bcrita, bcritb, &
                   zdcrit, zdcrlo, zdcrhi, zdcstar, pactd  )

  END SUBROUTINE ActCloudBase

  SUBROUTINE actInterst(kproma,kbdim,klev,prv,prs,temp)
    !
    ! Activate interstitial aerosols if they've reached their critical size
    !
    ! 1. Formulate the profiles of Dwet within bins
    !      - Get the slopes between adjacent bin mids (known)
    !      - Interpolate Dwet to bin edges
    ! 2. Based on Dcrit and Dwet slopes, estimate the Ddry where Dwet becomes > Dcrit
    ! 3. Formulate the slopes for number concentration
    ! 4. Use the Dry limits from (2) as the integration limits if they are defined
    !
    !
    USE classSection, ONLY : Section

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: kproma,kbdim,klev
    REAL, INTENT(IN) :: prv(kbdim,klev), prs(kbdim,klev)  ! Water vapour and saturation mixin ratios
    REAL, INTENT(in) :: temp(kbdim,klev)  ! Absolute temperature
    
    TYPE(Section), TARGET :: pactd(ncld) ! Local variable

    REAL, PARAMETER :: THvol = 1.e-28 ! m^3, less than the volume of a 1 nm particle

    REAL :: zkelvin               ! Coefficient for Kelvin effect

    REAL :: zdcstar, zvcstar      ! Critical diameter/volume corresponding to Smax_LES
    REAL :: vol_sol, vol_insol  ! Total mass of soluble and insoluble material within a bin

    REAL :: zactvol               ! Total volume of the activated particles
    
    REAL :: Nact, Vact(8)         ! Helper variables for transferring the activated particles
                                  ! The size is hardcoded as 8, since it's the same in the mass arrays. Not all of them are used 
                                  ! and it should be fixed at some point so that there's no unnecessary space.
    
    REAL :: Nmid, Nim1, Nip1      ! Bin number concentrations in current and adjacent bins
    REAL :: dNmid, dNim1, dNip1   ! Density function value of the number distribution for current and adjacent bins
    
    REAL :: Vmid, Vim1, Vip1      ! Dry particle volume in the middle of the bin
    REAL :: Vlo, Vhi              ! Dry particle volume scaled to bin edges
    REAL :: Vlom1, Vhim1          ! - '' - For adjacent bins
    REAL :: Vlop1, Vhip1          !
    
    REAL :: Vwmid, Vwim1, Vwip1   ! Wet particle volume in the middle of the bin
    REAL :: Vwlo, Vwhi            ! Wet particle volume at bin edges
    
    REAL :: zs1, zs2              ! Slopes for number concetration distributions within bins
    
    REAL :: N01, N02              ! Origin values for number distribution slopes
    REAL :: V01, V02              ! Origin values for wet particle volume slopes
    REAL :: Nnorm, Vnorm          ! Normalization factors for number and volume integrals
    
    REAL    :: vcut, vint1, vint2 ! cut volume, integration limit volumes
    LOGICAL :: intrange(4)        ! Logical table for integration ranges depending on the shape of the wet size profile, i.e. masks the ranges that can activate:
                                  ! [Vlo -- vint1][vint1 -- Vmid][Vmid -- vint2][vint1 -- Vhi]
    INTEGER :: cb, ab, ii, jj, ss
    INTEGER :: isol
    INTEGER :: ndry, nwet, iwa
    
    ndry = spec%getNSpec(type="dry")
    nwet = spec%getNSpec(type="wet")
    iwa = spec%getIndex("H2O")

    DO jj = 1, klev
       DO ii = 1, kbdim
          IF ( prv(ii,jj)/prs(ii,jj) <= 1.000 ) CYCLE  
          
          zkelvin = 4.*spec%mwa*surfw0/(rg*spec%rhowa*temp(ii,jj)) ! Kelvin effect [m]
            
          ! Determine Dstar == critical diameter corresponding to the host model S
          zdcstar = 2.*zkelvin/( 3.*( (prv(ii,jj)/prs(ii,jj))-1. ) )
          zvcstar = pi6*zdcstar**3
          
          ! Loop over cloud droplet (and aerosol) bins
          DO cb = ica%cur, fcb%cur
             IF (cb <= fca%cur) THEN
                ! a-bins
                ab = ica%par + (cb-ica%cur)
             ELSE
                ! b-bins
                ab = icb%par + (cb-icb%cur)
             END IF
             pactd(cb)%numc = 0.
             pactd(cb)%volc(:) =0.

             ! Determine the total mass of soluble and insoluble material in the current aerosol bin
             vol_sol = 0.
             vol_insol = 0.
             IF ( spec%Nsoluble >= 1 ) THEN ! If this isn't true, we should not be here in the first place...               
                DO isol = 1,spec%Nsoluble
                   vol_sol = vol_sol + aero(ii,jj,ab)%volc( spec%ind_soluble(isol) )
                END DO
             END IF
             IF ( spec%Ninsoluble >= 1 ) THEN
                DO isol = 1,spec%Ninsoluble
                   vol_insol = vol_insol + aero(ii,jj,ab)%volc( spec%ind_insoluble(isol) )
                END DO
             END IF

             IF ( aero(ii,jj,ab)%numc < aero(ii,jj,ab)%nlim) CYCLE                            ! Must have a reasonable number of particles
             IF ( vol_sol < aero(ii,jj,ab)%numc*THvol) CYCLE                   ! Must have at least some soluble material in the particles
             IF ( aero(ii,jj,ab)%volc(iwa) < aero(ii,jj,ab)%numc*THvol ) CYCLE ! Must not be totally dry
             
             ! Initialize the integration range mask for current bin
             intrange = .FALSE.
             
             ! Define some parameters
             Nmid = aero(ii,jj,ab)%numc                     ! Number concentration at the current bin center
             Vwmid = SUM(aero(ii,jj,ab)%volc(1:nwet))/Nmid  ! Wet volume at the current bin center
             Vmid = SUM(aero(ii,jj,ab)%volc(1:ndry))/Nmid   ! Dry volume at the current bin center
             Vlo = Vmid*aero(ii,jj,ab)%vratiolo             ! Dry vol at low limit
             Vhi = Vmid*aero(ii,jj,ab)%vratiohi             ! Dry vol at high limit
             
             ! Number concentrations and volumes at adjacent bins (check for sizedistribution boundaries)
             IF (ab == in1a .OR. ab == in2b) THEN
                Nim1 = aero(ii,jj,ab)%nlim
                Vim1 = Vlo/2.
                Vlom1 = 0.
                Vhim1 = Vlo
                Vwim1 = Vwmid/3.
             ELSE
                Nim1 = aero(ii,jj,ab-1)%numc
                IF (Nim1 > aero(ii,jj,ab-1)%nlim) THEN
                   Vim1 = SUM(aero(ii,jj,ab-1)%volc(1:ndry))/Nim1
                   Vwim1 = SUM(aero(ii,jj,ab-1)%volc(1:nwet))/Nim1
                ELSE
                   Vim1 = pi6*aero(ii,jj,ab-1)%dmid**3
                   Vwim1 = pi6*aero(ii,jj,ab-1)%dmid**3
                END IF
                Vlom1 = Vim1*aero(ii,jj,ab-1)%vratiolo
                Vhim1 = Vim1*aero(ii,jj,ab-1)%vratiohi
             END IF
             IF (ab == fn2a .OR. ab == fn2b ) THEN
                Nip1 = aero(ii,jj,ab)%nlim
                Vip1 = Vhi + 0.5*(Vhi-Vlo)
                Vlop1 = Vhi
                Vhip1 = Vhi + (Vhi-Vlo)
                Vwip1 = Vhip1
             ELSE
                Nip1 = aero(ii,jj,ab+1)%numc
                IF (Nip1 > aero(ii,jj,ab+1)%nlim) THEN
                   Vip1 = SUM(aero(ii,jj,ab+1)%volc(1:ndry))/Nip1
                   Vwip1 = SUM(aero(ii,jj,ab+1)%volc(1:nwet))/Nip1
                ELSE
                   Vip1 = pi6*aero(ii,jj,ab+1)%dmid**3
                   Vwip1 = pi6*aero(ii,jj,ab+1)%dmid**3
                END IF
                Vlop1 = Vip1*aero(ii,jj,ab+1)%vratiolo
                Vhip1 = Vip1*aero(ii,jj,ab+1)%vratiohi
             END IF
             
             ! Keeping things smooth...
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
             Nact    = 0.
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
             
             DO ss = 1, nwet
                Vact(ss) = zactvol*( aero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
             END DO
             
             ! Store the number concentration and mass of activated particles for current bins
             pactd(cb)%numc = MIN(Nact,Nmid)
             pactd(cb)%volc(1:nwet) = MIN(Vact(1:nwet),aero(ii,jj,ab)%volc(1:nwet))
             
          END DO ! cb
          
          ! Make things cleaner
          ASSOCIATE(zaer => aero(ii,jj,ica%par:fcb%par),  &
                    zcld => cloud(ii,jj,ica%cur:fcb%cur), &
                    zact => pactd(ica%cur:fcb%cur)        )
            
            ! Apply the number and mass activated to aerosol and cloud bins 
            DO cb = 1, fcb%cur - ica%cur + 1
               zaer(cb)%numc = MAX(0., zaer(cb)%numc - zact(cb)%numc)
               zcld(cb)%numc = zcld(cb)%numc + zact(cb)%numc
               DO ss = 1,nwet
                  zaer(cb)%volc(ss) = MAX(0., zaer(cb)%volc(ss) - zact(cb)%volc(ss))
                  zcld(cb)%volc(ss) = zcld(cb)%volc(ss) + zact(cb)%volc(ss)
               END DO
            END DO
            
          END ASSOCIATE
          
       END DO ! ii
       
    END DO ! jj
    
  END SUBROUTINE actInterst
  
  ! ----------------------------------------------
  
  SUBROUTINE activate3(kproma,kbdim,klev,pbcrita,pbcritb, &
                       pdcrit, pdcrlo, pdcrhi, pdcstar, pactd   )
    !
    ! Gets the number and mass activated in the critical aerosol size bin
    !
    USE classSection, ONLY : Section

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: kproma,kbdim,klev
    INTEGER, INTENT(IN) :: pbcrita(kbdim,klev),     &   ! Index of the critical aerosol bin in regime a
         pbcritb(kbdim,klev)            ! Index of the critical aerosol bin in regime b
    REAL, INTENT(IN) :: pdcrit(kbdim,klev,fn2b),    & ! Bin middle critical diameter
         pdcrlo(kbdim,klev,fn2b),    & ! Critical diameter at low limit
         pdcrhi(kbdim,klev,fn2b)       ! Critical diameter at high limit
    REAL, INTENT(IN) :: pdcstar(kbdim,klev)           ! Critical diameter corresponding to Smax
    TYPE(Section), INTENT(OUT) :: pactd(kbdim,klev,ncld) ! Properties of the maximum amount of newly activated droplets
    
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
    
    INTEGER :: ii,jj,ss,cb,ab
    
    INTEGER :: ndry,nwet,iwa
    
    ndry = spec%getNSpec(type="dry")
    nwet = spec%getNSpec(type="wet")
    iwa = spec%getIndex("H2O")
    
    DO jj = 1, klev
       DO ii = 1, kbdim
          
          ! This means in practice that vertical velocity is <= 0 or Ntot == 0
          IF ( ALL(pdcrit(ii,jj,:) < epsilon(1.0)) ) CYCLE
          
          zvcstar = 0.
          
          IF ( aero(ii,jj,pbcrita(ii,jj))%numc < aero(ii,jj,pbcrita(ii,jj))%nlim ) THEN
             Vmid = pi6*aero(ii,jj,pbcrita(ii,jj))%dmid**3
          ELSE
             Vmid = SUM( aero(ii,jj,pbcrita(ii,jj))%volc(1:ndry) )/MAX(aero(ii,jj,pbcrita(ii,jj))%nlim, &
                                                                       aero(ii,jj,pbcrita(ii,jj))%numc)
          END IF
          Vhi = aero(ii,jj,pbcrita(ii,jj))%vhilim
          Vlo = aero(ii,jj,pbcrita(ii,jj))%vlolim
          
          IF ( pdcstar(ii,jj) >= pdcrit(ii,jj,pbcrita(ii,jj)) ) THEN
             zhlp = ( pi6*pdcstar(ii,jj)**3 - pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 ) / &
                  MAX( epsilon(1.0), pi6*pdcrhi(ii,jj,pbcrita(ii,jj))**3 - pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 )
             
             zvcstar = Vmid + zhlp*(Vhi-Vmid)
          ELSE IF (pdcstar(ii,jj) < pdcrit(ii,jj,pbcrita(ii,jj)) ) THEN
             
             zhlp = ( pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 - pi6*pdcstar(ii,jj)**3 ) / &
                  MAX( epsilon(1.0), pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 - pi6*pdcrlo(ii,jj,pbcrita(ii,jj))**3 )
             
             zvcstar = Vmid - zhlp*(Vmid-Vlo)
          END IF
          
          zvcstar = MAX( zvcstar, aero(ii,jj,pbcrita(ii,jj))%vlolim )
          zvcstar = MIN( zvcstar, aero(ii,jj,pbcrita(ii,jj))%vhilim )
          
          ! Loop over cloud droplet (and aerosol) bins
          DO cb = ica%cur, fcb%cur
             IF (cb <= fca%cur) THEN
                ab = ica%par + (cb-ica%cur)
             ELSE
                ab = icb%par + (cb-icb%cur)
             END IF
             
             IF ( aero(ii,jj,ab)%numc < aero(ii,jj,ab)%nlim) CYCLE
             
             ! Formulate a slope for Wet particle size within bins and integrate over
             ! the particles larger than zvcstar
             
             Nmid = MAX(aero(ii,jj,ab)%numc, aero(ii,jj,ab)%nlim)
             Vmid = SUM(aero(ii,jj,ab)%volc(1:ndry))/Nmid ! Dry bin mid volume
             Vlo = aero(ii,jj,ab)%vlolim      ! Mid dry volume scaled to bin low limit (this is mostly an educated guess... )
             Vhi = aero(ii,jj,ab)%vhilim      ! Same for high limit
             
             IF (ab == in1a .OR. ab == in2b) THEN
                Nim1 = aero(ii,jj,ab)%nlim
                Vim1 = Vlo/2.
                Vlom1 = 0.
                Vhim1 = Vlo
             ELSE
                Nim1 = MAX(aero(ii,jj,ab-1)%numc, aero(ii,jj,ab-1)%nlim)
                IF (Nim1 > aero(ii,jj,ab-1)%nlim) THEN
                   Vim1 = SUM(aero(ii,jj,ab-1)%volc(1:ndry))/Nim1
                ELSE
                   Vim1 = pi6*aero(ii,jj,ab-1)%dmid**3
                END IF
                Vlom1 = aero(ii,jj,ab-1)%vlolim
                Vhim1 = aero(ii,jj,ab-1)%vhilim
             END IF
             IF (ab == fn2a .OR. ab == fn2b) THEN
                Nip1 = aero(ii,jj,ab)%nlim
                Vip1 = Vhi + 0.5*(Vhi-Vlo)
                Vlop1 = Vhi
                Vhip1 = Vhi + (Vhi-Vlo)
             ELSE
                Nip1 = MAX(aero(ii,jj,ab+1)%numc, aero(ii,jj,ab+1)%nlim)
                IF (Nip1 > aero(ii,jj,ab+1)%nlim) THEN
                   Vip1 = SUM(aero(ii,jj,ab+1)%volc(1:ndry))/Nip1
                ELSE
                   Vip1 = pi6*aero(ii,jj,ab+1)%dmid**3
                END IF
                Vlop1 = aero(ii,jj,ab+1)%vlolim
                Vhip1 = aero(ii,jj,ab+1)%vhilim
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
                DO ss = 1, nwet
                   pactd(ii,jj,cb)%volc(ss) = zactvol*( aero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
                END DO
                
             ELSE IF (zvcstar >= Vmid) THEN
                
                ! Use actual critical volume only in the critical bin, otherwise current bin limits
                zvcint = MIN(zvcstar,Vhi)
                
                pactd(ii,jj,cb)%numc = (Nmid/Nnorm) * ( intgN(zs2,N02,zvcint,Vhi) )
                zactvol = (Nmid*Vmid/Vnorm) * ( intgV(zs2,N02,zvcint,Vhi) )
                DO ss = 1, nwet
                   pactd(ii,jj,cb)%volc(ss) = zactvol*( aero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
                END DO
                
             END IF
             
             pactd(ii,jj,cb)%numc = MAX(0., pactd(ii,jj,cb)%numc)
             DO ss = 1, nwet
                pactd(ii,jj,cb)%volc(ss) = MAX(0., pactd(ii,jj,cb)%volc(ss))
             END DO
             
             ! "Artificially" adjust the wet size of newly activated a little bit to prevent them from being
             ! evaporated immediately
             pactd(ii,jj,cb)%volc(iwa) = pactd(ii,jj,cb)%numc*pi6*(pdcrit(ii,jj,ab)**3) *  &
                  MIN(2.,(3.e-6/max(epsilon(1.0),pdcrit(ii,jj,ab)))**2)
             
          END DO ! cb

       END DO ! ii
    END DO ! jj
    
  END SUBROUTINE activate3
  
  ! ------------------------------------------------
  REAL FUNCTION intgN(ikk,icc,ilow,ihigh)
    ! Gets the integral over a (linear) number concentration distribution
    
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
  SUBROUTINE autoconv2(kproma,kbdim,klev,ptstep  )
    !
    ! Uses a more straightforward method for converting cloud droplets to drizzle.
    ! Assume a lognormal cloud droplet distribution for each bin. Sigma_g is an adjustable
    ! parameter and is set to 1.2 by default
    !

    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: kproma,kbdim,klev
    REAL, INTENT(in)    :: ptstep
    
    REAL :: Vrem, Nrem, Vtot, Ntot
    REAL :: dvg,dg
    REAL :: tot
    
    REAL, PARAMETER :: zd0 = 50.e-6
    REAL, PARAMETER :: sigmag = 1.2
    REAL, PARAMETER :: max_rate_autoc=1.0e10 ! Maximum autoconversion rate (#/m^3/s)
    
    INTEGER :: ii,jj,cc,ss
    INTEGER :: nwet, ndry, iwa
    
    nwet = spec%getNSpec(type='wet')
    ndry = spec%getNSpec(type='dry')
    iwa = spec%getIndex("H2O")
        
    ! Find the cloud bins where the mean droplet diameter is above 50 um
    ! Do some fitting...
    DO jj = 1, klev
       DO ii = 1, kbdim
          DO cc = ncld,1,-1 ! Start from the largest drops
             ! Autoconversion rate can be limited
             tot = 0.
             
             Ntot = cloud(ii,jj,cc)%numc
             Vtot = SUM(cloud(ii,jj,cc)%volc(1:nwet))
             
             IF ( Ntot > cloud(ii,jj,cc)%nlim .AND. Vtot > 0. ) THEN
                
                ! Volume geometric mean diameter
                dvg = ((Vtot/Ntot/pi6)**(1./3.))*EXP( (3.*LOG(sigmag)**2)/2. )
                dg = dvg*EXP( -3.*LOG(sigmag)**2 )
                
                Vrem = Vtot*( 1. - cumlognorm(dvg,sigmag,zd0) )
                Nrem = Ntot*( 1. - cumlognorm(dg,sigmag,zd0) )
                
                IF ( Vrem > 0. .AND. Nrem > precp(ii,jj,1)%nlim) THEN
                   
                   ! Put the mass and number to the first precipitation bin and remove from
                   ! cloud droplets
                   DO ss = 1, ndry
                      precp(ii,jj,1)%volc(ss) = precp(ii,jj,1)%volc(ss) + cloud(ii,jj,cc)%volc(ss)*(Nrem/Ntot)
                      cloud(ii,jj,cc)%volc(ss) = cloud(ii,jj,cc)%volc(ss)*(1. - (Nrem/Ntot))
                   END DO
                   
                   precp(ii,jj,1)%volc(iwa) = precp(ii,jj,1)%volc(iwa) + cloud(ii,jj,cc)%volc(iwa)*(Vrem/Vtot)
                   cloud(ii,jj,cc)%volc(iwa) = cloud(ii,jj,cc)%volc(iwa)*(1. - (Vrem/Vtot))
                   
                   precp(ii,jj,1)%numc = precp(ii,jj,1)%numc + Nrem
                   cloud(ii,jj,cc)%numc = cloud(ii,jj,cc)%numc - Nrem
                   
                   tot = tot + Nrem
                   IF (tot > max_rate_autoc*ptstep) THEN
                      WRITE(*,*) tot, max_rate_autoc
                      EXIT
                   END IF
                   
                END IF ! Nrem Vrem
                
             END IF ! Ntot Vtot
             
          END DO ! cc
       END DO ! ii
    END DO ! jj
    
  END SUBROUTINE autoconv2
  

  
END MODULE mo_salsa_cloud
