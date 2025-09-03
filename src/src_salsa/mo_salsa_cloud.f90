MODULE mo_salsa_cloud

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

  SUBROUTINE cloud_activation(kbdim, klev,           &
                              temp, rv, rs, paero, pcloud)

    USE mo_submctl, ONLY : t_section, nbins, ncld

    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) :: kbdim, klev      ! Dimensions
    REAL, INTENT(in)    :: temp(kbdim,klev) ! Temperature
    REAL, INTENT(inout) :: rv(kbdim,klev) ! Water vapor mixing ratio
    REAL, INTENT(in)    :: rs(kbdim,klev) ! Saturation vapor mixing ratio
    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld),  &
                                      paero(kbdim,klev,nbins)

    ! -------------------------------------
    ! Interstitial activation
    ! -------------------------------------
    CALL actInterst(kbdim,klev,paero,pcloud,rv,rs,temp)

    ! -----------------------------------
    ! Activation at cloud base
    ! -----------------------------------
    ! Cloud base activation not supported!

  END SUBROUTINE cloud_activation
! -----------------------------------------------------------------


  SUBROUTINE actInterst(kbdim,klev,paero,pcloud,prv,prs,temp)
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
    ! Note: insoluble species are not properly accounted for
    !
    USE mo_submctl, ONLY :     &
         rg,surfw0,nlim,pi6,mwa,rhowa,t_section, & ! constants
         in1a,fn1a,in2a,fn2a,in2b,fn2b,inp2a,fnp2b, & ! bin indexes
         nbins,ncld,nspec
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kbdim,klev
    TYPE(t_section), INTENT(INOUT) :: paero(kbdim,klev,nbins),  &
                                      pcloud(kbdim,klev,ncld)
    REAL, INTENT(IN) :: prv(kbdim,klev),prs(kbdim,klev)  ! Water vapour and saturation mixin ratios
    REAL, INTENT(in) :: temp(kbdim,klev)  ! Absolute temperature

    REAL :: paa        ! Coefficient for Kelvin effect

    REAL :: zdcstar,zvcstar   ! Critical diameter/volume corresponding to S_LES
    REAL :: zactvol,zactn     ! Total volume of the activated particles

    REAL :: Nact(ncld), Vact(ncld,nspec+1) ! Helper variables for transferring the activated particles

    REAL :: Nmid, Nim1, Nip1     ! Bin number concentrations in current and adjacent bins
    REAL :: dNmid, dNim1, dNip1  ! Density function value of the number distribution for current and adjacent bins

    REAL :: Vmid, Vim1, Vip1     ! Dry particle volume in the middle of the bin
    REAL :: Vlo, Vhi             ! Dry particle volume scaled to bin edges
    REAL :: Vlom1,Vhim1          ! - '' - For adjacent bins
    REAL :: Vlop1,Vhip1          !

    REAL :: Vwmid, Vwim1, Vwip1  ! Wet particle volume in the middle of the bin
    REAL :: Vwlo,Vwhi            ! Wet particle volume at bin edges

    REAL :: zs1,zs2           ! Slopes for number concetration distributions within bins

    REAL :: N01,N02           ! Origin values for number distribution slopes
    REAL :: V01,V02           ! Origin values for wet particle volume slopes
    REAL :: Nnorm, Vnorm      ! Normalization factors for number and volume integrals

    REAL :: vint1,vint2       ! integration limit volumes
    INTEGER :: cb,ab, ii,jj,ss,nn

    nn = nspec+1 ! Aerosol species + water

    DO jj = 1,klev
       DO ii = 1,kbdim
          IF ( prv(ii,jj)/prs(ii,jj) <= 1.000 ) CYCLE

          paa = 4.*mwa*surfw0/(rg*rhowa*temp(ii,jj)) ! Kelvin effect [m]

          ! Determine Dstar == critical droplet diameter corresponding to the host model S
          zdcstar = 2.*paa/( 3.*( (prv(ii,jj)/prs(ii,jj))-1. ) )
          zvcstar = pi6*zdcstar**3

          ! Note: this is valid for dilute droplets, which include all particles containing soluble
          ! substances. Critical droplet diameter for a completely insoluble particle (DU or BC) is
          !   zdstar = paa/(prv(ii,jj)/prs(ii,jj)-1.)
          ! but dry particles should not be moved to cloud bins. Otherwise just assume that there
          ! is enough soluble material so that the current approach is valid.

          Nact(:) = 0.
          Vact(:,:) = 0.

          ! Loop over cloud droplet (and aerosol) bins
          DO cb = inp2a, fnp2b
             ab = fn1a + cb
             ! Dry particles are not activated (volume 1e-28 m^3 is less than that in a 1 nm droplet)
             IF ( paero(ii,jj,ab)%numc < nlim .OR. paero(ii,jj,ab)%volc(1)<paero(ii,jj,ab)%numc*1e-28 ) CYCLE

             ! Define some parameters
             Nmid = paero(ii,jj,ab)%numc     ! Number concentration at the current bin center
             Vwmid = SUM(paero(ii,jj,ab)%volc(1:nn))/Nmid  ! Wet volume at the current bin center
             Vmid = SUM(paero(ii,jj,ab)%volc(2:nn))/Nmid ! Dry volume at the current bin center
             Vlo = Vmid*paero(ii,jj,ab)%vlolim/paero(ii,jj,ab)%vmid ! Dry vol at low limit
             Vhi = Vmid*paero(ii,jj,ab)%vhilim/paero(ii,jj,ab)%vmid ! Dry vol at high limit

             ! Number concentrations and volumes at adjacent bins (check for sizedistribution boundaries)
             IF (ab==in1a .OR. ab==in2b) THEN
                   Nim1 = nlim
                   Vim1 = Vlo/2.
                   Vlom1 = 0.
                   Vhim1 = Vlo
                   Vwim1 = Vwmid/3.
             ELSE
                   Nim1 = paero(ii,jj,ab-1)%numc
                   IF (Nim1 > nlim) THEN
                      Vim1 = SUM(paero(ii,jj,ab-1)%volc(2:nn))/Nim1
                      Vwim1 = SUM(paero(ii,jj,ab-1)%volc(1:nn))/Nim1
                   ELSE
                      Vim1 = paero(ii,jj,ab-1)%vmid
                      Vwim1 = Vim1
                   END IF
                   Vlom1 = Vim1*paero(ii,jj,ab-1)%vlolim/paero(ii,jj,ab-1)%vmid
                   Vhim1 = Vim1*paero(ii,jj,ab-1)%vhilim/paero(ii,jj,ab-1)%vmid
             END IF
             IF (ab==fn2a .OR. ab==fn2b ) THEN
                   Nip1 = nlim
                   Vip1 = Vhi + 0.5*(Vhi-Vlo)
                   Vlop1 = Vhi
                   Vhip1 = Vhi + (Vhi-Vlo)
                   Vwip1 = Vhip1
             ELSE
                   Nip1 = paero(ii,jj,ab+1)%numc
                   IF (Nip1 > nlim) THEN
                      Vip1 = SUM(paero(ii,jj,ab+1)%volc(2:nn))/Nip1
                      Vwip1 = SUM(paero(ii,jj,ab+1)%volc(1:nn))/Nip1
                   ELSE
                      Vip1 = paero(ii,jj,ab+1)%vmid
                      Vwip1 = Vip1
                   END IF
                   Vlop1 = Vip1*paero(ii,jj,ab+1)%vlolim/paero(ii,jj,ab+1)%vmid
                   Vhip1 = Vip1*paero(ii,jj,ab+1)%vhilim/paero(ii,jj,ab+1)%vmid
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
                Nact(cb) = paero(ii,jj,ab)%numc
                Vact(cb,1:nn) = paero(ii,jj,ab)%volc(1:nn)
             ELSE IF ( zvcstar > Vwlo .AND. zvcstar > Vwmid .AND. zvcstar > Vwhi) THEN
                ! None activates
             ELSE
                ! Partial activation:
                 ! Slope1
                vint1 = (zvcstar - V01)/zs1  ! Where the wet size profile intersects the critical size (slope 1)
                IF (vint1 < Vlo .OR. vint1 > Vmid) THEN
                   ! intersection volume outside the current size range -> set as the lower limit
                   vint1 = Vlo
                END IF

                ! Slope2
                vint2 = (zvcstar - V02)/zs2  ! Where the wet size profile intersects the critical size (slope 2)
                IF (vint2 < Vmid .OR. vint2 > Vhi) THEN
                   ! Intersection volume outside the current size range -> set as the lower limit
                   vint2 = Vmid
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
                zactn = 0.
                zactvol = 0.

                ! Integration over each size range within a bin
                IF ( Vwlo > zvcstar ) THEN
                   zactn = zactn + intgN(zs1,N01,Vlo,vint1)
                   zactvol = zactvol + intgV(zs1,N01,Vlo,vint1)
                END IF

                IF ( Vwmid > zvcstar ) THEN
                   zactn = zactn + intgN(zs1,N01,vint1,Vmid)
                   zactvol = zactvol + intgV(zs1,N01,vint1,Vmid)
                   zactn = zactn + intgN(zs2,N02,Vmid,vint2)
                   zactvol = zactvol + intgV(zs2,N02,Vmid,vint2)
                END IF

                IF ( Vwhi > zvcstar ) THEN
                   zactn = zactn + intgN(zs2,N02,vint2,Vhi)
                   zactvol = zactvol + intgV(zs2,N02,vint2,Vhi)
                END IF

                ! Store the number concentration and mass of activated particles for current bins
                Nact(cb) = MIN(zactn/Nnorm*paero(ii,jj,ab)%numc,paero(ii,jj,ab)%numc)
                Vact(cb,1:nn) = MIN(zactvol/Vnorm*paero(ii,jj,ab)%volc(1:nn),paero(ii,jj,ab)%volc(1:nn))
             END IF

          END DO ! cb

          ! Apply the number and mass activated to aerosol and cloud bins
          paero(ii,jj,in2a:)%numc = MAX(0., paero(ii,jj,in2a:)%numc - Nact(1:))
          pcloud(ii,jj,:)%numc = pcloud(ii,jj,:)%numc + Nact(:)
          DO ss = 1,nn
             paero(ii,jj,in2a:)%volc(ss) = MAX(0., paero(ii,jj,in2a:)%volc(ss) - Vact(1:,ss))
             pcloud(ii,jj,:)%volc(ss) = pcloud(ii,jj,:)%volc(ss) + Vact(:,ss)
          END DO

       END DO ! ii

    END DO ! jj

  END SUBROUTINE actInterst


  ! ----------------------------------------------
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


  ! ------------------------------------------------
  SUBROUTINE ReleaseDrops(kbdim,klev,paero,pcloud,pprecp,prv,prs,ptemp)
    !
    ! Release cloud and rain drops back to aerosol when they have become small enough
    !
    USE mo_submctl, ONLY :  t_section,nbins,ncld,nprc,rg,surfw0,pi6,pi,nlim,prlim, &
                diss,dens,mws,fn1a,in2a,fn2a,in2b,nspec,calc_correlation
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kbdim,klev
    TYPE(t_section), INTENT(INOUT) :: paero(kbdim,klev,nbins), &
                pcloud(kbdim,klev,ncld),pprecp(kbdim,klev,nprc)
    REAL, INTENT(IN) :: prv(kbdim,klev),prs(kbdim,klev),ptemp(kbdim,klev)
    REAL :: paa, ns, cd, zvol, ra, rb, vol
    INTEGER :: ii,jj,kk,ab,bb,nn

    nn = nspec+1 ! Aerosol species + water

    DO jj = 1,klev
      DO ii = 1,kbdim
        IF ( prv(ii,jj)/prs(ii,jj) >= 0.999 ) CYCLE
        ! For calculating critical droplet diameter
        paa =rg*ptemp(ii,jj)/(2.*pi*surfw0)

        ! Loop over cloud droplet (and aerosol) bins
        DO kk = 1,ncld
            IF ( pcloud(ii,jj,kk)%numc>nlim .AND. pcloud(ii,jj,kk)%volc(1)<1e-8 ) THEN
                ! Critical droplet diameter
                ns = SUM( pcloud(ii,jj,kk)%volc(2:nn)*diss(2:nn)*dens(2:nn)/mws(2:nn) )/pcloud(ii,jj,kk)%numc
                cd = 3.*SQRT(ns*paa)
                ! Wet diameter
                zvol = (SUM(pcloud(ii,jj,kk)%volc(1:nn))/pcloud(ii,jj,kk)%numc/pi6)**(1./3.)
                ! Lose the droplets if smaller than 0.2*critical diameter or 2 um or if there is no water
                IF ( zvol < MAX(0.2*cd,2.e-6) .OR. pcloud(ii,jj,kk)%volc(1)<1e-28*pcloud(ii,jj,kk)%numc ) THEN
                    ab = fn1a + kk ! Index for parallel aerosol bin
                    ! Move the number of particles from cloud to aerosol bins
                    paero(ii,jj,ab)%numc = paero(ii,jj,ab)%numc + pcloud(ii,jj,kk)%numc
                    pcloud(ii,jj,kk)%numc = 0.0
                    ! Move ccn material back to aerosol regime (including water)
                    paero(ii,jj,ab)%volc(1:nn) = paero(ii,jj,ab)%volc(1:nn) + pcloud(ii,jj,kk)%volc(1:nn)
                    pcloud(ii,jj,kk)%volc(1:nn) = 0.0
                END IF
            END IF
        END DO

        ! Loop over precipitation bins
        DO kk = 1,nprc
            IF ( pprecp(ii,jj,kk)%numc>prlim .AND. pprecp(ii,jj,kk)%volc(1)<1e-9 ) THEN
                ! Critical droplet diameter
                ns = SUM( pprecp(ii,jj,kk)%volc(2:nn)*diss(2:nn)*dens(2:nn)/mws(2:nn) )/pprecp(ii,jj,kk)%numc
                cd = 3.*SQRT(ns*paa)
                ! Wet diameter
                zvol = (SUM(pprecp(ii,jj,kk)%volc(1:nn))/pprecp(ii,jj,kk)%numc/pi6)**(1./3.)

               ! Lose the droplets if smaller than 0.02*critical diameter or 2 um or if there is no water
               IF ( zvol < MAX(0.02*cd,2.e-6) .OR. pprecp(ii,jj,kk)%volc(1)<1e-28*pprecp(ii,jj,kk)%numc ) THEN
                    ! Move evaporating precipitation to aerosol bin based on dry radius and chemical composition

                    ! 1) Matching a-bin
                    vol=SUM(pprecp(ii,jj,kk)%volc(2:nn))/pprecp(ii,jj,kk)%numc
                    ab=MAX(1,COUNT(vol>paero(ii,jj,1:fn2a)%vlolim)) ! Aerosol a-bin
                    ! Corresponding b bin
                    bb=ab-in2a+in2b
                    ! 2) Select a or b bin
                    IF (ab<in2a .OR. paero(ii,jj,bb)%numc<=nlim) THEN
                        ! Empty b bin so select a
                        !ab = ab
                    ELSEIF (paero(ii,jj,ab)%numc<=nlim) THEN
                        ! Empty a bin so select b
                        ab = bb
                    ELSE
                        ! Both are present - find bin based on compositional similarity
                        ra = calc_correlation(paero(ii,jj,ab)%volc(2:nn),pprecp(ii,jj,kk)%volc(2:nn),nspec)
                        rb = calc_correlation(paero(ii,jj,bb)%volc(2:nn),pprecp(ii,jj,kk)%volc(2:nn),nspec)
                        IF (ra<rb) ab = bb
                    ENDIF

                    ! Move the number of particles from rain to aerosol bins
                    paero(ii,jj,ab)%numc = paero(ii,jj,ab)%numc + pprecp(ii,jj,kk)%numc
                    pprecp(ii,jj,kk)%numc = 0.0
                    ! Move ccn material back to aerosol regime (including water)
                    paero(ii,jj,ab)%volc(1:nn) = paero(ii,jj,ab)%volc(1:nn) + pprecp(ii,jj,kk)%volc(1:nn)
                    pprecp(ii,jj,kk)%volc(1:nn) = 0.0
                END IF
            END IF
        END DO
      END DO
    END DO

  END SUBROUTINE ReleaseDrops

  ! ------------------------------------------------
  SUBROUTINE ReleaseIce(kbdim,klev,paero,pice,psnow,prv,prsi)
    !
    ! Release ice and snow back to aerosol when they have become small enough
    !
    USE mo_submctl, ONLY : t_section,nbins,nice,nsnw,pi6,nlim,prlim, &
                fn1a,in2a,fn2a,in2b,nspec,calc_correlation
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kbdim,klev
    TYPE(t_section), INTENT(INOUT) :: paero(kbdim,klev,nbins), &
                pice(kbdim,klev,nice),psnow(kbdim,klev,nsnw)
    REAL, INTENT(IN) :: prv(kbdim,klev),prsi(kbdim,klev)
    REAL :: diam, vrat, zvol, ra, rb
    INTEGER :: ii,jj,kk,ab,bb,nn

    nn = nspec+1 ! Aerosol species + water

    DO jj = 1,klev
      DO ii = 1,kbdim
        IF ( prv(ii,jj)/prsi(ii,jj) >= 0.999 ) CYCLE

        ! Loop over ice (and aerosol) bins
        DO kk = 1, nice
            IF ( pice(ii,jj,kk)%numc>prlim .AND. pice(ii,jj,kk)%volc(1)<1e-11 ) THEN
                ! Diameter (assuming water density for ice)
                diam = (SUM(pice(ii,jj,kk)%volc(1:nn))/pice(ii,jj,kk)%numc/pi6)**(1./3.)
                ! Dry to total volume ratio
                vrat = SUM(pice(ii,jj,kk)%volc(2:nn))/SUM(pice(ii,jj,kk)%volc(1:nn))
                ! Ice and snow don't have a critical size, but lose particles smaller than 2e-6 m
                ! and particles which dry to total mass ratio is more than 0.5
                IF ( vrat>0.5 .OR. diam<2e-6 ) THEN
                    ab = fn1a + kk ! Index for parallel aerosol bin
                    ! Move the number of particles from ice to aerosol bins
                    paero(ii,jj,ab)%numc = paero(ii,jj,ab)%numc + pice(ii,jj,kk)%numc
                    pice(ii,jj,kk)%numc = 0.0
                    ! Move ccn material back to aerosol regime (including water)
                    paero(ii,jj,ab)%volc(1:nn) = paero(ii,jj,ab)%volc(1:nn) + pice(ii,jj,kk)%volc(1:nn)
                    pice(ii,jj,kk)%volc(1:nn) = 0.0
                END IF
            END IF
        END DO

        ! Loop over snow bins
        DO kk = 1,nsnw
            IF ( psnow(ii,jj,kk)%numc > prlim .AND. psnow(ii,jj,kk)%volc(1)<1e-11 ) THEN
                ! Diameter (assuming water density for snow)
                diam =(SUM(psnow(ii,jj,kk)%volc(1:nn))/psnow(ii,jj,kk)%numc/pi6)**(1./3.)
                ! Dry to total volume ratio
                vrat = SUM(psnow(ii,jj,kk)%volc(2:nn))/SUM(psnow(ii,jj,kk)%volc(1:nn))
                ! Lose particles smaller than 2e-6 m and particles which dry to total
                ! mass ratio is more than 0.5
                IF ( vrat>0.5  .OR. diam<2.e-6 ) THEN
                    ! Move snow to aerosol bin based on dry radius and chemical composition

                    ! 1) Matching a-bin
                    zvol=SUM(psnow(ii,jj,kk)%volc(2:nn))/psnow(ii,jj,kk)%numc
                    ab=MAX(1,COUNT(zvol>paero(ii,jj,1:fn2a)%vlolim)) ! Aerosol a-bin
                    ! Corresponding b bin
                    bb=ab-in2a+in2b
                    ! 2) Select a or b bin
                    IF (ab<in2a .OR. paero(ii,jj,bb)%numc<=nlim) THEN
                        ! Empty b bin or 1a, so select a
                        !ab = ab
                    ELSEIF (paero(ii,jj,ab)%numc<=nlim) THEN
                        ! Empty a bin so select b
                        ab = bb
                    ELSE
                        ! Both are present - find bin based on compositional similarity
                        ra = calc_correlation(paero(ii,jj,ab)%volc(2:nn),psnow(ii,jj,kk)%volc(2:nn),nspec)
                        rb = calc_correlation(paero(ii,jj,bb)%volc(2:nn),psnow(ii,jj,kk)%volc(2:nn),nspec)
                        IF (ra<rb) ab = bb
                    ENDIF

                    ! Move the number of particles from rain to aerosol bins
                    paero(ii,jj,ab)%numc = paero(ii,jj,ab)%numc + psnow(ii,jj,kk)%numc
                    psnow(ii,jj,kk)%numc = 0.0
                    ! Move ccn material back to aerosol regime (including water)
                    paero(ii,jj,ab)%volc(1:nn) = paero(ii,jj,ab)%volc(1:nn) + psnow(ii,jj,kk)%volc(1:nn)
                    psnow(ii,jj,kk)%volc(1:nn) = 0.0
                END IF
            END IF
        END DO
      END DO
    END DO

  END SUBROUTINE ReleaseIce

  ! ------------------------------------------------
  SUBROUTINE ReleaseAerosol(kbdim,klev,paero,pc_gas,ngas)
    !
    ! Remove or release aerosol back to gas phase when they have become small enough
    !
    USE mo_submctl, ONLY : t_section, nbins, nspec, pi6, nlim, &
            lscndgas, part_h2so4, isog, iso, part_ocnv, iocg, ioc
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kbdim,klev,ngas
    TYPE(t_section), INTENT(INOUT) :: paero(kbdim,klev,nbins)
    REAL, INTENT(INOUT) :: pc_gas(kbdim,klev,ngas)
    REAL ::  zvol
    INTEGER :: ii,jj,kk,nn

    nn = nspec+1 ! Aerosol species + water

    DO jj = 1,klev
      DO ii = 1,kbdim
        ! Loop over aerosol bins
        DO kk = 1,nbins
            IF (paero(ii,jj,kk)%numc > nlim) THEN
                ! Dry volume
                zvol = SUM(paero(ii,jj,kk)%volc(2:nn))/paero(ii,jj,kk)%numc

                ! Particles smaller than 0.1 nm diameter are set to zero
                IF ( zvol < pi6*1.e-30 ) THEN
                    ! Volatile species to the gas phase
                    IF (lscndgas .AND. part_h2so4 .AND. isog>0 .AND. iso>0) THEN
                        pc_gas(kbdim,klev,isog) = pc_gas(kbdim,klev,isog) + paero(ii,jj,kk)%volc(iso)
                    END IF
                    IF (lscndgas .AND. part_ocnv .AND. iocg>0 .AND. ioc>0) THEN
                        pc_gas(kbdim,klev,iocg) = pc_gas(kbdim,klev,iocg) + paero(ii,jj,kk)%volc(ioc)
                    END IF

                    ! Mass and number to zero (insolube species and water are lost)
                    paero(ii,jj,kk)%volc(1:nn) = 0.
                    paero(ii,jj,kk)%numc = 0.
                END IF
            END IF
        END DO
      END DO
    END DO

  END SUBROUTINE ReleaseAerosol


  SUBROUTINE clean_missing(kbdim,klev,n,psect)
    !
    ! Remove negative number concentration values
    ! and particles that have number but no volume (or mass).
    !
    USE mo_submctl, ONLY : t_section
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kbdim,klev,n
    TYPE(t_section), INTENT(INOUT) :: psect(kbdim,klev,n)
    INTEGER :: ii,jj,kk

    DO jj = 1,klev
      DO ii = 1,kbdim
        DO kk = 1,n
            psect(ii,jj,kk)%numc = MAX(0.0,psect(ii,jj,kk)%numc)
            IF (psect(ii,jj,kk)%numc > 0. .AND. SUM(psect(ii,jj,kk)%volc(:)) <= 0.) THEN
                psect(ii,jj,kk)%numc = 0.0
                psect(ii,jj,kk)%volc(:) = 0.0
            END IF
        END DO
      END DO
    END DO

  END SUBROUTINE clean_missing



  !-----------------------------------------
  SUBROUTINE autoconv2(kbdim,klev,   &
                      pcloud,pprecp)
  !
  ! Uses a more straightforward method for converting cloud droplets to drizzle.
  ! Assume a lognormal cloud droplet distribution for each bin. Sigma_g is an adjustable
  ! parameter and is set to 1.2 by default
  !
    
    USE mo_submctl, ONLY : t_section,   &
                               ncld,        &
                               nprc,        &
                               nspec,       &
                               pi6,         &
                               nlim, prlim, &
                               autoc_rain_zd0, autoc_rain_sigmag
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,klev
    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld)
    TYPE(t_section), INTENT(inout) :: pprecp(kbdim,klev,nprc)

    REAL :: Vrem, Nrem, Vtot, Ntot
    REAL :: dvg,dg
    INTEGER :: ii,jj,cc,ss

    DO jj = 1,klev
       DO ii = 1,kbdim
          DO cc = 1,ncld

             Ntot = pcloud(ii,jj,cc)%numc
             Vtot = SUM(pcloud(ii,jj,cc)%volc(:))

             IF ( Ntot > nlim .AND. Vtot > 0. ) THEN

                ! Volume geometric mean diameter
                dvg = ((Vtot/Ntot/pi6)**(1./3.))*EXP( (3.*LOG(autoc_rain_sigmag)**2)/2. )
                dg = dvg*EXP( -3.*LOG(autoc_rain_sigmag)**2 )

                Vrem = Vtot*( 1. - cumlognorm(dvg,autoc_rain_sigmag,autoc_rain_zd0) )
                Nrem = Ntot*( 1. - cumlognorm(dg,autoc_rain_sigmag,autoc_rain_zd0) )

                IF ( Vrem > 0. .AND. Nrem > prlim) THEN

                   ! Put the mass and number to the first precipitation bin and remove from
                   ! cloud droplets
                   DO ss = 2,nspec+1
                      pprecp(ii,jj,1)%volc(ss) = pprecp(ii,jj,1)%volc(ss) + pcloud(ii,jj,cc)%volc(ss)*(Nrem/Ntot)
                      pcloud(ii,jj,cc)%volc(ss) = pcloud(ii,jj,cc)%volc(ss)*(1. - (Nrem/Ntot))
                   END DO
                   
                   pprecp(ii,jj,1)%volc(1) = pprecp(ii,jj,1)%volc(1) + pcloud(ii,jj,cc)%volc(1)*(Vrem/Vtot)
                   pcloud(ii,jj,cc)%volc(1) = pcloud(ii,jj,cc)%volc(1)*(1. - (Vrem/Vtot))

                   pprecp(ii,jj,1)%numc = pprecp(ii,jj,1)%numc + Nrem
                   pcloud(ii,jj,cc)%numc = pcloud(ii,jj,cc)%numc - Nrem

                END IF ! Nrem Vrem

             END IF ! Ntot Vtot

          END DO ! cc
       END DO ! ii
    END DO ! jj

  END SUBROUTINE autoconv2


  !-----------------------------------------
  SUBROUTINE autoconv_sb(kbdim,klev,ptstep,pcloud,pprecp)
    !
    ! Autoconversion based on the Seifert & Beheng (2001) microphysics
    !
    ! Seifert, A., and Beheng, K.D.: A two-moment cloud microphysics parameterization
    ! for mixed-phase clouds. Part 1: Model description, Meteorol. Atmos. Phys., 92, 45-66, 2006.
    ! DOI: 10.1007/s00703-005-0112-4
    !
    USE mo_submctl, ONLY : t_section,   &
                               ncld,        &
                               nprc,        &
                               nlim, prlim, &
                               rhowa
    IMPLICIT NONE

    ! Parameters
    real, parameter :: &
            nu_c  = 0, &       ! Width parameter of cloud DSD
            k_c = 9.44e+9, &   ! Long-Kernel (m^3/kg^2/s)
            k_1 = 6.e+2, k_2 = 0.68, &  ! Parameters for phi function
            X_min = 4.2e-15, & ! Minimum cloud droplet mass (kg), D_min=2e-6 m
            X_bnd = 2.6e-10 ! Droplet mass that separates cloud droplets from rain drops (kg), D_bnd=80e-6 m
    INTEGER, SAVE :: iout=-1  ! Specified bin for new raindrops: <=0 find the bin matching with X_bnd; >0 select one bin (hard coded)

    INTEGER, INTENT(in) :: kbdim,klev
    REAL, INTENT(IN) :: ptstep
    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld)
    TYPE(t_section), INTENT(inout) :: pprecp(kbdim,klev,nprc)

    REAL :: Vrem, Nrem, Vtot, Ntot, fact, scaling
    REAL :: k_au, Xc, tau, phi, au
    INTEGER :: ii,jj,cc,io

    IF (iout<=0) THEN
        ! Find the bin matching with X_bnd
        io=1
        DO WHILE (io<nprc .AND. X_bnd/rhowa .GT. pprecp(1,1,io)%vhilim)
            io=io+1
        END DO
        iout=io
    ENDIF

    ! Constant for the autoconversion equation (m^3/kg^3/s)
    k_au  = k_c / (20.*X_bnd) * (nu_c+2.)*(nu_c+4.)/(nu_c+1.)**2

    DO jj = 1,klev
      DO ii = 1,kbdim
        ! Total cloud droplet number and water volume
        Ntot = SUM(pcloud(ii,jj,:)%numc)
        Vtot = SUM(pcloud(ii,jj,:)%volc(1))
        IF (Ntot <= nlim .OR. Vtot*rhowa <= Ntot*X_min) CYCLE
        !
        ! Parameters
        if (Vtot > 1.e-10) then ! Limit set to 1e-10 m^3/m^3
            !   Dimensionless internal time scale - based on water only
            tau = 1.0-Vtot/( Vtot+SUM(pprecp(ii,jj,:)%volc(1)) )
            tau = MIN(MAX(tau,epsilon(1.0)),0.9)
            !   Universal function for autoconversion
            phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
        ELSE
            tau = 0.
            phi = 0.
        ENDIF
        !
        ! Mean cloud droplet mass (kg) - assuming it is just water
        Xc = MIN(MAX(Vtot/Ntot*rhowa,X_min),X_bnd)
        !
        ! Autoconversion rate (kg/m^3/s)
        au = k_au*(Vtot*rhowa)**2*Xc**2*(1.0+phi/(1.0-tau)**2)
        !
        ! Convert to absolute change in rain water volume mixing ratio (m^3/m^3)
        Vrem = MAX(au/rhowa*ptstep,0.0)
        !
        ! Detemine the change in rain droplet number concentration
        !   Note: Nrem=dNrain=fact*dNcloud where 0<fact<=1 represents number reduction due to cloud droplet collisions
        IF (Vrem>=Vtot) THEN
            ! If all volume should go to the rain bins, then just move all cloud droplets
            fact=0.5
            Nrem=Ntot*fact
            Vrem=Vtot
            ! Find bin
            io=1
            DO WHILE (io<nprc .AND. Vrem .GT. Nrem*pprecp(ii,jj,io)%vhilim)
               io=io+1
            END DO
            ! Note: must move all to this bin
        ELSE
            ! Add rain drops to the specified bin (e.g. the first bin, the bin with 65-100 um rain drops, or any other bin) and
            ! calculate number based on that: N*Vmid=V => N=V/Vmid
            io=iout ! The default bin
            Nrem=Vrem/pprecp(ii,jj,io)%vmid
            ! Determine fact: Nrem=SUM(fact*Nc(i)*Vrem/Vtot)=fact*Ntot*Vrem/Vtot
            fact=MIN(1.,Nrem/Ntot*Vtot/Vrem) ! Limited to 1
        END IF
        !
        IF (Vrem <= 0. .OR. Nrem <= prlim) CYCLE
        !
        ! How to split the change between bins?
        !   - Don't change cloud droplet dry size (bins), so volc(2:) and numc must change with the same fraction
        !   - Should not change aerosol-water ratios, so volc(:) and numc must change with the same fraction
        !   - Number need not to be conserved in the coagulation-based autoconversion (adjusting rain drop number)
        !   - Linear scaling, dV(i)=V(i)*dVtot/Vtot, is robust and does not produce unphysical results
        !
        scaling = Vrem/Vtot
        !
        DO cc = 1,ncld
            IF ( pcloud(ii,jj,cc)%numc > nlim ) THEN
                pprecp(ii,jj,io)%volc(:) = pprecp(ii,jj,io)%volc(:) + pcloud(ii,jj,cc)%volc(:)*scaling
                pcloud(ii,jj,cc)%volc(:) = pcloud(ii,jj,cc)%volc(:)*(1. - scaling)
                pprecp(ii,jj,io)%numc = pprecp(ii,jj,io)%numc + pcloud(ii,jj,cc)%numc*scaling*fact
                pcloud(ii,jj,cc)%numc = pcloud(ii,jj,cc)%numc*(1. - scaling)
            END IF
        END DO ! cc
      END DO ! ii
    END DO ! jj

  END SUBROUTINE autoconv_sb


  !***********************************************
  ! Ice nucleation
  !***********************************************

  SUBROUTINE ice_nucl_driver(kbdim,klev,   &
                      paero,pcloud,pprecp,pice,psnow, &
                      ptemp,prv,prs,prsi,ptstep)

    USE mo_submctl, ONLY : t_section,   &
                               in2a, fn2b, &
                               ncld, nprc, nice, nsnw,  &
                               pi, nlim, prlim, &
                               calc_Sw_eq, &
                               ice_hom, ice_imm, ice_dep, &
                               calc_correlation, &
                               idu, nspec

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,klev
    REAL, INTENT(in) :: ptstep
    REAL, INTENT(in) :: ptemp(kbdim,klev),  &
                            prv(kbdim,klev),    &
                            prs(kbdim,klev),    &
                            prsi(kbdim,klev)

    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,fn2b), &
                                      pcloud(kbdim,klev,ncld), &
                                      pprecp(kbdim,klev,nprc), &
                                      pice(kbdim,klev,nice), &
                                      psnow(kbdim,klev,nsnw)

    ! Which species are allowed to freeze
    LOGICAL, PARAMETER :: ice_aerosol = .TRUE., ice_cloud = .TRUE., ice_precip = .TRUE.
    ! Maximum temperature (K) for homogeneous nucleation
    REAL, PARAMETER :: tmax_homog=243.15

    INTEGER :: ii,jj,kk,nn
    REAL :: pf_imm, pf_dep, pf_hom, jf
    REAL :: Sw_eq, Si, rn, rw, frac, dN

    ! Dust is the only ice nucleation active species
    IF (idu<0) RETURN

    ! Total number of active species
    nn = nspec+1

    ! Ice nucleation modes
    ! 1) Homogeneous freezing: possible for any aqueous droplet with or without insoluble species (DU or BC)
    ! 2) Immersion freezing: possible for aqueous droplets that have an insoluble core (DU or BC)
    ! 3) Deposition freezing: possible for dry insoluble aerosol at sub-saturated conditions (RH < 100%)
    ! 4) Contact freezing: not implemented, because most INPs are immersed in liquid droplets
    ! 5) Condensation freezing: not implemented, because these droplets can freeze during the modelled condensational growth

    DO ii = 1,kbdim
    DO jj = 1,klev
        if (ptemp(ii,jj) > 273.15) cycle

        ! Precipitation
        ! ------------
        ! Always aqueous, so homogeneous and immersion freezing modes are possible
        DO kk = 1,nprc
            IF (pprecp(ii,jj,kk)%numc<prlim .OR. .NOT.ice_precip) CYCLE

            ! Radius of the insoluble portion of the droplet
            rn = MAX(0., (3.*pprecp(ii,jj,kk)%volc(idu)/pprecp(ii,jj,kk)%numc/4./pi)**(1./3.) )
            ! Droplet radius
            rw = (3.*sum(pprecp(ii,jj,kk)%volc(1:nn))/pprecp(ii,jj,kk)%numc/4./pi)**(1./3.)
            ! Equilibrium saturation ratio
            Sw_eq = calc_Sw_eq(pprecp(ii,jj,kk),ptemp(ii,jj))

            ! Immersion freezing
            pf_imm = 0.
            IF (rn>1e-10 .AND. ice_imm) THEN
                jf = calc_Jhet(rn,ptemp(ii,jj),Sw_eq)
                pf_imm = 1. - exp( -jf*ptstep )
            ENDIF

            ! Homogeneous freezing
            pf_hom = 0.
            IF (rw-rn>1e-10 .AND. ptemp(ii,jj)<tmax_homog .AND. ice_hom) THEN
                jf = calc_Jhf(ptemp(ii,jj),Sw_eq)
                pf_hom = 1. - exp( -jf*4./3.*pi*(rw**3-rn**3)*ptstep )
            ENDIF

            frac = MIN(1.,pf_imm+pf_hom-pf_imm*pf_hom)
            dN = pprecp(ii,jj,kk)%numc*frac

            ! Move to the parallel ice bin or to a snow bin
            IF (dN>prlim) CALL rain2ice_snow(kbdim,klev,pprecp,pice,psnow,ii,jj,kk,dN)

        end do


        ! Cloud droplets
        ! --------------
        ! Always aqueous, so homogeneous and immersion freezing modes are possible
        DO kk = 1,ncld
            IF (pcloud(ii,jj,kk)%numc<nlim .OR. .NOT.ice_cloud) CYCLE

            ! Radius of the insoluble portion of the droplet
            rn = MAX(0., (3.*pcloud(ii,jj,kk)%volc(idu)/pcloud(ii,jj,kk)%numc/4./pi)**(1./3.) )
            ! Droplet radius
            rw = (3.*sum(pcloud(ii,jj,kk)%volc(1:nn))/pcloud(ii,jj,kk)%numc/4./pi)**(1./3.)
            ! Equilibrium saturation ratio
            Sw_eq = calc_Sw_eq(pcloud(ii,jj,kk),ptemp(ii,jj))

            ! Immersion freezing
            pf_imm = 0.
            IF (rn>1e-10 .AND. ice_imm) THEN
                jf = calc_Jhet(rn,ptemp(ii,jj),Sw_eq)
                pf_imm = 1. - exp( -jf*ptstep )
            ENDIF

            ! Homogeneous freezing
            pf_hom = 0.
            IF (rw-rn>1e-10 .AND. ptemp(ii,jj)<tmax_homog .AND. ice_hom) THEN
                jf = calc_Jhf(ptemp(ii,jj),Sw_eq)
                pf_hom = 1. - exp( -jf*4./3.*pi*(rw**3-rn**3)*ptstep )
            ENDIF

            frac = MIN(1.,pf_imm+pf_hom-pf_imm*pf_hom)
            dN = pcloud(ii,jj,kk)%numc*frac

            ! Move to the parallel ice bin or to a snow bin
            IF (dN>prlim) CALL cloud2ice_snow(kbdim,klev,pcloud,pice,psnow,ii,jj,kk,dN)

        END DO


        ! Aerosol
        ! -------
        ! Insoluble particles are possible, so all freezing modes are possible
        DO kk = in2a, fn2b
            IF (paero(ii,jj,kk)%numc<nlim .OR. .NOT.ice_aerosol) CYCLE

            ! Radius of the insoluble portion of the droplet
            rn = MAX(0., (3.*paero(ii,jj,kk)%volc(idu)/paero(ii,jj,kk)%numc/4./pi)**(1./3.) )
            ! Droplet radius
            rw = (3.*sum(paero(ii,jj,kk)%volc(1:nn))/paero(ii,jj,kk)%numc/4./pi)**(1./3.)
            ! Equilibrium saturation ratio
            Sw_eq = calc_Sw_eq(paero(ii,jj,kk),ptemp(ii,jj))

            ! Immersion (aqueous aerosol) or deposition (dry insoluble aerosol) freezing
            pf_imm = 0.
            pf_dep = 0.
            IF (rn>1e-10 .AND. rw-rn>1e-10 .AND. ice_imm) THEN
                jf = calc_Jhet(rn,ptemp(ii,jj),Sw_eq)
                pf_imm = 1. - exp( -jf*ptstep )
            ELSEIF (rn>1e-10 .AND. rw-rn<1e-10 .AND. prv(ii,jj)/prs(ii,jj)<1.0 .AND. ice_dep) THEN
                Si=prv(ii,jj)/prsi(ii,jj) ! Water vapor saturation ratio over ice
                jf = calc_Jdep(rn,ptemp(ii,jj),Si)
                pf_dep = 1. - exp( -jf*ptstep )
            ENDIF

            ! Homogeneous freezing
            pf_hom = 0.
            IF (rw-rn>1e-10 .AND. ptemp(ii,jj)<tmax_homog .AND. ice_hom) THEN
                ! Homogeneous freezing
                jf = calc_Jhf(ptemp(ii,jj),Sw_eq)
                pf_hom = 1. - exp( -jf*4./3.*pi*(rw**3-rn**3)*ptstep )
            ENDIF

            frac = MIN(1.,pf_imm+pf_hom+pf_dep-(pf_imm+pf_dep)*pf_hom)
            dN = paero(ii,jj,kk)%numc*frac

            ! Move to the parallel ice bin or to a snow bin
            IF (dN>prlim) CALL aero2ice_snow(kbdim,klev,paero,pice,psnow,ii,jj,kk,dN)

        END DO

    END DO
    END DO

  END SUBROUTINE ice_nucl_driver


  !
  ! Ice nucleation based on Khvorostyanov and Curry, Geophys. Res. Lett., 27, 4081-4084, 2000 [KC00] and
  ! Khvorostyanov and Sassen, Geophys. Res. Lett., 25, 3155-3158, 1998 [KS98]
  !

  REAL FUNCTION calc_Jhet(rn,temp,Sw)
    ! The rate of germ formation (#/s) through heterogeneous freezing following
    ! Khvorostyanov and Curry, Geophys. Res. Lett., 27, 4081-4084, 2000 [KC00]
    ! Additional parameters from
    ! Khvorostyanov and Curry, J. Atmos. Sci., 61, 2676-2691, 2004 [KC04]
    ! Khvorostyanov and Sassen, Geophys. Res. Lett., 25, 3155-3158, 1998 [KS98]
    ! Jeffery and Austin, J. Geophys. Res., 102, 25269-25279, 1997
    ! Li et al., J. Geophys. Res., 118, 11213-11227, 2013
    ! - Here heterogeneos freezing includes just immersion freezing mode
    USE mo_submctl, ONLY : boltz, planck, pi, rg, mwa, avog
    IMPLICIT NONE
    REAL, INTENT(in) :: rn,temp,Sw

    REAL :: Tc, act_energy, Lefm, GG, sigma_is, r_g, sf, crit_energy
    REAL, PARAMETER :: & ! Constants
            C = 1.7e10, & ! Constant (1.7e11 dyn cm^-2 = 1.7e11*1e-5/1e-4 N m^-2 = 1.7e10 N m^-2)
            rho_ice = 900., & ! Density of ice (kg/m^3)
            c_1s = 1e19, & ! The concentration of water molecules adsorbed on 1 m^2 a surface (#/m^2)
            T0 = 273.15 ! 0 C in Kelvins
    REAL, PARAMETER :: & ! Case-dependent parameters
            epsi = 0., & ! Elastic strain produced in ice embryo by the insoluble substrate
            alpha = 0.0, & ! Relative area of active sites
            mis = 0.5 ! Cosine of the contact angle

    calc_Jhet = 0.

    ! Must have a core
    IF (rn<1e-10) RETURN

    Tc = temp-T0 ! Temperature in Celsius

    ! Activation energy (case-dependent)
    !   Eq. 15 in Jeffery and Austin (1997) and parameters for p=1 bar (Table 2) - used in KC04
    act_energy = rg*temp*(347./(temp-177.)-log(4.14/349.))/avog
    !   Khvorostyanov and Sassen (1998) for T < -30 C - used in KC00
    !act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))
    !   Eq. 2 in Li et al. (2013)
    !IF (Tc<=-30.0) THEN
    !    act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))
    !ELSE
    !    act_energy = 5.55*exp(-8.423e-3*Tc+6.384e-4*Tc**2+7.891e-6*Tc**3)/avog*4.1868e3
    !ENDIF

    ! Critical energy of germ formation
    ! a) Ice germ radius (eq. 2.6 in KC00)
    Lefm = (79.7+0.708*Tc-2.5e-3*Tc**2)*4.1868e3 ! Effective latent heat of melting (eq. 6 in KS98)
    GG = rg*temp/mwa/Lefm ! Eq 2.7 in KC00
    IF ( (T0/temp)*Sw**GG<1.0001 ) RETURN
    sigma_is = 28e-3+0.25e-3*Tc ! Surface tension between ice and solution (from KS98)
    r_g = 2.*sigma_is/( rho_ice*Lefm*log((T0/temp)*Sw**GG)-C*epsi**2)
    IF (r_g<=1e-10) RETURN
    ! b) Shape factor (eq. 2.9 in KC00)
    sf=calc_het_sf(rn/r_g,mis)
    ! c) Critical energy (eq. 2.10 in KC00)
    crit_energy = 4.*pi/3.*sigma_is*r_g**2*sf-alpha*(1-mis)*rn**2

    ! Eq 2.1 in KC00
    calc_Jhet= boltz*temp/planck*c_1s*4.*pi*rn**2*exp((-act_energy-crit_energy)/(boltz*temp))

  END FUNCTION calc_Jhet


  REAL FUNCTION calc_Jdep(rn,temp,Si)
    ! The rate of germ formation (#/s) through deposition freezing following
    ! Khvorostyanov and Curry, Geophys. Res. Lett., 27, 4081-4084, 2000 [KC00]
    ! Additional parameters from
    ! Hoose et al., J. Atmos. Sci., 67, 2483-2503, 2010 [Ho10]
    USE mo_submctl, ONLY : boltz, pi, rg, mwa
    IMPLICIT NONE
    REAL, INTENT(in) :: rn,temp,Si

    REAL :: Tc, act_energy, sigma_iv, r_g, sf, crit_energy
    REAL, PARAMETER :: & ! Constants
            C = 1.7e10, & ! Constant (1.7e11 dyn cm^-2 = 1.7e11*1e-5/1e-4 N m^-2 = 1.7e10 N m^-2)
            rho_ice = 900., & ! Density of ice (kg/m^3)
            T0 = 273.15 ! 0 C in Kelvins
    REAL, PARAMETER :: & ! Case-dependent parameters
            epsi = 0., & ! Elastic strain produced in ice embryo by the insoluble substrate
            mis = 0.5 ! Cosine of the contact angle

    calc_Jdep = 0.

    ! Must have a core and supersaturation over ice
    IF (rn<1e-10 .OR. Si<1.0001) RETURN

    Tc = temp-T0 ! Temperature in Celsius

    ! Activation energy (case-dependent)
    !   Set to zero in KC00
    act_energy = 0.

    ! Critical energy of germ formation
    ! a) Ice germ radius (based on eq. 2.12 in KC00)
    sigma_iv = ( (76.1-0.155*Tc) + (28.5+0.25*Tc) )*1e-3 ! Surface tension between ice and vapor (from Ho10)
    r_g = 2.*sigma_iv/( rg*rho_ice/mwa*temp*log(Si)-C*epsi**2)   ! R_v=R*rho_ice/M_ice
    IF (r_g<=1e-10) RETURN
    ! b) Shape factor (eq. 2.9 in KC00)
    sf=calc_het_sf(rn/r_g,mis)
    ! c) Critical energy (eq. 2.12 in KC00)
    crit_energy = 4.*pi/3.*sigma_iv*r_g**2*sf

    ! Eq 2.13 in KC00
    !   The pre-exponential factor (kineticc oefficient) is about (1e26 cm^-2)*rn**2
    calc_Jdep = 1e30*rn**2*exp((-act_energy-crit_energy)/(boltz*temp))

  END FUNCTION calc_Jdep


  REAL FUNCTION calc_Jhf(temp,Sw)
    ! Homogeneous freezing based on Khvorostyanov and Sassen, Geophys. Res. Lett., 25, 3155-3158, 1998 [KS98]
    USE mo_submctl, ONLY : boltz, planck, pi, mwa, rhowa, rg
    IMPLICIT NONE
    REAL, intent(in) :: temp, Sw ! Temperature (K) and water vapor saturation ratio

    REAL :: Tc, act_energy, Lefm, GG, sigma_is, r_g, crit_energy
    REAL, PARAMETER :: & ! Constants
            Nc=5.85e16, & ! The number of water molecules contacting unit area of ice germ (#/m^2) [from KC00]
            rho_ice = 900., & ! Density of ice (kg/m^3)
            T0 = 273.15 ! 0 C in Kelvins

    calc_Jhf = 0.

    Tc = temp-T0 ! Temperature in Celsius

    ! Activation energy (case-dependent)
    !   Khvorostyanov and Sassen (1998) for T < -30 C
    act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))

    ! Critical energy of germ formation
    ! a) Ice germ radius (eq. 9a)
    Lefm = (79.7+0.708*Tc-2.5e-3*Tc**2)*4.1868e3 ! Effective latent heat of melting (eq. 6)
    GG = rg*temp/Lefm/mwa
    IF ( (T0/temp)*Sw**GG<1.0001 ) RETURN
    sigma_is = 28e-3+0.25e-3*Tc ! Surface tension between ice and solution
    r_g = 2.*sigma_is/( rho_ice*Lefm*log((T0/temp)*Sw**GG) )
    IF (r_g<=1e-10) RETURN
    ! c) Critical energy (eq. 9b)
    crit_energy = 4.*pi/3.*sigma_is*r_g**2

    ! Eq. 1
    calc_Jhf = 2.0*Nc*(rhowa*boltz*temp/rho_ice/planck)*sqrt(sigma_is/boltz/temp)*exp((-crit_energy-act_energy)/(boltz*temp))

  END FUNCTION calc_Jhf


  REAL function calc_het_sf(x,mis)
    ! Calculate shape factor for heterogeneous ice nucleation.
    !   Khvorostyanov and Curry, J. Atmos. Sci., 61, 2676-2691, 2004
    REAL :: x ! x=r_core/r_crit
    REAL :: mis ! Cosine of the contact angle
    REAL :: fii, fpsi
    !
    IF (x>100.) THEN
        ! Problems with numerical accuracy when x>>1; x=100 seems to be good limit for using the limiting value
        calc_het_sf = (mis**3-3*mis+2)/4
    ELSE
        fii = sqrt(1.-2.*mis*x+x**2)
        fpsi = (x-mis)/fii
        calc_het_sf = 0.5* (1. + ( ( 1.-mis*x)/fii)**3 + (2.-3.*fpsi+fpsi**3)*x**3 + 3.*mis*(fpsi-1.)*x**2)
    ENDIF
  END function calc_het_sf


  !***********************************************
  !
  ! Prescribed (concentration fixinc>=0.0 #/kg) or parameterized (integer option ice_diag<0) ice
  ! number concentration: if current ice concentartion is lower than the prescribed/parameterized
  ! target ice number concentration, then this target concentartion is reached by converting cloud
  ! droplets to ice or snow.
  !
  !***********************************************

  SUBROUTINE fixed_ice_driver(kbdim, klev,   &
            pcloud, pice, psnow, ptemp, ppres, prv, prs, prsi)

    USE mo_submctl, ONLY : t_section, ncld, nice, nsnw, &
                    rhowa, rda, nlim, prlim, fixinc, fixinc_slope, ice_diag, &
                    fixed_ice_min_Si, fixed_ice_min_rc, fixed_ice_max_T, &
                    ice_source_opt
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,klev
    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), &
                    pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev), &
                    prv(kbdim,klev), prs(kbdim,klev), prsi(kbdim,klev)

    INTEGER :: ii,jj,kk
    REAL :: pdn, S_ice, rc, Ni0, sumICE, dnice, frac

    DO ii = 1,kbdim
    DO jj = 1,klev
        pdn=ppres(ii,jj)/(rda*ptemp(ii,jj)) ! Air density (kg/m^3)

        ! Conditions for ice nucleation
        S_ice = prv(ii,jj)/prsi(ii,jj) ! Saturation with respect to ice
        rc = sum( pcloud(ii,jj,:)%volc(1) )*rhowa/pdn ! Cloud water mixing ratio (kg/kg)
        if ( S_ice < fixed_ice_min_Si .OR. rc < fixed_ice_min_rc .OR. ptemp(ii,jj)>fixed_ice_max_T ) cycle

        IF (ice_diag<0) THEN
            ! Temperature and/or saturation dependent parameterizations for ice concentration (#/m3)
            Ni0 = n_ice_diagnostic(ABS(ice_diag),ptemp(ii,jj),prv(ii,jj),prs(ii,jj),prsi(ii,jj))
        ELSE
            ! Target number concentration of ice, converted to #/m^3
            Ni0 = fixinc*EXP(fixinc_slope*(273.15-ptemp(ii,jj))) * pdn
            Ni0 = fixinc ! COMBLE: concentration in #/m3
        ENDIF

        ! Current ice number concentration (#/m^3)
        sumICE = SUM( pice(ii,jj,:)%numc ) + SUM( psnow(ii,jj,:)%numc )

        if ( sumICE > Ni0 ) cycle

        IF (ice_source_opt>0) THEN
            ! Activate ice starting from the largest cloud bin
            DO kk = nice,1,-1 ! ncld=nice
                IF( Ni0 - sumICE > prlim .AND. pcloud(ii,jj,kk)%numc > nlim) THEN
                    dnice = MAX( 0.0, MIN( Ni0 - sumICE , pcloud(ii,jj,kk)%numc ) )
                    sumICE = sumICE + dnice
                    CALL cloud2ice_snow(kbdim,klev,pcloud,pice,psnow,ii,jj,kk,dnice)
                END IF
            END DO
        ELSEIF (ice_source_opt<0) THEN
            ! Activate ice starting from the smallest cloud bin
            DO kk = 1,nice ! ncld=nice
                IF( Ni0 - sumICE > prlim .AND. pcloud(ii,jj,kk)%numc > nlim) THEN
                    dnice = MAX( 0.0, MIN( Ni0 - sumICE , pcloud(ii,jj,kk)%numc ) )
                    sumICE = sumICE + dnice
                    CALL cloud2ice_snow(kbdim,klev,pcloud,pice,psnow,ii,jj,kk,dnice)
                END IF
            END DO
        ELSE
            ! Activate ice from all cloud bins
            frac = max(0.0, min(1.0, (Ni0-sumICE)/SUM(pcloud(ii,jj,:)%numc)))
            DO kk = 1,nice ! ncld=nice
                IF(pcloud(ii,jj,kk)%numc > nlim .AND. frac*pcloud(ii,jj,kk)%numc > prlim) THEN
                    dnice = frac*pcloud(ii,jj,kk)%numc
                    CALL cloud2ice_snow(kbdim,klev,pcloud,pice,psnow,ii,jj,kk,dnice)
                END IF
            END DO
        ENDIF
    END DO
    END DO

  END SUBROUTINE fixed_ice_driver
  ! ------------------------------------------------------------


  !***********************************************
  !
  ! INP/ice parameterizations based on temperature and saturation ratio
  !
  ! Cooper, W. A.: Ice initiation in natural clouds. Meteor. Monogr., 21, 29-32,.
  ! https://doi.org/10.1175/0065-9401-21.43.29, 1986.
  !
  ! DeMott , P.J., et al.: Predicting global atmospheric ice nuclei distributions and their
  ! impacts on climate, PNAS, 107, 11217-11222, https://doi.org/10.1073/pnas.0910818107, 2010
  !
  ! Fletcher, N.H.: The physics of rainclouds, Cambridge University Press, 1962.
  !
  ! Meyers, M.P.,  DeMott, P.J., and Cotton, W.R.: New primary ice-nucleation parameterizations
  ! in an explicit cloud model,  J. Appl. Meteor., 31, 708-721, 1992.
  !
  ! Murakami, M.: Numerical modeling of dynamical and microphysical evolution of an isolated
  ! convective cloud - The 19 July 1981 CCOPE cloud, J. Met. Soc. Japan. Ser. II, 1990
  !
  ! Reisner, J., Rasmussen, R.M., and Bruintjes, R.T.: Explicit forecasting of supercooled
  ! liquid water in winter storms using the MM5 mesoscale model, Q. J. R. Meteorol. Soc.,
  ! 124, 1071-1107, 1998.
  !
  ! Thompson, G., Rasmussen, R.M., and Manning, K.: Explicit Forecasts of Winter Precipitation
  ! Using an Improved Bulk Microphysics Scheme. Part I: Description and Sensitivity Analysis,
  ! Mon. Wea. Rev., 132, 519-542,
  ! https://doi.org/10.1175/1520-0493(2004)132<0519:EFOWPU>2.0.CO;2, 2004
  !
  REAL FUNCTION n_ice_diagnostic(iopt,ptemp,prv,prs,prsi) RESULT(ni)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: iopt
    REAL, INTENT(in) :: ptemp, prv,prs,prsi

    ! INP parameterizations
    IF (iopt==1) THEN
        ! Fletcher (1962) deposition/condensation freezing parameterization with
        ! temperature limit by Reisner et al. (1998)
        ni = 1.0e-2*EXP(0.6*(273.15-MAX(246.0,ptemp))) ! 1/m3
    ELSEIF (iopt==2) THEN
        ! Meyers et al. (1992) deposition/condensation freezing
        ni = 1.0e3*EXP(-0.639+12.960*MIN(prv/prsi-1.0,0.25)) ! 1/m3
    ELSEIF (iopt==3) THEN
        ! Cooper (1986) parameterization with temperature limit by Thompson et al. (2004)
        ni = 5.0*EXP(0.304*(273.15-MAX(233.0,ptemp))) ! 1/m3
    ELSEIF (iopt==4) THEN
        ! Murakami (1990):  Fletcher (1962) and Reisner et al. (1998) equation with
        ! supersaturation dependency of given by Huffmann and Vali (1973)
        ni = ((prv/prsi-1.0)/(prs/prsi-1.0))**4.5 * &
             1.0e-2*EXP(0.6*(273.15-MAX(246.0,ptemp))) ! 1/m3
    ELSEIF (iopt==5) THEN
        ! Fletcher (1962) and Reisner et al. (1998) equation  + contact freezing from Meyers et al. (1992)
        ni = 1.0e-2*EXP(0.6*(273.15-MAX(246.0,ptemp))) + &
             1.0e3*EXP(-2.8+0.262*(273.15-ptemp))
    ELSEIF (iopt==6) THEN
        ! Meyers et al. (1992) deposition/condensation freezing + contact freezing
        ni = 1.0e3*EXP(-0.639+12.960*MIN(prv/prsi-1.0,0.25)) + &
             1.0e3*EXP(-2.8+0.262*(273.15-ptemp))
    ELSE
        STOP 'INP parameterization not supported!'
    ENDIF

  END FUNCTION n_ice_diagnostic

  !***********************************************



  !***********************************************
  !
  ! Ice formation based on various parameterizations
  !
  SUBROUTINE param_ice_driver(kbdim, klev, paero, pcloud, pprecp, pice, psnow, &
            iopt, ptemp, ppres, prv, prs, prsi)

    USE mo_submctl, ONLY : t_section, nbins, ncld, nprc, nice, nsnw
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim, klev, iopt
    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,nbins),  &
                    pcloud(kbdim,klev,ncld), pprecp(kbdim,klev,nprc), &
                    pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev), &
                    prv(kbdim,klev), prs(kbdim,klev), prsi(kbdim,klev)

    IF (0<iopt .AND. iopt<=10) THEN
        ! 1-10: INAS parameterizations
        CALL ice_inas_driver(kbdim, klev, paero, pcloud, pprecp, pice, psnow, &
                iopt, ptemp, ppres, prv, prsi)
    ELSEIF (10<iopt .AND. iopt<=20) THEN
        ! 11-20: Other parameterizations
        CALL cloud_freeze_driver(kbdim, klev, pcloud, pprecp, pice, psnow, &
                iopt-10, ptemp, ppres, prv, prsi)
    ELSE
        STOP 'Ice parameterization not supported!'
    ENDIF

  END SUBROUTINE param_ice_driver

  ! ------------------------------------------------------------
  !
  ! Ice-nucleation active site (INAS) parameterizations: these give the number of INPs
  ! depending on current temperature, saturation, and aerosol, cloud or rain drop size
  ! distribution. The INP excess will freeze.
  !
  SUBROUTINE ice_inas_driver(kbdim, klev, paero, pcloud, pprecp, pice, psnow, &
            iopt, ptemp, ppres, prv, prsi)

    USE mo_submctl, ONLY : t_section, in2a, nbins, ncld, nprc, nice, nsnw, &
                    rhowa, rda, pi, nlim, prlim, &
                    fixed_ice_min_Si, fixed_ice_min_rc, fixed_ice_max_T
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim, klev, iopt
    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,nbins),  &
                    pcloud(kbdim,klev,ncld), pprecp(kbdim,klev,nprc), &
                    pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev), &
                    prv(kbdim,klev), prsi(kbdim,klev)

    INTEGER :: ii,jj,kk
    REAL :: pdn, S_ice, rc, ns, sumICE, dnice, inp_a(nbins), inp_c(ncld), inp_r(nprc)

    DO ii = 1,kbdim
    DO jj = 1,klev
        pdn=ppres(ii,jj)/(rda*ptemp(ii,jj)) ! Air density (kg/m^3)

        ! Conditions for ice nucleation
        S_ice = prv(ii,jj)/prsi(ii,jj) ! Saturation with respect to ice
        rc = sum( pcloud(ii,jj,:)%volc(1) )*rhowa/pdn ! Cloud water mixing ratio (kg/kg)
        if ( S_ice < fixed_ice_min_Si .OR. rc < fixed_ice_min_rc .OR. ptemp(ii,jj)>fixed_ice_max_T ) cycle

        ! INAS
        IF (iopt==1) THEN
            ! McCluskey, C. S., et al. (2018), Marine and Terrestrial Organic Ice-Nucleating
            ! Particles in Pristine Marine to Continentally Influenced Northeast Atlantic Air
            ! Masses, J. Geophys. Res.-Atmos., 123(11), 6196-6212, doi:10.1029/2017jd028033.
            ns = exp(-0.545*(ptemp(ii,jj)-273.15) + 1.0125) ! m-2
        ELSE
            STOP 'INAS parameterization not supported!'
        ENDIF

        ! Calculate INP concentration based on ...
        ! a) aerosol - currently not included
        inp_a=0.
        DO kk = in2a,nbins
            IF(paero(ii,jj,kk)%numc > nlim) THEN
                ! Here using dry surface area calculated from dry volume:
                !   V=pi/6*D^3=pi/6*(A/pi)^(3/2) => A=pi*(V*6/pi)^(2/3)
                inp_a(kk) = paero(ii,jj,kk)%numc*ns*pi* &
                        (6.0/pi*SUM(paero(ii,jj,kk)%volc(2:))/paero(ii,jj,kk)%numc)**(2./3.)
            ENDIF
        ENDDO
        ! b) cloud droplets
        inp_c=0.
        DO kk = 1,ncld
            IF(pcloud(ii,jj,kk)%numc > nlim) THEN
                ! Dry surface area
                inp_c(kk) = pcloud(ii,jj,kk)%numc*ns*pi* &
                        (6.0/pi*SUM(pcloud(ii,jj,kk)%volc(2:))/pcloud(ii,jj,kk)%numc)**(2./3.)
            ENDIF
        ENDDO
        ! c) rain drops
        inp_r=0.
        DO kk = 1,nprc
            IF(pprecp(ii,jj,kk)%numc > prlim) THEN
                ! Dry surface area
                inp_r(kk) = pprecp(ii,jj,kk)%numc*ns*pi* &
                        (6.0/pi*SUM(pprecp(ii,jj,kk)%volc(2:))/pprecp(ii,jj,kk)%numc)**(2./3.)
            ENDIF
        ENDDO

        ! Calculate the INP excess = the number of new ice particles
        sumICE = SUM(inp_a) + SUM(inp_c) + SUM(inp_r)
        dnice = sumICE - SUM(pice(ii,jj,:)%numc) - SUM(psnow(ii,jj,:)%numc)

        IF (dnice<prlim) CYCLE

        ! New ice based on normalized INP concentrations
        ! a) Aerosol
        inp_a=inp_a/sumICE*dnice
        DO kk = in2a,nbins
            IF(paero(ii,jj,kk)%numc > nlim .AND. inp_a(kk) > prlim) THEN
                dnice = MIN(inp_a(kk),paero(ii,jj,kk)%numc)
                CALL aero2ice_snow(kbdim,klev,paero,pice,psnow,ii,jj,kk,dnice)
            ENDIF
        END DO
        ! b) Cloud
        inp_c=inp_c/sumICE*dnice
        DO kk = 1,ncld
            IF(pcloud(ii,jj,kk)%numc > nlim .AND. inp_c(kk) > prlim) THEN
                dnice = MIN(inp_c(kk),pcloud(ii,jj,kk)%numc)
                CALL cloud2ice_snow(kbdim,klev,pcloud,pice,psnow,ii,jj,kk,dnice)
            ENDIF
        END DO
        ! c) Rain
        inp_r=inp_r/sumICE*dnice
        DO kk = 1,nprc
            IF(pprecp(ii,jj,kk)%numc > nlim .AND. inp_r(kk) > prlim) THEN
                dnice = MIN(inp_r(kk),pprecp(ii,jj,kk)%numc)
                CALL rain2ice_snow(kbdim,klev,pprecp,pice,psnow,ii,jj,kk,dnice)
            ENDIF
        END DO
    END DO
    END DO

  END SUBROUTINE ice_inas_driver


  ! ------------------------------------------------------------
  !
  ! Aerosol/cloud/rain freezing parameterizations.
  !
  SUBROUTINE cloud_freeze_driver(kbdim, klev, pcloud, pprecp, pice, psnow, &
             iopt, ptemp, ppres, prv, prsi)

    USE mo_submctl, ONLY : t_section, ncld, nprc, nice, nsnw, &
                    rhowa, rda, nlim, prlim, &
                    fixed_ice_min_Si, fixed_ice_min_rc, fixed_ice_max_T
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim, klev, iopt
    TYPE(t_section), INTENT(inout) ::  &
                    pcloud(kbdim,klev,ncld), pprecp(kbdim,klev,nprc), &
                    pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev), &
                    prv(kbdim,klev), prsi(kbdim,klev)

    INTEGER :: ii,jj,kk
    REAL :: pdn, S_ice, rc, ns, dnice

    DO ii = 1,kbdim
    DO jj = 1,klev
        pdn=ppres(ii,jj)/(rda*ptemp(ii,jj)) ! Air density (kg/m^3)

        ! Conditions for ice nucleation
        S_ice = prv(ii,jj)/prsi(ii,jj) ! Saturation with respect to ice
        rc = sum( pcloud(ii,jj,:)%volc(1) )*rhowa/pdn ! Cloud water mixing ratio (kg/kg)
        if ( S_ice < fixed_ice_min_Si .OR. rc < fixed_ice_min_rc .OR. ptemp(ii,jj)>fixed_ice_max_T ) cycle

        IF (iopt==1) THEN
            ! Droplet freezing parameterization by Bigg (1953)
            !
            ! Thompson, G., P. R. Field, R. M. Rasmussen, and W. D. Hall, 2008:
            ! Explicit Forecasts of Winter Precipitation Using an Improved Bulk Microphysics Scheme.
            ! Part II: Implementation of a New Snow Parameterization. Mon. Wea. Rev., 136, 5095-5115,
            ! https://doi.org/10.1175/2008MWR2387.1.
            !
            ! Probability of freezing for water drops with volume vol:
            !   P = 1 - exp(ns*vol)
            ! where constant ns is
            ns = -120.0*5.2e-4*( EXP(ptemp(ii,jj)-273.15)-1.0 )
            !
            ! Cloud droplet freezing
            DO kk = 1,ncld
                IF(pcloud(ii,jj,kk)%numc > nlim) THEN
                    dnice = pcloud(ii,jj,kk)%numc * &
                            (1.0-EXP(ns*pcloud(ii,jj,kk)%volc(1)/pcloud(ii,jj,kk)%numc))
                    !
                    ! Apply to bin kk
                    IF (dnice>prlim) CALL cloud2ice_snow(kbdim,klev,pcloud,pice,psnow,ii,jj,kk,dnice)
                ENDIF
            ENDDO
            !
            ! Rain drop freezing
            DO kk = 1,nprc
                IF(pprecp(ii,jj,kk)%numc > prlim) THEN
                    dnice = pprecp(ii,jj,kk)%numc * &
                            (1.0-EXP(ns*pprecp(ii,jj,kk)%volc(1)/pprecp(ii,jj,kk)%numc))
                    !
                    ! Apply to bin kk
                    IF (dnice>prlim) CALL rain2ice_snow(kbdim,klev,pprecp,pice,psnow,ii,jj,kk,dnice)
                ENDIF
            ENDDO
        ELSE
            STOP 'Freezing parameterization not supported!'
        ENDIF

    END DO
    END DO

  END SUBROUTINE cloud_freeze_driver
  !***********************************************


  ! Move frozen aerosol to ice or snow bins
  SUBROUTINE aero2ice_snow(kbdim,klev,paero,pice,psnow,ii,jj,kk,dN)
    USE mo_submctl, ONLY : t_section, nbins, nice, nsnw, &
            fn1a, nspec, ice_target_opt
    IMPLICIT NONE
    ! Inputs/outputs
    INTEGER, INTENT(in) :: kbdim,klev
    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,nbins), &
                                      pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    INTEGER, INTENT(in) :: ii, jj, kk ! Aerosol bin
    REAL, INTENT(inout) :: dN ! Change in number concentration
    ! Local
    REAL :: frac, vol
    INTEGER :: nn, ss
    !
    nn = nspec + 1 ! Aerosol species + water
    !
    dN = MIN(dN,paero(ii,jj,kk)%numc)
    frac = dN/paero(ii,jj,kk)%numc
    !
    IF (ice_target_opt<0) THEN
        ! Add to the parallel ice bin
        ss = MAX(1,kk-fn1a) ! subtract 1a aerosol bins
        pice(ii,jj,ss)%volc(1:nn) = pice(ii,jj,ss)%volc(1:nn) + paero(ii,jj,kk)%volc(1:nn)*frac
        pice(ii,jj,ss)%numc   = pice(ii,jj,ss)%numc + dN
    ELSEIF (ice_target_opt==0) THEN
        ! Add to the matching snow bin - based on wet volume
        vol=SUM(paero(ii,jj,kk)%volc(1:nn))/paero(ii,jj,kk)%numc
        ss =MAX(1,COUNT(vol>psnow(ii,jj,:)%vlolim))

        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + paero(ii,jj,kk)%volc(1:nn)*frac
        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dN
    ELSE
        ! Add to the ss:th snow bin
        ss=MIN(nsnw,ice_target_opt)
        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + paero(ii,jj,kk)%volc(1:nn)*frac
        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dN
    ENDIF
    !
    paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc - dN
    paero(ii,jj,kk)%volc(1:nn) = paero(ii,jj,kk)%volc(1:nn) - paero(ii,jj,kk)%volc(1:nn)*frac
    !
  END SUBROUTINE aero2ice_snow


  ! Move frozen cloud droplets to ice or snow bins
  SUBROUTINE cloud2ice_snow(kbdim,klev,pcloud,pice,psnow,ii,jj,kk,dN)
    USE mo_submctl, ONLY : t_section, ncld, nice, nsnw, &
            nspec, ice_target_opt
    IMPLICIT NONE
    ! Inputs/outputs
    INTEGER, INTENT(in) :: kbdim,klev
    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), &
                                      pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    INTEGER, INTENT(in) :: ii, jj, kk ! Cloud bin
    REAL, INTENT(inout) :: dN ! Change in number concentration
    ! Local
    REAL :: frac, vol
    INTEGER :: nn, ss
    !
    nn = nspec + 1 ! Aerosol species + water
    !
    dN = MIN(dN,pcloud(ii,jj,kk)%numc)
    frac = dN/pcloud(ii,jj,kk)%numc
    !
    IF (ice_target_opt<0) THEN
        ! Add to the parallel ice bin
        pice(ii,jj,kk)%volc(1:nn) = pice(ii,jj,kk)%volc(1:nn) + pcloud(ii,jj,kk)%volc(1:nn)*frac
        pice(ii,jj,kk)%numc   = pice(ii,jj,kk)%numc + dN
    ELSEIF (ice_target_opt==0) THEN
        ! Add to the matching snow bin - based on wet volume
        vol=SUM(pcloud(ii,jj,kk)%volc(1:nn))/pcloud(ii,jj,kk)%numc
        ss =MAX(1,COUNT(vol>psnow(ii,jj,:)%vlolim))
        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + pcloud(ii,jj,kk)%volc(1:nn)*frac
        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dN
    ELSE
        ! Add to the ss:th snow bin
        ss=MIN(nsnw,ice_target_opt)
        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + pcloud(ii,jj,kk)%volc(1:nn)*frac
        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dN
    ENDIF
    !
    pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc - dN
    pcloud(ii,jj,kk)%volc(1:nn) = pcloud(ii,jj,kk)%volc(1:nn) - pcloud(ii,jj,kk)%volc(1:nn)*frac
    !
  END SUBROUTINE cloud2ice_snow


  ! Move frozen rain droplets to ice or snow bins
  SUBROUTINE rain2ice_snow(kbdim,klev,pprecp,pice,psnow,ii,jj,kk,dN)
    USE mo_submctl, ONLY : t_section, nprc, nice, nsnw, &
            fnp2a, nspec, prlim, ice_target_opt, calc_correlation
    IMPLICIT NONE
    ! Inputs/outputs
    INTEGER, INTENT(in) :: kbdim,klev
    TYPE(t_section), INTENT(inout) :: pprecp(kbdim,klev,nprc), &
                                      pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    INTEGER, INTENT(in) :: ii, jj, kk ! Rain bin
    REAL, INTENT(inout) :: dN ! Change in number concentration
    ! Local
    REAL :: frac, vol, ra, rb
    INTEGER :: nn, ss
    !
    nn = nspec + 1 ! Aerosol species + water
    !
    dN = MIN(dN,pprecp(ii,jj,kk)%numc)
    frac = dN/pprecp(ii,jj,kk)%numc
    !
    IF (ice_target_opt<0) THEN
        ! Add to the matching ice bin - based on dry volume
        ! 1) Matching a-bin
        vol=SUM(pprecp(ii,jj,kk)%volc(2:nn))/pprecp(ii,jj,kk)%numc
        ss=MAX(1,COUNT(vol>pice(ii,jj,1:fnp2a)%vlolim)) ! Ice a-bin
        ! 2) Select a or b bin
        IF (pice(ii,jj,ss+fnp2a)%numc<=prlim) THEN
            ! Empty b bin so select a
            !ss = ss
        ELSEIF (pice(ii,jj,ss)%numc<=prlim) THEN
            ! Empty a bin so select b
            ss = ss + fnp2a
        ELSE
            ! Both are present - find bin based on compositional similarity
            ra = calc_correlation(pice(ii,jj,ss)%volc(2:nn),pprecp(ii,jj,kk)%volc(2:nn),nspec)
            rb = calc_correlation(pice(ii,jj,ss+fnp2a)%volc(2:nn),pprecp(ii,jj,kk)%volc(2:nn),nspec)
            IF (ra<rb) ss = ss + fnp2a
        ENDIF
        pice(ii,jj,ss)%volc(1:nn) = pice(ii,jj,ss)%volc(1:nn) + pprecp(ii,jj,kk)%volc(1:nn)*frac
        pice(ii,jj,ss)%numc   = pice(ii,jj,ss)%numc + dN
    ELSEIF (ice_target_opt==0) THEN
        ! Add to the matching snow bin - based on wet volume
        vol=SUM(pprecp(ii,jj,kk)%volc(1:nn))/pprecp(ii,jj,kk)%numc
        ss=MAX(1,COUNT(vol>psnow(ii,jj,:)%vlolim))
        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + pprecp(ii,jj,kk)%volc(1:nn)*frac
        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dN
    ELSE
        ! Add to the ss:th snow bin
        ss=MIN(nsnw,ice_target_opt)
        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + pprecp(ii,jj,kk)%volc(1:nn)*frac
        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dN
    ENDIF
    !
    pprecp(ii,jj,kk)%numc = pprecp(ii,jj,kk)%numc - dN
    pprecp(ii,jj,kk)%volc(1:nn) = pprecp(ii,jj,kk)%volc(1:nn) - pprecp(ii,jj,kk)%volc(1:nn)*frac
    !
  END SUBROUTINE rain2ice_snow


  SUBROUTINE ice_melt(kbdim,klev,   &
                      pcloud,pice,pprecp,psnow, &
                      ptemp )

    USE mo_submctl, ONLY : t_section,   &
                               ncld,        &
                               nice,        &
                               nsnw,        &
                               nprc,        &
                               nspec,       &
                               prlim

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,klev
    REAL, INTENT(in) :: ptemp(kbdim,klev)

    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), &
                                      pice(kbdim,klev,nice),   &
                                      psnow(kbdim,klev,nsnw),  &
                                      pprecp(kbdim,klev,nprc)

    INTEGER :: ii,jj,kk,ss

    DO ii = 1,kbdim
       DO jj = 1,klev
          ! Ice and snow melt when temperature above 273.15 K
          ! => should add the effect of freezing point depression
          if (ptemp(ii,jj) <= 273.15 ) cycle

          DO kk = 1,nice
              ! Ice => cloud water (parallel bin)
              IF (pice(ii,jj,kk)%numc<prlim) CYCLE
              DO ss = 1,nspec+1
                  pcloud(ii,jj,kk)%volc(ss) = pcloud(ii,jj,kk)%volc(ss) + pice(ii,jj,kk)%volc(ss)
                  pice(ii,jj,kk)%volc(ss) = 0.
              END DO
              pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc + pice(ii,jj,kk)%numc
              pice(ii,jj,kk)%numc = 0.
          END DO

          DO kk =1,nsnw
              ! Snow => precipitation (bin 1)
              IF (psnow(ii,jj,kk)%numc<prlim) CYCLE
              DO ss = 1,nspec+1
                  pprecp(ii,jj,1)%volc(ss) = pprecp(ii,jj,1)%volc(ss) + psnow(ii,jj,kk)%volc(ss)
                  psnow(ii,jj,kk)%volc(ss) = 0.
              END DO
              pprecp(ii,jj,1)%numc = pprecp(ii,jj,1)%numc + psnow(ii,jj,kk)%numc
              psnow(ii,jj,kk)%numc = 0.
            END DO
       END DO
    END DO

  END SUBROUTINE ice_melt


  SUBROUTINE autosnow(kbdim,klev,   &
                      pice,psnow)
  !
  ! Uses a more straightforward method for converting ice to snow.
  ! Assume a lognormal ice distribution for each bin. Sigma_g is an adjustable
  ! parameter and is set to 1.2 by default
  !
    
    USE mo_submctl, ONLY : t_section,   &
                               nice,        &
                               nsnw,        &
                               nspec,       &
                               pi6,         &
                               prlim,       &
                               autoc_snow_zd0, autoc_snow_sigmag
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,klev
    TYPE(t_section), INTENT(inout) :: pice(kbdim,klev,nice)
    TYPE(t_section), INTENT(inout) :: psnow(kbdim,klev,nsnw)

    REAL :: Vrem, Nrem, Vtot, Ntot
    REAL :: dvg,dg
    INTEGER :: ii,jj,cc,ss

    DO jj = 1,klev
       DO ii = 1,kbdim
          DO cc = 1,nice

             Ntot = pice(ii,jj,cc)%numc
             Vtot = SUM(pice(ii,jj,cc)%volc(:))

             IF ( Ntot > prlim .AND. Vtot > 0. ) THEN

                ! Volume geometric mean diameter
                dvg = ((Vtot/Ntot/pi6)**(1./3.))*EXP( (3.*LOG(autoc_snow_sigmag)**2)/2. )
                dg = dvg*EXP( -3.*LOG(autoc_snow_sigmag)**2 )

                Vrem = Vtot*( 1. - cumlognorm(dvg,autoc_snow_sigmag,autoc_snow_zd0) )
                Nrem = Ntot*( 1. - cumlognorm(dg,autoc_snow_sigmag,autoc_snow_zd0) )

                IF ( Vrem > 0. .AND. Nrem > prlim) THEN

                   ! Put the mass and number to the first snow bin and remove from ice
                   DO ss = 1,nspec+1
                      psnow(ii,jj,1)%volc(ss) = psnow(ii,jj,1)%volc(ss) + pice(ii,jj,cc)%volc(ss)*(Nrem/Ntot)
                      pice(ii,jj,cc)%volc(ss) = pice(ii,jj,cc)%volc(ss)*(1. - (Nrem/Ntot))
                   END DO

                   psnow(ii,jj,1)%numc = psnow(ii,jj,1)%numc + Nrem
                   pice(ii,jj,cc)%numc = pice(ii,jj,cc)%numc - Nrem

                END IF ! Nrem Vrem

             END IF ! Ntot Vtot

          END DO ! cc
       END DO ! ii
    END DO ! jj

  END SUBROUTINE autosnow

  !
  ! -----------------------------------------------------------------
  !
  REAL FUNCTION cumlognorm(dg,sigmag,dpart)
    
    IMPLICIT NONE
    ! Cumulative lognormal function
    REAL, INTENT(in) :: dg
    REAL, INTENT(in) :: sigmag
    REAL, INTENT(in) :: dpart

    REAL :: hlp1,hlp2

    hlp1 = ( LOG(dpart) - LOG(dg) )
    hlp2 = SQRT(2.)*LOG(sigmag)
    cumlognorm = 0.5 + 0.5*ERF( hlp1/hlp2 )

  END FUNCTION cumlognorm


  !***********************************************
  ! Secondary ice production
  !***********************************************

  ! Hallett-Mossop or splintering during riming
  !   Cotton, W. R., Tripoli,  G. J., Rauber, R. M., and Mulvihill, E. A.: Numerical Simulation of
  !   the Effects of Varying Ice Crystal Nucleation Rates and Aggregation Processes on Orographic
  !   Snowfall, J. Appl. Meteor. Climatol., 25, 1658-1680,
  !   https://doi.org/10.1175/1520-0450(1986)025<1658:NSOTEO>2.0.CO;2, 1986.
  !
  !   Sotiropoulou, G., Vignon, E., Young, G., Morrison, H., O'Shea, S. J., Lachlan-Cope, T.,
  !   Berne, A., and Nenes, A.: Secondary ice production in summer clouds over the Antarctic coast:
  !   an underappreciated process in atmospheric models, Atmos. Chem. Phys., 21, 755-771,
  !   https://doi.org/10.5194/acp-21-755-2021, 2021.
  SUBROUTINE sip_hm(kbdim,klev,pice,psnow,ptemp)
    USE mo_submctl, ONLY : t_section, nice, nsnw, fnp2a, prlim, rhowa, &
        hm_c_mult, hm_frag_vfrac, & ! Parameters
        rime_volc_ice, rime_volc_snw ! Accumulated rime (rime water volume concentration, m3/m3)
    IMPLICIT NONE
    ! Inputs/outputs
    INTEGER, INTENT(in) :: kbdim,klev
    TYPE(t_section), INTENT(inout) :: pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    REAL, INTENT(in) :: ptemp(kbdim,klev)
    ! Constants for Hallett-Mossop (following SB microphysics)
    real, parameter :: t_mult_min = 265.0  !..min temp. splintering
    real, parameter :: t_mult_max = 270.0  !..max temp. splintering
    real, parameter :: t_mult_opt = 268.0  !..opt temp. splintering
    ! Local parameters
    INTEGER :: ii, jj, cc, bb
    REAL :: fact, dN, vol
    !
    DO jj = 1,klev
    DO ii = 1,kbdim
        IF (t_mult_min<ptemp(ii,jj) .AND. ptemp(ii,jj)<t_mult_max) THEN
            ! The number of splinters depends on temperature
            IF (ptemp(ii,jj) > t_mult_opt) THEN
                fact = hm_c_mult * (t_mult_max - ptemp(ii,jj))/(t_mult_max-t_mult_opt)
            ELSE
                fact = hm_c_mult * (ptemp(ii,jj) - t_mult_min)/(t_mult_opt-t_mult_min)
            ENDIF
            !
            ! Ice collecting rime
            DO cc = 1,nice
                ! Splinters
                dN=fact*rime_volc_ice(ii,jj,cc)*rhowa
                IF (dN>prlim .AND. pice(ii,jj,cc)%numc>prlim) THEN
                    ! Fragment dry volume
                    vol=hm_frag_vfrac*SUM(pice(ii,jj,cc)%volc(2:))/pice(ii,jj,cc)%numc
                    ! Ice a-bin
                    bb=MAX(1,COUNT(vol>pice(ii,jj,1:fnp2a)%vlolim))
                    IF (cc>fnp2a) bb=bb+fnp2a
                    !
                    CALL ice2ice(ii,jj,nice,pice,dN,cc,bb)
                ENDIF
            ENDDO
            !
            ! Snow collecting rime
            DO cc = 1,nsnw
                ! Splinters
                dN=fact*rime_volc_snw(ii,jj,cc)*rhowa
                IF (dN>prlim .AND. psnow(ii,jj,cc)%numc>prlim) THEN
                    ! Fragment wet volume
                    vol=hm_frag_vfrac*SUM(psnow(ii,jj,cc)%volc(:))/psnow(ii,jj,cc)%numc
                    ! Snow bin
                    bb=MAX(1,COUNT(vol>psnow(ii,jj,:)%vlolim))
                    !
                    CALL snow2snow(ii,jj,nsnw,psnow,dN,cc,bb)
                ENDIF
            ENDDO
        ENDIF
    END DO
    END DO
  END SUBROUTINE sip_hm

  ! Droplet fragmentation during freezing
  !   Sullivan, S. C., Hoose, C., Kiselev, A., Leisner, T., and Nenes, A.: Initiation of secondary
  !   ice production in clouds, Atmos. Chem. Phys., 18, 1593-1610,
  !   https://doi.org/10.5194/acp-18-1593-2018, 2018.
  SUBROUTINE sip_df(kbdim,klev,pcloud,pprecp,pice,psnow,ptemp)
    USE mo_submctl, ONLY : t_section, ncld, nprc, nice, nsnw, fnp2a, prlim, nlim, pi, &
        df_c_mult, df_tmin, df_tmax, df_frag_vfrac, & ! Parameters
        coll_rate_ic, coll_rate_ir, coll_rate_sc, coll_rate_sr ! Accumulated collisions (#/m3)
    IMPLICIT NONE
    ! Inputs/outputs
    INTEGER, INTENT(in) :: kbdim,klev
    TYPE(t_section), INTENT(inout) :: pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    TYPE(t_section), INTENT(in) ::  pcloud(kbdim,klev,ncld), pprecp(kbdim,klev,nprc)
    REAL, INTENT(in) :: ptemp(kbdim,klev)
    ! Local parameters
    INTEGER :: ii, jj, cc, bb, aa
    REAL :: fact, dN, vol
    !
    DO jj = 1,klev
    DO ii = 1,kbdim
        IF (df_tmin<ptemp(ii,jj) .AND. ptemp(ii,jj)<df_tmax) THEN
            ! The number of fragments depends on temperature: 0.2*f(T;m=258 K,s=10 K)
            fact=0.2/(10.*sqrt(2.*pi))*exp(-0.5*((ptemp(ii,jj)-258.)/10.)**2)
            !
            ! Ice-cloud/rain collisions producing ice
            DO cc = 1,nice
                ! New particles are taken from snow bin cc
                IF (pice(ii,jj,cc)%numc<prlim) CYCLE
                !
                ! Ice-cloud
                DO aa = 1,ncld
                    ! New particles
                    dN=df_c_mult*(pcloud(ii,jj,aa)%dwet**4)*fact*coll_rate_ic(ii,jj,cc,aa)
                    IF (dN<prlim .OR. pcloud(ii,jj,aa)%numc<nlim) CYCLE
                    !
                    ! Fragment dry volume
                    vol=df_frag_vfrac*SUM(pcloud(ii,jj,aa)%volc(2:))/pcloud(ii,jj,aa)%numc
                    ! Ice a-bin
                    bb=MAX(1,COUNT(vol>pice(ii,jj,1:fnp2a)%vlolim))
                    IF (cc>fnp2a) bb=bb+fnp2a
                    !
                    CALL ice2ice(ii,jj,nice,pice,dN,cc,bb)
                ENDDO
                !
                ! Ice-rain
                DO aa = 1,nprc
                    ! New particles
                    dN=df_c_mult*(pprecp(ii,jj,aa)%dwet**4)*fact*coll_rate_ir(ii,jj,cc,aa)
                    IF (dN<prlim .OR. pprecp(ii,jj,aa)%numc<prlim) CYCLE
                    !
                    ! Fragment dry volume
                    vol=df_frag_vfrac*SUM(pprecp(ii,jj,aa)%volc(2:))/pprecp(ii,jj,aa)%numc
                    ! Ice a-bin
                    bb=MAX(1,COUNT(vol>pice(ii,jj,1:fnp2a)%vlolim))
                    IF (cc>fnp2a) bb=bb+fnp2a
                    !
                    CALL ice2ice(ii,jj,nice,pice,dN,cc,bb)
                ENDDO
            ENDDO
            !
            ! Snow-cloud/rain collisions producing snow
            DO cc = 1,nsnw
                ! New particles are taken from snow bin cc
                IF (psnow(ii,jj,cc)%numc<prlim) CYCLE
                !
                ! Snow-cloud
                DO aa = 1,ncld
                    ! New particles
                    dN=df_c_mult*(pcloud(ii,jj,aa)%dwet**4)*fact*coll_rate_sc(ii,jj,cc,aa)
                    IF (dN<prlim .OR. pcloud(ii,jj,aa)%numc<nlim) CYCLE
                    !
                    ! Fragment wet volume
                    vol=df_frag_vfrac*SUM(pcloud(ii,jj,aa)%volc(:))/pcloud(ii,jj,aa)%numc
                    ! Snow bin
                    bb=MAX(1,COUNT(vol>psnow(ii,jj,:)%vlolim))
                    !
                    CALL snow2snow(ii,jj,nsnw,psnow,dN,cc,bb)
                ENDDO
                !
                ! Snow-rain
                DO aa = 1,nprc
                    ! New particles
                    dN=df_c_mult*(pprecp(ii,jj,aa)%dwet**4)*fact*coll_rate_sr(ii,jj,cc,aa)
                    IF (dN<prlim .OR. pprecp(ii,jj,aa)%numc<prlim) CYCLE
                    !
                    ! Fragment wet volume
                    vol=df_frag_vfrac*SUM(pprecp(ii,jj,aa)%volc(:))/pprecp(ii,jj,aa)%numc
                    ! Snow bin
                    bb=MAX(1,COUNT(vol>psnow(ii,jj,:)%vlolim))
                    !
                    CALL snow2snow(ii,jj,nsnw,psnow,dN,cc,bb)
                ENDDO
            ENDDO
        ENDIF
    END DO
    END DO
  END SUBROUTINE sip_df

  ! Ice-ice collisional breakup
  !   Sullivan, S. C., Hoose, C., Kiselev, A., Leisner, T., and Nenes, A.: Initiation of secondary
  !   ice production in clouds, Atmos. Chem. Phys., 18, 1593-1610,
  !   https://doi.org/10.5194/acp-18-1593-2018, 2018.
  !
  !   Sullivan, S. C., Barthlott, C., Crosier, J., Zhukov, I., Nenes, A., and Hoose, C.: The effect
  !   of secondary ice production parameterization on the simulation of a cold frontal rainband,
  !   Atmos. Chem. Phys., 18, 16461-16480, https://doi.org/10.5194/acp-18-16461-2018, 2018.
  !
  !   Sotiropoulou, G., Sullivan, S., Savre, J., Lloyd, G., Lachlan-Cope, T., Ekman, A. M. L., and
  !   Nenes, A.: The impact of secondary ice production on Arctic stratocumulus, Atmos. Chem.
  !   Phys., 20, 1301-1316, https://doi.org/10.5194/acp-20-1301-2020, 2020.
  SUBROUTINE sip_iibr(kbdim,klev,pice,psnow,ptemp)
    USE mo_submctl, ONLY : t_section, nice, nsnw, fnp2a, prlim, &
        iibr_fbr, iibr_tmin, iibr_tmax, iibr_dref, iibr_frag_vfrac, & ! Parameters
        coll_rate_ii, coll_rate_si, coll_rate_ss ! Accumulated collisions (#/m3)
    IMPLICIT NONE
    ! Inputs/outputs
    INTEGER, INTENT(in) :: kbdim,klev
    TYPE(t_section), INTENT(inout) :: pice(kbdim,klev,nice), psnow(kbdim,klev,nsnw)
    REAL, INTENT(in) :: ptemp(kbdim,klev)
    ! Local parameters
    INTEGER :: ii, jj, cc, bb, aa
    REAL :: fact, scaler, dN, vol
    !
    ! Size-scaling from Sotiropuolou et al. (2021): d/d0, where d0=0.02 m is the size of hail
    ! balls in the experiments and d is the size of the ice particle that undergoes fracturing.
    ! Here d is the the size of the collecting particle and d0=iibr_dref is adjustable parameter.
    ! Zero or negative iibr_dref means that this scaling is not used.
    scaler=1.
    !
    DO jj = 1,klev
    DO ii = 1,kbdim
        IF (iibr_tmin<ptemp(ii,jj) .AND. ptemp(ii,jj)<iibr_tmax .AND. iibr_fbr>0.) THEN
            ! The number of splinters depends on temperature
            fact= iibr_fbr*(ptemp(ii,jj)-iibr_tmin)**1.2*exp((iibr_tmin-ptemp(ii,jj))*0.2)
            !
            ! Ice-ice (smaller and equal) collisions producing ice
            DO cc = 1,nice
                ! New particles are taken from ice bin cc
                IF (pice(ii,jj,cc)%numc<prlim) CYCLE
                !
                ! Size scaling based on collecting ice particle size
                IF (iibr_dref>0.) scaler=pice(ii,jj,cc)%dwet/iibr_dref
                !
                DO aa = 1,nice
                    ! New particles
                    dN=fact*coll_rate_ii(ii,jj,cc,aa)*scaler
                    IF (dN<prlim) CYCLE
                    !
                    ! Fragment dry volume
                    vol=iibr_frag_vfrac*SUM(pice(ii,jj,aa)%volc(2:))/pice(ii,jj,aa)%numc
                    ! Ice a-bin
                    bb=MAX(1,COUNT(vol>pice(ii,jj,1:fnp2a)%vlolim))
                    IF (aa>fnp2a) bb=bb+fnp2a
                    !
                    CALL ice2ice(ii,jj,nice,pice,dN,cc,aa)
                ENDDO
            ENDDO
            !
            ! Snow-ice/snow collisions
            DO cc = 1,nsnw
                ! New particles are taken from snow bin cc
                IF (psnow(ii,jj,cc)%numc<prlim) CYCLE
                !
                ! Size scaling based on collecting snow size
                IF (iibr_dref>0.) scaler=psnow(ii,jj,cc)%dwet/iibr_dref
                !
                ! Snow-ice producing ice
                DO aa = 1,nice
                    ! New particles
                    dN=fact*coll_rate_si(ii,jj,cc,aa)*scaler
                    IF (dN<prlim) CYCLE
                    !
                    ! Fragment dry volume
                    vol=iibr_frag_vfrac*SUM(pice(ii,jj,aa)%volc(2:))/pice(ii,jj,aa)%numc
                    ! Ice a-bin
                    bb=MAX(1,COUNT(vol>pice(ii,jj,1:fnp2a)%vlolim))
                    IF (aa>fnp2a) bb=bb+fnp2a
                    !
                    CALL snow2ice(ii,jj,nice,nsnw,pice,psnow,dN,cc,aa)
                ENDDO
                !
                ! Snow-snow producing snow
                DO aa = 1,cc
                    ! New particles
                    dN=fact*coll_rate_ss(ii,jj,cc,aa)*scaler
                    IF (dN<prlim) CYCLE
                    !
                    ! Fragment wet volume
                    vol=iibr_frag_vfrac*SUM(psnow(ii,jj,aa)%volc(:))/psnow(ii,jj,aa)%numc
                    ! Snow bin
                    bb=MAX(1,COUNT(vol>psnow(ii,jj,:)%vlolim))
                    !
                    CALL snow2snow(ii,jj,nsnw,psnow,dN,cc,aa)
                ENDDO
            ENDDO
        ENDIF
    END DO
    END DO
  END SUBROUTINE sip_iibr

  SUBROUTINE ice2ice(ii,jj,nice,pice,dN,cc,bb)
    USE mo_submctl, ONLY : t_section
    ! Inputs/outputs
    INTEGER, INTENT(in) :: ii, jj, nice ! Dimensions
    INTEGER, INTENT(in) :: cc, bb ! Bin indices for the source (cc) and target (bb) ice
    TYPE(t_section), INTENT(inout) :: pice(ii,jj,nice)
    REAL, INTENT(inout) :: dN ! Change in number concentration
    ! Local
    REAL :: vfrac
    !
    ! Limit dN to 10% of the parent number to avoid large changes in mean dry size
    dN=MIN(dN,0.1*pice(ii,jj,cc)%numc)
    !
    IF (bb==cc) THEN
        ! Just increase the number concentration
        pice(ii,jj,bb)%numc = pice(ii,jj,bb)%numc + dN
    ELSE
        ! Volume fraction to be removed (based on dry size)
        vfrac = dN/pice(ii,jj,cc)%numc*(pice(ii,jj,bb)%Vmid/pice(ii,jj,cc)%Vmid)
        ! Move dN splinters from ice bin cc to ice bin bb
        pice(ii,jj,bb)%numc = pice(ii,jj,bb)%numc + dN
        pice(ii,jj,bb)%volc(:) = pice(ii,jj,bb)%volc(:) + vfrac*pice(ii,jj,cc)%volc(:)
        pice(ii,jj,cc)%volc(:) = (1.-vfrac)*pice(ii,jj,cc)%volc(:)
    ENDIF
    !
  END SUBROUTINE ice2ice

  SUBROUTINE snow2snow(ii,jj,nsnw,psnow,dN,cc,bb)
    USE mo_submctl, ONLY : t_section
    ! Inputs/outputs
    INTEGER, INTENT(in) :: ii, jj, nsnw ! Dimensions
    INTEGER, INTENT(in) :: cc, bb ! Bin indices for the source (cc) and target (bb) snow
    TYPE(t_section), INTENT(inout) :: psnow(ii,jj,nsnw)
    REAL, INTENT(inout) :: dN ! Change in number concentration
    ! Local
    REAL :: vfrac
    !
    ! Limit dN to 10% of the parent number to avoid large changes in mean dry size
    dN=MIN(dN,0.1*psnow(ii,jj,cc)%numc)
    !
    IF (bb==cc) THEN
        ! Just increase the number concentration
        psnow(ii,jj,bb)%numc = psnow(ii,jj,bb)%numc + dN
    ELSE
        ! Volume fraction to be removed (based on wet size)
        vfrac = dN/psnow(ii,jj,cc)%numc*(psnow(ii,jj,bb)%vmid/psnow(ii,jj,cc)%vmid)
        ! Move dN splinters from snow bin cc to snow bin bb
        psnow(ii,jj,bb)%numc = psnow(ii,jj,bb)%numc + dN
        psnow(ii,jj,bb)%volc(:) = psnow(ii,jj,bb)%volc(:) + vfrac*psnow(ii,jj,cc)%volc(:)
        psnow(ii,jj,cc)%volc(:) = (1.-vfrac)*psnow(ii,jj,cc)%volc(:)
    ENDIF
    !
  END SUBROUTINE snow2snow

  SUBROUTINE snow2ice(ii,jj,nice,nsnw,pice,psnow,dN,cc,bb)
    USE mo_submctl, ONLY : t_section
    ! Inputs/outputs
    INTEGER, INTENT(in) :: ii, jj, nice, nsnw ! Dimensions
    INTEGER, INTENT(in) :: cc, bb ! Bin indices for the source snow (cc) and target ice (bb)
    TYPE(t_section), INTENT(inout) :: pice(ii,jj,nice), psnow(ii,jj,nsnw)
    REAL, INTENT(inout) :: dN ! Change in number concentration
    ! Local
    REAL :: vfrac
    !
    ! Limit dN to 10% of the parent number to avoid large changes in mean dry size
    dN=MIN(dN,0.1*psnow(ii,jj,cc)%numc)
    !
    ! Volume fraction to be removed (based on dry size)
    vfrac = dN*pice(ii,jj,bb)%vmid/SUM(psnow(ii,jj,cc)%volc(2:))
    ! Move dN splinters from snow bin cc to ice bin bb
    pice(ii,jj,bb)%numc = pice(ii,jj,bb)%numc + dN
    pice(ii,jj,bb)%volc(:) = pice(ii,jj,bb)%volc(:) + vfrac*psnow(ii,jj,cc)%volc(:)
    psnow(ii,jj,cc)%volc(:) = (1.-vfrac)*psnow(ii,jj,cc)%volc(:)
    !
  END SUBROUTINE snow2ice

END MODULE mo_salsa_cloud
