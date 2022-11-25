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
                              temp,   pres,  rv,     &
                              rs,     paero, pcloud  )

    USE mo_submctl, ONLY : t_section, nbins, ncld, &
              lsactintst, lsactbase

    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::              &
             kbdim,                     & ! dimension for arrays
             klev                       ! number of vertical levels

    REAL, INTENT(in) ::             &
             pres(kbdim,klev),          &
             temp(kbdim,klev)

    REAL, INTENT(inout) :: rv(kbdim,klev) ! Water vapor mixing ratio
    REAL, INTENT(in)    :: rs(kbdim,klev) ! Saturation vapor mixing ratio

    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld),  &
                                      paero(kbdim,klev,nbins)

    ! -------------------------------------
    ! Interstitial activation
    ! -------------------------------------
    IF ( lsactintst ) THEN

       CALL actInterst(kbdim,klev,paero,pcloud,rv,rs,temp)

    END IF

    ! -----------------------------------
    ! Activation at cloud base
    ! -----------------------------------
    IF ( lsactbase ) THEN

        STOP 'Cloud base activation not supported!'

    END IF

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

    TYPE(t_section) :: pactd(ncld) ! Local variable

    REAL :: paa        ! Coefficient for Kelvin effect

    REAL :: zdcstar,zvcstar   ! Critical diameter/volume corresponding to S_LES
    REAL :: zactvol           ! Total volume of the activated particles

    REAL :: Nact, Vact(nspec+1)  ! Helper variables for transferring the activated particles

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

    REAL :: vcut,vint1,vint2  ! cut volume, integration limit volumes
    LOGICAL :: intrange(4)    ! Logical table for integration ranges depending on the shape of the wet size profile:
                                  ! [Vlo -- vint1][vint1 -- Vmid][Vmid -- vint2][vint1 -- Vhi]
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

          ! Loop over cloud droplet (and aerosol) bins
          DO cb = inp2a, fnp2b
             ab = fn1a + cb
             pactd(cb)%numc = 0.d0
             pactd(cb)%volc(:) =0.d0
             ! Dry particles are not activated (volume 1e-28 m^3 is less than that in a 1 nm droplet)
             IF ( paero(ii,jj,ab)%numc < nlim .OR. paero(ii,jj,ab)%volc(1)<paero(ii,jj,ab)%numc*1e-28 ) CYCLE
             intrange = .FALSE.

             ! Define some parameters
             Nmid = paero(ii,jj,ab)%numc     ! Number concentration at the current bin center
             Vwmid = SUM(paero(ii,jj,ab)%volc(1:nn))/Nmid  ! Wet volume at the current bin center
             Vmid = SUM(paero(ii,jj,ab)%volc(2:nn))/Nmid ! Dry volume at the current bin center
             Vlo = Vmid*paero(ii,jj,ab)%vratiolo        ! Dry vol at low limit
             Vhi = Vmid*paero(ii,jj,ab)%vratiohi        ! Dry vol at high limit

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
                      Vim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                      Vwim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                   END IF
                   Vlom1 = Vim1*paero(ii,jj,ab-1)%vratiolo
                   Vhim1 = Vim1*paero(ii,jj,ab-1)%vratiohi
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
                      Vip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                      Vwip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                   END IF
                   Vlop1 = Vip1*paero(ii,jj,ab+1)%vratiolo
                   Vhip1 = Vip1*paero(ii,jj,ab+1)%vratiohi
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

             DO ss = 1,nn
                Vact(ss) = zactvol*( paero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
             END DO

             ! Store the number concentration and mass of activated particles for current bins
             pactd(cb)%numc = MIN(Nact,Nmid)
             pactd(cb)%volc(1:nn) = MIN(Vact(1:nn),paero(ii,jj,ab)%volc(1:nn))

          END DO ! cb

          ! Apply the number and mass activated to aerosol and cloud bins
          paero(ii,jj,in2a:fn2b)%numc = MAX(0., paero(ii,jj,in2a:fn2b)%numc - pactd(inp2a:fnp2b)%numc)
          pcloud(ii,jj,inp2a:fnp2b)%numc = pcloud(ii,jj,inp2a:fnp2b)%numc + pactd(inp2a:fnp2b)%numc
          DO ss = 1,nn
             paero(ii,jj,in2a:fn2b)%volc(ss) = MAX(0., paero(ii,jj,in2a:fn2b)%volc(ss) - pactd(inp2a:fnp2b)%volc(ss))
             pcloud(ii,jj,inp2a:fnp2b)%volc(ss) = pcloud(ii,jj,inp2a:fnp2b)%volc(ss) + pactd(inp2a:fnp2b)%volc(ss)
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
                               pi6,         &
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
            ! calculate number based on that: N*pi/6*Dmid**3=V => N=V/(pi/6*Dmid**3)
            io=iout ! The default bin
            Nrem=Vrem/(pi6*pprecp(ii,jj,io)%dmid**3)
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
                               in2a, fn2b, fnp2a,  &
                               ncld, nprc, nice, nsnw,  &
                               pi, nlim, prlim, &
                               calc_Sw_eq, &
                               ice_hom, ice_imm, ice_dep, &
                               calc_correlation, &
                               ice_target_opt, idu, nspec

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

    INTEGER :: ii,jj,kk,ss,bb,nn
    REAL :: pf_imm, pf_dep, pf_hom, jf
    REAL :: Sw_eq, Si, rn, rw, frac, zvol, ra, rb

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
            IF (pprecp(ii,jj,kk)%numc*frac <prlim) CYCLE

            ! Move to the parallel ice bin or to a snow bin
            IF (ice_target_opt<0) THEN
                ! Move to ice bin with matching dry volume (%vhilim) and select a or b bin based on composition (%volc(2:)).
                zvol = SUM( pprecp(ii,jj,kk)%volc(2:nn) )/pprecp(ii,jj,kk)%numc ! Dry volume

                ! 1) Find the matching bin
                bb=1
                DO WHILE (zvol>pice(ii,jj,bb)%vhilim .AND. bb<fnp2a)
                    bb=bb+1
                ENDDO

                ! 2) Select a or b bin
                IF (pice(ii,jj,bb+fnp2a)%numc<=prlim) THEN
                    ! Empty b bin so select a
                    !bb = bb
                ELSEIF (pice(ii,jj,bb)%numc<=prlim) THEN
                    ! Empty a bin so select b
                    bb = bb + fnp2a
                ELSE
                    ! Both are present - find bin based on compositional similarity
                    ra = calc_correlation(pice(ii,jj,bb)%volc(2:nn),pprecp(ii,jj,kk)%volc(2:nn),nspec)
                    rb = calc_correlation(pice(ii,jj,bb+fnp2a)%volc(2:nn),pprecp(ii,jj,kk)%volc(2:nn),nspec)
                    IF (ra<rb) bb = bb + fnp2a
                ENDIF
                ! Add to the matching ice bin
                pice(ii,jj,bb)%volc(1:nn) = pice(ii,jj,bb)%volc(1:nn) + max(0., pprecp(ii,jj,kk)%volc(1:nn)*frac )
                pice(ii,jj,bb)%numc   = pice(ii,jj,bb)%numc + pprecp(ii,jj,kk)%numc*frac
            ELSEIF (ice_target_opt==0) THEN
                ! Add to the matching snow bin
                ss=1
                zvol=SUM(pprecp(ii,jj,kk)%volc(1:nn))/pprecp(ii,jj,kk)%numc
                DO WHILE (zvol>psnow(ii,jj,ss)%vhilim .AND. ss<nsnw)
                    ss=ss+1
                ENDDO
                psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., pprecp(ii,jj,kk)%volc(1:nn)*frac )
                psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + pprecp(ii,jj,kk)%numc*frac
             ELSE
                ! Add to the ss:th snow bin
                ss=MIN(nsnw,ice_target_opt)
                psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., pprecp(ii,jj,kk)%volc(1:nn)*frac )
                psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + pprecp(ii,jj,kk)%numc*frac
            ENDIF

            pprecp(ii,jj,kk)%numc = pprecp(ii,jj,kk)%numc - pprecp(ii,jj,kk)%numc*frac
            pprecp(ii,jj,kk)%volc(1:nn) = pprecp(ii,jj,kk)%volc(1:nn) - max(0., pprecp(ii,jj,kk)%volc(1:nn)*frac )

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
            IF (pcloud(ii,jj,kk)%numc*frac <prlim) CYCLE

            ! Move to the parallel ice bin or to a snow bin
            IF (ice_target_opt<0) THEN
                ! Add to the matching ice bin
                pice(ii,jj,kk)%volc(1:nn) = pice(ii,jj,kk)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                pice(ii,jj,kk)%numc   = pice(ii,jj,kk)%numc + pcloud(ii,jj,kk)%numc*frac
            ELSEIF (ice_target_opt==0) THEN
                ! Add to the matching snow bin
                ss=1
                zvol=SUM(pcloud(ii,jj,kk)%volc(1:nn))/pcloud(ii,jj,kk)%numc
                DO WHILE (zvol>psnow(ii,jj,ss)%vhilim .AND. ss<nsnw)
                    ss=ss+1
                ENDDO
                psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + pcloud(ii,jj,kk)%numc*frac
             ELSE
                ! Add to the ss:th snow bin
                ss=MIN(nsnw,ice_target_opt)
                psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + pcloud(ii,jj,kk)%numc*frac
            ENDIF

            pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc - pcloud(ii,jj,kk)%numc*frac
            pcloud(ii,jj,kk)%volc(1:nn) = pcloud(ii,jj,kk)%volc(1:nn) - max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )

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
            IF (paero(ii,jj,kk)%numc*frac <prlim) CYCLE

            ! Move to the parallel ice bin or to a snow bin
            IF (ice_target_opt<0) THEN
                ! Add to the matching ice bin
                bb = kk-in2a+1
                pice(ii,jj,bb)%volc(1:nn) = pice(ii,jj,bb)%volc(1:nn) + max(0., paero(ii,jj,kk)%volc(1:nn)*frac )
                pice(ii,jj,bb)%numc   = pice(ii,jj,bb)%numc + max(0., paero(ii,jj,kk)%numc*frac )
            ELSEIF (ice_target_opt==0) THEN
                ! Add to the matching snow bin
                ss=1
                zvol=SUM(paero(ii,jj,kk)%volc(1:nn))/paero(ii,jj,kk)%numc
                DO WHILE (zvol>psnow(ii,jj,ss)%vhilim .AND. ss<nsnw)
                    ss=ss+1
                ENDDO
                psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., paero(ii,jj,kk)%volc(1:nn)*frac )
                psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + max(0., paero(ii,jj,kk)%numc*frac )
            ELSE
                ! Add to the ss:th snow bin
                ss=MIN(nsnw,ice_target_opt)
                psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., paero(ii,jj,kk)%volc(1:nn)*frac )
                psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + max(0., paero(ii,jj,kk)%numc*frac )
            ENDIF

             paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc - max(0., paero(ii,jj,kk)%numc*frac )
             paero(ii,jj,kk)%volc(1:nn) = paero(ii,jj,kk)%volc(1:nn) - max(0., paero(ii,jj,kk)%volc(1:nn)*frac )

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
  ! Prescribed ice number concentration: ice number concentration is increased to the
  ! target concentration (fixinc, #/kg) by converting cloud droplets to ice or snow.
  !
  !***********************************************

  SUBROUTINE fixed_ice_driver(kbdim, klev,   &
            pcloud,  pice,   psnow, ptemp,  ppres,  prv,  prsi)

    USE mo_submctl, ONLY : t_section,      &
                    ncld, nice, nsnw, nspec, &
                    rhowa,  &
                    rda, nlim, prlim, &
                    fixinc, ice_source_opt, ice_target_opt
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,klev
    REAL, INTENT(in) :: &
                    ptemp(kbdim,klev), &
                    ppres(kbdim,klev), &
                    prv(kbdim,klev),   &
                    prsi(kbdim,klev)
    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), &
                    pice(kbdim,klev,nsnw), psnow(kbdim,klev,nsnw)

    ! Limits for ice formation
    !   a) Minimum  water vapor satuturation ratio ove ice (-)
    !   b) Minimum cloud water mixing ratio (kg/kg)
    REAL, PARAMETER :: min_S_ice=1.05, min_rc=1e-6

    INTEGER :: ii,jj,kk,ss,nn
    REAL :: pdn, S_ice, rc, Ni0, vol, sumICE, dnice, frac

    nn = nspec + 1 ! Aerosol species + water

    ! Target snow bin when ice_target_opt>0
    ss=MIN(nsnw,ice_target_opt)

    DO ii = 1,kbdim
    DO jj = 1,klev
        pdn=ppres(ii,jj)/(rda*ptemp(ii,jj)) ! Air density (kg/m^3)

        ! Conditions for ice nucleation
        S_ice = prv(ii,jj)/prsi(ii,jj) ! Saturation with respect to ice
        rc = sum( pcloud(ii,jj,:)%volc(1) )*rhowa/pdn ! Cloud water mixing ratio (kg/kg)
        if ( S_ice < min_S_ice .OR. rc < min_rc ) cycle

        ! Target number concentration of ice, converted to #/m^3
        Ni0 = fixinc * pdn

        ! Current ice number concentration (#/m^3)
        IF (ice_target_opt<0) THEN
            sumICE = SUM( pice(ii,jj,:)%numc )
        ELSE
            sumICE = SUM( psnow(ii,jj,:)%numc )
        ENDIF

        if ( sumICE > Ni0 ) cycle

        IF (ice_source_opt>0) THEN
            ! Activate ice starting from the largest cloud bin
            DO kk = nice,1,-1 ! ncld=nice
                IF( Ni0 - sumICE > prlim .AND. pcloud(ii,jj,kk)%numc > nlim) THEN
                    dnice = MAX( 0.0, MIN( Ni0 - sumICE , pcloud(ii,jj,kk)%numc ) )
                    frac = MAX( 0.0, MIN( 1.0, dnice/pcloud(ii,jj,kk)%numc ) )
                    sumICE = sumICE + dnice

                    IF (ice_target_opt<0) THEN
                        ! Add to the matching ice bin
                        pice(ii,jj,kk)%volc(1:nn) = pice(ii,jj,kk)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                        pice(ii,jj,kk)%numc   = pice(ii,jj,kk)%numc + dnice
                    ELSEIF (ice_target_opt==0) THEN
                        ! Add to the matching snow bin
                        ss=1
                        vol=SUM(pcloud(ii,jj,kk)%volc(1:nn))/pcloud(ii,jj,kk)%numc
                        DO WHILE (vol>psnow(ii,jj,ss)%vhilim .AND. ss<nsnw)
                            ss=ss+1
                        ENDDO
                        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dnice
                    ELSE
                        ! Add to the ss:th snow bin
                        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dnice
                    ENDIF

                    pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc - dnice
                    pcloud(ii,jj,kk)%volc(1:nn) = pcloud(ii,jj,kk)%volc(1:nn) - max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                END IF
            END DO
        ELSEIF (ice_source_opt<0) THEN
            ! Activate ice starting from the smallest cloud bin
            DO kk = 1,nice ! ncld=nice
                IF( Ni0 - sumICE > prlim .AND. pcloud(ii,jj,kk)%numc > nlim) THEN
                    dnice = MAX( 0.0, MIN( Ni0 - sumICE , pcloud(ii,jj,kk)%numc ) )
                    frac = MAX( 0.0, MIN( 1.0, dnice/pcloud(ii,jj,kk)%numc ) )
                    sumICE = sumICE + dnice

                    IF (ice_target_opt<0) THEN
                        ! Add to the matching ice bin
                        pice(ii,jj,kk)%volc(1:nn) = pice(ii,jj,kk)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                        pice(ii,jj,kk)%numc   = pice(ii,jj,kk)%numc + dnice
                    ELSEIF (ice_target_opt==0) THEN
                        ! Add to the matching snow bin
                        ss=1
                        vol=SUM(pcloud(ii,jj,kk)%volc(1:nn))/pcloud(ii,jj,kk)%numc
                        DO WHILE (vol>psnow(ii,jj,ss)%vhilim .AND. ss<nsnw)
                            ss=ss+1
                        ENDDO
                        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dnice
                    ELSE
                        ! Add to the ss:th snow bin
                        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dnice
                    ENDIF

                    pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc - dnice
                    pcloud(ii,jj,kk)%volc(1:nn) = pcloud(ii,jj,kk)%volc(1:nn) - max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                END IF
            END DO
        ELSE
            ! Activate ice from all cloud bins
            frac = max(0.0, min(1.0, (Ni0-sumICE)/SUM(pcloud(ii,jj,:)%numc)))
            DO kk = 1,nice ! ncld=nice
                IF(pcloud(ii,jj,kk)%numc > nlim .AND. frac*pcloud(ii,jj,kk)%numc > prlim) THEN
                    dnice = frac*pcloud(ii,jj,kk)%numc

                    IF (ice_target_opt<0) THEN
                        ! Add to the matching ice bin
                        pice(ii,jj,kk)%volc(1:nn) = pice(ii,jj,kk)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                        pice(ii,jj,kk)%numc   = pice(ii,jj,kk)%numc + dnice
                    ELSEIF (ice_target_opt==0) THEN
                        ! Add to the matching snow bin
                        ss=1
                        vol=SUM(pcloud(ii,jj,kk)%volc(1:nn))/pcloud(ii,jj,kk)%numc
                        DO WHILE (vol>psnow(ii,jj,ss)%vhilim .AND. ss<nsnw)
                            ss=ss+1
                        ENDDO
                        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dnice
                    ELSE
                        ! Add to the ss:th snow bin
                        psnow(ii,jj,ss)%volc(1:nn) = psnow(ii,jj,ss)%volc(1:nn) + max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                        psnow(ii,jj,ss)%numc   = psnow(ii,jj,ss)%numc + dnice
                    ENDIF

                    pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc - dnice
                    pcloud(ii,jj,kk)%volc(1:nn) = pcloud(ii,jj,kk)%volc(1:nn) - max(0., pcloud(ii,jj,kk)%volc(1:nn)*frac )
                END IF
            END DO
        ENDIF
    END DO
    END DO

  END SUBROUTINE fixed_ice_driver
  ! ------------------------------------------------------------


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



END MODULE mo_salsa_cloud
