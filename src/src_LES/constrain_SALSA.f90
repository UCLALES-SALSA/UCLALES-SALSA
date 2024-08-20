MODULE constrain_SALSA
  USE mo_progn_state, ONLY : a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,   &
                             a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,   &
                             a_nicep,  a_nicet,  a_micep,   a_micet, a_gaerop, a_indefp, a_rp
  USE mo_diag_state, ONLY : a_rtot, a_rc, a_srp, a_snrp, a_rh, a_temp, a_ri, a_riri, a_rhi
  USE mo_aux_state, ONLY : aetot
  USE mo_submctl, ONLY : spec, nlim, prlim, ice_theta_dist
  USE grid, ONLY : level
  USE mo_structured_datatypes
  IMPLICIT NONE


  CONTAINS


   SUBROUTINE tend_constrain2
     USE grid, ONLY : dtlt
     USE mo_field_state, ONLY : SALSA_tracers_4d
     IMPLICIT NONE

      TYPE(FloatArray4d), POINTER :: varp => NULL(),vart => NULL()
      INTEGER :: nv
      
      DO nv = 1,SALSA_tracers_4d%count
         CALL SALSA_tracers_4d%getData(1,varp,nv)
         CALL SALSA_tracers_4d%getData(2,vart,nv)
         vart%d(:,:,:,:) = MAX( (-1.+1.e-8)*varp%d(:,:,:,:)/dtlt, vart%d(:,:,:,:)  )
      END DO
      varp => NULL(); vart => NULL()

      IF (ice_theta_dist .AND. level > 4)   &
           a_indefp%d(:,:,:,:) = MIN( MAX( a_indefp%d(:,:,:,:),0. ),1. )
      
   END SUBROUTINE tend_constrain2
   
   !
   ! ---------------------------------------------------------------------
   ! SALSA_diagnostics: Update properties for the current timestep:
   !                    E.g. if enough water has evaporated from droplets,
   !                    deplete the cloud droplet bins and move CCN material
   !                    back to the aerosol regime.
   !                    In addition, update the diagnostic scalars for total grid-cell
   !                    liquid water contents.
   !
   ! Juha Tonttila, FMI, 2014
   ! Tomi Raatikainen, FMI, 2016

   SUBROUTINE SALSA_diagnostics(onlyDiag)
     USE grid, ONLY : nxp,nyp,nzp,    &
                      level
     USE mo_derived_procedures, ONLY : binMixrat ! Maybe at some point one could use the actual variables instead of this subroutine directly?
     USE mo_submctl, ONLY : nbins,ncld,nprc,nice,ica,fca,icb,fcb,ira,fra,              &
                            in1a,fn2b,                        &
                            nice,iia,fia,                    &
                            surfw0, rg, pi, &
                            lscndgas, pi6, avog
     USE util, ONLY : getMassIndex
     IMPLICIT NONE

     LOGICAL, INTENT(in), OPTIONAL :: onlyDiag ! If true, only update diagnostic concentrations and don't do anything else. Default value is FALSE.

     REAL, PARAMETER :: massTH = 1.e-25! Minimum mass threshold; corresponds to about a 1 nm particle == does not make sense
     
     INTEGER :: i,j,k,bc,ba,s,sc,sa,str,end,nc,nspec
     INTEGER :: mi,mi2 !"mass index" variables

     REAL :: zvol
     REAL :: zdh2o
     REAL :: ns, zbb, zaa ! Number of moles, Raoult effect, Kelvin effect; For calculating the critical radius
     REAL :: cdcld,cdprc, cdice! Critical diameter for cloud droplets and precipitation
     REAL :: zdrms, zwams
     REAL :: mice_tot   ! Variable for total ice mass
     REAL :: vice_tot   ! Variable for total ice vol

     LOGICAL :: l_onlyDiag
     
     nspec = spec%getNSpec(type="wet") ! Aerosol species + water. For rime add +1

     IF (PRESENT(onlyDiag)) THEN
        l_onlyDiag = onlyDiag
     ELSE
        l_onlyDiag = .FALSE.
     END IF

     ! Remove negative values
     a_naerop%d = MAX(0.,a_naerop%d)
     a_ncloudp%d = MAX(0.,a_ncloudp%d)
     a_nprecpp%d = MAX(0.,a_nprecpp%d)
     a_maerop%d = MAX(0.,a_maerop%d)
     a_mcloudp%d = MAX(0.,a_mcloudp%d)
     a_mprecpp%d = MAX(0.,a_mprecpp%d)
     
     IF (level == 5) THEN
        a_nicep%d = MAX(0.,a_nicep%d)
        a_micep%d = MAX(0.,a_micep%d)
     END IF
     
     IF ( .NOT. l_onlyDiag) THEN
                   
        ! Remove particles that have number but not mass
        DO j = 3,nyp-2
           DO i = 3,nxp-2
              DO k = 1,nzp
                 ! Aerosols
                 DO bc = 1,nbins
                    mi = getMassIndex(nbins,bc,nspec)
                    IF (a_naerop%d(k,i,j,bc) > 0. .AND. SUM(a_maerop%d(k,i,j,bc:mi:nbins)) <= 0.) THEN
                       a_naerop%d(k,i,j,bc) = 0.
                       a_maerop%d(k,i,j,bc:mi:nbins) = 0.
                    END IF
                 END DO
                 
                 ! Clouds
                 DO bc = 1,ncld
                    mi = getMassIndex(ncld,bc,nspec)
                    IF (a_ncloudp%d(k,i,j,bc) > 0. .AND. SUM(a_mcloudp%d(k,i,j,bc:mi:ncld)) <= 0.) THEN
                       a_ncloudp%d(k,i,j,bc) = 0.
                       a_mcloudp%d(k,i,j,bc:mi:ncld) = 0.
                    END IF
                 END DO ! ncld
                 
                 ! Precipitation
                 DO bc = 1,nprc
                    mi = getMassIndex(nprc,bc,nspec)
                    IF (a_nprecpp%d(k,i,j,bc) > 0. .AND. a_mprecpp%d(k,i,j,mi) <= 0.) THEN
                       a_nprecpp%d(k,i,j,bc) = 0.
                       a_mprecpp%d(k,i,j,bc:mi:nprc) = 0.
                    END IF
                 END DO ! nprc
                 
                 ! Ice
                 IF (level<5) CYCLE
                 DO bc = 1,nice
                    mi = getMassIndex(nice,bc,nspec)     ! Pristine
                    mi2 = getMassIndex(nice,bc,nspec+1)  ! Rime
                    mice_tot = a_micep%d(k,i,j,mi) + a_micep%d(k,i,j,mi2)
                    IF (a_nicep%d(k,i,j,bc) > 0. .AND. mice_tot < massTH) THEN
                       a_nicep%d(k,i,j,bc) = 0.
                       a_micep%d(k,i,j,bc:mi:nice) = 0.
                       a_micep%d(k,i,j,bc:mi2:nice) = 0.
                    END IF
                 END DO ! nice
                 
              END DO !k
           END DO !i
        END DO !j
        
        ! Check evaporation of particles
        ! --------------------------------------
        ! Loop over cloud droplet bins
        DO bc = ica%cur,fcb%cur    
           DO j = 3,nyp-2
              DO i = 3,nxp-2
                 DO k = 1,nzp
                    
                    mi = getMassIndex(ncld,bc,nspec)    ! all compounds
                    mi2 = getMassIndex(ncld,bc,nspec-1) ! Dry aerosol part
                    IF ( a_ncloudp%d(k,i,j,bc) > nlim .AND. a_rh%d(k,i,j)<0.999 .AND.   &
                         a_mcloudp%d(k,i,j,mi) < 1.e-5  ) THEN
                       
                       ! Critical diameter
                       ns = SUM( spec%diss(1:nspec-1)*a_mcloudp%d(k,i,j,bc:mi2:ncld)/spec%MM(1:nspec-1) ) / &
                            a_ncloudp%d(k,i,j,bc)
                       zbb = 6.*spec%mwa*ns/(pi*spec%rhowa)
                       zaa = 4.*spec%mwa*surfw0/(rg*spec%rhowa*a_temp%d(k,i,j))
                       cdcld = SQRT(3.*zbb/zaa)
                       
                       ! Wet diameter
                       zvol = SUM( a_mcloudp%d(k,i,j,bc:mi:ncld)/spec%rholiq(1:nspec) )/a_ncloudp%d(k,i,j,bc)
                       zdh2o = (zvol/pi6)**(1./3.)
                       
                       ! Lose the droplets if small
                       IF ( zdh2o < MIN(0.7*cdcld,1.3e-5) .OR.     &
                            a_mcloudp%d(k,i,j,mi) < massTH*a_ncloudp%d(k,i,j,bc) )  THEN
                          
                          IF (bc<=fca%cur) THEN
                             ba = ica%par + (bc-ica%cur) ! Index for parallel aerosol bin
                          ELSE
                             ba = icb%par + (bc-icb%cur) ! Index for parallel aerosol bin
                          END IF
                          ! Move the number of particles from cloud to aerosol bins
                          a_naerop%d(k,i,j,ba) = a_naerop%d(k,i,j,ba) + a_ncloudp%d(k,i,j,bc)
                          a_ncloudp%d(k,i,j,bc) = 0.
                          
                          ! Move ccn material back to aerosol regime (including water)
                          DO s = 1,nspec
                             sc = getMassIndex(ncld,bc,s)
                             sa = getMassIndex(nbins,ba,s)
                             a_maerop%d(k,i,j,sa) = a_maerop%d(k,i,j,sa) + a_mcloudp%d(k,i,j,sc)
                             a_mcloudp%d(k,i,j,sc) = 0.
                          END DO
                          
                       END IF ! critical diameter
                    END IF  ! nlim
                    
                 END DO
              END DO
           END DO
        END DO ! bc
        
        ! Loop over precipitation bins
        DO bc = 1,nprc
           DO j = 3,nyp-2
              DO i = 3,nxp-2
                 DO k = 1,nzp
                    
                    mi = getMassIndex(nprc,bc,nspec)    ! all
                    mi2 = getMassIndex(nprc,bc,nspec-1) ! Dry aerosol part
                    IF ( a_nprecpp%d(k,i,j,bc) > prlim .AND. a_rh%d(k,i,j)<0.999 .AND.  &
                         a_mprecpp%d(k,i,j,mi) < 1.e-6 ) THEN
                       
                       ! Critical diameter
                       ns = SUM( spec%diss(1:nspec-1)*a_mprecpp%d(k,i,j,bc:mi2:nprc)/spec%MM(1:nspec-1) ) / &
                            a_nprecpp%d(k,i,j,bc)
                       zbb = 6.*spec%mwa*ns/(pi*spec%rhowa)
                       zaa = 4.*spec%mwa*surfw0/(rg*spec%rhowa*a_temp%d(k,i,j))
                       cdprc = SQRT(3.*zbb/zaa)
                       
                       ! Wet diameter
                       zvol = SUM( a_mprecpp%d(k,i,j,bc:mi:nprc)/spec%rholiq(1:nspec) )/a_nprecpp%d(k,i,j,bc)
                       zdh2o = (zvol/pi6)**(1./3.)
                       
                       ! Lose the droplets if small
                       IF ( zdh2o < MIN(0.2*cdprc,10.e-6) .OR.   &
                            a_mprecpp%d(k,i,j,mi) < massTH*a_nprecpp%d(k,i,j,bc) ) THEN
                          
                          ba = findDry4Wet(a_nprecpp,a_mprecpp,nprc,bc,k,i,j)
                          
                          ! Move the number of particles from cloud to aerosol bins
                          a_naerop%d(k,i,j,ba) = a_naerop%d(k,i,j,ba) + a_nprecpp%d(k,i,j,bc)
                          a_nprecpp%d(k,i,j,bc) = 0.
                          
                          ! Move ccn material back to aerosol regime (including water)
                          DO s = 1,nspec
                             sc = getMassIndex(nprc,bc,s)
                             sa = getMassIndex(nbins,ba,s)
                             a_maerop%d(k,i,j,sa) = a_maerop%d(k,i,j,sa) + a_mprecpp%d(k,i,j,sc)
                             a_mprecpp%d(k,i,j,sc) = 0.
                          END DO
                          
                       END IF ! Critical diameter
                       
                    END IF ! prlim
                    
                 END DO
              END DO
           END DO
        END DO ! bc
        
        ! Loop over ice bins
        IF (level == 5) THEN
           DO bc = 1,nice
              DO j = 3,nyp-2
                 DO i = 3,nxp-2
                    DO k = 1,nzp
                       
                       mi = getMassIndex(nice,bc,nspec)     ! Pristine ice
                       mi2 = getMassIndex(nice,bc,nspec+1)  ! rimed ice
                       mice_tot = a_micep%d(k,i,j,mi) + a_micep%d(k,i,j,mi2)
                       vice_tot = SUM(a_micep%d(k,i,j,bc:mi2:nice)/spec%rhoice(1:nspec+1))
                       
                       IF ( a_nicep%d(k,i,j,bc) > prlim .AND. a_rhi%d(k,i,j)<0.999 .AND.   &
                            mice_tot < 1.e-15 ) THEN
                          
                          ! Ice and snow don't have a critical size, but lose particles when water content becomes low enough or size is < 2e-6m
                          
                          ! Diameter          
                          cdice = vice_tot / a_nicep%d(k,i,j,bc)
                          cdice = (cdice/pi6)**(1./3.)
                          
                          ! Lose ice when dry to total mass ratio is more than 0.5
                          CALL binMixrat("ice","dry",bc,i,j,k,zdrms)
                          CALL binMixrat("ice","wet",bc,i,j,k,zwams)
                          zvol = zdrms/zwams
                          
                          IF ( zvol>0.5 .OR. cdice < 2.e-6 ) THEN
                             
                             ! Find the aerosol bin corresponding to the composition of the IN the current bin
                             ba = findDry4Wet(a_nicep,a_micep,nice,bc,k,i,j)                     
                             
                             ! Move the number of particles from ice to aerosol bins
                             a_naerop%d(k,i,j,ba) = a_naerop%d(k,i,j,ba) + a_nicep%d(k,i,j,bc)
                             a_nicep%d(k,i,j,bc) = 0.
                             
                             ! Move mass material back to aerosol regime (including water)
                             DO s = 1,nspec
                                sc = getMassIndex(nice,bc,s)
                                sa = getMassIndex(nbins,ba,s)
                                a_maerop%d(k,i,j,sa) = a_maerop%d(k,i,j,sa) + a_micep%d(k,i,j,sc)
                                a_micep%d(k,i,j,sc) = 0.
                             END DO
                             ! Rimed ice
                             sa = getMassIndex(nbins,ba,nspec)
                             sc = getMassIndex(nice,bc,nspec+1)
                             a_maerop%d(k,i,j,sa) = a_maerop%d(k,i,j,sa) + a_micep%d(k,i,j,sc)
                             a_micep%d(k,i,j,sc) = 0.
                          END IF
                          
                       END IF  ! prlim
                       
                    END DO
                 END DO
              END DO
           END DO ! bc
        END IF ! level = 5
           
        ! Loop over aerosol bins
        DO ba = 1,nbins
           DO j = 3,nyp-2
              DO i = 3,nxp-2
                 DO k = 1,nzp
                    
                    mi = getMassIndex(nbins,ba,nspec)    ! all
                    mi2 = getMassIndex(nbins,ba,nspec-1) ! dry part
                    IF (a_naerop%d(k,i,j,ba) > nlim) THEN
                       zvol = SUM( a_maerop%d(k,i,j,ba:mi2:nbins)/spec%rholiq(1:nspec-1) )/ &
                            a_naerop%d(k,i,j,ba) ! Dry volume
                       
                       ! Particles smaller then 0.1 nm diameter are set to zero 
                       IF ( zvol < pi6*1.e-10**3 ) THEN
                          ! Volatile species to the gas phase
                          IF (spec%isUsed('SO4') .AND. lscndgas) THEN
                             nc = spec%getIndex('SO4')
                             s = getMassIndex(nbins,ba,nc)
                             a_gaerop%d(k,i,j,1) = a_gaerop%d(k,i,j,1) + a_maerop%d(k,i,j,s) / spec%msu * avog
                          END IF
                          IF (spec%isUsed('OC') .AND. lscndgas) THEN
                             nc = spec%getIndex('OC')
                             s = getMassIndex(nbins,ba,nc)
                             a_gaerop%d(k,i,j,5) = a_gaerop%d(k,i,j,5) + a_maerop%d(k,i,j,s) / spec%moc * avog
                          END IF
                          IF (spec%isUsed('NO') .AND. lscndgas) THEN
                             nc = spec%getIndex('NO')
                             s = getMassIndex(nbins,ba,nc)
                             a_gaerop%d(k,i,j,2) = a_gaerop%d(k,i,j,2) + a_maerop%d(k,i,j,s) / spec%mno * avog
                          END IF
                          IF (spec%isUsed('NH') .AND. lscndgas) THEN
                             nc = spec%getIndex('NH')
                             s = getMassIndex(nbins,ba,nc)
                             a_gaerop%d(k,i,j,3) = a_gaerop%d(k,i,j,3) + a_maerop%d(k,i,j,s) / spec%mnh * avog
                          END IF
                          
                          ! Mass and number to zero (insolube species and water are lost)
                          a_maerop%d(k,i,j,ba:mi:nbins) = 0.
                          a_naerop%d(k,i,j,ba) = 0.
                       END IF
                    END IF
                    
                 END DO
              END DO
           END DO
        END DO
        
     END IF ! onlyDiag

    !!!!!!!!!!!!!!!!!!!!!!!
    ! Update diagnostic tracers
    !!!!!!!!!!!!!!!!!!!!!!!

    ! Liquid water content
    nc = spec%getIndex('H2O')
    ! Aerosols, regimes a and b
    str = getMassIndex(nbins,in1a,nc)
    end = getMassIndex(nbins,fn2b,nc)
    a_rc%d(:,:,:) = SUM(a_maerop%d(:,:,:,str:end),DIM=4)
    ! Clouds, regime a and b
    str = getMassIndex(ncld,ica%cur,nc)
    end = getMassIndex(ncld,fcb%cur,nc)
    a_rc%d(:,:,:) = a_rc%d(:,:,:) + SUM(a_mcloudp%d(:,:,:,str:end),DIM=4)
    ! Precipitation
    str = getMassIndex(nprc,ira,nc)
    end = getMassIndex(nprc,fra,nc)
    a_srp%d(:,:,:) = SUM(a_mprecpp%d(:,:,:,str:end),DIM=4)
    a_snrp%d(:,:,:) = SUM(a_nprecpp%d(:,:,:,ira:fra),DIM=4)
    ! Total water mix rat
    a_rtot%d = a_rc%d + a_srp%d + a_rp%d

    ! ice
    IF (level == 5) THEN 
       str = getMassIndex(nice,iia,nc)
       end = getMassIndex(nice,fia,nc)
       a_ri%d(:,:,:) = SUM(a_micep%d(:,:,:,str:end),DIM=4)
       ! Rimed ice
       nc = spec%getIndex("rime")
       str = getMassIndex(nice,iia,nc)
       end = getMassIndex(nice,fia,nc)
       a_riri%d(:,:,:) = SUM(a_micep%d(:,:,:,str:end),DIM=4)
       ! Update total water
       a_rtot%d = a_rtot%d + a_ri%d + a_riri%d       
    END IF
        
  END SUBROUTINE SALSA_diagnostics

  !
  ! ----------------------------------------------------------------------
  !

  INTEGER FUNCTION findDry4Wet(nevap,mevap,nb,ib,iz,ix,iy)
    USE mo_submctl, ONLY : in2a,fn2a,nbins,pi6
    USE util, ONLY :  calc_correlation, getMassIndex
    IMPLICIT NONE
    
    INTEGER, INTENT(in)            :: nb  ! Number of bins and number of species in the evaporating particle class 
    INTEGER, INTENT(in)            :: ib,iz,ix,iy  ! Current bin and grid indices
    TYPE(FloatArray4d), INTENT(in) :: mevap
    TYPE(FloatArray4d), INTENT(in) :: nevap

    REAL :: cd        ! dry particle diameter
    INTEGER :: mi,mi2 ! Mass indices
    INTEGER :: ba,bb  ! Corresponding aerosol indices for regimes a and b
    REAL :: ra, rb ! Correlation coefficients for a and b aerosol bins

    INTEGER :: ndry
    
    ! This function finds a suitable aerosol bin for evaporating
    ! precipitation or ice bins
    
    ndry = spec%getNSpec(type="dry")
    
    mi = getMassIndex(nb,ib,ndry)  ! Final dry mass index
                           
    ! 1) Find the closest matching bin based on dry particle diameter (a and b bins)
    !    -- note that spec%rholiq for all evaporating categories is ok here,
    !    -- since the aerosol densities are always the same
    cd = SUM( mevap%d(iz,ix,iy,ib:mi:nb)/spec%rholiq(1:ndry) ) / &
              nevap%d(iz,ix,iy,ib)
    cd = (cd/pi6)**(1./3.) ! Dry diameter
    
    ba=in2a !Ignore 1a and note that "aerobins" contains the lower limit of bin dry diameter
    ba = MIN( MAX(in2a + COUNT( cd > aetot%d(ba:fn2a) )-1, in2a), fn2a )

    ! Corresponding b bin is ba+(fn2a-fn1a)
    bb = fn2a + (ba - in2a)
    
    ! 2) Select a or b bin
    IF (a_naerop%d(iz,ix,iy,bb) <= nlim) THEN
       ! Empty b bin so select a, i.e. do nothing here
       findDry4Wet = ba
    ELSE IF (a_naerop%d(iz,ix,iy,ba) <= nlim) THEN
       ! Empty a bin so select b
       findDry4Wet = bb
    ELSE
       mi = getMassIndex(nbins,ba,ndry) ! Index of the last dry species for current bin in aerosol regime a
       mi2 = getMassIndex(nb,ib,ndry)   ! The same for the evaporating particle
       
       ! Both are present - find bin based on compositional similarity
       ra = calc_correlation(a_maerop%d(iz,ix,iy,ba:mi:nbins),  &
                             mevap%d(iz,ix,iy,ib:mi2:nb),   &
                             ndry                       )
                         
       mi = getMassIndex(nbins,bb,ndry) ! Index of the last dry species for current bin in aerosol regime b
       rb = calc_correlation(a_maerop%d(iz,ix,iy,bb:mi:nbins),  &
                             mevap%d(iz,ix,iy,ib:mi2:nb),   &
                             ndry                       )

       findDry4Wet = ba
       IF (ra<rb) findDry4Wet = bb
    END IF

  END FUNCTION findDry4Wet
  

END MODULE constrain_SALSA
