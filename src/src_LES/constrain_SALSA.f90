MODULE constrain_SALSA
  USE grid, ONLY : a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,   &
                   a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,   &
                   a_nicep,  a_nicet,  a_micep,   a_micet, a_salsat,a_salsap
  USE mo_submctl, ONLY : spec, aerobins, nlim, prlim
  
  IMPLICIT NONE


  CONTAINS

   !
   !----------------------------------------------------------------------
   ! In case of negative tendencies to SALSA arrays, put some constrains
   ! in order to avoid concentrations going negative. This will possibly
   ! slightly affect the conservation of mass - needs testing/revision
   ! Juha Tonttila, FMI, 2014
   !
   SUBROUTINE tend_constrain(nn)

      USE grid, ONLY : dtlt, nxp,nyp,nzp, level
      USE mo_submctl, ONLY : nbins, ncld, nprc, &
                             nice       !ice'n'snow

      INTEGER, INTENT(in) :: nn

      INTEGER :: cc, ii,jj,kk,ni

      DO jj = 3, nyp-2

         DO ii = 3, nxp-2

            DO kk = 1, nzp

               ! Aerosols
               DO cc = 1, nbins

                  IF ( a_naerop(kk,ii,jj,cc)+a_naerot(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_naerot(kk,ii,jj,cc) = MAX( (((1.e-10)-1.0)*a_naerop(kk,ii,jj,cc))/dtlt,a_naerot(kk,ii,jj,cc) )
                     DO ni = 1, nn
                        a_maerot(kk,ii,jj,(ni-1)*nbins+cc) = MAX( (((1.e-10)-1.0)*a_maerop(kk,ii,jj,(ni-1)*nbins+cc))/dtlt,  &
                                                                 a_maerot(kk,ii,jj,(ni-1)*nbins+cc) )
                     END DO

                  END IF

               END DO

               ! Cloud droplets
               DO cc = 1, ncld

                  IF ( a_ncloudp(kk,ii,jj,cc)+a_ncloudt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_ncloudt(kk,ii,jj,cc) = MAX( (((1.e-10)-1.0)*a_ncloudp(kk,ii,jj,cc))/dtlt,a_ncloudt(kk,ii,jj,cc) )
                     DO ni = 1, nn
                        a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) = MAX( (((1.e-10)-1.0)*a_mcloudp(kk,ii,jj,(ni-1)*ncld+cc))/dtlt,  &
                                                                 a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) )
                     END DO

                  END IF

               END DO

               ! Precipitation
               DO cc = 1, nprc

                  IF ( a_nprecpp(kk,ii,jj,cc)+a_nprecpt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nprecpt(kk,ii,jj,cc) = MAX( (((1.e-10)-1.0)*a_nprecpp(kk,ii,jj,cc))/dtlt,a_nprecpt(kk,ii,jj,cc) )
                     DO ni = 1, nn
                        a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) = MAX( (((1.e-10)-1.0)*a_mprecpp(kk,ii,jj,(ni-1)*nprc+cc))/dtlt,  &
                                                                 a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) )
                     END DO

                  END IF

               END DO

               ! ice particles
               IF (level<5) CYCLE
               DO cc = 1,nice

                  IF ( a_nicep(kk,ii,jj,cc)+a_nicet(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nicet(kk,ii,jj,cc) = MAX( (((1.e-10)-1.0)*a_nicep(kk,ii,jj,cc))/dtlt,a_nicet(kk,ii,jj,cc) )
                     DO ni = 1, nn
                        a_micet(kk,ii,jj,(ni-1)*nice+cc) = MAX( (((1.e-10)-1.0)*a_micep(kk,ii,jj,(ni-1)*nice+cc))/dtlt,  &
                                                               a_micet(kk,ii,jj,(ni-1)*nice+cc) )
                     END DO
                     ! Rime
                     a_micet(kk,ii,jj,nn*nice+cc) = MAX( (((1.e-10)-1.0)*a_micep(kk,ii,jj,nn*nice+cc))/dtlt,  &
                                                         a_micet(kk,ii,jj,nn*nice+cc) )
                  END IF
                  
               END DO

            END DO ! kk

         END DO ! ii

      END DO ! jj

   END SUBROUTINE tend_constrain


   SUBROUTINE tend_constrain2
      USE grid, ONLY : dtlt, nxp,nyp,nzp,nsalsa,level
      USE mo_submctl, ONLY : nbins, ncld, nprc, &
                             nice    
      IMPLICIT NONE

      INTEGER :: nb,k,i,j
            
      DO nb = 1,nsalsa
         DO j = 1,nyp
            DO i = 1,nxp
               DO k = 1,nzp
                  a_salsat(k,i,j,nb) = MAX( ((1.e-10)-1.)*a_salsap(k,i,j,nb)/dtlt, a_salsat(k,i,j,nb) )
               END DO
            END DO
         END DO
      END DO
          
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

   SUBROUTINE SALSA_diagnostics(callno,onlyDiag)
     USE grid, ONLY : nxp,nyp,nzp,    &
                      a_gaerop, &
                      a_rc, a_srp,a_snrp,  &
                      a_rh, a_temp, a_ri, a_riri,a_rhi,                                      &
                      a_nicep,a_micep, level,   &
                      binMixrat
     USE mo_submctl, ONLY : nbins,ncld,nprc,nice,ica,fca,icb,fcb,ira,fra,              &
                            in1a,in2a,fn2a,fn2b,                        &
                            nice,iia,fia,                    &
                            surfw0, rg, pi, &
                            lscndgas, pi6, avog
     USE util, ONLY : getMassIndex
     IMPLICIT NONE

     INTEGER, INTENT(in) :: callno
     LOGICAL, INTENT(in), OPTIONAL :: onlyDiag ! If true, only update diagnostic concentrations and don't do anything else. Default value is FALSE.

     REAL, PARAMETER :: massTH = 1.e-25! Minimum mass threshold; corresponds to about a 1 nm particle == does not make sense
     
     INTEGER :: i,j,k,bb,bc,ba,s,sc,sa,str,end,nc,nspec
     INTEGER :: mi,mi2 !"mass index" variables

     REAL :: ra, rb

     REAL :: zvol
     REAL :: zdh2o
     REAL :: ns, zbb, zaa ! Number of moles, Raoult effect, Kelvin effect; For calculating the critical radius
     REAL :: cdcld,cdprc, cdice! Critical diameter for cloud droplets and precipitation
     REAL :: zdrms, zwams
     REAL :: mice_tot   ! Variable for total ice mass
     REAL :: vice_tot   ! Variable for total ice vol

     LOGICAL :: l_onlyDiag
     
     nspec = spec%getNSpec() ! total number of species

     IF (PRESENT(onlyDiag)) THEN
        l_onlyDiag = onlyDiag
     ELSE
        l_onlyDiag = .FALSE.
     END IF
     

     IF ( .NOT. l_onlyDiag) THEN
        
        ! Remove negative values
        a_naerop = MAX(0.,a_naerop)
        a_ncloudp = MAX(0.,a_ncloudp)
        a_nprecpp = MAX(0.,a_nprecpp)
        a_maerop = MAX(0.,a_maerop)
        a_mcloudp = MAX(0.,a_mcloudp)
        a_mprecpp = MAX(0.,a_mprecpp)
        
        a_nicep = MAX(0.,a_nicep)
        a_micep = MAX(0.,a_micep)
        
        ! Remove particles that have number but not mass
        DO j = 3,nyp-2
           DO i = 3,nxp-2
              DO k = 1,nzp
                 ! Aerosols
                 DO bc = 1,nbins
                    mi = getMassIndex(nbins,bc,nspec)
                    IF (a_naerop(k,i,j,bc) > 0. .AND. SUM(a_maerop(k,i,j,bc:mi:nbins)) <= 0.) THEN
                       a_naerop(k,i,j,bc) = 0.
                       a_maerop(k,i,j,bc:mi:nbins) = 0.
                    END IF
                 END DO
                 
                 ! Clouds
                 DO bc = 1,ncld
                    mi = getMassIndex(ncld,bc,nspec)
                    IF (a_ncloudp(k,i,j,bc) > 0. .AND. SUM(a_mcloudp(k,i,j,bc:mi:ncld)) <= 0.) THEN
                       a_ncloudp(k,i,j,bc) = 0.
                       a_mcloudp(k,i,j,bc:mi:ncld) = 0.
                    END IF
                 END DO ! ncld
                 
                 ! Precipitation
                 DO bc = 1,nprc
                    mi = getMassIndex(nprc,bc,nspec)
                    IF (a_nprecpp(k,i,j,bc) > 0. .AND. a_mprecpp(k,i,j,mi) <= 0.) THEN
                       a_nprecpp(k,i,j,bc) = 0.
                       a_mprecpp(k,i,j,bc:mi:nprc) = 0.
                    END IF
                 END DO ! nprc
                 
                 ! Ice
                 IF (level<5) CYCLE
                 DO bc = 1,nice
                    mi = getMassIndex(nice,bc,nspec)     ! Pristine
                    mi2 = getMassIndex(nice,bc,nspec+1)  ! Rime
                    mice_tot = a_micep(k,i,j,mi) + a_micep(k,i,j,mi2)
                    IF (a_nicep(k,i,j,bc) > 0. .AND. mice_tot < massTH) THEN
                       a_nicep(k,i,j,bc) = 0.
                       a_micep(k,i,j,bc:mi:nice) = 0.
                       a_micep(k,i,j,bc:mi2:nice) = 0.
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
                    IF ( a_ncloudp(k,i,j,bc) > nlim .AND. a_rh(k,i,j)<0.999 .AND.   &
                         a_mcloudp(k,i,j,mi) < 1.e-5  ) THEN
                       
                       ! Critical diameter
                       ns = SUM( spec%diss(1:nspec-1)*a_mcloudp(k,i,j,bc:mi2:ncld)/spec%MM(1:nspec-1) ) / &
                            a_ncloudp(k,i,j,bc)
                       zbb = 6.*spec%mwa*ns/(pi*spec%rhowa)
                       zaa = 4.*spec%mwa*surfw0/(rg*spec%rhowa*a_temp(k,i,j))
                       cdcld = SQRT(3.*zbb/zaa)
                       
                       ! Wet diameter
                       zvol = SUM( a_mcloudp(k,i,j,bc:mi:ncld)/spec%rholiq(1:nspec) )/a_ncloudp(k,i,j,bc)
                       zdh2o = (zvol/pi6)**(1./3.)
                       
                       ! Lose the droplets if small
                       IF ( zdh2o < MIN(0.7*cdcld,1.3e-5) .OR.     &
                            a_mcloudp(k,i,j,mi) < massTH*a_ncloudp(k,i,j,bc) )  THEN
                          
                          IF (bc<=fca%cur) THEN
                             ba = ica%par + (bc-ica%cur) ! Index for parallel aerosol bin
                          ELSE
                             ba = icb%par + (bc-icb%cur) ! Index for parallel aerosol bin
                          END IF
                          ! Move the number of particles from cloud to aerosol bins
                          a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_ncloudp(k,i,j,bc)
                          a_ncloudp(k,i,j,bc) = 0.
                          
                          ! Move ccn material back to aerosol regime (including water)
                          DO s = 1,nspec
                             sc = getMassIndex(ncld,bc,s)
                             sa = getMassIndex(nbins,ba,s)
                             a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mcloudp(k,i,j,sc)
                             a_mcloudp(k,i,j,sc) = 0.
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
                    IF ( a_nprecpp(k,i,j,bc) > prlim .AND. a_rh(k,i,j)<0.999 .AND.  &
                         a_mprecpp(k,i,j,mi) < 1.e-6 ) THEN
                       
                       ! Critical diameter
                       ns = SUM( spec%diss(1:nspec-1)*a_mprecpp(k,i,j,bc:mi2:nprc)/spec%MM(1:nspec-1) ) / &
                            a_nprecpp(k,i,j,bc)
                       zbb = 6.*spec%mwa*ns/(pi*spec%rhowa)
                       zaa = 4.*spec%mwa*surfw0/(rg*spec%rhowa*a_temp(k,i,j))
                       cdprc = SQRT(3.*zbb/zaa)
                       
                       ! Wet diameter
                       zvol = SUM( a_mprecpp(k,i,j,bc:mi:nprc)/spec%rholiq(1:nspec) )/a_nprecpp(k,i,j,bc)
                       zdh2o = (zvol/pi6)**(1./3.)
                       
                       ! Lose the droplets if small
                       IF ( zdh2o < MAX(0.1*cdprc,1.3e-5) .OR.   &
                            a_mprecpp(k,i,j,mi) < massTH*a_nprecpp(k,i,j,bc) ) THEN
                          
                          ba = findDry4Wet(a_nprecpp,a_mprecpp,nprc,nspec,bc,k,i,j,3)
                          
                          ! Move the number of particles from cloud to aerosol bins
                          a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nprecpp(k,i,j,bc)
                          a_nprecpp(k,i,j,bc) = 0.
                          
                          ! Move ccn material back to aerosol regime (including water)
                          DO s = 1,nspec
                             sc = getMassIndex(nprc,bc,s)
                             sa = getMassIndex(nbins,ba,s)
                             a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mprecpp(k,i,j,sc)
                             a_mprecpp(k,i,j,sc) = 0.
                          END DO
                          
                       END IF ! Critical diameter
                       
                    END IF ! prlim
                    
                 END DO
              END DO
           END DO
        END DO ! bc
        
        ! Loop over ice bins
        DO bc = 1,nice
           DO j = 3,nyp-2
              DO i = 3,nxp-2
                 DO k = 1,nzp
                    
                    mi = getMassIndex(nice,bc,nspec)     ! Pristine ice
                    mi2 = getMassIndex(nice,bc,nspec+1)  ! rimed ice
                    mice_tot = a_micep(k,i,j,mi) + a_micep(k,i,j,mi2)
                    vice_tot = SUM(a_micep(k,i,j,bc:mi:nice)/spec%rhoice(1:nspec)) + &
                         a_micep(k,i,j,mi2)/spec%rhori
                    
                    IF ( a_nicep(k,i,j,bc) > prlim .AND. a_rhi(k,i,j)<0.999 .AND.   &
                         mice_tot < 1.e-15 ) THEN
                       
                       ! Ice and snow don't have a critical size, but lose particles when water content becomes low enough or size is < 2e-6m
                       
                       ! Diameter          
                       cdice = vice_tot / a_nicep(k,i,j,bc)
                       cdice = (cdice/pi6)**(1./3.)
                       
                       ! Lose ice when dry to total mass ratio is more than 0.5
                       CALL binMixrat("ice","dry",bc,i,j,k,zdrms)
                       CALL binMixrat("ice","wet",bc,i,j,k,zwams)
                       zvol = zdrms/zwams
                       
                       IF ( zvol>0.5 .OR. cdice < 2.e-6 ) THEN
                          
                          ! Find the aerosol bin corresponding to the composition of the IN the current bin
                          ba = findDry4Wet(a_nicep,a_micep,nice,nspec+1,bc,k,i,j,4)                     
                          
                          ! Move the number of particles from ice to aerosol bins
                          a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nicep(k,i,j,bc)
                          a_nicep(k,i,j,bc) = 0.
                          
                          ! Move mass material back to aerosol regime (including water)
                          DO s = 1,nspec
                             sc = getMassIndex(nice,bc,s)
                             sa = getMassIndex(nbins,ba,s)
                             a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_micep(k,i,j,sc)
                             a_micep(k,i,j,sc) = 0.
                          END DO
                          ! Rimed ice
                          sa = getMassIndex(nbins,ba,nspec)
                          sc = getMassIndex(nice,bc,nspec+1)
                          a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_micep(k,i,j,sc)
                          a_micep(k,i,j,sc) = 0.
                       END IF
                       
                    END IF  ! prlim
                    
                 END DO
              END DO
           END DO
        END DO ! bc
        
        ! Loop over aerosol bins
        DO ba = 1,nbins
           DO j = 3,nyp-2
              DO i = 3,nxp-2
                 DO k = 1,nzp
                    
                    mi = getMassIndex(nbins,ba,nspec)    ! all
                    mi2 = getMassIndex(nbins,ba,nspec-1) ! dry part
                    IF (a_naerop(k,i,j,ba) > nlim) THEN
                       zvol = SUM( a_maerop(k,i,j,ba:mi2:nbins)/spec%rholiq(1:nspec-1) )/ &
                            a_naerop(k,i,j,ba) ! Dry volume
                       
                       ! Particles smaller then 0.1 nm diameter are set to zero 
                       IF ( zvol < pi6*1.e-10**3 ) THEN
                          ! Volatile species to the gas phase
                          IF (spec%isUsed('SO4') .AND. lscndgas) THEN
                             nc = spec%getIndex('SO4')
                             s = getMassIndex(nbins,ba,nc)
                             a_gaerop(k,i,j,1) = a_gaerop(k,i,j,1) + a_maerop(k,i,j,s) / spec%msu * avog
                          END IF
                          IF (spec%isUsed('OC') .AND. lscndgas) THEN
                             nc = spec%getIndex('OC')
                             s = getMassIndex(nbins,ba,nc)
                             a_gaerop(k,i,j,5) = a_gaerop(k,i,j,5) + a_maerop(k,i,j,s) / spec%moc * avog
                          END IF
                          IF (spec%isUsed('NO') .AND. lscndgas) THEN
                             nc = spec%getIndex('NO')
                             s = getMassIndex(nbins,ba,nc)
                             a_gaerop(k,i,j,2) = a_gaerop(k,i,j,2) + a_maerop(k,i,j,s) / spec%mno * avog
                          END IF
                          IF (spec%isUsed('NH') .AND. lscndgas) THEN
                             nc = spec%getIndex('NH')
                             s = getMassIndex(nbins,ba,nc)
                             a_gaerop(k,i,j,3) = a_gaerop(k,i,j,3) + a_maerop(k,i,j,s) / spec%mnh * avog
                          END IF
                          
                          ! Mass and number to zero (insolube species and water are lost)
                          a_maerop(k,i,j,ba:mi:nbins) = 0.
                          a_naerop(k,i,j,ba) = 0.
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
    a_rc(:,:,:) = SUM(a_maerop(:,:,:,str:end),DIM=4)
    ! Clouds, regime a and b
    str = getMassIndex(ncld,ica%cur,nc)
    end = getMassIndex(ncld,fcb%cur,nc)
    a_rc(:,:,:) = a_rc(:,:,:) + SUM(a_mcloudp(:,:,:,str:end),DIM=4)
    ! Precipitation
    str = getMassIndex(nprc,ira,nc)
    end = getMassIndex(nprc,fra,nc)
    a_srp(:,:,:) = SUM(a_mprecpp(:,:,:,str:end),DIM=4)
    a_snrp(:,:,:) = SUM(a_nprecpp(:,:,:,ira:fra),DIM=4)
    ! ice
    str = getMassIndex(nice,iia,nc)
    end = getMassIndex(nice,fia,nc)
    a_ri(:,:,:) = SUM(a_micep(:,:,:,str:end),DIM=4)
    ! Rimed ice
    str = getMassIndex(nice,iia,nc+1)
    end = getMassIndex(nice,fia,nc+1)
    a_riri(:,:,:) = SUM(a_micep(:,:,:,str:end),DIM=4)

  END SUBROUTINE SALSA_diagnostics


  !
  ! ----------------------------------------------------------------------
  !
  INTEGER FUNCTION findDry4Wet(nevap,mevap,nb,nsp,ib,iz,ix,iy,iphase)
    USE grid, ONLY : nzp,nxp,nyp
    USE mo_submctl, ONLY : in2a,fn2a,fn1a,nbins,pi6
    USE util, ONLY :  calc_correlation, getMassIndex
    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: iphase    ! Evaporating particle phase: 3: precip, 4: ice 
    INTEGER, INTENT(in) :: nb,nsp  ! Number of bins and number of species in the evaporating particle class 
    INTEGER, INTENT(in) :: ib,iz,ix,iy  ! Current bin and grid indices
    REAL, INTENT(in) :: mevap(nzp,nxp,nyp,nb*nsp)
    REAL, INTENT(in) :: nevap(nzp,nxp,nyp,nb)

    REAL :: cd        ! dry particle diameter
    INTEGER :: mi,mi2 ! Mass indices
    INTEGER :: ba,bb  ! Corresponding aerosol indices for regimes a and b
    REAL :: ra, rb ! Correlation coefficients for a and b aerosol bins

    INTEGER :: nspec
    
    ! This function finds a suitable aerosol bin for evaporating
    ! precipitation or ice bins

    ! First, for the rest of the function, get the true nspec, i.e. without rime if evaporating species is ice. No change for precip
    
    IF (iphase == 4) THEN
       nspec = nsp-1
    ELSE
       nspec = nsp
    END IF
    
    mi = getMassIndex(nb,ib,nspec)  ! Wet final mass index
    mi2 = getMassIndex(nb,ib,nspec-1) ! Dry mass index
                           
    ! 1) Find the closest matching bin based on dry particle diameter (a and b bins)
    !    -- note that for ice the index is already corrected for rime. Also, spec%rholiq
    !    -- for all evaporating categories is ok here, since the aerosol densities
    !    -- are always the same
    cd = SUM( mevap(iz,ix,iy,ib:mi2:nb)/spec%rholiq(1:nspec-1) ) / &
              nevap(iz,ix,iy,ib)
    cd = (cd/pi6)**(1./3.) ! Dry diameter
    
    ba=in2a !Ignore 1a and note that "aerobins" contains the lower limit of bin dry diameter
    ba = MIN( MAX(in2a + COUNT( cd > aerobins(ba:fn2a) )-1, in2a), fn2a )

    ! Corresponding b bin is ba+(fn2a-fn1a)
    bb = fn2a + (ba - in2a)
    
    ! 2) Select a or b bin
    IF (a_naerop(iz,ix,iy,bb) <= nlim) THEN
       ! Empty b bin so select a, i.e. do nothing here
       findDry4Wet = ba
    ELSE IF (a_naerop(iz,ix,iy,ba) <= nlim) THEN
       ! Empty a bin so select b
       findDry4Wet = bb
    ELSE
       mi = getMassIndex(nbins,ba,nspec-1) ! Index of the last dry species for current bin in aerosol regime a
       mi2 = getMassIndex(nb,ib,nspec-1)   ! The same for the evaporating particle
       
       ! Both are present - find bin based on compositional similarity
       ra = calc_correlation(a_maerop(iz,ix,iy,ba:mi:nbins),  &
                             mevap(iz,ix,iy,ib:mi2:nb),   &
                             nspec-1                       )
                         
       mi = getMassIndex(nbins,bb,nspec-1) ! Index of the last dry species for current bin in aerosol regime b
       rb = calc_correlation(a_maerop(iz,ix,iy,bb:mi:nbins),  &
                             mevap(iz,ix,iy,ib:mi2:nb),   &
                             nspec-1                       )

       findDry4Wet = ba
       IF (ra<rb) findDry4Wet = bb
    END IF

  END FUNCTION findDry4Wet
  

END MODULE constrain_SALSA
