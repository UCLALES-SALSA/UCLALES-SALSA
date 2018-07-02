MODULE constrain_SALSA
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

      USE grid, ONLY : a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,   &
                       a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,   &
                       a_nicep,  a_nicet, a_nsnowp, a_nsnowt,                            & ! ice'n'snow
                       a_micep,  a_micet, a_msnowp, a_msnowt,                            & ! ice'n'snow
                       dtlt, nxp,nyp,nzp, level
      USE mo_submctl, ONLY : nbins, ncld, nprc, &
                             nice, nsnw          !ice'n'snow

      INTEGER, INTENT(in) :: nn

      INTEGER :: cc, ii,jj,kk,ni

      DO jj = 3, nyp-2

         DO ii = 3, nxp-2

            DO kk = 1, nzp

               ! Aerosols
               DO cc = 1, nbins

                  IF ( a_naerop(kk,ii,jj,cc)+a_naerot(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_naerot(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_naerop(kk,ii,jj,cc))/dtlt,a_naerot(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_maerot(kk,ii,jj,(ni-1)*nbins+cc) = MAX( ((1.e-10-1.0)*a_maerop(kk,ii,jj,(ni-1)*nbins+cc))/dtlt,  &
                                                                 a_maerot(kk,ii,jj,(ni-1)*nbins+cc) )
                     END DO

                  END IF

               END DO

               ! Cloud droplets
               DO cc = 1, ncld

                  IF ( a_ncloudp(kk,ii,jj,cc)+a_ncloudt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_ncloudt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_ncloudp(kk,ii,jj,cc))/dtlt,a_ncloudt(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) = MAX( ((1.e-10-1.0)*a_mcloudp(kk,ii,jj,(ni-1)*ncld+cc))/dtlt,  &
                                                                 a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) )
                     END DO

                  END IF

               END DO

               ! Precipitation
               DO cc = 1, nprc

                  IF ( a_nprecpp(kk,ii,jj,cc)+a_nprecpt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nprecpt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nprecpp(kk,ii,jj,cc))/dtlt,a_nprecpt(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) = MAX( ((1.e-10-1.0)*a_mprecpp(kk,ii,jj,(ni-1)*nprc+cc))/dtlt,  &
                                                                 a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) )
                     END DO

                  END IF

               END DO

             ! ice particles
             IF (level<5) CYCLE
             DO cc = 1,nice

                  IF ( a_nicep(kk,ii,jj,cc)+a_nicet(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nicet(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nicep(kk,ii,jj,cc))/dtlt,a_nicet(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_micet(kk,ii,jj,(ni-1)*ncld+cc) = MAX( ((1.e-10-1.0)*a_micep(kk,ii,jj,(ni-1)*nice+cc))/dtlt,  &
                                                               a_micet(kk,ii,jj,(ni-1)*nice+cc) )
                     END DO

                  END IF

               END DO

               ! Snow
               DO cc = 1, nsnw

                  IF ( a_nsnowp(kk,ii,jj,cc)+a_nsnowt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nsnowt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nsnowp(kk,ii,jj,cc))/dtlt,a_nsnowt(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_msnowt(kk,ii,jj,(ni-1)*nsnw+cc) = MAX( ((1.e-10-1.0)*a_msnowp(kk,ii,jj,(ni-1)*nsnw+cc))/dtlt,  &
                                                                a_msnowt(kk,ii,jj,(ni-1)*nsnw+cc) )
                     END DO

                  END IF

               END DO

            END DO ! kk

         END DO ! ii

      END DO ! jj

   END SUBROUTINE tend_constrain



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

   SUBROUTINE SALSA_diagnostics
      USE grid, ONLY : nxp,nyp,nzp,    &
                       a_naerop,a_maerop,a_ncloudp,a_mcloudp,a_nprecpp,a_mprecpp,a_gaerop, &
                       a_rc, a_srp,a_snrp,  &
                       a_rh, a_temp, a_ri,a_srs,a_snrs,a_rhi,                                      &
                       a_nicep,a_micep,a_nsnowp,a_msnowp, level,   &
                       binMixrat
      USE mo_submctl, ONLY : nbins,ncld,nprc,ica,fca,icb,fcb,ira,fra,              &
                             in1a,in2a,fn2a,fn2b,                        &
                             nice,nsnw,iia,fia,iib,fib,isa,fsa,                    &
                             spec, surfw0, rg, nlim, prlim, pi, &
                             lscndgas, pi6, avog,                                  &
                             aerobins
      USE util, ONLY : getMassIndex, calc_correlation
      IMPLICIT NONE

      INTEGER :: i,j,k,bb,bc,ba,s,sc,sa,str,end,nc,nn

      REAL :: cd

      REAL :: ra, rb

    REAL :: zvol
    REAL :: zdh2o
    REAL :: ns, zbb, zaa ! Number of moles, Raoult effect, Kelvin effect; For calculating the critical radius
    REAL :: cdcld,cdprc ! Critical diameter for cloud droplets and precipitation
    REAL :: zdrms, zwams

    ! Remove negative values
    a_naerop = MAX(0.,a_naerop)
    a_ncloudp = MAX(0.,a_ncloudp)
    a_nprecpp = MAX(0.,a_nprecpp)
    a_maerop = MAX(0.,a_maerop)
    a_mcloudp = MAX(0.,a_mcloudp)
    a_mprecpp = MAX(0.,a_mprecpp)
    
    a_nicep = MAX(0.,a_nicep)
    a_nsnowp = MAX(0.,a_nsnowp)
    a_micep = MAX(0.,a_micep)
    a_msnowp = MAX(0.,a_msnowp)
    
    nn = spec%getNSpec() ! total number of species
    
    ! Remove particles that have number but not mass
    DO j = 3,nyp-2
       DO i = 3,nxp-2
          DO k = 1,nzp
             ! Aerosols
             DO bc = 1,nbins    
                IF (a_naerop(k,i,j,bc) > 0. .AND. SUM(a_maerop(k,i,j,bc:getMassIndex(nbins,bc,nn):nbins)) <= 0.) THEN
                   a_naerop(k,i,j,bc) = 0.
                   a_maerop(k,i,j,bc:getMassIndex(nbins,bc,nn):nbins) = 0.
                END IF
             END DO

             ! Clouds
             DO bc = 1,ncld
                IF (a_ncloudp(k,i,j,bc) > 0. .AND. SUM(a_mcloudp(k,i,j,bc:getMassIndex(ncld,bc,nn):ncld)) <= 0.) THEN
                   a_ncloudp(k,i,j,bc) = 0.
                   a_mcloudp(k,i,j,bc:getMassIndex(ncld,bc,nn):ncld) = 0.
                END IF
             END DO ! ncld

             ! Precipitation
             DO bc = 1,nprc
                IF (a_nprecpp(k,i,j,bc) > 0. .AND. a_mprecpp(k,i,j,getMassIndex(nprc,bc,nn)) <= 0.) THEN
                   a_nprecpp(k,i,j,bc) = 0.
                   a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn):nprc) = 0.
                END IF
             END DO ! nprc

             ! Ice
             IF (level<5) CYCLE
             DO bc = 1,nice
                IF (a_nicep(k,i,j,bc) > 0. .AND. SUM(a_micep(k,i,j,bc:getMassIndex(nice,bc,nn):nice)) <= 0.) THEN
                   a_nicep(k,i,j,bc) = 0.
                   a_micep(k,i,j,bc:getMassIndex(nice,bc,nn):nice) = 0.
                END IF
             END DO ! ncld

             ! Snow
             DO bc = 1,nsnw
                IF (a_nsnowp(k,i,j,bc) > 0. .AND. a_msnowp(k,i,j,getMassIndex(nsnw,bc,nn)) <= 0.) THEN
                   a_nsnowp(k,i,j,bc) = 0.
                   a_msnowp(k,i,j,bc:getMassIndex(nsnw,bc,nn):nsnw) = 0.
                END IF
             END DO ! nsnw

          END DO !k
       END DO !i
    END DO !j

   ! Ghost species
    DO j = 3,nyp-2
       DO i = 3,nxp-2
          DO k = 1,nzp

             ! Loop over cloud droplet bins
             DO bc = ica%cur,fcb%cur

                IF ( a_ncloudp(k,i,j,bc) > nlim .AND. a_rh(k,i,j)<0.999 .AND.   &
                     a_mcloudp(k,i,j,getMassIndex(ncld,bc,nn)) < 1.e-5  ) THEN

                   ! Critical diameter
                   ns = SUM( spec%diss(1:nn-1)*a_mcloudp(k,i,j,bc:getMassIndex(ncld,bc,nn-1):ncld)/spec%MM(1:nn-1) ) / &
                        a_ncloudp(k,i,j,bc)
                   zbb = 6.*spec%mwa*ns/(pi*spec%rhowa)
                   zaa = 4.*spec%mwa*surfw0/(rg*spec%rhowa*a_temp(k,i,j))
                   cdcld = SQRT(3.*zbb/zaa)

                   ! Wet diameter
                   zvol = SUM( a_mcloudp(k,i,j,bc:getMassIndex(ncld,bc,nn):ncld)/spec%rholiq(1:nn) )/a_ncloudp(k,i,j,bc)
                   zdh2o = (zvol/pi6)**(1./3.)

                   ! Lose the droplets if smaller than 0.2*(critical size) or 2 um
                   IF ( zdh2o < MAX(0.2*cdcld,2.e-6) .OR.     &
                        a_mcloudp(k,i,j,getMassIndex(ncld,bc,nn)) < 1.e-25*a_ncloudp(k,i,j,bc) )  THEN

                      IF (bc<=fca%cur) THEN
                          ba = ica%par + (bc-ica%cur) ! Index for parallel aerosol bin
                      ELSE
                          ba = icb%par + (bc-icb%cur) ! Index for parallel aerosol bin
                      END IF
                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_ncloudp(k,i,j,bc)
                      a_ncloudp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = getMassIndex(ncld,bc,s)
                         sa = getMassIndex(nbins,ba,s)
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mcloudp(k,i,j,sc)
                         a_mcloudp(k,i,j,sc) = 0.
                      END DO

                   END IF ! critical diameter

                END IF  ! blim

             END DO ! bc

             ! Loop over precipitation bins
             DO bc = ira,fra

                IF ( a_nprecpp(k,i,j,bc) > prlim .AND. a_rh(k,i,j)<0.999 .AND.  &
                     a_mprecpp(k,i,j,getMassIndex(nprc,bc,nn)) < 1.e-6 ) THEN

                   ! Critical radius
                   ns = SUM( spec%diss(1:nn-1)*a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn-1):nprc)/spec%MM(1:nn-1) ) / &
                        a_nprecpp(k,i,j,bc)
                   zbb = 6.*spec%mwa*ns/(pi*spec%rhowa)
                   zaa = 4.*spec%mwa*surfw0/(rg*spec%rhowa*a_temp(k,i,j))
                   cdprc = SQRT(3.*zbb/zaa)

                   ! Wet radius
                   zvol = SUM( a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn):nprc)/spec%rholiq(1:nn) )/a_nprecpp(k,i,j,bc)
                   zdh2o = (zvol/pi6)**(1./3.)

                   ! Lose the droplets if smaller than 0.02*critical radius or 2 um
                   IF ( zdh2o < MAX(0.02*cdprc,2.e-6) .OR.   &
                        a_mprecpp(k,i,j,getMassIndex(nprc,bc,nn))<1e-25*a_nprecpp(k,i,j,bc) ) THEN
                      ! Move evaporating rain drops to a soluble aerosol bin with
                      ! the closest match in dry particle mass. Ain't perfect but
                      ! the bin update subroutine in SALSA will take care of the rest.
                      zvol = SUM( a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn-1):nprc) )/a_nprecpp(k,i,j,bc) ! Dry mass

                      ! 1) Find the closest matching bin based on dry particle diameter (a and b bins)
                      cd = SUM( a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn-1):nprc)/spec%rholiq(1:nn-1) ) / &
                           a_nprecpp(k,i,j,bc) ! Dry radius
                      cd = (cd/pi6)**(1./3.)
                      ba=in2a !Ignore 1a and note that "aerobins" contains the lower limit of bin dry diameter
                      DO WHILE (cd>=aerobins(ba+1) .AND. ba<fn2a)
                         ba=ba+1
                      END DO
                      ! Corresponding b bin is ba+(fn2a-fn1a)=ba+fn2a-(in2a-1)=ba+fn2a-in2a+1
                      bb=ba+fn2a-in2a+1
                      ! 2) Select a or b bin
                      IF (a_naerop(k,i,j,bb) <= nlim) THEN
                         ! Empty b bin so select a
                         !ba = ba
                      ELSE IF (a_naerop(k,i,j,ba) <= nlim) THEN
                         ! Empty a bin so select b
                         ba = bb
                      ELSE
                         ! Both are present - find bin based on compositional similarity
                         ra = calc_correlation(a_maerop(k,i,j,ba:getMassIndex(nbins,ba,nn-1):nbins),  &
                                               a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn-1):nprc),   &
                                               nn-1)

                         rb = calc_correlation(a_maerop(k,i,j,bb:getMassIndex(nbins,bb,nn-1):nbins),  &
                                               a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn-1):nprc),   &
                                               nn-1)
                         IF (ra<rb) ba = bb
                      END IF                     

                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nprecpp(k,i,j,bc)
                      a_nprecpp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = getMassIndex(nprc,bc,s)
                         sa = getMassIndex(nbins,ba,s)
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mprecpp(k,i,j,sc)
                         a_mprecpp(k,i,j,sc) = 0.
                      END DO

                   END IF ! Critical radius

                END IF ! prlim

             END DO ! bc

             ! Loop over ice bins
             DO bc = iia%cur,fib%cur

                IF ( a_nicep(k,i,j,bc) > prlim .AND. a_rhi(k,i,j)<0.999 .AND.   &
                     a_micep(k,i,j,getMassIndex(nice,bc,nn)) < 1.e-15 ) THEN

                   ! Ice and snow don't have a critical size, but lose particles when water content becomes low enough or size is < 2e-6m
                   
                   ! Diameter
                   cd = SUM( a_micep(k,i,j,bc:getMassIndex(nice,bc,nn):nice)/spec%rhoice(1:nn) ) / &
                        a_nicep(k,i,j,bc)
                   cd = (cd/pi6)**(1./3.)

                   ! Lose ice when dry to total mass ratio is more than 0.5
                   CALL binMixrat("ice","dry",bc,i,j,k,zdrms)
                   CALL binMixrat("ice","wet",bc,i,j,k,zwams)
                   zvol = zdrms/zwams

                   IF ( zvol>0.5 .OR. cd < 2.e-6 ) THEN
                      IF (bc<=fia%cur) THEN
                         ba = iia%par + (bc-iia%cur) ! Index for parallel aerosol bin
                      ELSE
                         ba = iib%par + (bc-iib%cur) ! Index for parallel aerosol bin
                      END IF

                      ! Move the number of particles from ice to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nicep(k,i,j,bc)
                      a_nicep(k,i,j,bc) = 0.

                      ! Move mass material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = getMassIndex(nice,bc,s)
                         sa = getMassIndex(nbins,ba,s)
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_micep(k,i,j,sc)
                         a_micep(k,i,j,sc) = 0.
                      END DO
                   END IF

                END IF  ! prlim

             END DO ! bc

             ! Loop over snow bins
             DO bc = isa,fsa

                IF ( a_nsnowp(k,i,j,bc) > prlim .AND. a_rhi(k,i,j)<0.999 .AND.  &
                     a_msnowp(k,i,j,getMassIndex(nsnw,bc,nn)) < 1.e-20 ) THEN

                   cd = SUM( a_msnowp(k,i,j,bc:getMassIndex(nsnw,bc,nn):nsnw)/spec%rhosnow(1:nn) ) / &
                        a_nsnowp(k,i,j,bc)
                   cd = (cd/pi6)**(1./3.)

                   ! Lose snow when dry to total mass ratio is more than 0.5 or diameter < 2.e-6m
                   CALL binMixrat("snow","dry",bc,i,j,k,zdrms)
                   CALL binMixrat("snow","wet",bc,i,j,k,zwams)
                   zvol = zdrms/zwams

                   IF ( zvol>0.5 .OR. cd < 2.e-6) THEN
 
                      ! Move evaporating snow to aerosol bin based on dry radius and chemical composition

                      ! 1) Find the closest matching bin based on dry particle diameter (a and b bins)
                      cd = SUM( a_msnowp(k,i,j,bc:getMassIndex(nsnw,bc,nn-1):nsnw)/spec%rhosnow(1:nn-1) ) / &
                           a_nsnowp(k,i,j,bc)
                      cd = (cd/pi6)**(1./3.)
                      
                      ba=in2a ! Ignore 1a and note that aerobins contains the lower limit of bin dry radius
                      DO WHILE (cd>=aerobins(ba+1) .AND. ba<fn2a)
                         ba=ba+1
                      END DO
                      ! Corresponding b bin is ba+(fn2a-fn1a)=ba+fn2a-(in2a-1)=ba+fn2a-in2a+1
                      bb=ba+fn2a-in2a+1
                      ! 2) Select a or b bin
                      IF (a_naerop(k,i,j,bb) <= nlim) THEN
                         ! Empty b bin so select a
                         !ba = ba
                      ELSE IF (a_naerop(k,i,j,ba) <= nlim) THEN
                         ! Empty a bin so select b
                         ba = bb
                      ELSE
                         ! Both are present - find bin based on compositional similarity
                         ra = calc_correlation(a_maerop(k,i,j,ba:getMassIndex(nbins,ba,nn-1):nbins),  &
                                               a_msnowp(k,i,j,bc:getMassIndex(nsnw,bc,nn-1):nsnw),    &
                                               nn-1)
                         rb = calc_correlation(a_maerop(k,i,j,bb:getMassIndex(nbins,bb,nn-1):nbins),  &
                                               a_msnowp(k,i,j,bc:getMassIndex(nsnw,bc,nn-1):nsnw),    &
                                               nn-1)
                         IF (ra<rb) ba = bb
                      END IF


                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nsnowp(k,i,j,bc)
                      a_nsnowp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = getMassIndex(nsnw,bc,s)
                         sa = getMassIndex(nbins,ba,s)
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_msnowp(k,i,j,sc)
                         a_msnowp(k,i,j,sc) = 0.
                      END DO
                   END IF

                END IF ! prlim

             END DO ! bc

             ! Loop over aerosol bins
             DO ba = 1,nbins
                IF (a_naerop(k,i,j,ba) > nlim) THEN
                   zvol = SUM( a_maerop(k,i,j,ba:getMassIndex(nbins,ba,nn-1):nbins)/spec%rholiq(1:nn-1) )/ &
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
                      a_maerop(k,i,j,ba:getMassIndex(nbins,ba,nn):nbins) = 0.
                      a_naerop(k,i,j,ba) = 0.
                   END IF
                END IF
             END DO

          END DO   ! k
       END DO   ! i
    END DO   ! j

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

      ! ice, regimes a and b
      str = getMassIndex(nice,iia%cur,nc)
      end = getMassIndex(nice,fib%cur,nc)
      a_ri(:,:,:) = SUM(a_micep(:,:,:,str:end),DIM=4)
      ! Snow
      str = getMassIndex(nsnw,isa,nc)
      end = getMassIndex(nsnw,fsa,nc)
      a_srs(:,:,:) = SUM(a_msnowp(:,:,:,str:end),DIM=4)
      a_snrs(:,:,:) = SUM(a_nsnowp(:,:,:,isa:fsa),DIM=4)

   END SUBROUTINE SALSA_diagnostics



END MODULE constrain_SALSA
