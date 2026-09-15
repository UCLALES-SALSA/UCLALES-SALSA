MODULE mo_derived_procedures
  USE mo_submctl, ONLY : ica,fca,icb,fcb,ira,fra,          &  
                         iia,fia,in1a,in2b,fn2a,fn2b,      &
                         nbins,ncld,nprc,nice,             &
                         nlim, prlim,                      &
                         spec, pi6,                        & 
                         mean_theta_imm, sigma_theta_imm,  &
                         mean_theta_dep, sigma_theta_dep
  USE util, ONLY : getBinMassArray, getMassIndex
  USE defs, ONLY : cp, alvl
  USE mo_particle_external_properties, ONLY : calcDiamLES
  USE grid, ONLY : nzp,nxp,nyp,level
  USE mo_progn_state, ONLY : a_naerop, a_ncloudp, a_nprecpp, a_nicep,  &
                             a_maerop, a_mcloudp, a_mprecpp, a_micep,  &
                             a_indefaba, a_indefabb, a_indefcba, a_indefcbb, a_indefpba, &
                             a_rp, a_rpp, a_gaerop
  USE mo_diag_state, ONLY : a_rc, a_ri, a_riri, a_srp, a_dn, wt_sfc, wq_sfc, &
                            a_temp, a_press, a_dn
  USE mo_aux_state, ONLY : dzt,dn0
  USE mo_structured_datatypes
  USE math_functions, ONLY : erfm1
  IMPLICIT NONE

  PRIVATE

  ! Procedures for derived output diagnostics currently implemented
  PUBLIC :: bulkNumc,        &  ! Total number concentration of particles of given type (for level >= 4)
            totalWater,      &  ! Total water mixing ratio (makes sense only for level >= 4, since with level <= 3 this is given by rp)
            bulkDiameter,    &  ! Mean diameter of particles of given type (level >= 4)
            bulkMixrat,      &  ! Total mixing ratio of given aerosol constituent (level >= 4)
            binMixrat,       &  ! Binned mixing ratio of given aerosol constituent (level >= 4)
            getBinDiameter,  &  ! Binned particle (wet) diameter (level >= 4)
            binIceDensities, &  ! Get the binned ice particle effective densities (level = 5)
            waterPaths,      &  ! Get waterpaths (lwp,iwp or rpw; level >= 2)
            surfaceFluxes,   &  ! Diagnose surface fluxes in W/m2
            getCDNC,         &  ! Diagnose the "real" CDNC ( 2 um < D < 80 um; level >= 4 )
            getCNC,          &  ! Diagnose the "cloud number concentration" (D > 2 um; level >= 4)
            getReff,         &  ! Get the effective radius using all liquid hydrometeor > 2 um (level >= 4)
            getBinTotMass,   &  ! Get the binned total mass
            getGasConc,      &  ! Get the concentration of specific precursor gas
            initContactAngle,&  ! Convert the IN nucleated fraction into the initial value for contact angle integration
            binSpecMixrat,   &  ! Get binned mass of given aerosol constituent (level >= 4)
            getTerminalVelocity, & ! Calculate the terminal velocity of hydrometeors
            getIceArea,      &  ! Calculate the cross sectional area A=gamma*D**sigma of ice particles   
            getIceAspRatio,  &  ! Calculate the aspect ratio of ice particles 
            getExtinctionCoeffSW, & ! Calculate the extinction coefficient per model layer at 550 nm (you can change wavelength)
            getExtinctionCoeffLW, & ! Calculate the extinction coefficient per model layer at 2100 nm
            getOpticalDepthSW, &    ! Calculate optical depth at 550 nm per model layer as cumulative values viewer at surface
            getOpticalDepthLW       ! Calculate optical depth at 2100 nm per model layer as cumulative values viewer at surface   
  
  CONTAINS

   !
   ! ----------------------------------------------
   ! Subroutine bulkNumc: Calculate the total number
   ! concentration of particles of given type
   !
   ! Juha Tonttila, FMI, 2015
   !
   SUBROUTINE bulkNumc(name,output)
      CHARACTER(LEN=*), INTENT(in) :: name
      REAL, INTENT(out) :: output(nzp,nxp,nyp)
      
      INTEGER :: istr,iend

      istr = 0
      iend = 0

      ! Outputs #/kg
      ! No concentration limits (nlim or prlim) for number

      SELECT CASE(name)
      CASE('Natot','Naa','Nab')
         SELECT CASE(name)
         CASE('Natot')
            istr = in1a
            iend = fn2b
         CASE('Naa')
            istr = in1a
            iend = fn2a
         CASE('Nab')
            istr = in2b
            iend = fn2b
         END SELECT
         output(:,:,:) = SUM(a_naerop%d(:,:,:,istr:iend),DIM=4)
      CASE('Nctot','Nca','Ncb')
         SELECT CASE(name)
         CASE('Nctot')  ! Note: 1a and 2a, 2b combined
            istr = ica%cur
            iend = fcb%cur
         CASE('Nca')
            istr = ica%cur
            iend = fca%cur
         CASE('Ncb')
            istr = icb%cur
            iend = fcb%cur
         END SELECT
         output(:,:,:) = SUM(a_ncloudp%d(:,:,:,istr:iend),DIM=4)
      CASE('Np')
         istr = ira
         iend = fra
         output(:,:,:) = SUM(a_nprecpp%d(:,:,:,istr:iend),DIM=4)
      CASE('Ni')
         istr = iia
         iend = fia
         output(:,:,:) = SUM(a_nicep%d(:,:,:,istr:iend),DIM=4)
      END SELECT

   END SUBROUTINE bulkNumc

   ! --------------------------------------------------------------

   SUBROUTINE totalWater(name,output)
     CHARACTER(len=*), INTENT(in) :: name
     REAL, INTENT(out) :: output(nzp,nxp,nyp)

     output = a_rp%d + a_rc%d + a_srp%d
     IF (level == 5) output = output + a_ri%d + a_riri%d
         
   END SUBROUTINE totalWater
   
   ! --------------------------------------------------------------

   SUBROUTINE bulkDiameter(name,output)
     CHARACTER(len=*), INTENT(in) :: name
     REAL, INTENT(out) :: output(nzp,nxp,nyp)
     INTEGER :: istr,iend
     INTEGER :: nspec
     
     nspec = spec%getNSpec(type="wet")
     
     output = 0.
     
     SELECT CASE(name)
     CASE('Dwatot','Dwaa','Dwab')
        SELECT CASE(name)
        CASE('Dwatot')
           istr = in1a
           iend = fn2b
        CASE('Dwaa')
           istr = in1a
           iend = fn2a
        CASE('Dwab')
           istr = in2b
           iend = fn2b
        END SELECT       
        CALL getMeanDiameter(istr,iend,nbins,nspec,a_naerop,a_maerop,nlim,output,1)
                
     CASE('Dwctot','Dwca','Dwcb')
        SELECT CASE(name)
        CASE('Dwctot')
           istr = ica%cur
           iend = fcb%cur
        CASE('Dwca')
           istr = ica%cur
           iend = fca%cur
        CASE('Dwcb')
           istr = icb%cur
           iend = fcb%cur
        END SELECT        
        CALL getMeanDiameter(istr,iend,ncld,nspec,a_ncloudp,a_mcloudp,nlim,output,2)
        
     CASE('Dwpa')        
        istr = ira
        iend = fra        
        CALL getMeanDiameter(istr,iend,nprc,nspec,a_nprecpp,a_mprecpp,prlim,output,3)
        
     CASE('Dwia')
        istr = iia
        iend = fia         
        CALL getMeanDiameter(istr,iend,nice,nspec+1,a_nicep,a_micep,prlim,output,4)
               
     END SELECT

   END SUBROUTINE bulkDiameter

   !------------------------------------------

   SUBROUTINE getMeanDiameter(zstr,zend,nb,ns,numc,mass,numlim,zdiam,flag)

     IMPLICIT NONE
     
     INTEGER, INTENT(in) :: nb, ns ! Number of bins (nb) and compounds (ns)
     INTEGER, INTENT(in) :: zstr,zend  ! Start and end index for averaging
     TYPE(FloatArray4d), INTENT(in) :: numc
     TYPE(FloatArray4d), INTENT(in) :: mass
     REAL, INTENT(in) :: numlim
     INTEGER, INTENT(IN) :: flag
     REAL, INTENT(out) :: zdiam(nzp,nxp,nyp)
     
     INTEGER :: k,i,j,bin
     REAL :: tot, dwet, tmp(ns)
     REAL :: zlm(nb*ns),zln(nb) ! Local grid point binned mass and number concentrations 
     
     zdiam(:,:,:)=0.
     DO j = 3,nyp-2
        DO i = 3,nxp-2
           DO k = 1,nzp
              zlm(:) = mass%d(k,i,j,:)
              zln(:) = numc%d(k,i,j,:)
              tot=0.
              dwet=0.
              DO bin = zstr,zend                  
                 IF (zln(bin)>numlim) THEN
                    tot=tot+zln(bin)
                    tmp(:) = 0.
                    CALL getBinMassArray(nb,ns,bin,zlm,tmp)
                    dwet=dwet+calcDiamLES(ns,zln(bin),tmp,flag,sph=.FALSE.)*zln(bin)
                 ENDIF
              ENDDO
              IF (tot>numlim) THEN
                 zdiam(k,i,j) = dwet/tot
              ENDIF
           END DO
        END DO
     END DO
     
   END SUBROUTINE getMeanDiameter

   ! -------------------------------------------------

   !
   ! ---------------------------------------------------
   ! SUBROUTINE getBinDiameter
   ! Calculates wet diameter for each bin in the whole domain - this function is for outputs only
   SUBROUTINE getBinDiameter(name, output, nstr, nend)
     USE util, ONLY : getBinMassArray
     USE mo_particle_external_properties, ONLY : calcDiamLES
     USE mo_submctl, ONLY : pi6
     IMPLICIT NONE
     
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr, nend 
     REAL, INTENT(out) :: output(nzp,nxp,nyp,nend-nstr+1)

     INTEGER :: flag, k,i,j,bin, nb, ntot
     INTEGER :: nspec
     TYPE(FloatArray4d), POINTER :: numc
     TYPE(FloatArray4d), POINTER :: mass
     REAL, ALLOCATABLE :: tmp(:), zlm(:), zln(:)
     REAL :: numlim
     
     nspec = spec%getNSpec(type="wet")
     numlim = 0.
     numc => NULL(); mass => NULL(); nb = 0.

     
     SELECT CASE(name)
     CASE('Dwaba')
        flag = 1
        numlim = nlim
        numc => a_naerop
        mass => a_maerop
        nb = nbins
     CASE('Dwabb')
        flag = 1
        numlim = nlim
        numc => a_naerop
        mass => a_maerop
        nb = nbins          
     CASE('Dwcba')
        flag = 2
        numlim = nlim
        numc => a_ncloudp
        mass => a_mcloudp     
        nb = ncld
     CASE('Dwcbb')
        flag = 2
        numlim = nlim
        numc => a_ncloudp
        mass => a_mcloudp
        nb = ncld
     CASE('Dwpba')        
        flag = 3   
        numlim = prlim
        numc => a_nprecpp
        mass => a_mprecpp
        nb = nprc
     CASE('Dwiba')
        flag = 4
        numlim = prlim
        numc => a_nicep
        mass => a_micep 
        nspec = nspec + 1 ! For rime
        nb = nice        
     END SELECT    

     ! Could this be done with pointers or even just straight up? -Juha
     ALLOCATE(tmp(nspec), zlm(nb*nspec), zln(nb))      
      
     output(:,:,:,:)=0.
     DO j = 3,nyp-2
        DO i = 3,nxp-2
           DO k = 1,nzp
              zlm(:) = mass%d(k,i,j,:)
              zln(:) = numc%d(k,i,j,:)
              DO bin = nstr,nend
                 IF (zln(bin)>numlim) THEN
                    tmp(:) = 0.
                    CALL getBinMassArray(nb,nspec,bin,zlm,tmp)
                    ! sph=false enables calculation of non-spherical ice diameter. This will only affect ice.
                    output(k,i,j,bin-nstr+1)=calcDiamLES(nspec,zln(bin),tmp,flag,sph=.FALSE.)  
                 ENDIF
              END DO
           END DO
        END DO
     END DO

     DEALLOCATE(tmp, zlm, zln)
     
   END SUBROUTINE getBinDiameter

   
   !
   ! -----------------------------------
   ! Subroutine bulkMixrat: Find and calculate
   ! the total mixing ratio of a given compound
   ! in aerosol particles or hydrometeors
   !
   ! Juha Tonttila, FMI, 2015
   ! Jaakko Ahola, FMI, 2015
   SUBROUTINE bulkMixrat(name,output)
     CHARACTER(len=*), INTENT(in) :: name
     REAL, INTENT(out) :: output(nzp,nxp,nyp)

     CHARACTER(len=4) :: icomp
     
     INTEGER :: istr,iend,mm,cmax

     iend = 0
     istr = 0
     
     cmax = LEN_TRIM(name)
     output = 0.
     icomp = name(2:cmax-1)
     mm = spec%getIndex(icomp)

     ! Given in kg/kg
     SELECT CASE(name(1:1))
     CASE('a') ! aerosol
        SELECT CASE(name(cmax:cmax))
        CASE('t') ! total
           istr = getMassIndex(nbins,in1a,mm)   
           iend = getMassIndex(nbins,fn2b,mm)    
        CASE('a')
           istr = getMassIndex(nbins,in1a,mm)
           iend = getMassIndex(nbins,fn2a,mm)
        CASE('b')
           istr = getMassIndex(nbins,in2b,mm)
           iend = getMassIndex(nbins,fn2b,mm)
        END SELECT
        output(:,:,:) = SUM(a_maerop%d(:,:,:,istr:iend),DIM=4)

     CASE('c') ! cloud
        SELECT CASE(name(cmax:cmax))
        CASE('t')
           istr = getMassIndex(ncld,ica%cur,mm)
           iend = getMassIndex(ncld,fcb%cur,mm)
        CASE('a')
           istr = getMassIndex(ncld,ica%cur,mm)  
           iend = getMassIndex(ncld,fca%cur,mm)
        CASE('b')
           istr = getMassIndex(ncld,icb%cur,mm)
           iend = getMassIndex(ncld,fcb%cur,mm)
        END SELECT
        output(:,:,:) = SUM(a_mcloudp%d(:,:,:,istr:iend),DIM=4)

     CASE('p')
        istr = getMassIndex(nprc,ira,mm)
        iend = getMassIndex(nprc,fra,mm)
        output(:,:,:) = SUM(a_mprecpp%d(:,:,:,istr:iend),DIM=4)

     CASE('i')
        istr = getMassIndex(nice,iia,mm)
        iend = getMassIndex(nice,fia,mm)
        output(:,:,:) = SUM(a_micep%d(:,:,:,istr:iend),DIM=4)
     END SELECT
     
   END SUBROUTINE bulkMixrat
   !
   
   SUBROUTINE waterPaths(name, output)
     CHARACTER(len=*), INTENT(in) :: name
     REAL, INTENT(out) :: output(nxp,nyp)
     INTEGER :: kk
     
     output = 0.
     
     SELECT CASE(name)
     CASE("lwp")
        DO kk = 1,nzp
           output(:,:) = output(:,:) + a_rc%d(kk,:,:)*a_dn%d(kk,:,:)/dzt%d(kk)   !dzt = 1/dz
        END DO
     CASE("iwp")
        DO kk = 1,nzp
           output(:,:) = output(:,:) +    &
                (a_ri%d(kk,:,:)+a_riri%d(kk,:,:))*a_dn%d(kk,:,:)/dzt%d(kk)
        END DO                
     CASE("rwp")
        IF (level < 4) THEN
           DO kk = 1,nzp
              output(:,:) = output(:,:) + a_rpp%d(kk,:,:)*a_dn%d(kk,:,:)/dzt%d(kk)               
           END DO
        ELSE
           DO kk = 1,nzp
              output(:,:) = output(:,:) + a_srp%d(kk,:,:)*a_dn%d(kk,:,:)/dzt%d(kk)  ! For level 4+ should make some more robust
           END DO                                                             ! since the smallest precipitation bins are not
        END IF                                                                ! actually precipitation.
              
     END SELECT
                      
   END SUBROUTINE waterPaths

   ! -----------------------------

   SUBROUTINE binIceDensities(name,output,nstr,nend)
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr, nend 
     REAL, INTENT(out) :: output(nzp,nxp,nyp,nend-nstr+1)

     INTEGER :: nspec,i,j,k,b
     REAL, ALLOCATABLE :: pmass(:)
     REAL :: diam
     LOGICAL :: sphtype

     output = 0.
     
     nspec = spec%getNSpec(type="total")
     ALLOCATE(pmass(nspec))
     pmass = 0.
     ! sph=false enables calculation of non-spherical ice diameter. This will only affect ice.
     sphtype = .FALSE.
     IF (name == "irhob") sphtype = .TRUE.
     
     DO b = nstr,nend
        DO j = 1,nyp
           DO i = 1,nxp
              DO k = 1,nzp
                 IF (a_nicep%d(k,i,j,b) > prlim) THEN
                    CALL getBinMassArray(nice,nspec,b,a_micep%d(k,i,j,:),pmass)
                    diam = calcDiamLES(nspec,a_nicep%d(k,i,j,b),pmass,4,sph=sphtype)
                    output(k,i,j,b-nstr+1) = SUM(pmass)/a_nicep%d(k,i,j,b)/(pi6*diam**3)
                 END IF
              END DO
           END DO
        END DO
     END DO
        
     DEALLOCATE(pmass)
     
   END SUBROUTINE binIceDensities
   
   ! ------------------------------------------

   SUBROUTINE surfaceFluxes(name,output)
     CHARACTER(len=*), INTENT(in) :: name
     REAL, INTENT(out) :: output(nxp,nyp)
     INTEGER :: kk

     output = 0.
     SELECT CASE(name)
     CASE("lhf")
        output = wq_sfc%d * alvl*(dn0%d(1) + dn0%d(2))*0.5
     CASE("shf")
        output = wt_sfc%d * cp*(dn0%d(1) + dn0%d(2))*0.5
     END SELECT
       
   END SUBROUTINE surfaceFluxes

   ! -------------------------------------------------

   SUBROUTINE getCDNC(name,output)
     CHARACTER(len=*), INTENT(in) :: name
     REAL, INTENT(out) :: output(nzp,nxp,nyp)

     REAL :: Dwcba(nzp,nxp,nyp,fca%cur),Dwcbb(nzp,nxp,nyp,fcb%cur-fca%cur),  &
             Dwpba(nzp,nxp,nyp,nprc)

     REAL, PARAMETER :: lowlim = 2.e-6
     REAL, PARAMETER :: highlim = 80.e-6
     
     Dwcba = 0.; Dwcbb = 0.; Dwpba = 0.
     output = 0.
     
     ! Get binned diameters for clouds and drizzle modes
     CALL getBinDiameter("Dwcba",Dwcba,ica%cur,fca%cur)
     CALL getBinDiameter("Dwcbb",Dwcbb,icb%cur,fcb%cur)
     CALL getBinDiameter("Dwpba",Dwpba,ira,fra)

     output = output + &
              a_dn%d(:,:,:) * SUM( a_ncloudp%d(:,:,:,ica%cur:fca%cur), &
                                   DIM=4, MASK=(Dwcba > lowlim .AND. Dwcba < highlim) )  
     
     output = output + &
              a_dn%d(:,:,:) * SUM( a_ncloudp%d(:,:,:,icb%cur:fcb%cur), &
                                   DIM=4, MASK=(Dwcbb > lowlim .AND. Dwcbb < highlim) )  

     output = output + &
              a_dn%d(:,:,:) * SUM( a_nprecpp%d(:,:,:,ira:fra), &
                                   DIM=4, MASK=(Dwpba > lowlim .AND. Dwpba < highlim) )       
   END SUBROUTINE getCDNC
     
   ! -----------------------------------------------------------

   SUBROUTINE getCNC(name,output)
     CHARACTER(len=*), INTENT(in) :: name
     REAL, INTENT(out) :: output(nzp,nxp,nyp)

     REAL :: Dwcba(nzp,nxp,nyp,fca%cur),Dwcbb(nzp,nxp,nyp,fcb%cur-fca%cur),  &
             Dwpba(nzp,nxp,nyp,nprc)

     REAL, PARAMETER :: lowlim = 2.e-6
     
     Dwcba = 0.; Dwcbb = 0.; Dwpba = 0.
     output = 0.
     
     ! Get binned diameters for clouds and drizzle modes
     CALL getBinDiameter("Dwcba",Dwcba,ica%cur,fca%cur)
     CALL getBinDiameter("Dwcbb",Dwcbb,icb%cur,fcb%cur)
     CALL getBinDiameter("Dwpba",Dwpba,ira,fra)

     output = output +   &
              a_dn%d(:,:,:) * SUM( a_ncloudp%d(:,:,:,ica%cur:fca%cur), &
                                   DIM=4, MASK=(Dwcba > lowlim) )  
     
     output = output +   &
              a_dn%d(:,:,:) * SUM( a_ncloudp%d(:,:,:,icb%cur:fcb%cur), &
                                   DIM=4, MASK=(Dwcbb > lowlim) )  

     output = output +   &
              a_dn%d(:,:,:) * SUM( a_nprecpp%d(:,:,:,ira:fra), &
                                   DIM=4, MASK=(Dwpba > lowlim) )       
   END SUBROUTINE getCNC

   ! -----------------------------------------------------------

   SUBROUTINE getReff(name,output)
     CHARACTER(len=*), INTENT(in) :: name
     REAL, INTENT(out) :: output(nzp,nxp,nyp)

     REAL :: Dwcba(nzp,nxp,nyp,fca%cur),Dwcbb(nzp,nxp,nyp,fcb%cur-fca%cur),  &
             Dwpba(nzp,nxp,nyp,nprc)

     REAL :: third(nzp,nxp,nyp), second(nzp,nxp,nyp) ! Third and second moments of the size distribution
     REAL :: numsum(nzp,nxp,nyp)
     
     REAL, PARAMETER :: lowlim = 2.e-6
     
     Dwcba = 0.; Dwcbb = 0.; Dwpba = 0.
     output = 0.; third=0.; second=0.

     numsum = SUM(a_ncloudp%d,DIM=4) + SUM(a_nprecpp%d,DIM=4)
     
     ! Get binned diameters for clouds and drizzle modes
     CALL getBinDiameter("Dwcba",Dwcba,ica%cur,fca%cur)
     CALL getBinDiameter("Dwcbb",Dwcbb,icb%cur,fcb%cur)
     CALL getBinDiameter("Dwpba",Dwpba,ira,fra)

     third = third + SUM( a_ncloudp%d(:,:,:,ica%cur:fca%cur)*Dwcba**3, DIM=4, MASK=(Dwcba > lowlim) ) 
     third = third + SUM( a_ncloudp%d(:,:,:,icb%cur:fcb%cur)*Dwcbb**3, DIM=4, MASK=(Dwcbb > lowlim) )
     third = third + SUM( a_nprecpp%d(:,:,:,ira:fra)*Dwpba**3, DIM=4, MASK=(Dwpba > lowlim) )

     second = second + SUM( a_ncloudp%d(:,:,:,ica%cur:fca%cur)*Dwcba**2, DIM=4, MASK=(Dwcba > lowlim)  )
     second = second + SUM( a_ncloudp%d(:,:,:,icb%cur:fcb%cur)*Dwcbb**2, DIM=4, MASK=(Dwcbb > lowlim)  )
     second = second + SUM( a_nprecpp%d(:,:,:,ira:fra)*Dwpba**2, DIM=4, MASK=(Dwpba > lowlim)  )
     
     output = 0.5*MERGE(third/MAX(second,1.e-20), 0., numsum > 1.)
     
     
   END SUBROUTINE getReff

   ! -------------------------------------------------------------------------------

   SUBROUTINE getBinTotMass(name,output,nstr,nend)
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr,nend
     REAL, INTENT(out) :: output(nzp,nxp,nyp,nend-nstr+1)
     REAL, POINTER :: parr(:,:,:,:)
     REAL, ALLOCATABLE :: tmp(:)    
     INTEGER :: nspec, i,j,k,bin
     INTEGER :: nb

     parr => NULL()
     output = 0.
     nspec = spec%getNSpec(type="wet") ! Update below if ice
     
     SELECT CASE(name)
     CASE("Maba","Mabb")
        parr => a_maerop%d(:,:,:,:)
        nb = nbins
     CASE("Mcba","Mcbb")
        parr => a_mcloudp%d(:,:,:,:)
        nb = ncld
     CASE("Mpba")
        parr => a_mprecpp%d(:,:,:,:)
        nb = nprc
     CASE("Miba")
        parr => a_micep%d(:,:,:,:)
        nspec = spec%getNSpec(type="total")
        nb = nice
     END SELECT
     
     ALLOCATE(tmp(nspec))
     tmp = 0.
     
     DO bin = nstr,nend
        DO j = 1,nyp
           DO i = 1,nxp
              DO k = 1,nzp
                 CALL getBinMassArray(nb,nspec,bin,parr(k,i,j,:),tmp)
                 output(k,i,j,bin-nstr+1) = SUM(tmp)
              END DO
           END DO
        END DO
     END DO
     
     DEALLOCATE(tmp)
     parr => NULL()
     
   END SUBROUTINE getBinTotMass

! -------------------------------------------------------------------------------

   ! ----------------------------------------------------------------------------
   ! SUBROUTINE GETGASCONC: Gets the sseparate gas concentration for output
   !
   SUBROUTINE getGasConc(name,output)
     CHARACTER(len=*), INTENT(in) :: name
     REAL, INTENT(out) :: output(nzp,nxp,nyp)
     INTEGER :: i
     output = 0. 

     SELECT CASE(name)
     CASE("gSO4")
        i = 1
     CASE("gNO3")
        i = 2
     CASE("gNH4")
        i = 3
     CASE("gOCNV")
        i = 4
     CASE("gOCSV")
        i = 5
     END SELECT

     output(:,:,:) = a_gaerop%d(:,:,:,i)
     
   END SUBROUTINE getGasConc

   ! ----------------------------------------------------------------------------
   ! SUBROUTINE INITCONTACTANGLE: Get the initial value for contact angle
   ! integration in ice nucleation
   !
   SUBROUTINE initContactAngle(name,output,nstr,nend)
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr,nend
     REAL, INTENT(out) :: output(nzp,nxp,nyp,nend-nstr+1)
     REAL, POINTER :: var(:,:,:,:) 
     REAL, PARAMETER :: sqrt2 = SQRT(2.)
     REAL :: mean,sigma
     INTEGER :: bin,i,j,k
     
     var => NULL()

     ! Associate to the target values. Note that some tweaking with the indices is necessary
     ! because of the indices the output routine assumes for binned variables...
     SELECT CASE(name)
     CASE("immThetaaba","depThetaaba")
        var => a_indefaba%d(:,:,:,nstr:nend)
     CASE("immThetaabb","depThetaabb")
        var => a_indefabb%d(:,:,:,1:nend-nstr+1)    
     CASE("immThetacba","depThetacba")
        var => a_indefcba%d(:,:,:,nstr:nend)
     CASE("immThetacbb","depThetacbb")
        var => a_indefcbb%d(:,:,:,1:nend-nstr+1)
     CASE("immThetapba","depThetapba")
        var => a_indefpba%d(:,:,:,nstr:nend)
     END SELECT
        
     IF (name(1:3) == "imm") THEN
        mean = mean_theta_imm
        sigma = sigma_theta_imm
     ELSE IF (name(1:3) == "dep") THEN
        mean = mean_theta_dep
        sigma = sigma_theta_dep
     END IF
         
     DO bin = nstr,nend
        DO j = 1,nyp
           DO i = 1,nxp
              DO k = 1,nzp
                 output(k,i,j,bin-nstr+1) = mean +     &
                      sqrt2*sigma*erfm1( 2.*var(k,i,j,bin-nstr+1) - 1. )
              END DO
           END DO
        END DO
     END DO

     output = MAX(0.,output)
     
     var => NULL()
     
   END SUBROUTINE initContactAngle
   !
   ! ----------------------------------------------
   ! Subroutine binMixrat: Calculate the total dry or wet
   ! Mass concentration for individual bins
   !
   ! Juha Tonttila, FMI, 2015
   ! Tomi Raatikainen, FMI, 2016
   SUBROUTINE binMixrat(ipart,itype,ibin,ii,jj,kk,sumc)
      USE util, ONLY : getBinTotalMass
      USE mo_submctl, ONLY : ncld,nbins,nprc,nice
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(in) :: ipart
      CHARACTER(len=*), INTENT(in) :: itype
      INTEGER, INTENT(in) :: ibin,ii,jj,kk
      REAL, INTENT(out) :: sumc

      CHARACTER(len=20), PARAMETER :: name = "binMixrat"

      INTEGER :: iend

      REAL, POINTER :: tmp(:) => NULL()

      IF (itype == 'dry') THEN
         iend = spec%getNSpec(type="dry")   ! dry CASE
      ELSE IF (itype == 'wet') THEN
         iend = spec%getNSpec(type="wet")   ! wet CASE
         IF (ipart == 'ice') &
              iend = iend+1   ! For ice, take also rime
      ELSE
         STOP 'Error in binMixrat!'
      END IF

      SELECT CASE(ipart)
         CASE('aerosol')
            tmp => a_maerop%d(kk,ii,jj,1:iend*nbins)
            CALL getBinTotalMass(nbins,iend,ibin,tmp,sumc)
         CASE('cloud')
            tmp => a_mcloudp%d(kk,ii,jj,1:iend*ncld)
            CALL getBinTotalMass(ncld,iend,ibin,tmp,sumc)
         CASE('precp')
            tmp => a_mprecpp%d(kk,ii,jj,1:iend*nprc)
            CALL getBinTotalMass(nprc,iend,ibin,tmp,sumc)
         CASE('ice')
            tmp => a_micep%d(kk,ii,jj,1:iend*nice)
            CALL getBinTotalMass(nice,iend,ibin,tmp,sumc) 
         CASE DEFAULT
            STOP 'bin mixrat error'
      END SELECT

      tmp => NULL()
      
   END SUBROUTINE binMixrat


   ! -----------------------------------
   ! Subroutine BinSpecMixrat: Find and calculate
   ! the mass mixing ratio per bin of a given compound
   ! in aerosol particles or hydrometeors
   !
   ! Juha Tonttila, FMI, 2015
   ! Jaakko Ahola, FMI, 2015
   SUBROUTINE binSpecMixrat(name,output,nstr,nend)
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr,nend
     REAL, INTENT(out) :: output(nzp,nxp,nyp,nend-nstr+1)
     CHARACTER(len=4) :: icomp
     
     INTEGER :: istr,iend,mm,cmax

     iend = 0
     istr = 0
     
     cmax = LEN_TRIM(name)
     output = 0.
     icomp = name(3:cmax-1)
     mm = spec%getIndex(icomp)

     ! Given in kg/kg
     SELECT CASE(name(1:2))
     CASE('ma') ! aerosol
        SELECT CASE(name(cmax:cmax))
        CASE('t') ! total
           istr = getMassIndex(nbins,in1a,mm)   
           iend = getMassIndex(nbins,fn2b,mm)    
        CASE('a')
           istr = getMassIndex(nbins,in1a,mm)
           iend = getMassIndex(nbins,fn2a,mm)
        CASE('b')
           istr = getMassIndex(nbins,in2b,mm)
           iend = getMassIndex(nbins,fn2b,mm)
        END SELECT
        output(:,:,:,:) = a_maerop%d(:,:,:,istr:iend)

     CASE('mc') ! cloud
        SELECT CASE(name(cmax:cmax))
        CASE('t')
           istr = getMassIndex(ncld,ica%cur,mm)
           iend = getMassIndex(ncld,fcb%cur,mm)
        CASE('a')
           istr = getMassIndex(ncld,ica%cur,mm)  
           iend = getMassIndex(ncld,fca%cur,mm)
        CASE('b')
           istr = getMassIndex(ncld,icb%cur,mm)
           iend = getMassIndex(ncld,fcb%cur,mm)
        END SELECT
        output(:,:,:,:) = a_mcloudp%d(:,:,:,istr:iend)

     CASE('mp')
        istr = getMassIndex(nprc,ira,mm)
        iend = getMassIndex(nprc,fra,mm)
        output(:,:,:,:) = a_mprecpp%d(:,:,:,istr:iend)

     CASE('mi')
        istr = getMassIndex(nice,iia,mm)
        iend = getMassIndex(nice,fia,mm)
        output(:,:,:,:) = a_micep%d(:,:,:,istr:iend)
     END SELECT

     
   END SUBROUTINE binSpecMixRat
   
 ! ---------------------------------------------------
 ! Silvia Calderon FMI-Kuopio 11.09.2026
 ! All based on derived procedures made by Juha Tonttila
   ! SUBROUTINE getTerminalVelocity(name,output,nstr,nend)
   ! Calculates terminal velocity for hydrometeors for 
   ! each bin in the whole domain - this function is for outputs only
   ! It uses the function terminal_vel(diam,rhop,rhoa,visc,beta,flag,shape,dnsp)
   ! Inside /src/src_shared/
   !
   SUBROUTINE getTerminalVelocity(name,output,nstr,nend)
     USE util, ONLY : getBinMassArray
     USE mo_particle_external_properties, ONLY : calcDiamLES,terminal_vel
     USE mo_submctl, ONLY : pi6,rhowa,pstand,pi
     USE mo_ice_shape, ONLY : t_shape_coeffs, getShapeCoefficients
     IMPLICIT NONE
     
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr, nend 
     REAL, INTENT(out) :: output(nzp,nxp,nyp,nend-nstr+1)

     INTEGER :: flag, k,i,j,bin, nb, ntot
     INTEGER :: nspec
     TYPE(FloatArray4d), POINTER :: numc
     TYPE(FloatArray4d), POINTER :: mass
     TYPE(FloatArray3d), POINTER :: temp  ! ambient temperature [K]
     TYPE(FloatArray3d), POINTER :: pres  ! ambient pressure [Pa]
     TYPE(FloatArray3d), POINTER :: rhoa  ! air density [kg/m3]
     REAL, ALLOCATABLE :: tmp(:), zlm(:), zln(:)
     REAL :: numlim
     REAL :: knud,   &   ! particle knudsen number [1]
             beta,   &   ! Cunningham correction factor [1]
             visc,   &   ! viscosity of air [kg/(m s)]
             mfp,    &   ! mean free path of air molecules [m]
             diam,   &   ! Particle diameters; for ice this is spherical equivalent diameter
             dnsp,   &   ! Nonspherical diameter for ice particles or effective max dimension
             rhop,   &   ! Particle density; For ice this is the effective density (low for non-spherical)
             masspri, massrim ! mass of pristine ice and rimed ice in particle    
     
     TYPE(t_shape_coeffs) :: shape 
     
     nspec = spec%getNSpec(type="wet")
     numlim = 0.
     numc => NULL(); mass => NULL(); nb = 0.
     temp=> a_temp; pres => a_press; 
     rhoa => a_dn! Air density [kg/m3]
     
     !ns = spec%getNSpec(type="total") ! includes rime
     ! iwa = ns-1 irim=ns

     SELECT CASE(name)
     CASE('Vtaba')
        flag = 1
        numlim = nlim
        numc => a_naerop
        mass => a_maerop
        nb = nbins        
     CASE('Vtabb')
        flag = 1
        numlim = nlim
        numc => a_naerop
        mass => a_maerop
        nb = nbins          
     CASE('Vtcba')
        flag = 2
        numlim = nlim
        numc => a_ncloudp
        mass => a_mcloudp     
        nb = ncld
     CASE('Vtcbb')
        flag = 2
        numlim = nlim
        numc => a_ncloudp
        mass => a_mcloudp
        nb = ncld
     CASE('Vtpba')        
        flag = 3   
        numlim = prlim
        numc => a_nprecpp
        mass => a_mprecpp
        nb = nprc
     CASE('Vtiba')
        flag = 4
        numlim = prlim
        numc => a_nicep
        mass => a_micep 
        nspec = nspec + 1 ! For rime
        nb = nice        
     END SELECT    

     ! Could this be done with pointers or even just straight up? -Juha
     ALLOCATE(tmp(nspec), zlm(nb*nspec), zln(nb))      
      
     output(:,:,:,:)=0.
     DO j = 3,nyp-2
        DO i = 3,nxp-2
           DO k = 1,nzp
              zlm(:) = mass%d(k,i,j,:)
              zln(:) = numc%d(k,i,j,:)
              visc = (7.44523e-3*SQRT(temp%d(k,i,j)**3))/(5093.*(temp%d(k,i,j)+110.4))  ! viscosity of air [kg/(m s)] Hinds,p.25 ~ Jacobson FAM eq.4-54
              mfp = (1.656e-10*temp%d(k,i,j)+1.828e-8)*pstand/pres%d(k,i,j) ! mean free path of air [m]
              DO bin = nstr,nend
                 IF (zln(bin)>numlim) THEN
                    tmp(:) = 0.
                    CALL getBinMassArray(nb,nspec,bin,zlm,tmp)
                    ! sph=false enables calculation of non-spherical ice diameter. This will only affect ice.   
                    IF (flag < 4) THEN                
                    	diam=calcDiamLES(nspec,zln(bin),tmp,flag,sph=.TRUE.) 
                    	! Shape coefficients not needed but kept here for consistency 
                    	shape%alpha = pi6*rhowa
      			shape%beta = 3.
      			shape%gamma = pi/4.
      			shape%sigma = 2. 
      			rhop=zlm(bin)/zln(bin)/(pi6*diam**3)  
      			dnsp = diam ! not applicable but just for consistency	
      			knud = 2 * mfp / diam        
                    ELSE 
                        diam=calcDiamLES(nspec,zln(bin),tmp,flag,sph=.TRUE.)  
                    	dnsp=calcDiamLES(nspec,zln(bin),tmp,flag,sph=.FALSE.)  
                    	masspri=SUM(tmp(1:nspec-1))
                    	massrim=tmp(nspec)
                    	CALL getShapeCoefficients(shape,masspri,massrim,zln(bin))
                    	rhop = (masspri*spec%rhoic + massrim*spec%rhori) / &
                               (masspri + massrim)
                        knud = 2 * mfp / dnsp
                    END IF
                    beta = 1.+knud*(1.142+0.558*exp(-0.999/knud))! Cunningham correction factor
		           ! (Allen and Raabe, Aerosol Sci. Tech. 4, 269)
                    output(k,i,j,bin-nstr+1) = terminal_vel(diam,rhop,rhoa%d(k,i,j),visc,beta,flag,shape,dnsp)
                END IF
              END DO
           END DO
        END DO
     END DO

     DEALLOCATE(tmp, zlm, zln)
     
   END SUBROUTINE  getTerminalVelocity   
   
   ! ---------------------------------------------------
 ! Silvia Calderon FMI-Kuopio 11.09.2026
   ! SUBROUTINE getIceArea(name,output,nstr,nend)
   ! Calculates the cross sectional area of ice particles 
   ! A = gamma*D**sigma, gamma and sigma rime-fraction and size dependent
   ! each bin in the whole domain - this function is for outputs only
   ! It uses the same function in cross_sec_area(D,flag,shape) 
   ! Inside /src/src_shared/
   !
   SUBROUTINE getIceArea(name,output,nstr,nend)
     USE util, ONLY : getBinMassArray
     USE mo_particle_external_properties, ONLY : calcDiamLES
     USE mo_ice_shape, ONLY : t_shape_coeffs, getShapeCoefficients
     IMPLICIT NONE
     
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr, nend 
     REAL, INTENT(out) :: output(nzp,nxp,nyp,nend-nstr+1)

     INTEGER :: flag, k,i,j,bin, nb, ntot
     INTEGER :: nspec
     TYPE(FloatArray4d), POINTER :: numc
     TYPE(FloatArray4d), POINTER :: mass
     REAL, ALLOCATABLE :: tmp(:), zlm(:), zln(:)
     REAL :: numlim
     REAL :: diam,   &        ! Particle diameters; for ice this should be the effective max dimension
             masspri, massrim ! mass of pristine ice and rimed ice in particle    
     
     TYPE(t_shape_coeffs) :: shape
     
     nspec = spec%getNSpec(type="total") ! includes rime
     ! iwa = ns-1 irim=ns
     numlim = 0.
     numc => NULL(); mass => NULL(); nb = 0.
       
     flag = 4
     numlim = prlim
     numc => a_nicep
     mass => a_micep 
     nb = nice        
     ALLOCATE(tmp(nspec), zlm(nb*nspec), zln(nb))     
     output(:,:,:,:)=0.
     
     DO j = 3,nyp-2
        DO i = 3,nxp-2
           DO k = 1,nzp
              zlm(:) = mass%d(k,i,j,:)
              zln(:) = numc%d(k,i,j,:)
             DO bin = nstr,nend
                 IF (zln(bin)>numlim) THEN
                    tmp(:) = 0.
                    CALL getBinMassArray(nb,nspec,bin,zlm,tmp)
		    diam=calcDiamLES(nspec,zln(bin),tmp,flag,sph=.FALSE.)  
		    masspri=SUM(tmp(1:nspec-1))
		    massrim=tmp(nspec)
		    CALL getShapeCoefficients(shape,masspri,massrim,zln(bin))
                    output(k,i,j,bin-nstr+1) = shape%gamma*diam**shape%sigma
                END IF
              END DO
           END DO
        END DO
     END DO

     DEALLOCATE(tmp, zlm, zln)
  END SUBROUTINE getIceArea
   
! ---------------------------------------------------
 ! Silvia Calderon FMI-Kuopio 11.09.2026
   ! SUBROUTINE getIceAspRatio(name,output,nstr,nend)
   ! Calculates the aspect ratio of ice particles 
   ! each bin in the whole domain - this function is for outputs only
   ! It applies the same calculation inside 
   ! src/src_salsa/mo_salsa_dynamics.f90
   
   SUBROUTINE getIceAspRatio(name,output,nstr,nend)
     USE util, ONLY : getBinMassArray
     USE mo_particle_external_properties, ONLY : calcDiamLES
     USE mo_submctl, ONLY : pi
     IMPLICIT NONE
     
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr, nend 
     REAL, INTENT(out) :: output(nzp,nxp,nyp,nend-nstr+1)

     INTEGER :: flag, k,i,j,bin, nb, ntot
     INTEGER :: nspec
     TYPE(FloatArray4d), POINTER :: numc
     TYPE(FloatArray4d), POINTER :: mass
     REAL, ALLOCATABLE :: tmp(:), zlm(:), zln(:)
     REAL :: numlim
     REAL :: diam,   &          ! Particle diameters; for ice this should be the effective max dimension
             masspri, massrim,& ! mass of pristine ice and rimed ice in particle    
             rhoice             ! The mean ice density for frozen particles. 
                                ! Takes into account only the bulk ice composition 
     
     nspec = spec%getNSpec(type="total") ! includes rime
     ! iwa = ns-1 irim=ns
     numlim = 0.
     numc => NULL(); mass => NULL(); nb = 0.
            
     flag = 4
     numlim = prlim
     numc => a_nicep
     mass => a_micep 
     nb = nice        
     ALLOCATE(tmp(nspec), zlm(nb*nspec), zln(nb))     
     output(:,:,:,:)=0.
     
     DO j = 3,nyp-2
        DO i = 3,nxp-2
           DO k = 1,nzp
              zlm(:) = mass%d(k,i,j,:)
              zln(:) = numc%d(k,i,j,:)
             DO bin = nstr,nend
                 IF (zln(bin)>numlim) THEN
                    tmp(:) = 0.
                    CALL getBinMassArray(nb,nspec,bin,zlm,tmp)
		    diam=calcDiamLES(nspec,zln(bin),tmp,flag,sph=.FALSE.) 
		    masspri=SUM(tmp(1:nspec-1))
                    massrim=tmp(nspec) 
                    rhoice = (masspri*spec%rhoic + massrim*spec%rhori) / &
                             (masspri + massrim) ! equivalent to pp%rhomean
		    ! Ice particle aspect ratio  
	            ! D for non-spherical ice is defined as the maximum particle length 
	            ! D is related to particle projected area then dnsp~a
		    ! aspect_ratio = c/a = polar radius /equatorial radius
		    ! c should correspond to the volume of the spheroid with eq.radius a
		    ! c = Vspheroid / (4/3*pi*a**2)  Vspheroid = (massice/numc)/rhoice 
		    ! aspect_ratio = 2 * (massice/ice(ii,jj,cc)%numc)/ rhoice / (pi/3*dnsp**3)        
                    output(k,i,j,bin-nstr+1) = 2 * (SUM(tmp)/zln(bin))/ rhoice / (pi/3*diam**3) 
                END IF
              END DO
           END DO
        END DO
     END DO
     DEALLOCATE(tmp, zlm, zln)
   END SUBROUTINE getIceAspRatio
  
   ! ---------------------------------------------------
   ! Silvia Calderon FMI-Kuopio 11.09.2026
   ! All based on derived procedures made by Juha Tonttila
   ! SUBROUTINE getExtinctionCoeffSW(name,output,nstr,nend)
   ! Calculates extinction coefficient of an extinction element
   ! per vertical layer in the whole domain - this function is for outputs only
   ! It uses the calculation approach already employed in 
   ! Inside /src/src_rad/rad_cldwtr.f90 --> aero_rad
   ! bext as SUM(Qext(alpha,r)*pi*r**2*N(r)dr) 
   !
   SUBROUTINE getExtinctionCoeffSW(name,output,nstr,nend)
     USE util, ONLY : getMassIndex,closest, getBinMassArray
     USE mo_salsa_optical_properties, ONLY : aerRefrIBands_SW, &
                                             riReSW, riImSW
     USE mo_submctl, ONLY : pi,pi6,nlim,spec
     USE cldwtr, ONLY : init_aerorad_lookuptables 
    
    IMPLICIT NONE
     
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr, nend 
     REAL, INTENT(out) :: output(nzp,nxp,nyp)

     INTEGER :: flag, k,i,j,bb, nb, ntot, nspec,ss,istr,iend
     TYPE(FloatArray4d), POINTER :: numc
     TYPE(FloatArray4d), POINTER :: mass
     REAL, ALLOCATABLE :: tmp(:), zlm(:), zln(:)
     REAL, ALLOCATABLE :: volc(:,:), voltot(:)
     ! Refractive index for each chemical (closest to current band from tables in submctl)
     REAL, ALLOCATABLE :: refrRe_all(:), refrIm_all(:), volspec(:)
     REAL :: numlim

     CHARACTER(len=34) :: filename
     REAL, ALLOCATABLE, TARGET :: aer_nre_SW(:), aer_nim_SW(:), aer_alpha_SW(:),   &
                       aer_sigma_SW(:,:,:), aer_asym_SW(:,:,:), aer_omega_SW(:,:,:)
     REAL, POINTER :: aer_nre(:) => NULL(), aer_nim(:) => NULL(),          &
                      aer_alpha(:) => NULL(), aer_sigma(:,:,:) => NULL(),  &
                      aer_asym(:,:,:) => NULL(), aer_omega(:,:,:) => NULL()
     
     ! Size parameter for given wavelength and size bin
     REAL :: sizeparam
     ! Selected wavelength 
     REAL :: lambda_r
     ! Volume mean refractive index for single bin
     REAL :: volmean_refrRe, volmean_refrIm 
     ! Lookup table: Indices in refractive index vectors and size parameter vector
     INTEGER :: i_re, i_im, i_alpha
     INTEGER :: refi_ind  ! index for the vector with refractive indices for each wavelength
     
     ! Getting the number of chemical species 
     nspec = spec%getNSpec(type="wet")
     numlim = 0.
     numc => NULL(); mass => NULL(); nb = 0.
     !ns = spec%getNSpec(type="total") ! includes rime
     ! iwa = ns-1 irim=ns

     ! Juha:
     ! Lookup table variables for aerosol optical properties for radiation calculations
     ! Real and imaginary parts of refractive indices, size parameter, extinction crossection, asymmetry parameter and omega
     filename = "datafiles/lut_uclales_salsa_sw.nc"
     CALL init_aerorad_lookuptables(filename, aer_nre_SW, aer_nim_SW, aer_alpha_SW,  &
                                   aer_sigma_SW, aer_asym_SW, aer_omega_SW          )  
     aer_nre => aer_nre_SW(:)
     aer_nim => aer_nim_SW(:)
     aer_alpha => aer_alpha_SW(:)
     aer_sigma => aer_sigma_SW(:,:,:)
     aer_asym => aer_asym_SW(:,:,:)
     aer_omega => aer_omega_SW(:,:,:) 

     SELECT CASE(name)
     CASE('swbextaa')
        flag = 1
        numlim = nlim
        numc => a_naerop
        mass => a_maerop
        nb = nbins 
     CASE('swbextab')
        flag = 1
        numlim = nlim
        numc => a_naerop
        mass => a_maerop
        nb = nbins          
     CASE('swbextca')
        flag = 2
        numlim = nlim
        numc => a_ncloudp
        mass => a_mcloudp     
        nb = ncld
     CASE('swbextcb')
        flag = 2
        numlim = nlim
        numc => a_ncloudp
        mass => a_mcloudp
        nb = ncld
     CASE('swbextpa')        
        flag = 3   
        numlim = prlim
        numc => a_nprecpp
        mass => a_mprecpp
        nb = nprc     
     END SELECT    

     ALLOCATE(refrRe_all(nspec), refrIm_all(nspec), volspec(nspec))
     ALLOCATE(tmp(nb), zlm(nb*nspec), zln(nb))     
     ALLOCATE(volc(nspec,nb))! Corresponding particle volume concentrations for each bin (0 if not used)
     ALLOCATE(voltot(nb))! Total particle volume for each bin  
     
     ! Getting extinction efficiency at 550 nm, representative of visible light
     !  Band:   1:   619.60 Wm^-2, between 50000. and 14500. cm^-1
     !  1 gase(s): and  10 g-points
     !  200 nm - 689.7 nm
     lambda_r = 1/5.5E-05 ! wavenumber in cm-1 for lambda=550nm
     ! Get the refractive indices from the LUT-SW for the current band
     refi_ind = closest(aerRefrIbands_SW,1./lambda_r)
     refrRe_all(:) = riReSW(:,refi_ind)
     refrIm_all(:) = riImSW(:,refi_ind) 
     
     output(:,:,:)=0.
     
     DO j = 3,nyp-2
        DO i = 3,nxp-2
           DO k = 1,nzp
              zlm(:) = mass%d(k,i,j,:)
              zln(:) = numc%d(k,i,j,:)
              tmp(:) = 0.
              ! Loop over chemical species
       	      DO ss = 1,nspec
          	! Mass bin indices
          	istr = getMassIndex(nb,1,ss)
          	iend = getMassIndex(nb,nb,ss) ! ss is index of the aerosol species       
          	!WRITE(*,*) 'Species, n+ki', spec%names(ss), refrRe_all(ss),refrIm_all(ss)
          	! Volumes for each species, 0 if not used or if nothing present
          	IF (flag < 4) THEN                
          	   volc(ss,1:nb) = MERGE( zlm(istr:iend)/spec%rholiq(ss), 0., &
                                        (zln(1:nb)> nlim )            )
                ELSE
                   volc(ss,1:nb) = MERGE( zlm(istr:iend)/spec%rhoice(ss), 0., &
                                        (zln(1:nb)> nlim )            )
                END IF                               
              END DO
       	      voltot(1:nb) = SUM(volc(1:nspec,1:nb),DIM=1)
       	      !
              DO bb = 1,nb!nstr,nend
                 IF (zln(bb)<numlim .OR. voltot(bb) < 1.e-30) CYCLE
                 ! Volume mean refractive indices in current bin
          	 volmean_refrRe = SUM(volc(1:nspec,bb) * refrRe_all(1:nspec)) &
          	                 / voltot(bb)
          	 volmean_refrIm = SUM(volc(1:nspec,bb) * refrIm_all(1:nspec)) &
          	                / voltot(bb)
                 !WRITE(*,*) 'Re+Im', volmean_refrRe, volmean_refrIm
          	 ! Size parameter in current bin sizeparam=alpha=x=pi*D/lambda= pi*D*lambda_r
          	 ! Since lambda_r is in cm-1  volume mean diameter of the bin from m to cm 
          	 sizeparam = 1.e2*lambda_r*pi*(((voltot(bb)/zln(bb))/pi6)**(1./3.))        
          	 ! Corresponding lookup table indices      	 
          	 i_re = closest(aer_nre,volmean_refrRe)
          	 i_im = closest(aer_nim,volmean_refrIm)
          	 i_alpha = closest(aer_alpha,sizeparam)
          	 ! Binned optical properties
          
		 ! LUT tables were build using the size parameter alpha=x=2*pi*radius/lambda as independent variable
		 ! Internally sigma = pi*x**2*Qext and it was already corrected by *1/(4pi**2) but additional
		 ! renormalization is needed because the program is written in x
		 ! aer_sigma(i_re, i_im, i_alpha) still 
		 ! needs to be renormalized by multiplying by lambda**2 = (1./lambda_r)**2
		 ! and transformed from cm to m
		 
		 ! Bin contribution to the extinction coefficient of the model layer 
		 ! bext_aer(kk,bb) = (aer_sigma(i_re,i_im,i_alpha)* 1.0e-4*(1./lambda_r)**2)*naerobin(kk,bb)
		 tmp(bb) = aer_sigma_SW(i_re,i_im,i_alpha)* 1.0e-4*(1./lambda_r)**2*zln(bb)            
              END DO
              output(k,i,j) = SUM(tmp(1:nb))
           END DO
        END DO
     END DO
     
    aer_nre => NULL()
    aer_nim => NULL()
    aer_alpha => NULL()
    aer_sigma => NULL()
    aer_asym => NULL()
    aer_omega => NULL()

     DEALLOCATE(tmp, zlm, zln,volc,voltot,refrRe_all, refrIm_all, volspec)
     
   END SUBROUTINE getExtinctionCoeffSW
   
   ! ---------------------------------------------------
   ! Silvia Calderon FMI-Kuopio 11.09.2026
   ! All based on derived procedures made by Juha Tonttila
   ! SUBROUTINE getExtinctionCoeffSW(name,output,nstr,nend)
   ! Calculates extinction coefficient of an extinction element
   ! per vertical layer in the whole domain - this function is for outputs only
   ! It uses the calculation approach already employed in 
   ! Inside /src/src_rad/rad_cldwtr.f90 --> aero_rad
   ! bext as SUM(Qext(alpha,r)*pi*r**2*N(r)dr) 
   !
   SUBROUTINE getExtinctionCoeffLW(name,output,nstr,nend)
     USE util, ONLY : getMassIndex,closest, getBinMassArray
     USE mo_salsa_optical_properties, ONLY : aerRefrIBands_LW,  &
                                             riReLW, riImLW
     USE mo_submctl, ONLY : pi,pi6,nlim,spec
     USE cldwtr, ONLY : init_aerorad_lookuptables 
    
    IMPLICIT NONE
     
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr, nend 
     REAL, INTENT(out) :: output(nzp,nxp,nyp)

     INTEGER :: flag, k,i,j,bb, nb, ntot, nspec,ss, istr,iend
     TYPE(FloatArray4d), POINTER :: numc
     TYPE(FloatArray4d), POINTER :: mass
     REAL, ALLOCATABLE :: tmp(:), zlm(:), zln(:)
     REAL, ALLOCATABLE :: volc(:,:), voltot(:)
     ! Refractive index for each chemical (closest to current band from tables in submctl)
     REAL,ALLOCATABLE :: refrRe_all(:), refrIm_all(:), volspec(:)   
     REAL :: numlim
     
     CHARACTER(len=34) :: filename
     REAL, ALLOCATABLE, TARGET :: aer_nre_LW(:), aer_nim_LW(:), aer_alpha_LW(:),   &
                       aer_sigma_LW(:,:,:), aer_asym_LW(:,:,:), aer_omega_LW(:,:,:)
     REAL, POINTER :: aer_nre(:) => NULL(), aer_nim(:) => NULL(),          &
                      aer_alpha(:) => NULL(), aer_sigma(:,:,:) => NULL(),  &
                      aer_asym(:,:,:) => NULL(), aer_omega(:,:,:) => NULL()
     
     ! Size parameter for given wavelength and size bin
     REAL :: sizeparam
     ! Selected wavelength 
     REAL :: lambda_r
     ! Volume mean refractive index for single bin
     REAL :: volmean_refrRe, volmean_refrIm 
     ! Lookup table: Indices in refractive index vectors and size parameter vector
     INTEGER :: i_re, i_im, i_alpha
     INTEGER :: refi_ind  ! index for the vector with refractive indices for each wavelength
     
     ! Getting the number of chemical species 
     nspec = spec%getNSpec(type="wet")
     numlim = 0.
     numc => NULL(); mass => NULL(); nb = 0.
     !ns = spec%getNSpec(type="total") ! includes rime
     ! iwa = ns-1 irim=ns

     ! Juha:
     ! Lookup table variables for aerosol optical properties for radiation calculations
     ! Real and imaginary parts of refractive indices, size parameter, extinction crossection, asymmetry parameter and omega
     filename = "datafiles/lut_uclales_salsa_lw.nc"
     CALL init_aerorad_lookuptables(filename, aer_nre_LW, aer_nim_LW, aer_alpha_LW,  &
                                   aer_sigma_LW, aer_asym_LW, aer_omega_LW          )  
     aer_nre => aer_nre_LW(:)
     aer_nim => aer_nim_LW(:)
     aer_alpha => aer_alpha_LW(:)
     aer_sigma => aer_sigma_LW(:,:,:)
     aer_asym => aer_asym_LW(:,:,:)
     aer_omega => aer_omega_LW(:,:,:)
       
     SELECT CASE(name)
     CASE('lwbextaa')
        flag = 1
        numlim = nlim
        numc => a_naerop
        mass => a_maerop
        nb = nbins        
     CASE('lwbextab')
        flag = 1
        numlim = nlim
        numc => a_naerop
        mass => a_maerop
        nb = nbins          
     CASE('lwbextca')
        flag = 2
        numlim = nlim
        numc => a_ncloudp
        mass => a_mcloudp     
        nb = ncld
     CASE('lwbextcb')
        flag = 2
        numlim = nlim
        numc => a_ncloudp
        mass => a_mcloudp
        nb = ncld
     CASE('lwbextpa')        
        flag = 3   
        numlim = prlim
        numc => a_nprecpp
        mass => a_mprecpp
        nb = nprc      
     END SELECT   
                                     
     ALLOCATE(refrRe_all(nspec), refrIm_all(nspec), volspec(nspec))
     ALLOCATE(tmp(nb), zlm(nb*nspec), zln(nb))    
     ALLOCATE(volc(nspec,nb))! Corresponding particle volume concentrations for each bin (0 if not used)
     ALLOCATE(voltot(nb))! Total particle volume for each bin   
     
     ! Getting extinction efficiency in the near IR 
     ! --------------------------------------------------------------------------
     ! You can change the wavelength if needed
     lambda_r = 1/(2100.*1.E-4) ! wavenumber in cm-1 for lambda=2100 nm
     ! Get the refractive indices from the LUT-SW for the current band
     refi_ind = closest(aerRefrIbands_LW,1./lambda_r)
     refrRe_all(:) = riReLW(:,refi_ind)
     refrIm_all(:) = riImLW(:,refi_ind)
     
     output(:,:,:)=0.
     DO j = 3,nyp-2
        DO i = 3,nxp-2
           DO k = 1,nzp
              zlm(:) = mass%d(k,i,j,:)
              zln(:) = numc%d(k,i,j,:)
              tmp(:) = 0.
              ! Loop over chemical species
       	      DO ss = 1,nspec
          	! Mass bin indices
          	istr = getMassIndex(nb,1,ss)
          	iend = getMassIndex(nb,nb,ss) ! ss is index of the aerosol species       
          	!WRITE(*,*) 'Species, n+ki', spec%names(ss), refrRe_all(ss),refrIm_all(ss)
          	! Volumes for each species, 0 if not used or if nothing present
          	IF (flag < 4) THEN                
          	   volc(ss,1:nb) = MERGE( zlm(istr:iend)/spec%rholiq(ss), 0., &
                                        (zln(1:nb)> nlim )            )
                ELSE
                   volc(ss,1:nb) = MERGE( zlm(istr:iend)/spec%rhoice(ss), 0., &
                                        (zln(1:nb)> nlim )            )
                END IF                               
              END DO
       	      voltot(1:nb) = SUM(volc(1:nspec,1:nb),DIM=1)
       	      !
              DO bb = 1, nb ! nstr,nend
                 IF (zln(bb)<numlim .OR. voltot(bb) < 1.e-30) CYCLE
                 ! Volume mean refractive indices in current bin
          	 volmean_refrRe = SUM(volc(1:nspec,bb) * refrRe_all(1:nspec)) &
          	                  / voltot(bb)
          	 volmean_refrIm = SUM(volc(1:nspec,bb) * refrIm_all(1:nspec)) &
          	                  / voltot(bb)
             
          	 ! Size parameter in current bin sizeparam=alpha=x=pi*D/lambda= pi*D*lambda_r
          	 ! Since lambda_r is in cm-1  volume mean diameter of the bin from m to cm 
          	 sizeparam = 1.e2*lambda_r*pi*(((voltot(bb)/zln(bb))/pi6)**(1./3.))        
          	 ! Corresponding lookup table indices
          	 i_re = closest(aer_nre,volmean_refrRe)
          	 i_im = closest(aer_nim,volmean_refrIm)
          	 i_alpha = closest(aer_alpha,sizeparam)
          	 ! Binned optical properties
          
		 ! LUT tables were build using the size parameter alpha=x=2*pi*radius/lambda as independent variable
		 ! Internally sigma = pi*x**2*Qext and it was already corrected by *1/(4pi**2) but additional
		 ! renormalization is needed because the program is written in x
		 ! aer_sigma(i_re, i_im, i_alpha) still 
		 ! needs to be renormalized by multiplying by lambda**2 = (1./lambda_r)**2
		 ! and transformed from cm to m
		  
		 ! Bin contribution to the extinction coefficient of the model layer 
		 ! bext_aer(kk,bb) = (aer_sigma(i_re,i_im,i_alpha)* 1.0e-4*(1./lambda_r)**2)*naerobin(kk,bb)
		 tmp(bb) = aer_sigma_LW(i_re,i_im,i_alpha)* 1.0e-4*(1./lambda_r)**2*zln(bb)              
              END DO
              output(k,i,j) = SUM(tmp(1:nb))
           END DO
        END DO
     END DO

     DEALLOCATE(tmp, zlm, zln, volc,voltot,refrRe_all,refrIm_all,volspec)
     
    aer_nre => NULL()
    aer_nim => NULL()
    aer_alpha => NULL()
    aer_sigma => NULL()
    aer_asym => NULL()
    aer_omega => NULL()
     
   END SUBROUTINE getExtinctionCoeffLW
   
      ! ---------------------------------------------------
   ! Silvia Calderon FMI-Kuopio 11.09.2026
   ! All based on derived procedures made by Juha Tonttila
   ! SUBROUTINE getOpticalDepthSW(name,output,nstr,nend)
   ! Calculates optical depth for the model domain
   ! as a cumulative value along the vertical axis
   ! viewer at surface bottom -> top approach (high -> low)
   ! this function is for outputs only
   ! It uses the calculation approach already employed in 
   ! Inside /src/src_rad/rad_cldwtr.f90 
   ! IMPORTANT: 
   ! For ice particles you must use iodsw and/or iodlw
   ! For total values including gases (i.e. water and ozone) 
   ! you must use todsw and/or todlw
   !
   SUBROUTINE getOpticalDepthSW(name,output,nstr,nend)

     IMPLICIT NONE
     
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr, nend 
     REAL, INTENT(out) :: output(nzp,nxp,nyp)
     REAL :: bext(nzp,nxp,nyp)
     INTEGER :: k,i,j
     
     bext(:,:,:) = 0.
      
     SELECT CASE(name)
      CASE('swAODaa')
        CALL getExtinctionCoeffSW('swbextaa',bext,nstr,nend)
      CASE('swAODab')
        CALL getExtinctionCoeffSW('swbextab',bext,nstr,nend)         
      CASE('swCODca')
        CALL getExtinctionCoeffSW('swbextca',bext,nstr,nend)
      CASE('swCODcb')
        CALL getExtinctionCoeffSW('swbextcb',bext,nstr,nend)       
      CASE('swCODpa')
        CALL getExtinctionCoeffSW('swbextpa',bext,nstr,nend)
     END SELECT 
     
     output(:,:,:)=0.
    
     DO j = 3,nyp-2
        DO i = 3,nxp-2
           output(nzp,i,j) = bext(nzp,i,j)/dzt%d(nzp)   !dzt = 1/dz 
           DO k = nzp-1,1,-1
              output(k,i,j) = output(k,i,j) + bext(k,i,j)/dzt%d(k)   !dzt = 1/dz           
           END DO
        END DO
     END DO
     
   END SUBROUTINE getOpticalDepthSW
   
      ! ---------------------------------------------------
   ! Silvia Calderon FMI-Kuopio 11.09.2026
   ! All based on derived procedures made by Juha Tonttila
   ! SUBROUTINE getOpticalDepthLW(name,output,nstr,nend)
   ! Calculates optical depth for the model domain
   ! as a cumulative value along the vertical axis
   ! viewer at surface bottom -> top approach (high -> low)
   ! this function is for outputs only
   ! It uses the calculation approach already employed in 
   ! Inside /src/src_rad/rad_cldwtr.f90 
   ! IMPORTANT: 
   ! For ice particles you must use iodsw and/or iodlw
   ! For total values including gases (i.e. water and ozone) 
   ! you must use todsw and/or todlw
   !
   SUBROUTINE getOpticalDepthLW(name,output,nstr,nend)
   
     IMPLICIT NONE
     
     CHARACTER(len=*), INTENT(in) :: name
     INTEGER, INTENT(in) :: nstr, nend 
     REAL, INTENT(out) :: output(nzp,nxp,nyp)
     REAL :: bext(nzp,nxp,nyp)
     INTEGER :: k,i,j
     
     bext(:,:,:) = 0.
     
     SELECT CASE(name)
      CASE('lwAODab')
        CALL getExtinctionCoeffLW('lwbextaa',bext,nstr,nend)
      CASE('lwAODbb')
        CALL getExtinctionCoeffLW('lwbextab',bext,nstr,nend)         
      CASE('lwCODca')
        CALL getExtinctionCoeffLW('lwbextca',bext,nstr,nend)
      CASE('lwCODcb')
        CALL getExtinctionCoeffLW('lwbextcb',bext,nstr,nend)       
      CASE('lwCODpa')
        CALL getExtinctionCoeffLW('lwbextpa',bext,nstr,nend)
     END SELECT 
     
     output(:,:,:)=0.
    
     DO j = 3,nyp-2
        DO i = 3,nxp-2
           output(nzp,i,j) = bext(nzp,i,j)/dzt%d(nzp)   !dzt = 1/dz 
           DO k = nzp-1,1,-1
              output(k,i,j) = output(k,i,j) + bext(k,i,j)/dzt%d(k)   !dzt = 1/dz           
           END DO
        END DO
     END DO
     
   END SUBROUTINE getOpticalDepthLW
  
END MODULE mo_derived_procedures
