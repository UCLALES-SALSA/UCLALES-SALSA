MODULE mo_derived_procedures
  USE mo_submctl, ONLY : ica,fca,icb,fcb,ira,fra,     & 
                         iia,fia,in1a,in2b,fn2a,fn2b, &
                         nbins,ncld,nprc,nice,        &
                         nlim, prlim,                 &
                         spec, pi6
  USE util, ONLY : getBinMassArray, getMassIndex
  USE mo_particle_external_properties, ONLY : calcDiamLES
  USE grid, ONLY : nzp,nxp,nyp,level
  USE mo_progn_state, ONLY : a_naerop, a_ncloudp, a_nprecpp, a_nicep,  &
                             a_maerop, a_mcloudp, a_mprecpp, a_micep,  &
                             a_rp
  USE mo_diag_state, ONLY : a_rc, a_ri, a_riri, a_srp
  USE mo_structured_datatypes
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: bulkNumc, totalWater, bulkDiameter, bulkMixrat
    
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
      CASE('precp')
         istr = ira
         iend = fra
         output(:,:,:) = SUM(a_nprecpp%d(:,:,:,istr:iend),DIM=4)
      CASE('ice')
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
         
   END SUBROUTINE
   
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
        CALL getRadius(istr,iend,nbins,nspec,a_naerop,a_maerop,nlim,output,1)
                
     CASE('Dwctot','Dwca','Dwcb')
        SELECT CASE(name)
        CASE('Dwctot')
           istr = ica%cur
           iend = fcb%cur
        CASE('Dwca')
           istr = ica%cur
           iend = fca%cur
        CASE('dwcb')
           istr = icb%cur
           iend = fcb%cur
        END SELECT        
        CALL getRadius(istr,iend,ncld,nspec,a_ncloudp,a_mcloudp,nlim,output,2)
        
     CASE('precp')        
        istr = ira
        iend = fra        
        CALL getRadius(istr,iend,nprc,nspec,a_nprecpp,a_mprecpp,prlim,output,3)
        
     CASE('ice')
        istr = iia
        iend = fia         
        CALL getRadius(istr,iend,nice,nspec+1,a_nicep,a_micep,prlim,output,4)
               
     END SELECT

   END SUBROUTINE bulkDiameter

   !------------------------------------------

   SUBROUTINE getRadius(zstr,zend,nb,ns,numc,mass,numlim,zdiam,flag)

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
                    dwet=dwet+calcDiamLES(ns,zln(bin),tmp,flag)*zln(bin)
                 ENDIF
              ENDDO
              IF (tot>numlim) THEN
                 zdiam(k,i,j) = dwet/tot
              ENDIF
           END DO
        END DO
     END DO
     
   END SUBROUTINE getRadius

   ! -------------------------------------------------

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
   

   
   
END MODULE mo_derived_procedures
