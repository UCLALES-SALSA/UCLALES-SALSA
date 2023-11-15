MODULE mo_salsa_secondary_ice

 
  CONTAINS

   SUBROUTINE secondary_ice(kbdim,kproma,klev,ppres,ptemp,ptstep)
     USE mo_submctl, ONLY : lssecice, lssiprimespln, lssipdropfrac, &
          lssipicecollbreak, spec
      USE mo_salsa_SIP_DF, ONLY: dropfracturing
      USE mo_salsa_SIP_RS, ONLY: rimesplintering
      USE mo_salsa_SIP_IIBR, ONLY: iceicecollbreak
      
      INTEGER, INTENT(in) :: kbdim,kproma,klev
      REAL, INTENT(in)    :: ppres(kbdim,klev), ptemp(kbdim,klev)
      REAL, INTENT(in)    :: ptstep 

      INTEGER :: nspec 

      ! Secondary ice processes.
      ! These need information about the ice collected liquid drops; Therefore it is imperative, that
      ! the bin redistribution routine is NOT called between coagulation and the secondary ice
      ! parameterizations.
      IF (lssecice%state) THEN

         nspec = spec%getNSpec(type='total')
         
         ! Rime splintering
         IF (lssiprimespln%state) &
             CALL rimesplintering(kbdim,kproma,klev,nspec,ppres,ptemp,ptstep)
         ! Drop fracturing
         IF (lssipdropfrac%state) &
            CALL dropfracturing(kbdim,kproma,klev,nspec,ppres,ptemp,ptstep)
         ! Ice-ice collisional breakup
         IF (lssipicecollbreak%state) &
            CALL iceicecollbreak(kbdim,kproma,klev,nspec,ppres,ptemp,ptstep)
         
      END IF
         
   END SUBROUTINE secondary_ice 

    
END MODULE mo_salsa_secondary_ice
