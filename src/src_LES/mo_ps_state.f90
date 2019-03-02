MODULE mo_ps_state
  USE classFieldArray
  USE mo_structured_datatypes
  USE mo_ps_procedures
  USE mo_check_state, ONLY : checkOutputs
  IMPLICIT NONE

  ! Variables for profile statistics output.
  ! Note that for some of these variables (so far all) no data is stored, but the onDemand function is associated to get
  ! the necessary diagnostics when writing output.

  TYPE(FloatArray1d), TARGET :: ps_theta, ps_temp, ps_press,     &
                                ps_rp, ps_rc, ps_rr, ps_ri,      &
                                ps_riri, ps_Naa, ps_Nab,         &
                                ps_Nca, ps_Ncb, ps_Np, ps_Ni,    &
                                ps_RH, ps_rsl, ps_RHI, ps_rsi,   &

                                ps_rrate, ps_irate,              &

                                ps_rflx, ps_sflx, ps_lwup,       &
                                ps_lwdn, ps_swup, ps_swdn,       &

                                ps_Dwaa, ps_Dwab, ps_Dwca,       &
                                ps_Dwcb, ps_Dwpa, ps_Dwia,       &

                                ps_aSO4a, ps_aSO4b, ps_cSO4a, ps_cSO4b, ps_pSO4a, ps_iSO4a,  &
                                ps_aOCa, ps_aOCb, ps_cOCa, ps_cOCb, ps_pOCa, ps_iOCa,  &
                                ps_aBCa, ps_aBCb, ps_cBCa, ps_cBCb, ps_pBCa, ps_iBCa,  &
                                ps_aDUa, ps_aDUb, ps_cDUa, ps_cDUb, ps_pDUa, ps_iDUa,  &
                                ps_aSSa, ps_aSSb, ps_cSSa, ps_cSSb, ps_pSSa, ps_iSSa,  &
                                ps_aNOa, ps_aNOb, ps_cNOa, ps_cNOb, ps_pNOa, ps_iNOa,  &
                                ps_aNHa, ps_aNHb, ps_cNHa, ps_cNHb, ps_pNHa, ps_iNHa

 
  CONTAINS

    SUBROUTINE setPSVariables(PS,outputlist,level)
      TYPE(FieldArray), INTENT(inout) :: PS
      CHARACTER(len=*), INTENT(in)   :: outputlist(:)
      INTEGER, INTENT(in)             :: level
      CLASS(*), POINTER :: pipeline => NULL()
      

      IF (ANY(outputlist == "theta")) THEN
         pipeline => NULL()
         ps_theta = FloatArray1d()
         ps_theta%onDemand => globalAvgProfile
         pipeline => ps_theta
         CALL PS%newField("theta", "Potential temperature", "K", "ztt",   &
                          ANY(outputlist == "theta"), pipeline)
      END IF

      IF (ANY(outputlist == "temp")) THEN
         pipeline => NULL()
         ps_temp = FloatArray1d()
         ps_temp%onDemand => globalAvgProfile
         pipeline => ps_temp
         CALL PS%newField("temp", "Abs temperature", "K", "ztt",   &
                          ANY(outputlist == "temp"), pipeline)
      END IF

      IF (ANY(outputlist == "press")) THEN
         pipeline => NULL()
         ps_press = FloatArray1d()
         ps_press%onDemand => globalAvgProfile
         pipeline => ps_press
         CALL PS%newField("press", "Pressure", "Pa", "ztt",   &
                          ANY(outputlist == "press"), pipeline)
      END IF

      IF (ANY(outputlist == "rp")) THEN
         pipeline => NULL()
         ps_rp = FloatArray1d()
         ps_rp%onDemand => globalAvgProfile
         pipeline => ps_rp
         CALL PS%newField("rp", "Water vapor mixing ratio", "kg/kg", "ztt",   &
                          ANY(outputlist == "rp"), pipeline)
      END IF

      IF (ANY(outputlist == "rc")) THEN
         pipeline => NULL()
         ps_rc = FloatArray1d()
         ps_rc%onDemand => globalAvgProfile
         pipeline => ps_rc
         CALL PS%newField("rc", "Cloud water mixing ratio", "kg/kg", "ztt",   &
                          ANY(outputlist == "rc"), pipeline)
      END IF

      IF (ANY(outputlist == "rr")) THEN
         pipeline => NULL()
         ps_rr = FloatArray1d()
         ps_rr%onDemand => globalAvgProfile
         pipeline => ps_rr
         CALL PS%newField("rr", "Precipitation mixing ratio", "kg/kg", "ztt",   &
                          ANY(outputlist == "rr"), pipeline)
      END IF
         
      IF (level == 5) THEN
         IF (ANY(outputlist == "ri")) THEN
            pipeline => NULL()
            ps_ri = FloatArray1d()
            ps_ri%onDemand => globalAvgProfile
            pipeline => ps_ri
            CALL PS%newField("ri", "Unrimed ice mixing ratio", "kg/kg", "ztt",   &
                             ANY(outputlist == "ri"), pipeline)
         END IF

         IF (ANY(outputlist == "riri")) THEN
            pipeline => NULL()
            ps_riri = FloatArray1d()
            ps_riri%onDemand => globalAvgProfile
            pipeline => ps_riri
            CALL PS%newField("riri", "Rimed ice mixing ratio", "kg/kg", "ztt",   &
                             ANY(outputlist == "riri"), pipeline)
         END IF
      END IF

      IF (level >= 4) THEN
         IF (ANY(outputlist == "Naa")) THEN
            pipeline => NULL()
            ps_Naa = FloatArray1d()
            ps_Naa%onDemand => globalAvgProfile
            pipeline => ps_Naa
            CALL PS%newField("Naa", "Aerosol A, bulk number concentration", "m-3", "ztt",   &
                             ANY(outputlist == "Naa"), pipeline)
         END IF

         IF (ANY(outputlist == "Nab")) THEN
            pipeline => NULL()
            ps_Nab = FloatArray1d()
            ps_Nab%onDemand => globalAvgProfile
            pipeline => ps_Nab
            CALL PS%newField("Nab", "Aerosol B, bulk number concentration", "m-3", "ztt",   &
                             ANY(outputlist == "Nab"), pipeline)
         END IF

         IF (ANY(outputlist == "Nca")) THEN
            pipeline => NULL()
            ps_Nca = FloatArray1d()
            ps_Nca%onDemand => globalAvgProfile
            pipeline => ps_Nca
            CALL PS%newField("Nca", "Cloud A, bulk number concentration", "m-3", "ztt",   &
                             ANY(outputlist == "Nca"), pipeline)
         END IF

         IF (ANY(outputlist == "Ncb")) THEN
            pipeline => NULL()
            ps_Ncb = FloatArray1d()
            ps_Ncb%onDemand => globalAvgProfile
            pipeline => ps_Ncb
            CALL PS%newField("Ncb", "Cloud B, bulk number concentration", "m-3", "ztt",   &
                             ANY(outputlist == "Ncb"), pipeline)
         END IF

         IF (ANY(outputlist == "Np")) THEN
            pipeline => NULL()
            ps_Np = FloatArray1d()
            ps_Np%onDemand => globalAvgProfile
            pipeline => ps_Np
            CALL PS%newField("Np", "Precipitation bulk number concentration", "m-3", "ztt",   &
                             ANY(outputlist == "Np"), pipeline)
         END IF
      END IF

      IF (level == 5 .AND. ANY(outputlist == "Ni")) THEN
         pipeline => NULL()
         ps_Ni = FloatArray1d()
         ps_Ni%onDemand => globalAvgProfile
         pipeline => ps_Ni
         CALL PS%newField("Ni", "Ice bulk number concentration", "m-3", "ztt",   &
                          ANY(outputlist == "Ni"), pipeline)
      END IF

      IF (ANY(outputlist == "RH")) THEN
         pipeline => NULL()
         ps_RH = FloatArray1d()
         ps_RH%onDemand => globalAvgProfile
         pipeline => ps_RH
         CALL PS%newField("RH", "Relative humidity", "1", "ztt",   &
                          ANY(outputlist == "RH"), pipeline)
      END IF

      IF (ANY(outputlist == "rsl")) THEN
         pipeline => NULL()
         ps_rsl = FloatArray1d()
         ps_rsl%onDemand => globalAvgProfile
         pipeline => ps_rsl
         CALL PS%newField("rsl", "Saturation mixing ratio, liquid", "kg/kg", "ztt",   &
                          ANY(outputlist == "rsl"), pipeline)
      END IF
         
      IF (level == 5) THEN
         IF (ANY(outputlist == "RHI")) THEN
            pipeline => NULL()
            ps_RHI = FloatArray1d()
            ps_RHI%onDemand => globalAvgProfile
            pipeline => ps_RHI
            CALL PS%newField("RHI", "Relative humidity over ice", "1", "ztt",   &
                             ANY(outputlist == "RHI"), pipeline)
         END IF

         IF (ANY(outputlist == "rsi")) THEN
            pipeline => NULL()
            ps_rsi = FloatArray1d()
            ps_rsi%onDemand => globalAvgProfile
            pipeline => ps_rsi
            CALL PS%newField("rsi", "Saturation mixing ratio, ice", "kg/kg", "ztt",   &
                             ANY(outputlist == "rsi"), pipeline)
         END IF         
      END IF

      IF (ANY(outputlist == "rrate")) THEN
         pipeline => NULL()
         ps_rrate  = FloatArray1d()
         ps_rrate%onDemand => globalAvgProfile
         pipeline => ps_rrate
         CALL PS%newField("rrate", "Precipitation flux", "kg/m2s", "ztt",   &
                          ANY(outputlist == "rrate"), pipeline)
      END IF
         
      IF (level == 5 .AND. ANY(outputlist == "irate")) THEN
         pipeline => NULL()
         ps_irate = FloatArray1d()
         ps_irate%onDemand => globalAvgProfile
         pipeline => ps_irate
         CALL PS%newField("irate", "Ice precipitation flux", "kg/m2s", "ztt",   &
                          ANY(outputlist == "irate"), pipeline)
      END IF

      IF (ANY(outputlist == "rflx")) THEN
         pipeline => NULL()
         ps_rflx = FloatArray1d()
         ps_rflx%onDemand => globalAvgProfile
         pipeline => ps_rflx
         CALL PS%newField("rflx", "Net LW flux", "W/m2", "ztt",   &
                          ANY(outputlist == "rflx"), pipeline)
      END IF

      IF (ANY(outputlist == "sflx")) THEN
         pipeline => NULL()
         ps_sflx = FloatArray1d()
         ps_sflx%onDemand => globalAvgProfile
         pipeline => ps_sflx
         CALL PS%newField("sflx", "Net SW flux", "W/m2", "ztt",   &
                          ANY(outputlist == "sflx"), pipeline)
      END IF

      IF (ANY(outputlist == "lwup")) THEN
         pipeline => NULL()
         ps_lwup = FloatArray1d()
         ps_lwup%onDemand => globalAvgProfile
         pipeline => ps_lwup
         CALL PS%newField("lwup", "Upward LW flux", "W/m2", "ztt",   &
                          ANY(outputlist == "lwup"), pipeline)
      END IF

      IF (ANY(outputlist == "lwdn")) THEN
         pipeline => NULL()
         ps_lwdn = FloatArray1d()
         ps_lwdn%onDemand => globalAvgProfile
         pipeline => ps_lwdn
         CALL PS%newField("lwdn", "Downward LW flux", "W/m2", "ztt",   &
                          ANY(outputlist == "lwdn"), pipeline)
      END IF

      IF (ANY(outputlist == "swup")) THEN
         pipeline => NULL()
         ps_swup = FloatArray1d()
         ps_swup%onDemand => globalAvgProfile
         pipeline => ps_swup
         CALL PS%newField("swup", "Upward SW flux", "W/m2", "ztt",   &
                          ANY(outputlist == "swup"), pipeline)
      END IF

      IF (ANY(outputlist == "swdn")) THEN
         pipeline => NULL()
         ps_swdn = FloatArray1d()
         ps_swdn%onDemand => globalAvgProfile
         pipeline => ps_swdn
         CALL PS%newField("swdn", "Downward SW flux", "W/m2", "ztt",   &
                          ANY(outputlist == "swdn"), pipeline)
      END IF
         
      IF (level >= 4) THEN
         IF (ANY(outputlist == "Dwaa")) THEN
            pipeline => NULL()
            ps_Dwaa = FloatArray1d()
            ps_Dwaa%onDemand => globalAvgProfile
            pipeline => ps_Dwaa
            CALL PS%newField("Dwaa", "Aerosol A bulk wet diameter", "m", "ztt",   &
                             ANY(outputlist == "Dwaa"), pipeline)
         END IF

         IF (ANY(outputlist == "Dwab")) THEN
            pipeline => NULL()
            ps_Dwab = FloatArray1d()
            ps_Dwab%onDemand => globalAvgProfile
            pipeline => ps_Dwab
            CALL PS%newField("Dwab", "Aerosol B bulk wet diameter", "m", "ztt",   &
                             ANY(outputlist == "Dwab"), pipeline)
         END IF

         IF (ANY(outputlist == "Dwca")) THEN
            pipeline => NULL()
            ps_Dwca = FloatArray1d()
            ps_Dwca%onDemand => globalAvgProfile
            pipeline => ps_Dwca
            CALL PS%newField("Dwca", "Cloud A bulk wet diameter", "m", "ztt",   &
                             ANY(outputlist == "Dwca"), pipeline)
         END IF

         IF (ANY(outputlist == "Dwcb")) THEN
            pipeline => NULL()
            ps_Dwcb = FloatArray1d()
            ps_Dwcb%onDemand => globalAvgProfile
            pipeline => ps_Dwcb
            CALL PS%newField("Dwcb", "Cloud B bulk wet diameter", "m", "ztt",   &
                             ANY(outputlist == "Dwcb"), pipeline)
         END IF

         IF (ANY(outputlist == "Dwpa")) THEN
            pipeline => NULL()
            ps_Dwpa = FloatArray1d()
            ps_Dwpa%onDemand => globalAvgProfile
            pipeline => ps_Dwpa
            CALL PS%newField("Dwpa", "Precipitation bulk wet diameter", "m", "ztt",   &
                             ANY(outputlist == "Dwpa"), pipeline)
         END IF
      END IF

      IF (level == 5 .AND. ANY(outputlist == "Dwia")) THEN
         pipeline => NULL()
         ps_Dwia = FloatArray1d()
         ps_Dwia%onDemand => globalAvgProfile
         pipeline => ps_Dwia
         CALL PS%newField("Dwia", "Ice bulk wet diameter (spherical)", "m", "ztt",   &
                          ANY(outputlist == "Dwia"), pipeline)
      END IF

      IF (level >= 4) THEN
         IF (ANY(outputlist == "aSO4a")) THEN
            pipeline => NULL()
            ps_aSO4a = FloatArray1d()
            ps_aSO4a%onDemand => globalAvgProfile
            pipeline => ps_aSO4a
            CALL PS%newField("aSO4a", "Aerosol A bulk SO4 mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aSO4a"), pipeline)
         END IF

         IF (ANY(outputlist == "aSO4b")) THEN
            pipeline => NULL()
            ps_aSO4b = FloatArray1d()
            ps_aSO4b%onDemand => globalAvgProfile
            pipeline => ps_aSO4b
            CALL PS%newField("aSO4b", "Aerosol B bulk SO4 mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aSO4b"), pipeline)
         END IF

         IF (ANY(outputlist == "cSO4a")) THEN
            pipeline => NULL()
            ps_cSO4a = FloatArray1d()
            ps_cSO4a%onDemand => globalAvgProfile
            pipeline => ps_cSO4a
            CALL PS%newField("cSO4a", "Cloud A bulk SO4 mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cSO4a"), pipeline)
         END IF

         IF (ANY(outputlist == "cSO4b")) THEN
            pipeline => NULL()
            ps_cSO4b = FloatArray1d()
            ps_cSO4b%onDemand => globalAvgProfile
            pipeline => ps_cSO4b
            CALL PS%newField("cSO4b", "Cloud B bulk SO4 mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cSO4b"), pipeline)
         END IF

         IF (ANY(outputlist == "pSO4a")) THEN
            pipeline => NULL()
            ps_pSO4a = FloatArray1d()
            ps_pSO4a%onDemand => globalAvgProfile
            pipeline => ps_pSO4a
            CALL PS%newField("pSO4a", "Precipitation bulk SO4 mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "pSO4a"), pipeline)
         END IF
      END IF

      IF (level == 5 .AND. ANY(outputlist == "iSO4a")) THEN
         pipeline => NULL()
         ps_iSO4a = FloatArray1d()
         ps_iSO4a%onDemand => globalAvgProfile
         pipeline => ps_iSO4a
         CALL PS%newField("iSO4a", "Ice bulk SO4 mass", "kg/kg", "ztt",   &
                          ANY(outputlist == "iSO4a"), pipeline)
      END IF

      IF (level >= 4) THEN
         IF (ANY(outputlist == "aOCa")) THEN
            pipeline => NULL()
            ps_aOCa = FloatArray1d()
            ps_aOCa%onDemand => globalAvgProfile
            pipeline => ps_aOCa
            CALL PS%newField("aOCa", "Aerosol A bulk OC mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aOCa"), pipeline)
         END IF

         IF (ANY(outputlist == "aOCb")) THEN
            pipeline => NULL()
            ps_aOCb = FloatArray1d()
            ps_aOCb%onDemand => globalAvgProfile
            pipeline => ps_aOCb
            CALL PS%newField("aOCb", "Aerosol B bulk OC mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aOCb"), pipeline)
         END IF

         IF (ANY(outputlist == "cOCa")) THEN
            pipeline => NULL()
            ps_cOCa = FloatArray1d()
            ps_cOCa%onDemand => globalAvgProfile
            pipeline => ps_cOCa
            CALL PS%newField("cOCa", "Cloud A bulk OC mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cOCa"), pipeline)
         END IF

         IF (ANY(outputlist == "cOCb")) THEN
            pipeline => NULL()
            ps_cOCb = FloatArray1d()
            ps_cOCb%onDemand => globalAvgProfile
            pipeline => ps_cOCb
            CALL PS%newField("cOCb", "Cloud B bulk OC mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cOCb"), pipeline)
         END IF

         IF (ANY(outputlist == "pOCa")) THEN
            pipeline => NULL()
            ps_pOCa = FloatArray1d()
            ps_pOCa%onDemand => globalAvgProfile
            pipeline => ps_pOCa
            CALL PS%newField("pOCa", "Precipitation bulk OC mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "pOCa"), pipeline)
         END IF
      END IF
         
      IF (level == 5 .AND. ANY(outputlist == "iOCa")) THEN
         pipeline => NULL()
         ps_iOCa = FloatArray1d()
         ps_iOCa%onDemand => globalAvgProfile
         pipeline => ps_iOCa
         CALL PS%newField("iOCa", "Ice bulk OC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == "iOCa"), pipeline)
      END IF

      IF (level >= 4) THEN
         IF (ANY(outputlist == "aBCa")) THEN
            pipeline => NULL()
            ps_aBCa = FloatArray1d()
            ps_aBCa%onDemand => globalAvgProfile
            pipeline => ps_aBCa
            CALL PS%newField("aBCa", "Aerosol A bulk BC mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aBCa"), pipeline)
         END IF
         
         IF (ANY(outputlist == "aBCb")) THEN
            pipeline => NULL()
            ps_aBCb = FloatArray1d()
            ps_aBCb%onDemand => globalAvgProfile
            pipeline => ps_aBCb
            CALL PS%newField("aBCb", "Aerosol B bulk BC mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aBCb"), pipeline)
         END IF

         IF (ANY(outputlist == "cBCa")) THEN
            pipeline => NULL()
            ps_cBCa = FloatArray1d()
            ps_cBCa%onDemand => globalAvgProfile
            pipeline => ps_cBCa
            CALL PS%newField("cBCa", "Cloud A bulk BC mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cBCa"), pipeline)
         END IF

         IF (ANY(outputlist == "cBCb")) THEN
            pipeline => NULL()
            ps_cBCb = FloatArray1d()
            ps_cBCb%onDemand => globalAvgProfile
            pipeline => ps_cBCb
            CALL PS%newField("cBCb", "Cloud B bulk BC mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cBCb"), pipeline)
         END IF

         IF (ANY(outputlist == "pBCa")) THEN
            pipeline => NULL()
            ps_pBCa = FloatArray1d()
            ps_pBCa%onDemand => globalAvgProfile
            pipeline => ps_pBCa
            CALL PS%newField("pBCa", "Precipitation bulk BC mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "pBCa"), pipeline)
         END IF
      END IF

      IF (level == 5 .AND. ANY(outputlist == "iBCa")) THEN
         pipeline => NULL()
         ps_iBCa = FloatArray1d()
         ps_iBCa%onDemand => globalAvgProfile
         pipeline => ps_iBCa
         CALL PS%newField("iBCa", "Ice bulk BC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == "iBCa"), pipeline)
      END IF

      IF (level >= 4) THEN
         IF (ANY(outputlist == "aDUa")) THEN
            pipeline => NULL()
            ps_aDUa = FloatArray1d()
            ps_aDUa%onDemand => globalAvgProfile
            pipeline => ps_aDUa
            CALL PS%newField("aDUa", "Aerosol A bulk DU mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aDUa"), pipeline)
         END IF

         IF (ANY(outputlist == "aDUb")) THEN
            pipeline => NULL()
            ps_aDUb = FloatArray1d()
            ps_aDUb%onDemand => globalAvgProfile
            pipeline => ps_aDUb
            CALL PS%newField("aDUb", "Aerosol B bulk DU ", "kg/kg", "ztt",   &
                             ANY(outputlist == "aDUb"), pipeline)
         END IF

         IF (ANY(outputlist == "cDUa")) THEN
            ps_cDUa = FloatArray1d()
            ps_cDUa%onDemand => globalAvgProfile
            pipeline => ps_cDUa
            CALL PS%newField("cDUa", "Cloud A bulk DU mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cDUa"), pipeline)
         END IF

         IF (ANY(outputlist == "cDUb")) THEN
            pipeline => NULL()
            ps_cDUb = FloatArray1d()
            ps_cDUb%onDemand => globalAvgProfile
            pipeline => ps_cDUb
            CALL PS%newField("cDUb", "Cloud B bulk DU mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cDUb"), pipeline)
         END IF

         IF (ANY(outputlist == "pDUa")) THEN
            pipeline => NULL()
            ps_pDUa = FloatArray1d()
            ps_pDUa%onDemand => globalAvgProfile
            pipeline => ps_pDUa
            CALL PS%newField("pDUa", "Precipitation bulk DU mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "pDUa"), pipeline)
         END IF
      END IF

      IF (level == 5 .AND. ANY(outputlist == "iDUa")) THEN
         pipeline => NULL()
         ps_iDUa = FloatArray1d()
         ps_iDUa%onDemand => globalAvgProfile
         pipeline => ps_iDUa
         CALL PS%newField("iDUa", "Ice bulk DU mass", "kg/kg", "ztt",   &
                          ANY(outputlist == "iDUa"), pipeline)
      END IF

      IF (level >= 4) THEN
         IF (ANY(outputlist == "aSSa")) THEN
            pipeline => NULL()
            ps_aSSa = FloatArray1d()
            ps_aSSa%onDemand => globalAvgProfile
            pipeline => ps_aSSA
            CALL PS%newField("aSSa", "Aerosol A bulk SS mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aSSa"), pipeline)
         END IF

         IF (ANY(outputlist == "aSSb")) THEN
            pipeline => NULL()
            ps_aSSb = FloatArray1d()
            ps_aSSb%onDemand => globalAvgProfile
            pipeline => ps_aSSb
            CALL PS%newField("aSSb", "Aerosol B bulk SS mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aSSb"), pipeline)
         END IF

         IF (ANY(outputlist == "cSSa")) THEN
            pipeline => NULL()
            ps_cSSa = FloatArray1d()
            ps_cSSa%onDemand => globalAvgProfile
            pipeline => ps_cSSa
            CALL PS%newField("cSSa", "Cloud A bulk SS mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cSSa"), pipeline)
         END IF

         IF (ANY(outputlist == "cSSb")) THEN
            pipeline => NULL()
            ps_cSSb = FloatArray1d()
            ps_cSSb%onDemand => globalAvgProfile
            pipeline => ps_cSSb
            CALL PS%newField("cSSb", "Cloud B bulk SS mass", "K", "ztt",   &
                             ANY(outputlist == "cSSb"), pipeline)
         END IF

         IF (ANY(outputlist == "pSSa")) THEN
            pipeline => NULL()
            ps_pSSa = FloatArray1d()
            ps_pSSa%onDemand => globalAvgProfile
            pipeline => ps_pSSa
            CALL PS%newField("pSSa", "Precipitation bulk SS mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "pSSa"), pipeline)
         END IF
      END IF

      IF (level == 5 .AND. ANY(outputlist == "iSSa")) THEN
         pipeline => NULL()
         ps_iSSa = FloatArray1d()
         ps_iSSa%onDemand => globalAvgProfile
         pipeline => ps_iSSa
         CALL PS%newField("iSSa", "Ice bulk SS mass", "kg/kg", "ztt",   &
                          ANY(outputlist == "iSSa"), pipeline)
      END IF

      IF (level >= 4) THEN
         IF (ANY(outputlist == "aNOa")) THEN
            pipeline => NULL()
            ps_aNOa = FloatArray1d()
            ps_aNOa%onDemand => globalAvgProfile
            pipeline => ps_aNOa
            CALL PS%newField("aNOa", "Aerosol A bulk NO mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aNOa"), pipeline)
         END IF

         IF (ANY(outputlist == "aNOb")) THEN
            pipeline => NULL()
            ps_aNOb = FloatArray1d()
            ps_aNOb%onDemand => globalAvgProfile
            pipeline => ps_aNOb
            CALL PS%newField("aNOb", "Aerosol B bulk NO mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aNOb"), pipeline)
         END IF

         IF (ANY(outputlist == "cNOa")) THEN
            pipeline => NULL()
            ps_cNOa = FloatArray1d()
            ps_cNOa%onDemand => globalAvgProfile
            pipeline => ps_cNOa
            CALL PS%newField("cNOa", "Cloud A bulk NO mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cNOa"), pipeline)
         END IF
            
         IF (ANY(outputlist == "cNOb")) THEN
            pipeline => NULL()
            ps_cNOb = FloatArray1d()
            ps_cNOb%onDemand => globalAvgProfile
            pipeline => ps_cNOb
            CALL PS%newField("cNOb", "Cloud B bulk NO mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cNOb"), pipeline)
         END IF
            
         IF (ANY(outputlist == "pNOa")) THEN
            pipeline => NULL()
            ps_pNOa = FloatArray1d()
            ps_pNOa%onDemand => globalAvgProfile
            pipeline => ps_pNOa
            CALL PS%newField("pNOa", "Precipitation bulk NO mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "pNOa"), pipeline)
         END IF
      END IF

      IF (level == 5 .AND. ANY(outputlist == "iNOa")) THEN
         pipeline => NULL()
         ps_iNOa = FloatArray1d()
         ps_iNOa%onDemand => globalAvgProfile
         pipeline => ps_iNOa
         CALL PS%newField("iNOa", "Ice bulk NO mass", "kg/kg", "ztt",   &
                          ANY(outputlist == "iNOa"), pipeline)
      END IF

      IF (level >= 4) THEN
         IF (ANY(outputlist == "aNHa")) THEN
            pipeline => NULL()
            ps_aNHa = FloatArray1d()
            ps_aNHa%onDemand => globalAvgProfile
            pipeline => ps_aNHa
            CALL PS%newField("aNHa", "Aerosol A bulk NH mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aNHa"), pipeline)
         END IF
            
         IF (ANY(outputlist == "aNHb")) THEN
            pipeline => NULL()
            ps_aNHb = FloatArray1d()
            ps_aNHb%onDemand => globalAvgProfile
            pipeline => ps_aNHb
            CALL PS%newField("aNHb", "Aerosol B bulk NH mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "aNHb"), pipeline)
         END IF
            
         IF (ANY(outputlist == "cNHa")) THEN
            pipeline => NULL()
            ps_cNHa = FloatArray1d()
            ps_cNHa%onDemand => globalAvgProfile
            pipeline => ps_cNHa
            CALL PS%newField("cNHa", "Cloud A bulk NH mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cNHa"), pipeline)
         END IF
            
         IF (ANY(outputlist == "cNHb")) THEN
            pipeline => NULL()
            ps_cNHb = FloatArray1d()
            ps_cNHb%onDemand => globalAvgProfile
            pipeline => ps_cNHb
            CALL PS%newField("cNHb", "Cloud B bulk NH mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "cNHb"), pipeline)
         END IF
            
         IF (ANY(outputlist == "pNHa")) THEN
            pipeline => NULL()
            ps_pNHa = FloatArray1d()
            ps_pNHa%onDemand => globalAvgProfile
            pipeline => ps_pNHa
            CALL PS%newField("pNHa", "Precipitation bulk NH mass", "kg/kg", "ztt",   &
                             ANY(outputlist == "pNHa"), pipeline)
         END IF
      END IF

      IF (level == 5 .AND. ANY(outputlist == "iNHa")) THEN
         pipeline => NULL()
         ps_iNHa = FloatArray1d()
         ps_iNHa%onDemand => globalAvgProfile
         pipeline => ps_iNHa
         CALL PS%newField("iNHa", "Ice bulk NH mass", "kg/kg", "ztt",   &
                          ANY(outputlist == "iNHa"), pipeline)
      END IF

      pipeline => NULL()

      ! Check the user specified output variable list for bad entries
      CALL checkOutputs(outputlist,PS)
      
    END SUBROUTINE setPSVariables


END MODULE mo_ps_state
