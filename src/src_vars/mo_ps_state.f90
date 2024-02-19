MODULE mo_ps_state
  USE classFieldArray
  USE mo_structured_datatypes
  USE mo_ps_procedures
  IMPLICIT NONE

  ! Contains variable definitions for profile statistics outputs.
  ! STATISTICAL OUTPUT IS ONLY WRITTEN BY THE ROOT PROCESS, I.E. RANK=MPI_ROOT. THEREFORE MAKE SURE TO USE
  ! MPI REDUCTION OPERATIONS TO DEFINE THE FINAL RESULT IN THE ROOT PROCESS.
  
  ! Note some of these have to be stored in memory in order to calculate them in the physics routines,
  ! others are associated with and onDemand procedure to calculate the values upon output. The onDemand
  ! subroutines should be implemented in mo_ps_procedures and they must have identical IO interface (*name,*output).

  ! There are 3 ways to calculate the statistics:
  ! 1. For stats variables associated with an onDemand procedure but no storage, assume the stats variable to
  !    have IDENTICAL SHORT NAME WITH THE CORRESPONDING FULL FIELD from where they are calculated (assuming there is one).
  !    This allows you to use automatic fetching of data using mo_stats_finder in the subroutines found in mo_ps_procedures. 
  ! 2. Stats variables associated with an onDemand procedure, but a UNIQUE short name: An explicit paring of the stats variable
  !    and the full field source variable as well as the calculation of the statistic must be implemented in mo_ps_procedures
  !    and associated with the onDemand procedure.
  ! 3. Stats variables associated with storage space and no onDemand procedure: calculate the values during the timestep
  !    flagged for statistics output (psflg in mo_output.f90) within the appropriate physics subroutine.
  !
  ! Regarding 3), the likely case would be accumulation statistics, which require memory. These are yet to be implemented.
  ! Note that for all cases 1)-3), MPI reduction operations should be applied before writing output, since the
  ! statistics io will only be performed by rank=mpi_root (from mpi_interface.f90)

  ! Adding new variables:
  ! Declare your new TYPE(FloatArray1d) variable. Add initialization and registration of the variable to the TYPE(FieldArray) :: PS
  ! in the subroutine setPSVariables found below. Just check how other variables are defined and do the same. Only thing to consider
  ! are the 3 ways to calculate the statistics listed above.
  !   - If you want to use onDemand, first see if mo_ps_procedures already has a suitable subroutine to calculate your variable
  !     (e.g. simple domain mean is already implemented).
  !   - If so, and your variable fits point 1) for calculating the statistic, just associate the onDemand variable in the initialization
  !     and you're good to go. Nothing else needed.
  !   - If instead you see a suitable subroutine but your variable fits point 2), you must add a new paring for your variable and a
  !     source variable name, and the calculation of the statistic in mo_ps_procedures.
  !   - If you don't want to use onDemand, associate the variable with the storage instead. Now you need to calculate the value during
  !     a statistics timestep, or accumulate the value over consecutive timesteps using the TYPE(FloatArray1d) variable you declared.
  
  ! Regular 3d variables
  ! Simple mean profiles
  TYPE(FloatArray1d), TARGET :: ps_theta, ps_temp, ps_press, ps_pexnr,  &
                                ps_rp, ps_rc, ps_srp, ps_rpp, ps_npp,   &
                                ps_ri, ps_riri, ps_Naa, ps_Nab,         &
                                ps_Nca, ps_Ncb, ps_Np, ps_Ni,           &
                                ps_CDNC, ps_CNC, ps_Reff,               &                                
                                ps_rh, ps_rsl, ps_rhi, ps_rsi,          &
                                ps_dn,                                  &
                                                                
                                ps_rrate, ps_irate,              &

                                ps_wmin, ps_wmax,                &
                                ps_ubar, ps_vbar, ps_wbar,       &
                                ps_u_2, ps_v_2, ps_w_2,          &

                                ps_tw_res, ps_qw_res, ps_lw_res,     &
                                ps_uw_res, ps_vw_res, ps_ww_res,     &  
                                ps_tke_res, ps_tke_sgs,                        &  ! sgs subroutine missing

                                ps_rflx, ps_sflx, ps_lwup,       &
                                ps_lwdn, ps_swup, ps_swdn,       &

                                ps_Dwaa, ps_Dwab, ps_Dwca,       &
                                ps_Dwcb, ps_Dwpa, ps_Dwia,       &

                                ps_aSO4a, ps_aSO4b, ps_cSO4a, ps_cSO4b, ps_pSO4a, ps_iSO4a,  &
                                ps_aOCa, ps_aOCb, ps_cOCa, ps_cOCb, ps_pOCa, ps_iOCa,        &
                                ps_aBCa, ps_aBCb, ps_cBCa, ps_cBCb, ps_pBCa, ps_iBCa,        &
                                ps_aDUa, ps_aDUb, ps_cDUa, ps_cDUb, ps_pDUa, ps_iDUa,        &
                                ps_aSSa, ps_aSSb, ps_cSSa, ps_cSSb, ps_pSSa, ps_iSSa,        &
                                ps_aNOa, ps_aNOb, ps_cNOa, ps_cNOb, ps_pNOa, ps_iNOa,        &
                                ps_aNHa, ps_aNHb, ps_cNHa, ps_cNHb, ps_pNHa, ps_iNHa

  ! Conditionally sampled profiles
  TYPE(FloatArray1d), TARGET :: psic_rc, psic_Naa, psic_Nab, psic_Nca, psic_Ncb, psic_CDNC, psic_CNC,   &
                                psic_Reff, psic_rh   
  
  ! -- mean across precipitating grid boxes
  TYPE(FloatArray1d), TARGET :: pspr_rrate, pspr_irate, pspr_srp, pspr_rpp, pspr_npp, pspr_Np 
       
  ! binned variables
  TYPE(FloatArray2d), TARGET :: ps_Dwaba, ps_Dwabb, ps_Dwcba, ps_Dwcbb, ps_Dwpba, ps_Dwiba
  TYPE(FloatArray2d), TARGET :: ps_Naba, ps_Nabb, ps_Ncba, ps_Ncbb, ps_Npba, ps_Niba
  TYPE(FloatArray2d), TARGET :: ps_sipdrfr, ps_sipiibr, ps_siprmspl
  
  ! Piggybacking slave microphysics diagnostics
  TYPE(FloatArray1d), TARGET :: ps_pb_rrate, ps_pb_rc,         &
                                pspr_pb_rrate, psic_pb_rc

  ! Process rate diagnostics (in-cloud)
  TYPE(FloatArray1d), TARGET :: psic_s_m_autoc,           &
                                psic_s_m_autoc50,         &
                                psic_s_m_autoc80,         &
                                psic_s_m_accr,            &
                                psic_s_m_accr50,          &
                                psic_s_m_accr80,          &
                                psic_s_m_ACcoll_dry,      &
                                psic_s_m_APcoll_dry,      &
                                psic_s_m_AIcoll_dry,      &
                                psic_s_n_activ,           &
                                psic_s_n_icehom,          &
                                psic_s_n_icedep,          &
                                psic_s_n_iceimm,          &
                                psic_s_m_conda,           &
                                psic_s_m_condc,           &
                                psic_s_m_condp,           &
                                psic_s_m_condi,           &
                                  
                                psic_b_m_autoc,           &
                                psic_b_n_autoc,           &
                                psic_b_m_accr       
    
  CONTAINS

    SUBROUTINE setPSVariables(PS,outputlist,level,lpback)
      TYPE(FieldArray), INTENT(inout) :: PS
      CHARACTER(len=*), INTENT(in)   :: outputlist(:)
      INTEGER, INTENT(in)             :: level
      LOGICAL, INTENT(in)             :: lpback
      CLASS(*), POINTER :: pipeline => NULL()
      

      !---------------------------------------
      ! Simple domain mean profiles
      !---------------------------------------
      
      pipeline => NULL()
      ps_theta = FloatArray1d("theta")
      ps_theta%onDemand => globalMeanProfile
      pipeline => ps_theta
      CALL PS%newField(ps_theta%shortName, "Potential temperature", "K", "ztt",   &
                       ANY(outputlist == ps_theta%shortName), pipeline)

      pipeline => NULL()
      ps_temp = FloatArray1d("temp")
      ps_temp%onDemand => globalMeanProfile
      pipeline => ps_temp
      CALL PS%newField(ps_temp%shortName, "Abs temperature", "K", "ztt",   &
                       ANY(outputlist == ps_temp%shortName), pipeline)

      pipeline => NULL()
      ps_press = FloatArray1d("press")
      ps_press%onDemand => globalMeanProfile
      pipeline => ps_press
      CALL PS%newField(ps_press%shortName, "Pressure", "Pa", "ztt",   &
                       ANY(outputlist == ps_press%shortName), pipeline)

      pipeline => NULL()
      ps_pexnr = FloatArray1d("pexnr")
      ps_pexnr%onDemand => globalMeanProfile
      pipeline => ps_pexnr
      CALL PS%newField(ps_pexnr%shortName, "Exner function", "", "ztt",   &
                       ANY(outputlist == ps_pexnr%shortName), pipeline)      

      pipeline => NULL()
      ps_rp = FloatArray1d("rp")
      ps_rp%onDemand => globalMeanProfile
      pipeline => ps_rp
      CALL PS%newField(ps_rp%shortName, "Water vapor mixing ratio", "kg/kg", "ztt",   &
                       ANY(outputlist == ps_rp%shortName), pipeline)

      pipeline => NULL()
      ps_rc = FloatArray1d("rc")
      ps_rc%onDemand => globalMeanProfile
      pipeline => ps_rc
      CALL PS%newField(ps_rc%shortName, "Cloud water mixing ratio", "kg/kg", "ztt",   &
                       ANY(outputlist == ps_rc%shortName), pipeline)

      IF (level < 4) THEN
         pipeline => NULL()
         ps_rpp = FloatArray1d("rpp")
         ps_rpp%onDemand => globalMeanProfile
         pipeline => ps_rpp
         CALL PS%newField(ps_rpp%shortName, "Precipitation mixing ratio", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_rpp%shortName), pipeline)

         pipeline => NULL()
         ps_npp = FloatArray1d("npp")
         ps_npp%onDemand => globalMeanProfile
         pipeline => ps_npp
         CALL PS%newField(ps_npp%shortName, "Precipitation number conc", "#/kg", "ztt",   &
                          ANY(outputlist == ps_npp%shortName), pipeline)
         
      END IF
      
      IF (level >= 4) THEN
         pipeline => NULL()
         ps_srp = FloatArray1d("srp")
         ps_srp%onDemand => globalMeanProfile
         pipeline => ps_srp
         CALL PS%newField(ps_srp%shortName, "Precipitation mixing ratio", "kg/kg", "ztt",   &
              ANY(outputlist == ps_srp%shortName), pipeline)
      END IF
         
      IF (level == 5) THEN
         pipeline => NULL()
         ps_ri = FloatArray1d("ri")
         ps_ri%onDemand => globalMeanProfile
         pipeline => ps_ri
         CALL PS%newField(ps_ri%shortName, "Unrimed ice mixing ratio", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_ri%shortName), pipeline)

         pipeline => NULL()
         ps_riri = FloatArray1d("riri")
         ps_riri%onDemand => globalMeanProfile
         pipeline => ps_riri
         CALL PS%newField(ps_riri%shortName, "Rimed ice mixing ratio", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_riri%shortName), pipeline)
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         ps_Naa = FloatArray1d("Naa")
         ps_Naa%onDemand => globalMeanProfile
         pipeline => ps_Naa
         CALL PS%newField(ps_Naa%shortName, "Aerosol A, bulk number concentration", "#/kg", "ztt",   &
                          ANY(outputlist == ps_Naa%shortName), pipeline)

         pipeline => NULL()
         ps_Nab = FloatArray1d("Nab")
         ps_Nab%onDemand => globalMeanProfile
         pipeline => ps_Nab
         CALL PS%newField(ps_Nab%shortName, "Aerosol B, bulk number concentration", "#/kg", "ztt",   &
                          ANY(outputlist == ps_Nab%shortName), pipeline)

         pipeline => NULL()
         ps_Nca = FloatArray1d("Nca")
         ps_Nca%onDemand => globalMeanProfile
         pipeline => ps_Nca
         CALL PS%newField(ps_Nca%shortName, "Cloud A, bulk number concentration", "#/kg", "ztt",   &
                          ANY(outputlist == ps_Nca%shortName), pipeline)

         pipeline => NULL()
         ps_Ncb = FloatArray1d("Ncb")
         ps_Ncb%onDemand => globalMeanProfile
         pipeline => ps_Ncb
         CALL PS%newField(ps_Ncb%shortName, "Cloud B, bulk number concentration", "#/kg", "ztt",   &
                          ANY(outputlist == ps_Ncb%shortName), pipeline)

         pipeline => NULL()
         ps_Np = FloatArray1d("Np")
         ps_Np%onDemand => globalMeanProfile
         pipeline => ps_Np
         CALL PS%newField(ps_Np%shortName, "Precipitation bulk number concentration", "#/kg", "ztt",   &
                          ANY(outputlist == ps_Np%shortName), pipeline)

         pipeline => NULL()
         ps_Naba = FloatArray2d()
         ps_Naba%onDemand => globalMeanProfileBinned
         pipeline => ps_Naba
         CALL PS%newField("Naba", "Aerosol A binned number concentration", "#/kg", "zttaea",   &
                          ANY(outputlist == "Naba"), pipeline)

         pipeline => NULL()
         ps_Nabb = FloatArray2d()
         ps_Nabb%onDemand => globalMeanProfileBinned
         pipeline => ps_Nabb
         CALL PS%newfield("Nabb", "Aerosol B binned number concentration", "#/kg", "zttaeb",   &
                          ANY(outputlist == "Nabb"), pipeline)

         pipeline => NULL()
         ps_Ncba = FloatArray2d()
         ps_Ncba%onDemand => globalMeanProfileBinned
         pipeline => ps_Ncba
         CALL PS%newfield("Ncba", "Cloud A binned number concentration", "#/kg", "zttcla",   &
                          ANY(outputlist == "Ncba"), pipeline)         

         pipeline => NULL()
         ps_Ncbb = FloatArray2d()
         ps_Ncbb%onDemand => globalMeanProfileBinned
         pipeline => ps_Ncbb
         CALL PS%newfield("Ncbb", "Cloud B binned number concentration", "#/kg", "zttclb",   &
                          ANY(outputlist == "Ncbb"), pipeline)

         pipeline => NULL()
         ps_Npba = FloatArray2d()
         ps_Npba%onDemand => globalMeanProfileBinned
         pipeline => ps_Npba
         CALL PS%newfield("Npba", "Precipitation binned number concentration", "#/kg", "zttprc",   &
                          ANY(outputlist == "Npba"), pipeline)

         pipeline => NULL()
         ps_CDNC = FloatArray1d("CDNC")
         ps_CDNC%onDemand => globalMeanProfile
         pipeline => ps_CDNC
         CALL PS%newfield(ps_CDNC%shortName, "Cloud droplet number concentration", "m-3", "ztt",   &
                          ANY(outputlist == ps_CDNC%shortName), pipeline)

         pipeline => NULL()
         ps_CNC = FloatArray1d("CNC")
         ps_CNC%onDemand => globalMeanProfile
         pipeline => ps_CNC
         CALL PS%newfield(ps_CNC%shortName, "Hydrometeor number concentration", "m-3", "ztt",   &
                          ANY(outputlist == ps_CNC%shortName), pipeline)         

         pipeline => NULL()
         ps_Reff = FloatArray1d("Reff")
         ps_Reff%onDemand => globalMeanProfile
         pipeline => ps_Reff
         CALL PS%newfield(ps_Reff%shortName, "Droplet effective radius", "m", "ztt",   &
                          ANY(outputlist == ps_Reff%shortName), pipeline)         
         
      END IF

      IF (level == 5) THEN
         pipeline => NULL()
         ps_Ni = FloatArray1d("Ni")
         ps_Ni%onDemand => globalMeanProfile
         pipeline => ps_Ni
         CALL PS%newField(ps_Ni%shortName, "Ice bulk number concentration", "#/kg", "ztt",   &
                          ANY(outputlist == ps_Ni%shortName), pipeline)

         pipeline => NULL()
         ps_Niba = FloatArray2d()
         ps_Niba%onDemand => globalMeanProfileBinned
         pipeline => ps_Niba
         CALL PS%newField("Niba", "Ice binned number concentration", "#/kg", "zttice",   &
                          ANY(outputlist == "Niba"), pipeline)

         pipeline => NULL()
         ps_sipdrfr = FloatArray2d()
         ps_sipdrfr%onDemand => globalMeanProfileBinned
         pipeline => ps_sipdrfr
         CALL PS%newField("sipdrfr","Drop fracturing SIP concentration", "#/kg", "zttice",  &
              ANY(outputlist == "sipdrfr"), pipeline)

         pipeline => NULL()
         ps_sipiibr = FloatArray2d()
         ps_sipiibr%onDemand => globalMeanProfileBinned
         pipeline => ps_sipiibr
         CALL PS%newField("sipiibr","Ice-ice collisional breakup SIP concentration", "#/kg", "zttice",  &
              ANY(outputlist == "sipiibr"), pipeline)

         pipeline => NULL()
         ps_siprmspl = FloatArray2d()
         ps_siprmspl%onDemand => globalMeanProfileBinned
         pipeline => ps_siprmspl
         CALL PS%newField("siprmspl","Rime splintering SIP concentration", "#/kg", "zttice",  &
              ANY(outputlist == "siprmspl"), pipeline)
        
         
      END IF

      pipeline => NULL()
      ps_rh = FloatArray1d("rh")
      ps_rh%onDemand => globalMeanProfile
      pipeline => ps_rh
      CALL PS%newField(ps_rh%shortName, "Relative humidity", "1", "ztt",   &
                       ANY(outputlist == ps_rh%shortName), pipeline)
      
      pipeline => NULL()
      ps_rsl = FloatArray1d("rsl")
      ps_rsl%onDemand => globalMeanProfile
      pipeline => ps_rsl
      CALL PS%newField(ps_rsl%shortName, "Saturation mixing ratio, liquid", "kg/kg", "ztt",   &
                       ANY(outputlist == ps_rsl%shortName), pipeline)

      IF (level == 5) THEN
         pipeline => NULL()
         ps_rhi = FloatArray1d("rhi")
         ps_rhi%onDemand => globalMeanProfile
         pipeline => ps_rhi
         CALL PS%newField(ps_rhi%shortName, "Relative humidity over ice", "1", "ztt",   &
                          ANY(outputlist == ps_rhi%shortName), pipeline)

         pipeline => NULL()
         ps_rsi = FloatArray1d("rsi")
         ps_rsi%onDemand => globalMeanProfile
         pipeline => ps_rsi
         CALL PS%newField(ps_rsi%shortName, "Saturation mixing ratio, ice", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_rsi%shortName), pipeline)
      END IF

      pipeline => NULL()
      ps_dn = FloatArray1d("dn")
      ps_dn%onDemand => globalMeanProfile
      pipeline => ps_dn
      CALL PS%newField(ps_dn%shortName, "Air density", "kg/m3", "ztt",   &
                       ANY(outputlist == ps_dn%shortName), pipeline)      
      

      pipeline => NULL()
      ps_rrate  = FloatArray1d("rrate")
      ps_rrate%onDemand => globalMeanProfile
      pipeline => ps_rrate
      CALL PS%newField(ps_rrate%shortName, "Precipitation flux", "W m-2", "ztt",   &
           ANY(outputlist == ps_rrate%shortName), pipeline)

      pipeline => NULL()
      pspr_rrate  = FloatArray1d("pr_rrate",srcName="rrate")
      pspr_rrate%onDemand => globalMeanProfile 
      pipeline => pspr_rrate
      CALL PS%newField(pspr_rrate%shortName, "Conditionally avg precipitation flux", "W m-2", "ztt",   &
                       ANY(outputlist == pspr_rrate%shortName), pipeline)     
         
      IF (level == 5) THEN
         pipeline => NULL()
         ps_irate = FloatArray1d("irate")
         ps_irate%onDemand => globalMeanProfile
         pipeline => ps_irate
         CALL PS%newField(ps_irate%shortName, "Ice precipitation flux", "W m-2", "ztt",   &
              ANY(outputlist == ps_irate%shortName), pipeline)

         pipeline => NULL()
         pspr_irate = FloatArray1d("pr_irate",srcName="irate")
         pspr_irate%onDemand => globalMeanProfile 
         pipeline => pspr_irate
         CALL PS%newField(pspr_irate%shortName, "Conditionally avg ice precipitation flux", "W m-2", "ztt",   &
                          ANY(outputlist == pspr_irate%shortName), pipeline)
         
      END IF

      pipeline => NULL()
      ps_rflx = FloatArray1d("rflx")
      ps_rflx%onDemand => globalMeanProfile
      pipeline => ps_rflx
      CALL PS%newField(ps_rflx%shortName, "Net rad flux", "W/m2", "ztt",   &
                       ANY(outputlist == ps_rflx%shortName), pipeline)

      pipeline => NULL()
      ps_sflx = FloatArray1d("sflx")
      ps_sflx%onDemand => globalMeanProfile
      pipeline => ps_sflx
      CALL PS%newField(ps_sflx%shortName, "Net SW flux", "W/m2", "ztt",   &
                       ANY(outputlist == ps_sflx%shortName), pipeline)

      pipeline => NULL()
      ps_lwup = FloatArray1d("lwup",srcName="fuir")
      ps_lwup%onDemand => globalMeanProfile
      pipeline => ps_lwup
      CALL PS%newField(ps_lwup%shortName, "Upward LW flux", "W/m2", "ztt",   &
                       ANY(outputlist == ps_lwup%shortName), pipeline)

      pipeline => NULL()
      ps_lwdn = FloatArray1d("lwdn",srcName="fdir")
      ps_lwdn%onDemand => globalMeanProfile
      pipeline => ps_lwdn
      CALL PS%newField(ps_lwdn%shortName, "Downward LW flux", "W/m2", "ztt",   &
                       ANY(outputlist == ps_lwdn%shortName), pipeline)

      pipeline => NULL()
      ps_swup = FloatArray1d("swup",srcName="fus")
      ps_swup%onDemand => globalMeanProfile
      pipeline => ps_swup
      CALL PS%newField(ps_swup%shortName, "Upward SW flux", "W/m2", "ztt",   &
                       ANY(outputlist == ps_swup%shortName), pipeline)

      pipeline => NULL()
      ps_swdn = FloatArray1d("swdn",srcName="fds")
      ps_swdn%onDemand => globalMeanProfile
      pipeline => ps_swdn
      CALL PS%newField(ps_swdn%shortName, "Downward SW flux", "W/m2", "ztt",   &
                       ANY(outputlist == ps_swdn%shortName), pipeline)

         
      IF (level >= 4) THEN
         pipeline => NULL()
         ps_Dwaa = FloatArray1d("Dwaa")
         ps_Dwaa%onDemand => globalMeanProfile
         pipeline => ps_Dwaa
         CALL PS%newField(ps_Dwaa%shortName, "Aerosol A bulk wet diameter", "m", "ztt",   &
                          ANY(outputlist == ps_Dwaa%shortName), pipeline)

         pipeline => NULL()
         ps_Dwab = FloatArray1d("Dwab")
         ps_Dwab%onDemand => globalMeanProfile
         pipeline => ps_Dwab
         CALL PS%newField(ps_Dwab%shortName, "Aerosol B bulk wet diameter", "m", "ztt",   &
                          ANY(outputlist == ps_Dwab%shortName), pipeline)

         pipeline => NULL()
         ps_Dwca = FloatArray1d("Dwca")
         ps_Dwca%onDemand => globalMeanProfile
         pipeline => ps_Dwca
         CALL PS%newField(ps_Dwca%shortName, "Cloud A bulk wet diameter", "m", "ztt",   &
                          ANY(outputlist == ps_Dwca%shortName), pipeline)

         pipeline => NULL()
         ps_Dwcb = FloatArray1d("Dwcb")
         ps_Dwcb%onDemand => globalMeanProfile
         pipeline => ps_Dwcb
         CALL PS%newField(ps_Dwcb%shortName, "Cloud B bulk wet diameter", "m", "ztt",   &
                          ANY(outputlist == ps_Dwcb%shortName), pipeline)

         pipeline => NULL()
         ps_Dwpa = FloatArray1d("Dwpa")
         ps_Dwpa%onDemand => globalMeanProfile
         pipeline => ps_Dwpa
         CALL PS%newField(ps_Dwpa%shortName, "Precipitation bulk wet diameter", "m", "ztt",   &
                          ANY(outputlist == ps_Dwpa%shortName), pipeline)

         pipeline => NULL()
         ps_Dwaba = FloatArray2d()
         ps_Dwaba%onDemand => globalMeanProfileBinned
         pipeline => ps_Dwaba
         CALL PS%newField("Dwaba", "Aerosol A binned diameter", "m", "zttaea",    &
                          ANY(outputlist == "Dwaba"), pipeline)

         pipeline => NULL()
         ps_Dwabb = FloatArray2d()
         ps_Dwabb%onDemand => globalMeanProfileBinned
         pipeline => ps_Dwabb
         CALL PS%newField("Dwabb", "Aerosol B binned diameter", "m", "zttaeb",    &
                          ANY(outputlist == "Dwabb"), pipeline)

         pipeline => NULL()
         ps_Dwcba = FloatArray2d()
         ps_Dwcba%onDemand => globalMeanProfileBinned
         pipeline => ps_Dwcba
         CALL PS%newField("Dwcba", "Cloud A binned diameter", "m", "zttcla",    &
                          ANY(outputlist == "Dwcba"), pipeline)

         pipeline => NULL()
         ps_Dwcbb = FloatArray2d()
         ps_Dwcbb%onDemand => globalMeanProfileBinned
         pipeline => ps_Dwcbb
         CALL PS%newField("Dwcbb", "Cloud B binned diameter", "m", "zttclb",    &
                          ANY(outputlist == "Dwcbb"), pipeline)

         pipeline => NULL()
         ps_Dwpba = FloatArray2d()
         ps_Dwpba%onDemand => globalMeanProfileBinned
         pipeline => ps_Dwpba
         CALL PS%newField("Dwpba", "Precipitation binned diameter", "m", "zttprc",    &
                          ANY(outputlist == "Dwpba"), pipeline)
         
      END IF

      IF (level == 5) THEN
         pipeline => NULL()
         ps_Dwia = FloatArray1d("Dwia")
         ps_Dwia%onDemand => globalMeanProfile
         pipeline => ps_Dwia
         CALL PS%newField(ps_Dwia%shortName, "Ice bulk wet diameter (spherical)", "m", "ztt",   &
                          ANY(outputlist == ps_Dwia%shortName), pipeline)

         pipeline => NULL()
         ps_Dwiba = FloatArray2d()
         ps_Dwiba%onDemand => globalMeanProfileBinned
         pipeline => ps_Dwiba
         CALL PS%newField("Dwiba", "Ice binned diameter", "m", "zttice",    &
                          ANY(outputlist == "Dwiba"), pipeline)         
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         ps_aSO4a = FloatArray1d("aSO4a")
         ps_aSO4a%onDemand => globalMeanProfile
         pipeline => ps_aSO4a
         CALL PS%newField(ps_aSO4a%shortName, "Aerosol A bulk SO4 mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aSO4a%shortName), pipeline)

         pipeline => NULL()
         ps_aSO4b = FloatArray1d("aSO4b")
         ps_aSO4b%onDemand => globalMeanProfile
         pipeline => ps_aSO4b
         CALL PS%newField(ps_aSO4b%shortName, "Aerosol B bulk SO4 mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aSO4b%shortName), pipeline)

         pipeline => NULL()
         ps_cSO4a = FloatArray1d("cSO4a")
         ps_cSO4a%onDemand => globalMeanProfile
         pipeline => ps_cSO4a
         CALL PS%newField(ps_cSO4a%shortName, "Cloud A bulk SO4 mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cSO4a%shortName), pipeline)

         pipeline => NULL()
         ps_cSO4b = FloatArray1d("cSO4b")
         ps_cSO4b%onDemand => globalMeanProfile
         pipeline => ps_cSO4b
         CALL PS%newField(ps_cSO4b%shortName, "Cloud B bulk SO4 mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cSO4b%shortName), pipeline)

         pipeline => NULL()
         ps_pSO4a = FloatArray1d("pSO4a")
         ps_pSO4a%onDemand => globalMeanProfile
         pipeline => ps_pSO4a
         CALL PS%newField(ps_pSO4a%shortName, "Precipitation bulk SO4 mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_pSO4a%shortName), pipeline)
      END IF

      IF (level == 5) THEN
         pipeline => NULL()
         ps_iSO4a = FloatArray1d("iSO4a")
         ps_iSO4a%onDemand => globalMeanProfile
         pipeline => ps_iSO4a
         CALL PS%newField(ps_iSO4a%shortName, "Ice bulk SO4 mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_iSO4a%shortName), pipeline)
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         ps_aOCa = FloatArray1d("aOCa")
         ps_aOCa%onDemand => globalMeanProfile
         pipeline => ps_aOCa
         CALL PS%newField(ps_aOCa%shortName, "Aerosol A bulk OC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aOCa%shortName), pipeline)

         pipeline => NULL()
         ps_aOCb = FloatArray1d("aOCb")
         ps_aOCb%onDemand => globalMeanProfile
         pipeline => ps_aOCb
         CALL PS%newField(ps_aOCb%shortName, "Aerosol B bulk OC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aOCb%shortName), pipeline)

         pipeline => NULL()
         ps_cOCa = FloatArray1d("cOCa")
         ps_cOCa%onDemand => globalMeanProfile
         pipeline => ps_cOCa
         CALL PS%newField(ps_cOCa%shortName, "Cloud A bulk OC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cOCa%shortName), pipeline)

         pipeline => NULL()
         ps_cOCb = FloatArray1d("cOCb")
         ps_cOCb%onDemand => globalMeanProfile
         pipeline => ps_cOCb
         CALL PS%newField(ps_cOCb%shortName, "Cloud B bulk OC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cOCb%shortName), pipeline)

         pipeline => NULL()
         ps_pOCa = FloatArray1d("pOCa")
         ps_pOCa%onDemand => globalMeanProfile
         pipeline => ps_pOCa
         CALL PS%newField(ps_pOCa%shortName, "Precipitation bulk OC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_pOCa%shortName), pipeline)
      END IF
         
      IF (level == 5) THEN
         pipeline => NULL()
         ps_iOCa = FloatArray1d("iOCa")
         ps_iOCa%onDemand => globalMeanProfile
         pipeline => ps_iOCa
         CALL PS%newField(ps_iOCa%shortName, "Ice bulk OC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_iOCa%shortName), pipeline)
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         ps_aBCa = FloatArray1d("aBCa")
         ps_aBCa%onDemand => globalMeanProfile
         pipeline => ps_aBCa
         CALL PS%newField(ps_aBCa%shortName, "Aerosol A bulk BC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aBCa%shortName), pipeline)

         pipeline => NULL()
         ps_aBCb = FloatArray1d("aBCb")
         ps_aBCb%onDemand => globalMeanProfile
         pipeline => ps_aBCb
         CALL PS%newField(ps_aBCb%shortName, "Aerosol B bulk BC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aBCb%shortName), pipeline)

         pipeline => NULL()
         ps_cBCa = FloatArray1d("cBCa")
         ps_cBCa%onDemand => globalMeanProfile
         pipeline => ps_cBCa
         CALL PS%newField(ps_cBCa%shortName, "Cloud A bulk BC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cBCa%shortName), pipeline)

         pipeline => NULL()
         ps_cBCb = FloatArray1d("cBCb")
         ps_cBCb%onDemand => globalMeanProfile
         pipeline => ps_cBCb
         CALL PS%newField(ps_cBCb%shortName, "Cloud B bulk BC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cBCb%shortName), pipeline)

         pipeline => NULL()
         ps_pBCa = FloatArray1d("pBCa")
         ps_pBCa%onDemand => globalMeanProfile
         pipeline => ps_pBCa
         CALL PS%newField(ps_pBCa%shortName, "Precipitation bulk BC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_pBCa%shortName), pipeline)
      END IF

      IF (level == 5) THEN
         pipeline => NULL()
         ps_iBCa = FloatArray1d("iBCa")
         ps_iBCa%onDemand => globalMeanProfile
         pipeline => ps_iBCa
         CALL PS%newField(ps_iBCa%shortName, "Ice bulk BC mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_iBCa%shortName), pipeline)
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         ps_aDUa = FloatArray1d("aDUa")
         ps_aDUa%onDemand => globalMeanProfile
         pipeline => ps_aDUa
         CALL PS%newField(ps_aDUa%shortName, "Aerosol A bulk DU mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aDUa%shortName), pipeline)

         pipeline => NULL()
         ps_aDUb = FloatArray1d("aDUb")
         ps_aDUb%onDemand => globalMeanProfile
         pipeline => ps_aDUb
         CALL PS%newField(ps_aDUb%shortName, "Aerosol B bulk DU ", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aDUb%shortName), pipeline)

         pipeline => NULL()
         ps_cDUa = FloatArray1d("cDUa")
         ps_cDUa%onDemand => globalMeanProfile
         pipeline => ps_cDUa
         CALL PS%newField(ps_cDUa%shortName, "Cloud A bulk DU mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cDUa%shortName), pipeline)

         pipeline => NULL()
         ps_cDUb = FloatArray1d("cDUb")
         ps_cDUb%onDemand => globalMeanProfile
         pipeline => ps_cDUb
         CALL PS%newField(ps_cDUb%shortName, "Cloud B bulk DU mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cDUb%shortName), pipeline)

         pipeline => NULL()
         ps_pDUa = FloatArray1d("pDUa")
         ps_pDUa%onDemand => globalMeanProfile
         pipeline => ps_pDUa
         CALL PS%newField(ps_pDUa%shortName, "Precipitation bulk DU mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_pDUa%shortName), pipeline)
      END IF

      IF (level == 5) THEN
         pipeline => NULL()
         ps_iDUa = FloatArray1d("iDUa")
         ps_iDUa%onDemand => globalMeanProfile
         pipeline => ps_iDUa
         CALL PS%newField(ps_iDUa%shortName, "Ice bulk DU mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_iDUa%shortName), pipeline)
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         ps_aSSa = FloatArray1d("aSSa")
         ps_aSSa%onDemand => globalMeanProfile
         pipeline => ps_aSSA
         CALL PS%newField(ps_aSSa%shortName, "Aerosol A bulk SS mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aSSa%shortName), pipeline)

         pipeline => NULL()
         ps_aSSb = FloatArray1d("aSSb")
         ps_aSSb%onDemand => globalMeanProfile
         pipeline => ps_aSSb
         CALL PS%newField(ps_aSSb%shortName, "Aerosol B bulk SS mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aSSb%shortName), pipeline)

         pipeline => NULL()
         ps_cSSa = FloatArray1d("cSSa")
         ps_cSSa%onDemand => globalMeanProfile
         pipeline => ps_cSSa
         CALL PS%newField(ps_cSSa%shortName, "Cloud A bulk SS mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cSSa%shortName), pipeline)

         pipeline => NULL()
         ps_cSSb = FloatArray1d("cSSb")
         ps_cSSb%onDemand => globalMeanProfile
         pipeline => ps_cSSb
         CALL PS%newField(ps_cSSb%shortName, "Cloud B bulk SS mass", "K", "ztt",   &
                          ANY(outputlist == ps_cSSb%shortName), pipeline)

         pipeline => NULL()
         ps_pSSa = FloatArray1d("pSSa")
         ps_pSSa%onDemand => globalMeanProfile
         pipeline => ps_pSSa
         CALL PS%newField(ps_pSSa%shortName, "Precipitation bulk SS mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_pSSa%shortName), pipeline)

      END IF

      IF (level == 5) THEN
         pipeline => NULL()
         ps_iSSa = FloatArray1d("iSSa")
         ps_iSSa%onDemand => globalMeanProfile
         pipeline => ps_iSSa
         CALL PS%newField(ps_iSSa%shortName, "Ice bulk SS mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_iSSa%shortName), pipeline)
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         ps_aNOa = FloatArray1d("aNOa")
         ps_aNOa%onDemand => globalMeanProfile
         pipeline => ps_aNOa
         CALL PS%newField(ps_aNOa%shortName, "Aerosol A bulk NO mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aNOa%shortName), pipeline)

         pipeline => NULL()
         ps_aNOb = FloatArray1d("aNOb")
         ps_aNOb%onDemand => globalMeanProfile
         pipeline => ps_aNOb
         CALL PS%newField(ps_aNOb%shortName, "Aerosol B bulk NO mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aNOb%shortName), pipeline)

         pipeline => NULL()
         ps_cNOa = FloatArray1d("cNOa")
         ps_cNOa%onDemand => globalMeanProfile
         pipeline => ps_cNOa
         CALL PS%newField(ps_cNOa%shortName, "Cloud A bulk NO mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cNOa%shortName), pipeline)

         pipeline => NULL()
         ps_cNOb = FloatArray1d("cNOb")
         ps_cNOb%onDemand => globalMeanProfile
         pipeline => ps_cNOb
         CALL PS%newField(ps_cNOb%shortName, "Cloud B bulk NO mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cNOb%shortName), pipeline)
            
         pipeline => NULL()
         ps_pNOa = FloatArray1d("pNOa")
         ps_pNOa%onDemand => globalMeanProfile
         pipeline => ps_pNOa
         CALL PS%newField(ps_pNOa%shortName, "Precipitation bulk NO mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_pNOa%shortName), pipeline)
      END IF

      IF (level == 5) THEN
         pipeline => NULL()
         ps_iNOa = FloatArray1d("iNOa")
         ps_iNOa%onDemand => globalMeanProfile
         pipeline => ps_iNOa
         CALL PS%newField(ps_iNOa%shortName, "Ice bulk NO mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_iNOa%shortName), pipeline)
      END IF

      IF (level >= 4) THEN

         pipeline => NULL()
         ps_aNHa = FloatArray1d("aNHa")
         ps_aNHa%onDemand => globalMeanProfile
         pipeline => ps_aNHa
         CALL PS%newField(ps_aNHa%shortName, "Aerosol A bulk NH mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aNHa%shortName), pipeline)
           
         pipeline => NULL()
         ps_aNHb = FloatArray1d("aNHb")
         ps_aNHb%onDemand => globalMeanProfile
         pipeline => ps_aNHb
         CALL PS%newField(ps_aNHb%shortName, "Aerosol B bulk NH mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_aNHb%shortName), pipeline)
                  
         pipeline => NULL()
         ps_cNHa = FloatArray1d("cNHa")
         ps_cNHa%onDemand => globalMeanProfile
         pipeline => ps_cNHa
         CALL PS%newField(ps_cNHa%shortName, "Cloud A bulk NH mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cNHa%shortName), pipeline)
                 
         pipeline => NULL()
         ps_cNHb = FloatArray1d("cNHb")
         ps_cNHb%onDemand => globalMeanProfile
         pipeline => ps_cNHb
         CALL PS%newField(ps_cNHb%shortName, "Cloud B bulk NH mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_cNHb%shortName), pipeline)
                 
         pipeline => NULL()
         ps_pNHa = FloatArray1d("pNHa")
         ps_pNHa%onDemand => globalMeanProfile
         pipeline => ps_pNHa
         CALL PS%newField(ps_pNHa%shortName, "Precipitation bulk NH mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_pNHa%shortName), pipeline)
      END IF

      IF (level == 5) THEN
         pipeline => NULL()
         ps_iNHa = FloatArray1d("iNHa")
         ps_iNHa%onDemand => globalMeanProfile
         pipeline => ps_iNHa
         CALL PS%newField(ps_iNHa%shortName, "Ice bulk NH mass", "kg/kg", "ztt",   &
                          ANY(outputlist == ps_iNHa%shortName), pipeline)
      END IF

      !---------------------------------------
      ! In-cloud mean profiles
      !---------------------------------------

      pipeline => NULL()
      psic_rc = FloatArray1d("ic_rc",srcName="rc")
      psic_rc%onDemand => inLiqMeanProfile
      pipeline => psic_rc
      CALL PS%newField(psic_rc%shortName, "In-cloud cloud water mix rat", "kg/kg", "ztt",   &
                       ANY(outputlist == psic_rc%shortName), pipeline)

      pipeline => NULL()
      psic_rh = FloatArray1d("ic_rh",srcName="rh")
      psic_rh%onDemand => inCloudMeanProfile
      pipeline => psic_rh
      CALL PS%newField(psic_rh%shortName, "In-cloud relative humidity", "1", "ztt",   &
                       ANY(outputlist == psic_rh%shortName), pipeline)
      
      IF (level >= 4) THEN
         pipeline => NULL()
         psic_Naa = FloatArray1d("ic_Naa",srcName="Naa")
         psic_Naa%onDemand => inCloudMeanProfile 
         pipeline => psic_Naa
         CALL PS%newField(psic_Naa%shortName, "Interstitial aerosol A", "kg-1", "ztt",   &
                          ANY(outputlist == psic_Naa%shortName), pipeline)

         pipeline => NULL()
         psic_Nab = FloatArray1d("ic_Nab",srcName="Nab")
         psic_Nab%onDemand => inCloudMeanProfile
         pipeline => psic_Nab
         CALL PS%newField(psic_Nab%shortName, "Interstitial aerosol B", "kg-1", "ztt",   &
                          ANY(outputlist == psic_Nab%shortName), pipeline)
         
         pipeline => NULL()
         psic_Nca = FloatArray1d("ic_Nca",srcName="Nca")
         psic_Nca%onDemand => incloudMeanProfile
         pipeline => psic_Nca
         CALL PS%newField(psic_Nca%shortName, "Cloud droplets A in-cloud", "kg-1", "ztt",   &
                          ANY(outputlist == psic_Nca%shortName), pipeline)

         pipeline => NULL()
         psic_Ncb = FloatArray1d("ic_Ncb",srcName="Ncb")
         psic_Ncb%onDemand => incloudMeanProfile
         pipeline => psic_Ncb
         CALL PS%newField(psic_Ncb%shortName, "Cloud droplets B in-cloud", "kg-1", "ztt",   &
                          ANY(outputlist == psic_Ncb%shortName), pipeline)

         pipeline => NULL()
         psic_CDNC = FloatArray1d("ic_CDNC",srcName="CDNC")
         psic_CDNC%onDemand => inLiqMeanProfile
         pipeline => psic_CDNC
         CALL PS%newField(psic_CDNC%shortName, "Cloud droplet number concentration in-cloud", "m-3", "ztt",   &
                          ANY(outputlist == psic_CDNC%shortName), pipeline)

         pipeline => NULL()
         psic_CNC = FloatArray1d("ic_CNC",srcName="CNC")
         psic_CNC%onDemand => inLiqMeanProfile
         pipeline => psic_CNC
         CALL PS%newField(psic_CNC%shortName, "Hydrometeor number concentration in-cloud", "m-3", "ztt",   &
                          ANY(outputlist == psic_CNC%shortName), pipeline)        

         pipeline => NULL()
         psic_Reff = FloatArray1d("ic_Reff",srcName="Reff")
         psic_Reff%onDemand => inLiqMeanProfile
         pipeline => psic_Reff
         CALL PS%newField(psic_Reff%shortName, "Droplet effective radius, in-cloud", "m", "ztt",   &
                          ANY(outputlist == psic_Reff%shortName), pipeline)          
         
      END IF

      !---------------------------------------
      ! precipitating grid-box mean profiles
      !---------------------------------------

      pipeline => NULL()
      pspr_rrate = FloatArray1d("pr_rrate",srcName="rrate")
      pspr_rrate%onDemand => precipMeanProfile  
      pipeline => pspr_rrate
      CALL PS%newField(pspr_rrate%shortName, "Conditionally avg precipitation flux", "W m-2", "ztt",   &
                       ANY(outputlist == pspr_rrate%shortName), pipeline)   

      IF (level == 3) THEN
         pipeline => NULL()
         pspr_rpp = FloatArray1d("pr_rpp",srcName="rpp")
         pspr_rpp%onDemand => precipMeanProfile
         pipeline => pspr_rpp
         CALL PS%newField(pspr_rpp%shortName, "Conditionally avg precip mix rat", "kg kg-1", "ztt",   &
                          ANY(outputlist == pspr_rpp%shortName), pipeline)   
         
         pipeline => NULL()
         pspr_npp = FloatArray1d("pr_npp",srcName="npp")
         pspr_npp%onDemand => precipMeanProfile
         pipeline => pspr_npp
         CALL PS%newField(pspr_npp%shortName, "Conditionally avg precip number", "kg-1", "ztt",   &
                          ANY(outputlist == pspr_npp%shortName), pipeline)   

      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         pspr_srp = FloatArray1d("pr_srp",srcName="srp")
         pspr_srp%onDemand => precipMeanProfile
         pipeline => pspr_srp
         CALL PS%newField(pspr_srp%shortName, "Conditionally avg precip mix rat", "kg kg-1", "ztt",   &
                          ANY(outputlist == pspr_srp%shortName), pipeline)

         pipeline => NULL()
         pspr_Np = FloatArray1d("pr_Np",srcName="Np")
         pspr_Np%onDemand => precipMeanProfile
         pipeline => pspr_Np
         CALL PS%newField(pspr_Np%shortName, "Conditionally avg precip number", "kg-1", "ztt",   &
                          ANY(outputlist == pspr_Np%shortName), pipeline)           

      END IF


      ! ---------------------------------------------
      ! Scalar and momentum fluxes
      ! ---------------------------------------------
      pipeline => NULL()
      ps_tw_res = FloatArray1d("tw_res")
      ps_tw_res%onDemand => meanScalarFlux
      pipeline => ps_tw_res
      CALL PS%newField(ps_tw_res%shortName, "Resolved vertical flux of theta", "W m-2", "ztt",   &
                       ANY(outputlist == ps_tw_res%shortName), pipeline)

      pipeline => NULL()
      ps_qw_res = FloatArray1d("qw_res")
      ps_qw_res%onDemand => meanScalarFlux
      pipeline => ps_qw_res
      CALL PS%newField(ps_qw_res%shortName, "Resolved vertical flux of moisture", "W m-2", "ztt",  &
                       ANY(outputlist == ps_qw_res%shortName), pipeline)

      pipeline => NULL()
      ps_lw_res = FloatArray1d("lw_res")
      ps_lw_res%onDemand => meanScalarFlux
      pipeline => ps_lw_res
      CALL PS%newField(ps_lw_res%shortName, "Resolved vertical flux of liquid water", "W m-2", "ztt",  &
                       ANY(outputlist == ps_lw_res%shortName), pipeline)
      
      pipeline => NULL()
      ps_uw_res = FloatArray1d("uw_res")
      ps_uw_res%onDemand => meanMomFlux
      pipeline => ps_uw_res
      CALL PS%newField(ps_uw_res%shortName, "Resolved flux of u momentum", "m2 s-2", "ztt",   &
                       ANY(outputlist == ps_uw_res%shortName), pipeline)

      pipeline => NULL()
      ps_vw_res = FloatArray1d("vw_res")
      ps_vw_res%onDemand => meanMomFlux
      pipeline => ps_vw_res
      CALL PS%newField(ps_vw_res%shortName, "Resolved flux of v momentum", "m2 s-2", "ztt",   &
                       ANY(outputlist == ps_vw_res%shortName), pipeline)

      pipeline => NULL()
      ps_ww_res = FloatArray1d("ww_res")
      ps_ww_res%onDemand => meanMomFlux
      pipeline => ps_ww_res
      CALL PS%newField(ps_ww_res%shortName, "Resolved flux of w momentum", "m2 s-2", "ztt",   &
                       ANY(outputlist == ps_ww_res%shortName), pipeline)

      
      ! --------------------------------------------------------
      ! Turbulence statistics
      ! --------------------------------------------------------
      pipeline => NULL()
      ps_wmin = FloatArray1d("wmin",srcname="wwind")
      ps_wmin%onDemand => globalMinProfile
      pipeline => ps_wmin
      CALL PS%newField(ps_wmin%shortName, "Minimum vertical velocity", "m s-1", "ztt",    &
                       ANY(outputlist == ps_wmin%shortName), pipeline)

      pipeline => NULL()
      ps_wmax = FloatArray1d("wmax",srcname="wwind")
      ps_wmax%onDemand => globalMaxProfile
      pipeline => ps_wmax
      CALL PS%newField(ps_wmax%shortName, "Maximum vertical velocity", "m s-1", "ztt",    &
                       ANY(outputlist == ps_wmax%shortName), pipeline)

      pipeline => NULL()
      ps_wbar = FloatArray1d("wbar",srcname="wwind")
      ps_wbar%onDemand => globalMeanProfile
      pipeline => ps_wbar
      CALL PS%newField(ps_wbar%shortName, "Mean vertical velocity", "m s-1", "ztt",     &
                       ANY(outputlist == ps_wbar%shortName), pipeline)

      pipeline => NULL()
      ps_ubar = FloatArray1d("ubar",srcname="uwin")
      ps_ubar%onDemand => globalMeanProfile
      pipeline => ps_ubar
      CALL PS%newField(ps_ubar%shortName, "Mean U wind", "m s-1", "ztt",     &
                       ANY(outputlist == ps_ubar%shortName), pipeline)

      pipeline => NULL()
      ps_vbar = FloatArray1d("vbar",srcname="vwind")
      ps_vbar%onDemand => globalMeanProfile
      pipeline => ps_vbar
      CALL PS%newField(ps_vbar%shortName, "Mean V wind", "m s-1", "ztt",     &
                       ANY(outputlist == ps_vbar%shortName), pipeline)      
      
      pipeline => NULL()
      ps_w_2 = FloatArray1d("w_2",srcname="wwind")
      ps_w_2%onDemand => globalVarProfile
      pipeline => ps_w_2
      CALL PS%newField(ps_w_2%shortName, "Vertical velocity variance", "m s-1", "ztt",    &
                       ANY(outputlist == ps_w_2%shortName), pipeline)

      pipeline => NULL()
      ps_u_2 = FloatArray1d("u_2",srcname="uwind")
      ps_u_2%onDemand => globalVarProfile
      pipeline => ps_u_2
      CALL PS%newField(ps_u_2%shortName, "U wind variance", "m s-1", "ztt",    &
                       ANY(outputlist == ps_u_2%shortName), pipeline)      

      pipeline => NULL()
      ps_v_2 = FloatArray1d("v_2",srcname="vwind")
      ps_v_2%onDemand => globalVarProfile
      pipeline => ps_v_2
      CALL PS%newField(ps_v_2%shortName, "V wind variance", "m s-1", "ztt",    &
                       ANY(outputlist == ps_v_2%shortName), pipeline)      

      pipeline => NULL()
      ps_tke_res = FloatArray1d("tke_res")
      ps_tke_res%onDemand => meanTKEres
      pipeline => ps_tke_res
      CALL PS%newField(ps_tke_res%shortName, "TKE from resolved flow", "m2 s-2", "ztt",  &
                       ANY(outputlist == ps_tke_res%shortName), pipeline)

      pipeline => NULL()

      ! ---------------------------------------------------------------
      ! Piggybacking slave microphysical diagnostics
      ! ---------------------------------------------------------------
      IF (lpback) THEN
         pipeline => NULL()
         ps_pb_rrate = FloatArray1d("pb_rrate")
         ps_pb_rrate%onDemand => globalMeanProfile
         pipeline => ps_pb_rrate
         CALL PS%newField(ps_pb_rrate%shortName, "Precip flux (bulk slave)", "W m-2", "ztt",   &
                          ANY(outputlist == ps_pb_rrate%shortName), pipeline)

         pipeline => NULL()
         pspr_pb_rrate = FloatArray1d("pr_pb_rrate",srcname="pb_rrate")
         pspr_pb_rrate%onDemand => precipMeanProfile
         pipeline => pspr_pb_rrate
         CALL PS%newField(pspr_pb_rrate%shortName, "Conditional precip flux (bulk slave)", "W m-2", "ztt",  &
                          ANY(outputlist == pspr_pb_rrate%shortName), pipeline)

         pipeline => NULL()
         ps_pb_rc = FloatArray1d("pb_rc")
         ps_pb_rc%onDemand => globalMeanProfile
         pipeline => ps_pb_rc
         CALL PS%newField(ps_pb_rc%shortName, "Cloud condensate mixing ratio (bulk slave)", "W m-2", "ztt",  &
                          ANY(outputlist == ps_pb_rc%shortName), pipeline)

         pipeline => NULL()
         psic_pb_rc = FloatArray1d("ic_pb_rc",srcname="pb_rc")
         psic_pb_rc%onDemand => inLiqMeanProfile
         pipeline => psic_pb_rc
         CALL PS%newField(psic_pb_rc%shortName, "Conditional in-cloud condensate mix rat (bulk slave)",      &
                          "W m-2", "ztt", ANY(outputlist == psic_pb_rc%shortName), pipeline)                  

      END IF

      ! --------------------------------------------------------------
      ! Process rate diagnostics
      ! --------------------------------------------------------------
      IF (level <= 3 .OR. lpback) THEN
         pipeline => NULL()
         psic_b_m_autoc = FloatArray1d("ic_b_m_autoc",srcname="b_m_autoc")
         psic_b_m_autoc%onDemand => inCloudMeanProfile
         pipeline => psic_b_m_autoc
         CALL PS%newField(psic_b_m_autoc%shortName, "Conditional avg autoconversion rate (bulk scheme, mass)",  &
                          "kg/kg s", "ztt", ANY(outputlist == psic_b_m_autoc%shortName), pipeline)

         pipeline => NULL()
         psic_b_n_autoc = FloatArray1d("ic_b_n_autoc",srcname="b_n_autoc")
         psic_b_n_autoc%onDemand => inCloudMeanProfile
         pipeline => psic_b_n_autoc
         CALL PS%newField(psic_b_n_autoc%shortName, "Conditional avg autoconversion rate (bulk scheme, number)", &
                          "#/kg s", "ztt", ANY(outputlist == psic_b_n_autoc%shortName), pipeline)

         pipeline => NULL()
         psic_b_m_accr = FloatArray1d("ic_b_m_accr",srcname="b_m_accr")
         psic_b_m_accr%onDemand => inCloudMeanProfile
         pipeline => psic_b_m_accr
         CALL PS%newField(psic_b_m_accr%shortName, "Conditional avg accretion rate (bulk scheme, mass)",    &
                          "kg/kg s", "ztt", ANY(outputlist == psic_b_m_accr%shortName), pipeline)
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         psic_s_m_autoc = FloatArray1d("ic_s_m_autoc",srcname="s_m_autoc")
         psic_s_m_autoc%onDemand => inCloudMeanProfile
         pipeline => psic_s_m_autoc
         CALL PS%newField(psic_s_m_autoc%shortName, "Conditional avg autoconversion rate (SALSA, mass)",   &
                          "kg/kg s", "ztt", ANY(outputlist == psic_s_m_autoc%shortName), pipeline)

         pipeline => NULL()
         psic_s_m_autoc50 = FloatArray1d("ic_s_m_autoc50",srcname="s_m_autoc50")
         psic_s_m_autoc50%onDemand => inCloudMeanProfile
         pipeline => psic_s_m_autoc50
         CALL PS%newField(psic_s_m_autoc50%shortName,      &
                          "Conditional avg autoconversion rate (SALSA, mass, over 50um)",   &
                          "kg/kg s", "ztt", ANY(outputlist == psic_s_m_autoc50%shortName), pipeline)  
         
         pipeline => NULL()
         psic_s_m_autoc80 = FloatArray1d("ic_s_m_autoc80",srcname="s_m_autoc80")
         psic_s_m_autoc80%onDemand => inCloudMeanProfile
         pipeline => psic_s_m_autoc80
         CALL PS%newField(psic_s_m_autoc80%shortName,      &
                          "Conditional avg autoconversion rate (SALSA, mass, over 80um)",   &
                          "kg/kg s", "ztt", ANY(outputlist == psic_s_m_autoc80%shortName), pipeline)         
         
         pipeline => NULL()         
         psic_s_m_accr = FloatArray1d("ic_s_m_accr",srcname="s_m_accr")
         psic_s_m_accr%onDemand => inCloudMeanProfile
         pipeline => psic_s_m_accr
         CALL PS%newField(psic_s_m_accr%shortName, "Conditional avg accretion rate (SALSA, mass)",    &
                          "kg/kg s", "ztt", ANY(outputlist == psic_s_m_accr%shortName), pipeline)

         pipeline => NULL()         
         psic_s_m_accr50 = FloatArray1d("ic_s_m_accr50",srcname="s_m_accr50")
         psic_s_m_accr50%onDemand => inCloudMeanProfile
         pipeline => psic_s_m_accr50
         CALL PS%newField(psic_s_m_accr50%shortName,       &
                          "Conditional avg accretion rate (SALSA, mass, over 50um)",    &
                          "kg/kg s", "ztt", ANY(outputlist == psic_s_m_accr50%shortName), pipeline)
         
         pipeline => NULL()         
         psic_s_m_accr80 = FloatArray1d("ic_s_m_accr80",srcname="s_m_accr80")
         psic_s_m_accr80%onDemand => inCloudMeanProfile
         pipeline => psic_s_m_accr80
         CALL PS%newField(psic_s_m_accr80%shortName,       &
                          "Conditional avg accretion rate (SALSA, mass, over 80um)",    &
                          "kg/kg s", "ztt", ANY(outputlist == psic_s_m_accr80%shortName), pipeline)
         
         pipeline => NULL()
         psic_s_m_ACcoll_dry = FloatArray1d("ic_s_m_ACcoll_dry",srcname="s_m_ACcoll_dry")
         psic_s_m_ACcoll_dry%onDemand => inCloudMeanProfile
         pipeline => psic_s_m_ACcoll_dry
         CALL PS%newField(psic_s_m_ACcoll_dry%shortName,   &
                          "Conditional avg cloud-aerosol collection rate (SALSA, mass)",    &
                          "kg/kg s", "ztt", ANY(outputlist == psic_s_m_ACcoll_dry%shortName), pipeline)
                  
         pipeline => NULL()
         psic_s_m_APcoll_dry = FloatArray1d("ic_s_m_APcoll_dry",srcname="s_m_APcoll_dry")
         psic_s_m_APcoll_dry%onDemand => inCloudMeanProfile
         pipeline => psic_s_m_APcoll_dry
         CALL PS%newField(psic_s_m_APcoll_dry%shortName,   &
                          "Conditional avg precip-aerosol collection rate (SALSA, mass)", "kg/kg s",  &
                          "ztt", ANY(outputlist == psic_s_m_APcoll_dry%shortName), pipeline)
                  
         pipeline => NULL()
         psic_s_m_AIcoll_dry = FloatArray1d("ic_s_m_AIcoll_dry",srcname="s_m_AIcoll_dry")
         psic_s_m_AIcoll_dry%onDemand => inCloudMeanProfile   ! Should implement in-ice averaging
         pipeline => psic_s_m_AIcoll_dry
         CALL PS%newField(psic_s_m_AIcoll_dry%shortName,   &
                          "Conditional avg ice-aerosol collection rate (SALSA, mass)", "kg/kg s",  &
                          "ztt", ANY(outputlist == psic_s_m_AIcoll_dry%shortName), pipeline)
                  
         pipeline => NULL()
         psic_s_n_activ = FloatArray1d("ic_s_n_activ",srcname="s_n_activ")
         psic_s_n_activ%onDemand => inCloudMeanProfile
         pipeline => psic_s_n_activ
         CALL PS%newField(psic_s_n_activ%shortName,   &
                          "Conditional avg activation rate (SALSA, num)", "#/kg s",  &
                          "ztt", ANY(outputlist == psic_s_n_activ%shortName), pipeline)
                  
         pipeline => NULL()
         psic_s_n_icehom = FloatArray1d("ic_s_n_icehom",srcname="s_n_icehom")
         psic_s_n_icehom%onDemand => inCloudMeanProfile  ! Should implement in-ice averaging
         pipeline => psic_s_n_icehom
         CALL PS%newField(psic_s_n_icehom%shortName,  &
                          "Conditional avg homog ice nucl rate (SALSA, num)", "#/kg s",  &
                          "ztt", ANY(outputlist == psic_s_n_icehom%shortName), pipeline)
                  
         pipeline => NULL()
         psic_s_n_icedep = FloatArray1d("ic_s_n_icedep",srcname="s_n_icedep")
         psic_s_n_icedep%onDemand => inCloudMeanProfile  ! Should implement in-ice averaging
         pipeline => psic_s_n_icedep
         CALL PS%newField(psic_s_n_icedep%shortName,  &
                          "Conditional avg deposition ice nucl rate (SALSA, num)", "#/kg s",  &
                          "ztt", ANY(outputlist == psic_s_n_icedep%shortName), pipeline)
                  
         pipeline => NULL()
         psic_s_n_iceimm = FloatArray1d("ic_s_n_iceimm",srcname="s_n_iceimm")
         psic_s_n_iceimm%onDemand => inCloudMeanProfile  ! Should implement in-ice averaging
         pipeline => psic_s_n_iceimm
         CALL PS%newField(psic_s_n_iceimm%shortName,  &
                          "Conditional avg immersion ice nucl rate (SALSA, num)", "#/kg s",  &
                          "ztt", ANY(outputlist == psic_s_n_iceimm%shortName), pipeline)
                  
         pipeline => NULL()
         psic_s_m_conda = FloatArray1d("ic_s_m_conda",srcname="s_m_conda")
         psic_s_m_conda%onDemand => inCloudMeanProfile 
         pipeline => psic_s_m_conda
         CALL PS%newField(psic_s_m_conda%shortName,   &
                          "Conditional avg condensation on aero (SALSA, mass)", "kg/kg s",  &
                          "ztt", ANY(outputlist == psic_s_m_conda%shortName), pipeline)
                  
         pipeline => NULL()
         psic_s_m_condc = FloatArray1d("ic_s_m_condc",srcname="s_m_condc")
         psic_s_m_condc%onDemand => inCloudMeanProfile
         pipeline => psic_s_m_condc
         CALL PS%newField(psic_s_m_condc%shortName,   &
                          "Conditional avg condensation on cloud (SALSA, mass)", "kg/kg s",  &
                          "ztt", ANY(outputlist == psic_s_m_condc%shortName), pipeline)
                  
         pipeline => NULL()
         psic_s_m_condp = FloatArray1d("ic_s_m_condp",srcname="s_m_condp")
         psic_s_m_condp%onDemand => inCloudMeanProfile
         pipeline => psic_s_m_condp
         CALL PS%newField(psic_s_m_condp%shortName,   &
                          "Conditional avg condensation on precip (SALSA, mass)", "kg/kg s",  &
                          "ztt", ANY(outputlist == psic_s_m_condp%shortName), pipeline)
                  
         pipeline => NULL()
         psic_s_m_condi = FloatArray1d("ic_s_m_condi",srcname="s_m_condi")
         psic_s_m_condi%onDemand => inCloudMeanProfile  ! Should implement in-ice averaging
         pipeline => psic_s_m_condi
         CALL PS%newField(psic_s_m_condi%shortName,   &
                          "Conditional avg condensation on ice (SALSA, mass)", "kg/kg s",  &
                          "ztt", ANY(outputlist == psic_s_m_condi%shortName), pipeline)
         
      END IF
      
    END SUBROUTINE setPSVariables


END MODULE mo_ps_state
