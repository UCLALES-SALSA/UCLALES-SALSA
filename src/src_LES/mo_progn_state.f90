MODULE mo_progn_state
  USE classFieldArray
  USE mo_structured_datatypes, ONLY : FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d
  USE mo_submctl, ONLY : spec, nbins, ncld, nprc, nice, in1a,fn2a, in2b,fn2b, ica,fca, icb,fcb, ice_theta_dist
  IMPLICIT NONE

  SAVE
  
  ! prognostic scalar variables
  ! ===============================
  !
  ! General
  ! --------------------------------------------------------------------------------------------------------------------------
  TYPE(FloatArray3D), TARGET :: a_tp,a_tt   ! Liquid potential temperature - defined as deviation from a mean state
  TYPE(FloatArray3D), TARGET :: a_rp,a_rt   ! With level < 4 this is the TOTAL water content.
                                                          ! With level >= 4 this is the water VAPOUR content
  TYPE(FloatArray3D), TARGET :: a_rpp,a_rpt ! Precipitation mixing ratio for level < 4
  TYPE(FloatArray3D), TARGET :: a_npp,a_npt ! Precipitation mnumber for level < 4
  TYPE(FloatArray3D), TARGET :: a_qp,a_qt   ! 
  
  ! SALSA tracers
  !---------------------------------------------------------------------------
  ! -- Masses given in kg/kg, number concentrations in #/kg
  ! -- Each size bin/species will be treated as a separate tracer.
  ! -- NOTE: Mass mixing ratio  arrays are reduced to 4 dims.
  !          The 4th dim contains all the size bins sequentially for
  !          each aerosol species  + water + for ice also rimed ice
  !
  !          Gas tracers are contained sequentially in dimension
  !          4 as: 1. SO4, 2. HNO3, 3. NH3, 4. OCNV, 5. OCSV

  ! =================================
  ! NOTE that the SALSA tracers here cannot be used directly as output. Use the corresponding variables
  ! separated for A and B regimes in mo_derived_state.f90
  ! =================================
  
  ! -- Number concentrations
  TYPE(FloatArray4D), TARGET :: a_naerop,  a_naerot,  &    ! Number of aerosol
                                a_ncloudp, a_ncloudt, &    ! Number of cloud
                                a_nprecpp, a_nprecpt, &    ! Number of precip
                                a_nicep,   a_nicet         ! Number of ice
  
  !                      ! A and B bins for output
  TYPE(FloatArray4d), TARGET :: a_Naba, a_Nabb, a_Ncba, a_Ncbb  ! Aerosol A, B, Cloud droplets A,B
  TYPE(FloatArray4d), TARGET :: a_Npba, a_Niba                  ! Precipitation and ice; these are not really neede but added
                                                                ! anyway to keep consistent output variable naming convention
  
  ! -- Mass concentrations
  TYPE(FloatArray4D), TARGET :: a_maerop,  a_maerot,  &    ! Aerosol
                                a_mcloudp, a_mcloudt, &    ! Cloud
                                a_mprecpp, a_mprecpt, &    ! Precip
                                a_micep,   a_micet         ! Ice

  TYPE(FloatArray4D), TARGET :: a_indefp, a_indeft
  
  ! -- Gas compound tracers
  TYPE(FloatArray4D), TARGET :: a_gaerop, a_gaerot
  
  CONTAINS
  
    SUBROUTINE setPrognosticVariables(a_sclrp,a_sclrt,Prog,outputlist,level,isgstyp,lbinanl,nzp,nxp,nyp,nscl)
      INTEGER, INTENT(in) :: level,isgstyp,nzp,nxp,nyp,nscl
      LOGICAL, INTENT(in) :: lbinanl
      REAL, INTENT(in) :: a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl)
      CHARACTER(len=*), INTENT(in) :: outputlist(:)
      TYPE(FieldArray), INTENT(inout) :: Prog
 
      CLASS(*), POINTER :: pipeline_p => NULL(), pipeline_t => NULL()

      INTEGER :: iscl, nspec
      
      nspec = 0
      IF (level >= 4) nspec = spec%getNSpec(type="wet")  ! Number of aerosol compounds used including water (but not rimed ice!!)
      
      iscl = 0

      iscl = iscl + 1
      a_tp = FloatArray3D(a_sclrp(:,:,:,iscl))
      a_tt = FloatArray3D(a_sclrt(:,:,:,iscl))
      pipeline_p => a_tp
      pipeline_t => a_tt
      CALL Prog%newField( "tp","Liquid water potential temperature, deviation from the mean",    &
                          "K", "tttt", ANY(outputlist == "theta_l"),                                  &
                          pipeline_p, in_t_data = pipeline_t                                          &
                        )
      
      iscl = iscl + 1
      pipeline_p => NULL(); pipeline_t => NULL()
      a_rp = FloatArray3D(a_sclrp(:,:,:,iscl))
      a_rt = FloatArray3D(a_sclrt(:,:,:,iscl))
      pipeline_p => a_rp
      pipeline_t => a_rt
      CALL Prog%newField( "rp","Progostic water vapor (lev4+) or total water content(lev3-)",   &
                          "kg/kg", "tttt", ANY(outputlist == "rp"),                             &
                          pipeline_p, in_t_data = pipeline_t                                    &
                        )
      

      IF (level < 4) THEN
         iscl = iscl + 1
         pipeline_p => NULL(); pipeline_t => NULL()
         a_rpp = FloatArray3D(a_sclrp(:,:,:,iscl))
         a_rpt = FloatArray3D(a_sclrt(:,:,:,iscl))
         pipeline_p => a_rpp
         pipeline_t => a_rpt
         CALL Prog%newField( "rpp","Precipitation mixing ratio",               &
                             "kg/kg", "tttt", ANY(outputlist == "rpp"),        &
                             pipeline_p, in_t_data = pipeline_t                &
                           )
         
         iscl = iscl + 1
         pipeline_p => NULL(); pipeline_t => NULL()
         a_npp = FloatArray3D(a_sclrp(:,:,:,iscl))
         a_npt = FloatArray3D(a_sclrt(:,:,:,iscl))
         pipeline_p => a_npp
         pipeline_t => a_npt
         CALL Prog%newField( "npp","Precipitation number",             &
                             "#/kg", "tttt", ANY(outputlist == "npp"), &
                             pipeline_p, in_t_data = pipeline_t        &
                           )                              
      END IF
                         
      IF (isgstyp > 1) THEN
         iscl = iscl + 1
         pipeline_p => NULL(); pipeline_t => NULL()
         a_qp = FloatArray3D(a_sclrp(:,:,:,iscl))
         a_qt = FloatArray3D(a_sclrt(:,:,:,iscl))
         pipeline_p => a_qp
         pipeline_t => a_qt 
         CALL Prog%newField( "qp","CHECK",                            &
                             "CHECK", "tttt", ANY(outputlist == "qp"), &
                             pipeline_p, in_t_data = pipeline_t       &
                           )    
      END IF

      
      ! SALSA Variables
      ! ------------------------------------------------------
      IF (level >= 4) THEN
         
         iscl = iscl + 1
         pipeline_p => NULL(); pipeline_t => NULL()
         a_naerop = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+nbins-1))
         a_naerot = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+nbins-1))
         pipeline_p => a_naerop
         pipeline_t => a_naerot
         CALL Prog%newField( "naero","Binned aerosol number concentration",    &
                             "#/kg", "N/A", .FALSE.,                           &
                             pipeline_p, in_t_data = pipeline_t,               &
                             in_group = ["SALSA_4d"]                           &
                             )
         pipeline_p => NULL()
         a_Naba = FloatArray4d(a_naerop%d(:,:,:,in1a:fn2a))
         pipeline_p => a_Naba
         CALL Prog%newField( "Naba", "Binned aerosol number, A",                           &
                             "#/kg", "ttttaea", ANY(outputlist == "Naba" .AND. lbinanl),   &
                             pipeline_p                                                    &
                             )
         pipeline_p => NULL()
         a_Nabb = FloatArray4d(a_naerop%d(:,:,:,in2b:fn2b))
         pipeline_p => a_Nabb
         CALL Prog%newField( "Nabb", "Binned aerosol number, B",                           & 
                             "#/kg", "ttttaeb", ANY(outputlist == "Nabb" .AND. lbinanl),   &
                             pipeline_p                                                    &
                           )

         iscl = iscl + nbins
         pipeline_p => NULL(); pipeline_t => NULL()
         a_maerop = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+(nbins*nspec)-1))
         a_maerot = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+(nbins*nspec)-1))
         pipeline_p => a_maerop
         pipeline_t => a_maerot
         CALL Prog%newField( "maero","Binned aerosol mass concentration, sequentially for each compound",    &
                             "kg/kg", "N/A", .FALSE.,                                                        &
                             pipeline_p, in_t_data = pipeline_t,                                             &
                             in_group = ["SALSA_4d"]                                                         &
                           )    
         
         iscl = iscl + nbins*nspec
         pipeline_p => NULL(); pipeline_t => NULL()
         a_ncloudp = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+ncld-1))
         a_ncloudt = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+ncld-1))
         pipeline_p => a_ncloudp
         pipeline_t => a_ncloudt
         CALL Prog%newField( "ncloud","Binned cloud number concentration",    &
                             "#/kg", "N/A", .FALSE.,                          &
                             pipeline_p, in_t_data = pipeline_t,              &
                             in_group = ["SALSA_4d"]                          &
                             )
         pipeline_p => NULL()
         a_Ncba = FloatArray4d(a_ncloudp%d(:,:,:,ica%cur:fca%cur))
         pipeline_p => a_Ncba
         CALL Prog%newField( "Ncba", "Binned cloud number, A",                             & 
                             "#/kg", "ttttcla", ANY(outputlist == "Ncba" .AND. lbinanl),   &
                             pipeline_p                                                    &
                             )
         pipeline_p => NULL()
         a_Ncbb = FloatArray4d(a_ncloudp%d(:,:,:,icb%cur:fcb%cur))
         pipeline_p => a_Ncbb
         CALL Prog%newField( "Ncbb", "Binned cloud number, B",                             &
                             "#/kg", "ttttclb", ANY(outputlist == "Ncbb" .AND. lbinanl),   &
                             pipeline_p                                                    &
                           )

         iscl = iscl + ncld
         pipeline_p => NULL(); pipeline_t => NULL()
         a_mcloudp = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+(ncld*nspec)-1))
         a_mcloudt = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+(ncld*nspec)-1))
         pipeline_p => a_mcloudp
         pipeline_t => a_mcloudt
         CALL Prog%newField( "mcloud","Binned cloud mass concentration, sequentially for each compound",    &
                             "kg/kg", "N/A", .FALSE.,                                                       &
                             pipeline_p, in_t_data = pipeline_t,                                            &
                             in_group = ["SALSA_4d"]                                                        &
                           )    
         
         iscl = iscl + ncld*nspec
         pipeline_p => NULL(); pipeline_t => NULL()
         a_nprecpp = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+nprc-1))
         a_nprecpt = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+nprc-1))
         pipeline_p => a_nprecpp
         pipeline_t => a_nprecpt
         CALL Prog%newField( "nprecp","Binned precip number concentration",     &
                             "#/kg", "N/A", .FALSE.,                            &  
                             pipeline_p, in_t_data = pipeline_t,                &
                             in_group = ["SALSA_4d"]                            &
                             )
         pipeline_p => NULL()
         a_Npba = FloatArray4d(a_nprecpp%d(:,:,:,1:nprc))
         pipeline_p => a_Npba
         CALL Prog%newField( "Npba", "Binned precipitation number",                         &
                             "#/kg", "ttttprc", ANY(outputlist == "Npba" .AND. lbinanl),    &
                             pipeline_p                                                     &
                           )
         
         iscl = iscl + nprc
         pipeline_p => NULL(); pipeline_t => NULL()
         a_mprecpp = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+(nprc*nspec)-1))
         a_mprecpt = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+(nprc*nspec)-1))
         pipeline_p => a_mprecpp
         pipeline_t => a_mprecpt
         CALL Prog%newField( "mprecp","Binned precip mass concentration, sequentially for each compound",    &
                             "kg/kg", "N/A", .FALSE.,                                                        &
                             pipeline_p, in_t_data = pipeline_t,                                             &
                             in_group = ["SALSA_4d"]                                                         &
                           )                            
         
         iscl = iscl + nprc*nspec
         pipeline_p => NULL(); pipeline_t => NULL()
         a_gaerop = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+5-1))
         a_gaerot = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+5-1))
         pipeline_p => a_gaerop
         pipeline_t => a_gaerot
         CALL Prog%newField( "gaero","Binned ice number concentration",    &
                             "#/kg", "N/A", .FALSE.,                       &
                             pipeline_p, in_t_data = pipeline_t,           &
                             in_group = ["SALSA_4d"]                       &
                           )                            
         iscl = iscl + 5-1
         
      END IF

      IF (level == 5) THEN
         iscl = iscl + 1
         pipeline_p => NULL(); pipeline_t => NULL()
         a_nicep = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+nice-1))
         a_nicet = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+nice-1))
         pipeline_p => a_nicep
         pipeline_t => a_nicet
         CALL Prog%newField( "nice","Binned ice number concentration",      &
                             "#/kg", "N/A", .FALSE.,                        &
                             pipeline_p, in_t_data = pipeline_t,            &
                             in_group = ["SALSA_4d"]                        &
                             )
         pipeline_p => NULL()
         a_Niba = FloatArray4d(a_nicep%d(:,:,:,1:nice))
         pipeline_p => a_Niba
         CALL Prog%newField( "Niba", "Binned ice number",                                 &
                             "#/kg", "ttttice", ANY(outputlist == "Niba" .AND. lbinanl),  &
                             pipeline_p                                                   &
                           )
                  
         iscl = iscl + nice
         pipeline_p => NULL(); pipeline_t => NULL()
         a_micep = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+(nice*(nspec+1))-1)) ! nspec+1 for rimed ice!
         a_micet = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+(nice*(nspec+1))-1))
         pipeline_p => a_micep
         pipeline_t => a_micet
         CALL Prog%newField( "mice","Binned ice mass concentration, sequentially for each compound",    &
                             "kg/kg", "N/A", .FALSE.,                                                   &
                             pipeline_p, in_t_data = pipeline_t,                                        &
                             in_group = ["SALSA_4d"]                                                    &
                           )

         IF (ice_theta_dist) THEN
            iscl = iscl + nice*(nspec+1)
            pipeline_p => NULL(); pipeline_t => NULL()
            a_indefp = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+nbins+ncld+nprc-1)) 
            a_indeft = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+nbins+ncld+nprc-1))
            pipeline_p => a_indefp
            pipeline_t => a_indeft
            CALL Prog%newField( "indef","IN deficit fraction for contact angle distributions",    &
                                "kg/kg", "N/A", .FALSE.,     &
                                pipeline_p,in_t_data = pipeline_t                                 &
                              )
            
            
         END IF

            
      END IF
      

      pipeline_p => NULL()
      pipeline_t => NULL()
            
    END SUBROUTINE setPrognosticVariables

    


    
END MODULE mo_progn_state
