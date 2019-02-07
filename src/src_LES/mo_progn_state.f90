MODULE mo_progn_state
  USE classFieldArray
  USE mo_structured_datatypes, ONLY : FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d
  USE mo_submctl, ONLY : spec, nbins, ncld, nprc, nice
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
  TYPE(FloatArray3D), TARGET :: a_qp,a_qt   ! T
  
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
  
  ! Prognostic tracers
  ! -- Number concentrations
  TYPE(FloatArray4D), TARGET :: a_naerop,  a_naerot,  &    ! Number of aerosol
                                a_ncloudp, a_ncloudt, &    ! Number of cloud
                                a_nprecpp, a_nprecpt, &    ! Number of precip
                                a_nicep,   a_nicet         ! Number of ice

  ! -- Mass concentrations
  TYPE(FloatArray4D), TARGET :: a_maerop,  a_maerot,  &    ! Aerosol
                                a_mcloudp, a_mcloudt, &    ! Cloud
                                a_mprecpp, a_mprecpt, &    ! Precip
                                a_micep,   a_micet         ! Ice

  ! -- Gas compound tracers
  TYPE(FloatArray4D), TARGET :: a_gaerop, a_gaerot

  CONTAINS
  
    SUBROUTINE setPrognosticVariables(a_sclrp,a_sclrt,Prog,outputlist,memsize,level,isgstyp,lbinanl,nzp,nxp,nyp,nscl)
      INTEGER, INTENT(in) :: level,isgstyp,nzp,nxp,nyp,nscl
      LOGICAL, INTENT(in) :: lbinanl
      REAL, INTENT(in) :: a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl)
      CHARACTER(len=10), INTENT(in) :: outputlist(:)
      TYPE(FieldArray), INTENT(inout) :: Prog
      INTEGER, INTENT(inout) :: memsize
 
      CLASS(*), POINTER :: pipeline_p, pipeline_t

      INTEGER :: iscl, nspec
      
      IF (level >= 4) nspec = spec%getNSpec(type="wet")  ! Number of aerosol compounds used including water (but not rimed ice!!)
      
      iscl = 0

      iscl = iscl + 1
      a_tp = FloatArray3D(a_sclrp(:,:,:,iscl),store=.FALSE.)
      a_tt = FloatArray3D(a_sclrt(:,:,:,iscl),store=.FALSE.)
      pipeline_p => a_tp
      pipeline_t => a_tt
      CALL Prog%NewField( "theta_l","Liquid water potential temperature, deviation from the mean",    &
                          "K", "tttt", ANY(outputlist == "theta_l"),                                  &
                          pipeline_p, in_t_data = pipeline_t,                                         &
                          in_group = "LES"                                                            &
                        )

      iscl = iscl + 1
      a_rp = FloatArray3D(a_sclrp(:,:,:,iscl),store=.FALSE.)
      a_rt = FloatArray3D(a_sclrt(:,:,:,iscl),store=.FALSE.)
      pipeline_p => a_rp
      pipeline_t => a_rt
      CALL Prog%NewField( "rp","Progostic water vapor (lev4+) or total water content(lev3-)",   &
                          "kg/kg", "tttt", ANY(outputlist == "rp"),                             &
                          pipeline_p, in_t_data = pipeline_t,                                   &
                          in_group = "LES"                                                      &
                        )
      

      IF (level < 4) THEN
         iscl = iscl + 1
         a_rpp = FloatArray3D(a_sclrp(:,:,:,iscl),store=.FALSE.)
         a_rpt = FloatArray3D(a_sclrt(:,:,:,iscl),store=.FALSE.)
         pipeline_p => a_rpp
         pipeline_t => a_rpt
         CALL Prog%NewField( "rpp","Precipitation mixing ratio",               &
                             "kg/kg", "tttt", ANY(outputlist == "rpp"),        &
                             pipeline_p, in_t_data = pipeline_t,               &
                             in_group = "LES"                                  &
                           )

         iscl = iscl + 1
         a_npp = FloatArray3D(a_sclrp(:,:,:,iscl),store=.FALSE.)
         a_npt = FloatArray3D(a_sclrt(:,:,:,iscl),store=.FALSE.)
         pipeline_p => a_npp
         pipeline_t => a_npt
         CALL Prog%NewField( "npp","Precipitation number",             &
                             "#/kg", "tttt", ANY(outputlist == "npp"), &
                             pipeline_p, in_t_data = pipeline_t,       &
                             in_group = "LES"                          &
                           )                              
      END IF
                         
      IF (isgstyp > 1) THEN
         iscl = iscl + 1
         a_qp = FloatArray3D(a_sclrp(:,:,:,iscl),store=.FALSE.)
         a_qt = FloatArray3D(a_sclrt(:,:,:,iscl),store=.FALSE.)
         pipeline_p => a_qp
         pipeline_t => a_qt 
         CALL Prog%NewField( "qp","CHECK",                            &
                             "CHECK", "tttt", ANY(outputlist == "qp"), &
                             pipeline_p, in_t_data = pipeline_t,      &
                             in_group = "LES"                         &
                           )    

      END IF

      
      ! SALSA Variables
      ! ------------------------------------------------------
      IF (level >= 4) THEN
         
         iscl = iscl + 1
         a_naerop = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+nbins-1),store=.FALSE.)
         a_naerot = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+nbins-1),store=.FALSE.)
         pipeline_p => a_naerop
         pipeline_t => a_naerot
         CALL Prog%NewField( "naero","Binned aerosol number concentration",                  &
                             "#/kg", "ttttaea", ( ANY(outputlist == "naero") .AND. lbinanl),  &
                             pipeline_p, in_t_data = pipeline_t,                             &
                             in_group = "SALSA_4d"                                           &
                           )    
         iscl = iscl + nbins-1

         iscl = iscl + 1
         a_maerop = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+nbins*nspec-1),store=.FALSE.)
         a_maerot = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+nbins*nspec-1),store=.FALSE.)
         pipeline_p => a_maerop
         pipeline_t => a_maerot
         CALL Prog%NewField( "maero","Binned aerosol mass concentration, sequentially for each compound",    &
                             "kg/kg", "N/A", .FALSE.,                                                        &
                             pipeline_p, in_t_data = pipeline_t,                                             &
                             in_group = "SALSA_4d"                                                           &
                           )    
         iscl = iscl + nbins*nspec-1
         
         iscl = iscl + 1
         a_ncloudp = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+ncld-1),store=.FALSE.)
         a_ncloudt = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+ncld-1),store=.FALSE.)
         pipeline_p => a_ncloudp
         pipeline_t => a_ncloudt
         CALL Prog%NewField( "ncloud","Binned cloud number concentration",                    &
                             "#/kg", "ttttcld", ( ANY(outputlist == "ncloud") .AND. lbinanl),  &
                             pipeline_p, in_t_data = pipeline_t,                              &
                             in_group = "SALSA_4d"                                            &
                           )    
         iscl = iscl + ncld-1

         iscl = iscl + 1
         a_mcloudp = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+ncld*nspec-1),store=.FALSE.)
         a_mcloudt = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+ncld*nspec-1),store=.FALSE.)
         pipeline_p => a_mcloudp
         pipeline_t => a_mcloudt
         CALL Prog%NewField( "mcloud","Binned cloud mass concentration, sequentially for each compound",    &
                             "kg/kg", "N/A", .FALSE.,                                                       &
                             pipeline_p, in_t_data = pipeline_t,                                            &
                             in_group = "SALSA_4d"                                                          &
                           )    
         iscl = iscl + ncld*nspec-1
         
         iscl = iscl + 1
         a_nprecpp = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+nprc-1),store=.FALSE.)
         a_nprecpt = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+nprc-1),store=.FALSE.)
         pipeline_p => a_nprecpp
         pipeline_t => a_nprecpt
         CALL Prog%NewField( "nprecp","Binned precip number concentration",                  &
                             "#/kg", "ttttprc", ( ANY(outputlist == "nprecp") .AND. lbinanl), &
                             pipeline_p, in_t_data = pipeline_t,                             &
                             in_group = "SALSA_4d"                                           &
                           )                            
         iscl = iscl + nprc-1

         iscl = iscl + 1
         a_mprecpp = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+nprc*nspec-1),store=.FALSE.)
         a_mprecpt = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+nprc*nspec-1),store=.FALSE.)
         pipeline_p => a_mprecpp
         pipeline_t => a_mprecpt
         CALL Prog%NewField( "mprecp","Binned precip mass concentration, sequentially for each compound",    &
                             "kg/kg", "N/A", .FALSE.,                                                        &
                             pipeline_p, in_t_data = pipeline_t,                                             &
                             in_group = "SALSA_4d"                                                           &
                           )                            
         iscl = iscl + nprc*nspec-1
         
         iscl = iscl + 1
         a_gaerop = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+5-1),store=.FALSE.)
         a_gaerot = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+5-1),store=.FALSE.)
         pipeline_p => a_gaerop
         pipeline_t => a_gaerop
         CALL Prog%NewField( "gaero","Binned ice number concentration",    &
                             "#/kg", "ttttice", .FALSE.,                   &
                             pipeline_p, in_t_data = pipeline_t,           &
                             in_group = "SALSA_4d"                         &
                           )                            
         iscl = iscl + 5-1
         
      END IF

      IF (level == 5) THEN
         iscl = iscl + 1
         a_nicep = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+nice-1),store=.FALSE.)
         a_nicet = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+nice-1),store=.FALSE.)
         pipeline_p => a_nicep
         pipeline_t => a_nicet
         CALL Prog%NewField( "nice","Binned ice number concentration",                     &
                             "#/kg", "ttttice", ( ANY(outputlist == "nice") .AND. lbinanl), &
                             pipeline_p, in_t_data = pipeline_t,                           &
                             in_group = "SALSA_4d"                                         &
                           )                            
         iscl = iscl + nice-1

         iscl = iscl + 1
         a_micep = FloatArray4D(a_sclrp(:,:,:,iscl:iscl+nice*(nspec+1)-1),store=.FALSE.) ! nspec+1 for rimed ice!
         a_micet = FloatArray4D(a_sclrt(:,:,:,iscl:iscl+nice*(nspec+1)-1),store=.FALSE.)
         pipeline_p => a_micep
         pipeline_t => a_micet
         CALL Prog%NewField( "mice","Binned ice mass concentration, sequentially for each compound",    &
                             "kg/kg", "N/A", .FALSE.,                                                   &
                             pipeline_p, in_t_data = pipeline_t,                                        &
                             in_group = "SALSA_4d"                                                      &
                           )                            
         iscl = iscl + (nice+1)-1
         
      END IF
      

      pipeline_p => NULL()
      pipeline_t => NULL()
      
      
      
    END SUBROUTINE setPrognosticVariables

    


    
END MODULE mo_progn_state
