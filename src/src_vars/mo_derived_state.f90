MODULE mo_derived_state
  USE classFieldArray
  USE mo_structured_datatypes
  USE mo_derived_procedures
  IMPLICIT NONE

  SAVE

  !
  ! Derived diagnostic quantities
  ! Note that no data is stored for these variables, but the onDemand function is associated to get
  ! the necessary diagnostics when writing output.
  !

  ! General LES variables
  TYPE(FloatArray3d), TARGET :: qtot                     ! Total water content
  TYPE(FloatArray2d), TARGET :: lwp, iwp, rwp            ! Liquid water path, ice water path, rain water path.
                                                         ! LWP contains cloud droplets and aerosol, IWP both pristine and rimed ice,
                                                         ! and rwp the water in the "precipitation" category
  TYPE(FloatArray2d), TARGET :: shf, lhf                 ! Surface sensible and latent heat fluxes
  
  ! SALSA related variables
  TYPE(FloatArray3d), TARGET :: Naa, Nab, Nca, Ncb, Np, Ni,                &  ! Bulk number concentrations, aerosol, cloud, precip, ice
                                Dwaa, Dwab, Dwca, Dwcb, Dwpa, Dwia,        &  ! Bulk mean diameters
                                aSO4a, aSO4b, cSO4a, cSO4b, pSO4a, iSO4a,  &  ! Sulfate mass
                                aOCa, aOCb, cOCa, cOCb, pOCa, iOCa,        &  ! Organic carbon
                                aBCa, aBCb, cBCa, cBCb, pBCa, iBCa,        &  ! Black carbon
                                aDUa, aDUb, cDUa, cDUb, pDUa, iDUa,        &  ! Dust
                                aSSa, aSSb, cSSa, cSSb, pSSa, iSSa,        &  ! Sea salt
                                aNOa, aNOb, cNOa, cNOb, pNOa, iNOa,        &  ! Nitrate
                                aNHa, aNHb, cNHa, cNHb, pNHa, iNHa            ! Ammonia

  TYPE(FloatArray3d), TARGET :: gSO4, gNO3, gNH4, gOCNV, gOCSV

  ! A bit more derived SALSA variables
  TYPE(FloatArray3d), TARGET :: CDNC,                                      &  ! CDNC - number of all droplets for whom 2 um < D < 80 um. Note units in #/m3
                                CNC,                                       &  ! Cloud number concentration - number of all droplet for whom D > 2 um. In #/m3
                                Reff                                          ! Droplet effective radius for all D > 2 um
                                
  ! SALSA related variables
  TYPE(FloatArray4d), TARGET :: Dwaba, Dwabb, Dwcba, Dwcbb, Dwpba, Dwiba     ! Binned wet diameters  
  TYPE(FloatArray4d), TARGET :: Maba, Mabb, Mcba, Mcbb, Mpba, Miba           ! Binned particle/droplet total masses
  
  
  ! Some binned diagnostics
  TYPE(FloatArray4d), TARGET :: irhob, irhoe ! Bulk mean and effective ice densities

  
  CONTAINS

    SUBROUTINE setDerivedVariables(Derived,outputlist,level)
      TYPE(FieldArray), INTENT(inout) :: Derived
      CHARACTER(len=*), INTENT(in) :: outputlist(:)
      INTEGER, INTENT(in) :: level
      CLASS(*), POINTER :: pipeline => NULL()

      ! Since these are only for output, mask their allocation with the output status.
      ! Most of these variables are also used for profile mean statistical outputs,
      ! so allocate them also if they appear in outputlist_ps! However the outputstatus of
      ! the variables in this module is still determined solely by the main outputlist.

      pipeline => NULL()
      qtot = FloatArray3d()
      qtot%onDemand => totalWater
      pipeline => qtot
      CALL Derived%newField("qtot", "Total water content", "kg/kg", "tttt",   &
                            ANY(outputlist == "qtot"), pipeline               )

      pipeline => NULL()
      lwp = FloatArray2d()
      lwp%onDemand => waterPaths
      pipeline => lwp
      CALL Derived%newField("lwp", "Liquid water path", "kg/m2", "xtytt",   &
                            ANY(outputlist == "lwp"), pipeline              )

      pipeline => NULL()
      lhf = FloatArray2d()
      lhf%onDemand => surfaceFluxes
      pipeline => lhf
      CALL Derived%newField("lhf", "Latent heat flux", "W m-2", "xtytt",   &
                            ANY(outputlist == "lhf"), pipeline             )

      pipeline => NULL()
      shf = FloatArray2d()
      shf%onDemand => surfaceFluxes
      pipeline => shf
      CALL Derived%newField("shf", "Sensible heat flux", "W m-2", "xtytt",   &
                            ANY(outputlist == "shf"), pipeline                )
      
      
      IF (level > 4) THEN
         pipeline => NULL()
         iwp = FloatArray2d()
         iwp%onDemand => waterPaths
         pipeline => iwp
         CALL Derived%newField("iwp", "Ice water path", "kg/m2", "xtytt",    &
                               ANY(outputlist == "iwp"), pipeline           )
      END IF

      IF (level >= 3) THEN
         pipeline => NULL()
         rwp = FloatArray2d()
         rwp%onDemand => waterPaths
         pipeline => rwp
         CALL Derived%newField("rwp", "Rain water path", "kg/m2", "xtytt",   &
                               ANY(outputlist == "rwp"), pipeline            )
      END IF
            
      IF (level >= 4) THEN
         pipeline => NULL()
         Naa = FloatArray3d()
         Naa%onDemand => bulkNumc
         pipeline => Naa
         CALL Derived%newField("Naa", "Bulk number of aerosol, A", "kg-1", 'tttt',   &
                               ANY(outputlist == "Naa"), pipeline                   )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Nab = FloatArray3d()
         Nab%onDemand => bulkNumc
         pipeline => Nab
         CALL Derived%newField("Nab", "Bulk number of aerosol, B", "kg-1", 'tttt',   &
                               ANY(outputlist == "Nab"), pipeline                   )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Nca = FloatArray3d()
         Nca%onDemand => bulkNumc
         pipeline => Nca
         CALL Derived%newField("Nca", "Bulk number of cloud droplets, A", "kg-1", 'tttt',   &
                               ANY(outputlist == "Nca"), pipeline                          )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Ncb = FloatArray3d()
         Ncb%onDemand => bulkNumc
         pipeline => Ncb
         CALL Derived%newField("Ncb", "Bulk number of cloud droplets, B", "kg-1", 'tttt',   &
                               ANY(outputlist == "Ncb"), pipeline                          )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Np = FloatArray3d()
         Np%onDemand => bulkNumc
         pipeline => Np
         CALL Derived%newField("Np", "Bulk number of precip", "kg-1", 'tttt',   &
                               ANY(outputlist == "Np"), pipeline               )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Ni = FloatArray3d()
         Ni%onDemand => bulkNumc
         pipeline => Ni
         CALL Derived%newField("Ni", "Bulk number of ice", "kg-1", 'tttt',   &
                               ANY(outputlist == "Ni"), pipeline            )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         CDNC = FloatArray3d()
         CDNC%onDemand => getCDNC
         pipeline => CDNC
         CALL Derived%newField("CDNC", "Cloud droplet number concentration",      &
                               "m-3", 'tttt', ANY(outputlist == "CDNC"), pipeline)
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         CNC = FloatArray3d()
         CNC%onDemand => getCNC
         pipeline => CNC
         CALL Derived%newField("CNC", "Cloud number concentration", "m-3", 'tttt',  &
                               ANY(outputlist == "CNC"), pipeline                    )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Reff = FloatArray3d()
         Reff%onDemand => getReff
         pipeline => Reff
         CALL Derived%newField("Reff", "Droplet effective radius", "m", 'tttt',   &
                               ANY(outputlist == "Reff"), pipeline                )
      END IF
      
      IF (level >= 4) THEN
         pipeline => NULL()
         Dwaa = FloatArray3d()
         Dwaa%onDemand => bulkDiameter
         pipeline => Dwaa
         CALL Derived%newField("Dwaa", "Bulk mean diameter, aerosol A", "m", 'tttt',   &
                               ANY(outputlist == "Dwaa"), pipeline                     )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Dwab = FloatArray3d()
         Dwab%onDemand => bulkDiameter
         pipeline => Dwab
         CALL Derived%newField("Dwab", "Bulk mean diameter, aerosol B", "m", 'tttt',   &
                               ANY(outputlist == "Dwab"), pipeline                     )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Dwca = FloatArray3d()
         Dwca%onDemand => bulkDiameter
         pipeline => Dwca
         CALL Derived%newField("Dwca", "Bulk mean diameter, clouds A", "m", 'tttt',   &
                               ANY(outputlist == "Dwca"), pipeline                    )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Dwcb = FloatArray3d()
         Dwcb%onDemand => bulkDiameter
         pipeline => Dwcb
         CALL Derived%newField("Dwcb", "Bulk mean diameter, clouds B", "m", 'tttt',   &
                               ANY(outputlist == "Dwcb"), pipeline                    )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Dwpa = FloatArray3d()
         Dwpa%onDemand => bulkDiameter
         pipeline => Dwpa
         CALL Derived%newField("Dwpa", "Bulk mean diameter, precip", "m", 'tttt',   &
                               ANY(outputlist == "Dwpa"), pipeline                  )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Dwia = FloatArray3d()
         Dwia%onDemand => bulkDiameter
         pipeline => Dwia
         CALL Derived%newField("Dwia", "Bulk mean diameter, ice", "m", 'tttt',   &
                               ANY(outputlist == "Dwia"), pipeline               )
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         Dwaba = FloatArray4d()
         Dwaba%onDemand => getBinDiameter
         pipeline => Dwaba
         CALL Derived%newField("Dwaba", "Bin diameter, aerosol A", "m", 'ttttaea',   &
                               ANY(outputlist == "Dwaba"), pipeline                  )

         pipeline => NULL()
         Maba = FloatArray4d()
         Maba%onDemand => getBinTotMass
         pipeline => Maba
         CALL Derived%newField("Maba", "Bin total mass, aerosol A", "kg/kg", "ttttaea",   &
                               ANY(outputlist == "Maba"), pipeline                        )
         
         pipeline => NULL()
         Dwabb = FloatArray4d()
         Dwabb%onDemand => getBinDiameter
         pipeline => Dwabb
         CALL Derived%newField("Dwabb", "Bin diameter, aerosol B", "m", 'ttttaeb',   &
                               ANY(outputlist == "Dwabb"), pipeline                  )

         pipeline => NULL()
         Mabb = FloatArray4d()
         Mabb%onDEmand => getBinTotMass
         pipeline => Mabb
         CALL Derived%newField("Mabb", "Bin total mass, aerosol B", "kg/kg", "ttttaeb",   &
                               ANY(outputlist == "Mabb"), pipeline                        )
                  
         pipeline => NULL()
         Dwcba = FloatArray4d()
         Dwcba%onDemand => getBinDiameter
         pipeline => Dwcba
         CALL Derived%newField("Dwcba", "Bin diameter, clouds A", "m", 'ttttcla',   &
                               ANY(outputlist == "Dwcba"), pipeline                 )

         pipeline => NULL()
         Mcba = FloatArray4d()
         Mcba%onDemand => getBinTotMass
         pipeline => Mcba
         CALL Derived%newField("Mcba", "Bin total mass, clouds A", "kg/kg", "ttttcla",  &
                               ANY(outputlist == "Mcba"), pipeline                      )
         
         pipeline => NULL()
         Dwcbb = FloatArray4d()
         Dwcbb%onDemand => getBinDiameter
         pipeline => Dwcbb
         CALL Derived%newField("Dwcbb", "Bin diameter, clouds B", "m", 'ttttclb',   &
                               ANY(outputlist == "Dwcbb"), pipeline                 )

         pipeline => NULL()
         Mcbb = FloatArray4d()
         Mcbb%onDemand => getBinTotMass
         pipeline => Mcbb
         CALL Derived%newField("Mcbb", "Bin total mass, clouds B", "kg/kg", "ttttclb",   &
                               ANY(outputlist == "Mcbb"), pipeline                       )
         
         pipeline => NULL()
         Dwpba = FloatArray4d()
         Dwpba%onDemand => getBinDiameter
         pipeline => Dwpba
         CALL Derived%newField("Dwpba", "Bin diameter, precip", "m", 'ttttprc',   &
                               ANY(outputlist == "Dwpba"), pipeline               )

         pipeline => NULL()
         Mpba = FloatArray4d()
         Mpba%onDemand => getBinTotMass
         pipeline => Mpba
         CALL Derived%newField("Mpba", "Bin total mass, precip", "kg/kg", "ttttprc",   &
                               ANY(outputlist == "Mpba"), pipeline                     )

      END IF

      IF (level == 5) THEN
         pipeline => NULL()
         Dwiba = FloatArray4d()
         Dwiba%onDemand => getBinDiameter
         pipeline => Dwiba
         CALL Derived%newField("Dwiba", "Bin diameter, ice", "m", 'ttttice',   &
                               ANY(outputlist == "Dwiba"), pipeline            )

         pipeline => NULL()
         Miba = FloatArray4d()
         Miba%onDemand => getBinTotMass
         pipeline => Miba
         CALL Derived%newField("Miba", "Bin total mass, ice", "m", "ttttice",       &
                               ANY(outputlist == "Miba"), pipeline                  )

      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         aSO4a = FloatArray3d()
         aSO4a%onDemand => bulkMixrat
         pipeline => aSO4a
         CALL Derived%newField("aSO4a", "Bulk SO4 in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aSO4a"), pipeline                 )

         pipeline => NULL()
         aSO4b = FloatArray3d()
         aSO4b%onDemand => bulkMixrat
         pipeline => aSO4b
         CALL Derived%newField("aSO4b", "Bulk SO4 in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aSO4b"), pipeline                 )

         pipeline => NULL()
         cSO4a = FloatArray3d()
         cSO4a%onDemand => bulkMixrat
         pipeline => cSO4a
         CALL Derived%newField("cSO4a", "Bulk SO4 in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cSO4a"), pipeline                )

         pipeline => NULL()
         cSO4b = FloatArray3d()
         cSO4b%onDemand => bulkMixrat
         pipeline => cSO4b
         CALL Derived%newField("cSO4b", "Bulk SO4 in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cSO4b"), pipeline                )

         pipeline => NULL()
         pSO4a = FloatArray3d()
         pSO4a%onDemand => bulkMixrat
         pipeline => pSO4a
         CALL Derived%newField("pSO4a", "Bulk SO4 in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pSO4a"), pipeline              )

         pipeline => NULL()
         aOCa = FloatArray3d()
         aOCa%onDemand => bulkMixrat
         pipeline => aOCa
         CALL Derived%newField("aOCa", "Bulk OC in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aOCa"), pipeline                )

         pipeline => NULL()
         aOCb = FloatArray3d()
         aOCb%onDemand => bulkMixrat
         pipeline => aOCb
         CALL Derived%newField("aOCb", "Bulk OC in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aOCb"), pipeline                )

         pipeline => NULL()
         cOCa = FloatArray3d()
         cOCa%onDemand => bulkMixrat
         pipeline => cOCa
         CALL Derived%newField("cOCa", "Bulk OC in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cOCa"), pipeline               )

         pipeline => NULL()
         cOCb = FloatArray3d()
         cOCb%onDemand => bulkMixrat
         pipeline => cOCb
         CALL Derived%newField("cOCb", "Bulk OC in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cOCb"), pipeline               )

         pipeline => NULL()
         pOCa = FloatArray3d()
         pOCa%onDemand => bulkMixrat
         pipeline => pOCa
         CALL Derived%newField("pOCa", "Bulk OC in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pOCa"), pipeline             )

         pipeline => NULL()
         aBCa = FloatArray3d()
         aBCa%onDemand => bulkMixrat
         pipeline => aBCa
         CALL Derived%newField("aBCa", "Bulk BC in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aBCa"), pipeline                )

         pipeline => NULL()
         aBCb = FloatArray3d()
         aBCb%onDemand => bulkMixrat
         pipeline => aBCb
         CALL Derived%newField("aBCb", "Bulk BC in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aBCb"), pipeline                )

         pipeline => NULL()
         cBCa = FloatArray3d()
         cBCa%onDemand => bulkMixrat
         pipeline => cBCa
         CALL Derived%newField("cBCa", "Bulk BC in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cBCa"), pipeline               )

         pipeline => NULL()
         cBCb = FloatArray3d()
         cBCb%onDemand => bulkMixrat
         pipeline => cBCb
         CALL Derived%newField("cBCb", "Bulk BC in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cBCb"), pipeline               )

         pipeline => NULL()
         pBCa = FloatArray3d()
         pBCa%onDemand => bulkMixrat
         pipeline => pBCa
         CALL Derived%newField("pBCa", "Bulk BC in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pBCa"), pipeline             ) 

         pipeline => NULL()
         aDUa = FloatArray3d()
         aDUa%onDemand => bulkMixrat
         pipeline => aDUa
         CALL Derived%newField("aDUa", "Bulk DU in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aDUa"), pipeline                )

         pipeline => NULL()
         aDUb = FloatArray3d()
         aDUb%onDemand => bulkMixrat
         pipeline => aDUb
         CALL Derived%newField("aDUb", "Bulk DU in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aDUb"), pipeline                )

         pipeline => NULL()
         cDUa = FloatArray3d()
         cDUa%onDemand => bulkMixrat
         pipeline => cDUa
         CALL Derived%newField("cDUa", "Bulk DU in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cDUa"), pipeline               )

         pipeline => NULL()
         cDUb = FloatArray3d()
         cDUb%onDemand => bulkMixrat
         pipeline => cDUb
         CALL Derived%newField("cDUb", "Bulk DU in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cDUb"), pipeline               )

         pipeline => NULL()
         pDUa = FloatArray3d()
         pDUa%onDemand => bulkMixrat
         pipeline => pDUa
         CALL Derived%newField("pDUa", "Bulk DU in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pDUa"), pipeline             )

         pipeline => NULL()
         aSSa = FloatArray3d()
         aSSa%onDemand => bulkMixrat
         pipeline => aSSa
         CALL Derived%newField("aSSa", "Bulk SS in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aSSa"), pipeline                )

         pipeline => NULL()
         aSSb = FloatArray3d()
         aSSb%onDemand => bulkMixrat
         pipeline => aSSb
         CALL Derived%newField("aSSb", "Bulk SS in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aSSb"), pipeline                )

         pipeline => NULL()
         cSSa = FloatArray3d()
         cSSa%onDemand => bulkMixrat
         pipeline => cSSa
         CALL Derived%newField("cSSa", "Bulk SS in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cSSa"), pipeline               ) 

         pipeline => NULL()
         cSSb = FloatArray3d()
         cSSb%onDemand => bulkMixrat
         pipeline => cSSb
         CALL Derived%newField("cSSb", "Bulk SS in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cSSb"), pipeline               )

         pipeline => NULL()
         pSSa = FloatArray3d()
         pSSa%onDemand => bulkMixrat
         pipeline => pSSa
         CALL Derived%newField("pSSa", "Bulk SS in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pSSa"), pipeline             )

         pipeline => NULL()
         aNOa = FloatArray3d()
         aNOa%onDemand => bulkMixrat
         pipeline => aNOa
         CALL Derived%newField("aNOa", "Bulk NO in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aNOa"), pipeline                )

         pipeline => NULL()
         aNOb = FloatArray3d()
         aNOb%onDemand => bulkMixrat
         pipeline => aNOb
         CALL Derived%newField("aNOb", "Bulk NO in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aNOb"), pipeline                )

         pipeline => NULL()
         cNOa = FloatArray3d()
         cNOa%onDemand => bulkMixrat
         pipeline => cNOa
         CALL Derived%newField("cNOa", "Bulk NO in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cNOa"), pipeline               )

         pipeline => NULL()
         cNOb = FloatArray3d()
         cNOb%onDemand => bulkMixrat
         pipeline => cNOb
         CALL Derived%newField("cNOb", "Bulk NO in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cNOb"), pipeline               )

         pipeline => NULL()
         pNOa = FloatArray3d()
         pNOa%onDemand => bulkMixrat
         pipeline => pNOa
         CALL Derived%newField("pNOa", "Bulk NO in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pNOa"), pipeline             )

         pipeline => NULL()
         aNHa = FloatArray3d()
         aNHa%onDemand => bulkMixrat
         pipeline => aNHa
         CALL Derived%newField("aNHa", "Bulk NH in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aNHa"), pipeline                )

         pipeline => NULL()
         aNHb = FloatArray3d()
         aNHb%onDemand => bulkMixrat
         pipeline => aNHb
         CALL Derived%newField("aNHb", "Bulk NH in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aNHb"), pipeline                )

         pipeline => NULL()
         cNHa = FloatArray3d()
         cNHa%onDemand => bulkMixrat
         pipeline => cNHa
         CALL Derived%newField("cNHa", "Bulk NH in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cNHa"), pipeline               )

         pipeline => NULL()
         cNHb = FloatArray3d()
         cNHb%onDemand => bulkMixrat
         pipeline => cNHb
         CALL Derived%newField("cNHb", "Bulk NH in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cNHb"), pipeline               )

         pipeline => NULL()
         pNHa = FloatArray3d()
         pNHa%onDemand => bulkMixrat
         pipeline => pNHa
         CALL Derived%newField("pNHa", "Bulk NH in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pNHa"), pipeline             )
      END IF

      IF ( level == 5 ) THEN

         pipeline => NULL()
         iSO4a = FloatArray3d()
         iSO4a%onDemand => bulkMixrat
         pipeline => iSO4a
         CALL Derived%newField("iSO4a", "Bulk SO4 in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iSO4a"), pipeline           ) 

         pipeline => NULL()
         iOCa = FloatArray3d()
         iOCa%onDemand => bulkMixrat
         pipeline => iOCa
         CALL Derived%newField("iOCa", "Bulk OC in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iOCa"), pipeline          )

         pipeline => NULL()
         iBCa = FloatArray3d()
         iBCa%onDemand => bulkMixrat
         pipeline => iBCa
         CALL Derived%newField("iBCa", "Bulk BC in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iBCa"), pipeline          )

         pipeline => NULL()
         iDUa = FloatArray3d()
         iDUa%onDemand => bulkMixrat
         pipeline => iDUa
         CALL Derived%newField("iDUa", "Bulk DU in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iDUa"), pipeline          )      
         
         pipeline => NULL()
         iSSa = FloatArray3d()
         iSSa%onDemand => bulkMixrat
         pipeline => iSSa
         CALL Derived%newField("iSSa", "Bulk SS in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iSSa"), pipeline          )

         pipeline => NULL()
         iNOa = FloatArray3d()
         iNOa%onDemand => bulkMixrat
         pipeline => iNOa
         CALL Derived%newField("iNOa", "Bulk NO in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iNOa"), pipeline          )      

         pipeline => NULL()
         iNHa = FloatArray3d()
         iNHa%onDemand => bulkMixrat
         pipeline => iNHa
         CALL Derived%newField("iNHa", "Bulk NH in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iNHa"), pipeline          )

         pipeline => NULL()
         irhob = FloatArray4d()
         irhob%onDemand => binIceDensities 
         pipeline => irhob
         CALL Derived%newField("irhob", "Bulk mean density of ice", "kg/m3", "ttttice",   &
                            ANY(outputlist == "irhob"), pipeline)

         pipeline => NULL()
         irhoe = FloatArray4d()
         irhoe%onDemand => binIceDensities
         pipeline => irhoe
         CALL Derived%newField("irhoe", "Effective mean density of ice", "kg/m3", "ttttice",   &
                            ANY(outputlist == "irhoe"), pipeline)
      END IF

      IF (level >= 4) THEN
         pipeline => NULL()
         gSO4 = FloatArray3d()
         gSO4%onDemand => getGasConc
         pipeline => gSO4
         CALL Derived%newField("gSO4", "Gas concentration SO4", "#/kg", "tttt",   &
                               ANY(outputlist == "gSO4"), pipeline)

         pipeline => NULL()
         gNO3 = FloatArray3d()
         gNO3%onDemand => getGasConc
         pipeline => gNO3
         CALL Derived%newField("gNO3", "Gas concentration NO3", "#/kg", "tttt",   &
              ANY(outputlist == "gNO3"), pipeline)

         pipeline => NULL()
         gNH4 = FloatArray3d()
         gNH4%onDemand => getGasConc
         pipeline => gSO4
         CALL Derived%newField("gSO4", "Gas concentration SO4", "#/kg", "tttt",   &
              ANY(outputlist == "gSO4"), pipeline)

         pipeline => NULL()
         gOCNV = FloatArray3d()
         gOCNV%onDemand => getGasConc
         pipeline => gOCNV
         CALL Derived%newField("gOCNV", "Gas concentration OCNV", "#/kg", "tttt",   &
                               ANY(outputlist == "gOCNV"), pipeline)                  

         pipeline => NULL()
         gOCSV = FloatArray3d()
         gOCSV%onDemand => getGasConc
         pipeline => gOCSV
         CALL Derived%newField("gOCSV", "Gas concentration OCSV", "#/kg", "tttt",   &
                               ANY(outputlist == "gOCSV"), pipeline)                  

         
      END IF

      
      pipeline => NULL()
      
    END SUBROUTINE setDerivedVariables        

END MODULE mo_derived_state
