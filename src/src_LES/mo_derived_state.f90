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
  TYPE(FLoatArray3d), TARGET :: qtot     ! Total water content
  
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

  
  
  CONTAINS

    SUBROUTINE setDerivedVariables(Derived,outputlist,outputlist_ps,lsalsabbins,level)
      TYPE(FieldArray), INTENT(inout) :: Derived
      CHARACTER(len=*), INTENT(in) :: outputlist(:), outputlist_ps(:)
      LOGICAL, INTENT(in) :: lsalsabbins
      INTEGER, INTENT(in) :: level
      CLASS(*), POINTER :: pipeline => NULL()

      ! Since these are only for output, mask their allocation with the output status.
      ! Most of these variables are also used for profile mean statistical outputs,
      ! so allocate them also if they appear in outputlist_ps! However the outputstatus of
      ! the variables in this module is still determined solely by the main outputlist.

      IF (ANY(outputlist == "qtot" .OR. outputlist_ps == "qtot")) THEN
         qtot = FloatArray3d()
         qtot%onDemand => totalWater
         pipeline => qtot
         CALL Derived%newField("qtot", "Total water content", "kg/kg", "tttt",   &
                               ANY(outputlist == "qtot"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Naa" .OR. outputlist_ps == "Naa")) THEN
         Naa = FloatArray3d()
         Naa%onDemand => bulkNumc
         pipeline => Naa
         CALL Derived%newField("Naa", "Bulk number of aerosol, A", "m-3", 'tttt',   &
                               ANY(outputlist == "Naa"), pipeline                   )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Nab" .OR. outputlist_ps == "Nab")  &
          .AND. lsalsabbins) THEN
         Nab = FloatArray3d()
         Nab%onDemand => bulkNumc
         pipeline => Nab
         CALL Derived%newField("Nab", "Bulk number of aerosol, B", "m-3", 'tttt',   &
                               ANY(outputlist == "Nab"), pipeline                   )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Nca" .OR. outputlist_ps == "Nca")) THEN
         Nca = FloatArray3d()
         Nca%onDemand => bulkNumc
         pipeline => Nca
         CALL Derived%newField("Nca", "Bulk number of cloud droplets, A", "m-3", 'tttt',   &
                               ANY(outputlist == "Nca"), pipeline                          )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Ncb" .OR. outputlist_ps == "Ncb") &
          .AND. lsalsabbins) THEN         
         Ncb = FloatArray3d()
         Ncb%onDemand => bulkNumc
         pipeline => Ncb
         CALL Derived%newField("Ncb", "Bulk number of cloud droplets, B", "m-3", 'tttt',   &
                               ANY(outputlist == "Ncb"), pipeline                          )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Np" .OR. outputlist_ps == "Np")) THEN      
         Np = FloatArray3d()
         Np%onDemand => bulkNumc
         pipeline => Np
         CALL Derived%newField("Np", "Bulk number of precip", "m-3", 'tttt',   &
                               ANY(outputlist == "Np"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Ni" .OR. outputlist_ps == "Ni")) THEN
         Ni = FloatArray3d()
         Ni%onDemand => bulkNumc
         pipeline => Ni
         CALL Derived%newField("Ni", "Bulk number of ice", "m-3", 'tttt',   &
                               ANY(outputlist == "Ni"), pipeline            )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwaa" .OR. outputlist_ps == "Dwaa")) THEN
         Dwaa = FloatArray3d()
         Dwaa%onDemand => bulkDiameter
         pipeline => Dwaa
         CALL Derived%newField("Dwaa", "Bulk mean diameter, aerosol A", "m", 'tttt',   &
                               ANY(outputlist == "Dwaa"), pipeline                     )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwab" .OR. outputlist_ps == "Dwab") &
          .AND. lsalsabbins) THEN
         Dwab = FloatArray3d()
         Dwab%onDemand => bulkDiameter
         pipeline => Dwab
         CALL Derived%newField("Dwab", "Bulk mean diameter, aerosol B", "m", 'tttt',   &
                               ANY(outputlist == "Dwab"), pipeline                     )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwca" .OR. outputlist_ps == "Dwca")) THEN
         Dwca = FloatArray3d()
         Dwca%onDemand => bulkDiameter
         pipeline => Dwca
         CALL Derived%newField("Dwca", "Bulk mean diameter, clouds A", "m", 'tttt',   &
                               ANY(outputlist == "Dwca"), pipeline                    )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwcb" .OR. outputlist_ps == "Dwcb")) THEN
         Dwcb = FloatArray3d()
         Dwcb%onDemand => bulkDiameter
         pipeline => Dwcb
         CALL Derived%newField("Dwcb", "Bulk mean diameter, clouds B", "m", 'tttt',   &
                               ANY(outputlist == "Dwcb"), pipeline                    )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwpa" .OR. outputlist_ps == "Dwpa")) THEN
         Dwpa = FloatArray3d()
         Dwpa%onDemand => bulkDiameter
         pipeline => Dwpa
         CALL Derived%newField("Dwpa", "Bulk mean diameter, precip", "m", 'tttt',   &
                               ANY(outputlist == "Dwpa"), pipeline                  )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwia" .OR. outputlist_ps == "Dwia")) THEN
         Dwia = FloatArray3d()
         Dwia%onDemand => bulkDiameter
         pipeline => Dwia
         CALL Derived%newField("Dwia", "Bulk mean diameter, ice", "m", 'tttt',   &
                               ANY(outputlist == "Dwia"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aSO4a" .OR. outputlist_ps == "aSO4a")) THEN
         aSO4a = FloatArray3d()
         aSO4a%onDemand => bulkMixrat
         pipeline => aSO4a
         CALL Derived%newField("aSO4a", "Bulk SO4 in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aSO4a"), pipeline                 )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aSO4b" .OR. outputlist_ps == "aSO4b") &
          .AND. lsalsabbins) THEN
         aSO4b = FloatArray3d()
         aSO4b%onDemand => bulkMixrat
         pipeline => aSO4b
         CALL Derived%newField("aSO4b", "Bulk SO4 in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aSO4b"), pipeline                 )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cSO4a" .OR. outputlist_ps == "cSO4a")) THEN
         cSO4a = FloatArray3d()
         cSO4a%onDemand => bulkMixrat
         pipeline => cSO4a
         CALL Derived%newField("cSO4a", "Bulk SO4 in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cSO4a"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cSO4b" .OR. outputlist_ps == "cSO4b")  &
          .AND. lsalsabbins) THEN
         cSO4b = FloatArray3d()
         cSO4b%onDemand => bulkMixrat
         pipeline => cSO4b
         CALL Derived%newField("cSO4b", "Bulk SO4 in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cSO4b"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pSO4a" .OR. outputlist_ps == "pSO4a")) THEN
         pSO4a = FloatArray3d()
         pSO4a%onDemand => bulkMixrat
         pipeline => pSO4a
         CALL Derived%newField("pSO4a", "Bulk SO4 in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pSO4a"), pipeline              )
      END IF
      
      IF (level >= 4 .AND. ANY(outputlist == "iSO4a" .OR. outputlist_ps == "iSO4a")) THEN      
         iSO4a = FloatArray3d()
         iSO4a%onDemand => bulkMixrat
         pipeline => iSO4a
         CALL Derived%newField("iSO4a", "Bulk SO4 in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iSO4a"), pipeline           ) 
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aOCa" .OR. outputlist_ps == "aOCa")) THEN
         aOCa = FloatArray3d()
         aOCa%onDemand => bulkMixrat
         pipeline => aOCa
         CALL Derived%newField("aOCa", "Bulk OC in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aOCa"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aOCb" .OR. outputlist_ps == "aOCb")  &
          .AND. lsalsabbins) THEN
         aOCb = FloatArray3d()
         aOCb%onDemand => bulkMixrat
         pipeline => aOCb
         CALL Derived%newField("aOCb", "Bulk OC in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aOCb"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cOCa" .OR. outputlist_ps == "cOCa")) THEN
         cOCa = FloatArray3d()
         cOCa%onDemand => bulkMixrat
         pipeline => cOCa
         CALL Derived%newField("cOCa", "Bulk OC in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cOCa"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cOCb" .OR. outputlist_ps == "cOCb")   &
          .AND. lsalsabbins) THEN
         cOCb = FloatArray3d()
         cOCb%onDemand => bulkMixrat
         pipeline => cOCb
         CALL Derived%newField("cOCb", "Bulk OC in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cOCb"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pOCa" .OR. outputlist_ps == "pOCa")) THEN
         pOCa = FloatArray3d()
         pOCa%onDemand => bulkMixrat
         pipeline => pOCa
         CALL Derived%newField("pOCa", "Bulk OC in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pOCa"), pipeline             )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iOCa" .OR. outputlist_ps == "iOCa")) THEN
         iOCa = FloatArray3d()
         iOCa%onDemand => bulkMixrat
         pipeline => iOCa
         CALL Derived%newField("iOCa", "Bulk OC in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iOCa"), pipeline          )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aBCa" .OR. outputlist_ps == "aBCa")) THEN
         aBCa = FloatArray3d()
         aBCa%onDemand => bulkMixrat
         pipeline => aBCa
         CALL Derived%newField("aBCa", "Bulk BC in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aBCa"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aBCb" .OR. outputlist_ps == "aBCb")  &
          .AND. lsalsabbins) THEN
         aBCb = FloatArray3d()
         aBCb%onDemand => bulkMixrat
         pipeline => aBCb
         CALL Derived%newField("aBCb", "Bulk BC in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aBCb"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cBCa" .OR. outputlist_ps == "cBCa")) THEN
         cBCa = FloatArray3d()
         cBCa%onDemand => bulkMixrat
         pipeline => cBCa
         CALL Derived%newField("cBCa", "Bulk BC in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cBCa"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cBCb" .OR. outputlist_ps == "cBCb") &
          .AND. lsalsabbins) THEN
         cBCb = FloatArray3d()
         cBCb%onDemand => bulkMixrat
         pipeline => cBCb
         CALL Derived%newField("cBCb", "Bulk BC in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cBCb"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pBCa" .OR. outputlist_ps == "pBCa")) THEN
         pBCa = FloatArray3d()
         pBCa%onDemand => bulkMixrat
         pipeline => pBCa
         CALL Derived%newField("pBCa", "Bulk BC in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pBCa"), pipeline             ) 
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iBCa" .OR. outputlist_ps == "iBCa")) THEN
         iBCa = FloatArray3d()
         iBCa%onDemand => bulkMixrat
         pipeline => iBCa
         CALL Derived%newField("iBCa", "Bulk BC in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iBCa"), pipeline          )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aDUa" .OR. outputlist_ps == "aDUa")) THEN
         aDUa = FloatArray3d()
         aDUa%onDemand => bulkMixrat
         pipeline => aDUa
         CALL Derived%newField("aDUa", "Bulk DU in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aDUa"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aDUb" .OR. outputlist_ps == "aDUb")   &
          .AND. lsalsabbins) THEN
         aDUb = FloatArray3d()
         aDUb%onDemand => bulkMixrat
         pipeline => aDUb
         CALL Derived%newField("aDUb", "Bulk DU in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aDUb"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cDUa" .OR. outputlist_ps == "cDUa")) THEN
         cDUa = FloatArray3d()
         cDUa%onDemand => bulkMixrat
         pipeline => cDUa
         CALL Derived%newField("cDUa", "Bulk DU in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cDUa"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cDUb" .OR. outputlist_ps == "cDUb")) THEN
         cDUb = FloatArray3d()
         cDUb%onDemand => bulkMixrat
         pipeline => cDUb
         CALL Derived%newField("cDUb", "Bulk DU in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cDUb"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pDUa" .OR. outputlist_ps == "pDUa")) THEN
         pDUa = FloatArray3d()
         pDUa%onDemand => bulkMixrat
         pipeline => pDUa
         CALL Derived%newField("pDUa", "Bulk DU in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pDUa"), pipeline             )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iDUa" .OR. outputlist_ps == "iDUa")) THEN
         iDUa = FloatArray3d()
         iDUa%onDemand => bulkMixrat
         pipeline => iDUa
         CALL Derived%newField("iDUa", "Bulk DU in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iDUa"), pipeline          )      
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aSSa" .OR. outputlist_ps == "aSSa")) THEN
         aSSa = FloatArray3d()
         aSSa%onDemand => bulkMixrat
         pipeline => aSSa
         CALL Derived%newField("aSSa", "Bulk SS in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aSSa"), pipeline                )
      END IF
      
      IF (level >= 4 .AND. ANY(outputlist == "aSSb" .OR. outputlist_ps == "aSSb")   &
          .AND. lsalsabbins) THEN  
         aSSb = FloatArray3d()
         aSSb%onDemand => bulkMixrat
         pipeline => aSSb
         CALL Derived%newField("aSSb", "Bulk SS in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aSSb"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cSSa" .OR. outputlist_ps == "cSSa")) THEN
         cSSa = FloatArray3d()
         cSSa%onDemand => bulkMixrat
         pipeline => cSSa
         CALL Derived%newField("cSSa", "Bulk SS in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cSSa"), pipeline               ) 
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cSSb" .OR. outputlist_ps == "cSSb")  &
          .AND. lsalsabbins) THEN
         cSSb = FloatArray3d()
         cSSb%onDemand => bulkMixrat
         pipeline => cSSb
         CALL Derived%newField("cSSb", "Bulk SS in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cSSb"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pSSa" .OR. outputlist_ps == "pSSa")) THEN
         pSSa = FloatArray3d()
         pSSa%onDemand => bulkMixrat
         pipeline => pSSa
         CALL Derived%newField("pSSa", "Bulk SS in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pSSa"), pipeline             )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iSSa" .OR. outputlist_ps == "iSSa")) THEN
         iSSa = FloatArray3d()
         iSSa%onDemand => bulkMixrat
         pipeline => iSSa
         CALL Derived%newField("iSSa", "Bulk SS in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iSSa"), pipeline          )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aNOa" .OR. outputlist_ps == "aNOa")) THEN
         aNOa = FloatArray3d()
         aNOa%onDemand => bulkMixrat
         pipeline => aNOa
         CALL Derived%newField("aNOa", "Bulk NO in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aNOa"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aNOb" .OR. outputlist_ps == "aNOb")   &
          .AND. lsalsabbins) THEN
         aNOb = FloatArray3d()
         aNOb%onDemand => bulkMixrat
         pipeline => aNOb
         CALL Derived%newField("aNOb", "Bulk NO in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aNOb"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cNOa" .OR. outputlist_ps == "cNOa")) THEN
         cNOa = FloatArray3d()
         cNOa%onDemand => bulkMixrat
         pipeline => cNOa
         CALL Derived%newField("cNOa", "Bulk NO in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cNOa"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cNOb" .OR. outputlist_ps == "cNOb")  &
          .AND. lsalsabbins) THEN
         cNOb = FloatArray3d()
         cNOb%onDemand => bulkMixrat
         pipeline => cNOb
         CALL Derived%newField("cNOb", "Bulk NO in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cNOb"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pNOa" .OR. outputlist_ps == "pNOa")) THEN
         pNOa = FloatArray3d()
         pNOa%onDemand => bulkMixrat
         pipeline => pNOa
         CALL Derived%newField("pNOa", "Bulk NO in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pNOa"), pipeline             )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iNOa" .OR. outputlist_ps == "iNOa")) THEN
         iNOa = FloatArray3d()
         iNOa%onDemand => bulkMixrat
         pipeline => iNOa
         CALL Derived%newField("iNOa", "Bulk NO in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iNOa"), pipeline          )      
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aNHa" .OR. outputlist_ps == "aNHa")) THEN
         aNHa = FloatArray3d()
         aNHa%onDemand => bulkMixrat
         pipeline => aNHa
         CALL Derived%newField("aNHa", "Bulk NH in aerosol A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aNHa"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aNHb" .OR. outputlist_ps == "aNHb")  &
          .AND. lsalsabbins) THEN
         aNHb = FloatArray3d()
         aNHb%onDemand => bulkMixrat
         pipeline => aNHb
         CALL Derived%newField("aNHb", "Bulk NH in aerosol B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "aNHb"), pipeline                )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cNHa" .OR. outputlist_ps == "cNHa")) THEN
         cNHa = FloatArray3d()
         cNHa%onDemand => bulkMixrat
         pipeline => cNHa
         CALL Derived%newField("cNHa", "Bulk NH in clouds A", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cNHa"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cNHb" .OR. outputlist_ps == "cNHb") &
          .AND. lsalsabbins) THEN
         cNHb = FloatArray3d()
         cNHb%onDemand => bulkMixrat
         pipeline => cNHb
         CALL Derived%newField("cNHb", "Bulk NH in clouds B", "kg/kg", 'tttt',   &
                               ANY(outputlist == "cNHb"), pipeline               )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pNHa" .OR. outputlist_ps == "pNHa")) THEN
         pNHa = FloatArray3d()
         pNHa%onDemand => bulkMixrat
         pipeline => pNHa
         CALL Derived%newField("pNHa", "Bulk NH in precip", "kg/kg", 'tttt',   &
                               ANY(outputlist == "pNHa"), pipeline             )
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iNHa" .OR. outputlist_ps == "iNHa")) THEN
         iNHa = FloatArray3d()
         iNHa%onDemand => bulkMixrat
         pipeline => iNHa
         CALL Derived%newField("iNHa", "Bulk NH in ice", "kg/kg", 'tttt',   &
                               ANY(outputlist == "iNHa"), pipeline          )
      END IF
         
      pipeline => NULL()
      
    END SUBROUTINE setDerivedVariables
      

  
  

END MODULE mo_derived_state
