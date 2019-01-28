MODULE mo_derived_state
  USE classFieldArray
  USE mo_structured_datatypes
  USE mo_derived_procedures
  IMPLICIT NONE

  SAVE

  !
  ! Derived diagnostic quantities
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

  ! Since the quantities in this module are only calculated for output, we don't really need to save any of the data.
  ! Therefore, make all the arrays to point to this dummy array. The subroutines in mo_derived_procedures will provide
  ! the output arrays on on-demand basis. However, if some of these need to be saved in the future, that's an easy change.
  ! This is just to reduce the amount of allocated memory
  REAL, ALLOCATABLE, TARGET :: zeros3d(:,:,:)
  
  CONTAINS

    SUBROUTINE setDerivedVariables(Derived,outputlist,level,nzp,nxp,nyp)
      TYPE(FieldArray), INTENT(inout) :: Derived
      CHARACTER(len=10), INTENT(in) :: outputlist(:)
      INTEGER, INTENT(in) :: level, nzp,nxp,nyp
      CLASS(*), POINTER :: pipeline

      ALLOCATE(zeros3d(nzp,nxp,nyp))
      zeros3d = 0.

      ! Since these are only for output, mask their allocation with the output status

      IF (ANY(outputlist == "qtot")) THEN
         qtot = FloatArray3d(zeros3d, store=.FALSE.)
         qtot%onDemand => totalWater
         pipeline => qtot
         CALL Derived%newField("qtot", "Total water content", "kg/kg", "tttt",   &
              ANY(outputlist == "qtot"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Naa")) THEN
         Naa = FloatArray3d(zeros3d,store=.FALSE.)
         Naa%onDemand => bulkNumc
         pipeline => Naa
         CALL Derived%newField("Naa", "Bulk number of aerosol, A", "m-3", 'tttt',   &
              ANY(outputlist == "Naa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Nab")) THEN
         Nab = FloatArray3d(zeros3d,store=.FALSE.)
         Nab%onDemand => bulkNumc
         pipeline => Nab
         CALL Derived%newField("Nab", "Bulk number of aerosol, B", "m-3", 'tttt',   &
              ANY(outputlist == "Nab"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Nca")) THEN
         Nca = FloatArray3d(zeros3d,store=.FALSE.)
         Nca%onDemand => bulkNumc
         pipeline => Nca
         CALL Derived%newField("Nca", "Bulk number of cloud droplets, A", "m-3", 'tttt',   &
              ANY(outputlist == "Nca"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Ncb")) THEN         
         Ncb = FloatArray3d(zeros3d,store=.FALSE.)
         Ncb%onDemand => bulkNumc
         pipeline => Ncb
         CALL Derived%newField("Ncb", "Bulk number of cloud droplets, B", "m-3", 'tttt',   &
              ANY(outputlist == "Ncb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Np")) THEN      
         Np = FloatArray3d(zeros3d,store=.FALSE.)
         Np%onDemand => bulkNumc
         pipeline => Np
         CALL Derived%newField("Np", "Bulk number of precip", "m-3", 'tttt',   &
              ANY(outputlist == "Np"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Ni")) THEN
         Ni = FloatArray3d(zeros3d,store=.FALSE.)
         Ni%onDemand => bulkNumc
         pipeline => Ni
         CALL Derived%newField("Ni", "Bulk number of ice", "m-3", 'tttt',   &
              ANY(outputlist == "Ni"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwaa")) THEN
         Dwaa = FloatArray3d(zeros3d,store=.FALSE.)
         Dwaa%onDemand => bulkDiameter
         pipeline => Dwaa
         CALL Derived%newField("Dwaa", "Bulk mean diameter, aerosol A", "m", 'tttt',   &
              ANY(outputlist == "Dwaa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwab")) THEN
         Dwab = FloatArray3d(zeros3d,store=.FALSE.)
         Dwab%onDemand => bulkDiameter
         pipeline => Dwab
         CALL Derived%newField("Dwab", "Bulk mean diameter, aerosol B", "m", 'tttt',   &
              ANY(outputlist == "Dwab"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwca")) THEN
         Dwca = FloatArray3d(zeros3d,store=.FALSE.)
         Dwca%onDemand => bulkDiameter
         pipeline => Dwca
         CALL Derived%newField("Dwca", "Bulk mean diameter, clouds A", "m", 'tttt',   &
              ANY(outputlist == "Dwca"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwcb")) THEN
         Dwcb = FloatArray3d(zeros3d,store=.FALSE.)
         Dwcb%onDemand => bulkDiameter
         pipeline => Dwcb
         CALL Derived%newField("Dwcb", "Bulk mean diameter, clouds B", "m", 'tttt',   &
              ANY(outputlist == "Dwcb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwpa")) THEN
         Dwpa = FloatArray3d(zeros3d,store=.FALSE.)
         Dwpa%onDemand => bulkDiameter
         pipeline => Dwpa
         CALL Derived%newField("Dwpa", "Bulk mean diameter, precip", "m", 'tttt',   &
              ANY(outputlist == "Naa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "Dwia")) THEN
         Dwia = FloatArray3d(zeros3d,store=.FALSE.)
         Dwia%onDemand => bulkDiameter
         pipeline => Dwia
         CALL Derived%newField("Dwia", "Bulk mean diameter, ice", "m", 'tttt',   &
              ANY(outputlist == "Dwia"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aSO4a")) THEN
         aSO4a = FloatArray3d(zeros3d,store=.FALSE.)
         aSO4a%onDemand => bulkMixrat
         pipeline => aSO4a
         CALL Derived%newField("aSO4a", "Bulk SO4 in aerosol A", "kg/kg", 'tttt',   &
              ANY(outputlist == "aSO4a"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aSO4b")) THEN
         aSO4b = FloatArray3d(zeros3d,store=.FALSE.)
         aSO4b%onDemand => bulkMixrat
         pipeline => aSO4b
         CALL Derived%newField("aSO4b", "Bulk SO4 in aerosol B", "kg/kg", 'tttt',   &
              ANY(outputlist == "aSO4b"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cSO4a")) THEN
         cSO4a = FloatArray3d(zeros3d,store=.FALSE.)
         cSO4a%onDemand => bulkMixrat
         pipeline => cSO4a
         CALL Derived%newField("cSO4a", "Bulk SO4 in clouds A", "kg/kg", 'tttt',   &
              ANY(outputlist == "cSO4a"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cSO4b")) THEN
         cSO4b = FloatArray3d(zeros3d,store=.FALSE.)
         cSO4b%onDemand => bulkMixrat
         pipeline => cSO4b
         CALL Derived%newField("cSO4b", "Bulk SO4 in clouds B", "kg/kg", 'tttt',   &
              ANY(outputlist == "cSO4b"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pSO4a")) THEN
         pSO4a = FloatArray3d(zeros3d,store=.FALSE.)
         pSO4a%onDemand => bulkMixrat
         pipeline => pSO4a
         CALL Derived%newField("pSO4a", "Bulk SO4 in precip", "kg/kg", 'tttt',   &
              ANY(outputlist == "pSO4a"), pipeline)
      END IF
      
      IF (level >= 4 .AND. ANY(outputlist == "iSO4a")) THEN      
         iSO4a = FloatArray3d(zeros3d,store=.FALSE.)
         iSO4a%onDemand => bulkMixrat
         pipeline => iSO4a
         CALL Derived%newField("iSO4a", "Bulk SO4 in ice", "kg/kg", 'tttt',   &
              ANY(outputlist == "iSO4a"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aOCa")) THEN
         aOCa = FloatArray3d(zeros3d,store=.FALSE.)
         aOCa%onDemand => bulkMixrat
         pipeline => aOCa
         CALL Derived%newField("aOCa", "Bulk OC in aerosol A", "kg/kg", 'tttt',   &
              ANY(outputlist == "aOCa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aOCb")) THEN
         aOCb = FloatArray3d(zeros3d,store=.FALSE.)
         aOCb%onDemand => bulkMixrat
         pipeline => aOCb
         CALL Derived%newField("aOCb", "Bulk OC in aerosol B", "kg/kg", 'tttt',   &
              ANY(outputlist == "aOCb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cOCa")) THEN
         cOCa = FloatArray3d(zeros3d,store=.FALSE.)
         cOCa%onDemand => bulkMixrat
         pipeline => cOCa
         CALL Derived%newField("cOCa", "Bulk OC in clouds A", "kg/kg", 'tttt',   &
              ANY(outputlist == "cOCa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cOCb")) THEN
         cOCb = FloatArray3d(zeros3d,store=.FALSE.)
         cOCb%onDemand => bulkMixrat
         pipeline => cOCb
         CALL Derived%newField("cOCb", "Bulk OC in clouds B", "kg/kg", 'tttt',   &
              ANY(outputlist == "cOCb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pOCa")) THEN
         pOCa = FloatArray3d(zeros3d,store=.FALSE.)
         pOCa%onDemand => bulkMixrat
         pipeline => pOCa
         CALL Derived%newField("pOCa", "Bulk OC in precip", "kg/kg", 'tttt',   &
              ANY(outputlist == "pOCa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iOCa")) THEN
         iOCa = FloatArray3d(zeros3d,store=.FALSE.)
         iOCa%onDemand => bulkMixrat
         pipeline => iOCa
         CALL Derived%newField("iOCa", "Bulk OC in ice", "kg/kg", 'tttt',   &
              ANY(outputlist == "iOCa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aBCa")) THEN
         aBCa = FloatArray3d(zeros3d,store=.FALSE.)
         aBCa%onDemand => bulkMixrat
         pipeline => aBCa
         CALL Derived%newField("aBCa", "Bulk BC in aerosol A", "kg/kg", 'tttt',   &
              ANY(outputlist == "aBCa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aBCb")) THEN
         aBCb = FloatArray3d(zeros3d,store=.FALSE.)
         aBCb%onDemand => bulkMixrat
         pipeline => aBCb
         CALL Derived%newField("aBCb", "Bulk BC in aerosol B", "kg/kg", 'tttt',   &
              ANY(outputlist == "aBCb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cBCa")) THEN
         cBCa = FloatArray3d(zeros3d,store=.FALSE.)
         cBCa%onDemand => bulkMixrat
         pipeline => cBCa
         CALL Derived%newField("cBCa", "Bulk BC in clouds A", "kg/kg", 'tttt',   &
              ANY(outputlist == "cBCa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cBCb")) THEN
         cBCb = FloatArray3d(zeros3d,store=.FALSE.)
         cBCb%onDemand => bulkMixrat
         pipeline => cBCb
         CALL Derived%newField("cBCb", "Bulk BC in clouds B", "kg/kg", 'tttt',   &
              ANY(outputlist == "cBCb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pBCa")) THEN
         pBCa = FloatArray3d(zeros3d,store=.FALSE.)
         pBCa%onDemand => bulkMixrat
         pipeline => pBCa
         CALL Derived%newField("pBCa", "Bulk BC in precip", "kg/kg", 'tttt',   &
              ANY(outputlist == "pBCa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iBCa")) THEN
         iBCa = FloatArray3d(zeros3d,store=.FALSE.)
         iBCa%onDemand => bulkMixrat
         pipeline => iBCa
         CALL Derived%newField("iBCa", "Bulk BC in ice", "kg/kg", 'tttt',   &
              ANY(outputlist == "iBCa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aDUa")) THEN
         aDUa = FloatArray3d(zeros3d,store=.FALSE.)
         aDUa%onDemand => bulkMixrat
         pipeline => aDUa
         CALL Derived%newField("aDUa", "Bulk DU in aerosol A", "kg/kg", 'tttt',   &
              ANY(outputlist == "aDUa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aDUb")) THEN
         aDUb = FloatArray3d(zeros3d,store=.FALSE.)
         aDUb%onDemand => bulkMixrat
         pipeline => aDUb
         CALL Derived%newField("aDUb", "Bulk DU in aerosol B", "kg/kg", 'tttt',   &
              ANY(outputlist == "aDUb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cDUa")) THEN
         cDUa = FloatArray3d(zeros3d,store=.FALSE.)
         cDUa%onDemand => bulkMixrat
         pipeline => cDUa
         CALL Derived%newField("cDUa", "Bulk DU in clouds A", "kg/kg", 'tttt',   &
              ANY(outputlist == "cDUa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cDUb")) THEN
         cDUb = FloatArray3d(zeros3d,store=.FALSE.)
         cDUb%onDemand => bulkMixrat
         pipeline => cDUb
         CALL Derived%newField("cDUb", "Bulk DU in clouds B", "kg/kg", 'tttt',   &
              ANY(outputlist == "cDUb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pDUa")) THEN
         pDUa = FloatArray3d(zeros3d,store=.FALSE.)
         pDUa%onDemand => bulkMixrat
         pipeline => pDUa
         CALL Derived%newField("pDUa", "Bulk DU in precip", "kg/kg", 'tttt',   &
              ANY(outputlist == "pDUa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iDUa")) THEN
         iDUa = FloatArray3d(zeros3d,store=.FALSE.)
         iDUa%onDemand => bulkMixrat
         pipeline => iDUa
         CALL Derived%newField("iDUa", "Bulk DU in ice", "kg/kg", 'tttt',   &
              ANY(outputlist == "iDUa"), pipeline)      
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aSSa")) THEN
         aSSa = FloatArray3d(zeros3d,store=.FALSE.)
         aSSa%onDemand => bulkMixrat
         pipeline => aSSa
         CALL Derived%newField("aSSa", "Bulk SS in aerosol A", "kg/kg", 'tttt',   &
              ANY(outputlist == "aSSa"), pipeline)
      END IF
      
      IF (level >= 4 .AND. ANY(outputlist == "aSSb")) THEN  
         aSSb = FloatArray3d(zeros3d,store=.FALSE.)
         aSSb%onDemand => bulkMixrat
         pipeline => aSSb
         CALL Derived%newField("aSSb", "Bulk SS in aerosol B", "kg/kg", 'tttt',   &
              ANY(outputlist == "aSSb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cSSa")) THEN
         cSSa = FloatArray3d(zeros3d,store=.FALSE.)
         cSSa%onDemand => bulkMixrat
         pipeline => cSSa
         CALL Derived%newField("cSSa", "Bulk SS in clouds A", "kg/kg", 'tttt',   &
              ANY(outputlist == "cSSa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cSSb")) THEN
         cSSb = FloatArray3d(zeros3d,store=.FALSE.)
         cSSb%onDemand => bulkMixrat
         pipeline => cSSb
         CALL Derived%newField("cSSb", "Bulk SS in clouds B", "kg/kg", 'tttt',   &
              ANY(outputlist == "cSSb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pSSa")) THEN
         pSSa = FloatArray3d(zeros3d,store=.FALSE.)
         pSSa%onDemand => bulkMixrat
         pipeline => pSSa
         CALL Derived%newField("pSSa", "Bulk SS in precip", "kg/kg", 'tttt',   &
              ANY(outputlist == "pSSa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iSSa")) THEN
         iSSa = FloatArray3d(zeros3d,store=.FALSE.)
         iSSa%onDemand => bulkMixrat
         pipeline => iSSa
         CALL Derived%newField("iSSa", "Bulk SS in ice", "kg/kg", 'tttt',   &
              ANY(outputlist == "iSSa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aNOa")) THEN
         aNOa = FloatArray3d(zeros3d,store=.FALSE.)
         aNOa%onDemand => bulkMixrat
         pipeline => aNOa
         CALL Derived%newField("aNOa", "Bulk NO in aerosol A", "kg/kg", 'tttt',   &
              ANY(outputlist == "aNOa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aNOb")) THEN
         aNOb = FloatArray3d(zeros3d,store=.FALSE.)
         aNOb%onDemand => bulkMixrat
         pipeline => aNOb
         CALL Derived%newField("aNOb", "Bulk NO in aerosol B", "kg/kg", 'tttt',   &
              ANY(outputlist == "aNOb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cNOa")) THEN
         cNOa = FloatArray3d(zeros3d,store=.FALSE.)
         cNOa%onDemand => bulkMixrat
         pipeline => cNOa
         CALL Derived%newField("cNOa", "Bulk NO in clouds A", "kg/kg", 'tttt',   &
              ANY(outputlist == "cNOa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cNOb")) THEN
         cNOb = FloatArray3d(zeros3d,store=.FALSE.)
         cNOb%onDemand => bulkMixrat
         pipeline => cNOb
         CALL Derived%newField("cNOb", "Bulk NO in clouds B", "kg/kg", 'tttt',   &
              ANY(outputlist == "cNOb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pNOa")) THEN
         pNOa = FloatArray3d(zeros3d,store=.FALSE.)
         pNOa%onDemand => bulkMixrat
         pipeline => pNOa
         CALL Derived%newField("pNOa", "Bulk NO in precip", "kg/kg", 'tttt',   &
              ANY(outputlist == "pNOa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iNOa")) THEN
         iNOa = FloatArray3d(zeros3d,store=.FALSE.)
         iNOa%onDemand => bulkMixrat
         pipeline => iNOa
         CALL Derived%newField("iNOa", "Bulk NO in ice", "kg/kg", 'tttt',   &
              ANY(outputlist == "iNOa"), pipeline)      
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aNHa")) THEN
         aNHa = FloatArray3d(zeros3d,store=.FALSE.)
         aNHa%onDemand => bulkMixrat
         pipeline => aNHa
         CALL Derived%newField("aNHa", "Bulk NH in aerosol A", "kg/kg", 'tttt',   &
              ANY(outputlist == "aNHa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "aNHb")) THEN
         aNHb = FloatArray3d(zeros3d,store=.FALSE.)
         aNHb%onDemand => bulkMixrat
         pipeline => aNHb
         CALL Derived%newField("aNHb", "Bulk NH in aerosol B", "kg/kg", 'tttt',   &
              ANY(outputlist == "aNHb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cNHa")) THEN
         cNHa = FloatArray3d(zeros3d,store=.FALSE.)
         cNHa%onDemand => bulkMixrat
         pipeline => cNHa
         CALL Derived%newField("cNHa", "Bulk NH in clouds A", "kg/kg", 'tttt',   &
              ANY(outputlist == "cNHa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "cNHb")) THEN
         cNHb = FloatArray3d(zeros3d,store=.FALSE.)
         cNHb%onDemand => bulkMixrat
         pipeline => cNHb
         CALL Derived%newField("cNHb", "Bulk NH in clouds B", "kg/kg", 'tttt',   &
              ANY(outputlist == "cNHb"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "pNHa")) THEN
         pNHa = FloatArray3d(zeros3d,store=.FALSE.)
         pNHa%onDemand => bulkMixrat
         pipeline => pNHa
         CALL Derived%newField("pNHa", "Bulk NH in precip", "kg/kg", 'tttt',   &
              ANY(outputlist == "pNHa"), pipeline)
      END IF

      IF (level >= 4 .AND. ANY(outputlist == "iNHa")) THEN
         iNHa = FloatArray3d(zeros3d,store=.FALSE.)
         iNHa%onDemand => bulkMixrat
         pipeline => iNHa
         CALL Derived%newField("iNHa", "Bulk NH in ice", "kg/kg", 'tttt',   &
              ANY(outputlist == "iNHa"), pipeline)
      END IF
         
      
    END SUBROUTINE setDerivedVariables
      

  
  

END MODULE mo_derived_state
