MODULE classBinLayout
  ! This class is a container for parameters used to define 
  ! the bin setup for the particle categories in SALSA.
  ! 
  ! FOR NOW (@ 20.06.2018) ONLY USED FOR PRECIPITATION BINS
  ! BUT WITH SOME ADJUSTMENT OF THE BASIC LOGIC IN AEROSOL AND 
  ! CLOUDS, CAN BE EASILY USED FOR THE REST AS WELL.

  IMPLICIT NONE

  TYPE BinLayout
     INTEGER :: nbins       ! Number of bins
     REAL    :: dlo         ! Diameter of the low limit of the first bin
     REAL    :: vol_ratio   ! Spacing of the bins given as the ratio of the volume (~mass for hydrometeors)
                            ! between consecutive bins. Assuming spherical particles!!
  END TYPE BinLayout
  INTERFACE binLayout
     PROCEDURE :: cnstr
  END INTERFACE binLayout

  CONTAINS

    FUNCTION cnstr(iN,idlo,ivr)
      TYPE(BinLayout) :: cnstr
      INTEGER, INTENT(in) :: iN
      REAL, INTENT(in)    :: idlo
      REAL, INTENT(in)    :: ivr

      cnstr%nbins      = iN
      cnstr%dlo        = idlo
      cnstr%vol_ratio = ivr

    END FUNCTION cnstr



END MODULE classBinLayout
