MODULE classProcessRates
  IMPLICIT NONE

  SAVE
    
  TYPE ProcessRates
     TYPE(Rate), POINTER :: Autoconversion
     TYPE(Rate), POINTER :: Autoconversion50 ! For drops growing past 50um
     TYPE(Rate), POINTER :: Autoconversion80 ! For drops growing past 80um
     TYPE(Rate), POINTER :: Accretion
     TYPE(Rate), POINTER :: Accretion50      ! By drops L.T. 50 um
     TYPE(Rate), POINTER :: Accretion80      ! By drops L.T. 80 um
     TYPE(Rate), POINTER :: ACcoll
     TYPE(Rate), POINTER :: APcoll
     TYPE(Rate), POINTER :: AIcoll     
     TYPE(Rate), POINTER :: Activation
     TYPE(Rate), POINTER :: Ice_hom
     TYPE(Rate), POINTER :: Ice_imm
     TYPE(Rate), POINTER :: Ice_dep
     TYPE(Rate), POINTER :: Cond_a
     TYPE(Rate), POINTER :: Cond_c
     TYPE(Rate), POINTER :: Cond_p
     TYPE(rate), POINTER :: Cond_i
     CONTAINS
       PROCEDURE :: Reset
  END TYPE ProcessRates
  INTERFACE ProcessRates
     PROCEDURE :: ProcessRates_cnstr
  END INTERFACE ProcessRates
  
  TYPE Rate
     ! Holds a process rate for volume and number concentrations, i.e. m3/m3s and #/m3s
     REAL :: volc(9)
     REAL :: numc
     CONTAINS
       PROCEDURE :: Reset_R
       PROCEDURE :: Accumulate
  END TYPE Rate
  INTERFACE Rate
     PROCEDURE :: Rate_cnstr
  END INTERFACE Rate

  INTEGER, PARAMETER :: nproc = 17
  TYPE(Rate), TARGET :: Rates(nproc)

  
  ! ------------
  
  CONTAINS

    SUBROUTINE Initialize_processrates()
      INTEGER :: i
      DO i = 1,nproc
         Rates(i) = Rate()
      END DO
    END SUBROUTINE Initialize_processrates
    
    FUNCTION ProcessRates_cnstr()
      TYPE(ProcessRates), TARGET :: ProcessRates_cnstr
      
      ProcessRates_cnstr%Autoconversion => Rates(1)

      ProcessRates_cnstr%Autoconversion50 => Rates(2)

      ProcessRates_cnstr%Autoconversion80 => Rates(3)

      ProcessRates_cnstr%Accretion => Rates(4)

      ProcessRates_cnstr%Accretion50 => Rates(5)
      
      ProcessRates_cnstr%Accretion80 => Rates(6)

      ProcessRates_cnstr%ACcoll => Rates(7)

      ProcessRates_cnstr%APcoll => Rates(8)

      ProcessRates_cnstr%AIcoll => Rates(9)

      ProcessRates_cnstr%Activation => Rates(10)

      ProcessRates_cnstr%Ice_hom => Rates(11)

      ProcessRates_cnstr%Ice_imm => Rates(12)

      ProcessRates_cnstr%Ice_dep => Rates(13)

      ProcessRates_cnstr%Cond_a => Rates(14)

      ProcessRates_cnstr%Cond_c => Rates(15)
      
      ProcessRates_cnstr%Cond_p => Rates(16)
      
      ProcessRates_cnstr%Cond_i => Rates(17)
      
    END FUNCTION ProcessRates_cnstr

    ! --------------------
    
    FUNCTION Rate_cnstr()
      TYPE(Rate) :: Rate_cnstr
      Rate_cnstr%volc(:) = 0.
      Rate_cnstr%numc = 0.
    END FUNCTION Rate_cnstr

    ! ====================
    
    SUBROUTINE Reset(SELF)
      CLASS(ProcessRates), INTENT(inout) :: SELF
      INTEGER :: i
      DO i = 1,nproc
         CALL Rates(i)%Reset_R()
      END DO      
    END SUBROUTINE Reset

    ! -------------------
    
    SUBROUTINE Reset_R(SELF)
      CLASS(Rate), INTENT(inout) :: SELF
      SELF%volc(:) = 0.
      SELF%numc = 0.
    END SUBROUTINE Reset_R

    ! --------------------------------

    SUBROUTINE Accumulate(SELF,n,v)
      CLASS(Rate), INTENT(inout) :: SELF
      REAL, INTENT(in), OPTIONAL :: n
      REAL, INTENT(in), OPTIONAL :: v(:)

      INTEGER :: nv

      IF (PRESENT(v)) THEN
         nv = SIZE(v)
         IF (nv > SIZE(SELF%volc)) STOP "classProcess RATE accumulate ERROR"
         SELF%volc(1:nv) = SELF%volc(1:nv) + v(1:nv)        
      END IF
      IF (PRESENT(n)) &
           SELF%numc = SELF%numc + n

    END SUBROUTINE Accumulate

      
  
END MODULE classProcessRates
