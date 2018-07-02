!****************************************************************
!*                                                              *
!*   MODULE nudg_defs                                           *
!*                                                              *
!*   Ali                                                        *
!*    variables and function of module nudg are moved to        *
!*    this module, to make them available for module grid       *
!*                                                              *
!****************************************************************


MODULE nudg_defs

  IMPLICIT NONE

  TYPE t_nudge
    REAL :: tau_min = 300.      ! Minimum and maximum tau for changing relaxation time schemes (tau_type = 1-3)
    REAL :: tau_max = 300.
    LOGICAL :: tau_max_continue = .FALSE. ! Continue nudging with tau_max after the actual nudging time period?
                                            ! Valid for tau_type = 1-3
    INTEGER :: tau_type = 0     ! 0: constant, 1: linear increase,   
                                  ! 2: negative exponential increase, 3: positive exponential increase
    INTEGER :: nudgetype = 0     ! Nudging options (nudge_*: 0=disabled, 1=soft, 2=hard)
  CONTAINS

    PROCEDURE :: f_tau

  END TYPE t_nudge

  REAL, ALLOCATABLE :: theta_ref(:), rv_ref(:), u_ref(:), v_ref(:), aero_ref(:,:), aero_target(:,:)
  TYPE(t_nudge) :: ndg_theta, ndg_rv, ndg_u, ndg_v, ndg_aero

   REAL, PARAMETER :: global_tau_min = 300.

   REAL  :: nudge_time = 3600., nudge_zmin = -1.e10, nudge_zmax = 1.e10

 CONTAINS

  !
  ! Procedures bound to t_nudge
  ! ------------------------------
  !  
  FUNCTION f_tau(SELF,t)
    IMPLICIT NONE
    CLASS(t_nudge) :: SELF
    REAL, INTENT(in) :: t
    REAL :: f_tau

    CHARACTER(len=50), PARAMETER :: name = "f_tau"

    SELECT CASE(SELF%tau_type)
       CASE(0)
          ! For constant tau use always TAU_MIN
          f_tau = SELF%tau_min
       CASE(1)
          f_tau = linear(t)
       CASE(2)
          f_tau = negative_exponential(t)
       CASE(3)
          f_tau = positive_exponential(t)
        END SELECT
        
  CONTAINS
      
      FUNCTION linear(tt)
        REAL, INTENT(in) :: tt
        CHARACTER(len=50), PARAMETER :: name = "f_tau/linear"
        REAL :: linear
        REAL :: hlp,ttloc
        ttloc = MIN(tt,nudge_time) ! If tt > nudge_time and you end here, then tau_max_continue == TRUE
        hlp = ttloc*(SELF%tau_max - SELF%tau_min)/nudge_time
        linear = MAX(SELF%tau_min + hlp, global_tau_min)
      END FUNCTION linear

      FUNCTION negative_exponential(tt)
        REAL, INTENT(in) :: tt
        CHARACTER(len=50), PARAMETER :: name = "f_tau/negative_exponential"
        REAL :: negative_exponential
        REAL :: hlp, hlp2
        REAL, PARAMETER :: z = -1.0067837
        REAL :: ttloc
        ttloc = MIN(tt,nudge_time) ! If tt > nudge_time and you end here, then tau_max_continue == TRUE
        hlp = (SELF%tau_max - SELF%tau_min) * z
        hlp2 = EXP(-5.*ttloc/nudge_time) - 1.
        negative_exponential = MAX(SELF%tau_min + hlp*hlp2, global_tau_min)
      END FUNCTION negative_exponential

      FUNCTION positive_exponential(tt)
        REAL, INTENT(in) :: tt
        CHARACTER(len=50), PARAMETER :: name = "f_tau/positive_exponential"
        REAL :: positive_exponential
        REAL :: hlp, hlp2
        REAL, PARAMETER :: z = 0.0067837
        REAL :: ttloc
        ttloc = MIN(tt,nudge_time) ! If tt > nudge_time and you end here, then tau_max_continue == TRUE
        hlp = (SELF%tau_max - SELF%tau_min) * z
        hlp2 = EXP(5.*ttloc/nudge_time) - 1.
        positive_exponential = MAX(SELF%tau_min + hlp*hlp2, global_tau_min)
      END FUNCTION positive_exponential 
        
  END FUNCTION f_tau        
        
END MODULE nudg_defs
