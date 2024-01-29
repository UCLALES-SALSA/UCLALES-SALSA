MODULE perturbation_forc
  USE defs, ONLY : pi, cp
  USE mo_aux_state, ONLY : xt,yt,zt,pi0,pi1
  USE grid, ONLY : nxp,nyp,nzp,level

  IMPLICIT NONE

  TYPE t_warm_bubble
     ! Set up a 3d sinusoidal temperature perturbation. This is intended to be only
     ! called on a single timestep for the affected area.

     ! NAMELIST parameters
     ! -------------------------------------------------------
     LOGICAL :: switch = .FALSE.              ! Whether to use or not
     REAL :: tstart = 1.                      ! time in seconds when the bubble is created
     REAL :: center(3) = [1000.,0.,0.]        ! z,x,y, meters, determined from 
                                              ! grid point displacements
     REAL :: diameter(3) = [500.,1000.,1000.] ! z,x,y, meters
     REAL :: temp_ampl = 1.5                  ! Kelvins
     ! ------------------------------------------------------------------------------
     ! Other parameters
     LOGICAL :: state                                ! runtime process switch dpending on the logical switch and time specs
     
     CONTAINS
       PROCEDURE :: run => run_warm_bubble
  END TYPE t_warm_bubble

  TYPE t_gaussian_flux_perturbation

     ! RAJOTA ENSIN KOSKEMAAN VAKIOVUO-OLETUSTA!!
     
     ! Set up gaussian surface flux perturbation (sensible and latent). This is intended to
     ! be called over consecutive timesteps over a suitable period.

     ! NAMELIST parameters
     ! ----------------------------------------------------------
     LOGICAL :: switch = .FALSE.        ! Whether to use or not
     REAL :: t_start = 3600.          ! Time in seconds when flux perturbation is started
     REAL :: t_end = 7200.            ! Time in seconds when flux perturbation is stopped
     REAL :: amp_lhf                    ! Latent heat flux amplitude Wm-2. Added to background
     REAL :: amp_shf                    ! Sensible heat flux amplitude Wm-2. Added to background
     REAL :: center(2) = [0.,0.]    ! Center location of the perturbation in meters from domain center
     REAL :: sigma = 2000.          ! Standard deviation of the area distribution of perturbation
     ! -------------------------------------------------------------------------------------------------
     ! Other parameters
     LOGICAL :: state                   ! runtime process switch depending on the logical switch and time specs. 
     
     
     CONTAINS
       PROCEDURE :: run => run_gaussian_flux_perturbation       
  END TYPE t_gaussian_flux_perturbation

  ! Set these from namelist
  TYPE(t_warm_bubble) :: warm_bubble
  TYPE(t_gaussian_flux_perturbation) :: gaussian_flux_perturbation
    
  ! -----------------------------------------------------------------
  ! Contains routines and parameters for setting up a warm bubble in 
  ! the beginning of simulation to initialize deep convection. the
  ! warm bubble initialization will be called instead of the random
  ! perturbations routine.
  ! -----------------------------------------
  
  CONTAINS

    SUBROUTINE run_warm_bubble(SELF)
      USE mo_progn_state, ONLY : a_tp
      IMPLICIT NONE
      CLASS(t_warm_bubble), INTENT(in) :: SELF
      INTEGER :: i,j,k
      REAL :: total_dist_norm
      REAL :: coord(3)

      WRITE(*,*) '------------------------------'
      WRITE(*,*) 'INITIALIZING WARM BUBBLE'
      WRITE(*,*) '------------------------------'

      coord = 0.

      DO j = 3,nyp-2
         DO i = 3,nxp-2
            DO k = 2,nzp
               coord = [zt%d(k),xt%d(i),yt%d(j)]
               total_dist_norm = SQRT( SUM( ((coord-SELF%center)/(0.5*SELF%diameter))**2 ) )
               IF ( total_dist_norm < 1. ) THEN                  
                  ! Add the temperature perturbation to the liquid potential temperature
                  ! (Just a straightforward addition, don't worry about the resulting decrease in RH?)                                    
                  !exner = (pi0%d(k) + pi1%d(k))/cp
                  a_tp%d(k,i,j) = a_tp%d(k,i,j) +   &
                       SIN( (1.-total_dist_norm)*0.5*pi )*SELF%temp_ampl
               END IF                
            END DO
         END DO
      END DO
    END SUBROUTINE run_warm_bubble


    SUBROUTINE run_gaussian_flux_perturbation(SELF)
      !
      ! ----------------------------------------------------------------------------------
      ! Set horizontally gaussian perturbation to sensible and latent surface heat fluxes
      ! centered at specified location with specified amplitudes over a specified period in
      ! time. The perturbations are ADDED to the existing background surface fluxes. Note
      ! that this should only be used with isfctyp=0, i.e. with fixed surface fluxes defined
      ! in the NAMELIST. 
      ! -----------------------------------------------x
      !
      USE math_functions, ONLY : f_gauss_NN
      USE srfc, ONLY : lh_flx, sh_flx
      IMPLICIT NONE
      CLASS(t_gaussian_flux_perturbation), INTENT(in) :: SELF
      INTEGER :: i,j,k
      REAL :: coord(2)
      REAL :: dist_from_center
            
      coord = 0.
      DO j=3,nyp-2
         DO i = 3,nxp-2
            coord = [xt%d(i),yt%d(j)]
            dist_from_center = SQRT( SUM( (coord-SELF%center)**2 ) )
            lh_flx(i,j) = lh_flx(i,j) + SELF%amp_lhf * f_gauss_NN(dist_from_center,SELF%sigma,0.) 
            sh_flx(i,j) = sh_flx(i,j) + SELF%amp_shf * f_gauss_NN(dist_from_center,SELF%sigma,0.)                        
         END DO
      END DO                  
    END SUBROUTINE run_gaussian_flux_perturbation

    
END MODULE perturbation_forc
