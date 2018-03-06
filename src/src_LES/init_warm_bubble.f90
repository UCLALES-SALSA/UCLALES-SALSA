MODULE init_warm_bubble
  USE grid, ONLY :
  USE defs, ONLY : pi
  IMPLICIT NONE

  ! -----------------------------------------------------------------
  ! Contains routines and parameters for setting up a warm bubble in 
  ! the beginning of simulation to initialize deep convection. the
  ! warm bubble initialization will be called instead of the random
  ! perturbations routine. Note that this module will not take care
  ! of providing a suitable input sounding, that's left for the user
  ! to do offline. Input sounding or model grid configuration not designed
  ! for a warm bubble case will most likely just crash the model.
  ! -----------------------------------------------------------------
  
  ! Warm bubble parameters: Define the center location (meters, (0.,0.) is the domain center), 
  ! amplitude (in K) and bubble diameter in the horizontal and vertical. Assume a sinusoidal bubble. 
  REAL :: wb_center(3) = [500.,0.,0.] ! Order: Z,X,Y, unit: meters
  REAL :: wb_ampl      = 1.5          ! unit: K
  REAL :: wb_dia_hor   = 1000.        ! unit: meters
  REAL :: wb_dia_ver   = 500.         ! unit: meters


  CONTAINS

    SUBROUTINE warm_bubble
      IMPLICIT NONE

      REAL :: x_min,x_max,y_min,y_max,z_min,z_max
      REAL :: d2r_hor, d2r_ver
      REAL :: z_r,z_hr2,z_vr2


      ! Diameter to radian scaling factor, taking into account the 0.1 RAD cutoff (from both ends)
      d2r_hor = wb_dia_hor/(pi-0.2)
      d2r_ver = wb_dia_ver/(pi-0.2)


      ! Check if it's a 2d run somehow

      ! The bubble will be essentially half wavelength of a sin-wave in each direction. I.e.
      ! the diameter is normalized to cover 0...pi in radian space. 

      DO j = 3,nyp-2
         DO i = 3,nxp-2
            DO k = 1,nzp

               ! Get the distance
               z_hr2 = ((xt(i) - wb_center(2))/wb_dia_hor)**2 +  &
                       ((yt(j) - wb_center(3))/wb_dia_hor)**2
               z_vr2 = ((zt(i) - wb_center(1))/wb_dia_ver)**2

               z_r = SQRT(z_hr2 + z_vr2)

               IF (z_r < 1.) THEN

                  


            END DO
         END DO
      END DO




    END SUBROUTINE warm_bubble

    



END MODULE init_warm_bubble
