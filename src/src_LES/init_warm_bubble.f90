MODULE init_warm_bubble
  USE defs, ONLY : pi, cp
  USE grid, ONLY : nxp,nyp,nzp, xt,yt,zt,  &
                   a_tp, level, pi0, pi1


  IMPLICIT NONE

  ! -----------------------------------------------------------------
  ! Contains routines and parameters for setting up a warm bubble in 
  ! the beginning of simulation to initialize deep convection. the
  ! warm bubble initialization will be called instead of the random
  ! perturbations routine.
  ! -----------------------------------------
  
  REAL :: bubble_center(3) = [1000.,0.,0.]        ! z,x,y, meters, determined from 
                                                  ! grid point displacements
  REAL :: bubble_diameter(3) = [500.,1000.,1000.] ! z,x,y, meters
  REAL :: bubble_temp_ampl = 1.5                  ! Kelvins

  CONTAINS

    SUBROUTINE warm_bubble()
      IMPLICIT NONE

      INTEGER :: i,j,k
      REAL :: total_dist_norm
      REAL :: coord(3)
      REAL :: exner

      WRITE(*,*) '------------------------------'
      WRITE(*,*) 'INITIALIZING WARM BUBBLE'
      WRITE(*,*) '------------------------------'

      coord = 0.

      DO j = 3,nyp-2
         DO i = 3,nxp-2
            DO k = 2,nzp

               coord = [zt(k),xt(i),yt(j)]
               total_dist_norm = SQRT( SUM( ((coord-bubble_center)/(0.5*bubble_diameter))**2 ) )

               !IF ( k<10 ) WRITE(*,*) k,i,j,total_dist_norm, coord-bubble_center

               IF ( total_dist_norm < 1. ) THEN
                  
                  ! Add the temperature perturbation to the liquid potential temperature
                  ! (Just a straightforward addition, don't worry about the resulting decrease in RH?)
                                    
                  exner = (pi0(k) + pi1(k))/cp
                  !WRITE(*,*) k,i,j, exner, SIN( (1.-total_dist_norm)*0.5*pi ),zt(k)
                  a_tp(k,i,j) = a_tp(k,i,j) +   &
                       SIN( (1.-total_dist_norm)*0.5*pi )*bubble_temp_ampl!/exner

               END IF 
               
            END DO
         END DO
      END DO

    END SUBROUTINE warm_bubble
    
END MODULE init_warm_bubble
