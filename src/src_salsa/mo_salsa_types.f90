MODULE mo_salsa_types
  USE classProcessRates, ONLY : ProcessRates
  USE classSection, ONLY : Section
  IMPLICIT NONE

  SAVE
  
  ! This module is a container for SALSA datatypes.
  ! Previously most of the stuff was placed in mo_submctl.
  ! Replacing them here help avoiding cyclic dependencies and
  ! improves the code structure.
  
  ! All particle properties for SALSA. All the setup and pointer associations will be done in mo_aero_init
  TYPE(Section), ALLOCATABLE, TARGET :: allSALSA(:,:,:)   ! Parent array holding all particle and hydrometeor types consecutively 
  
  ! Particle type specific pointers to "allSALSA" master array defined in mo_salsa_driver.
  ! Pointer association is done in mo_salsa_init. These should be accessed by importin mo_submctl,
  ! not by dummy arguments.
  TYPE(Section), POINTER :: aero(:,:,:)  => NULL(),   &
                            cloud(:,:,:) => NULL(),  &
                            precp(:,:,:) => NULL(),  &
                            ice(:,:,:)   => NULL(),    &
                            liquid(:,:,:) => NULL(),  &
                            frozen(:,:,:) => NULL()

  ! Star and end indices for different particle types in the allSALSA array
  INTEGER :: iaero, faero, icloud, fcloud, iprecp, fprecp, iice, fice
  
  ! Process specific pointers for process rate containers. Initialize in salsa_init, remember to reset for each cycle in salsa_driver
  TYPE(ProcessRates) :: rateDiag     

  
END MODULE mo_salsa_types
