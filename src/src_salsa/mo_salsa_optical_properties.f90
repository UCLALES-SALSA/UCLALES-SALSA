
MODULE mo_salsa_optical_properties
  USE classSpecies, ONLY : maxspec
  USE mo_submctl, ONLY : spec
  IMPLICIT NONE

  ! Refractive indices and the corresponding wavelengths
  ! Shortwave
  INTEGER, PARAMETER :: NaerRadPropsSW = 13
  REAL, PARAMETER :: aerRefrIBands_SW(NaerRadPropsSW) =            &
                                    (/3.46, 2.79, 2.33, 2.05,      &
                                    1.78, 1.46, 1.27, 1.01,        &
                                    0.70, 0.53, 0.39, 0.30,        &
                                    0.23/)*1.e-4  ! cm
  REAL, PARAMETER :: riReSWsu(NaerRadPropsSW) =                    &
                                    (/1.361, 1.295, 1.364, 1.382,  &
                                     1.393, 1.406, 1.413, 1.422,   &
                                     1.427, 1.432, 1.445, 1.450,   &
                                     1.450/),                      &
                     riImSWsu(NaerRadPropsSW) =                    &
                                    (/1.400E-01, 5.500E-02, 2.100E-03, 1.300E-03,   &
                                      5.100E-04, 9.000E-05, 7.900E-06, 1.300E-06,   &
                                      5.200E-08, 1.000E-09, 1.000E-09, 1.000E-09,   &
                                      1.000E-09/)
  REAL, PARAMETER :: riReSWbc(NaerRadPropsSW) =                      &
                                    (/1.984, 1.936, 1.917, 1.905,    &
                                      1.894, 1.869, 1.861, 1.861,    &
                                      1.850, 1.850, 1.839, 1.839,    &
                                      1.713/),                       &
                     riImSWbc(NaerRadPropsSW) =                      &
                                    (/8.975E-01, 8.510E-01, 8.120E-01, 7.939E-01,   &
                                      7.765E-01, 7.397E-01, 7.274E-01, 7.106E-01,   &
                                      6.939E-01, 7.213E-01, 7.294E-01, 7.584E-01,   &
                                      7.261E-01/)
  REAL, PARAMETER :: riReSWoc(NaerRadPropsSW) =                      &
                                    (/1.530, 1.510, 1.510, 1.420,    &
                                      1.464, 1.520, 1.420, 1.420,    &
                                      1.530, 1.530, 1.530, 1.443,    &
                                      1.530/),                    &
                     riImSWoc(NaerRadPropsSW) =                      &
                                    (/2.75E-02, 7.33E-03, 7.33E-03, 4.58E-03,    &
                                      6.42E-03, 1.43E-02, 1.77E-02, 2.01E-02,    &
                                      1.50E-02, 7.70E-03, 9.75E-03, 1.63E-02,    &
                                      5.27E-03/)
  REAL, PARAMETER :: riReSWss(NaerRadPropsSW) =                      &
                                    (/1.480, 1.400, 1.440, 1.450,    &
                                      1.450, 1.460, 1.470, 1.470,    &
                                      1.480, 1.490, 1.500, 1.510,    &          
                                      1.510/),                       &
                     riImSWss(NaerRadPropsSW) =                      &
                                    (/1.300E-02, 8.000E-03, 2.500E-03, 1.500E-03, &
                                      1.000E-03, 5.500E-04, 3.300E-04, 1.000E-04, &
                                      1.000E-07, 1.000E-08, 2.000E-08, 1.000E-06, &
                                      1.000E-05/)
  REAL, PARAMETER :: riReSWdu(NaerRadPropsSW) =                      &
                                    (/1.460, 1.460, 1.460, 1.450,    &
                                      1.450, 1.450, 1.450, 1.450,    &
                                      1.450, 1.450, 1.450, 1.450,    &
                                      1.450/),                       &
                     riImSWdu(NaerRadPropsSW) =                      &
                                    (/1.180E-02, 6.000E-03, 2.500E-03, 1.500E-03, &
                                      1.000E-03, 8.000E-04, 6.000E-04, 7.500E-04, &
                                      9.500E-04, 1.000E-03, 2.500E-03, 2.000E-02, &
                                      2.500E-02/)
  ! VALUES FOR HNO3 AND NH3 NOT SET YET!!!
  ! --------------------------------------------------------------------
  REAL, PARAMETER :: riReSWhno(NaerRadPropsSW) =                       &
                                    (/0.0, 0.0, 0.0, 0.0,     &
                                      0.0, 0.0, 0.0, 0.0,     &
                                      0.0, 0.0, 0.0, 0.0,     &
                                      0.0/),                  &
                     riImSWhno(NaerRadPropsSW) =                       &
                                    (/0.0, 0.0, 0.0, 0.0,                 &
                                      0.0, 0.0, 0.0, 0.0,                 &
                                      0.0, 0.0, 0.0, 0.0,                 &
                                      0.0/)

  REAL, PARAMETER :: riReSWnh(NaerRadPropsSW) =                       &
                                    (/0.0, 0.0, 0.0, 0.0,     &
                                      0.0, 0.0, 0.0, 0.0,     &
                                      0.0, 0.0, 0.0, 0.0,     &
                                      0.0/),                  &
                     riImSWnh(NaerRadPropsSW) =                       &
                                    (/0.0, 0.0, 0.0, 0.0,                 &
                                      0.0, 0.0, 0.0, 0.0,                 &
                                      0.0, 0.0, 0.0, 0.0,                 &
                                      0.0/)
  ! ----------------------------------------------------------------------

  REAL, PARAMETER :: riReSWh2o(NaerRadPropsSW) =                     &
                                    (/1.423, 1.244, 1.283, 1.300,    &
                                      1.312, 1.319, 1.324, 1.328,    &
                                      1.331, 1.335, 1.341, 1.350,    &
                                      1.377/),                       &
                     riImSWh2o(NaerRadPropsSW) =                     &
                                    (/5.000E-02, 1.300E-01, 6.500E-04, 6.700E-04, &
                                      1.200E-04, 1.100E-04, 1.200E-05, 2.100E-06, &
                                      6.800E-08, 2.800E-09, 3.900E-09, 1.700E-08, &
                                      6.400E-08/)

  ! Longwave
  INTEGER, PARAMETER :: NaerRadPropsLW = 16
  REAL, PARAMETER :: aerRefrIBands_LW(NaerRadPropsLW) =               &
                                    (/55.56, 23.53, 17.70, 15.04,     &
                                      13.16, 11.11, 9.71,  8.85,      &
                                      7.78,  6.97,  6.10,  5.15,      &
                                      4.62,  4.32,  4.02,  3.42/)*1.e-4  !cm

  REAL, PARAMETER :: riReLWsu(NaerRadPropsLW) =                       &
                                    (/1.889, 1.588, 1.804, 1.537,     &
                                      1.709, 1.879, 2.469, 0.685,     &
                                      1.427, 0.956, 1.336, 1.450,     &
                                      1.489, 1.512, 1.541, 1.602/),   &
                     riImLWsu(NaerRadPropsLW) =                       &
                                    (/0.967E-01, 0.380E-01, 0.287E-01, 0.225E-01,    &
                                      0.200E-01, 0.396E-01, 0.269E+00, 0.111E+01,    &
                                      0.705E-01, 0.678E+00, 0.143E-01, 0.664E-02,    &
                                      0.657E-02, 0.944E-02, 0.148E-01, 0.156E+00/)
  REAL, PARAMETER :: riReLWbc(NaerRadPropsLW) =                       &
                                    (/2.84, 2.63, 2.53, 2.46,         &
                                      2.42, 2.36, 2.33, 2.30,         &
                                      2.23, 2.17, 2.14, 2.09,         &
                                      2.06, 2.04, 2.03, 1.98/),       &
                     riImLWbc(NaerRadPropsLW) =                       &
                                    (/1.61E+00, 1.42E+00, 1.33E+00, 1.28E+00,         &
                                      1.23E+00, 1.18E+00, 1.16E+00, 1.14E+00,         &
                                      1.08E+00, 1.04E+00, 1.00E+00, 9.73E-01,         &
                                      9.56E-01, 9.46E-01, 9.37E-01, 8.91E-01/)
  REAL, PARAMETER :: riReLWoc(NaerRadPropsLW) =                       &
                                    (/1.86, 1.95, 2.02, 1.43,         &
                                      1.61, 1.71, 1.81, 2.64,         &
                                      1.23, 1.42, 1.42, 1.45,         &
                                      1.46, 1.46, 1.46, 1.44/),       &
                     riImLWoc(NaerRadPropsLW) =                       &
                                    (/4.58E-01, 2.35E-01, 1.86E-01, 1.82E-01,         &
                                      5.31E-02, 4.52E-02, 4.54E-02, 3.76E-01,         &
                                      6.04E-02, 5.30E-02, 2.29E-02, 1.27E-02,         &
                                      1.17E-02, 9.28E-03, 4.88E-03, 5.95E-03/)
  REAL, PARAMETER :: riReLWss(NaerRadPropsLW) =                       &
                                    (/1.668, 1.749, 1.763, 1.447,     &
                                      1.408, 1.485, 1.563, 1.638,     &
                                      1.401, 1.450, 1.505, 1.459,     &
                                      1.483, 1.488, 1.478, 1.484/),   &
                     riImLWss(NaerRadPropsLW) =                       &
                                    (/0.981E+00, 0.193E+00, 0.111E+00, 0.344E-01,     &
                                      0.192E-01, 0.140E-01, 0.179E-01, 0.293E-01,     &
                                      0.138E-01, 0.543E-02, 0.180E-01, 0.288E-02,     &
                                      0.251E-02, 0.246E-02, 0.175E-02, 0.206E-02/)
  REAL, PARAMETER :: riReLWdu(NaerRadPropsLW) =                       &
                                    (/2.552, 2.552, 1.865, 1.518,     &
                                      1.697, 1.816, 2.739, 1.613,     &
                                      1.248, 1.439, 1.423, 1.526,     &
                                      1.502, 1.487, 1.480, 1.468/),   &
                     riImLWdu(NaerRadPropsLW) =                       &
                                    (/0.7412, 0.7412, 0.5456, 0.2309,                 &
                                      0.1885, 0.2993, 0.7829, 0.4393,                 &
                                      0.1050, 0.0976, 0.0540, 0.0228,                 &
                                      0.0092, 0.0053, 0.0044, 0.0101/)

  ! VALUES FOR HNO3 AND NH3 NOT SET YET!!!
  ! --------------------------------------------------------------------
  REAL, PARAMETER :: riReLWhno(NaerRadPropsLW) =                       &
                                    (/0.0, 0.0, 0.0, 0.0,     &
                                      0.0, 0.0, 0.0, 0.0,     &
                                      0.0, 0.0, 0.0, 0.0,     &
                                      0.0, 0.0, 0.0, 0.0/),   &
                     riImLWhno(NaerRadPropsLW) =                       &
                                    (/0.0, 0.0, 0.0, 0.0,                 &
                                      0.0, 0.0, 0.0, 0.0,                 &
                                      0.0, 0.0, 0.0, 0.0,                 &
                                      0.0, 0.0, 0.0, 0.0/)

  REAL, PARAMETER :: riReLWnh(NaerRadPropsLW) =                       &
                                    (/0.0, 0.0, 0.0, 0.0,     &
                                      0.0, 0.0, 0.0, 0.0,     &
                                      0.0, 0.0, 0.0, 0.0,     &
                                      0.0, 0.0, 0.0, 0.0/),   &
                     riImLWnh(NaerRadPropsLW) =                       &
                                    (/0.0, 0.0, 0.0, 0.0,                 &
                                      0.0, 0.0, 0.0, 0.0,                 &
                                      0.0, 0.0, 0.0, 0.0,                 &
                                      0.0, 0.0, 0.0, 0.0/)
  ! ----------------------------------------------------------------------

  REAL, PARAMETER :: riReLWh2o(NaerRadPropsLW) =                      &
                                    (/1.689, 1.524, 1.401, 1.283,     &
                                      1.171, 1.149, 1.230, 1.264,     &
                                      1.295, 1.314, 1.312, 1.316,     &
                                      1.327, 1.333, 1.348, 1.416/),   &
                     riImLWh2o(NaerRadPropsLW) =                      &
                                    (/0.618E+00, 0.392E+00, 0.428E+00, 0.395E+00,     &
                                      0.317E+00, 0.107E+00, 0.481E-01, 0.392E-01,     &
                                      0.347E-01, 0.348E-01, 0.132E+00, 0.106E-01,     &
                                      0.151E-01, 0.881E-02, 0.483E-02, 0.169E-01/)

  ! allocatable arrays for the chemicals configured active in SALSA (order as in the NAMELIST and mass arrays)
  ! Shape = (number of species, number of spectral bands)
  REAL, ALLOCATABLE :: riReSW(:,:), riImSW(:,:), riReLW(:,:), riImLW(:,:)

  CONTAINS

    SUBROUTINE initialize_optical_properties()
      IMPLICIT NONE
      
      INTEGER :: nspec, ii, jj

      ! Refractive indices for all available chemicals (order as in classSpecies/allNames) 
      REAL :: riReSW_all(maxspec+1,NaerRadPropsSW), riImSW_all(maxspec+1,NaerRadPropsSW),   &
              riReLW_all(maxspec+1,NaerRadPropsLW), riImLW_all(maxspec+1,NaerRadPropsLW)
      
      ! The number of active compounds defined in the NAMELIST (+water)
      nspec = spec%getNSpec()

      ! This gives the optical properties for all possible SALSA compounds in single arrays, hard coded in the same order
      ! as classSpecies/allNames. This is mostly needed because it makes it easy to use the indexing provided by the class object "spec" to
      ! construct the trucated arrays below
      riReSW_all(:,:) = RESHAPE( [riReSWsu,riReSWoc,riReSWbc,riReSWdu,riReSWss,riReSWhno,riReSWnh,riReSWh2o], &
                                 SHAPE(riReSW_all), ORDER=[2,1]                                               )
      riReLW_all(:,:) = RESHAPE( [riReLWsu,riReLWoc,riReLWbc,riReLWdu,riReLWss,riReLWhno,riReLWnh,riReLWh2o], &
                                 SHAPE(riReLW_all), ORDER=[2,1]                                               )
      riImSW_all(:,:) = RESHAPE( [riImSWsu,riImSWoc,riImSWbc,riImSWdu,riImSWss,riImSWhno,riImSWnh,riImSWh2o], &
                                 SHAPE(riImSW_all), ORDER=[2,1]                                               )
      riImLW_all(:,:) = RESHAPE( [riImLWsu,riImLWoc,riImLWbc,riImLWdu,riImLWss,riImLWhno,riImLWnh,riImLWh2o], &
                                 SHAPE(riImLW_all), ORDER=[2,1]                                               )
      
      ! Allocate and construct the truncated optical property arrays (contain only active compounds and given in the order in which they
      ! are given in the NAMELIST and the LES/SALSA mass arrays)
      ALLOCATE( riReSW(nspec,NaerRadPropsSW), riImSW(nspec,NaerRadPropsSW),  &
                riReLW(nspec,NaerRadPropsLW), riImLW(nspec,NaerRadPropsLW)   )
      
      riReSW = 0.; riImSW = 0.; riReLW = 0.; riImLW = 0.

      ! Sort the properties
      DO ii = 1,maxspec+1
         jj = spec%allInd(ii)
         IF (jj == 0) CYCLE  ! Species not used
         ! DEBUGGING
         IF (jj > nspec) STOP "optical props vaarin meni!!!"
         riReSW(jj,:) = riReSW_all(ii,:)
         riImSW(jj,:) = riImSW_all(ii,:)
         riReLW(jj,:) = riReLW_all(ii,:)
         riImLW(jj,:) = riImLW_all(ii,:)
      END DO
         

    END SUBROUTINE initialize_optical_properties


END MODULE mo_salsa_optical_properties
