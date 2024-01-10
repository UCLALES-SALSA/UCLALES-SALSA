!===============================================================================!
! UCLALES code from https://github.com/uclales/uclales/blob/master/src/ice_sb.F90
!
!===============================================================================!
!
module mcrp_ice_sb
    use defs, only : pi, rho_w=>rowt, rho_ice=>roice, L_wd=>alvl, L_ed=>alvi, R_d=>Rm, R_l=>Rd
    use thrm, only : esl, esi

    IMPLICIT NONE

    PRIVATE
!
!===============================================================================!
!
! Two-moment bulk microphysics by Axel Seifert
!
! (with contributions from Uli Blahak and Heike Noppel)
!
! Description:
! Provides various modules and subroutines for two-moment bulk microphysics
!
! Current Code Owner: Axel Seifert, DWD
!                     axel.seifert@dwd.de
!
! Language: Fortran 90.
!
! Some code standards or recommendations, at least:
!
! - All new variables/subroutines should be named in English
! - In the future also comments should be written in English, 
!   but temporary use of some German is OK, too.
! - Length of names of subroutines should be <= 20
! - Length of names of variables should be <= 15
! - Length of lines has to be <= 100 including comments,
!   recommended is <= 80 for code without comments. 
! - Temporary modifications for experiments should be marked by
!
!     AS_YYYYMMDD>
!         ! Change for / Bugfix ....
!     <AS_YYYYMMDD
!
!   until they are validated as improvments and made official
!   (with AS, or whatever, being the initials of the guy doing the stuff
!   and YYYYMMDD=year/month/day).
!
!===============================================================================!

!MODULE wolken_konstanten

  TYPE PARTICLE
    CHARACTER(10) :: name
    REAL :: nu    !..width parameter of the distribution
    REAL :: mu    !..exponential parameter of the distribution
    REAL :: x_max !..maximum particle mass
    REAL :: x_min !..minimum particle mass
    REAL :: a_geo !..coeff. geometry
    REAL :: b_geo !..coeff. geometry = 1/3
    REAL :: a_vel !..coeff. fall velocity
    REAL :: b_vel !..coeff. fall velocity
    REAL :: a_ven !..coeff. ventilation param.
    REAL :: b_ven !..coeff. ventilation param.
    REAL :: cap   !..coeff. capacity
  END TYPE PARTICLE

  ! .. allgemeine physikalische Konstanten und Parameter ..

  REAL, PARAMETER :: nu_l = 1.4086e-5    !..Kinem. Visc. von Luft
  REAL, PARAMETER :: D_v  = 3.000e-5     !..Diffusivitaet von Wasserdampf
  REAL, PARAMETER :: K_T  = 2.500e-2     !..Waermeleitfaehigkeit     
  REAL, PARAMETER :: L_ew = L_ed-L_wd    !..Schmelzwaerme
  REAL, PARAMETER :: rho0 = 1.225e+0     !..Norm-Luftdichte
  REAL, PARAMETER :: T_3  = 2.7315e+2    !..Tripelpunkt Wasser
  REAL, PARAMETER :: T_f  = 2.330e+2     !..Bei T < T_f kein Fl.wasser
  REAL, PARAMETER :: c_w  = 4.187e+3     !..spezifische Waerme von Wasser

  ! .. wolkenphysikalische Konstanten und Parameter ..

  REAL, PARAMETER :: N_sc = 0.710        !..Schmidt-Zahl (PK, S.541)
  REAL, PARAMETER :: n_f  = 0.333        !..Exponent von N_sc im Vent-koeff. (PK, S.541)
  REAL, PARAMETER :: m_f  = 0.500        !..Exponent von N_re im Vent-koeff. (PK, S.541)

  REAL, PARAMETER :: e_3  = 6.10780000e2 !..Saettigungsdamppfdruck bei T = T_3

  ! .. Hallett-Mossop
  REAL            :: C_mult     = 3.5e8    !..Koeff. fuer Splintering
  REAL, PARAMETER :: T_mult_min = 265.0    !..Minimale Temp. Splintering
  REAL, PARAMETER :: T_mult_max = 270.0    !..Maximale Temp. Splintering
  REAL, PARAMETER :: T_mult_opt = 268.0    !..Optimale Temp. Splintering

  ! .. collisions
  REAL, PARAMETER :: e_ic  = 0.80              !..max. Eff. fuer ice_cloud_riming
  REAL, PARAMETER :: e_sc  = 0.80              !..max. Eff. fuer snow_cloud_riming
  REAL, PARAMETER :: e_gc  = 1.00              !..max. Eff. fuer graupel_cloud_riming
  REAL, PARAMETER :: e_hc  = 1.00              !..max. Eff. fuer hail_cloud_riming
  REAL, PARAMETER :: e_min = 0.01              !..min. Eff. fuer gc,ic,sc
  REAL, PARAMETER :: alpha_spacefilling = 0.1  !..Raumerfuellungskoeff (max. 0.68)
  REAL, PARAMETER :: ice_s_vel  = 0.00         !..Dispersion der Fallgeschw. 
  REAL, PARAMETER :: snow_s_vel = 0.10         !..Dispersion der Fallgeschw.
  REAL, PARAMETER :: r_shedding = 500.0e-6     !..mittlerer Radius Shedding
  REAL, PARAMETER :: T_shed = 263.2
  REAL, PARAMETER :: q_krit_ii = 1.000e-6 ! q-Schwellenwert fuer ice_selfcollection 
  REAL, PARAMETER :: D_krit_ii = 100.0e-6 ! D-Schwellenwert fuer ice_selfcollection
  REAL, PARAMETER :: D_conv_ii = 75.00e-6 ! D-Schwellenwert fuer ice_selfcollection
  ! .. riming mixing ratio (q) and diameter (D) limits
  REAL :: q_krit_ic = 1.000e-5, D_krit_ic = 150.0e-6 ! ice in ice_cloud_riming
  REAL :: q_krit_ir = 1.000e-5, D_krit_ir = 0.0      ! ice in ice_rain_riming
  REAL :: q_krit_sc = 1.000e-5, D_krit_sc = 150.0e-6 ! snow in snow_cloud_riming
  REAL :: q_krit_sr = 1.000e-5, D_krit_sr = 0.0      ! snow in snow_rain_riming
  REAL :: q_krit_gc = 1.000e-6, D_krit_gc = 100.0e-6 ! graupel in graupel_cloud_riming
  REAL :: q_krit_gr = 1.000e-9 ! graupel in graupel_rain_riming; no diameter limit
  REAL :: q_krit_hc = 1.000e-6, D_krit_hc = 100.0e-6 ! hail in hail_cloud_riming
  REAL :: q_krit_hr = 1.000e-9 ! hail in hail_rain_riming; no diameter limit
  REAL :: q_krit_c  = 1.000e-6 ! cloud in all above; diameter limit is D_krit_c
  REAL :: q_krit_r  = 1.000e-9 ! rain in all above; no diameter limit
  ! .. other riming constants
  REAL, PARAMETER :: q_krit    = 1.000e-9 ! q-Schwellenwert sonst
  REAL, PARAMETER :: D_conv_sg = 200.0e-6 ! D-Schwellenwert
  REAL, PARAMETER :: D_conv_ig = 200.0e-6 ! D-Schwellenwert
  REAL, PARAMETER :: x_conv_g  = 0.100e-9 ! minimale Graupel-/Hagelmasse riming
  REAL, PARAMETER :: D_shed_g  = 3.000e-3 ! D-Schwellenwert fuer graupel_shedding
  REAL, PARAMETER :: D_shed_h  = 5.000e-3 ! D-Schwellenwert fuer hagel_shedding
  REAL, PARAMETER :: D_krit_c  = 10.00e-6 ! D-Schwellenwert fuer cloud_collection
  REAL, PARAMETER :: D_coll_c  = 40.00e-6 ! oberer Wert fuer cloud_coll_eff


  ! Adjustable parameters
  ! =====================

  ! Default cloud, rain, ice, snow, graupel and hail from
  ! 1) Seifert, A. and Beheng, K. D.: A two-moment cloud microphysics parameterization
  !    for mixed-phase clouds. Part 1: Model description, Meteorol. Atmos. Phys., 92,
  !    45-66, DOI: 10.1007/s00703-005-0112-4, 2006.
  ! 2) Seifert, A., Kohler, C., and Beheng, K. D.: Aerosol-cloud-precipitation effects
  !    over Germany as simulated by a convective-scale numerical weather prediction model,
  !    Atmos. Chem. Phys., 12, 709-725, https://doi.org/10.5194/acp-12-709-2012, 2012.
  ! 3) Seifert, A., Blahak, U., and Buhr, R.: On the analytic approximation of bulk
  !    collision rates of non-spherical hydrometeors, Geosci. Model Dev., 7, 463-478,
  !    https://doi.org/10.5194/gmd-7-463-2014, 2014.
  TYPE(PARTICLE) :: cloud = PARTICLE('cloud', &
        1.000000, & !.nu
        1.000000, & !.mu
        2.60e-10, & !.x_max, D=80e-6m
        4.20e-15, & !.x_min, D=2.e-6m
        1.24e-01, & !.a_geo
        0.333333, & !.b_geo
        3.75e+05, & !.a_vel
        0.666667, & !.b_vel
        0.780000, & !.a_ven
        0.308000, & !.b_ven
        2.0)        !.cap

  TYPE(PARTICLE) :: rain = PARTICLE('rain', &
        0.000000, & !.nu
        0.333333, & !.mu
        3.00e-06, & !.x_max
        2.60e-10, & !.x_min, D=80.e-6m
        1.24e-01, & !.a_geo
        0.333333, & !.b_geo
        1.14e+02, & !.a_vel
        0.234000, & !.b_vel
        0.780000, & !.a_ven
        0.308000, & !.b_ven
        2.0)        !.cap

  TYPE(PARTICLE) :: ice =PARTICLE('ice', &
         0.000000, & !.nu
         0.333333, & !.mu
         1.00e-06, & !.x_max
         1.00e-12, & !.x_min, D=200e-6m
         0.835000, & !.a_geo
         0.390000, & !.b_geo
         27.70000, & !.a_vel
         0.216000, & !.b_vel
         0.780000, & !.a_ven
         0.308000, & !.b_ven
         2.0)        !.cap

  TYPE(PARTICLE) :: snow = PARTICLE('snow', &
         1.000000, & !.nu
         0.500000, & !.mu
         2.00e-05, & !.x_max
         1.00e-10, & !.x_min
         2.400000, & !.a_geo
         0.455000, & !.b_geo
         8.800000, & !.a_vel
         0.150000, & !.b_vel
         0.780000, & !.a_ven
         0.308000, & !.b_ven
         2.0)        !.cap

  TYPE(PARTICLE) :: graupel = PARTICLE('graupel' ,&
        1.000000, & !.nu
        0.333333, & !.mu
        5.00e-04, & !.x_max
        1.00e-09, & !.x_min
        0.142000, & !.a_geo
        0.314000, & !.b_geo
        3.30e+01, & !.a_vel
        0.187000, & !.b_vel
        0.780000, & !.a_ven
        0.308000, & !.b_ven
        2.0)        !.cap

  TYPE(PARTICLE) :: hail = PARTICLE('hail' ,&
        1.000000, & !.nu
        0.333333, & !.mu
        1.00e-04, & !.x_max
        2.60e-09, & !.x_min
        1.34e-01, & !.a_geo
        0.333333, & !.b_geo
        3.93e+01, & !.a_vel
        0.167000, & !.b_vel
        0.780000, & !.a_ven
        0.308000, & !.b_ven
        2.0)        !.cap


  ! An alternative method for ice/snow - cloud/rain riming
  LOGICAL :: use_ice_graupel_conv_uli = .TRUE.

  ! Hallett-Mossop ice multiplication: ice/snow/graupel/hail - cloud/rain riming -> ice
  LOGICAL :: ice_multiplication = .TRUE.

  ! Graupel/hail - cloud/rain riming produces rain at temperatures above T_shed (<T_3)
  LOGICAL :: graupel_shedding = .FALSE.
  LOGICAL :: hail_shedding = .FALSE.

  ! Graupel/hail - cloud/rain riming produces rain at temperature above T_3
  LOGICAL :: enhanced_melting   = .TRUE.

  ! Allow freezing of cloud and rain drops
  LOGICAL :: drop_freeze = .TRUE.

  ! Cloud settings
  INTEGER :: nuc_c_typ = 0 ! Currently the only possible parameter
  !  0    Fixed CDNC as an input (model/CCN)
  INTEGER :: cloud_typ = 3 ! Seifert and Beheng (2000) (2-moment-scheme)
  !  0    No warm rain
  !  1-8  Different autoconversion and accretion schemes
  INTEGER :: ice_typ = 3 ! All ice species
  !  0    No ice
  !  >0   Ice, snow and graupel
  !  >1   Ice, snow, graupel and hail
  !  2    Ice/snow/graupel+rain=>graupel
  !  3    Ice/snow/graupel+rain=>graupel (smaller rain drop) or hail (larger rain drop)

  ! Fixed INP concentration (micro/nin_set) applied to supercooled liquid clouds
  REAL :: nin_set = 0.

  REAL, DIMENSION (:,:,:), ALLOCATABLE :: rrho_04
  REAL, DIMENSION (:,:,:), ALLOCATABLE :: rrho_c

!MODULE wolken_driver

  ! ... Gitter ...
  INTEGER          :: loc_ix, loc_iy, loc_iz

  ! ... Zeitschritt ....
  REAL :: dt

  ! ... Felder ...
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: &
       & q_cloud, q_ice, q_rain, q_graupel, q_snow, q_hail, &
       & n_cloud, n_ice, n_rain, n_graupel, n_snow, n_hail, &
       & q, p_0, t_0, rho_0, S_i

  ! ub>>
  ! Speicherfelder fuer die Depositionsraten von Eis und Schnee fuer die
  ! Entscheidung, ob Eis durch Riming Eis bleibt oder zu Graupel werden kann:
  REAL, DIMENSION(:,:,:), ALLOCATABLE ::   &
       deprate_ice, deprate_snow, &
       rimeqcrate_ice, rimencrate_ice, rimeqrrate_ice, rimeqirate_ice, rimenrrate_ice, &
       rimeqcrate_snow, rimencrate_snow, rimeqrrate_snow, rimeqirate_snow, rimenrrate_snow, &
       d_id_sp, d_sd_sp, d_rd_sp_ice, d_rd_sp_snow
  ! ub<<

  ! Statistics
  LOGICAL :: sflg=.FALSE.
  INTEGER :: out_mcrp_nout = 0
  real, save, allocatable :: out_mcrp_data(:,:,:,:)
  CHARACTER(LEN=:), SAVE, ALLOCATABLE :: out_mcrp_list(:)
  real, save, allocatable :: tmp_rv(:,:,:),tmp_rc(:,:,:),tmp_nr(:,:,:),tmp_rr(:,:,:), &
    tmp_ni(:,:,:),tmp_ri(:,:,:),tmp_rs(:,:,:),tmp_ns(:,:,:),tmp_rg(:,:,:),tmp_ng(:,:,:), &
    tmp_rh(:,:,:),tmp_nh(:,:,:)

  PUBLIC loc_ix, loc_iy, loc_iz, dt, p_0, T_0, rho_0, S_i, q, &
        q_cloud, q_ice, q_rain, q_snow, q_graupel, q_hail, &
        n_cloud, n_ice, n_rain, n_snow, n_graupel, n_hail, &
        init_seifert, alloc_driver, dealloc_driver, &
        alloc_wolken, dealloc_wolken, clouds, &
        particle, cloud, rain, ice, snow, graupel, hail, &
        sflg, out_mcrp_nout, out_mcrp_data, out_mcrp_list

CONTAINS


 !MODULE wolken_konstanten

  ! This subroutine has to be called once at the start of the model run by
  ! the main program. It properly sets the parameters for the different hydrometeor
  ! classes according to predefined parameter sets (see above).

  SUBROUTINE init_seifert( )
    IMPLICIT NONE
    INTEGER :: io, istat
    TYPE(PARTICLE) :: cldw ! Temporary cloud to match with level 4
    namelist /micro/ &
        drop_freeze, graupel_shedding, hail_shedding, enhanced_melting, ice_multiplication, &
        use_ice_graupel_conv_uli, &
        ice_typ, cloud_typ, nin_set, &
        q_krit_ic, D_krit_ic,  q_krit_ir, D_krit_ir, q_krit_sc, D_krit_sc, q_krit_sr, D_krit_sr, &
        q_krit_gc, D_krit_gc, q_krit_gr, q_krit_hc, D_krit_hc, q_krit_hr, q_krit_c, q_krit_r, &
        c_mult, &
        cldw,rain,ice,snow,graupel,hail

    ! Copy default cloud to cldw
    cldw = cloud

    ! read the optional microphysics namelist - overwrite the default parameters
    open(NEWUNIT=io,status='old',file='NAMELIST')
    read(io,nml=micro,iostat=istat)
    close(io)
    IF (istat>0) THEN
        ! namelist exists, but reading failed: find the problem by not using iostat
        open(NEWUNIT=io,status='old',file='NAMELIST')
        read(io,nml=micro)
        close(io)
    ENDIF

    ! Copy cloud back from cldw
    cloud = cldw

  END SUBROUTINE init_seifert

  REAL FUNCTION vent_coeff_a(parti,n)
    IMPLICIT NONE

    INTEGER        :: n
    TYPE(PARTICLE) :: parti

    vent_coeff_a = parti%a_ven * gfct((parti%nu+n+parti%b_geo)/parti%mu)              &
         &                  / gfct((parti%nu+1.0)/parti%mu)                        & 
         &                * ( gfct((parti%nu+1.0)/parti%mu)                        & 
         &                  / gfct((parti%nu+2.0)/parti%mu) )**(parti%b_geo+n-1.0) 

  END FUNCTION vent_coeff_a

  REAL FUNCTION vent_coeff_b(parti,n)
    IMPLICIT NONE

    INTEGER        :: n
    TYPE(PARTICLE) :: parti

    vent_coeff_b = parti%b_ven                                                  & 
         & * gfct((parti%nu+n+(m_f+1.0)*parti%b_geo+m_f*parti%b_vel)/parti%mu)  &
         &             / gfct((parti%nu+1.0)/parti%mu)                          & 
         &           * ( gfct((parti%nu+1.0)/parti%mu)                          &
         &             / gfct((parti%nu+2.0)/parti%mu)                          &
         &             )**((m_f+1.0)*parti%b_geo+m_f*parti%b_vel+n-1.0)

  END FUNCTION vent_coeff_b

  REAL FUNCTION moment_gamma(p,n) 
    IMPLICIT NONE

    INTEGER          :: n
    TYPE(PARTICLE)   :: p

    moment_gamma  = gfct((n+p%nu+1.0)/p%mu) / gfct((p%nu+1.0)/p%mu)        &
         &     * ( gfct((  p%nu+1.0)/p%mu) / gfct((p%nu+2.0)/p%mu) )**n
  END FUNCTION moment_gamma

  REAL FUNCTION lambda_gamma(p,x) 
    IMPLICIT NONE

    REAL :: x
    TYPE(PARTICLE)   :: p

    lambda_gamma  = ( gfct((p%nu+1.0)/p%mu) / gfct((p%nu+2.0)/p%mu) * x)**(-p%mu)
  END FUNCTION lambda_gamma

  REAL FUNCTION gfct(x)
    !*******************************************************************************
    !                                                                              *
    !       Gammafunktion aus Numerical Recipes (F77)                              *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    REAL cof(6)
    REAL stp,half,one,x,xx,fpf,tmp,ser,gamma
    INTEGER j

    DATA cof,stp/76.18009173e0,-86.50532033e0,24.01409822e0,  &
         &     -1.231739516e0,.120858003e-2,-.536382e-5,2.50662827465e0/
    DATA half,one,fpf/0.5e0,1.0e0,5.5e0/

    xx  = x  - one
    tmp = xx + fpf
    tmp = (xx + half) * LOG(tmp) - tmp
    ser = one
    DO j = 1,6
      xx  = xx  + one
      ser = ser + cof(j) / xx
    ENDDO
    gamma = tmp + LOG(stp*ser)
    gamma = EXP(gamma)

    gfct = gamma
    RETURN
  END FUNCTION gfct

!END MODULE wolken_konstanten

!MODULE gamma_functions_mp_seifert

  REAL FUNCTION gammln(x)
  !*******************************************************************************
  !                                                                              *
  !       LOG(Gamma function) taken from Press et al.,  Numerical Recipes (F77)
  !
  !       Log(Gammafunktion) aus Numerical Recipes (F77)                         *
  !       (original)                                                             *
  !*******************************************************************************
    IMPLICIT NONE
      
    REAL, INTENT(in) :: x
    REAL, SAVE :: cof(6), stp
    REAL :: xx,tmp,ser
    INTEGER :: j
    DATA cof /76.18009172947146e0,-86.50532032941677e0, &
         24.01409824083091e0,-1.231739572450155e0,.1208650973866179e-2, &
         -.5395239384953e-5/
    DATA stp /2.5066282746310005e0/

    xx  = x
    tmp = xx + 5.5e0
    tmp = (xx + 0.5e0) * LOG(tmp) - tmp
    ser = 1.000000000190015e0
    DO j = 1,6
       xx  = xx  + 1.0
       ser = ser + cof(j) / xx
    ENDDO
    gammln = tmp + LOG(stp*ser/x)
    RETURN
  END FUNCTION gammln
  
  !*******************************************************************************
  !
  !       Incomplete Gamma function taken from Press et al.,  Numerical Recipes (F77)
  !
  !       Unvollstaendige Gammafunktion aus Numerical Recipes (F77)              *
  !       (etwas umformuliert, aber dieselben Ergebnisse wie obige Originalfunktion)
  !*******************************************************************************

  !*******************************************************************************
  ! 1) diverse Hilfsfunktionen:

  SUBROUTINE gcf(gammcf,a,x,gln)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ITMAX = 100
    REAL, PARAMETER :: EPS = 3.e-7, FPMIN = 1.e-30
    REAL, INTENT(in) :: a, x
    REAL, INTENT(out) :: gammcf, gln

    INTEGER :: i
    REAL :: an,b,c,d,del,h

    gln=gammln(a)
    b=x+1.-a
    c=1./FPMIN
    d=1./b
    h=d
    DO i=1,ITMAX
      an=-i*(i-a)
      b=b+2.0
      d=an*d+b
      IF (ABS(d).LT.FPMIN) d=FPMIN
      c=b+an/c
      IF (ABS(c).LT.FPMIN) c=FPMIN
      d=1./d
      del=d*c
      h=h*del
      IF (ABS(del-1.).LT.EPS) EXIT
    END DO

    IF (ABS(del-1.).GE.EPS) THEN
      WRITE (*,*) 'ERROR in GCF: a too large, ITMAX too small (MODULE gamma_functions, src_seifert.f90)'
      gammcf = 0.0
      RETURN
    END IF

    gammcf=EXP(-x+a*LOG(x)-gln)*h

    RETURN
  END SUBROUTINE gcf

  SUBROUTINE gser(gamser,a,x,gln)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ITMAX = 100
    REAL, PARAMETER :: EPS=3.e-7
    REAL, INTENT(in) :: a, x
    REAL, INTENT(out) :: gamser, gln
      
    INTEGER :: n
    REAL :: ap,del,sum
    
    gln=gammln(a)
    IF (x.LE.0.) THEN
      IF (x.LT.0.) THEN
        WRITE (*,*) 'ERROR in GSER: x < 0 (MODULE gamma_functions, src_seifert.f90)'
      END IF
      gamser=0.0
      RETURN
    ENDIF

    ap=a
    sum=1./a
    del=sum
    DO n=1,ITMAX
      ap=ap+1.
      del=del*x/ap
      sum=sum+del
      IF (ABS(del).LT.ABS(sum)*EPS) EXIT
    END DO

    IF (ABS(del).GE.ABS(sum)*EPS) THEN
      WRITE (*,*) 'ERROR in GSER: a too large, ITMAX too small'
      WRITE (*,*) '  (MODULE gamma_functions, src_seifert.f90)'
      gamser = 0.0
      RETURN
    END IF

    gamser = sum*EXP(-x+a*LOG(x)-gln)

    RETURN
  END SUBROUTINE gser

  REAL FUNCTION gammp(a,x,gln)
  
    IMPLICIT NONE

    REAL, INTENT(in) :: a, x
    REAL, INTENT(out) :: gln
   
    REAL :: gammcf, gamser
    
    IF (x.LT.0.0 .OR. a .LE. 0.0) THEN
      WRITE(*,*) 'ERROR in GAMMP: bad arguments'
      WRITE(*,*) '  (MODULE gamma_functions, src_seifert.f90)'
      gammp = 0.0
      RETURN
    END IF

    IF (x .LT. a+1.)THEN
      CALL gser(gamser,a,x,gln)
      gammp = gamser
    ELSE
      CALL gcf(gammcf,a,x,gln)
      gammp = 1.0 - gammcf
    ENDIF
    RETURN
  END FUNCTION gammp


  !*******************************************************************************
  !
  ! Lower incomplete gamma function
  !
  ! Eigentliche untere unvollstaendige Gamma-Funktion, hier direkt
  ! das Integral
  !              int(0)(x) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

  REAL FUNCTION incgfct_lower(a,x)

    IMPLICIT NONE

    REAL, INTENT(in) :: a, x
    REAL :: gam, gln
    
    gam = gammp(a,x,gln)
    incgfct_lower = EXP(gln) * gam

  END FUNCTION incgfct_lower

!END MODULE gamma_functions_mp_seifert


!=======================================================================

!MODULE wolken_driver

  SUBROUTINE alloc_driver ()

    IMPLICIT NONE

    ALLOCATE( &
         &        q(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      p_0(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      T_0(0:loc_ix,1:loc_iy,1:loc_iz), &
         &    rho_0(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      S_i(0:loc_ix,1:loc_iy,1:loc_iz), &
         &  rrho_04(0:loc_ix,1:loc_iy,1:loc_iz), &
         &   rrho_c(0:loc_ix,1:loc_iy,1:loc_iz))
    q     = 0.0
    p_0   = 1e5
    T_0   = 3e2
    rho_0 = 1.2
    S_i   = 0.0

    ! ub>>
    IF (use_ice_graupel_conv_uli) THEN
      ALLOCATE( deprate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           deprate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqcrate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimencrate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqirate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqrrate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimenrrate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqcrate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimencrate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqirate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqrrate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimenrrate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           d_id_sp(0:loc_ix,1:loc_iy,1:loc_iz), &
           d_sd_sp(0:loc_ix,1:loc_iy,1:loc_iz), &
           d_rd_sp_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           d_rd_sp_snow(0:loc_ix,1:loc_iy,1:loc_iz))
      deprate_ice = 0.0
      deprate_snow = 0.0
      rimeqcrate_ice = 0.0
      rimencrate_ice = 0.0
      rimeqirate_ice = 0.0
      rimeqrrate_ice = 0.0
      rimenrrate_ice = 0.0
      rimeqcrate_snow = 0.0
      rimencrate_snow = 0.0
      rimeqirate_snow = 0.0
      rimeqrrate_snow = 0.0
      rimenrrate_snow = 0.0
      d_id_sp = 0.0
      d_sd_sp = 0.0
      d_rd_sp_ice = 0.0
      d_rd_sp_snow = 0.0
    END IF
    ! ub<<

    IF (sflg) THEN
        ALLOCATE( out_mcrp_data(0:loc_ix,1:loc_iy,1:loc_iz,out_mcrp_nout), &
            tmp_rv(0:loc_ix,1:loc_iy,1:loc_iz),tmp_rc(0:loc_ix,1:loc_iy,1:loc_iz), &
            tmp_nr(0:loc_ix,1:loc_iy,1:loc_iz),tmp_rr(0:loc_ix,1:loc_iy,1:loc_iz), &
            tmp_ni(0:loc_ix,1:loc_iy,1:loc_iz),tmp_ri(0:loc_ix,1:loc_iy,1:loc_iz), &
            tmp_rs(0:loc_ix,1:loc_iy,1:loc_iz),tmp_ns(0:loc_ix,1:loc_iy,1:loc_iz), &
            tmp_rg(0:loc_ix,1:loc_iy,1:loc_iz),tmp_ng(0:loc_ix,1:loc_iy,1:loc_iz), &
            tmp_rh(0:loc_ix,1:loc_iy,1:loc_iz),tmp_nh(0:loc_ix,1:loc_iy,1:loc_iz) )
        out_mcrp_data=0.
        tmp_rv=0.; tmp_rc=0.; tmp_nr=0.; tmp_rr=0.; tmp_ni=0.; tmp_ri=0.
        tmp_rs=0.; tmp_ns=0.; tmp_rg=0.; tmp_ng=0.; tmp_rh=0.; tmp_nh=0.
    ENDIF

  END SUBROUTINE alloc_driver

  SUBROUTINE dealloc_driver ()
    IMPLICIT NONE

    DEALLOCATE(q,p_0,T_0,rho_0,S_i,rrho_04,rrho_c)
    ! ub>>
    IF (use_ice_graupel_conv_uli) THEN
      DEALLOCATE(deprate_ice,deprate_snow,&
           rimeqcrate_ice,rimencrate_ice,&
           rimeqirate_ice,rimeqrrate_ice,&
           rimenrrate_ice,&
           rimeqcrate_snow,rimencrate_snow,&
           rimeqirate_snow,rimeqrrate_snow,&
           rimenrrate_snow, &
           d_id_sp, d_sd_sp, d_rd_sp_ice, d_rd_sp_snow)
    END IF
    ! ub<<
    IF (sflg) DEALLOCATE( out_mcrp_data,tmp_rv,tmp_rc,tmp_nr,tmp_rr,&
            tmp_ni,tmp_ri,tmp_rs,tmp_ns,tmp_rg,tmp_ng,tmp_rh,tmp_nh )
  END SUBROUTINE dealloc_driver

!END MODULE wolken_driver


!==============================================================================

!MODULE wolken_eis
  !*******************************************************************************
  !                                                                              *
  !     Diverse Unterprogramme zur Berechnung der Wolkenphysik mit Eisphase!     *
  !                                                                              *
  !*******************************************************************************

  SUBROUTINE ice_nucleation()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Nukleation der Wolkeneispartikel                        *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    REAL, PARAMETER :: eps  = 1.e-20
    REAL    :: nuc_q, ndiag
    INTEGER :: i,j,k

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          IF (s_i(i,j,k)>0.0 .AND. q_cloud(i,j,k)>0.001e-3 .AND. nin_set*ice%x_min>eps) THEN
            ! Cloud droplet freezing with fixed INP concentration
            ndiag = MAX(nin_set*rho_0(i,j,k) - (n_ice(i,j,k)+n_snow(i,j,k)+n_graupel(i,j,k)+n_hail(i,j,k)),0.0)
            nuc_q = MIN(ndiag*ice%x_min, q_cloud(i,j,k))

            q_ice(i,j,k) = q_ice(i,j,k) + nuc_q
            n_ice(i,j,k) = n_ice(i,j,k) + nuc_q/ice%x_min
            q_cloud(i,j,k) = q_cloud(i,j,k) - nuc_q
            n_cloud(i,j,k) = n_cloud(i,j,k) - nuc_q/cloud%x_max
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE ice_nucleation

  SUBROUTINE ice_cloud_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Eispartikel                               *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a             !..Absolute Temperatur
    REAL            :: q_i,n_i,x_i,d_i,v_i
    REAL            :: q_c,n_c,x_c,d_c,v_c,e_coll
    REAL            :: rime_n,rime_q
    REAL            :: conv_n,conv_q
    REAL            :: mult_n,mult_q,mult_1,mult_2
    REAL, SAVE      :: delta_n_ii,delta_n_ic,delta_n_cc
    REAL, SAVE      :: delta_q_ii,delta_q_ic,delta_q_cc
    REAL, SAVE      :: theta_n_ii,theta_n_ic,theta_n_cc
    REAL, SAVE      :: theta_q_ii,theta_q_ic,theta_q_cc
    REAL            :: const1,const3,const4,const5

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_ii = coll_delta_11(ice,cloud,0)
      delta_n_ic = coll_delta_12(ice,cloud,0)
      delta_n_cc = coll_delta_22(ice,cloud,0)
      delta_q_ii = coll_delta_11(ice,cloud,0) 
      delta_q_ic = coll_delta_12(ice,cloud,1)
      delta_q_cc = coll_delta_22(ice,cloud,1)

      theta_n_ii = coll_theta_11(ice,cloud,0)
      theta_n_ic = coll_theta_12(ice,cloud,0)
      theta_n_cc = coll_theta_22(ice,cloud,0)
      theta_q_ii = coll_theta_11(ice,cloud,0)
      theta_q_ic = coll_theta_12(ice,cloud,1)
      theta_q_cc = coll_theta_22(ice,cloud,1)

      firstcall = 1
    ENDIF

    const1 = e_ic/(D_coll_c - D_krit_c)
    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)
    const5 = alpha_spacefilling * rho_w/rho_ice
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          n_c = n_cloud(i,j,k)                                   !..Anzahldichte
          q_c = q_cloud(i,j,k)                                   !..Fluessigwassergehalt
          n_i = n_ice(i,j,k)                                     !..Anzahldichte
          q_i = q_ice(i,j,k)                                     !..Massendichte

          x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max)  !..mittlere Masse
          D_c = cloud%a_geo * x_c**cloud%b_geo                   !..mittlerer Durchmesser

          x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)      !..mittlere Masse
          D_i = ice%a_geo * x_i**ice%b_geo                       !..mittlerer Durchmesser

          T_a = T_0(i,j,k)

          IF (q_c > q_krit_c .AND. q_i > q_krit_ic .AND. D_i > D_krit_ic .AND. D_c > D_krit_c) THEN

            v_c = cloud%a_vel * x_c**cloud%b_vel * rrho_c(i,j,k)  !..mittlere Sedimentationsgeschw.
            v_i = ice%a_vel   * x_i**ice%b_vel   * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            e_coll = MIN(e_ic, MAX(const1*(D_c - D_krit_c), e_min))

            rime_n = pi/4.0 * e_coll * n_i * n_c * dt & 
                 &   *     (delta_n_ii * D_i*D_i + delta_n_ic * D_i*D_c + delta_n_cc * D_c*D_c) &
                 &   * SQRT(theta_n_ii * v_i*v_i - theta_n_ic * v_i*v_c + theta_n_cc * v_c*v_c  &
                 &         +ice_s_vel**2)

            rime_q = pi/4.0 * e_coll * n_i * q_c * dt & 
                 &   *     (delta_q_ii * D_i*D_i + delta_q_ic * D_i*D_c + delta_q_cc * D_c*D_c) &
                 &   * SQRT(theta_q_ii * v_i*v_i - theta_q_ic * v_i*v_c + theta_q_cc * v_c*v_c  &
                 &          +ice_s_vel**2)

            IF (.NOT.use_ice_graupel_conv_uli) THEN

              rime_q = MIN(q_c,rime_q)
              rime_n = MIN(n_c,rime_n)

              q_ice(i,j,k)   = q_ice(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4 
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
              ENDIF

              ! Umwandlung ice -> graupel

              IF (D_i > D_conv_ig) THEN
                conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*rho_ice*d_i**3/x_i - 1.0) )
                ! D_i darf durch CONV nicht kleiner werden als D_conv_ig
                !conv_q = MIN(q_ice(i,j,k)-n_i*(D_conv_ig/ice%a_geo)**(1.0/ice%b_geo),conv_q)
                conv_q = MIN(q_ice(i,j,k),conv_q)
                ! ub >>
                x_i = MIN(MAX((q_ice(i,j,k))/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse incl. Riming
                ! ub <<
                conv_n = conv_q / MAX(x_i,x_conv_g)
                conv_n = MIN(n_ice(i,j,k),conv_n)

                q_ice(i,j,k)     = q_ice(i,j,k)     - conv_q
                q_graupel(i,j,k) = q_graupel(i,j,k) + conv_q
                n_ice(i,j,k)     = n_ice(i,j,k)     - conv_n
                n_graupel(i,j,k) = n_graupel(i,j,k) + conv_n
              ENDIF

            ELSE

              rimeqcrate_ice(i,j,k) = rimeqcrate_ice(i,j,k) + rime_q
              rimencrate_ice(i,j,k) = rimencrate_ice(i,j,k) + rime_n

            END IF
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE ice_cloud_riming

  SUBROUTINE snow_cloud_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Schneepartikel                            *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a             !..Absolute Temperatur
    REAL            :: q_s,n_s,x_s,d_s,v_s
    REAL            :: q_c,n_c,x_c,d_c,v_c,e_coll
    REAL            :: rime_n,rime_q
    REAL            :: conv_n,conv_q
    REAL            :: mult_n,mult_q,mult_1,mult_2
    REAL, SAVE      :: delta_n_ss,delta_n_sc,delta_n_cc
    REAL, SAVE      :: delta_q_ss,delta_q_sc,delta_q_cc
    REAL, SAVE      :: theta_n_ss,theta_n_sc,theta_n_cc
    REAL, SAVE      :: theta_q_ss,theta_q_sc,theta_q_cc
    REAL            :: const1,const3,const4,const5

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_ss = coll_delta_11(snow,cloud,0)
      delta_n_sc = coll_delta_12(snow,cloud,0)
      delta_n_cc = coll_delta_22(snow,cloud,0)
      delta_q_ss = coll_delta_11(snow,cloud,0) 
      delta_q_sc = coll_delta_12(snow,cloud,1)
      delta_q_cc = coll_delta_22(snow,cloud,1)

      theta_n_ss = coll_theta_11(snow,cloud,0)
      theta_n_sc = coll_theta_12(snow,cloud,0)
      theta_n_cc = coll_theta_22(snow,cloud,0)
      theta_q_ss = coll_theta_11(snow,cloud,0)
      theta_q_sc = coll_theta_12(snow,cloud,1)
      theta_q_cc = coll_theta_22(snow,cloud,1)

      firstcall = 1
    ENDIF

! UB_20090316: re-invent constx 
    const1 = e_ic/(D_coll_c - D_krit_c)
    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)
    const5 = alpha_spacefilling * rho_w/rho_ice

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          n_c = n_cloud(i,j,k)                                   !..Anzahldichte
          q_c = q_cloud(i,j,k)                                   !..Fluessigwassergehalt
          n_s = n_snow(i,j,k)                                    !..Anzahldichte
          q_s = q_snow(i,j,k)                                    !..Massendichte

          x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max)  !..mittlere Masse
          D_c = cloud%a_geo * x_c**cloud%b_geo                   !..mittlerer Durchmesser

          x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse
          D_s = snow%a_geo * x_s**snow%b_geo                     !..mittlerer Durchmesser

          T_a = T_0(i,j,k)

          IF (q_c > q_krit_c .AND. q_s > q_krit_sc .AND. D_s > D_krit_sc .AND. D_c > D_krit_c) THEN

            v_c = cloud%a_vel * x_c**cloud%b_vel * rrho_c(i,j,k)  !..mittlere Sedimentationsgeschw.
            v_s = snow%a_vel  * x_s**snow%b_vel  * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            e_coll = MIN(e_sc, MAX(const1*(D_c - D_krit_c), e_min))

            rime_n = pi/4.0 * e_coll * n_s * n_c * dt & 
                 &   * (delta_n_ss * D_s**2 + delta_n_sc * D_s*D_c + delta_n_cc * D_c**2) &
                 &   * (theta_n_ss * v_s**2 - theta_n_sc * v_s*v_c + theta_n_cc * v_c**2  &
                 &     +snow_s_vel**2)**0.5

            rime_q = pi/4.0 * e_coll * n_s * q_c * dt & 
                 &   * (delta_q_ss * D_s**2 + delta_q_sc * D_s*D_c + delta_q_cc * D_c**2) &
                 &   * (theta_q_ss * v_s**2 - theta_q_sc * v_s*v_c + theta_q_cc * v_c**2  &
                 &     +snow_s_vel**2)**0.5

            IF (.NOT.use_ice_graupel_conv_uli) THEN

              rime_q = MIN(q_c,rime_q)
              rime_n = MIN(n_c,rime_n)

              q_snow(i,j,k)  = q_snow(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

              ! Umwandlung snow -> graupel

              IF (D_s > D_conv_sg) THEN
                conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*rho_ice*d_s**3/x_s - 1.0) )
                !conv_q = MIN(q_snow(i,j,k)-n_s*(D_conv_sg/snow%a_geo)**(1.0/snow%b_geo),conv_q)
                !conv_q = MAX(0.,conv_q)
                conv_q = MIN(q_snow(i,j,k),conv_q)
                ! ub >>
                x_s = MIN(MAX((q_snow(i,j,k))/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse incl. Riming
                ! ub <<
                conv_n = conv_q / MAX(x_s,x_conv_g)
                conv_n = MIN(n_snow(i,j,k),conv_n)

                q_snow(i,j,k)    = q_snow(i,j,k)    - conv_q
                q_graupel(i,j,k) = q_graupel(i,j,k) + conv_q
                n_snow(i,j,k)    = n_snow(i,j,k)    - conv_n
                n_graupel(i,j,k) = n_graupel(i,j,k) + conv_n
              ENDIF

            ELSE

              rimeqcrate_snow(i,j,k) = rimeqcrate_snow(i,j,k) + rime_q
              rimencrate_snow(i,j,k) = rimencrate_snow(i,j,k) + rime_n

            END IF

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE snow_cloud_riming

  SUBROUTINE graupel_cloud_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Graupelpartikel                           *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0
    REAL            :: T_a
    REAL            :: q_g,n_g,x_g,d_g,v_g
    REAL            :: q_c,n_c,x_c,d_c,v_c
    REAL            :: rime_n,rime_q,e_coll_n
    REAL            :: melt_n,melt_q,e_coll_q
    REAL            :: shed_n,shed_q,x_shed
    REAL            :: mult_n,mult_q,mult_1,mult_2
    REAL, SAVE      :: delta_n_gg,delta_n_gc,delta_n_cc
    REAL, SAVE      :: delta_q_gg,delta_q_gc,delta_q_cc
    REAL, SAVE      :: theta_n_gg,theta_n_gc,theta_n_cc
    REAL, SAVE      :: theta_q_gg,theta_q_gc,theta_q_cc
    REAL            :: const1,const2,const3,const4

    REAL, PARAMETER :: eps = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,cloud,0)
      delta_n_gc = coll_delta_12(graupel,cloud,0)
      delta_n_cc = coll_delta_22(graupel,cloud,0)
      delta_q_gg = coll_delta_11(graupel,cloud,0) 
      delta_q_gc = coll_delta_12(graupel,cloud,1)
      delta_q_cc = coll_delta_22(graupel,cloud,1)

      theta_n_gg = coll_theta_11(graupel,cloud,0)
      theta_n_gc = coll_theta_12(graupel,cloud,0)
      theta_n_cc = coll_theta_22(graupel,cloud,0)
      theta_q_gg = coll_theta_11(graupel,cloud,0)
      theta_q_gc = coll_theta_12(graupel,cloud,1)
      theta_q_cc = coll_theta_22(graupel,cloud,1)

      firstcall = 1
    ENDIF

    x_shed = 4./3.*pi * rho_w * r_shedding**3     !..Mittlere Masse der Sheddingtropfen

    const1 = e_gc/(D_coll_c - D_krit_c)    
    const2 = 1.0/(T_mult_opt - T_mult_min)
    const3 = 1.0/(T_mult_opt - T_mult_max)
    const4 = c_w / L_ew
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix             
          q_c = q_cloud(i,j,k)                                      !..Massendichte
          q_g = q_graupel(i,j,k)                                    !..Massendichte
          n_c = n_cloud(i,j,k)                                      !..Anzahldichte
          n_g = n_graupel(i,j,k)                                    !..Anzahldichte
          x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse
          D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
          x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max)     !..mittlere Masse
          D_c = cloud%a_geo * x_c**cloud%b_geo                      !..mittlerer Durchmesser
          T_a = T_0(i,j,k)                                          !..abs. Temperatur

          IF (q_c > q_krit_c .AND. q_g > q_krit_gc .AND. D_g > D_krit_gc .AND. D_c > D_krit_c) THEN

            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.
            v_c = cloud%a_vel   * x_c**cloud%b_vel   * rrho_c(i,j,k)  !..mittlere Sedimentationsgeschw.

            e_coll_n = MIN(e_gc, MAX(const1*(D_c - D_krit_c),e_min))
            e_coll_q = e_coll_n

            rime_n = pi/4.0 * e_coll_n * n_g * n_c * dt & 
                 &   * (delta_n_gg * D_g*D_g + delta_n_gc * D_g*D_c + delta_n_cc * D_c*D_c) &
                 &   * SQRT(theta_n_gg * v_g*v_g - theta_n_gc * v_g*v_c + theta_n_cc * v_c*v_c)

            rime_q = pi/4.0 * e_coll_q * n_g * q_c * dt & 
                 &   * (delta_q_gg * D_g*D_g + delta_q_gc * D_g*D_c + delta_q_cc * D_c*D_c) &
                 &   * SQRT(theta_q_gg * v_g*v_g - theta_q_gc * v_g*v_c + theta_q_cc * v_c*v_c)

            rime_q = MIN(q_c,rime_q)
            rime_n = MIN(n_c,rime_n)

            q_graupel(i,j,k) = q_graupel(i,j,k) + rime_q
            q_cloud(i,j,k)   = q_cloud(i,j,k)   - rime_q
            n_cloud(i,j,k)   = n_cloud(i,j,k)   - rime_n

            ! Eismultiplikation nach Hallet und Mossop

            IF (T_a < T_3 .AND. ice_multiplication) THEN
              mult_1 = const2*(T_a - T_mult_min) 
              mult_2 = const3*(T_a - T_mult_max) 
              mult_1 = MAX(0.,MIN(mult_1,1.))
              mult_2 = MAX(0.,MIN(mult_2,1.))
              mult_n = C_mult * mult_1 * mult_2 * rime_q
              mult_q = mult_n * ice%x_min
              mult_q = MIN(rime_q,mult_q)

              n_ice(i,j,k)     = n_ice(i,j,k)     + mult_n
              q_ice(i,j,k)     = q_ice(i,j,k)     + mult_q
              q_graupel(i,j,k) = q_graupel(i,j,k) - mult_q

            ENDIF

            ! enhancement of melting of graupel

            IF (T_a > T_3 .AND. enhanced_melting) THEN
              melt_q = const4 * (T_a - T_3) * rime_q
              melt_n = melt_q / x_g

              melt_q = MIN(q_graupel(i,j,k),melt_q)
              melt_n = MIN(n_graupel(i,j,k),melt_n)

              q_graupel(i,j,k) = q_graupel(i,j,k) - melt_q
              q_rain(i,j,k)    = q_rain(i,j,k)    + melt_q

              n_graupel(i,j,k) = n_graupel(i,j,k) - melt_n
              n_rain(i,j,k)    = n_rain(i,j,k)    + melt_n
            ELSE
              melt_q = 0.0
            ENDIF

            ! Shedding

            IF ((graupel_shedding .AND. D_g > D_shed_g .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
              q_g = q_graupel(i,j,k)
              n_g = n_graupel(i,j,k)
              x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse in SI

              shed_q = MIN(q_g,rime_q)
              shed_n = shed_q / MIN(x_shed,x_g)

              q_graupel(i,j,k) = q_graupel(i,j,k) - shed_q
              q_rain(i,j,k)    = q_rain(i,j,k)    + shed_q                   
              n_rain(i,j,k)    = n_rain(i,j,k)    + shed_n
            ENDIF

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_cloud_riming

  ! hn >>
  SUBROUTINE hail_cloud_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Hagelpartikel                             *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0
    REAL            :: T_a
    REAL            :: q_h,n_h,x_h,d_h,v_h
    REAL            :: q_c,n_c,x_c,d_c,v_c
    REAL            :: rime_n,rime_q,e_coll_n
    REAL            :: melt_n,melt_q,e_coll_q
    REAL            :: shed_n,shed_q,x_shed
    REAL            :: mult_n,mult_q,mult_1,mult_2
    REAL, SAVE      :: delta_n_hh,delta_n_hc,delta_n_cc
    REAL, SAVE      :: delta_q_hh,delta_q_hc,delta_q_cc
    REAL, SAVE      :: theta_n_hh,theta_n_hc,theta_n_cc
    REAL, SAVE      :: theta_q_hh,theta_q_hc,theta_q_cc
    REAL            :: const1,const2,const3,const4

    REAL, PARAMETER :: eps = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,cloud,0)
      delta_n_hc = coll_delta_12(hail,cloud,0)
      delta_n_cc = coll_delta_22(hail,cloud,0)
      delta_q_hh = coll_delta_11(hail,cloud,0) 
      delta_q_hc = coll_delta_12(hail,cloud,1)
      delta_q_cc = coll_delta_22(hail,cloud,1)

      theta_n_hh = coll_theta_11(hail,cloud,0)
      theta_n_hc = coll_theta_12(hail,cloud,0)
      theta_n_cc = coll_theta_22(hail,cloud,0)
      theta_q_hh = coll_theta_11(hail,cloud,0)
      theta_q_hc = coll_theta_12(hail,cloud,1)
      theta_q_cc = coll_theta_22(hail,cloud,1)

      firstcall = 1
    ENDIF

    x_shed = 4./3.*pi * rho_w * r_shedding**3     !..Mittlere Masse der Sheddingtropfen

    const1 = e_hc/(D_coll_c - D_krit_c)    
    const2 = 1.0/(T_mult_opt - T_mult_min)
    const3 = 1.0/(T_mult_opt - T_mult_max)
    const4 = c_w / L_ew
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix             
          q_c = q_cloud(i,j,k)                                      !..Massendichte
          q_h = q_hail(i,j,k)                                    !..Massendichte
          n_c = n_cloud(i,j,k)                                      !..Anzahldichte
          n_h = n_hail(i,j,k)                                    !..Anzahldichte
          x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max) !..mittlere Masse
          D_h = hail%a_geo * x_h**hail%b_geo                  !..mittlerer Durchmesser
          x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max)     !..mittlere Masse
          D_c = cloud%a_geo * x_c**cloud%b_geo                      !..mittlerer Durchmesser
          T_a = T_0(i,j,k)                                          !..abs. Temperatur

          IF (q_c > q_krit_c .AND. q_h > q_krit_hc .AND. D_h > D_krit_hc .AND. D_c > D_krit_c) THEN

            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.
            v_c = cloud%a_vel   * x_c**cloud%b_vel   * rrho_c(i,j,k)  !..mittlere Sedimentationsgeschw.

            e_coll_n = MIN(e_hc, MAX(const1*(D_c - D_krit_c),e_min))
            e_coll_q = e_coll_n

            rime_n = pi/4.0 * e_coll_n * n_h * n_c * dt & 
                 &   * (delta_n_hh * D_h*D_h + delta_n_hc * D_h*D_c + delta_n_cc * D_c*D_c) &
                 &   * SQRT(theta_n_hh * v_h*v_h - theta_n_hc * v_h*v_c + theta_n_cc * v_c*v_c)

            rime_q = pi/4.0 * e_coll_q * n_h * q_c * dt & 
                 &   * (delta_q_hh * D_h*D_h + delta_q_hc * D_h*D_c + delta_q_cc * D_c*D_c) &
                 &   * SQRT(theta_q_hh * v_h*v_h - theta_q_hc * v_h*v_c + theta_q_cc * v_c*v_c)

            rime_q = MIN(q_c,rime_q)
            rime_n = MIN(n_c,rime_n)

            q_hail(i,j,k)  = q_hail(i,j,k) + rime_q
            q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
            n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

            ! Eismultiplikation nach Hallet und Mossop

            IF (T_a < T_3 .AND. ice_multiplication) THEN
              mult_1 = const2*(T_a - T_mult_min) 
              mult_2 = const3*(T_a - T_mult_max) 
              mult_1 = MAX(0.,MIN(mult_1,1.))
              mult_2 = MAX(0.,MIN(mult_2,1.))
              mult_n = C_mult * mult_1 * mult_2 * rime_q
              mult_q = mult_n * ice%x_min
              mult_q = MIN(rime_q,mult_q)

              n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
              q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
              q_hail(i,j,k) = q_hail(i,j,k) - mult_q

            ENDIF

            ! enhancement of melting of hail

            IF (T_a > T_3 .AND. enhanced_melting) THEN
              melt_q = const4 * (T_a - T_3) * rime_q
              melt_n = melt_q / x_h

              melt_q = MIN(q_hail(i,j,k),melt_q)
              melt_n = MIN(n_hail(i,j,k),melt_n)

              q_hail(i,j,k) = q_hail(i,j,k) - melt_q
              q_rain(i,j,k) = q_rain(i,j,k) + melt_q

              n_hail(i,j,k) = n_hail(i,j,k) - melt_n
              n_rain(i,j,k) = n_rain(i,j,k) + melt_n
            ELSE
              melt_q = 0.0
            ENDIF

            ! Shedding

            IF ((D_h > D_shed_h .AND. T_a > T_shed .AND. hail_shedding) .OR. T_a > T_3 ) THEN
              q_h = q_hail(i,j,k)
              n_h = n_hail(i,j,k)
              x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max) !..mittlere Masse in SI

              shed_q = MIN(q_h,rime_q)
              shed_n = shed_q / MIN(x_shed,x_h)

              q_hail(i,j,k) = q_hail(i,j,k) - shed_q
              q_rain(i,j,k) = q_rain(i,j,k) + shed_q                   
              n_rain(i,j,k) = n_rain(i,j,k) + shed_n
            ENDIF

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE hail_cloud_riming

  REAL FUNCTION coll_delta(p1,n)
    IMPLICIT NONE
    TYPE(PARTICLE) :: p1
    INTEGER        :: n

    coll_delta = gfct((2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)         &
         &                  / gfct((p1%nu+1.0  )/p1%mu)         &
         &        * gfct((p1%nu+1.0)/p1%mu)**(2.0*p1%b_geo+n)   &
         &        / gfct((p1%nu+2.0)/p1%mu)**(2.0*p1%b_geo+n)

    RETURN
  END FUNCTION coll_delta

  REAL FUNCTION coll_delta_11(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_delta_11 = coll_delta(p1,n)

    RETURN
  END FUNCTION coll_delta_11

  REAL FUNCTION coll_delta_22(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_delta_22 = coll_delta(p2,n)

    RETURN
  END FUNCTION coll_delta_22

  REAL FUNCTION coll_delta_12(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_delta_12 = 2.0 * gfct((p1%b_geo+p1%nu+1.0)/p1%mu)               &
         &                       / gfct((p1%nu+1.0)/p1%mu)               &
         &                * gfct((p1%nu+1.0)/p1%mu)**(p1%b_geo)          &
         &                / gfct((p1%nu+2.0)/p1%mu)**(p1%b_geo)          &
         &              * gfct((p2%b_geo+p2%nu+1.0+n)/p2%mu)             &
         &                        /gfct((p2%nu+1.0  )/p2%mu)             &
         &                * gfct((p2%nu+1.0)/p2%mu)**(p2%b_geo+n)        &
         &                / gfct((p2%nu+2.0)/p2%mu)**(p2%b_geo+n)

    RETURN
  END FUNCTION coll_delta_12

  REAL FUNCTION coll_theta(p1,n)

    TYPE(PARTICLE) :: p1
    INTEGER        :: n

    coll_theta = gfct((2.0*p1%b_vel+2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)    &
         &               / gfct((2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)    &
         &          * gfct((p1%nu+1.0)/p1%mu)**(2.0*p1%b_vel)        &
         &          / gfct((p1%nu+2.0)/p1%mu)**(2.0*p1%b_vel)

    RETURN
  END FUNCTION coll_theta

  REAL FUNCTION coll_theta_11(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_theta_11 = coll_theta(p1,n)

    RETURN
  END FUNCTION coll_theta_11

  REAL FUNCTION coll_theta_22(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_theta_22 = coll_theta(p2,n)

    RETURN
  END FUNCTION coll_theta_22

  REAL FUNCTION coll_theta_12(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_theta_12 = 2.0 * gfct((p1%b_vel+2.0*p1%b_geo+p1%nu+1.0)/p1%mu)       &  
         &                    / gfct((2.0*p1%b_geo+p1%nu+1.0)/p1%mu)       &
         &              * gfct((p1%nu+1.0)/p1%mu)**(p1%b_vel)              &
         &              / gfct((p1%nu+2.0)/p1%mu)**(p1%b_vel)              &
         &           * gfct((p2%b_vel+2.0*p2%b_geo+p2%nu+1.0+n)/p2%mu)     &  
         &                    / gfct((2.0*p2%b_geo+p2%nu+1.0+n)/p2%mu)     &
         &              * gfct((p2%nu+1.0)/p2%mu)**(p2%b_vel)              &
         &              / gfct((p2%nu+2.0)/p2%mu)**(p2%b_vel)

    RETURN
  END FUNCTION coll_theta_12

  ! Faktor zur Berechnung des mittleren Durchmessers von verallg. gammaverteilten Hydrometeoren:
  FUNCTION D_average_factor (parti)

    ! gueltig fuer D = a_geo * x^b_geo
    ! Berechnung des mittleren Durchmessers: D_average = parti%b_geo * D_average_factor * (q/qn)**parti%b_geo
  
    IMPLICIT NONE
    
    REAL :: D_average_factor
    TYPE(PARTICLE), INTENT(in) :: parti
    
    D_average_factor = &
         ( gfct( (parti%b_geo+parti%nu+1.0)/parti%mu ) / &
         gfct( (parti%nu+1.0)/parti%mu ) ) * &
         (gfct( (parti%nu+1.0)/parti%mu ) / gfct( (parti%nu+2.0)/parti%mu ) ) ** parti%b_geo
    
  END FUNCTION D_average_factor

  SUBROUTINE ice_rain_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Eispartikel                               *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a             !..Absolute Temperatur
    REAL            :: q_i,n_i,x_i,d_i,v_i,d_id
    REAL            :: q_r,n_r,x_r,d_r,v_r,d_rd
    REAL            :: rime_n,rime_qi,rime_qr
    REAL            :: mult_n,mult_q,mult_1,mult_2
    REAL, SAVE      :: delta_n_ii,delta_n_ir,           delta_n_rr
    REAL, SAVE      :: delta_q_ii,delta_q_ir,delta_q_ri,delta_q_rr
    REAL, SAVE      :: theta_n_ii,theta_n_ir,           theta_n_rr
    REAL, SAVE      :: theta_q_ii,theta_q_ir,theta_q_ri,theta_q_rr,D_av_fakt_i,D_av_fakt_r
    REAL            :: const3,const4

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_ii = coll_delta_11(ice,rain,0)
      delta_n_ir = coll_delta_12(ice,rain,0)
      delta_n_rr = coll_delta_22(ice,rain,0)
      delta_q_ii = coll_delta_11(ice,rain,1) 
      delta_q_ir = coll_delta_12(ice,rain,1)
      delta_q_ri = coll_delta_12(rain,ice,1)
      delta_q_rr = coll_delta_22(ice,rain,1)

      theta_n_ii = coll_theta_11(ice,rain,0)
      theta_n_ir = coll_theta_12(ice,rain,0)
      theta_n_rr = coll_theta_22(ice,rain,0)
      theta_q_ii = coll_theta_11(ice,rain,1)
      theta_q_ir = coll_theta_12(ice,rain,1)
      theta_q_ri = coll_theta_12(rain,ice,1)
      theta_q_rr = coll_theta_22(ice,rain,1)

      D_av_fakt_r = D_average_factor(rain)
      D_av_fakt_i = D_average_factor(ice)

      firstcall = 1
    ENDIF

    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)                              !..Fluessigwassergehalt in SI
          q_i = q_ice(i,j,k)                               !..Fluessigwassergehalt in SI
          n_r = n_rain(i,j,k)                              !..Anzahldichte in SI
          n_i = n_ice(i,j,k)                               !..Anzahldichte in SI

          T_a = T_0(i,j,k)

          x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse in SI     
          D_i = ice%a_geo * x_i**ice%b_geo                     !..mittlerer Durchmesser
          D_id = D_av_fakt_i * D_i 

          IF (q_r > q_krit_r .AND. q_i > q_krit_ir .AND. D_i > D_krit_ir) THEN

            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)  !..mittlere Masse in SI     

            v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)    !..mittlere Sedimentationsgeschw.

            D_r = rain%a_geo * x_r**rain%b_geo                   !..mittlerer Durchmesser
            v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.
            D_rd = D_av_fakt_r * D_r

            rime_n  = pi/4.0 * n_i * n_r * dt & 
                 &   * (delta_n_ii * D_i*D_i + delta_n_ir * D_i*D_r + delta_n_rr * D_r*D_r) &
                 &   * SQRT(theta_n_ii * v_i*v_i - theta_n_ir * v_i*v_r + theta_n_rr * v_r*v_r  &
                 &     +ice_s_vel**2)

            rime_qr = pi/4.0 * n_i * q_r * dt & 
                 &   * (delta_n_ii * D_i*D_i + delta_q_ir * D_i*D_r + delta_q_rr * D_r*D_r) &
                 &   * SQRT(theta_n_ii * v_i*v_i - theta_q_ir * v_i*v_r + theta_q_rr * v_r*v_r  &
                 &     +ice_s_vel**2)

            rime_qi = pi/4.0 * n_r * q_i * dt & 
                 &   * (delta_q_ii * D_i*D_i + delta_q_ri * D_i*D_r + delta_n_rr * D_r*D_r) &
                 &   * SQRT(theta_q_ii * v_i*v_i - theta_q_ri * v_i*v_r + theta_n_rr * v_r*v_r  &
                 &     +ice_s_vel**2)

            IF (.NOT.use_ice_graupel_conv_uli) THEN

              rime_n  = MIN(MIN(n_r,n_i),rime_n)
              rime_qr = MIN(q_r,rime_qr)
              rime_qi = MIN(q_i,rime_qi)

              n_ice(i,j,k)  = n_ice(i,j,k)  - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_ice(i,j,k)  = q_ice(i,j,k)  - rime_qi
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr


              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              mult_n = 0.0  ! UB_20081120
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)
              ENDIF

              IF (T_a >= T_3) THEN
                IF (D_id > D_rd) THEN
                  ! Regen wird wieder abgeshedded
                  n_ice(i,j,k) = n_ice(i,j,k) + rime_n
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_qr / x_r
                  q_ice(i,j,k) = q_ice(i,j,k) + rime_qi     ! UB_20081120 - mult_q
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr
                ELSE
                  ! Eisteilchen schmelzen instantan zu Regen.
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_n
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr + rime_qi   ! UB_20081120 - mult_q
                END IF
              ELSE
                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n  ! UB_20081120
                q_ice(i,j,k) = q_ice(i,j,k)  + mult_q  ! UB_20081120
                IF (ice_typ < 3) THEN
                  ! Eis + angefrorenes Regenwasser ergibt Graupel:
                  n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                  q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qi + rime_qr - mult_q
                ELSE
                  IF (D_id > D_rd) THEN
                    ! Eis + angefrorenes Regenwasser ergibt Graupel:
                    n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                    q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qi + rime_qr - mult_q
                  ELSE
                    ! Eis + angefrorenes Regenwasser ergibt Hagel:
                    n_hail(i,j,k) = n_hail(i,j,k) + rime_n
                    q_hail(i,j,k) = q_hail(i,j,k) + rime_qi + rime_qr - mult_q
                  END IF
                END IF
              END IF

            ELSE

              rimenrrate_ice(i,j,k) = rimenrrate_ice(i,j,k) + rime_n
              rimeqirate_ice(i,j,k) = rimeqirate_ice(i,j,k) + rime_qi
              rimeqrrate_ice(i,j,k) = rimeqrrate_ice(i,j,k) + rime_qr
              d_id_sp(i,j,k)        = D_id
              d_rd_sp_ice(i,j,k)    = D_rd

            END IF

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE ice_rain_riming

  SUBROUTINE snow_rain_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Schneepartikel                            *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a             !..Absolute Temperatur
    REAL            :: q_s,n_s,x_s,d_s,v_s,d_sd
    REAL            :: q_r,n_r,x_r,d_r,v_r,d_rd
    REAL            :: rime_n,rime_qr,rime_qs
    REAL            :: mult_n,mult_q,mult_1,mult_2
    REAL, SAVE      :: delta_n_ss,delta_n_sr,delta_n_rr
    REAL, SAVE      :: delta_q_ss,delta_q_sr,delta_q_rs,delta_q_rr
    REAL, SAVE      :: theta_n_ss,theta_n_sr,theta_n_rr
    REAL, SAVE      :: theta_q_ss,theta_q_sr,theta_q_rs,theta_q_rr,D_av_fakt_s,D_av_fakt_r
    REAL            :: const3,const4

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_ss = coll_delta_11(snow,rain,0)
      delta_n_sr = coll_delta_12(snow,rain,0)
      delta_n_rr = coll_delta_22(snow,rain,0)
      delta_q_ss = coll_delta_11(snow,rain,1) 
      delta_q_sr = coll_delta_12(snow,rain,1)
      delta_q_rs = coll_delta_12(rain,snow,1)
      delta_q_rr = coll_delta_22(snow,rain,1)

      theta_n_ss = coll_theta_11(snow,rain,0)
      theta_n_sr = coll_theta_12(snow,rain,0)
      theta_n_rr = coll_theta_22(snow,rain,0)
      theta_q_ss = coll_theta_11(snow,rain,1)
      theta_q_sr = coll_theta_12(snow,rain,1)
      theta_q_rs = coll_theta_12(rain,snow,1)
      theta_q_rr = coll_theta_22(snow,rain,1)

      D_av_fakt_r = D_average_factor(rain)
      D_av_fakt_s = D_average_factor(snow)

      firstcall = 1
    ENDIF

    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)                              !..Fluessigwassergehalt in SI
          q_s = q_snow(i,j,k)                              !..Fluessigwassergehalt in SI
          n_r = n_rain(i,j,k)                              !..Anzahldichte in SI
          n_s = n_snow(i,j,k)                              !..Anzahldichte in SI

          x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max) !..mittlere Masse in SI     
          D_s = snow%a_geo * x_s**snow%b_geo                  !..mittlerer Durchmesser
          D_sd = D_av_fakt_s * D_s

          T_a = T_0(i,j,k)

          IF (q_r > q_krit_r .AND. q_s > q_krit_sr .AND. D_s > D_krit_sr) THEN

            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)  !..mittlere Masse in SI     
            D_r = rain%a_geo * x_r**rain%b_geo                   !..mittlerer Durchmesser
            v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.
            D_rd = D_av_fakt_r * D_r

            rime_n = pi/4.0 * n_s * n_r * dt & 
                 &   * (delta_n_ss * D_s**2 + delta_n_sr * D_s*D_r + delta_n_rr * D_r**2) &
                 &   * (theta_n_ss * v_s**2 - theta_n_sr * v_s*v_r + theta_n_rr * v_r**2  &
                 &     +snow_s_vel**2)**0.5

            rime_qr = pi/4.0 * n_s * q_r * dt & 
                 &   * (delta_n_ss * D_s**2 + delta_q_sr * D_s*D_r + delta_q_rr * D_r**2) &
                 &   * (theta_n_ss * v_s**2 - theta_q_sr * v_s*v_r + theta_q_rr * v_r**2  &
                 &     +snow_s_vel**2)**0.5

            rime_qs = pi/4.0 * n_r * q_s * dt & 
                 &   * (delta_q_ss * D_s**2 + delta_q_rs * D_s*D_r + delta_n_rr * D_r**2) &
                 &   * (theta_q_ss * v_s**2 - theta_q_rs * v_s*v_r + theta_n_rr * v_r**2  &
                 &     +snow_s_vel**2)**0.5

            IF (.NOT.use_ice_graupel_conv_uli) THEN

              rime_qr = MIN(q_r,rime_qr)
              rime_qs = MIN(q_s,rime_qs)
              rime_n  = MIN(n_r,rime_n)
              rime_n  = MIN(n_s,rime_n)

              n_snow(i,j,k) = n_snow(i,j,k) - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_snow(i,j,k) = q_snow(i,j,k) - rime_qs
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr


              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              mult_n = 0.0 ! UB_20081120
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)

              ENDIF

              IF (T_a >= T_3) THEN
                IF (D_sd > D_rd) THEN
                  ! Regen wird wieder abgeshedded
                  n_snow(i,j,k) = n_snow(i,j,k) + rime_n
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_qr / x_r
                  q_snow(i,j,k) = q_snow(i,j,k) + rime_qs ! UB_20081120 - mult_q
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr
                ELSE
                  ! Schneeflocken werden instantan zu Regen.
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_n
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr + rime_qs ! UB_20081120 - mult_q
                END IF
              ELSE
                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n ! UB_20081120
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q ! UB_20081120
                IF (ice_typ < 3) THEN
                  ! Schnee + angefrorenes Regenwasser ergibt Graupel:
                  n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                  q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qr + rime_qs - mult_q
                ELSE
                  IF (D_sd > D_rd) THEN
                    ! Schnee + angefrorenes Regenwasser ergibt Graupel:
                    n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                    q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qr + rime_qs - mult_q
                  ELSE
                    ! Schnee + angefrorenes Regenwasser ergibt Hagel:
                    n_hail(i,j,k) = n_hail(i,j,k) + rime_n
                    q_hail(i,j,k) = q_hail(i,j,k) + rime_qr + rime_qs - mult_q
                  END IF
                END IF
              END IF


            ELSE

              rimenrrate_snow(i,j,k)  = rimenrrate_snow(i,j,k)  + rime_n
              rimeqirate_snow(i,j,k)  = rimeqirate_snow(i,j,k) + rime_qs
              rimeqrrate_snow(i,j,k)  = rimeqrrate_snow(i,j,k) + rime_qr
              d_sd_sp(i,j,k) = D_sd
              d_rd_sp_snow(i,j,k) = D_rd

            END IF

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE snow_rain_riming

  SUBROUTINE graupel_rain_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Graupelpartikel                           *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a
    REAL            :: q_g,n_g,x_g,d_g,v_g,d_h,x_h,d_gd
    REAL            :: q_r,n_r,x_r,d_r,v_r,x_shed,q_h,n_h,d_rd
    REAL            :: rime_n,rime_q,rime_qg,rime_qr,melt_n,melt_q,shed_n,shed_q
    REAL            :: mult_n,mult_q,mult_1,mult_2
    REAL, SAVE      :: delta_n_gg,delta_n_gr,delta_n_rr
    REAL, SAVE      :: delta_q_gg,delta_q_gr,delta_q_rg,delta_q_rr
    REAL, SAVE      :: theta_n_gg,theta_n_gr,theta_n_rr
    REAL, SAVE      :: theta_q_gg,theta_q_gr,theta_q_rg,theta_q_rr,D_av_fakt_g,D_av_fakt_r
    REAL            :: const3,const4

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,rain,0)
      delta_n_gr = coll_delta_12(graupel,rain,0)
      delta_n_rr = coll_delta_22(graupel,rain,0)
!      delta_q_gg = coll_delta_11(graupel,rain,0) ! UB_20081120
      delta_q_gg = coll_delta_11(graupel,rain,1)  ! UB_20081120
      delta_q_gr = coll_delta_12(graupel,rain,1)
      delta_q_rg = coll_delta_12(rain,graupel,1)
      delta_q_rr = coll_delta_22(graupel,rain,1)

      theta_n_gg = coll_theta_11(graupel,rain,0)
      theta_n_gr = coll_theta_12(graupel,rain,0)
      theta_n_rr = coll_theta_22(graupel,rain,0)
!      theta_q_gg = coll_theta_11(graupel,rain,0) ! UB_20081120
      theta_q_gg = coll_theta_11(graupel,rain,1)  ! UB_20081120
      theta_q_gr = coll_theta_12(graupel,rain,1)
      theta_q_rg = coll_theta_12(rain,graupel,1)
      theta_q_rr = coll_theta_22(graupel,rain,1)

      D_av_fakt_r = D_average_factor(rain)
      D_av_fakt_g = D_average_factor(graupel)

      firstcall = 1
    ENDIF

    x_shed = 4./3. * pi * rho_w * r_shedding**3

    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)                                    !..Fluessigwassergehalt in SI
          q_g = q_graupel(i,j,k)                                 !..Fluessigwassergehalt in SI

          T_a = T_0(i,j,k)
          IF (q_r > q_krit_r .AND. q_g > q_krit_gr) THEN

            n_r = n_rain(i,j,k)                                 !..Anzahldichte in SI
            n_g = n_graupel(i,j,k)                              !..Anzahldichte in SI

            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)       !..mittlere Masse in SI     
            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse in SI     

            D_r = rain%a_geo * x_r**rain%b_geo                        !..mittlerer Durchmesser
            v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)       !..mittlere Sedimentationsgeschw.
            D_rd = D_av_fakt_r * D_r

            D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.
            D_gd = D_av_fakt_g * D_g 

            IF (T_a >= T_3 .OR. ice_typ < 3 .OR. D_gd > D_rd) THEN

              rime_n = pi/4.0 * n_g * n_r * dt & 
                   &   * (delta_n_gg * D_g**2 + delta_n_gr * D_g*D_r + delta_n_rr * D_r**2) &
                   &   * (theta_n_gg * v_g**2 - theta_n_gr * v_g*v_r + theta_n_rr * v_r**2)**0.5

              rime_q = pi/4.0 * n_g * q_r * dt & 
                   &   * (delta_n_gg * D_g**2 + delta_q_gr * D_g*D_r + delta_q_rr * D_r**2) &
                   &   * (theta_n_gg * v_g**2 - theta_q_gr * v_g*v_r + theta_q_rr * v_r**2)**0.5

              rime_n = MIN(n_r,rime_n)
              rime_q = MIN(q_r,rime_q)

              q_graupel(i,j,k) = q_graupel(i,j,k) + rime_q
              q_rain(i,j,k)    = q_rain(i,j,k)    - rime_q
              n_rain(i,j,k)    = n_rain(i,j,k)    - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)     = n_ice(i,j,k)     + mult_n
                q_ice(i,j,k)     = q_ice(i,j,k)     + mult_q
                q_graupel(i,j,k) = q_graupel(i,j,k) - mult_q
              ENDIF

              ! enhancement of melting of graupel

              IF (T_a > T_3 .AND. enhanced_melting) THEN
                melt_q = c_w / L_ew * (T_a - T_3) * rime_q
                melt_n = melt_q / x_g

                melt_q = MIN(q_graupel(i,j,k),melt_q)
                melt_n = MIN(n_graupel(i,j,k),melt_n)

                q_graupel(i,j,k) = q_graupel(i,j,k) - melt_q
                q_rain(i,j,k)    = q_rain(i,j,k)    + melt_q

                n_graupel(i,j,k) = n_graupel(i,j,k) - melt_n
                n_rain(i,j,k)    = n_rain(i,j,k)    + melt_n
              ELSE
                melt_q = 0.0
              ENDIF

              ! Shedding

              IF ((graupel_shedding .AND. D_g > D_shed_g .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
                q_g = q_graupel(i,j,k)
                n_g = n_graupel(i,j,k)
                x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse in SI     

                shed_q = MIN(q_g,rime_q)
                IF (T_a <= T_3) THEN
                  shed_n = shed_q / MIN(x_shed,x_g)
                ELSE
                  shed_n = shed_q / MAX(x_r,x_g)
                ENDIF

                q_graupel(i,j,k) = q_graupel(i,j,k) - shed_q
                q_rain(i,j,k)    = q_rain(i,j,k)    + shed_q                   
                n_rain(i,j,k)    = n_rain(i,j,k)    + shed_n
              ENDIF

            ELSE   ! T_a < T_3 .and. ice_typ >= 3 .and. D_g <= D_r

              rime_n = pi/4.0 * n_g * n_r * dt & 
                   &   * (delta_n_gg * D_g**2 + delta_n_gr * D_g*D_r + delta_n_rr * D_r**2) &
                   &   * (theta_n_gg * v_g**2 - theta_n_gr * v_g*v_r + theta_n_rr * v_r**2)**0.5

              rime_qr = pi/4.0 * n_g * q_r * dt & 
                   &   * (delta_n_gg * D_g**2 + delta_q_gr * D_g*D_r + delta_q_rr * D_r**2) &
                   &   * (theta_n_gg * v_g**2 - theta_q_gr * v_g*v_r + theta_q_rr * v_r**2)**0.5

              rime_qg = pi/4.0 * n_r * q_g * dt & 
                   &   * (delta_q_gg * D_g**2 + delta_q_rg * D_g*D_r + delta_n_rr * D_r**2) &
                   &   * (theta_q_gg * v_g**2 - theta_q_rg * v_g*v_r + theta_n_rr * v_r**2)**0.5

              rime_n = MIN(MIN(n_r,n_g),rime_n)
              rime_qr = MIN(q_r,rime_qr)
              rime_qg = MIN(q_g,rime_qg)

              n_graupel(i,j,k)  = n_graupel(i,j,k)  - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_graupel(i,j,k)  = q_graupel(i,j,k)  - rime_qg
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr

              n_hail(i,j,k) = n_hail(i,j,k) + rime_n
              q_hail(i,j,k) = q_hail(i,j,k) + rime_qg + rime_qr

              ! Eismultiplikation nach Hallet und Mossop

              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)

                n_ice(i,j,k)     = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)     = q_ice(i,j,k)  + mult_q
                q_hail(i,j,k)    = q_hail(i,j,k) - mult_q
              ENDIF

              ! Shedding

              x_h = MIN(MAX(q_hail(i,j,k)/(n_hail(i,j,k)+eps),hail%x_min),hail%x_max)       !..mittlere Masse in SI     
              D_h = hail%a_geo * x_h**hail%b_geo                        !..mittlerer Durchmesser

              IF (hail_shedding .AND. D_h > D_shed_h .AND. T_a > T_shed) THEN
                q_h = q_hail(i,j,k)
                n_h = n_hail(i,j,k)

                ! Vorher wurde Graupel + Regen zu Hagel umgewandelt und n_h und q_h haben sich erhoeht.
                ! Durch Shedding soll der Hagel nun die durch den Regen hinzugekommene Masse wieder verlieren.
                ! graupel_rain_riming wirkt also hier als Katalysator zur Umwandlung von Graupel in Hagel und
                ! zur Umwandlung von Regen in kleine Regentropfen.
                ! Zu Ueberlegen waere, ob nicht in diesem Falle der Graupel Graupel bleiben soll!!!
                shed_q = MIN(q_h,rime_qr)
                shed_n = shed_q / MIN(x_shed,x_h)

                q_hail(i,j,k)    = q_hail(i,j,k)    - shed_q
                q_rain(i,j,k)    = q_rain(i,j,k)    + shed_q                   
                n_rain(i,j,k)    = n_rain(i,j,k)    + shed_n
              ENDIF

            ENDIF  ! T_a >= T_3 .or. ice_typ < 3
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_rain_riming

  ! hn >>
  SUBROUTINE hail_rain_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Hagelpartikel                             *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a
    REAL            :: q_h,n_h,x_h,d_h,v_h
    REAL            :: q_r,n_r,x_r,d_r,v_r,x_shed
    REAL            :: rime_n,rime_q,melt_n,melt_q,shed_n,shed_q
    REAL            :: mult_n,mult_q,mult_1,mult_2
    REAL, SAVE      :: delta_n_hh,delta_n_hr,delta_n_rr
    REAL, SAVE      :: delta_q_hh,delta_q_hr,delta_q_rr
    REAL, SAVE      :: theta_n_hh,theta_n_hr,theta_n_rr
    REAL, SAVE      :: theta_q_hh,theta_q_hr,theta_q_rr
    REAL            :: const3,const4

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,rain,0)
      delta_n_hr = coll_delta_12(hail,rain,0)
      delta_n_rr = coll_delta_22(hail,rain,0)
      delta_q_hh = coll_delta_11(hail,rain,0) 
      delta_q_hr = coll_delta_12(hail,rain,1)
      delta_q_rr = coll_delta_22(hail,rain,1)

      theta_n_hh = coll_theta_11(hail,rain,0)
      theta_n_hr = coll_theta_12(hail,rain,0)
      theta_n_rr = coll_theta_22(hail,rain,0)
      theta_q_hh = coll_theta_11(hail,rain,0)
      theta_q_hr = coll_theta_12(hail,rain,1)
      theta_q_rr = coll_theta_22(hail,rain,1)

      firstcall = 1
    ENDIF

    x_shed = 4./3. * pi * rho_w * r_shedding**3

    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)                                    !..Fluessigwassergehalt in SI
          q_h = q_hail(i,j,k)                                    !..Fluessigwassergehalt in SI

          T_a = T_0(i,j,k)
          IF (q_r > q_krit_r .AND. q_h > q_krit_hr) THEN
            n_r = n_rain(i,j,k)                                 !..Anzahldichte in SI
            n_h = n_hail(i,j,k)                                 !..Anzahldichte in SI

            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)       !..mittlere Masse in SI     
            x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max)       !..mittlere Masse in SI     

            D_r = rain%a_geo * x_r**rain%b_geo                        !..mittlerer Durchmesser
            v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)       !..mittlere Sedimentationsgeschw.

            D_h = hail%a_geo * x_h**hail%b_geo                  !..mittlerer Durchmesser
            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.


            rime_n = pi/4.0 * n_h * n_r * dt & 
                 &   * (delta_n_hh * D_h**2 + delta_n_hr * D_h*D_r + delta_n_rr * D_r**2) &
                 &   * (theta_n_hh * v_h**2 - theta_n_hr * v_h*v_r + theta_n_rr * v_r**2)**0.5

            rime_q = pi/4.0 * n_h * q_r * dt & 
                 &   * (delta_q_hh * D_h**2 + delta_q_hr * D_h*D_r + delta_q_rr * D_r**2) &
                 &   * (theta_q_hh * v_h**2 - theta_q_hr * v_h*v_r + theta_q_rr * v_r**2)**0.5

            rime_n = MIN(n_r,rime_n)
            rime_q = MIN(q_r,rime_q)

            q_hail(i,j,k) = q_hail(i,j,k) + rime_q
            q_rain(i,j,k) = q_rain(i,j,k) - rime_q
            n_rain(i,j,k) = n_rain(i,j,k) - rime_n

            ! Eismultiplikation nach Hallet und Mossop

            IF (T_a < T_3 .AND. ice_multiplication) THEN
              mult_1 = (T_a - T_mult_min) * const3
              mult_2 = (T_a - T_mult_max) * const4
              mult_1 = MAX(0.,MIN(mult_1,1.))
              mult_2 = MAX(0.,MIN(mult_2,1.))
              mult_n = C_mult * mult_1 * mult_2 * rime_q
              mult_q = mult_n * ice%x_min
              mult_q = MIN(rime_q,mult_q)

              n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
              q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
              q_hail(i,j,k) = q_hail(i,j,k) - mult_q
            ENDIF

            ! enhancement of melting of hail

            IF (T_a > T_3 .AND. enhanced_melting) THEN
              melt_q = c_w / L_ew * (T_a - T_3) * rime_q
              melt_n = melt_q / x_h

              melt_q = MIN(q_hail(i,j,k),melt_q)
              melt_n = MIN(n_hail(i,j,k),melt_n)

              q_hail(i,j,k) = q_hail(i,j,k) - melt_q
              q_rain(i,j,k) = q_rain(i,j,k) + melt_q

              n_hail(i,j,k) = n_hail(i,j,k) - melt_n
              n_rain(i,j,k) = n_rain(i,j,k) + melt_n
            ELSE
              melt_q = 0.0
            ENDIF

            ! Shedding

            IF ((hail_shedding .AND. D_h > D_shed_h .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
              q_h = q_hail(i,j,k)
              n_h = n_hail(i,j,k)
              x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max) !..mittlere Masse in SI     

              shed_q = MIN(q_h,rime_q)
              IF (T_a <= T_3) THEN
                shed_n = shed_q / MIN(x_shed,x_h)
              ELSE
                !                shed_n = shed_q / x_shed
                shed_n = shed_q / MAX(x_r,x_h)
              ENDIF

              q_hail(i,j,k) = q_hail(i,j,k) - shed_q
              q_rain(i,j,k) = q_rain(i,j,k)    + shed_q                   
              n_rain(i,j,k) = n_rain(i,j,k)    + shed_n
            ENDIF

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE hail_rain_riming
  ! << hn

  ! ub>>
  SUBROUTINE complete_ice_snow_riming ()
    !*******************************************************************************
    !                                                                              *
    !       Im Falle von use_ice_graupel_conv_uli = .true. :                       *
    !                                                                              *
    !       Summierung der vorher lediglich gespeicherten riming-Raten             *
    !       auf die entsprechenden Felder (ice, snow, cloud, rain)                 *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen
    INTEGER                     :: i,j,k

    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a             !..Absolute Temperatur
    REAL            :: n_i,x_i,d_i,D_id
    REAL            :: n_s,x_s,d_s,D_sd
    REAL            :: x_r
    REAL            :: D_rd
    REAL            :: rime_n,rime_q, rime_qr, rime_qs, rime_qi
    REAL            :: conv_n,conv_q
    REAL            :: mult_n,mult_q,mult_1,mult_2
    REAL            :: const3,const4,const5
    REAL, SAVE      :: D_av_fakt_i,D_av_fakt_s

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall .NE. 1) THEN
      D_av_fakt_i = D_average_factor(ice)
      D_av_fakt_s = D_average_factor(snow)
      firstcall = 1
    END IF

    const3 = 1.0/(T_mult_opt - T_mult_min)
    const4 = 1.0/(T_mult_opt - T_mult_max)
    const5 = alpha_spacefilling * rho_w/rho_ice

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          T_a = T_0(i,j,k)

          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !.. Vollenden des ice-rimings:
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          IF (deprate_ice(i,j,k) > 0.0 .AND. deprate_ice(i,j,k) >= rimeqcrate_ice(i,j,k)+rimeqrrate_ice(i,j,k)) THEN

            !.. Deposition ist groesser als gesamtes riming, deswegen bleibt Eis Eis:

            IF (rimeqcrate_ice(i,j,k) > 0.0 .OR. rimencrate_ice(i,j,k) > 0.0) THEN

              rime_q = rimeqcrate_ice(i,j,k)
              rime_n = rimencrate_ice(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_ice(i,j,k)   = q_ice(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4 
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n
              ENDIF

            END IF

            IF (rimeqrrate_ice(i,j,k) > 0.0 .OR. rimenrrate_ice(i,j,k) > 0.0) THEN

              rime_q = rimeqrrate_ice(i,j,k)
              rime_n = rimenrrate_ice(i,j,k)
              rime_q = MIN(q_rain(i,j,k),rime_q)
              rime_n  = MIN(n_rain(i,j,k),rime_n)
              
              q_ice(i,j,k)  = q_ice(i,j,k)  + rime_q
              q_rain(i,j,k) = q_rain(i,j,k) - rime_q
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              
              ! Eismultiplikation nach Hallet und Mossop
              
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                
                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n
              ENDIF
              
            END IF

          ELSE

            !.. Depositionsrate ist kleiner als riming. Jetzt findet ice -> graupel bzw. ice -> hail statt:

            !.. Operator-Splitting: Zuerst ice_cloud_riming behandeln:

            IF (rimeqcrate_ice(i,j,k) > 0.0 .OR. rimencrate_ice(i,j,k) > 0.0) THEN

              n_i = n_ice(i,j,k)
              x_i = MIN(MAX(q_ice(i,j,k)/(n_i+eps),ice%x_min),ice%x_max)
              D_i = ice%a_geo * x_i**ice%b_geo

              rime_q = rimeqcrate_ice(i,j,k)
              rime_n = rimencrate_ice(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_ice(i,j,k)   = q_ice(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4 
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n              
              ENDIF

              ! Umwandlung ice -> graupel

              IF (D_i > D_conv_ig) THEN
                conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*rho_ice*d_i**3/x_i - 1.0) )
                conv_q = MIN(q_ice(i,j,k),conv_q)
                x_i = MIN(MAX(q_ice(i,j,k)/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse incl. Riming
                conv_n = conv_q / MAX(x_i,x_conv_g)
                conv_n = MIN(n_ice(i,j,k),conv_n)

                q_ice(i,j,k)     = q_ice(i,j,k)     - conv_q
                q_graupel(i,j,k) = q_graupel(i,j,k) + conv_q
                n_ice(i,j,k)     = n_ice(i,j,k)     - conv_n
                n_graupel(i,j,k) = n_graupel(i,j,k) + conv_n
              ENDIF

            END IF

            !.. Operator-Splitting: Dann ice_rain_riming behandeln:

            IF (rimeqirate_ice(i,j,k) > 0.0 .OR. rimeqrrate_ice(i,j,k) > 0.0 .OR.  rimenrrate_ice(i,j,k) > 0.0) THEN
              D_id = d_id_sp(i,j,k)
              D_rd = d_rd_sp_ice(i,j,k)

              x_r = MIN(MAX(q_rain(i,j,k)/(n_rain(i,j,k)+eps),rain%x_min),rain%x_max)

              rime_qi = rimeqirate_ice(i,j,k)
              rime_qr = rimeqrrate_ice(i,j,k)
              rime_n = rimenrrate_ice(i,j,k)
              rime_n  = MIN(MIN(n_rain(i,j,k),n_ice(i,j,k)),rime_n)
              rime_qr = MIN(q_rain(i,j,k),rime_qr)
              rime_qi = MIN(q_ice(i,j,k),rime_qi)

              n_ice(i,j,k)  = n_ice(i,j,k)  - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_ice(i,j,k)  = q_ice(i,j,k)  - rime_qi
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              mult_n = 0.0  ! UB_20081120
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)
              ENDIF


              IF (T_a >= T_3) THEN
                IF (D_id > D_rd) THEN
                  ! Regen wird wieder abgeshedded
                  n_ice(i,j,k) = n_ice(i,j,k) + rime_n
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_qr / x_r
                  q_ice(i,j,k) = q_ice(i,j,k) + rime_qi     ! UB_20081120 - mult_q
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr
                ELSE
                  ! Eisteilchen schmelzen instantan zu Regen.
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_n
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr + rime_qi   ! UB_20081120 - mult_q
                END IF
              ELSE
                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n  ! UB_20081120
                q_ice(i,j,k) = q_ice(i,j,k)  + mult_q  ! UB_20081120
                IF (ice_typ < 3) THEN
                  ! Eis + angefrorenes Regenwasser ergibt Graupel:
                  n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                  q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qi + rime_qr - mult_q
                ELSE
                  IF (D_id > D_rd) THEN
                    ! Eis + angefrorenes Regenwasser ergibt Graupel:
                    n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                    q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qi + rime_qr - mult_q
                  ELSE
                    ! Eis + angefrorenes Regenwasser ergibt Hagel:
                    n_hail(i,j,k) = n_hail(i,j,k) + rime_n
                    q_hail(i,j,k) = q_hail(i,j,k) + rime_qi + rime_qr - mult_q
                  END IF
                END IF
              END IF

            END IF

          END IF

          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! Vollenden des snow-rimings:
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          IF (deprate_snow(i,j,k) > 0.0 .AND. deprate_snow(i,j,k) >= rimeqcrate_snow(i,j,k)+rimeqrrate_snow(i,j,k)) THEN

            !.. Deposition ist groesser als gesamtes riming, deswegen bleibt Schnee Schnee:

            IF (rimeqcrate_snow(i,j,k) > 0.0 .OR. rimencrate_snow(i,j,k) > 0.0) THEN

              rime_q = rimeqcrate_snow(i,j,k)
              rime_n = rimencrate_snow(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_snow(i,j,k)  = q_snow(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

            END IF

            IF (rimeqrrate_snow(i,j,k) > 0.0 .OR. rimenrrate_snow(i,j,k) > 0.0) THEN

              rime_q = rimeqrrate_snow(i,j,k)
              rime_n = rimenrrate_snow(i,j,k)
              rime_q = MIN(q_rain(i,j,k),rime_q)
              rime_n = MIN(n_rain(i,j,k),rime_n)

              q_snow(i,j,k) = q_snow(i,j,k) + rime_q
              q_rain(i,j,k) = q_rain(i,j,k) - rime_q
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n


              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

            END IF

          ELSE

            !.. Operator-Splitting: Zuerst snow_cloud_riming behandeln:

            IF (rimeqcrate_snow(i,j,k) > 0.0 .OR. rimencrate_snow(i,j,k) > 0.0) THEN

              n_s = n_snow(i,j,k)
              x_s = MIN(MAX(q_snow(i,j,k)/(n_s+eps),snow%x_min),snow%x_max)
              D_s = snow%a_geo * x_s**snow%b_geo

              rime_q = rimeqcrate_snow(i,j,k)
              rime_n = rimencrate_snow(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_snow(i,j,k)  = q_snow(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

              ! Umwandlung snow -> graupel

              IF (D_s > D_conv_sg) THEN
                conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*rho_ice*d_s**3/x_s - 1.0) )
                conv_q = MIN(q_snow(i,j,k),conv_q)
                x_s = MIN(MAX(q_snow(i,j,k)/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse incl. Riming
                conv_n = conv_q / MAX(x_s,x_conv_g)
                conv_n = MIN(n_snow(i,j,k),conv_n)

                q_snow(i,j,k)    = q_snow(i,j,k)    - conv_q
                q_graupel(i,j,k) = q_graupel(i,j,k) + conv_q
                n_snow(i,j,k)    = n_snow(i,j,k)    - conv_n
                n_graupel(i,j,k) = n_graupel(i,j,k) + conv_n
              ENDIF

            END IF

            !.. Operator-Splitting: Dann snow_rain_riming behandeln:

            IF (rimeqirate_snow(i,j,k) > 0.0 .OR. rimeqrrate_snow(i,j,k) > 0.0 .OR.  rimenrrate_snow(i,j,k) > 0.0) THEN

              D_sd = d_sd_sp(i,j,k)
              D_rd = d_rd_sp_snow(i,j,k)

              x_r = MIN(MAX(q_rain(i,j,k)/(n_rain(i,j,k)+eps),rain%x_min),rain%x_max)

              rime_qs = rimeqirate_snow(i,j,k)
              rime_qr = rimeqrrate_snow(i,j,k)
              rime_n  = rimenrrate_snow(i,j,k)
              rime_qr = MIN(q_rain(i,j,k),rime_qr)
              rime_qs = MIN(q_snow(i,j,k),rime_qs)
              rime_n  = MIN(n_rain(i,j,k),rime_n)
              rime_n  = MIN(n_snow(i,j,k),rime_n)

              n_snow(i,j,k) = n_snow(i,j,k) - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_snow(i,j,k) = q_snow(i,j,k) - rime_qs
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr


              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              mult_n = 0.0 ! UB_20081120
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.,MIN(mult_1,1.))
                mult_2 = MAX(0.,MIN(mult_2,1.))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)

              ENDIF

              IF (T_a >= T_3) THEN
                IF (D_sd > D_rd) THEN
                  ! Regen wird wieder abgeshedded
                  n_snow(i,j,k) = n_snow(i,j,k) + rime_n
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_qr / x_r
                  q_snow(i,j,k) = q_snow(i,j,k) + rime_qs ! UB_20081120 - mult_q
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr
                ELSE
                  ! Schneeflocken werden instantan zu Regen.
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_n
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr + rime_qs ! UB_20081120 - mult_q
                END IF
              ELSE
                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n ! UB_20081120
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q ! UB_20081120
                IF (ice_typ < 3) THEN
                  ! Schnee + angefrorenes Regenwasser ergibt Graupel:
                  n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                  q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qr + rime_qs - mult_q
                ELSE
                  IF (D_sd > D_rd) THEN
                    ! Schnee + angefrorenes Regenwasser ergibt Graupel:
                    n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                    q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qr + rime_qs - mult_q
                  ELSE
                    ! Schnee + angefrorenes Regenwasser ergibt Hagel:
                    n_hail(i,j,k) = n_hail(i,j,k) + rime_n
                    q_hail(i,j,k) = q_hail(i,j,k) + rime_qr + rime_qs - mult_q
                  END IF
                END IF
              END IF

            END IF

          END IF

        END DO
      END DO
    END DO

  END SUBROUTINE complete_ice_snow_riming

  SUBROUTINE graupel_melting()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Schmelzens der Graupelpartikel                          *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: q_g,n_g,x_g,d_g,v_g,T_a,N_re,e_a
    REAL            :: melt,melt_v,melt_h,melt_n,melt_q
    !REAL            :: fh_n,fv_n
    REAL            :: fh_q,fv_q
    REAL, SAVE      :: a_melt_n,b_melt_n
    REAL, SAVE      :: a_melt_q,b_melt_q

    REAL, PARAMETER :: eps   = 1.e-20

    IF (firstcall.NE.1) THEN
      a_melt_n = vent_coeff_a(graupel,0)
      b_melt_n = vent_coeff_b(graupel,0)
      a_melt_q = vent_coeff_a(graupel,1)
      b_melt_q = vent_coeff_b(graupel,1)
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k)
          q_g = q_graupel(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (T_a > T_3 .AND. q_g > 0.0) THEN
            e_a = esl(T_a)                                     !..Saettigungsdampfdruck
            n_g = n_graupel(i,j,k)                              !..Anzahldichte in SI

            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max)  !..mittlere Masse in SI     

            D_g = graupel%a_geo * x_g**graupel%b_geo                   !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

            N_re = v_g * D_g / nu_l                             !..mittlere Reynoldszahl
            fv_q = a_melt_q + b_melt_q * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf
  !          fv_n = a_melt_n + b_melt_n * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf

! UB_20081125:            D_T  = K_T / (cp * rho_0(i,j,k))
! UB_20081125:            fh_q = D_T / D_v * fv_q  ! <-- gives wrong values!
! UB_20081125: After formulas of Rasmussen and Heymsfield (1987) the ratio fh_q / fv_q is approx. 1.05
!              for a wide temperature- and pressure range:
            fh_q = 1.05 * fv_q
            !           fh_n = D_T / D_v * fv_n

            melt   = 2.0*pi / L_ew * D_g * n_g * dt

            melt_h = melt * K_T * (T_a - T_3)
            melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)

            melt_q = (melt_h * fh_q + melt_v * fv_q)
!            melt_n = (melt_h * fh_n + melt_v * fv_n) / x_g

! ub>> setzte melt_n so, dass x_h beim Schmelzvorgang erhalten bleibt:
            melt_n = MIN(MAX( (melt_q - q_g) / x_g + n_g, 0.0), n_g)

            melt_q = MIN(q_g,melt_q)
            melt_n = MIN(n_g,melt_n)

            melt_q = MAX(0.,melt_q)
            melt_n = MAX(0.,melt_n)

            q_graupel(i,j,k) = q_graupel(i,j,k) - melt_q
            q_rain(i,j,k)    = q_rain(i,j,k)    + melt_q

            n_graupel(i,j,k) = n_graupel(i,j,k) - melt_n
            n_rain(i,j,k)    = n_rain(i,j,k)    + melt_n
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_melting

  SUBROUTINE hail_melting()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Schmelzens der Hagelpartikel                            *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: q_h,n_h,x_h,d_h,v_h,T_a,N_re,e_a
    REAL            :: melt,melt_v,melt_h,melt_n,melt_q
    !REAL            :: fh_n,fv_n
    REAL            :: fh_q,fv_q
    REAL, SAVE      :: a_melt_n,b_melt_n
    REAL, SAVE      :: a_melt_q,b_melt_q

    REAL, PARAMETER :: eps   = 1.e-20

    IF (firstcall.NE.1) THEN
      a_melt_n = vent_coeff_a(hail,0)
      b_melt_n = vent_coeff_b(hail,0)
      a_melt_q = vent_coeff_a(hail,1)
      b_melt_q = vent_coeff_b(hail,1)
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k)
          q_h = q_hail(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (T_a > T_3 .AND. q_h > 0.0) THEN
            e_a = esl(T_a)                                     !..Saettigungsdampfdruck
            n_h = n_hail(i,j,k)                              !..Anzahldichte in SI

            x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max)  !..mittlere Masse in SI     

            D_h = hail%a_geo * x_h**hail%b_geo                   !..mittlerer Durchmesser
            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

            N_re = v_h * D_h / nu_l                             !..mittlere Reynoldszahl
            fv_q = a_melt_q + b_melt_q * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf
!            fv_n = a_melt_n + b_melt_n * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf

! UB_20081125:            D_T  = K_T / (cp * rho_0(i,j,k))
! UB_20081125:            fh_q = D_T / D_v * fv_q  ! <-- gives wrong values!
! UB_20081125: After formulas of Rasmussen and Heymsfield (1987) the ratio fh_q / fv_q is approx. 1.05
!              for a wide temperature- and pressure range:
            fh_q = 1.05 * fv_q
!            fh_n = D_T / D_v * fv_n

            melt   = 2.0*pi / L_ew * D_h * n_h * dt

            melt_h = melt * K_T * (T_a - T_3)
            melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)

            melt_q = (melt_h * fh_q + melt_v * fv_q)
!            melt_n = (melt_h * fh_n + melt_v * fv_n) / x_h

! ub>> setzte melt_n so, dass x_h beim Schmelzvorgang erhalten bleibt:
            melt_n = MIN(MAX( (melt_q - q_h) / x_h + n_h, 0.0), n_h)

            melt_q = MIN(q_h,melt_q)
            melt_n = MIN(n_h,melt_n)

            melt_q = MAX(0.,melt_q)
            melt_n = MAX(0.,melt_n)

            q_hail(i,j,k) = q_hail(i,j,k) - melt_q
            q_rain(i,j,k) = q_rain(i,j,k)    + melt_q

            n_hail(i,j,k) = n_hail(i,j,k) - melt_n
            n_rain(i,j,k) = n_rain(i,j,k)    + melt_n
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE hail_melting

  SUBROUTINE cloud_freeze ()
    !*******************************************************************************
    !                                                                              *
    ! Diese Routine behandelt das Gefrieren von Wolkentropfen                      *
    !                                                                              *
    !*******************************************************************************

    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER          :: i, j, k
    REAL :: fr_q, fr_n,T_a,q_c,x_c,n_c,j_het,j_hom,j_tot,T_c
    REAL, PARAMETER :: eps   = 1.e-20 
    REAL            :: facg

    facg = moment_gamma(cloud,2)    ! <hn

    !..Test auf Schmelzen oder Gefrieren von Wolkenteilchen
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a  = T_0(i,j,k)
          IF (T_a < T_3) THEN
            q_c = q_cloud(i,j,k)
            n_c = n_cloud(i,j,k)
            T_c = T_a - T_3
            IF (q_c > 0.0) THEN
              IF (T_c < -50.0) THEN            
                fr_q = q_c                                               !..Komplettes hom. Gefrieren
                fr_n = n_c                                               !..unterhalb -50 C
              ELSE                                      
                x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max)    !..mittlere Masse
                
                !..Hom. Gefrieren nach Jeffrey und Austin (1997), siehe auch Cotton und Field (2001)
                IF (T_c > -30.0) THEN            
                  j_hom = 1.0e6 * EXP(-7.63-2.996*(T_c+30.0))           !..J in 1/(m3 s) 
                ELSE
                  j_hom = 1.0e6 * EXP(-243.4-14.75*T_c-0.307*T_c**2-0.00287*T_c**3-0.0000102*T_c**4)
                ENDIF

                !..Het. Gefrieren: stochastisches Modell nach Bigg, Daten nach Barklie and Gokhale
                !j_het = b_HET * ( EXP( - a_HET * T_c) - 1.0 )            !..J in 1/(m3 s)
                j_het = 0.0 ! neglected for cloud droplets
                
                !..Umwandlungsraten fuer Anzahl- und Massendichte
                j_tot = (j_hom + j_het) / rho_w * dt                     !..J*dt in 1/kg
                fr_n  = j_tot * q_c
                fr_q  = j_tot * q_c * x_c * facg
              END IF
              fr_q  = MIN(fr_q,q_c)
              fr_n  = MIN(fr_n,n_c)

              !..Berechnung der H2O-Komponenten
              q_cloud(i,j,k) = q_cloud(i,j,k) - fr_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - fr_n

              fr_n  = MAX(fr_n,fr_q/cloud%x_max)
              
              IF (nuc_c_typ .EQ. 0) THEN
                fr_n = MAX(MIN(fr_n,n_cloud(i,j,k)-n_ice(i,j,k)),0.)
              ENDIF
              
              q_ice(i,j,k)   = q_ice(i,j,k)   + fr_q
              n_ice(i,j,k)   = n_ice(i,j,k)   + fr_n

              ENDIF

          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE cloud_freeze

  SUBROUTINE rain_freeze ()
    !*******************************************************************************
    !                                                                              *
    ! Diese Routine behandelt das Gefrieren von Regentropfen                       *
    !                                                                              *
    !*******************************************************************************

    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0
    REAL            :: fr_q,fr_n,T_a,q_r,x_r,n_r,j_het,&
         &  fr_q_i,fr_n_i,fr_q_g,fr_n_g,fr_q_h,fr_n_h,n_0,lam,xmax_ice,xmax_gr,fr_q_tmp,fr_n_tmp
    REAL, SAVE      :: coeff_z
    REAL, PARAMETER :: D_rainfrz_ig = 1.00e-3 ! rain --> ice oder graupel
    REAL, PARAMETER :: D_rainfrz_gh = 1.25e-3 ! rain --> graupel oder hail
    REAL, PARAMETER :: q_krit_fr = 1.000e-6 ! q-Schwellenwert fuer rain_freeze
    REAL, PARAMETER :: a_HET = 6.5e-1 ! Messung nach Barklie and Gokhale (PK S.350)
    REAL, PARAMETER :: b_HET = 2.0e+2 ! Messung nach Barklie and Gokhale (PK S.350)
    REAL, PARAMETER :: eps  = 1.e-20


    IF (firstcall.NE.1) THEN
      !..Koeff. fuer Reflektivitaet Z (2. Moment)
      coeff_z = moment_gamma(rain,2)
      firstcall = 1       
    ENDIF

    xmax_ice = (D_rainfrz_ig/rain%a_geo)**(1.0/rain%b_geo)
    xmax_gr  = (D_rainfrz_gh/rain%a_geo)**(1.0/rain%b_geo)

    !..Test auf Schmelzen oder Gefrieren von Regentropfen
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k)
          q_r = q_rain(i,j,k)
          n_r = n_rain(i,j,k)

          IF (T_a < T_3) THEN
            IF (q_r <= q_krit_fr) THEN
              IF (T_a < T_f) THEN
                fr_q = q_r                  !  Ausfrieren unterhalb T_f \approx -40 C
                fr_n = n_r
                fr_n_i= n_r
                fr_q_i= q_r
                fr_n_g= 0.0
                fr_q_g= 0.0
                fr_n_h= 0.0
                fr_q_h= 0.0
! UB_20080220>
                fr_n_tmp = 1.0
                fr_q_tmp = 1.0
! <UB_20080220
              ELSE
                fr_q = 0.0
                fr_n = 0.0
                fr_n_i= 0.0
                fr_q_i= 0.0
                fr_n_g= 0.0
                fr_q_g= 0.0
                fr_n_h= 0.0
                fr_q_h= 0.0
! UB_20080220>
                fr_n_tmp = 0.0
                fr_q_tmp = 0.0
! <UB_20080220
              END IF
            ELSE
              x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
              n_r = q_r / x_r
              IF (T_a < T_f) THEN            !..Nur Eis
                fr_q = q_r                  !  Ausfrieren unterhalb T_f \approx -40 C
                fr_n = n_r

                ! ub>> Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                !      oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                !      bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                !      und von xmax_gr bis unendlich (--> Hagel).
                
                lam = ( gfct((rain%nu+1.0)/rain%mu) / gfct((rain%nu+2.0)/rain%mu) * x_r)**(-rain%mu)
                n_0 = rain%mu * n_r * lam**((rain%nu+1.0)/rain%mu) / gfct((rain%nu+1.0)/rain%mu)
                fr_n_i = n_0/(rain%mu*lam**((rain%nu+1.0)/rain%mu))* &
                    incgfct_lower((rain%nu+1.0)/rain%mu, lam*xmax_ice**rain%mu)
                fr_q_i = n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                    incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_ice**rain%mu)
                fr_n_g = n_0/(rain%mu*lam**((rain%nu+1.0)/rain%mu))* &
                    incgfct_lower((rain%nu+1.0)/rain%mu, lam*xmax_gr**rain%mu)
                fr_q_g = n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                    incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_gr**rain%mu)
                
                fr_n_h = fr_n - fr_n_g
                fr_q_h = fr_q - fr_q_g
                fr_n_g = fr_n_g - fr_n_i
                fr_q_g = fr_q_g - fr_q_i
                fr_n_tmp = n_r/MAX(fr_n,n_r)
                fr_q_tmp = q_r/MAX(fr_q,q_r)

              ELSE                           !..Heterogenes Gefrieren
                j_het = MAX(b_HET * ( EXP( a_HET * (T_3 - T_a)) - 1.0 ),0.) / rho_w * dt

                ! ub>> Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                !      oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                !      bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                !      und von xmax_gr bis unendlich (--> Hagel).
                
                IF (j_het >= 1e-20) THEN
                  fr_n  = j_het * q_r
                  fr_q  = j_het * q_r * x_r * coeff_z
                  
                  lam = ( gfct((rain%nu+1.0)/rain%mu) / gfct((rain%nu+2.0)/rain%mu) * x_r)**(-rain%mu)
                  n_0 = rain%mu * n_r * lam**((rain%nu+1.0)/rain%mu) / gfct((rain%nu+1.0)/rain%mu)
                  fr_n_i = j_het * n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                      incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_ice**rain%mu)
                  fr_q_i = j_het * n_0/(rain%mu*lam**((rain%nu+3.0)/rain%mu))* &
                      incgfct_lower((rain%nu+3.0)/rain%mu, lam*xmax_ice**rain%mu)
                  fr_n_g = j_het * n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                      incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_gr**rain%mu)
                  fr_q_g = j_het * n_0/(rain%mu*lam**((rain%nu+3.0)/rain%mu))* &
                      incgfct_lower((rain%nu+3.0)/rain%mu, lam*xmax_gr**rain%mu)
                  
                  fr_n_h = fr_n - fr_n_g
                  fr_q_h = fr_q - fr_q_g
                  fr_n_g = fr_n_g - fr_n_i
                  fr_q_g = fr_q_g - fr_q_i
                  fr_n_tmp = n_r/MAX(fr_n,n_r)
                  fr_q_tmp = q_r/MAX(fr_q,q_r)
                ELSE
                  fr_n= 0.0
                  fr_q= 0.0
                  fr_n_i= 0.0
                  fr_q_i= 0.0
                  fr_n_g= 0.0
                  fr_q_g= 0.0
! UB_20080212>
                  fr_n_h= 0.0
                  fr_q_h= 0.0
! <UB_20080212
                  fr_n_tmp = 0.0
                  fr_q_tmp = 0.0
                END IF
                
              END IF

              fr_n = fr_n * fr_n_tmp
              fr_q = fr_q * fr_q_tmp

              fr_n_i = fr_n_i * fr_n_tmp
              fr_n_g = fr_n_g * fr_n_tmp
              fr_n_h = fr_n_h * fr_n_tmp
              fr_q_i = fr_q_i * fr_q_tmp
              fr_q_g = fr_q_g * fr_q_tmp
              fr_q_h = fr_q_h * fr_q_tmp

            END IF

            !..Berechnung der H2O-Komponenten

            q_rain(i,j,k) = q_rain(i,j,k) - fr_q
            n_rain(i,j,k) = n_r - fr_n

            IF (q_rain(i,j,k) < 0.0) THEN
              q_rain(i,j,k) = 0.0
            END IF
            IF (n_rain(i,j,k) < 0.0) THEN
              n_rain(i,j,k) = 0.0
            END IF            

            IF (ice_typ < 2) THEN 
              ! ohne Hagelklasse,  gefrierender Regen wird Eis oder Graupel
              q_ice(i,j,k) = q_ice(i,j,k)  + fr_q_i
              n_ice(i,j,k) = n_ice(i,j,k)  + fr_n_i
              q_graupel(i,j,k) = q_graupel(i,j,k)  + fr_q_h + fr_q_g
              n_graupel(i,j,k) = n_graupel(i,j,k)  + fr_n_h + fr_n_g
            ELSE
              ! mit Hagelklasse, gefrierender Regen wird Eis, Graupel oder Hagel
              q_ice(i,j,k) = q_ice(i,j,k)  + fr_q_i
              n_ice(i,j,k) = n_ice(i,j,k)  + fr_n_i
              q_graupel(i,j,k) = q_graupel(i,j,k)  + fr_q_g
              n_graupel(i,j,k) = n_graupel(i,j,k)  + fr_n_g
              q_hail(i,j,k) = q_hail(i,j,k)  + fr_q_h
              n_hail(i,j,k) = n_hail(i,j,k)  + fr_n_h
            ENDIF

            IF (q_ice(i,j,k) < 0.0) THEN
              q_ice(i,j,k) = 0.0
            END IF
            IF (n_ice(i,j,k) < 0.0) THEN
              n_ice(i,j,k) = 0.0
            END IF            

            IF (q_graupel(i,j,k) < 0.0) THEN
              q_graupel(i,j,k) = 0.0
            END IF
            IF (n_graupel(i,j,k) < 0.0) THEN
              n_graupel(i,j,k) = 0.0
            END IF            

            IF (q_hail(i,j,k) < 0.0) THEN
              q_hail(i,j,k) = 0.0
            END IF
            IF (n_hail(i,j,k) < 0.0) THEN
              n_hail(i,j,k) = 0.0
            END IF            

          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE rain_freeze

  SUBROUTINE rain_freeze_old ()
    !*******************************************************************************
    !                                                                              *
    ! Diese Routine behandelt das Gefrieren von Regentropfen                       *
    !                                                                              *
    !*******************************************************************************

    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0
    REAL            :: fr_q,fr_n,T_a,q_r,x_r,n_r,j_het
    REAL, SAVE      :: coeff_z
    REAL, PARAMETER :: q_krit_fr = 1.000e-6 ! q-Schwellenwert fuer rain_freeze
    REAL, PARAMETER :: a_HET = 6.5e-1 ! Messung nach Barklie and Gokhale (PK S.350)
    REAL, PARAMETER :: b_HET = 2.0e+2 ! Messung nach Barklie and Gokhale (PK S.350)
    REAL, PARAMETER :: eps  = 1.e-20

    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rate_q, rate_n

    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), &
         rate_n(0:loc_ix,1:loc_iy,1:loc_iz))
    rate_q = 0.0
    rate_n = 0.0

    IF (firstcall.NE.1) THEN
      !..Koeff. fuer Reflektivitaet Z (2. Moment)
      coeff_z = moment_gamma(rain,2)
      firstcall = 1       
    ENDIF

    !..Test auf Schmelzen oder Gefrieren von Regentropfen
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k)
          q_r = q_rain(i,j,k)
          n_r = n_rain(i,j,k)
          IF (T_a < T_3 .AND. q_r > q_krit_fr) THEN
            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
            IF (T_a < T_f) THEN            !..Nur Eis
              fr_q = q_r                  !  Ausfrieren unterhalb T_f \approx -40 C
              fr_n = n_r
            ELSE                           !..Heterogenes Gefrieren
              j_het = MAX(b_HET * ( EXP( a_HET * (T_3 - T_a)) - 1.0 ),0.) / rho_w * dt
              fr_n  = j_het * q_r
              fr_q  = j_het * q_r * x_r * coeff_z
            END IF
            fr_q  = MIN(fr_q,q_r)
            fr_n  = MIN(fr_n,n_r)

            rate_q(i,j,k) = fr_q
            rate_n(i,j,k) = fr_n

          END IF
        END DO
      END DO
    END DO

    !..Berechnung der H2O-Komponenten
    q_rain = q_rain - rate_q
    n_rain = n_rain - rate_n

    IF (ice_typ < 2) THEN 
      ! ohne Hagelklasse,  gefrierender Regen wird Graupel
      q_graupel = q_graupel  + rate_q
      n_graupel = n_graupel  + rate_n
    ELSE
      ! mit Hagelklasse, gefrierender Regen wird Hagel
      q_hail = q_hail  + rate_q
      n_hail = n_hail  + rate_n
    ENDIF

    DEALLOCATE(rate_q, rate_n)

  END SUBROUTINE rain_freeze_old

  SUBROUTINE ice_selfcollection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Selbsteinfangs der Eispartikel                          *
    !                                                                              *
    !*******************************************************************************

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a             !..absolute Temperatur
    REAL            :: q_i,n_i,x_i,d_i,v_i,e_coll,x_conv
    REAL            :: self_n,self_q
    REAL            :: delta_n_11,delta_n_12,delta_n_22
    REAL            :: delta_q_11,delta_q_12,delta_q_22
    REAL            :: theta_n_11,theta_n_12,theta_n_22
    REAL            :: theta_q_11,theta_q_12,theta_q_22
    REAL,SAVE       :: delta_n,delta_q
    REAL,SAVE       :: theta_n,theta_q

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_11 = coll_delta_11(ice,ice,0)
      delta_n_12 = coll_delta_12(ice,ice,0)
      delta_n_22 = coll_delta_22(ice,ice,0)
      delta_q_11 = coll_delta_11(ice,ice,0) 
      delta_q_12 = coll_delta_12(ice,ice,1)
      delta_q_22 = coll_delta_22(ice,ice,1)

      theta_n_11 = coll_theta_11(ice,ice,0)
      theta_n_12 = coll_theta_12(ice,ice,0)
      theta_n_22 = coll_theta_22(ice,ice,0)
      theta_q_11 = coll_theta_11(ice,ice,0)
      theta_q_12 = coll_theta_12(ice,ice,1)
      theta_q_22 = coll_theta_22(ice,ice,1)

      delta_n = delta_n_11 + delta_n_12 + delta_n_22
      delta_q = delta_q_11 + delta_q_12 + delta_q_22
      theta_n = theta_n_11 - theta_n_12 + theta_n_22
      theta_q = theta_q_11 - theta_q_12 + theta_q_22

      firstcall = 1
    ENDIF

    x_conv = (D_conv_ii/snow%a_geo)**(1./snow%b_geo)

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_i = q_ice(i,j,k)                                   !..Fluessigwassergehalt in SI
          n_i = n_ice(i,j,k)                                   !..Anzahldichte in SI
          x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse in SI     
          D_i = ice%a_geo * x_i**ice%b_geo                     !..mittlerer Durchmesser
          IF ( n_i > 0.0 .AND. q_i > q_krit_ii .AND. D_i > D_krit_ii ) THEN

            T_a = T_0(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige Efficiency nach Cotton et al. (1986) 
              !   (siehe auch Straka, 1989; S. 53)
              e_coll = MIN(10**(0.035*(T_a-T_3)-0.7),0.2)
            END IF


            v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)      !..mittlere Sedimentationsgeschw.

            self_n = pi * 0.25 * e_coll * delta_n * n_i * n_i * D_i * D_i &
                 & * ( theta_n * v_i * v_i + 2.0*ice_s_vel**2 )**0.5 * dt

            self_q = pi * 0.25 * e_coll * delta_q * n_i * q_i * D_i * D_i &
                 & * ( theta_q * v_i * v_i + 2.0*ice_s_vel**2 )**0.5 * dt

            self_q = MIN(self_q,q_i)
            self_n = MIN(MIN(self_n,self_q/x_conv),n_i)

            q_ice(i,j,k)  = q_ice(i,j,k)  - self_q
            q_snow(i,j,k) = q_snow(i,j,k) + self_q

            n_ice(i,j,k)  = n_ice(i,j,k)  - self_n
            n_snow(i,j,k) = n_snow(i,j,k) + self_n / 2.0

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE ice_selfcollection
 
  SUBROUTINE snow_selfcollection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Selbsteinfangs der Eispartikel                          *
    !                                                                              *
    !*******************************************************************************

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a             !..Absolute Temperatur
    REAL            :: q_s,n_s,x_s,d_s,v_s,e_coll
    REAL            :: self_n
    REAL            :: delta_n_11,delta_n_12
    REAL            :: theta_n_11,theta_n_12
    REAL,SAVE       :: delta_n
    REAL,SAVE       :: theta_n

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_11 = coll_delta_11(snow,snow,0)
      delta_n_12 = coll_delta_12(snow,snow,0)

      theta_n_11 = coll_theta_11(snow,snow,0)
      theta_n_12 = coll_theta_12(snow,snow,0)

      delta_n = (2.0*delta_n_11 + delta_n_12)
      theta_n = (2.0*theta_n_11 - theta_n_12)

      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_s = q_snow(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF ( q_s > q_krit ) THEN
            T_a = T_0(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige Efficiency nach Cotton et al. (1986) 
              !   (siehe auch Straka, 1989; S. 53)
              !e_coll = MIN(10**(0.035*(T_a-T_3)-0.7),0.2)
              !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
              !e_coll = MIN(exp(0.09*(T_a-T_3)),1.0)
              e_coll = MAX(0.1,MIN(EXP(0.09*(T_a-T_3)),1.0))
            ENDIF

            n_s = n_snow(i,j,k)                        !..Anzahldichte in SI
            x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse in SI     

            D_s = snow%a_geo * x_s**snow%b_geo                     !..mittlerer Durchmesser
            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k)    !..mittlere Sedimentationsgeschw.

            !self_n = pi/4.0 * e_coll * n_s**2 * delta_n * D_s**2 * &
            !     &                ( theta_n * v_s**2 + 2.0 * snow_s_vel**2 )**0.5 * dt

            ! >>hn
            ! fehlenden Faktor 1/2 ergaenzt                                
            self_n = pi * 0.125 * e_coll * n_s * n_s * delta_n * D_s * D_s * &
                 &                ( theta_n * v_s * v_s + 2.0 * snow_s_vel**2 )**0.5 * dt
            ! << hn

            self_n = MIN(self_n,n_s)

            n_snow(i,j,k) = n_snow(i,j,k) - self_n

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE snow_selfcollection

  SUBROUTINE snow_melting()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Schmelzens der Schneepartikel                           *
    !                                                                              *
    !*******************************************************************************

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: q_s,n_s,x_s,d_s,v_s,T_a,N_re,e_a
    REAL            :: melt,melt_v,melt_h,melt_n,melt_q
    !REAL            :: fh_n,fv_n
    REAL            :: fh_q,fv_q
    REAL, SAVE      :: a_melt_n,b_melt_n
    REAL, SAVE      :: a_melt_q,b_melt_q

    REAL, PARAMETER :: eps   = 1.e-20

    IF (firstcall.NE.1) THEN
      a_melt_n = vent_coeff_a(snow,0)
      b_melt_n = vent_coeff_b(snow,0)
      a_melt_q = vent_coeff_a(snow,1)
      b_melt_q = vent_coeff_b(snow,1)

      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k)
          q_s = q_snow(i,j,k)                                    !..Fluessigwassergehalt in SI
          IF (T_a > T_3 .AND. q_s > 0.0) THEN
            e_a = esl(T_a)                                     !..Saettigungsdampfdruck
            n_s = n_snow(i,j,k)                                 !..Anzahldichte in SI

            x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max) !..mittlere Masse in SI     

            D_s = snow%a_geo * x_s**snow%b_geo                  !..mittlerer Durchmesser
            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            N_re = v_s * D_s / nu_l                             !..mittlere Reynoldszahl
            fv_q = a_melt_q + b_melt_q * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf
!            fv_n = a_melt_n + b_melt_n * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf

! UB_20081111:           D_T  = K_T/(cp * (rho_0(i,j,k)+rho_g(i,j,k)))
! UB_20081125:            D_T  = K_T / (cp * rho_0(i,j,k))
! UB_20081125:            fh_q = D_T/D_v * fv_q  ! <-- gives wrong values!
! UB_20081125: After formulas of Rasmussen and Heymsfield (1987) the ratio fh_q / fv_q is approx. 1.05
!              for a wide temperature- and pressure range:
            fh_q = 1.05 * fv_q
!            fh_n = D_T/D_v * fv_n

            melt   = 2.0*pi/L_ew * D_s * n_s * dt

            melt_h = melt * K_T * (T_a - T_3)
            melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)

            melt_q = (melt_h * fh_q + melt_v * fv_q)
!            melt_n = (melt_h * fh_n + melt_v * fv_n) / x_s

! ub>> setzte melt_n so, dass x_s beim Schmelzvorgang erhalten bleibt:
            melt_n = MIN(MAX( (melt_q - q_s) / x_s + n_s, 0.0), n_s)
 
            melt_q = MIN(q_s,MAX(melt_q,0.))
            melt_n = MIN(n_s,MAX(melt_n,0.))

            IF (T_a - T_3 > 10.0) THEN
              melt_q = q_s
              melt_n = n_s
            ENDIF

            q_snow(i,j,k) = q_snow(i,j,k) - melt_q
            q_rain(i,j,k) = q_rain(i,j,k) + melt_q

            n_snow(i,j,k) = n_snow(i,j,k) - melt_n
            n_rain(i,j,k) = n_rain(i,j,k) + melt_n

            n_snow(i,j,k) = MAX(n_snow(i,j,k), q_snow(i,j,k)/snow%x_max)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE snow_melting

  SUBROUTINE graupel_snow_collection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Koagulation von Graupel und Schnee                      *
    !                                                                              *
    !*******************************************************************************

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a
    REAL            :: q_g,n_g,x_g,d_g,v_g
    REAL            :: q_s,n_s,x_s,d_s,v_s
    REAL            :: coll_n,coll_q,e_coll
    REAL, SAVE      :: delta_n_gg,delta_n_gs,delta_n_ss
    REAL, SAVE      :: delta_q_gg,delta_q_gs,delta_q_ss
    REAL, SAVE      :: theta_n_gg,theta_n_gs,theta_n_ss
    REAL, SAVE      :: theta_q_gg,theta_q_gs,theta_q_ss

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,snow,0)
      delta_n_gs = coll_delta_12(graupel,snow,0)
      delta_n_ss = coll_delta_22(graupel,snow,0)
      delta_q_gg = coll_delta_11(graupel,snow,0) 
      delta_q_gs = coll_delta_12(graupel,snow,1)
      delta_q_ss = coll_delta_22(graupel,snow,1)

      theta_n_gg = coll_theta_11(graupel,snow,0)
      theta_n_gs = coll_theta_12(graupel,snow,0)
      theta_n_ss = coll_theta_22(graupel,snow,0)
      theta_q_gg = coll_theta_11(graupel,snow,0)
      theta_q_gs = coll_theta_12(graupel,snow,1)
      theta_q_ss = coll_theta_22(graupel,snow,1)

      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_s = q_snow(i,j,k)                                    !..Fluessigwassergehalt in SI
          q_g = q_graupel(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (q_s > q_krit .AND. q_g > q_krit) THEN
            T_a = T_0(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
              e_coll = MIN(EXP(0.09*(T_a-T_3)),1.0)
            ENDIF

            n_s = n_snow(i,j,k)                                       !..Anzahldichte
            n_g = n_graupel(i,j,k)                                    !..Anzahldichte

            x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max)       !..mittlere Masse
            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse

            D_s = snow%a_geo * x_s**snow%b_geo                        !..mittlerer Durchmesser
            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k)       !..mittlere Sedimentationsgeschw.

            D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            coll_n = pi/4.0 * n_g * n_s * e_coll * dt & 
                 &   * (delta_n_gg * D_g**2 + delta_n_gs * D_g*D_s + delta_n_ss * D_s**2) &
                 &   * (theta_n_gg * v_g**2 - theta_n_gs * v_g*v_s + theta_n_ss * v_s**2  &
                 &     +snow_s_vel**2)**0.5

            coll_q = pi/4.0 * n_g * q_s * e_coll * dt & 
                 &   * (delta_q_gg * D_g**2 + delta_q_gs * D_g*D_s + delta_q_ss * D_s**2) &
                 &   * (theta_q_gg * v_g**2 - theta_q_gs * v_g*v_s + theta_q_ss * v_s**2  &
                 &     +snow_s_vel**2)**0.5

            coll_n = MIN(n_s,coll_n)
            coll_q = MIN(q_s,coll_q)

            q_graupel(i,j,k) = q_graupel(i,j,k) + coll_q
            q_snow(i,j,k)    = q_snow(i,j,k)    - coll_q
            n_snow(i,j,k)    = n_snow(i,j,k)    - coll_n

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_snow_collection


  ! hn>>
  SUBROUTINE hail_snow_collection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Koagulation von Hagel und Schnee                      *
    !                                                                              *
    !*******************************************************************************

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a
    REAL            :: q_h,n_h,x_h,d_h,v_h
    REAL            :: q_s,n_s,x_s,d_s,v_s
    REAL            :: coll_n,coll_q,e_coll
    REAL, SAVE      :: delta_n_hh,delta_n_hs,delta_n_ss
    REAL, SAVE      :: delta_q_hh,delta_q_hs,delta_q_ss
    REAL, SAVE      :: theta_n_hh,theta_n_hs,theta_n_ss
    REAL, SAVE      :: theta_q_hh,theta_q_hs,theta_q_ss

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,snow,0)
      delta_n_hs = coll_delta_12(hail,snow,0)
      delta_n_ss = coll_delta_22(hail,snow,0)
      delta_q_hh = coll_delta_11(hail,snow,0) 
      delta_q_hs = coll_delta_12(hail,snow,1)
      delta_q_ss = coll_delta_22(hail,snow,1)

      theta_n_hh = coll_theta_11(hail,snow,0)
      theta_n_hs = coll_theta_12(hail,snow,0)
      theta_n_ss = coll_theta_22(hail,snow,0)
      theta_q_hh = coll_theta_11(hail,snow,0)
      theta_q_hs = coll_theta_12(hail,snow,1)
      theta_q_ss = coll_theta_22(hail,snow,1)

      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_s = q_snow(i,j,k)                                 !..Fluessigwassergehalt in SI
          q_h = q_hail(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (q_s > q_krit .AND. q_h > q_krit) THEN
            T_a = T_0(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
              e_coll = MIN(EXP(0.09*(T_a-T_3)),1.0)
            ENDIF

            n_s = n_snow(i,j,k)                                    !..Anzahldichte
            n_h = n_hail(i,j,k)                                    !..Anzahldichte

            x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max) !..mittlere Masse
            x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max) !..mittlere Masse

            D_s = snow%a_geo * x_s**snow%b_geo                        !..mittlerer Durchmesser
            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k)       !..mittlere Sedimentationsgeschw.

            D_h = hail%a_geo * x_h**hail%b_geo                  !..mittlerer Durchmesser
            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            coll_n = pi/4.0 * n_h * n_s * e_coll * dt & 
                 &   * (delta_n_hh * D_h**2 + delta_n_hs * D_h*D_s + delta_n_ss * D_s**2) &
                 &   * (theta_n_hh * v_h**2 - theta_n_hs * v_h*v_s + theta_n_ss * v_s**2  &
                 &     +snow_s_vel**2)**0.5

            coll_q = pi/4.0 * n_h * q_s * e_coll * dt & 
                 &   * (delta_q_hh * D_h**2 + delta_q_hs * D_h*D_s + delta_q_ss * D_s**2) &
                 &   * (theta_q_hh * v_h**2 - theta_q_hs * v_h*v_s + theta_q_ss * v_s**2  &
                 &     +snow_s_vel**2)**0.5

            coll_n = MIN(n_s,coll_n)
            coll_q = MIN(q_s,coll_q)

            q_hail(i,j,k) = q_hail(i,j,k) + coll_q
            q_snow(i,j,k) = q_snow(i,j,k) - coll_q
            n_snow(i,j,k) = n_snow(i,j,k) - coll_n

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE hail_snow_collection
  ! <<hn


  SUBROUTINE graupel_ice_collection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Koagulation von Graupel und Schnee                      *
    !                                                                              *
    !*******************************************************************************

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a
    REAL            :: q_g,n_g,x_g,d_g,v_g
    REAL            :: q_i,n_i,x_i,d_i,v_i
    REAL            :: coll_n,coll_q,e_coll
    REAL, SAVE      :: delta_n_gg,delta_n_gi,delta_n_ii
    REAL, SAVE      :: delta_q_gg,delta_q_gi,delta_q_ii
    REAL, SAVE      :: theta_n_gg,theta_n_gi,theta_n_ii
    REAL, SAVE      :: theta_q_gg,theta_q_gi,theta_q_ii

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,ice,0)
      delta_n_gi = coll_delta_12(graupel,ice,0)
      delta_n_ii = coll_delta_22(graupel,ice,0)
      delta_q_gg = coll_delta_11(graupel,ice,0) 
      delta_q_gi = coll_delta_12(graupel,ice,1)
      delta_q_ii = coll_delta_22(graupel,ice,1)

      theta_n_gg = coll_theta_11(graupel,ice,0)
      theta_n_gi = coll_theta_12(graupel,ice,0)
      theta_n_ii = coll_theta_22(graupel,ice,0)
      theta_q_gg = coll_theta_11(graupel,ice,0)
      theta_q_gi = coll_theta_12(graupel,ice,1)
      theta_q_ii = coll_theta_22(graupel,ice,1)

      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_i = q_ice(i,j,k)                                     !..Fluessigwassergehalt in SI
          q_g = q_graupel(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (q_i > q_krit .AND. q_g > q_krit) THEN
            T_a = T_0(i,j,k)

            !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              e_coll = MIN(EXP(0.09*(T_a-T_3)),1.0)
            END IF

            n_i = n_ice(i,j,k)                                        !..Anzahldichte
            n_g = n_graupel(i,j,k)                                    !..Anzahldichte

            x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)         !..mittlere Masse
            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse

            D_i = ice%a_geo * x_i**ice%b_geo                          !..mittlerer Durchmesser
            v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)         !..mittlere Sedimentationsgeschw.

            D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            coll_n = pi/4.0 * n_g * n_i * e_coll * dt & 
                 &   * (delta_n_gg * D_g**2 + delta_n_gi * D_g*D_i + delta_n_ii * D_i**2) &
                 &   * (theta_n_gg * v_g**2 - theta_n_gi * v_g*v_i + theta_n_ii * v_i**2)**0.5

            coll_q = pi/4.0 * n_g * q_i * e_coll * dt & 
                 &   * (delta_q_gg * D_g**2 + delta_q_gi * D_g*D_i + delta_q_ii * D_i**2) &
                 &   * (theta_q_gg * v_g**2 - theta_q_gi * v_g*v_i + theta_q_ii * v_i**2)**0.5

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            q_graupel(i,j,k) = q_graupel(i,j,k) + coll_q
            q_ice(i,j,k)     = q_ice(i,j,k)     - coll_q
            n_ice(i,j,k)     = n_ice(i,j,k)     - coll_n

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_ice_collection

  ! hn>>
  SUBROUTINE hail_ice_collection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Koagulation von Hagel und Schnee                        *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a
    REAL            :: q_h,n_h,x_h,d_h,v_h
    REAL            :: q_i,n_i,x_i,d_i,v_i
    REAL            :: coll_n,coll_q,e_coll
    REAL, SAVE      :: delta_n_hh,delta_n_hi,delta_n_ii
    REAL, SAVE      :: delta_q_hh,delta_q_hi,delta_q_ii
    REAL, SAVE      :: theta_n_hh,theta_n_hi,theta_n_ii
    REAL, SAVE      :: theta_q_hh,theta_q_hi,theta_q_ii

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,ice,0)
      delta_n_hi = coll_delta_12(hail,ice,0)
      delta_n_ii = coll_delta_22(hail,ice,0)
      delta_q_hh = coll_delta_11(hail,ice,0) 
      delta_q_hi = coll_delta_12(hail,ice,1)
      delta_q_ii = coll_delta_22(hail,ice,1)

      theta_n_hh = coll_theta_11(hail,ice,0)
      theta_n_hi = coll_theta_12(hail,ice,0)
      theta_n_ii = coll_theta_22(hail,ice,0)
      theta_q_hh = coll_theta_11(hail,ice,0)
      theta_q_hi = coll_theta_12(hail,ice,1)
      theta_q_ii = coll_theta_22(hail,ice,1)

      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_i = q_ice(i,j,k)                                  !..Fluessigwassergehalt in SI
          q_h = q_hail(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (q_i > q_krit .AND. q_h > q_krit) THEN
            T_a = T_0(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
              e_coll = MIN(EXP(0.09*(T_a-T_3)),1.0)
            ENDIF

            n_i = n_ice(i,j,k)                                     !..Anzahldichte
            n_h = n_hail(i,j,k)                                    !..Anzahldichte

            x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)      !..mittlere Masse
            x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max)    !..mittlere Masse

            D_i = ice%a_geo * x_i**ice%b_geo                       !..mittlerer Durchmesser
            v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)      !..mittlere Sedimentationsgeschw.

            D_h = hail%a_geo * x_h**hail%b_geo                  !..mittlerer Durchmesser
            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            coll_n = pi/4.0 * n_h * n_i * e_coll * dt & 
                 &   * (delta_n_hh * D_h**2 + delta_n_hi * D_h*D_i + delta_n_ii * D_i**2) &
                 &   * (theta_n_hh * v_h**2 - theta_n_hi * v_h*v_i + theta_n_ii * v_i**2)**0.5

            coll_q = pi/4.0 * n_h * q_i * e_coll * dt & 
                 &   * (delta_q_hh * D_h**2 + delta_q_hi * D_h*D_i + delta_q_ii * D_i**2) &
                 &   * (theta_q_hh * v_h**2 - theta_q_hi * v_h*v_i + theta_q_ii * v_i**2)**0.5

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            q_hail(i,j,k) = q_hail(i,j,k) + coll_q
            q_ice(i,j,k)  = q_ice(i,j,k)  - coll_q
            n_ice(i,j,k)  = n_ice(i,j,k)  - coll_n

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE hail_ice_collection
  ! <<hn

  SUBROUTINE snow_ice_collection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Koagulation von Schnee und Eis                          *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: T_a
    REAL            :: q_g,n_g,x_g,d_g,v_g
    REAL            :: q_i,n_i,x_i,d_i,v_i
    REAL            :: coll_n,coll_q,e_coll
    REAL, SAVE      :: delta_n_ss,delta_n_si,delta_n_ii
    REAL, SAVE      :: delta_q_ss,delta_q_si,delta_q_ii
    REAL, SAVE      :: theta_n_ss,theta_n_si,theta_n_ii
    REAL, SAVE      :: theta_q_ss,theta_q_si,theta_q_ii

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      delta_n_ss = coll_delta_11(snow,ice,0)
      delta_n_si = coll_delta_12(snow,ice,0)
      delta_n_ii = coll_delta_22(snow,ice,0)
      delta_q_ss = coll_delta_11(snow,ice,0) 
      delta_q_si = coll_delta_12(snow,ice,1)
      delta_q_ii = coll_delta_22(snow,ice,1)

      theta_n_ss = coll_theta_11(snow,ice,0)
      theta_n_si = coll_theta_12(snow,ice,0)
      theta_n_ii = coll_theta_22(snow,ice,0)
      theta_q_ss = coll_theta_11(snow,ice,0)
      theta_q_si = coll_theta_12(snow,ice,1)
      theta_q_ii = coll_theta_22(snow,ice,1)

      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_i = q_ice(i,j,k)
          q_g = q_snow(i,j,k)
          IF (q_i > q_krit .AND. q_g > q_krit) THEN
            T_a = T_0(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
              e_coll = MAX(0.1,MIN(EXP(0.09*(T_a-T_3)),1.0))
              !.. Temperaturabhaengige Efficiency nach Cotton et al. (1986) 
              !   (siehe auch Straka, 1989; S. 53)
              !e_coll = MIN(10**(0.035*(T_a-T_3)-0.7),0.2)
            ENDIF

            n_i = n_ice(i,j,k)                       !..Anzahldichte in SI
            n_g = n_snow(i,j,k)                      !..Anzahldichte in SI

            x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)     !..mittlere Masse in SI     
            x_g = MIN(MAX(q_g/(n_g+eps),snow%x_min),snow%x_max)   !..mittlere Masse in SI     

            D_i = ice%a_geo * x_i**ice%b_geo                      !..mittlerer Durchmesser
            v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)     !..mittlere Sedimentationsgeschw.

            D_g = snow%a_geo * x_g**snow%b_geo                    !..mittlerer Durchmesser
            v_g = snow%a_vel * x_g**snow%b_vel * rrho_04(i,j,k)   !..mittlere Sedimentationsgeschw.

            coll_n = pi/4.0 * n_g * n_i * e_coll * dt & 
                 &   * (delta_n_ss * D_g**2 + delta_n_si * D_g*D_i + delta_n_ii * D_i**2) &
                 &   * (theta_n_ss * v_g**2 - theta_n_si * v_g*v_i + theta_n_ii * v_i**2  &
                 &     +snow_s_vel**2 + ice_s_vel**2)**0.5

            coll_q = pi/4.0 * n_g * q_i * e_coll * dt & 
                 &   * (delta_q_ss * D_g**2 + delta_q_si * D_g*D_i + delta_q_ii * D_i**2) &
                 &   * (theta_q_ss * v_g**2 - theta_q_si * v_g*v_i + theta_q_ii * v_i**2  &
                 &     +snow_s_vel**2 + ice_s_vel**2)**0.5

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            q_snow(i,j,k) = q_snow(i,j,k) + coll_q
            q_ice(i,j,k)  = q_ice(i,j,k)  - coll_q
            n_ice(i,j,k)  = n_ice(i,j,k)  - coll_n

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE snow_ice_collection

  SUBROUTINE graupel_selfcollection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Selbsteinfangs der Eispartikel                          *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0

    REAL            :: q_g,n_g,x_g,d_g,v_g
    REAL            :: self_n
    REAL,SAVE       :: delta_n_11,delta_n_12
    REAL,SAVE       :: theta_n_11,theta_n_12
    REAL,SAVE       :: delta_n
    REAL,SAVE       :: theta_n
    REAL,SAVE       :: coll_n

    REAL, PARAMETER :: eps  = 1.e-20
    REAL, PARAMETER :: e_gg = 0.1
    REAL, PARAMETER :: e_gg_wet = 0.4

    IF (firstcall.NE.1) THEN
      delta_n_11 = coll_delta_11(graupel,graupel,0)
      delta_n_12 = coll_delta_12(graupel,graupel,0)

      theta_n_11 = coll_theta_11(graupel,graupel,0)
      theta_n_12 = coll_theta_12(graupel,graupel,0)

      delta_n = (2.0*delta_n_11 + delta_n_12)
      theta_n = (2.0*theta_n_11 - theta_n_12)**0.5

      coll_n  = delta_n * theta_n

      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_g = q_graupel(i,j,k)
          n_g = n_graupel(i,j,k)                                       !..Anzahldichte in SI
          IF ( n_g > 0.0 ) THEN
            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse in SI
            D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.
            !self_n = pi/4.0 * e_gg * coll_n * n_g**2 * D_g**2 * v_g * dt
            ! hn>> Korrektur:
            ! Faktor 1/2 ergaenzt
            ! ub>> efficiency von nassem Graupel etwas erhoeht:
            IF (T_0(i,j,k) > T_3) THEN
              self_n = pi/8.0 * e_gg_wet * coll_n * n_g**2 * D_g**2 * v_g * dt
            ELSE
              self_n = pi/8.0 * e_gg * coll_n * n_g**2 * D_g**2 * v_g * dt
            END IF
            ! <<hn

            self_n = MIN(self_n,n_g)

            n_graupel(i,j,k) = n_graupel(i,j,k) - self_n

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_selfcollection

  SUBROUTINE ice_melting ()
    !*******************************************************************************
    !                                                                              *
    ! Diese Routine behandelt das Schmelzen der Eispartikel                        *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    REAL            :: T_a             !..absolute Temperatur
    REAL            :: q_i,x_i,n_i,weight
    REAL            :: melt_q, melt_n
    REAL, PARAMETER :: eps = 1.00e-20

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_i = q_ice(i,j,k)
          T_a = T_0(i,j,k)
          IF (T_a > T_3 .AND. q_i > 0.0) THEN            !..Nur Fluessigwasser
            n_i = n_ice(i,j,k)
            x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse in SI     

            melt_q = q_i                !  spontanes Schmelzen oberhalb des Tripelpunkts
            melt_n = n_i
            IF (x_i > cloud%x_max) THEN
              weight = 0.0
            ELSE
              weight = 1.0
            ENDIF

            q_ice(i, j, k)   = q_ice(i, j, k) - melt_q
            n_ice(i, j, k)   = n_ice(i, j, k) - melt_n

            q_cloud(i,j,k) = q_cloud(i,j,k) + melt_q * weight
            n_cloud(i,j,k) = n_cloud(i,j,k) + melt_n * weight
            q_rain(i,j,k)  = q_rain(i,j,k)  + melt_q * (1.0 - weight)
            n_rain(i,j,k)  = n_rain(i,j,k)  + melt_n * (1.0 - weight)
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE ice_melting

  SUBROUTINE rain_evaporation ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Verdunstungung von Regentropfen
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0
    REAL            :: T_a,p_a,e_sw,s_sw,g_d,eva,eva_q,eva_n
    REAL            :: q_d,q_r,n_r,x_r,d_r,v_r,f_v,N_re,e_d,f_q
    REAL            :: mue,d_m,gamma_eva,lam,n0r,d_vtp,gfak
    REAL, SAVE      :: c_r               !..Koeff. fuer mittlere Kapazitaet
    REAL, SAVE      :: a_q,b_q   !..Koeff. fuer mittleren Ventilationkoeff.

    LOGICAL, PARAMETER :: use_mu_Dm_rain_evap = .TRUE.
    REAL, PARAMETER :: & !..Axel's mu-Dm-relation for raindrops based on 1d-bin model
      rain_cmu0 = 6.0, &     ! Axel's 2007 relation
      rain_cmu1 = 30.0, &    !
      rain_cmu2 = 1.00e+3, & !
      rain_cmu3 = 1.10e-3, & ! D_eq,break
      rain_cmu4 = 1.0, &     !
      rain_gfak = 1.0
    INTEGER, PARAMETER :: &
      rain_cmu5 = 2          ! exponent

    REAL, PARAMETER :: aa = 9.65e+00     ! in SI [m/s]
    REAL, PARAMETER :: bb = 1.03e+01     ! in SI [m/s]
    REAL, PARAMETER :: cc = 6.00e+02     ! in SI [1/m]

    REAL, PARAMETER :: eps = 1.e-20

    IF (firstcall.NE.1) THEN
      a_q = vent_coeff_a(rain,1)
      b_q = vent_coeff_b(rain,1)
      c_r = 1.0 / rain%cap
      firstcall = 1
    END IF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)
          n_r  = n_rain(i,j,k)
          q_d  = q(i,j,k)
          p_a  = p_0(i,j,k)
          T_a  = T_0(i,j,k)
          e_d  = q(i,j,k) * R_d * T_a
          e_sw = esl(T_a)
          s_sw = e_d / e_sw - 1.0  !..Uebersaettigung bzgl. Wasser

          IF(s_sw < 0.0 .AND. q_r > 0.0 .AND. q_cloud(i,j,k) < q_krit)THEN

            D_vtp = 8.7602e-5 * T_a**(1.81) / p_a
            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_sw) )

            IF (use_mu_Dm_rain_evap) THEN
              x_r = q_r / (n_r+eps)
              x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)

              D_m = rain%a_geo * x_r**rain%b_geo                 

              IF (D_m.LE.rain_cmu3) THEN ! see Seifert (2008)            
                mue = rain_cmu0*TANH((4.*rain_cmu2*(D_m-rain_cmu3))**rain_cmu5) &
                  & + rain_cmu4
              ELSE
                mue = rain_cmu1*TANH((1.*rain_cmu2*(D_m-rain_cmu3))**rain_cmu5) &
                  & + rain_cmu4
              ENDIF

              lam = (pi/6.*rho_w &
                &      * (mue+3.0)*(mue+2.0)*(mue+1.0) / x_r)**(1./3.)

!AS20081207 optimierte Variante

              ! chebyshev approximation of gamma_fct(zrmue+5/2)/gamma_fct(zrmue+2)
              ! (for mue in [0,31], error smaller than 2%)
!             gfak =  0.1357940435E+01              &
!                  &   + mue * 0.3033273220E+00 * ( &
!                  &   - mue * 0.1299313363E-01 * ( &
!                  &   + mue * 0.4002257774E-03 * ( &
!                  &   - mue * 0.4856703981E-05 ) ) )

! Uli              
              gfak =  0.1357940435E+01 &
                &  + mue * ( +0.3033273220E+00  &
                &  + mue * ( -0.1299313363E-01  &
                &  + mue * ( +0.4002257774E-03  &
                &  - mue * 0.4856703981E-05 ) ) )

!              gfak = gfct(mue+5/2)/gfct(mue+2)

!               gfak = 0.1357940435E+01          &
!                  & + 0.3033273220E+00 * mue    &
!                  & - 0.1299313363E-01 * mue**2 &
!                  & + 0.4002257774E-03 * mue**3 & 
!                  & - 0.4856703981E-05 * mue**4

              f_q  = rain%a_ven                                          &
              &    + rain%b_ven * N_sc**n_f                              &
              &                 * (aa/nu_l*rrho_04(i,j,k))**m_f          &
              &    * gfak / SQRT(lam)                                    &
              &    * (1.0                                                &
              &      - 1./2.  * (bb/aa)**1 * (lam/(1.*cc+lam))**(mue+5.0/2.0) &
              &      - 1./8.  * (bb/aa)**2 * (lam/(2.*cc+lam))**(mue+5.0/2.0) &
              &      - 1./16. * (bb/aa)**3 * (lam/(3.*cc+lam))**(mue+5.0/2.0) &
              &      - 5./127.* (bb/aa)**4 * (lam/(4.*cc+lam))**(mue+5.0/2.0) &
              &      )

              
              IF (rain_gfak.GT.0) THEN
                gamma_eva = rain_gfak * (1.1e-3/D_m) * EXP(-0.2*mue)
              ELSE
                gamma_eva = 1.0
              END IF

              eva_q = g_d * c_r * n_r * (mue+1.0) / lam * f_q * s_sw * dt
              eva_n = gamma_eva * eva_q / x_r

            ELSE ! old evaporation

              IF (cloud_typ <= 1) THEN
                n0r = 1.00e+07
                n_r = n0r * (pi*rho_w*n0r/q_r)**(-1./4.)     !..Nur q_r  (n0 = const) 
              ENDIF
              x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
              x_r = MAX(1e3/2e7*q_r,MIN(x_r,1e4/2e5*q_r))

              d_r = rain%a_geo * x_r**rain%b_geo                 
              v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)
              
              N_re = v_r * d_r / nu_l                            
              f_v  = N_sc**n_f * N_re**m_f           
              
              f_q  = a_q + b_q * f_v
              
              eva   = g_d * n_r * c_r * d_r * s_sw * dt 
              eva_q = f_q * eva
              !eva_n = f_q * eva / x_r
              !eva_n = MAX(MIN( (eva_q + q_r) / x_r - n_r, 0.0), -n_r)
              eva_n = MAX(MIN( eva_q / x_r, 0.0), -n_r)
            ENDIF

            eva_q = MAX(-eva_q,0.) 
            eva_n = MAX(-eva_n,0.) 

            eva_q = MIN(eva_q,q_r) 
            eva_n = MIN(eva_n,n_r) 

            !..Berechnung der H_2O-Komponenten
            q(i,j,k)      = q(i,j,k)      + eva_q
            q_rain(i,j,k) = q_rain(i,j,k) - eva_q
            IF (cloud_typ > 1) THEN
              n_rain(i,j,k) = n_rain(i,j,k) - eva_n
            END IF

          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE rain_evaporation

  SUBROUTINE graupel_evaporation ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Verdunstung von schmelzendem Graupel                        *
    !                                                                              *
    !*******************************************************************************

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0
    REAL            :: T_a,e_sw,s_sw,g_d,eva
    REAL            :: q_d,q_g,n_g,x_g,d_g,v_g,f_v,N_re,e_d
    REAL, SAVE      :: c_g             !..Koeff. fuer mittlere Kapazitaet
    REAL, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      a_f = vent_coeff_a(graupel,1)
      b_f = vent_coeff_b(graupel,1)
      c_g = 1.0 / graupel%cap
      firstcall = 1
    END IF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_g = q_graupel(i,j,k)
          T_a  = T_0(i,j,k)
          IF(q_g > 0. .AND. T_a > T_3)THEN
            q_d = q(i,j,k)
            e_d  = q(i,j,k) * R_d * T_a 
            e_sw = esl(T_a)
            s_sw = e_d / e_sw - 1.0                   !..Uebersaettigung bzgl. Wasser

            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_g = n_graupel(i,j,k)                          !..Anzahldichte in SI
            q_g = q_graupel(i,j,k)                          !..Fluessigwassergehalt in SI

            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max)   !..mittlere Masse in SI

            d_g = graupel%a_geo * x_g**graupel%b_geo                   !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

            N_re = v_g * d_g / nu_l                         !..mittlere Reynoldszahl
            f_v  = a_f + b_f * N_sc**n_f * N_re**m_f        !..mittlerer Vent.Koeff.

            eva = g_d * n_g * c_g * d_g * f_v * s_sw * dt

            eva = MAX(-eva,0.) 
            eva = MIN(eva,q_g) 

            !..Berechnung der H_2O-Komponenten
            q(i,j,k)         = q(i,j,k)         + eva
            q_graupel(i,j,k) = q_graupel(i,j,k) - eva
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE graupel_evaporation

  SUBROUTINE hail_evaporation ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Verdunstung von schmelzendem Hagel                          *
    !                                                                              *
    !*******************************************************************************

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0
    REAL            :: T_a,e_sw,s_sw,g_d,eva
    REAL            :: q_d,q_h,n_h,x_h,d_h,v_h,f_v,N_re,e_d
    REAL, SAVE      :: c_h             !..Koeff. fuer mittlere Kapazitaet
    REAL, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      a_f = vent_coeff_a(hail,1)
      b_f = vent_coeff_b(hail,1)
      c_h = 1.0 / hail%cap
      firstcall = 1
    END IF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_h = q_hail(i,j,k)
          T_a  = T_0(i,j,k)
          IF(q_h > 0. .AND. T_a > T_3)THEN
            q_d = q(i,j,k)
            e_d  = q(i,j,k) * R_d * T_a 
            e_sw = esl(T_a)
            s_sw = e_d / e_sw - 1.0                   !..Uebersaettigung bzgl. Wasser

            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_h = n_hail(i,j,k)                          !..Anzahldichte in SI
            q_h = q_hail(i,j,k)                          !..Fluessigwassergehalt in SI

            x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max)   !..mittlere Masse in SI

            d_h = hail%a_geo * x_h**hail%b_geo                   !..mittlerer Durchmesser
            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

            N_re = v_h * d_h / nu_l                         !..mittlere Reynoldszahl
            f_v  = a_f + b_f * N_sc**n_f * N_re**m_f        !..mittlerer Vent.Koeff.

            eva = g_d * n_h * c_h * d_h * f_v * s_sw * dt

            eva = MAX(-eva,0.) 
            eva = MIN(eva,q_h) 

            !..Berechnung der H_2O-Komponenten
            q(i,j,k)      = q(i,j,k)      + eva
            q_hail(i,j,k) = q_hail(i,j,k) - eva
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE hail_evaporation

  SUBROUTINE snow_evaporation ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Verdunstung von schmelzendem Schnee                         *
    !                                                                              *
    !*******************************************************************************

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall = 0
    REAL            :: T_a,e_sw,s_sw,g_d,eva
    REAL            :: q_d,q_s,n_s,x_s,d_s,v_s,f_v,N_re,e_d
    REAL, SAVE      :: c_s             !..Koeff. fuer mittlere Kapazitaet
    REAL, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.

    REAL, PARAMETER :: eps  = 1.e-20

    IF (firstcall.NE.1) THEN
      a_f = vent_coeff_a(snow,1)
      b_f = vent_coeff_b(snow,1)
      c_s = 1.0/snow%cap
      firstcall = 1
    END IF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_s = q_snow(i,j,k)
          T_a  = T_0(i,j,k)
          IF(q_s > 0. .AND. T_a > T_3)THEN
            q_d = q(i,j,k)
            e_d  = q(i,j,k) * R_d * T_a 
            e_sw = esl(T_a)
            s_sw = e_d / e_sw - 1.0                   !..Uebersaettigung bzgl. Wasser


            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_s = n_snow(i,j,k)                             !..Anzahldichte in SI
            q_s = q_snow(i,j,k)                             !..Fluessigwassergehalt in SI

            x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse in SI     

            d_s = snow%a_geo * x_s**snow%b_geo                     !..mittlerer Durchmesser
            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k)    !..mittlere Sedimentationsgeschw.

            N_re = v_s * d_s / nu_l                          !..mittlere Reynoldszahl
            f_v  = a_f + b_f * N_sc**n_f * N_re**m_f   !..mittlerer Vent.Koeff.

            eva = g_d * n_s * c_s * d_s * f_v * s_sw * dt

            eva = MAX(-eva,0.) 
            eva = MIN(eva,q_s) 

            !..Berechnung der H_2O-Komponenten
            q(i,j,k)         = q(i,j,k)         + eva
            q_snow(i,j,k) = q_snow(i,j,k) - eva

          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE snow_evaporation

!END MODULE wolken_eis

!==============================================================================

!MODULE wolken


  SUBROUTINE clouds ()
    !*******************************************************************************

    rrho_04 = sqrt(rho0/rho_0)
    rrho_c  = rho0/rho_0

    IF (sflg) CALL sb_var_stat_reset() ! Reset

    WHERE (q         < 0.0) q         = 0.0
    WHERE (n_ice     < 0.0) n_ice     = 0.0
    WHERE (q_ice     < 0.0) q_ice     = 0.0
    WHERE (n_rain    < 0.0) n_rain    = 0.0
    WHERE (q_rain    < 0.0) q_rain    = 0.0
    WHERE (n_snow    < 0.0) n_snow    = 0.0
    WHERE (q_snow    < 0.0) q_snow    = 0.0
    WHERE (n_cloud   < 0.0) n_cloud   = 0.0
    WHERE (q_cloud   < 0.0) q_cloud   = 0.0
    WHERE (q_cloud   < 0.0) n_cloud   = 0.0
    WHERE (n_graupel < 0.0) n_graupel = 0.0
    WHERE (q_graupel < 0.0) q_graupel = 0.0
    WHERE (n_hail    < 0.0) n_hail    = 0.0
    WHERE (q_hail    < 0.0) q_hail    = 0.0

    IF (sflg) CALL sb_var_stat('diag') ! Diagnostics

    IF (cloud_typ > 1 ) THEN
      IF (nuc_c_typ .NE. 0) THEN
        WRITE(*,*) 'Fixed CDNC required!'
        STOP
      END IF
    END IF

    IF (nuc_c_typ.ne.0) THEN
       n_cloud = MAX(n_cloud, q_cloud / cloud%x_max)
       n_cloud = MIN(n_cloud, q_cloud / cloud%x_min)
    END IF

    IF (sflg) CALL sb_var_stat('nucl') ! Cloud nucleation

    IF (ice_typ .NE. 0) THEN

      ! Eisnukleation
      CALL ice_nucleation()
      n_ice     = MIN(n_ice, q_ice/ice%x_min)
      n_ice     = MAX(n_ice, q_ice/ice%x_max)
      IF (sflg) CALL sb_var_stat('nucl') ! Ice nucleation

      ! Gefrieren der Wolkentropfen:
      IF (drop_freeze) CALL cloud_freeze ()
      IF (sflg) CALL sb_var_stat('frez') ! Cloud freezing

      CALL vapor_dep_relaxation(dt)
      IF (sflg) CALL sb_var_stat('cond') ! Condensation

      CALL ice_selfcollection ()

      CALL snow_selfcollection ()

      CALL snow_ice_collection ()

      CALL graupel_selfcollection ()

      CALL graupel_ice_collection ()

      CALL graupel_snow_collection ()

      IF (ice_typ > 1) THEN

        ! Disabled (31.3.2023): CALL graupel_hail_conv_wet ()

        CALL hail_ice_collection ()

        CALL hail_snow_collection ()

      ENDIF

      IF (.NOT. use_ice_graupel_conv_uli) THEN 
        ! Old SB2006 scheme

        CALL ice_cloud_riming ()
        CALL snow_cloud_riming ()
        CALL graupel_cloud_riming ()
        IF (ice_typ > 1) CALL hail_cloud_riming ()
        
        ! Bereifen mit Regentropfen
        CALL ice_rain_riming ()
        CALL snow_rain_riming ()
        CALL graupel_rain_riming ()
        IF (ice_typ > 1) CALL hail_rain_riming ()

      ELSE
        ! new scheme

        ! Bereifen mit Wolkentropfen
        CALL ice_cloud_riming ()
        CALL snow_cloud_riming ()
        ! Bereifen mit Regentropfen
        CALL ice_rain_riming ()
        CALL snow_rain_riming ()
        ! an dieser Stelle werden nun die vorher lediglich gespeicherten 
        ! riming-Raten summiert:
        CALL complete_ice_snow_riming ()

        CALL graupel_cloud_riming ()
        IF (ice_typ > 1) CALL hail_cloud_riming ()
        CALL graupel_rain_riming ()
        IF (ice_typ > 1) CALL hail_rain_riming ()

      ENDIF

      IF (sflg) CALL sb_var_stat('coag') ! Collisions

      ! Gefrieren der Regentropfen:

      IF (drop_freeze) THEN
        CALL rain_freeze ()
        ! simpler SB2006 rain-to-graupel freezing scheme
        !CALL rain_freeze_old ()
      END IF
      IF (sflg) CALL sb_var_stat('frez') ! Rain freezing

      ! Schmelzen der Eispartikel
      CALL ice_melting ()
      CALL snow_melting ()
      CALL graupel_melting ()
      IF (ice_typ > 1) CALL hail_melting ()
      IF (sflg) CALL sb_var_stat('melt') ! Melting

      ! Verdunstung von schmelzenden Eispartikeln
      CALL graupel_evaporation ()
      IF (ice_typ > 1) CALL hail_evaporation ()
      CALL snow_evaporation ()
      IF (sflg) CALL sb_var_stat('cond') ! Condensation

      ! Groesse der Partikel beschraenken
      n_snow    = MAX(n_snow, q_snow / snow%x_max)
      n_graupel = MAX(n_graupel, q_graupel / graupel%x_max)
      n_hail    = MAX(n_hail, q_hail / hail%x_max)
      n_ice     = MIN(n_ice, q_ice/ice%x_min)
      n_ice     = MAX(n_ice, q_ice/ice%x_max)
      IF (sflg) CALL sb_var_stat('diag') ! Diagnostics

      ! Limit ice number conc. to 1000 per liter in any case:
      ! n_ice = MIN(n_ice, 1000e3)

    ENDIF

    ! Koagulation der Wolken- und Regentropfen
    IF (cloud_typ == 0) THEN
      !CALL autoconversionKS ()   ! Kessler (1-moment-scheme)
      !CALL accretionKS ()
    ELSE IF (cloud_typ == 1) THEN 
      CALL autoconversionSB ()   ! Seifert and Beheng (2000) (1-moment-scheme)
      CALL accretionSB ()
    ELSE IF (cloud_typ == 2) THEN 
      CALL autoconversionKS ()   ! Kessler (1969) (2-moment-scheme)
      CALL accretionKS ()
      CALL rain_selfcollectionSB ()
    ELSE IF (cloud_typ == 3 .OR. &
         &   cloud_typ == 6 .OR. &
         &   cloud_typ == 7 .OR. &
         &   cloud_typ == 8) THEN
      CALL autoconversionSB ()   ! Seifert and Beheng (2000) (2-moment-scheme)
      IF (sflg) CALL sb_var_stat('auto') ! Autoconversion
      CALL accretionSB ()
      CALL rain_selfcollectionSB ()
      IF (sflg) CALL sb_var_stat('coag') ! Collisions
    ELSE IF (cloud_typ == 4) THEN
      CALL autoconversionKK ()   ! Khai.. and Kogan (2000)
      CALL accretionKK ()
      CALL rain_selfcollectionSB ()
    ELSE IF (cloud_typ == 5) THEN
      CALL autoconversionKB ()   ! Beheng (1994)
      CALL accretionKB ()
      CALL rain_selfcollectionSB ()
    ENDIF

    ! Verdunstung von Regentropfen
    CALL rain_evaporation ()
    IF (sflg) CALL sb_var_stat('cond') ! Evaporation from rain

    IF (nuc_c_typ > 0) THEN
       n_cloud = MIN(n_cloud, q_cloud/cloud%x_min)
       n_cloud = MAX(n_cloud, q_cloud/cloud%x_max)

       ! Hard upper limit for cloud number conc.
       n_cloud = MIN(n_cloud, 5000e6)
    end if

    ! Hydrometeore in ihrer Groesse beschraenken:
    n_rain = MIN(n_rain, q_rain/rain%x_min)
    n_rain = MAX(n_rain, q_rain/rain%x_max)
    IF (ice_typ > 0) THEN
       n_ice = MIN(n_ice, q_ice/ice%x_min)
       n_ice = MAX(n_ice, q_ice/ice%x_max)
       n_snow = MIN(n_snow, q_snow/snow%x_min)
       n_snow = MAX(n_snow, q_snow/snow%x_max)
       n_graupel = MIN(n_graupel, q_graupel/graupel%x_min)
       n_graupel = MAX(n_graupel, q_graupel/graupel%x_max)
    END IF
    IF (ice_typ > 1) THEN
      n_hail = MIN(n_hail, q_hail/hail%x_min)
      n_hail = MAX(n_hail, q_hail/hail%x_max)
    END IF

    if (.false.) then
       WHERE (q         < 0.0) q         = 0.0
       WHERE (n_ice     < 0.0) n_ice     = 0.0
       WHERE (q_ice     < 0.0) q_ice     = 0.0
       WHERE (n_rain    < 0.0) n_rain    = 0.0
       WHERE (q_rain    < 0.0) q_rain    = 0.0
       WHERE (n_snow    < 0.0) n_snow    = 0.0
       WHERE (q_snow    < 0.0) q_snow    = 0.0
       WHERE (n_cloud   < 0.0) n_cloud   = 0.0
       WHERE (q_cloud   < 0.0) q_cloud   = 0.0
       WHERE (n_cloud   < 0.0) n_cloud   = 0.0
       WHERE (n_graupel < 0.0) n_graupel = 0.0
       WHERE (q_graupel < 0.0) q_graupel = 0.0
       WHERE (n_hail    < 0.0) n_hail    = 0.0
       WHERE (q_hail    < 0.0) q_hail    = 0.0
    end if

    IF (sflg) CALL sb_var_stat('diag') ! Diagnostics

  CONTAINS

      ! Collecting statistical outputs from the Seifert and Beheng micophysics
      !     Main program: out_mcrp_nout, out_mcrp_data, out_mcrp_list
      ! 1) Save the current state
      SUBROUTINE sb_var_stat_reset()
        tmp_rv = q ! Water vapor
        tmp_rc = q_cloud ! Cloud water mixing ratio
        tmp_nr = n_rain; tmp_rr = q_rain ! Rain drop number and mixing ratio
        tmp_ni = n_ice; tmp_ri = q_ice ! Ice
        tmp_ns = n_snow; tmp_rs = q_snow ! Snow
        tmp_ng = n_graupel; tmp_rg = q_graupel ! Graupel
        tmp_nh = n_hail; tmp_rh = q_hail ! Hail
      END SUBROUTINE sb_var_stat_reset
      ! 2) Collect statistics
      SUBROUTINE sb_var_stat(prefix) !imicro)
        IMPLICIT NONE
        character(len=4), INTENT(IN) :: prefix ! Process name
        ! Local
        INTEGER :: k
        !
        IF (out_mcrp_nout==0) RETURN
        !
        ! Generate the requested ouputs
        DO k=1,out_mcrp_nout
            IF ( prefix//'_rv' == out_mcrp_list(k) ) THEN
                ! Water vapor (diagnostic)
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (q - tmp_rv)/dt
            ELSEIF ( prefix//'_rc' == out_mcrp_list(k) ) THEN
                ! Cloud water (diagnostic)
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (q_cloud - tmp_rc)/dt
            ELSEIF ( prefix//'_nr' == out_mcrp_list(k) ) THEN
                ! Rain number
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (n_rain - tmp_nr)/dt
            ELSEIF ( prefix//'_rr' == out_mcrp_list(k) ) THEN
                ! Rain mixing ratio
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (q_rain - tmp_rr)/dt
            ELSEIF ( prefix//'_ni' == out_mcrp_list(k) ) THEN
                ! Ice number
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (n_ice - tmp_ni)/dt
            ELSEIF ( prefix//'_ri' == out_mcrp_list(k) ) THEN
                ! Ice mixing ratio
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (q_ice - tmp_ri)/dt
            ELSEIF ( prefix//'_ns' == out_mcrp_list(k) ) THEN
                ! Snow numbero
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (n_snow - tmp_ns)/dt
            ELSEIF ( prefix//'_rs' == out_mcrp_list(k) ) THEN
                ! Snow mixing ratio
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (q_snow - tmp_rs)/dt
            ELSEIF ( prefix//'_ng' == out_mcrp_list(k) ) THEN
                ! Graupel number
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (n_graupel- tmp_ng)/dt
            ELSEIF ( prefix//'_rg' == out_mcrp_list(k) ) THEN
                ! Graupel mixing ratio
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (q_graupel- tmp_rg)/dt
            ELSEIF ( prefix//'_nh' == out_mcrp_list(k) ) THEN
                ! Hail number
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (n_hail- tmp_nh)/dt
            ELSEIF ( prefix//'_hg' == out_mcrp_list(k) ) THEN
                ! Hail mixing ratio
                out_mcrp_data(:,:,:,k) = out_mcrp_data(:,:,:,k) + (q_hail- tmp_rh)/dt
            ENDIF
        ENDDO
        !
        ! .. and reset
        CALL sb_var_stat_reset()
        !
      END SUBROUTINE sb_var_stat

    SUBROUTINE autoconversionKS ()
      !*******************************************************************************
      !                                                                              *
      !   Berechnung der Autokonversion von Wolkenwasser zu Niederschlagswasser      *
      !                                                                              *
      !*******************************************************************************

      IMPLICIT NONE

      !..Parameter fuer Kessler-Ansatz
      REAL, PARAMETER :: C_au = 1.00e-3
      REAL, PARAMETER :: C_qc = 5.00e-4
      REAL, PARAMETER :: x_s  = 2.60e-10     ! Trennmasse Wolken-Regentropfen

      !..Locale Variablen
      INTEGER          :: i, j, k
      REAL :: au, dt_cau, q_c

      !..Skalar-Initialisierung
      dt_cau = dt * C_au

      !..Autoconversionrate nach Kessler
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            q_c = q_cloud(i,j,k)
            IF (q_c >= C_qc) THEN

              au = dt_cau * MAX(q_c - C_qc,0.)
              au = MIN(q_c,au)

              !..Berechnung der H_2O-Komponenten
              n_rain(i,j,k)  = n_rain(i,j,k)  + au / x_s
              q_rain(i,j,k)  = q_rain(i,j,k)  + au
              n_cloud(i,j,k) = n_cloud(i,j,k) - au / x_s * 2.0
              q_cloud(i,j,k) = q_cloud(i,j,k) - au
            ENDIF
          END DO
        END DO
      END DO

    END SUBROUTINE autoconversionKS

    SUBROUTINE accretionKS ()
      !*******************************************************************************
      !                                                                              *
      !       Berechnung der Akkretion von Wolkentroepfchen an Regentroepfchen       *
      !                                                                              *
      !*******************************************************************************

      IMPLICIT NONE

      !..Parameter fuer Kessler-Ansatz
      REAL, PARAMETER :: C_ac  = 0.94
      REAL, PARAMETER :: eps  = 1.00e-25

      !..Lokale Variablen
      INTEGER          :: i, j, k
      REAL :: ac,n_c,L_c,x_c,L_r

      !..Akkreszenzrate nach Kessler (siehe Dotzek, 1999, p. 39)
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            L_c = q_cloud(i,j,k)
            N_c = n_cloud(i,j,k)
            L_r = q_rain(i,j,k)
            IF (L_c > 0. .AND. L_r > 0.) THEN

              ac = C_ac * L_c * L_r**0.875 * rrho_04(i,j,k) * dt

              ac = MIN(L_c,ac)

              x_c = MIN(MAX(L_c/(N_c+eps),cloud%x_min),cloud%x_max) !..mittlere Masse in SI
              q_rain(i,j,k)  = q_rain(i,j,k)  + ac
              q_cloud(i,j,k) = q_cloud(i,j,k) - ac
              n_cloud(i,j,k) = n_cloud(i,j,k) - ac / x_c
            ENDIF
          END DO
        END DO
      END DO

    END SUBROUTINE accretionKS

  END SUBROUTINE clouds

  SUBROUTINE alloc_wolken ()
    !*******************************************************************************
    !           Allokierung der Wolken- und Niederschlagsfelder                    *
    !*******************************************************************************

    IMPLICIT NONE

    ALLOCATE (n_cloud(0:loc_ix, 1:loc_iy, 1:loc_iz), &
         q_cloud(0:loc_ix, 1:loc_iy, 1:loc_iz), &
         n_rain(0:loc_ix, 1:loc_iy, 1:loc_iz),  &
         q_rain(0:loc_ix, 1:loc_iy, 1:loc_iz),  &
         n_ice(0:loc_ix, 1:loc_iy, 1:loc_iz),   &
         q_ice(0:loc_ix, 1:loc_iy, 1:loc_iz),   &             
         n_snow(0:loc_ix, 1:loc_iy, 1:loc_iz),   &
         q_snow(0:loc_ix, 1:loc_iy, 1:loc_iz),   &             
         n_hail(0:loc_ix, 1:loc_iy, 1:loc_iz),   &  ! <hn
         q_hail(0:loc_ix, 1:loc_iy, 1:loc_iz),   &  ! <hn
         n_graupel(0:loc_ix, 1:loc_iy, 1:loc_iz),   &
         q_graupel(0:loc_ix, 1:loc_iy, 1:loc_iz))

    q_cloud   = 0.0
    q_rain    = 0.0
    q_ice     = 0.0
    q_snow    = 0.0
    q_graupel = 0.0
    q_hail    = 0.0   ! <hn
    n_cloud   = 0.0
    n_rain    = 0.0
    n_ice     = 0.0
    n_snow    = 0.0
    n_hail    = 0.0   ! <hn
    n_graupel = 0.0

  END SUBROUTINE alloc_wolken

  SUBROUTINE dealloc_wolken ()
    !****************************************************************************
    !           Allokierung der Wolken- und Niederschlagsfelder                 *
    !****************************************************************************
    IMPLICIT NONE

    DEALLOCATE (n_cloud,q_cloud,n_rain,q_rain,n_ice,q_ice, n_snow, &
         &  q_snow,n_graupel,q_graupel, q_hail, n_hail)

  END SUBROUTINE dealloc_wolken

  SUBROUTINE autoconversionSB ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Autokonversion und Selfcollection der Wolkentropfen         *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    REAL, PARAMETER :: eps  = 1.00e-25

    !..Locale Variablen
    INTEGER          :: i,j,k
    REAL :: q_c, q_r, n_c, x_c, nu, mu, &
         & tau, phi, k_au, k_sc, x_s, au, sc, k_c, k_1, k_2

    !..Skalar-Initialisierung

    IF (cloud_typ <= 3) THEN
      !..Parameter fuer Seifert & Beheng (2001)
      k_c  = 9.44e+9   !..Long-Kernel
      k_1  = 6.00e+2   !..Parameter fuer Phi-Fkt.
      k_2  = 0.68e+0   !..Parameter fuer Phi-Fkt.
    ELSEIF (cloud_typ == 6) THEN
      !..Parameter fuer Pinsky et al (2000) Kernel
      k_c  = 4.44e+9   !..CC-Kernel
      k_1  = 4.00e+2   !..Parameter fuer Phi-Fkt.
      k_2  = 0.70e+0   !..Parameter fuer Phi-Fkt.
    ELSEIF (cloud_typ == 7) THEN
      !..Parameter fuer Pinsky et al (2000) Kernel
      k_c  = 10.58e+9  !..CC-Kernel
      k_1  = 4.00e+2   !..Parameter fuer Phi-Fkt.
      k_2  = 0.70e+0   !..Parameter fuer Phi-Fkt.
    ELSEIF (cloud_typ == 8) THEN
      !..Parameter fuer HUCM Kernel
      k_c  = 30.00e+9  !..CC-Kernel
      k_1  = 4.00e+2   !..Parameter fuer Phi-Fkt.
      k_2  = 0.70e+0   !..Parameter fuer Phi-Fkt.
    ENDIF

    nu    = cloud%nu
    mu    = cloud%mu
    x_s   = cloud%x_max                     !..Trennmasse

    IF (mu == 1.0) THEN 
      k_au  = k_c / (2.0e1*x_s) * (nu+2.0)*(nu+4.0)/(nu+1.0)**2
      k_sc  = k_c * (nu+2.0)/(nu+1.0)
    ELSE
      k_au = k_c / (2.0e1*x_s)                                       &
        & * ( 2.0 * gfct((nu+4.0)/mu)**1                          &
        &         * gfct((nu+2.0)/mu)**1 * gfct((nu+1.0)/mu)**2   &
        &   - 1.0 * gfct((nu+3.0)/mu)**2 * gfct((nu+1.0)/mu)**2 ) &
        &   / gfct((nu+2.0)/mu)**4
      k_sc = k_c * moment_gamma(cloud,2)
    ENDIF

    !..Parametrisierung nach Seifert und Beheng (2000)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_c = q_cloud(i,j,k) !..Fluessigwassergehalt
          IF (q_c > eps) THEN
            n_c = n_cloud(i,j,k)                                  !..Anzahldichte
            q_r = q_rain(i,j,k)                                   !..Fluessigwassergehalt
            x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max) !..mittlere Masse in SI

            !..Berechnung der Autokonversionsrate nach SB2000
            au  =  k_au * q_c**2 * x_c**2 * dt
            !au  = k_au * q_c**2 * x_c**2 * dt * rrho_c(i,j,k)
            IF (q_c > 1.0e-6) THEN
              tau = MIN(MAX(1.0-q_c/(q_c+q_r+eps),eps),0.9)
              phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
              au  = au * (1.0 + phi/(1.0 - tau)**2)
            ENDIF
            au = MAX(MIN(q_c,au),0.)

            !sc = k_sc * q_c**2 * dt 
            sc = k_sc * q_c**2 * dt * rrho_c(i,j,k)

            q_rain(i,j,k) = q_rain(i,j,k) + au
            n_rain(i,j,k) = n_rain(i,j,k)  + au / x_s
            n_cloud(i,j,k) = n_cloud(i,j,k) - MIN(n_c,sc)
            q_cloud(i,j,k) = q_cloud(i,j,k) - au

          ENDIF
        END DO
      END DO
    END DO

  END SUBROUTINE autoconversionSB

  SUBROUTINE accretionSB ()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Akkretion von Wolkentroepfchen an Regentroepfchen       *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    !..Parameter fuer Seifert & Beheng (2001)
    REAL, PARAMETER :: k_r = 5.78e+0   ! Parameter Kernel
    REAL, PARAMETER :: k_1 = 5.00e-4   ! Parameter fuer Phi-Fkt.

    !..Parameter fuer Seifert (2002), not recommended
    !REAL, PARAMETER :: k_r = 5.25e+0   ! Parameter Kernel
    !REAL, PARAMETER :: k_1 = 5.00e-5   ! Parameter fuer Phi-Fkt.
    REAL, PARAMETER :: eps = 1.00e-25

    !..Lokale Variablen
    INTEGER          :: i,j,k
    REAL :: ac
    REAL :: L_c, L_r, tau, phi, n_c, x_c

    !..Parametrisierung nach Seifert und Beheng (2001)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          n_c = n_cloud(i,j,k) !..Anzahldichte
          L_c = q_cloud(i,j,k) !..Fluessigwassergehalt
          L_r = q_rain(i,j,k)  !..Fluessigwassergehalt

          IF (L_c > 0.0.AND.L_r > 0.0) THEN

            !..Berechnung der Akkreszenzrate nach SB2001
            tau = MIN(MAX(1.0-L_c/(L_c+L_r+eps),eps),1.0)
            phi = (tau/(tau+k_1))**4
            !ac  = k_r *  L_c * L_r * phi * dt
            ac  = k_r *  L_c * L_r * phi  * rrho_04(i,j,k) * dt

            ac = MIN(L_c,ac)

            x_c = MIN(MAX(L_c/(n_c+eps),cloud%x_min),cloud%x_max)

            q_rain(i,j,k)  = q_rain(i,j,k)  + ac
            q_cloud(i,j,k) = q_cloud(i,j,k) - ac
            n_cloud(i,j,k) = n_cloud(i,j,k) - MIN(n_c,ac/x_c)
          ENDIF
        END DO
      END DO
    END DO

  END SUBROUTINE accretionSB

  SUBROUTINE rain_selfcollectionSB ()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Selfcollection von Regentropfen                         *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    !..Parameter fuer Seifert & Beheng (2001)
    !REAL, PARAMETER :: k_r = 5.78e+00

    !..Parameter fuer Seifert (2002), DO NOT USE!
    !REAL, PARAMETER :: k_sc = 2.87e+00
    !REAL, PARAMETER :: k_rr = 7.12e+00
    !REAL, PARAMETER :: k_br = 1.00e+03
    !REAL, PARAMETER :: D_br = 0.90e-03
    !REAL, PARAMETER :: kap  = 6.07e+01

    !..Parameter fuer Seifert (2007)
    REAL, PARAMETER :: D_br = 1.10e-03
    REAL, PARAMETER :: k_rr = 4.33e+00
    REAL, PARAMETER :: k_br = 1.00e+03

    REAL, PARAMETER :: eps = 1.00e-25

    !..Lokale Variablen
    INTEGER          :: i, j, k
    REAL :: sc, br
    REAL :: q_r, n_r, x_r, d_r, phi1

    !..Parametrisierung nach Seifert und Beheng

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          n_r = n_rain(i,j,k)    !..Anzahldichte
          q_r = q_rain(i,j,k)    !..Fluessigwassergehalt

          IF (q_r > 0.0) THEN
            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
            D_r = rain%a_geo * x_r**rain%b_geo
            !lam = lambda_gamma(rain,x_r)

            !..Berechnung der Selfcollection nach S2002
            !sc = k_r *  n_r * q_r * rrho_04(i,j,k) * dt
            !sc = k_sc *  n_r * q_r * rrho_04(i,j,k) * dt
            !sc = k_rr *  n_r * q_r * (1.0 + kap/lam)**(-9) * rrho_04(i,j,k) * dt
            sc = k_rr *  n_r * q_r * rrho_04(i,j,k) * dt

            !..Berechnung Breakup nach S2002
            br = 0.0
            IF (D_r.GT.0.30e-3) THEN
              phi1 = (k_br * (D_r - D_br) + 1.0)
              !phi2 = 2.0 * exp(2.3e3 * (D_r - D_br)) - 1.0
              !br = max(phi1,phi2) * sc                   
              br = phi1 * sc                   
            ENDIF
            sc = MIN(n_r,sc-br)

            n_rain(i,j,k)  = n_rain(i,j,k) - sc

            ! Untere und obere Schranke fuer d_rain bzw. n_rain
            !x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max) !..mittlere Masse in SI
            !n_rain(i,j,k) = q_rain(i,j,k) / x_r
          ENDIF
        END DO
      END DO
    END DO

  END SUBROUTINE rain_selfcollectionSB

  SUBROUTINE autoconversionKB ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Autokonversion von Wolkenwasser zu Niederschlagswasser      *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    !..Parameter fuer Beheng (1994)
    REAL, PARAMETER :: r_c  = 12.0e-6
    REAL, PARAMETER :: eps  = 1.00e-25

    !..Locale Variablen
    INTEGER          :: i,j,k
    REAL :: q_c, x_c, k_a, x_s, au, nu_r

    !..Skalar-Initialisierung
    nu_r = 9.59
    x_s  = cloud%x_max                     !..Trennmasse
    x_c  = 4./3. * pi * rho_w * r_c**3     !..Mittlere Masse der Wolkentropfen
    k_a  = 6.0e+25 * nu_r**(-1.7)

    !..Parametrisierung nach Beheng (1994)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_c = q_cloud(i,j,k)!..Fluessigwassergehalt
          IF (q_c > eps) THEN

            !..Berechnung der Autokonversionsrate nach Beheng (1994)
            au = k_a * (x_c*1e3)**(3.3) * (q_c*1e-3)**(1.4) * dt * 1e3
            au = MIN(q_c,au)

            n_rain(i,j,k)  = n_rain(i,j,k)  + au / x_s * 2.0
            q_rain(i,j,k)  = q_rain(i,j,k)  + au
            q_cloud(i,j,k) = q_cloud(i,j,k) - au

          ENDIF
        END DO
      END DO
    END DO

  END SUBROUTINE autoconversionKB

  SUBROUTINE accretionKB ()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Akkretion von Wolkentroepfchen an Regentroepfchen       *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    !..Parameter fuer Beheng (1994)
    REAL, PARAMETER :: k_r = 6.00e+00   ! Parameter Kernel

    !..Lokale Variablen
    INTEGER          :: i,j,k
    REAL :: ac
    REAL :: L_c, L_r

    !..Parametrisierung nach Seifert und Beheng (2000)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          L_c = q_cloud(i,j,k) !..Fluessigwassergehalt
          L_r = q_rain(i,j,k)  !..Fluessigwassergehalt

          IF (L_c > 0.0.AND.L_r > 0.0) THEN
            !..Berechnung der Akkreszenzrate nach Beheng 1994

            ac = k_r *  L_c * L_r * dt

            ac = MIN(L_c,ac)

            q_rain(i,j,k)  = q_rain(i,j,k)  + ac
            q_cloud(i,j,k) = q_cloud(i,j,k) - ac
          ENDIF
        END DO
      END DO
    END DO

  END SUBROUTINE accretionKB

  SUBROUTINE autoconversionKK ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Autokonversion von Wolkenwasser zu Niederschlagswasser      *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    REAL, PARAMETER :: r_c  = 12.0e-6
    REAL, PARAMETER :: eps  = 1.00e-25
    REAL, PARAMETER :: k_a  = 3.47e+07

    !..Locale Variablen
    INTEGER          :: i,j,k
    REAL :: q_c, x_c, x_s, au

    !..Skalar-Initialisierung
    x_s  = cloud%x_max                     !..Trennmasse
    x_c  = 4./3. * pi * rho_w * r_c**3     !..Mittlere Masse der Wolkentropfen

    !..Parametrisierung nach  Khairoutdinov and Kogan (2000), MWR 128, 229-243
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_c = q_cloud(i,j,k) !..Fluessigwassergehalt
          IF (q_c > eps) THEN

            !..Berechnung der Autokonversionsrate nach KK2000
            au = k_a * (q_c*1e-3)**(0.68) * (x_c*1e3)**(1.79) * dt *1e3
            au = MIN(q_c,au)

            n_rain(i,j,k)  = n_rain(i,j,k)  + au / x_s * 2.0
            q_rain(i,j,k)  = q_rain(i,j,k)  + au
            q_cloud(i,j,k) = q_cloud(i,j,k) - au

          ENDIF
        END DO
      END DO
    END DO

  END SUBROUTINE autoconversionKK

  SUBROUTINE accretionKK ()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Akkretion von Wolkentroepfchen an Regentroepfchen       *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    REAL, PARAMETER :: k_a = 5.32e+05
    !REAL, PARAMETER :: k_a = 6.70e+01

    !..Lokale Variablen
    INTEGER          :: i,j,k
    REAL :: ac
    REAL :: L_c, L_r

    !..Parametrisierung nach Khairoutdinov and Kogan (2000), MWR 128, 229-243
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          L_c = q_cloud(i,j,k) !..Fluessigwassergehalt
          L_r = q_rain(i,j,k)  !..Fluessigwassergehalt

          IF (L_c > 0.0 .AND. L_r > 0.0) THEN
            ac  = k_a *  (L_c * L_r * 1e-6)**1.15 * dt * 1e3

            ac = MIN(L_c,ac)

            q_rain(i,j,k)  = q_rain(i,j,k)  + ac
            q_cloud(i,j,k) = q_cloud(i,j,k) - ac
          ENDIF
        END DO
      END DO
    END DO

  END SUBROUTINE accretionKK

  SUBROUTINE vapor_dep_relaxation(dt_local)
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Wachstums der Eispartikel durch Wasserdampfdiffusion    *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    REAL            :: dt_local
    REAL            :: D_vtp
    REAL            :: zdt,qvsidiff,Xi_i,Xfac
    REAL            :: tau_i_i,tau_s_i,tau_g_i,tau_h_i
    REAL, PARAMETER :: eps  = 1.e-20

    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: s_si,g_i
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dep_ice
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dep_snow
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dep_graupel, dep_graupel_n
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dep_hail, dep_hail_n

    ! Locale Variablen 
    REAL            :: T_a             !..Absolute Temperatur
    REAL            :: e_si            !..Wasserpartialdruck bei Eissaettigung
    REAL            :: e_sw            !..Wasserpartialdruck bei saettigung
    REAL            :: e_d,p_a,dep_sum
    INTEGER                     :: i,j,k

    ALLOCATE(s_si(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(g_i(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_ice(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_snow(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_graupel(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_graupel_n(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_hail(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_hail_n(0:loc_ix,1:loc_iy,1:loc_iz))

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          p_a  = p_0(i,j,k)
          T_a  = T_0(i,j,k)
          e_d  = q(i,j,k) * R_d * T_a
          e_si = esi(T_a)
          e_sw = esl(T_a)
          s_si(i,j,k) = e_d / e_si - 1.0                    !..Uebersaettigung bzgl. Eis
          D_vtp = 8.7602e-5 * T_a**(1.81) / p_a
          IF (T_a < T_3) THEN
            g_i(i,j,k) = 4.0*pi / ( L_ed**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_si) )
          ELSE
            g_i(i,j,k)  = 0.0
            s_si(i,j,k) = 0.0
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    dep_ice     = 0.0
    dep_snow    = 0.0
    dep_graupel = 0.0
    dep_hail    = 0.0

    ! ub>>
    ! Die Routinen vapor_deposition_ice(), ..., wurden zeitsparender programmiert, indem
    ! Ausdruecke wie "a ** b" durch exp(b*log(a)) ersetzt wurden und einige Berechnungen
    ! mit Konstanten vor die Schleifen gezogen wurde. Insbesondere durch die andere
    ! Potenzierungstechnik ergeben sich aber numerisch leicht unterschiedliche Ergebnisse.
    ! Im Rahmen der hier verwendeten Depositions-Approximation spielt das aber keine
    ! grosse Rolle, weil die Approximation an sich schon recht ungenau ist.
    ! ub<<

    CALL vapor_deposition_ice()
    CALL vapor_deposition_snow()
    CALL vapor_deposition_graupel()
    IF (ice_typ > 1) CALL vapor_deposition_hail()

    zdt = 1.0/dt_local

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          T_a  = T_0(i,j,k)

          ! Nur bei Temperaturen unterhalb des Tripelpunktes wird Deposition/Sublimation parametrisiert.
          ! Bei T > T_3 kann Eis nur dann verdunsten, wenn es vorher anschmilzt.

          IF (T_a < T_3) THEN

             ! Depositional growth with relaxation time-scale approach based on:
             ! "A New Double-Moment Microphysics Parameterization for Application in Cloud and
             ! Climate Models. Part 1: Description" by H. Morrison, J.A.Curry, V.I. Khvorostyanov
             ! (M05)
             
             qvsidiff  = q(i,j,k) - esi(T_a)/(R_d*T_a)

             if (abs(qvsidiff).lt.eps) then                

                dep_ice(i,j,k)     = 0.0
                dep_snow(i,j,k)    = 0.0
                dep_graupel(i,j,k) = 0.0               
                dep_hail(i,j,k)    = 0.0
                dep_sum            = 0.0

             else
             
                ! deposition rates are already multiplied with dt_local, therefore divide them here
                tau_i_i  = zdt/qvsidiff*dep_ice(i,j,k)
                tau_s_i  = zdt/qvsidiff*dep_snow(i,j,k)
                tau_g_i  = zdt/qvsidiff*dep_graupel(i,j,k)
                tau_h_i  = zdt/qvsidiff*dep_hail(i,j,k)
                
                Xi_i = ( tau_i_i + tau_s_i + tau_g_i + tau_h_i ) 
                
                if (Xi_i.lt.eps) then
                 Xfac=0.
                else
                 Xfac =  qvsidiff / Xi_i * (1.0 - EXP(- dt_local*Xi_i))
                end if
                
                dep_ice(i,j,k)     = Xfac * tau_i_i
                dep_snow(i,j,k)    = Xfac * tau_s_i
                dep_graupel(i,j,k) = Xfac * tau_g_i
                dep_hail(i,j,k)    = Xfac * tau_h_i
             
                IF (qvsidiff < 0.0) THEN
                   dep_ice(i,j,k)     = MAX(dep_ice(i,j,k),    -q_ice(i,j,k))
                   dep_snow(i,j,k)    = MAX(dep_snow(i,j,k),   -q_snow(i,j,k))
                   dep_graupel(i,j,k) = MAX(dep_graupel(i,j,k),-q_graupel(i,j,k)) 
                   dep_hail(i,j,k)    = MAX(dep_hail(i,j,k),   -q_hail(i,j,k))  
                END IF

                dep_sum = dep_ice(i,j,k) + dep_graupel(i,j,k) + dep_snow(i,j,k) + dep_hail(i,j,k)

             END IF
             
             q_graupel(i,j,k) = q_graupel(i,j,k) + dep_graupel(i,j,k)
             q_ice(i,j,k)     = q_ice(i,j,k)     + dep_ice(i,j,k)
             q_snow(i,j,k)    = q_snow(i,j,k)    + dep_snow(i,j,k)
             IF (ice_typ > 1) q_hail(i,j,k)    = q_hail(i,j,k)    + dep_hail(i,j,k)    ! <hn
             
             q(i,j,k) = q(i,j,k) - dep_sum
             
            IF (use_ice_graupel_conv_uli) THEN
               deprate_ice(i,j,k) = deprate_ice(i,j,k) + dep_ice(i,j,k)
               deprate_snow(i,j,k) = deprate_snow(i,j,k) + dep_snow(i,j,k)
            END IF
           
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(s_si,g_i,dep_ice,dep_snow,dep_graupel,dep_graupel_n, &
         &  dep_hail,dep_hail_n)

  CONTAINS

    SUBROUTINE vapor_deposition_ice()
      IMPLICIT NONE

      REAL            :: q_i,n_i,x_i,d_i,v_i,f_v,N_re,f_v_fakt
      REAL, SAVE      :: c_i             !..Koeff. fuer mittlere Kapazitaet
      REAL, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall = 0

      IF (firstcall.NE.1) THEN
        c_i = 1.0 / ice%cap
        a_f = vent_coeff_a(ice,1)
        b_f = vent_coeff_b(ice,1)
        firstcall = 1
      ENDIF

      f_v_fakt = N_sc**n_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            ! hn: in case q_ice=0, dep_ice has to be zero too

            IF (q_ice(i,j,k) <= 0.0) THEN
              dep_ice(i,j,k)   = 0.0
            ELSE  
              n_i = n_ice(i,j,k)                                   !..Anzahldichte
              q_i = q_ice(i,j,k)                                   !..Massendichte

              x_i = MIN(MAX(q_i/(n_i+eps),MAX(ice%x_min,eps)),ice%x_max)    !..mittlere Masse

              d_i = ice%a_geo * EXP(ice%b_geo*LOG(x_i))            !..mittlerer Durchmesser
              v_i = ice%a_vel * EXP(ice%b_vel*LOG(x_i)) * rrho_04(i,j,k)    !..mittlere Sedimentationsgeschw.

              N_re = v_i * d_i / nu_l                              !..mittlere Reynoldszahl
              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))     !..mittlerer Vent.Koeff.

              dep_ice(i,j,k) = g_i(i,j,k) * n_i * c_i * d_i * f_v * s_si(i,j,k) * dt_local
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_ice

    SUBROUTINE vapor_deposition_graupel()
      IMPLICIT NONE

      ! Locale Variablen 
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall = 0
      REAL            :: q_g,n_g,x_g,d_g,v_g,f_v,f_n,N_re,f_v_fakt,vent_fakt
      REAL, SAVE      :: c_g                 !..Koeff. fuer mittlere Kapazitaet
      REAL, SAVE      :: a_f,b_f,a_n,b_n     !..Koeff. fuer mittleren Ventilationkoeff.

      IF (firstcall.NE.1) THEN
        c_g = 1.0 / graupel%cap
        a_n = vent_coeff_a(graupel,0)
        b_n = vent_coeff_b(graupel,0)
        a_f = vent_coeff_a(graupel,1)
        b_f = vent_coeff_b(graupel,1)
        firstcall = 1
      ENDIF

      f_v_fakt = N_sc**n_f
      vent_fakt = b_n / b_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            ! hn: in case q_garupel=0, dep_graupel has to be zero too

            IF (q_graupel(i,j,k) <= 0.0) THEN
              dep_graupel(i,j,k)   = 0.0
              dep_graupel_n(i,j,k) = 0.0
            ELSE

              n_g = n_graupel(i,j,k)                                     !..Anzahldichte
              q_g = q_graupel(i,j,k)                                     !..Massendichte

              x_g = MIN(MAX(q_g/(n_g+eps),MAX(graupel%x_min,eps)),graupel%x_max)  !..mittlere Masse

              d_g = graupel%a_geo * EXP(graupel%b_geo*LOG(x_g))          !..mittlerer Durchmesser
              v_g = graupel%a_vel * EXP(graupel%b_vel*LOG(x_g)) * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

              N_re = v_g * d_g / nu_l                                    !..mittlere Reynoldszahl
              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))           !..mittlerer Vent.Koeff.
              f_n = a_n + vent_fakt * (f_v - a_f)                        !..mittlerer Vent.Koeff.
              f_v  = MAX(f_v,1.) !unnoetig??
              f_n  = MAX(f_n,1.) !unnoetig??

              dep_graupel(i,j,k) = g_i(i,j,k) * n_g * c_g * d_g * f_v * s_si(i,j,k) * dt_local
              dep_graupel_n(i,j,k) = dep_graupel(i,j,k) * f_n/f_v / x_g
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_graupel

    SUBROUTINE vapor_deposition_hail()
      IMPLICIT NONE

      ! Locale Variablen 
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall = 0
      REAL            :: q_h,n_h,x_h,d_h,v_h,f_v,f_n,N_re,f_v_fakt,vent_fakt
      REAL, SAVE      :: c_h                 !..Koeff. fuer mittlere Kapazitaet
      REAL, SAVE      :: a_f,b_f,a_n,b_n     !..Koeff. fuer mittleren Ventilationkoeff.

      IF (firstcall.NE.1) THEN
        c_h = 1.0 / hail%cap
        a_n = vent_coeff_a(hail,0)
        b_n = vent_coeff_b(hail,0)
        a_f = vent_coeff_a(hail,1)
        b_f = vent_coeff_b(hail,1)
        firstcall = 1
      ENDIF

      f_v_fakt = N_sc**n_f
      vent_fakt = b_n / b_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            ! hn: in case q_hail=0, dep_hail has to be zero too

            IF (q_hail(i,j,k) <= 0.0) THEN
              dep_hail(i,j,k)   = 0.0
              dep_hail_n(i,j,k) = 0.0
            ELSE
              n_h = n_hail(i,j,k)                                     !..Anzahldichte
              q_h = q_hail(i,j,k)                                     !..Massendichte

              x_h = MIN(MAX(q_h/(n_h+eps),MAX(hail%x_min,eps)),hail%x_max)  !..mittlere Masse

              d_h = hail%a_geo * EXP(hail%b_geo*LOG(x_h))          !..mittlerer Durchmesser
              v_h = hail%a_vel * EXP(hail%b_vel*LOG(x_h)) * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

              N_re = v_h * d_h / nu_l                                    !..mittlere Reynoldszahl

              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))

              f_n = a_n + vent_fakt * (f_v - a_f)                        !..mittlerer Vent.Koeff.

              f_v  = MAX(f_v,1.) !unnoetig??
              f_n  = MAX(f_n,1.) !unnoetig??

              dep_hail(i,j,k) = g_i(i,j,k) * n_h * c_h * d_h * f_v * s_si(i,j,k) * dt_local
              dep_hail_n(i,j,k) = dep_hail(i,j,k) * f_n/f_v / x_h
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_hail

    SUBROUTINE vapor_deposition_snow()
      IMPLICIT NONE

      ! Locale Variablen 
      REAL            :: q_s,n_s,x_s,d_s,v_s,f_v,N_re,f_v_fakt
      REAL, SAVE      :: c_s             !..Koeff. fuer mittlere Kapazitaet
      REAL, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall = 0

      IF (firstcall.NE.1) THEN
        c_s = 1.0 / snow%cap
        a_f = vent_coeff_a(snow,1)
        b_f = vent_coeff_b(snow,1)
        firstcall = 1
      ENDIF

      f_v_fakt = N_sc**n_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            IF (q_snow(i,j,k) <= 0.0) THEN
              dep_snow(i,j,k) = 0.
            ELSE
              n_s = n_snow(i,j,k)                        !..Anzahldichte in SI
              q_s = q_snow(i,j,k)                        !..Fluessigwassergehalt in SI

              x_s = MIN(MAX(q_s/(n_s+eps),MAX(snow%x_min,eps)),snow%x_max)   !..mittlere Masse in SI     

              d_s = snow%a_geo * EXP(snow%b_geo*LOG(x_s))           !..mittlerer Durchmesser
              v_s = snow%a_vel * EXP(snow%b_vel*LOG(x_s)) * rrho_04(i,j,k)   !..mittlere Sedimentationsgeschw.

              N_re = v_s * d_s / nu_l                         !..mittlere Reynoldszahl

              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))        !..mittlerer Vent.Koeff.
              f_v  = MAX(f_v,1.) !unnoetig??

              dep_snow(i,j,k) = g_i(i,j,k) * n_s * c_s * d_s * f_v * s_si(i,j,k) * dt_local
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_snow

  END SUBROUTINE vapor_dep_relaxation

!END MODULE wolken

END module mcrp_ice_sb