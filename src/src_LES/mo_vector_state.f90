MODULE mo_vector_state
  USE classFieldArray
  USE mo_structured_datatypes  
  IMPLICIT NONE

  SAVE

  TYPE(FloatArray3d), TARGET :: a_up, a_uc, a_ut
  TYPE(FloatArray3d), TARGET :: a_vp, a_vc, a_vt
  TYPE(FloatArray3d), TARGET :: a_wp, a_wc, a_wt
  
  REAL, ALLOCATABLE, TARGET :: a_vectort(:,:,:,:), a_vectorp(:,:,:,:), a_vectorc(:,:,:,:)

  
  CONTAINS

    SUBROUTINE setVectorVariables(Vector,outputlist,nzp,nxp,nyp)
      TYPE(FieldArray), INTENT(inout) :: Vector
      CHARACTER(len=*), INTENT(in) :: outputlist(:)
      INTEGER, INTENT(in) :: nzp,nxp,nyp

      INTEGER :: nvar, nn
      CLASS(*), POINTER :: pipeline_p => NULL(),  &
                           pipeline_c => NULL(),    &
                           pipeline_t => NULL()

      nvar = 3
      ALLOCATE(a_vectort(nzp,nxp,nyp,nvar),  &
               a_vectorp(nzp,nxp,nyp,nvar),  &
               a_vectorc(nzp,nxp,nyp,nvar))
      a_vectort(:,:,:,:) = 0.
      a_vectorp(:,:,:,:) = 0.
      a_vectorc(:,:,:,:) = 0.
      nn = 0
      
      nn = nn + 1
      a_up = FloatArray3d(a_vectorp(:,:,:,nn))
      a_uc = FloatArray3d(a_vectorc(:,:,:,nn))
      a_ut = FloatArray3d(a_vectort(:,:,:,nn))
      pipeline_p => a_up
      pipeline_t => a_ut
      pipeline_c => a_uc
      CALL Vector%newField("uwind","Zonal wind vector","m/s",'mttt',         &
                           ANY(outputlist == 'uwind'),in_p_data=pipeline_p,  &
                           in_t_data=pipeline_t,in_c_data=pipeline_c,        &
                           in_group='Vector')

      nn = nn + 1
      a_vp = FloatArray3d(a_vectorp(:,:,:,nn))
      a_vc = FloatArray3d(a_vectorc(:,:,:,nn))      
      a_vt = FloatArray3d(a_vectort(:,:,:,nn))
      pipeline_p => a_vp
      pipeline_t => a_vt
      pipeline_c => a_vc
      CALL Vector%newField("vwind","Meridional wind vector","m/s",'tmtt',         &
                           ANY(outputlist == 'vwind'),in_p_data=pipeline_p,  &
                           in_t_data=pipeline_t,in_c_data=pipeline_c,        &
                           in_group='Vector')

      nn = nn + 1
      a_wp = FloatArray3d(a_vectorp(:,:,:,nn))
      a_wc = FloatArray3d(a_vectorc(:,:,:,nn))
      a_wt = FloatArray3d(a_vectort(:,:,:,nn))
      pipeline_p => a_wp
      pipeline_t => a_wt
      pipeline_c => a_wc
      CALL Vector%newField("wwind","Vertical wind vector","m/s",'ttmt',         &
                           ANY(outputlist == 'wwind'),in_p_data=pipeline_p,  &
                           in_t_data=pipeline_t,in_c_data=pipeline_c,        &
                           in_group='Vector')      

      pipeline_p => NULL()
      pipeline_c => NULL()
      pipeline_t => NULL()
      
    END SUBROUTINE setVectorVariables


    
END MODULE mo_vector_state
