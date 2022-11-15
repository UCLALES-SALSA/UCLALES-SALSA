.SUFFIXES : .f .o .a .f90

VPATH = $(SRCUTIL):   \
	$(SRCLES):    \
	$(SRCSALSA):  \
	$(SRCRAD):    \
	$(SRCEMIS):   \
	$(SRCSHARED): \
	$(SRCIO):     \
        $(SRCVARS):

ifneq (,$(findstring mpi,$(LES_ARC)))
	seqmpi = mpi_interface.f90
	seqmpihis = mo_history_mpi.f90
        mpi_io = mo_mpi_io.f90
else
	seqmpi = seq_interface.f90
	seqmpihis = mo_history_seq.f90
	mpi_io = 
endif

SRCSF90 = mo_salsa_math.f90 mo_salsa_update.f90 mo_salsa_cloud.f90 mo_salsa_cloud_ice.f90 mo_salsa_cloud_ice_SE.f90   \
	  mo_salsa_secondary_ice.f90 mo_salsa_sizedist.f90 mo_salsa_dynamics.f90 mo_salsa_coagulation_kernels.f90     \
	  mo_salsa_coagulation_processes.f90 mo_salsa.f90 mo_salsa_driver.f90 mo_salsa_init.f90                       \
	  mo_salsa_nucleation.f90 mo_salsa_properties.f90 mo_salsa_types.f90 mo_submctl.f90                           \
	  aerosol_thermodynamics.f90 classSpecies.f90 classSection.f90 classProcessSwitch.f90 classBinLayout.f90      \
	  classProcessRates.f90 defs.f90 $(seqmpi) LES.f90 ncio.f90 util.f90 grid.f90 thrm.f90 forc.f90 advf.f90      \
	  advl.f90 srfc.f90 prss.f90 sgsm.f90 mcrp.f90 step.f90 init.f90 perturbation_forc.f90 nudg.f90                \
	  nudg_defs.f90 constrain_SALSA.f90 rad_solver.f90 rad_corkds.f90 rad_cldwtr.f90 rad_d4strm.f90      \
	  rad_driver.f90 radiation_main.f90 emission_main.f90 emission_init.f90 mo_lagrangian_tracker.f90             \
	  emission_types.f90 mo_seasalt_emission.f90 exceptionHandling.f90 mo_aerorad_lut.f90                         \
	  mo_salsa_optical_properties.f90 mo_particle_external_properties.f90 mo_ice_shape.f90                        \
	  mo_structured_datatypes.f90 mo_aux_state.f90 mo_diag_state.f90 mo_progn_state.f90 mo_vector_state.f90       \
	  mo_derived_state.f90 mo_derived_procedures.f90 mo_ps_state.f90 mo_ps_procedures.f90 mo_ts_state.f90         \
          mo_ts_procedures.f90 classFieldArray.f90 mo_output.f90 $(seqmpihis) mo_field_state.f90 mo_field_init.f90    \
          mo_check_state.f90 mo_stats_finder.f90 mo_stats_parameters.f90 $(mpi_io) math_functions.f90


OBJS =	$(SRCSF90:.f90=.o) $(SRCSF77:.f=.o)


all: $(OUT)


$(OUT): $(OBJS)
	$(F90) $(FFLAGS) -o $@ $(OBJS) rfft.o $(LIBS)

rfft.o: rfft.f90
	$(F90) -c $(F77FLAGS) $<

%.o: %.f90 
	$(F90) -c $(FFLAGS) $<

#
# Dependencies
# --------------------
#
LES.o: $(seqmpi:.f90=.o) grid.o init.o step.o mo_salsa_init.o \
       util.o sgsm.o srfc.o nudg_defs.o emission_types.o perturbation_forc.o   \
       forc.o radiation_main.o mcrp.o mo_output.o mo_stats_parameters.o

ncio.o: $(seqmpi:.f90=.o) classFieldArray.o
mo_output.o: $(seqmpi:.f90=.o) ncio.o grid.o classFieldArray.o mo_structured_datatypes.o  \
	     mo_submctl.o mo_field_state.o
$(seqmpihis:.f90=.o): $(seqmpi:.f90=.o) grid.o nudg_defs.o mo_progn_state.o   \
	              mo_vector_state.o mo_diag_state.o mo_aux_state.o $(mpi_io:.f90=.o)
mo_mpi_io.o: $(seqmpi:.f90=.o)
util.o: $(seqmpi:.f90=.o) rfft.o
math_functions.o: 
defs.o:
grid.o: $(seqmpi:.f90=.o) mo_submctl.o     \
	classFieldArray.o mo_aux_state.o
thrm.o: mo_structured_datatypes.o mo_diag_state.o mo_progn_state.o mo_aux_state.o \
        grid.o defs.o $(seqmpi:.f90=.o)
forc.o: grid.o mo_aux_state.o mo_diag_state.o mo_progn_state.o mo_vector_state.o \
	mo_field_state.o util.o defs.o radiation_main.o nudg.o emission_main.o \
	mo_structured_datatypes.o $(seqmpi:.f90=.o)
advf.o: mo_structured_datatypes.o mo_diag_state.o mo_progn_state.o mo_vector_state.o   \
	mo_aux_state.o grid.o util.o mo_submctl.o $(seqmpi:.f90=.o)	
advl.o: mo_structured_datatypes.o grid.o mo_vector_state.o mo_aux_state.o util.o 
srfc.o: mo_structured_datatypes.o defs.o mo_diag_state.o mo_aux_state.o     \
        mo_vector_state.o mo_progn_state.o grid.o thrm.o defs.o $(seqmpi:.f90=.o) 
prss.o: defs.o mo_structured_datatypes.o mo_aux_state.o mo_vector_state.o mo_diag_state.o    \
        grid.o util.o $(seqmpi:.f90=.o)
sgsm.o: util.o mo_structured_datatypes.o mo_aux_state.o mo_diag_state.o  \
        mo_vector_state.o mo_progn_state.o grid.o thrm.o defs.o $(seqmpi:.f90=.o)
mcrp.o: defs.o mo_aux_state.o mo_diag_state.o mo_progn_state.o grid.o  \
	thrm.o mo_submctl.o mo_particle_external_properties.o util.o   \
	classProcessSwitch.o mo_structured_datatypes.o mo_ice_shape.o
step.o: mo_submctl.o util.o mo_structured_datatypes.o mo_aux_state.o                 \
	mo_vector_state.o mo_progn_state.o mo_diag_state.o grid.o thrm.o mo_output.o \
	$(seqmpihis:.f90=.o) mcrp.o grid.o mo_field_state.o sgsm.o srfc.o                    \
	prss.o advf.o advl.o forc.o util.o mo_salsa_driver.o constrain_SALSA.o       \
	defs.o $(seqmpi:.f90=.o) perturbation_forc.o
init.o: grid.o mo_progn_state.o mo_vector_state.o mo_diag_state.o mo_derived_state.o \
	mo_aux_state.o mo_ts_state.o mo_ps_state.o mo_field_state.o mo_field_init.o  \
	step.o sgsm.o thrm.o mo_salsa_driver.o mo_submctl.o util.o nudg.o            \
	emission_init.o constrain_SALSA.o mo_structured_datatypes.o mo_output.o      \
	$(seqmpihis:.f90=.o) defs.o mo_submctl.o mo_salsa_sizedist.o      \
	mo_salsa_types.o ncio.o mo_check_state.o $(seqmpi:.f90=.o) srfc.o
perturbation_forc.o: defs.o grid.o mo_aux_state.o mo_progn_state.o math_functions.o srfc.o
nudg.o: mo_aux_state.o mo_diag_state.o mo_progn_state.o mo_vector_state.o grid.o mo_submctl.o nudg_defs.o util.o
nudg_defs.o:

mo_structured_datatypes.o:
classFieldArray.o: mo_structured_datatypes.o
mo_field_state.o: classFieldArray.o
mo_field_init.o:  $(seqmpi:.f90=.o) mo_progn_state.o mo_diag_state.o mo_derived_state.o mo_vector_state.o \
                  mo_aux_state.o mo_ts_state.o mo_ps_state.o mo_field_state.o mo_check_state.o grid.o
mo_aux_state.o: classFieldArray.o mo_structured_datatypes.o mo_submctl.o
mo_diag_state.o: classFieldArray.o mo_structured_datatypes.o util.o mo_submctl.o
mo_progn_state.o: classFieldArray.o mo_structured_datatypes.o mo_submctl.o 
mo_vector_state.o: classFieldArray.o mo_structured_datatypes.o 
mo_derived_state.o: classFieldArray.o mo_structured_datatypes.o mo_derived_procedures.o 
mo_derived_procedures.o: mo_submctl.o util.o defs.o mo_particle_external_properties.o grid.o \
                         mo_progn_state.o mo_diag_state.o mo_aux_state.o mo_structured_datatypes.o \
                         math_functions.o
mo_ps_state.o: classFieldArray.o mo_structured_datatypes.o mo_ps_procedures.o 
mo_ps_procedures.o: util.o mo_structured_datatypes.o grid.o mo_stats_finder.o mo_aux_state.o \
                    mo_vector_state.o advf.o advl.o defs.o $(seqmpi:.f90=.o)
mo_ts_state.o: classFieldArray.o mo_structured_datatypes.o mo_ts_procedures.o
mo_ts_procedures.o: grid.o mo_stats_finder.o util.o mo_progn_state.o mo_diag_state.o mo_aux_state.o  \
                    mo_structured_datatypes.o mo_stats_parameters.o
mo_check_state.o: $(seqmpi:.f90=.o) classFieldArray.o
mo_stats_finder.o: mo_structured_datatypes.o mo_field_state.o
mo_stats_parameters.o:

constrain_SALSA.o: mo_progn_state.o mo_diag_state.o mo_aux_state.o mo_submctl.o nudg_defs.o \
                   util.o mo_structured_datatypes.o grid.o mo_derived_procedures.o          \
	           mo_field_state.o
mo_particle_external_properties.o: mo_submctl.o classSection.o mo_ice_shape.o
mo_ice_shape.o: mo_submctl.o

rfft.o:
radiation_main.o: grid.o rad_driver.o mo_submctl.o mo_aux_state.o mo_diag_state.o mo_progn_state.o
rad_solver.o: defs.o
rad_corkds.o: defs.o $(seqmpi:.f90=.o)
rad_cldwtr.o: defs.o rad_corkds.o mo_aerorad_lut.o mo_salsa_optical_properties.o
rad_d4strm.o: defs.o rad_solver.o \
	      rad_corkds.o rad_cldwtr.o
rad_driver.o: defs.o rad_d4strm.o $(seqmpi:.f90=.o) 
mo_aerorad_lut.o:
$(seqmpi:.f90=.o):

emission_main.o: emission_types.o mo_seasalt_emission.o mo_submctl.o mo_salsa_types.o    \
	         mo_salsa_sizedist.o grid.o util.o exceptionHandling.o    \
	         $(seqmpi:.f90=.o) mo_aux_state.o mo_vector_state.o mo_progn_state.o mo_diag_state.o
emission_types.o: 
emission_init.o: mo_submctl.o emission_types.o mo_lagrangian_tracker.o mo_salsa_types.o  \
	         util.o mo_salsa_sizedist.o exceptionHandling.o          \
	         mo_aux_state.o grid.o $(seqmpi:.f90=.o)
mo_lagrangian_tracker.o: mo_aux_state.o grid.o util.o

exceptionHandling.o:

mo_salsa_driver.o: classSection.o util.o mo_salsa_types.o mo_submctl.o classFieldArray.o   \
	           mo_structured_datatypes.o mo_salsa.o mo_salsa_properties.o
mo_salsa_init.o: mo_salsa_types.o mo_submctl.o mo_salsa_driver.o mo_salsa_optical_properties.o \
	         classSpecies.o classSection.o classProcessSwitch.o classBinLayout.o          \
	         mo_salsa_coagulation_kernels.o mo_salsa_coagulation_processes.o classProcessRates.o  \
                 mo_ice_shape.o mo_salsa_cloud_ice_SE.o mo_salsa_secondary_ice.o
mo_submctl.o: classSpecies.o classProcessSwitch.o classBinLayout.o
mo_salsa_cloud.o: mo_salsa_types.o mo_submctl.o classSection.o mo_particle_external_properties.o   \
	          mo_salsa_math.o math_functions.o
mo_salsa_cloud_ice.o: mo_submctl.o mo_salsa_types.o mo_particle_external_properties.o mo_salsa_math.o math_functions.o
mo_salsa_cloud_ice_SE.o: mo_submctl.o mo_salsa_types.o mo_salsa_math.o mo_particle_external_properties.o  \
	                 classSection.o math_functions.o
mo_salsa_secondary_ice.o: mo_salsa_types.o mo_submctl.o
mo_salsa_sizedist.o: mo_salsa_types.o mo_submctl.o 
mo_salsa.o: mo_submctl.o mo_salsa_dynamics.o mo_salsa_update.o \
	    mo_salsa_cloud.o classSection.o mo_salsa_cloud_ice.o mo_salsa_cloud_ice_SE.o
mo_salsa_dynamics.o: mo_salsa_types.o mo_submctl.o mo_salsa_nucleation.o \
                     mo_salsa_properties.o aerosol_thermodynamics.o \
                     classSection.o mo_particle_external_properties.o \
	             mo_salsa_coagulation_kernels.o mo_salsa_coagulation_processes.o  \
	             mo_salsa_secondary_ice.o
mo_salsa_coagulation_kernels.o: mo_salsa_types.o mo_submctl.o classSection.o mo_particle_external_properties.o mo_ice_shape.o
mo_salsa_coagulation_processes.o: mo_salsa_types.o mo_submctl.o classSection.o mo_particle_external_properties.o  
mo_salsa_nucleation.o: mo_salsa_types.o mo_submctl.o 
mo_salsa_properties.o: mo_salsa_types.o mo_submctl.o
mo_salsa_optical_properties.o: classSpecies.o mo_submctl.o
mo_salsa_update.o: mo_salsa_types.o mo_submctl.o classSection.o
mo_salsa_types.o: classSection.o classProcessRates.o
mo_salsa_math.o: mo_submctl.o
aerosol_thermodynamics.o:
classSpecies.o:
classSection.o: mo_submctl.o mo_ice_shape.o
classProcessSwitch.o: 
classBinLayout.o:
classProcessRates.o: 


