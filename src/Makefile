.SUFFIXES : .f .o .a .f90

VPATH = $(SRCUTIL):   \
	$(SRCLES):    \
	$(SRCSALSA):  \
	$(SRCRAD):    \
	$(SRCEMIS):   \
	$(SRCSHARED)  

ifneq (,$(findstring mpi,$(LES_ARC)))
	seqmpi = mpi_interface.f90
else
	seqmpi = seq_interface.f90
endif

SRCSF90 = mo_salsa_update.f90 mo_salsa_cloud.f90 \
          mo_salsa_sizedist.f90 mo_salsa_dynamics.f90 \
          mo_salsa.f90 mo_salsa_driver.f90 mo_salsa_init.f90  \
          mo_salsa_nucleation.f90 mo_salsa_properties.f90 \
          mo_submctl.f90 aerosol_thermodynamics.f90       \
	  classSpecies.f90 classSection.f90 classProcessSwitch.f90 \
          defs.f90 $(seqmpi) LES.f90 ncio.f90 util.f90 grid.f90 stat.f90 \
	  thrm.f90 forc.f90 advf.f90 advl.f90 srfc.f90 prss.f90 sgsm.f90 \
	  mcrp.f90 step.f90 init.f90 nudg.f90 \
          rfft.f90 rad_solver.f90 rad_corkds.f90 rad_cldwtr.f90 \
	  rad_d4strm.f90 rad_driver.f90 radiation_main.f90 emission_main.f90 \
	  mo_seasalt_emission.f90 exceptionHandling.f90 mo_aerorad_lut.f90  \
	  mo_salsa_optical_properties.f90 mo_particle_external_properties.f90
SRCSF77 = 


OBJS =	$(SRCSF90:.f90=.o) $(SRCSF77:.f=.o)

all: $(OUT)


$(OUT): $(OBJS)
	$(F90) $(FFLAGS) $(SEQFFLAGS) -o $@ $(OBJS) $(LIBS)

%.o: %.f90 
	$(F90) -c $(FFLAGS) $<

%.o: %.f
	$(F90) -c $(F77FLAGS) $<

#
# Dependencies
# --------------------
#
LES.o: grid.o init.o step.o $(seqmpi:.f90=.o) util.o sgsm.o srfc.o \
       stat.o forc.o mo_salsa_init.o nudg.o emission_main.o radiation_main.o
ncio.o: $(seqmpi:.f90=.o)
util.o: $(seqmpi:.f90=.o) 
defs.o:
grid.o: ncio.o $(seqmpi:.f90=.o) mo_submctl.o util.o 
stat.o: ncio.o grid.o util.o defs.o mo_submctl.o $(seqmpi:.f90=.o)
thrm.o: grid.o defs.o 
forc.o: grid.o defs.o stat.o $(seqmpi:.f90=.o) nudg.o radiation_main.o \
	emission_main.o
advf.o: grid.o stat.o util.o mo_submctl.o $(seqmpi:.f90=.o) 
advl.o: grid.o stat.o util.o
srfc.o: grid.o stat.o defs.o thrm.o $(seqmpi:.f90=.o)
prss.o: grid.o stat.o defs.o util.o $(seqmpi:.f90=.o)
sgsm.o: grid.o stat.o defs.o util.o thrm.o $(seqmpi:.f90=.o)
mcrp.o: grid.o stat.o defs.o thrm.o mo_submctl.o util.o mo_particle_external_properties.o
step.o: grid.o stat.o util.o defs.o thrm.o sgsm.o srfc.o  \
        mcrp.o prss.o advf.o advl.o forc.o $(seqmpi:.f90=.o) \
        mo_salsa_driver.o mo_submctl.o 
init.o: grid.o stat.o util.o defs.o step.o sgsm.o thrm.o  \
	$(seqmpi:.f90=.o) mo_salsa_sizedist.o mo_submctl.o \
	mo_salsa_driver.o mo_submctl.o ncio.o nudg.o emission_main.o 
nudg.o: grid.o mo_submctl.o util.o
emission_main.o: mo_seasalt_emission.o mo_submctl.o grid.o util.o exceptionHandling.o \
	         mo_salsa_sizedist.o mo_salsa_driver.o 
mo_seasalt_emission.o: mo_submctl.o
mo_particle_external_properties.o: mo_submctl.o classSection.o
rfft.o:
radiation_main.o: grid.o rad_driver.o mo_submctl.o 
rad_solver.o: defs.o
rad_corkds.o: defs.o $(seqmpi:.f90=.o)
rad_cldwtr.o: defs.o rad_corkds.o mo_aerorad_lut.o mo_salsa_optical_properties.o
rad_d4strm.o: defs.o rad_solver.o \
	      rad_corkds.o rad_cldwtr.o
rad_driver.o: defs.o rad_d4strm.o $(seqmpi:.f90=.o) 
mo_aerorad_lut.o:
$(seqmpi:.f90=.o):

exceptionHandling.o:

mo_salsa_driver.o: mo_submctl.o mo_salsa.o mo_salsa_properties.o \
	           classSection.o util.o classProcessSwitch.o
mo_salsa_init.o: mo_submctl.o mo_salsa_driver.o mo_salsa_optical_properties.o \
	         classSpecies.o classSection.o classProcessSwitch.o
mo_submctl.o: classSpecies.o classSection.o classProcessSwitch.o
mo_salsa_cloud.o: mo_submctl.o classSection.o mo_particle_external_properties.o util.o
mo_salsa_sizedist.o: mo_salsa_driver.o mo_submctl.o 
mo_salsa.o: mo_submctl.o mo_salsa_dynamics.o mo_salsa_update.o \
	    mo_salsa_cloud.o classSection.o
mo_salsa_dynamics.o: mo_submctl.o mo_salsa_nucleation.o \
                     mo_salsa_properties.o aerosol_thermodynamics.o \
                     classSection.o mo_particle_external_properties.o 
mo_salsa_nucleation.o: mo_submctl.o 
mo_salsa_properties.o: mo_submctl.o
mo_salsa_optical_properties.o: classSpecies.o mo_submctl.o
mo_salsa_update.o: mo_submctl.o classSection.o
aerosol_thermodynamics.o:
classSpecies.o:
classSection.o:
classProcessSwitch.o: 


