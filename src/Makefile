.SUFFIXES : .f .o .a .f90



VPATH = $(SRCUTIL): \
	$(SRCLES):  \
	$(SRCSALSA)


ifneq (,$(findstring mpi,$(LES_ARC)))
	seqmpi = mpi_interface.f90
else
	seqmpi = seq_interface.f90
endif

SRCSF90 = mo_salsa_update.f90 mo_salsa_cloud.f90 \
          mo_salsa_dynamics.f90 \
          mo_salsa.f90 mo_salsa_driver.f90 mo_salsa_init.f90  \
          mo_salsa_properties.f90 \
          mo_submctl.f90 \
          defs.f90 $(seqmpi) rad_solver.f90 rad_corkds.f90 LES.f90 \
          rad_cldwtr.f90 rad_d4strm.f90 rad_driver.f90 \
          ncio.f90 util.f90 grid.f90 stat.f90 thrm.f90 forc.f90 advf.f90 \
          advl.f90 srfc.f90 prss.f90 sgsm.f90 mcrp.f90 mcrp_ice.f90 mcrp_ice_sb.f90 step.f90 init.f90 \
          rfft.f90 \
          mo_vbsctl.f90 mo_vbs.f90 mo_vbs_init.f90 mo_vbs_partition.f90
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
       stat.o forc.o mo_salsa_init.o
ncio.o: $(seqmpi:.f90=.o)
util.o: $(seqmpi:.f90=.o)
defs.o:
grid.o: ncio.o $(seqmpi:.f90=.o) mo_submctl.o
stat.o: ncio.o grid.o util.o defs.o mo_submctl.o $(seqmpi:.f90=.o)
thrm.o: grid.o defs.o $(seqmpi:.f90=.o) mo_submctl.o 
forc.o: grid.o defs.o stat.o util.o rad_driver.o mo_submctl.o $(seqmpi:.f90=.o)
advf.o: grid.o stat.o util.o mo_submctl.o $(seqmpi:.f90=.o)
advl.o: grid.o stat.o util.o
srfc.o: grid.o stat.o defs.o thrm.o util.o $(seqmpi:.f90=.o)
prss.o: grid.o stat.o defs.o util.o $(seqmpi:.f90=.o)
sgsm.o: grid.o stat.o defs.o util.o thrm.o $(seqmpi:.f90=.o)
mcrp.o: grid.o stat.o defs.o mo_submctl.o mcrp_ice.o
mcrp_ice_sb.o: defs.o thrm.o
mcrp_ice.o: grid.o stat.o defs.o thrm.o mcrp_ice_sb.o
step.o: grid.o stat.o util.o defs.o thrm.o sgsm.o srfc.o  \
        mcrp.o prss.o advf.o advl.o forc.o $(seqmpi:.f90=.o) \
        mo_salsa_driver.o mo_submctl.o
init.o: grid.o stat.o util.o defs.o step.o sgsm.o thrm.o  \
	$(seqmpi:.f90=.o) mo_submctl.o \
	mo_salsa_driver.o mo_submctl.o \
	ncio.o
rfft.o:
rad_solver.o: defs.o
rad_corkds.o: defs.o $(seqmpi:.f90=.o)
rad_cldwtr.o: defs.o rad_corkds.o
rad_d4strm.o: defs.o rad_solver.o \
	      rad_corkds.o rad_cldwtr.o
rad_driver.o: defs.o rad_d4strm.o $(seqmpi:.f90=.o)
$(seqmpi:.f90=.o):
mo_salsa_driver.o: mo_submctl.o mo_salsa.o mo_salsa_properties.o
mo_salsa_init.o: mo_submctl.o mo_salsa_driver.o mo_vbs_init.o step.o 
mo_submctl.o: 
mo_salsa_cloud.o: mo_submctl.o
mo_salsa.o: mo_submctl.o mo_salsa_dynamics.o \
          mo_salsa_update.o mo_salsa_cloud.o mo_vbs_partition.o
mo_salsa_dynamics.o: mo_submctl.o mo_salsa_properties.o
mo_salsa_properties.o: mo_submctl.o
mo_salsa_update.o: mo_submctl.o
mo_vbsctl.o:
mo_vbs_init.o: mo_submctl.o mo_vbs.o mo_vbsctl.o
mo_vbs.o: mo_vbsctl.o 
mo_vbs_partition.o: mo_submctl.o mo_vbsctl.o mo_vbs.o grid.o
