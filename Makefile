###############################################################
#
# Location of code (in $ROOT) and location where model is to be built $BIN
#
# when adding a new compiler environment these are needed parameters for compiler:
#  $(F90)  	compiler command
#  $(FFLAGS) 	flags for compiler F90
#  $(F77FLAGS) 	flags for compiler F77 (optional)
#  $(LIBFLAGS) 	links
#  $(LIBS)	libraries for netcdf4 and hdf5
#
# example use: make mpi COMP=taitointel RUNTYPE=debug
#	default values: COMP=intel RUNTYPE=fast

ROOT      :=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
BIN       = $(ROOT)/bin
ARCH      := $(shell uname)
VERS      := $(shell git symbolic-ref --short -q HEAD)

#
# Generic Variables
#
SRC         =$(ROOT)/src
SRC_UTIL    =$(SRC)/src_util
SRC_LES     =$(SRC)/src_LES
SRC_SALSA   =$(SRC)/src_salsa

VPATH = $(SRC_LES):$(SRC_SALSA):$(SRC_UTIL):$(SRC)

ECHO	= /bin/echo
RM      = /bin/rm -f

ifndef $(COMP)
	COMP=intel
endif

ifndef $(RUNTYPE)
	RUNTYPE=fast
endif

$(info $$COMP is [${COMP}])
$(info $$RUNTYPE is [${RUNTYPE}])
# Ubuntu -------------------------------------------------
ifeq ($(COMP),ubuntu)
	F90 = f95

	NCDF		= /usr
	NCDFLIB		= '-L$(NCDF)/lib -lnetcdff -lnetcdf'
	NCDFINC		= -I$(NCDF)/include

	# Libraries
	LIBS		= $(NCDFLIB)
	LIBFLAGS 	= -I$(SRC)

	ifeq ($(RUNTYPE),fast)
		# Optimized
		FFLAGS		= -O2 -fdefault-real-8 ${NCDFINC} 
		F77FLAGS	= -O2 
	else
		#Debug
		FFLAGS		= -O2 -fdefault-real-8 ${NCDFINC} -fbounds-check  -g -fcheck=all  -Wall -Wtabs -fbacktrace -ffpe-trap=invalid,zero,overflow
		F77FLAGS 	= -O2 -fbounds-check  -ffpe-trap=invalid,zero,overflow
	endif
endif
# Gnu -------------------------------------------------
ifeq ($(COMP),gnu)
	F90 = ftn

	NCDF		= /usr
	NCDFLIB		= '-L$(NCDF)/lib -lnetcdff -lnetcdf'
	NCDFINC		= -I$(NCDF)/include
	
	# Libraries
	LIBS		= $(NCDFLIB)
	LIBFLAGS	= -I$(SRC)

	ifeq ($(RUNTYPE),fast)
		# Optimized
		FFLAGS 		= -O2 -fdefault-real-8 -fbounds-check ${NCDFINC}
		F77FLAGS 	= -O2
	else
		# Debug
		FFLAGS		= -O0 -fdefault-real-8 -fbounds-check ${NCDFINC} -ggdb3 -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all
		F77FLAGS 	= -O0 -fbounds-check -ffpe-trap=invalid,zero,overflow
	endif
endif
# Cray -------------------------------------------------
ifeq ($(COMP),cray)
	F90 = ftn

	# Libraries
	LIBFLAGS = -I$(SRC)

	ifeq ($(RUNTYPE),fast)
		# Optimized
		FFLAGS 		= -O2 -s real64 -h fp0 # Note: flag "-h fp0" can be required!
		F77FLAGS 	= -O2
	else
		# Debug
		FFLAGS		= -O0 -s real64 -eD -G0
		F77FLAGS	= -G0
	endif
	
endif
# Taito Intel -------------------------------------------------
ifeq ($(COMP),taitointel)
	F90		= mpif90

	NETCDROOT	= /appl/opt/netcdf4/intel-16.0.0/intelmpi-5.1.1/4.3.3.1/
	NETCDF_LIB	= -L$(NETCDFROOT)lib -lnetcdff -lnetcdf 
	NETCDF_INCLUDE	= -I$(NETCDFROOT)/include

	HDF5ROOT	= /appl/opt/hdf5-par/intel-16.0.0/intelmpi-5.1.1/1.8.15
	HDF5_LIB	= -L$(HDF5ROOT)/lib -lhdf5_hl -lhdf5
	HDF5_INCLUDE       = -I$(HDF5ROOT)/include

	# Libraries
	LIBS		= '$(HDF5_LIB) $(NETCDF_LIB)'
	LIBFLAGS 	= -I$(SRC) $(HDF5_INCLUDE) $(NETCDF_INCLUDE)

	ifeq ($(RUNTYPE),fast)
		# Optimized
		FFLAGS		= -O2 -march=native -real-size 64 -convert big_endian -fpe0
		F77FLAGS	= -O2 -march=native -real-size 64 -convert big_endian -fpe0
	else
		# Debug
		FFLAGS		= -O2 -march=native -real-size 64 -convert big_endian -fpe0 -fp-model source -fp-model precise -g -traceback -integer-size 32 -check bounds
		F77FLAGS	= -O2 -march=native -real-size 64 -convert big_endian -fpe0 -fp-model source -fp-model precise -g -traceback -integer-size 32 -check bounds
	endif
endif
# Intel -------------------------------------------------
ifeq ($(COMP),intel)
	F90 = ftn

	NETCDFROOT	= /opt/cray/netcdf/4.3.0/intel/130
	NETCDF_LIB	= -L$(NETCDFROOT)/lib -lnetcdff -lnetcdf
	NETCDF_INCLUDE	= -I$(NETCDFROOT)/include

	HDF5ROOT	= /opt/cray/hdf5/1.8.11/intel/130
	HDF5_LIB	= -L$(HDF5ROOT)/lib -lhdf5_hl -lhdf5
	HDF5_INCLUDE	= -I$(HDF5ROOT)/include

	# Libraries
	LIBS		= '$(HDF5_LIB) $(NETCDF_LIB)'
	LIBFLAGS	= $(HDF5_INCLUDE) $(NETCDF_INCLUDE)

	ifeq ($(RUNTYPE),fast)
		# Optimized
		FFLAGS		= -O2 -msse2 -real-size 64 -fp-model precise -convert big_endian
		F77FLAGS	= -O2 -msse2 -real-size 64 -fp-model precise -convert big_endian
	else
		# Debug
		FFLAGS		= -O2 -msse2 -fp-model source -fp-model precise -g -traceback -convert big_endian -real-size 64 -check bounds -fpe0
		F77FLAGS	= -O2 -msse2 -fp-model source -fp-model precise -g -traceback -convert big_endian -real-size 64 -check bounds -fpe0
	endif
endif
# -------------------------------------------------------

LES_OUT_MPI=$(BIN)/les.mpi.$(VERS).$(COMP).$(RUNTYPE)

LES_OUT_SEQ=$(BIN)/les.seq.$(VERS).$(COMP).$(RUNTYPE)

default: mpi

all:  mpi seq

seq: $(LES_OUT_SEQ)

mpi: $(LES_OUT_MPI)

$(LES_OUT_SEQ): 
	cd $(SRC); $(MAKE) LES_ARC=seq \
	FFLAGS='$(FFLAGS) $(LIBFLAGS)' F90=$(F90) \
	F77FLAGS='$(F77FLAGS)' OUT=$(LES_OUT_SEQ) \
	LIBS=$(LIBS) SRCUTIL=$(SRC_UTIL) SRCLES=$(SRC_LES) \
	SRCSALSA=$(SRC_SALSA)

$(LES_OUT_MPI):
	cd $(SRC); $(MAKE) LES_ARC=mpi \
	FFLAGS='$(FFLAGS) $(LIBFLAGS)' F90=$(F90)  \
	F77FLAGS='$(F77FLAGS)' OUT=$(LES_OUT_MPI) \
	LIBS=$(LIBS) SRCUTIL=$(SRC_UTIL) SRCLES=$(SRC_LES) \
	SRCSALSA=$(SRC_SALSA)

.PHONY: $(LES_OUT_SEQ) 
.PHONY: $(LES_OUT_MPI)

#
# cleaning
# --------------------
#
clean: cleanmpi cleanseq 
	$(RM) $(SRC)/*mod $(SRC)/*.o

cleanmpi:
	$(ECHO) "cleaning mpi model"
	$(RM) core $(LES_OUT_MPI) $(SRC)/mpi/*mod $(LES_ARC_MPI)

cleanseq:
	$(ECHO) "clean sequential model"
	$(RM) core $(LES_OUT_SEQ) $(SRC)/seq/*mod $(LES_ARC_SEQ)

FORCE: 
.PRECIOUS: $(LIBS)
