###############################################################
#
# Llocation of code (in $ROOT) and location where model is to be built $BIN
#
ROOT      :=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
BIN       = $(ROOT)/bin
ARCH      := $(shell uname)
#
# Generic Variables
#
SRC         =$(ROOT)/src
SRC_UTIL    =$(SRC)/src_util
SRC_LES     =$(SRC)/src_LES
SRC_SALSA   =$(SRC)/src_salsa
SRC_RAD     =$(SRC)/src_rad
SRC_EMIS    =$(SRC)/src_emission
SRC_SHARED  =$(SRC)/src_shared

VPATH = $(SRC_LES):$(SRC_SALSA):$(SRC_UTIL):$(SRC_RAD):  \
	$(SRC_EMIS):$(SRC_SHARED):$(SRC)

ECHO    = /bin/echo
RM      = /bin/rm -f

ARCHIVE = ar rs
RANLIB =:
SEQFFLAGS = -I$(SRC)
MPIFFLAGS = -I$(SRC)
MPI = /usr
MPILIB = -L$(MPI)/lib -lmpich
MPIINC = -I$(MPI)/include/mpich
NCDF = /usr
NCDFLIB = -L$(NCDF)/lib -lnetcdf -lnetcdff
NCDFINC = -I$(NCDF)/include
LIBS = '$(NCDFLIB) $(MPILIB)'
F90 = gfortran
MPIF90 = gfortran
FFLAGS = -O2 -fdefault-real-8 ${NCDFINC} $(MPIINC) -std=f2008 -fbounds-check  -g -fcheck=all -Wall -Wtabs  -fbacktrace -ffpe-trap=invalid,zero,overflow
F77FLAGS = -O2 

LES_OUT_MPI=$(BIN)/les.mpi

LES_OUT_SEQ=$(BIN)/les.seq

default: mpi

all:  mpi seq

seq: $(LES_OUT_SEQ)

mpi: $(LES_OUT_MPI)

$(LES_OUT_SEQ): 
	cd $(SRC); $(MAKE) LES_ARC=seq                     \
	FFLAGS='$(FFLAGS) $(SEQFFLAGS)' F90=$(F90)         \
	F77FLAGS='$(F77FLAGS)' OUT=$(LES_OUT_SEQ)          \
	LIBS=$(LIBS) SRCUTIL=$(SRC_UTIL) SRCLES=$(SRC_LES) \
	SRCSALSA=$(SRC_SALSA) SRCRAD=$(SRC_RAD)            \
	SRCEMIS=$(SRC_EMIS) SRCSHARED=$(SRC_SHARED)
 
$(LES_OUT_MPI):
	cd $(SRC); $(MAKE) LES_ARC=mpi                     \
	FFLAGS='$(FFLAGS) $(MPIFFLAGS)' F90=$(MPIF90)      \
	F77FLAGS='$(F77FLAGS)' OUT=$(LES_OUT_MPI)          \
	LIBS=$(LIBS) SRCUTIL=$(SRC_UTIL) SRCLES=$(SRC_LES) \
	SRCSALSA=$(SRC_SALSA) SRCRAD=$(SRC_RAD)            \
	SRCEMIS=$(SRC_EMIS) SRCSHARED=$(SRC_SHARED)

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
