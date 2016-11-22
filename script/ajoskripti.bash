#!/bin/bash

# This script is an example how to run REMO REgional Model at FMI's
# VOIMA machine. It is quite simple and can be modified as
# pleased. The main idea is that this script is used in home directory
# and it creates all necessary folders/ files to lustre. Please check
# all paths and names! Also, please note that this script does not tar
# outputfiles and does not move them anywhere.
#
# Joni-Pekka Pietikäinen, FMI, 2013 #
#
# Modified for use of UCLA-LES_SALSA, Jaakko Ahola, FMI, 10/2015
#
#
#


# Exit on error
set -e

salsa=../src/src_salsa
les=../src/src_LES
bin=../bin/

if [ -z $1 ]; then
  echo "You didn't give any name for output directory"
  exit 1
fi

sed -i "/filprf = /c\  filprf = '"$1"'" ${bin}/NAMELIST
# sed -i "/hfilin = /c\  hfilin = '"$1".rst'" ${bin}/NAMELIST

alikansio=UCLALES-SALSA
#clock=`date +%Y-%m-%d-%H-%M`
#$ending=`echo $1`
# Main directory
wrkdir=/lustre/tmp/aholaj/${alikansio}
# Running directory
rundir=${wrkdir}/${1}
# if main directory exists -> clean it
# if not -> make it
if [ -d ${rundir} ] ; then
rm -fr ${rundir}/*
else
mkdir -p ${rundir}
fi



# datafile directory
datadir=${rundir}/datafiles

# Create all above directories to main directory
mkdir -p ${datadir} #${bdmdir} ${remdir} ${remfdir} ${othdir}

# Path to executable
exepath=/home/users/aholaj/${alikansio}/bin

# Name of the executable
exename=les.seq

# copy executable to running directory
cp ${exepath}/${exename} ${rundir}/

# Start time
#cat > ${rundir}/SIMTIME << EOF
#0
#EOF






# Link input files
cp ${exepath}/sound_in ${rundir}/
cp ${exepath}/NAMELIST ${rundir}/
cp ${exepath}/namelist.salsa ${rundir}/
cp ${exepath}/*.rst ${rundir}/

cp ${exepath}/datafiles/* ${datadir}/

length=$(( ${#1} < 6 ? ${#1} : 6))

# Create run script

cat > ${rundir}/runles.sh <<FINAL
#!/bin/sh
#PBS -N LES_${1:$((${#1}-${length})):${length}}
#PBS -l mppwidth=8
#PBS -l mppnppn=20
#PBS -l walltime=23:59:00
#PBS -j oe
#PBS -M jaakko.ahola@fmi.fi
#PBS -m ae

#export I_MPI_PLATFORM=auto
#export MPICH_ALLTOALLV_THROTTLE=2

export MPICH_ENV_DISPLAY=1

# Exit on error
set -e

cd ${rundir}

start=`date +%s`


#totalview aprun -a -n1 -N1 -d20 les.seq | tee ${PBS_JOBNAME:-interactive}.${PBS_JOBID:-help}
aprun -n 1 les.seq | tee ${PBS_JOBNAME:-interactive}.${PBS_JOBID:-help}

end=`date +%s`



cd ${rundir}

# Update ending time
#cat > SIMTIME <<EOF
#\${ENDTIME}
#EOF


# print date and time of simulation end
date
#
exit
FINAL

# Goto rundir
cd ${rundir}

# Make initial submit
qsub runles.sh

#
exit

