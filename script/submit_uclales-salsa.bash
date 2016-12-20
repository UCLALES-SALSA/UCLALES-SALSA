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
# Modified for use of UCLALES-SALSA, Jaakko Ahola, FMI, 10/2015
# Major update 20.12.2016
#
# input variables:
# $1 = name of output directory
# $2 = number of processors
# $3 = job flag of the job scheduling system (OPTIONAL) default value: PBS

# Exit on error
set -e
echo ' '
#################################
###			                  ###
### user spesific information ###
###			                  ###
#################################

username=aholaj
email=jaakko.ahola@fmi.fi
subfolder=UCLALES-SALSA
outputroot=/lustre/tmp/${username}/${subfolder}
root=/home/users/${username}/${subfolder}

# supercomputer related variable settings
WT=36:00:00 # walltime
nodeNPU=20  # number of processing units in a node  
JOBFLAG=PBS     # job flag of the job scheduling system ( e.g. PBS or SBATCH )

if [ -z $3 ]; then
  echo "You didn't give the optional JOB FLAG of the job scheduling system"
  echo "Using assumption: " $JOBFLAG
else
  JOBFLAG=$3
fi

if [ $JOBFLAG == 'PBS' ]; then
    echo 'using default values of Voima'

elif [ $JOBFLAG == 'SBATCH' ]; then ## CSC's Sisu machine values
    nodeNPU=24 # number of processing units in a node 
    QUEUE=small_long # name of the queue of Sisu machine
    outputroot=/wrk/${username}/${subfolder}
    root=/homeappl/home/${username}/appl_sisu/${subfolder}
fi

#################################
###			                  ###
### folders		              ###
###         			      ### 
#################################

salsa=${root}/src/src_salsa
les=${root}/src/src_LES
bin=${root}/bin/

################################
###			                 ###
### output directory	     ###
###			                 ###
################################

if [ -z $1 ]; then
  echo "You didn't give any name for output directory"
  exit 1
fi

################################
###			                 ###
### number of processors     ###
###			                 ###
################################
if [ -z $2 ]; then
  echo "You didn't give number of processors"
  exit 1
else
  nproc=$2
fi


################################
###			                 ###
### sequential or mpi        ###
###			                 ###
################################

if [ $nproc -gt 1 ]; then
   mode=mpi
else
   mode=seq
fi

################################
###			                 ###
### change output file names ###
###			                 ###
################################

sed -i "/filprf\s\{0,\}=\s\{0,\}/c\  filprf  = '"$1"'" ${bin}/NAMELIST
sed -i "/hfilin\s\{0,\}=\s\{0,\}/c\  hfilin  = '"$1".rst'" ${bin}/NAMELIST

################################
###			                 ###
### output directories       ###
### input files		         ###
###			                 ###
################################

rundir=${outputroot}/${1}
datadir=${rundir}/datafiles

# if main directory exists -> clean it
# if not -> make it
if [ -d ${rundir} ] ; then
   rm -fr ${rundir}/*
else
   mkdir -p ${rundir}
fi

mkdir -p ${datadir} 

# copy executables and input files to running directory
cp ${bin}/les.${mode} ${rundir}/
cp ${bin}/sound_in ${rundir}/
cp ${bin}/NAMELIST ${rundir}/
cp ${bin}/datafiles/* ${datadir}/

#########################
###			          ###
### Create run script ###
###		              ###
#########################

## modify the job name based on length: ###
length=$(( ${#1} < 6 ? ${#1} : 6))

if [ $JOBFLAG == 'PBS' ] ; then

cat > ${rundir}/runles.sh <<FINALPBS
#!/bin/sh
#PBS -N LES_${1:$((${#1}-${length})):${length}}
#PBS -l mppwidth=${nproc}
#PBS -l mppnppn=${nodeNPU}
#PBS -l walltime=${WT}
#PBS -j oe
#PBS -M ${email}
#PBS -m ae

#export I_MPI_PLATFORM=auto
#export MPICH_ALLTOALLV_THROTTLE=2

export MPICH_ENV_DISPLAY=1

# Exit on error
set -e

cd ${rundir}

aprun -n ${nproc} les.${mode} | tee ${PBS_JOBNAME:-interactive}.${PBS_JOBID:-help}

exit

FINALPBS

# Goto rundir
cd ${rundir}

# Make initial submit
qsub runles.sh

elif [ $JOBFLAG == 'SBATCH' ] ; then

cat > ${rundir}/runles.sh <<FINALSBATCH
#!/bin/sh
#SBATCH -J LES_${1:$((${#1}-${length})):${length}}
#SBATCH -n ${nproc}
#SBATCH --ntasks-per-node=${nodeNPU}
#SBATCH -t ${WT}
#SBATCH --output=LES_${1}-%j.out
#SBATCH --error=LES_${1}-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${email}
#SBATCH -p ${QUEUE}

#export I_MPI_PLATFORM=auto  
#export MPICH_ALLTOALLV_THROTTLE=2

export MPICH_ENV_DISPLAY=1

# Exit on error
set -e

cd ${rundir}

srun les.${mode}

exit

FINALSBATCH

# Goto rundir
cd ${rundir}

# Make initial submit
sbatch runles.sh

fi
###########################################    

exit