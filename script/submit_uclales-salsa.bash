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
# 
# it is recommended to give input folder value:
# e.g. input=/home/users/aholaj/UCLALES-SALSA/bin/case_emulator/emul01 ./submit_uclales-salsa.bash
#

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


# supercomputer related variable settings
WT=${WT:-24:00:00} # walltime

JOBFLAG=${JOBFLAG:-PBS}     # job flag of the job scheduling system ( e.g. PBS or SBATCH )
COPY=${COPY:-true}
clean=${clean:-true}

if [ -z $3 ]; then
  echo ' '
  echo "You didn't give the optional JOB FLAG of the job scheduling system"
  echo "Using assumption: " $JOBFLAG
else
  JOBFLAG=$3
fi

if [ $JOBFLAG == 'PBS' ]; then
    outputroot=/lustre/tmp/${username}/${subfolder}/${subsubfolder}
    root=/home/users/${username}/${subfolder}
    nodeNPU=20  # number of processing units in a node  

elif [ $JOBFLAG == 'SBATCH' ]; then ## CSC's Sisu machine values
    nodeNPU=24 # number of processing units in a node 
    QUEUE=small_long # name of the queue of Sisu machine
    outputroot=/wrk/${username}/${subfolder}
    root=/homeappl/home/${username}/appl_sisu/${subfolder}
fi
echo ' '
echo 'SUPERCOMPUTER VALUES'



#################################
###			                  ###
### folders		              ###
###         			      ### 
#################################


bin=${root}/bin
input=${input:-}


if [ -z $input ]; then
    input=$bin
fi

################################
###			                 ###
### output directory	     ###
###			                 ###
################################

if [ -z $1 ]; then
  echo "You didn't give any name for output directory"
  exit 1
fi
nimi=$1

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

mode=${mode:-mpi}

# nn=$(( ${#mode}- 4 ))
# if [[ $nn -gt 0 ]]; then
#     nimi=${1}_${mode:$((${#mode}-$nn)):$nn}
# else
#     nimi=$1
# fi


################################
###			                 ###
### output directories       ###
### input files		         ###
###			                 ###
################################

rundir=${outputroot}/${nimi}
datadir=${rundir}/datafiles

# if main directory exists -> clean it
# if not -> make it
if [ $clean == true ]; then
    if [ -d ${rundir} ] ; then
        rm -rf ${rundir}/*
    else
        mkdir -p ${rundir}
    fi
fi    

mkdir -p ${datadir} 

# copy executables and input files to running directory
if [ $COPY == 'true' ]; then
    
    cp ${bin}/les.${mode} ${rundir}/
    cp ${bin}/datafiles/* ${datadir}/
    
    cp ${input}/sound_in ${rundir}/
    cp ${input}/NAMELIST ${rundir}/
fi

################################
###			                 ###
### change output file names ###
###			                 ###
################################

modifyoutput=${modifyoutput:-true}
modifyoutputHistory=${modifyoutputHistory:-${modifyoutput}}

if [ $modifyoutput == 'true' ]; then
    sed -i "/filprf\s\{0,\}=\s\{0,\}/c\  filprf  = '"$nimi"'" ${rundir}/NAMELIST    
fi

if [ $modifyoutputHistory == 'true' ]; then
    sed -i "/hfilin\s\{0,\}=\s\{0,\}/c\  hfilin  = '"$nimi".rst'" ${rundir}/NAMELIST
fi

#########################
###			          ###
### Create run script ###
###		              ###
#########################
echo ' '
## modify the job name based on length: ###
length=$(( ${#nimi} < 6 ? ${#nimi} : 6))
if [[ -n $ownjobname ]]; then
    jobname=LES_$ownjobname
else
    jobname=LES_${nimi:$((${#nimi}-${length})):${length}}
fi    
echo 'Queuing system jobname' $jobname

if [ $JOBFLAG == 'PBS' ] ; then

cat > ${rundir}/runles.sh <<FINALPBS
#!/bin/sh
#PBS -N ${jobname}
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
echo 'Submit to job scheduler'
qsub runles.sh

elif [ $JOBFLAG == 'SBATCH' ] ; then

cat > ${rundir}/runles.sh <<FINALSBATCH
#!/bin/sh
#SBATCH -J ${jobname}
#SBATCH -n ${nproc}
#SBATCH --ntasks-per-node=${nodeNPU}
#SBATCH -t ${WT}
#SBATCH --output=LES_${nimi}-%j.out
#SBATCH --error=LES_${nimi}-%j.err
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
echo 'Submit to job scheduler'
sbatch runles.sh

fi
###########################################    

exit