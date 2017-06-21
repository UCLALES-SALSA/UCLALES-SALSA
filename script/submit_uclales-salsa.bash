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

if [[ -d ${SCRIPT} ]]; then
   scriptref=${SCRIPT}
else
   scriptref=.
fi
source ${scriptref}/subroutines_variables.bash

# supercomputer related variable settings
WT=${WT:-24:00:00} # walltime


COPY=${COPY:-true}
clean=${clean:-true}

#################################
###			                  ###
### folders		              ###
###         			      ### 
#################################

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
  echo "You didn't give any name for output subfolder of lustre"
  exit 1
fi
simulation=$1

# if outputname doesn't exist set it to be same as simulation
if [[ -z $outputname ]]; then
    outputname=${simulation}
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

if [[ -n $3 ]]; then
  jobflag=$3
fi
echo "job scheduling system" $jobflag

################################
###			                 ###
### sequential or mpi        ###
###			                 ###
################################

mode=${mode:-mpi}

# nn=$(( ${#mode}- 4 ))
# if [[ $nn -gt 0 ]]; then
#     outputname=${1}_${mode:$((${#mode}-$nn)):$nn}
# else
#     outputname=$1
# fi


################################
###			                 ###
### output directories       ###
### input files		         ###
###			                 ###
################################

rundir=${outputroot}/${simulation}
datadir=${rundir}/datafiles

# if main directory exists -> clean it
# if not -> make it
if [ $clean == true ]; then
    if [ -d ${rundir} ] ; then
        rm -rf ${rundir}/*
    else
        mkdir -p ${rundir} ${datadir}
    fi
fi    


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
    sed -i "/filprf\s\{0,\}=\s\{0,\}/c\  filprf  = '"$outputname"'" ${rundir}/NAMELIST    
fi

if [ $modifyoutputHistory == 'true' ]; then
    sed -i "/hfilin\s\{0,\}=\s\{0,\}/c\  hfilin  = '"$outputname".rst'" ${rundir}/NAMELIST
fi

#########################
###			          ###
### Create run script ###
###		              ###
#########################
echo ' '
## modify the job name based on length: ###
if [[ -n $ownjobname ]]; then
    apuNimi=$ownjobname
else
    apuNimi=$outputname
fi    
length=$(( ${#apuNimi} < 7 ? ${#apuNimi} : 7))
jobname=LES${apuNimi:$((${#apuNimi}-${length})):${length}}
echo 'Queuing system jobname' $jobname

if [ $jobflag == 'PBS' ] ; then

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

elif [ $jobflag == 'SBATCH' ] ; then

cat > ${rundir}/runles.sh <<FINALSBATCH
#!/bin/sh
#SBATCH -J ${jobname}
#SBATCH -n ${nproc}
#SBATCH --ntasks-per-node=${nodeNPU}
#SBATCH -t ${WT}
#SBATCH --output=LES_${outputname}-%j.out
#SBATCH --error=LES_${outputname}-%j.err
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
