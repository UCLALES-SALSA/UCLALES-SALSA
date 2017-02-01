#!/bin/bash

# Exit on error
set -e

#
#
# EXAMPLE USAGE
# A=28 B=30 ./submit_emulator_runs_parallel.bash 3
# where A is the case number that you want to start from
#       B is the case number that you want to end up
#    $1 = 3 is the number of how many runs you dare to submit at the same time
#    $2 = number of processors
#    $3 = job flag of the job scheduling system (OPTIONAL) default value: PBS

root=/home/users/aholaj/UCLALES-SALSA
salsa=${root}/src/src_salsa
les=${root}/src/src_LES
bin=${root}/bin
script=${root}/script


#### user parameters
username=aholaj
email=jaakko.ahola@fmi.fi

folder=case_emulator
subfolder=UCLALES-SALSA

inputrootfolder=${bin}/${folder}
outputrootfolder=/lustre/tmp/${username}/${subfolder}/${folder}
# ibrixrootfolder=/ibrix/arch/ClimRes/aholaj/${folder}

# supercomputer related variable settings
WT=24:00:00 # walltime
nodeNPU=20  # number of processing units in a node  
JOBFLAG=PBS     # job flag of the job scheduling system ( e.g. PBS or SBATCH )
scriptname=combine.py
nproc=${nproc:-100}

echo 'poista vanhat kansiot' $folder
rm -rf  ${outputrootfolder} 

if [ -z $1 ]; then
  echo "You didn't the number of simultaneous submits of emulator runs"
  exit 1
fi


Nro=$1

A=${A:-1}
B=${B:-90}
k=$A

if [ -z $2 ]; then
  echo "You didn't give the optional name of the post processing script"
  echo "Using assumption: " ${scriptname} 
else
  scriptname=$2
fi

if [ -z $3 ]; then
  echo "You didn't give the optional JOB FLAG of the job scheduling system"
  echo "Using assumption: " $JOBFLAG
else
  JOBFLAG=$3 # 
fi




if [ $JOBFLAG == 'PBS' ]; then
    echo 'using default values of Voima'
    echo ' '

elif [ $JOBFLAG == 'SBATCH' ]; then ## CSC's Sisu machine values
    nodeNPU=24 # number of processing units in a node 
    QUEUE=small_long # name of the queue of Sisu machine
    outputrootfolder=/wrk/${username}/${subfolder}/${folder}
    root=/homeappl/home/${username}/appl_sisu/${subfolder}
fi

mkdir -p ${outputrootfolder}
##########################################
###                                    ###
### copy all necessary files to lustre ###
###                                    ###
##########################################

cp ${script}/emulator_runs_parallel.bash ${outputrootfolder}/
cp ${script}/submit_uclales-salsa.bash   ${outputrootfolder}/
cp ${script}/submit_postpros.bash   ${outputrootfolder}/

cp ${script}/${scriptname} ${outputrootfolder}/



if [ $nproc -gt 1 ]; then
   mode=mpi
else
   mode=seq
fi

for i in $(seq -f"%02g" $A $B  )
do
    mkdir -p ${outputrootfolder}/emul${i}/datafiles
    cp ${inputrootfolder}/emul${i}/* ${outputrootfolder}/emul${i}/
    cp ${bin}/les.${mode} ${outputrootfolder}/emul${i}/
    cp ${bin}/datafiles/* ${outputrootfolder}/emul${i}/datafiles
done


########################################


if [ $JOBFLAG == 'PBS' ] ; then

cat > ${outputrootfolder}/control_multiple_emulator_run.sh <<FINALPBS
#!/bin/sh
#PBS -N EMULATOR
#PBS -l mppwidth=1
#PBS -l mppnppn=1
#PBS -l mppdepth=${nodeNPU}
#PBS -l walltime=${WT}
#PBS -j oe
#PBS -M ${email}
#PBS -m ae

cd ${outputrootfolder}

FINALPBS

elif [ $JOBFLAG == 'SBATCH' ] ; then

cat > ${outputrootfolder}/control_multiple_emulator_run.sh <<FINALSBATCH
#!/bin/sh
#SBATCH -J EMULATOR
#SBATCH -n 1
#SBATCH -t ${WT}
#SBATCH --output=emulator_${1}-%j.out
#SBATCH --error=emulator_${1}-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${email}
#SBATCH -p ${QUEUE}

cd ${outputrootfolder}

FINALSBATCH

fi


for n in $( seq ${Nro} )
do

echo "runNroBegin=$k simulNro=$Nro runNroEnd=$B threadNro=$n nproc=${nproc} jobflag=$JOBFLAG scriptname=$scriptname ${outputrootfolder}/emulator_runs_parallel.bash | tee ${outputrootfolder}/emulatoroutput${n} &" >> ${outputrootfolder}/control_multiple_emulator_run.sh

k=$((A+$n))
done

cat >> ${outputrootfolder}/control_multiple_emulator_run.sh <<EOF
wait
exit

EOF

cd ${outputrootfolder}

# Make initial submit
chmod +x  control_multiple_emulator_run.sh
echo 'Submit emulator controller to job scheduler'

if [ $JOBFLAG == 'PBS' ] ; then
    qsub control_multiple_emulator_run.sh
elif [ $JOBFLAG == 'SBATCH' ] ; then
    sbatch control_multiple_emulator_run.sh
fi


exit