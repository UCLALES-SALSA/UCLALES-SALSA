#!/bin/bash

# Submit postprosessing to queue, Jaakko Ahola, FMI, 09/2016
# Update 20.12.2016
#
# input variables:
# $1 = name of netcdf input file
# $2 = name of postprocessing script (OPTIONAL) default value: combine
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
root=/home/users/${username}/${subfolder}
scriptname=combine.py  # postp_uclales-salsa3.py

scriptfolder=${scriptfolder:-}

if [ -z $scriptfolder ]; then
    scriptfolder=${root}/script
fi    

# supercomputer related variable settings
WT=01:00:00 # walltime
nodeNPU=20  # number of processing units in a node  
JOBFLAG=PBS # job flag of the job scheduling system ( e.g. PBS or SBATCH )
postpostfix=pros

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

if [ -z $4 ]; then
  echo "You didn't give the optional and additional postfix"
  echo "Using assumption: " $postpostfix
else
  postpostfix=$4 # 
fi

if [ $JOBFLAG == 'PBS' ]; then
    echo 'using default values of Voima'

elif [ $JOBFLAG == 'SBATCH' ]; then ## CSC's Sisu machine values
    nodeNPU=24 # number of processing units in a node 
    QUEUE=serial # name of the queue of Sisu machine
    root=/homeappl/home/${username}/appl_sisu/${subfolder}
    scriptfolder=${root}/script
fi

################################
###			                 ###
### input directory	         ###
###			                 ###
################################
if [ -z $1 ]; then
  echo "You didn't give any name of netcdf input file"
  exit 1
fi

input=$1
dir=$(dirname ${input})

if [ $dir == '.' ]; then
  dir=$PWD
fi

base=$(basename ${input})
if  [ ${base: -3} == ".ts" ] || [ ${base: -3} == ".ps" ]
then
	postfix=${input:$((${#input}-3)):3}
else
	# "input argument is a regular .nc
	postfix=.nc
fi

echo 'postfix' $postfix
##########################
###			           ###
### Create run scripts ###
###		               ###
##########################
echo " "
## modify the job name based on length: ###
length=$(( ${#postpostfix} < 7 ? ${#postpostfix} : 7))
jobname=${postfix:$((${#postfix}-2)):2}_${postpostfix:$((${#postpostfix}-${length})):${length}}
echo 'Queuing system jobname' $jobname

rm -rf ${dir}/post_* ${dir}/*pros.sh ${dir}/${scriptname}

cp ${scriptfolder}/${scriptname} ${dir}/

### first script
cat > ${dir}/postpros${postfix}.sh <<EOF
#!/bin/bash

set -e

cd ${dir}

python ${scriptname} $input

EOF

### second script

if [ $JOBFLAG == 'PBS' ] ; then

cat > ${dir}/runpostpros${postfix}.sh <<FINALPBS
#!/bin/sh
#PBS -N ${jobname}
#PBS -l mppwidth=1
#PBS -l mppnppn=1
#PBS -l mppdepth=${nodeNPU}
#PBS -l walltime=${WT}
#PBS -j oe
#PBS -M ${email}
#PBS -m ae

source /etc/profile
module load Python

cd ${dir}

aprun -n1 -N1 -d${nodeNPU} ./postpros${postfix}.sh | tee ${PBS_JOBNAME:-post_interactive${postfix}}.${PBS_JOBID:-help}

exit
FINALPBS

cd ${dir}
# Make initial submit
chmod +x runpostpros${postfix}.sh ${scriptname} postpros${postfix}.sh
echo 'Submit to job scheduler'
qsub runpostpros${postfix}.sh

elif [ $JOBFLAG == 'SBATCH' ] ; then

cat > ${dir}/runpostpros${postfix}.sh <<FINALSBATCH
#!/bin/sh
#SBATCH -J ${jobname}
#SBATCH -n 1
#SBATCH -t ${WT}
#SBATCH --output=postpro_${input}-%j.out
#SBATCH --error=postpro_${input}-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${email}
#SBATCH -p ${QUEUE}

source /etc/profile
cd ${dir}

srun -n1 -N1 -d${nodeNPU} ./postpros${postfix}.sh

exit
FINALSBATCH

cd ${dir}
# Make initial submit
chmod +x runpostpros${postfix}.sh ${scriptname} postpros${postfix}.sh
echo 'Submit to job scheduler'
sbatch runpostpros${postfix}.sh

fi

exit
