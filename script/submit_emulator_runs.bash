#!/bin/bash

# Exit on error
set -e

# import subroutines & variables 
source ./subroutines_variables.bash

#
#
# EXAMPLE USAGE
# A=28 B=30 mode=mpi.v1.0.4 postfix=geggals ./submit_emulator_runs_parallel.bash 3
# where A is the case number that you want to start from
#       B is the case number that you want to end up
#    $1 = 3 is the number of how many runs you dare to submit at the same time
#    $2 = scriptname
#    $3 = job flag of the job scheduling system (OPTIONAL) default value: PBS
#
# if you want to run specific cases given them as list="01 02 09 25 63 82"


# essential input values
submit=${submit:-true}
restart=${restart:-false}

A=${A:-1}
B=${B:-90}


if [[ -z $list ]]
then
    array=($(seq -f"%02g" $A $B  ))
else
    u=0
    for kk in ${list[@]}
    do
        array[u]=$(printf %02d $kk)
        u=$((u+1))
    done
fi
echo "submit array" ${array[@]}
k=0

mode=${mode:-mpi}
# supercomputer related variable settings
WT=336:00:00 # walltime
nproc=${nproc:-100}


bashfolder="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [[ -n $postfix ]]; then
    echo "postfix annettu" $postfix
    postfix=_${postfix}
fi

# design version
cd $(dirname $DESIGN)
designV=$(git describe --tags)
cd $bashfolder

# LES version
nn=$(( ${#mode}- 4 ))
if [[ $nn -gt 0 ]]; then
    les=_${mode:$((${#mode}-$nn)):$nn}
fi

emulatorinput=case_emulator_DESIGN_${designV}
inputrootfolder=${bin}/${emulatorinput}


LVL=$(grep level ${inputrootfolder}/emul??/NAMELIST | tail -1)
LVL=${LVL: -1} # last char

emulatorname=${emulatorinput}_LES${les}_LVL${LVL}${postfix}

emulatoroutputroot=${outputroot}/${emulatorname}


if [[ $restart == 'false' ]]; then
    echo 'poista vanhat kansiot' $folder
    rm -rf  ${emulatoroutputroot} 
fi

if [[ -z $1 ]]; then
  echo "You didn't give the number of simultaneous submits of emulator runs  (= number of THREADS)"
  exit 1
fi
ThreadNro=$1


if [[ -n $2 ]]; then
  scriptname=$2
fi
echo "scriptname" ${scriptname}

if [[ -n $3 ]]; then
  jobflag=$3
fi
echo "job scheduling system" $jobflag


mkdir -p ${emulatoroutputroot}
echo 'emulator begin' >> ${emulatoroutputroot}/log
date '+%T %d-%m-%Y'   >> ${emulatoroutputroot}/log
date '+%s'            >> ${emulatoroutputroot}/log
echo ' '              >> ${emulatoroutputroot}/log
##########################################
###                                    ###
### copy all necessary files to lustre ###
###                                    ###
##########################################

cp ${script}/subroutines_variables.bash ${emulatoroutputroot}/
cp ${script}/emulator_runs_parallel.bash ${emulatoroutputroot}/
cp ${script}/submit_uclales-salsa.bash   ${emulatoroutputroot}/
cp ${script}/submit_postpros.bash        ${emulatoroutputroot}/

cp ${script}/${scriptname}               ${emulatoroutputroot}/

for i in ${array[@]}
do
    if [ $restart == 'false' ] ||  ([ $restart == 'true' ] && ([ ! -f ${emulatoroutputroot}/emul${i}/emul${i}.nc ] || [ ! -f ${emulatoroutputroot}/emul${i}/emul${i}.ts.nc ] || [ ! -f ${emulatoroutputroot}/emul${i}/emul${i}.ps.nc ] || [ ! -f ${emulatoroutputroot}/emul${i}/valmis${i} ]));
    then
        echo emul${i} "ei ole valmis"
        mkdir -p ${emulatoroutputroot}/emul${i}/
        if [[ -z $(ls -A ${emulatoroutputroot}/emul${i}/) ]]
        then
            echo 'folder is empty'
            
            cp ${inputrootfolder}/emul${i}/* ${emulatoroutputroot}/emul${i}/
            cp ${bin}/les.${mode} ${emulatoroutputroot}/emul${i}/
        else
            echo "poistetaan turhat"
            restart=$restart poistaturhat ${emulatorinput}/emul${i} emul${i}
        fi
        
        
        mkdir -p ${emulatoroutputroot}/emul${i}/datafiles
        cp ${bin}/datafiles/* ${emulatoroutputroot}/emul${i}/datafiles
    fi
done


########################################


if [ $jobflag == 'PBS' ] ; then

cat > ${emulatoroutputroot}/control_multiple_emulator_run.sh <<FINALPBS
#!/bin/sh
#PBS -N EMULATOR
#PBS -l mppwidth=1
#PBS -l mppnppn=1
#PBS -l mppdepth=${nodeNPU}
#PBS -l walltime=${WT}
#PBS -j oe
#PBS -M ${email}
#PBS -m ae

cd ${emulatoroutputroot}

FINALPBS

elif [ $jobflag == 'SBATCH' ] ; then

cat > ${emulatoroutputroot}/control_multiple_emulator_run.sh <<FINALSBATCH
#!/bin/sh
#SBATCH -J EMULATOR
#SBATCH -n 1
#SBATCH -t ${WT}
#SBATCH --output=emulator_${1}-%j.out
#SBATCH --error=emulator_${1}-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${email}
#SBATCH -p ${QUEUE}

cd ${emulatoroutputroot}

FINALSBATCH

fi

# echo ' '
aloitusindeksi=0
nroJobs=0
for n in $( seq 0 $((ThreadNro-1)) )
do
nroJobs=$(python -c "from math import ceil; print int( ceil( ( ${#array[@]}-$aloitusindeksi  )/float( $ThreadNro-$n ) ) )")
# echo apu $apu bang $bang
echo "submit->parallel" ${array[@]:$aloitusindeksi:$nroJobs}
echo "emulatorname=${emulatorname} list='"${array[@]:$aloitusindeksi:$nroJobs}"' threadNro=$n nproc=${nproc} jobflag=$jobflag mode=${mode} restart=${restart} scriptname=$scriptname ${emulatoroutputroot}/emulator_runs_parallel.bash | tee ${emulatoroutputroot}/emulatoroutput$(printf %02d $n) &" >> ${emulatoroutputroot}/control_multiple_emulator_run.sh
aloitusindeksi=$((aloitusindeksi+nroJobs))
done

cat >> ${emulatoroutputroot}/control_multiple_emulator_run.sh <<EOF
wait
exit

EOF

cd ${emulatoroutputroot}

# Make initial submit
chmod +x  control_multiple_emulator_run.sh

if [ $submit == 'true' ]; then
    echo 'Submit emulator controller to job scheduler'
    if [ $jobflag == 'PBS' ] ; then
        qsub control_multiple_emulator_run.sh
        qstat -u $USER
    elif [ $jobflag == 'SBATCH' ] ; then
        sbatch control_multiple_emulator_run.sh
    fi
else
    echo "NOT submitting"
fi

exit