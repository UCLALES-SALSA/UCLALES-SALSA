#!/bin/bash

# Exit on error
set -e

# import subroutines & variables 
source ./subroutines_variables.bash


emulatorname=${emulatorname:-case_emulator}
emulatoroutputroot=${outputroot}/${emulatorname}

runNroBegin=${runNroBegin:-1}
simulNro=${simulNro:-1}
runNroEnd=${runNroEnd:-90}
threadNro=${threadNro:-0}

nproc=${nproc:-100}    # number of processors

mode=${mode:-mpi}

restart=${restart:-false}

if [[ -z $list ]]
then
    array=($(seq -f"%02g" ${runNroBegin} ${simulNro} ${runNroEnd}  ))
else
    u=0
    for kk in ${list[@]}
    do
        array[u]=$(printf %02d $kk)
        u=$((u+1))
    done
fi  
echo "parallel array" ${array[@]}

function submitting {
	
	simulation=$1
    
	echo ' '
	echo ' '
	echo 'Nyt suoritetaan funktiota submitting '

	rundir=${emulatoroutputroot}/${simulation}

	echo 'simulation ' $simulation
	LVL=$(grep level ${rundir}/NAMELIST)
	LVL=${LVL: -1} # last char
	echo 'level' $LVL
    if [ $LVL -le 3 ]; then
        walltime=12:00:00
    else
        walltime=36:00:00
    fi
	## submit
	
    input=${rundir} outputname=${simulation} modifyoutput='false' COPY=false clean=false WT=$walltime mode=${mode} ${emulatoroutputroot}/submit_uclales-salsa.bash ${emulatorname}/${simulation} $nproc


}


for i in ${array[@]}
do
    if [ $restart == 'false' ] ||  ([ $restart == 'true' ] && ([ ! -f ${emulatoroutputroot}/emul${i}/emul${i}.nc ] || [ ! -f ${emulatoroutputroot}/emul${i}/emul${i}.ts.nc ] || [ ! -f ${emulatoroutputroot}/emul${i}/emul${i}.ps.nc ]))
    then

        echo ' '
        echo ' '
        echo -n 'Käynnistetään emulaattoriajo ' $i ' '; date '+%T %d-%m-%Y'
        echo ' '
        submitting emul${i} 
        odota LES_emul${i}
        echo -n 'Simulaatio on valmis: ' LES_emul${i}' '; date '+%T %d-%m-%Y'
        
        
        scriptname=${scriptname} scriptfolder=${emulatoroutputroot} postprosessoi ${emulatorname}/emul${i} emul${i}

        odota nc_emul${i}
        echo -n 'Postprosessointi on valmis: ' nc_emul${i}' '; date '+%T %d-%m-%Y'
        odota ps_emul${i} 5s
        echo -n 'Postprosessointi on valmis: ' ps_emul${i}' '; date '+%T %d-%m-%Y'
        odota ts_emul${i} 5s
        echo -n 'Postprosessointi on valmis: ' ts_emul${i}' '; date '+%T %d-%m-%Y'
        echo ' '
        echo 'Kaikki on postprosessoitu'
        
        poistaturhat ${emulatorname}/emul${i} emul${i}
        echo -n 'Turhat poistettu'' '; date '+%T %d-%m-%Y'
        echo 'Valmis' emul$i 
        echo 'Valmis' emul$i >>  ${emulatoroutputroot}/emul${i}/valmis${i}
        date '+%T %d-%m-%Y' >>   ${emulatoroutputroot}/emul${i}/valmis${i}
    fi
done

echo -n "Valmis threadNro" $threadNro ' '; date '+%T %d-%m-%Y'

echo "Valmis threadNro" $threadNro >> ${emulatoroutputroot}/log
date '+%T %d-%m-%Y'   >> ${emulatoroutputroot}/log
date '+%s'            >> ${emulatoroutputroot}/log
echo ' '              >> ${emulatoroutputroot}/log