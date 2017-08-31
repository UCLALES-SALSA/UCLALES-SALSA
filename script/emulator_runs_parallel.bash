#!/bin/bash

# Exit on error
set -e
shopt -s extglob
# import subroutines & variables 
if [[ -d ${SCRIPT} ]]; then
   scriptref=${SCRIPT}
else
   scriptref=.
fi
source ${scriptref}/subroutines_variables.bash


emulatorname=${emulatorname:-case_emulator}
emulatoroutputroot=${outputroot}/${emulatorname}

runNroBegin=${runNroBegin:-1}
simulNro=${simulNro:-1}
runNroEnd=${runNroEnd:-90}
threadNro=${threadNro:-0}

nproc=${nproc:-100}    # number of processors

exe=${exe:-les.mpi}

etunolla=${etunolla:-2}
echo "etunolla" $etunolla 
if [[ -z $list ]]
then
    array=($(seq -f"%0${etunolla}g" ${runNroBegin} ${simulNro} ${runNroEnd}  ))
else
    u=0
    for kk in ${list[@]}
    do
        array[u]=$(printf %0${etunolla}d ${kk##+(0)} )
        u=$((u+1))
    done
fi  
echo "parallel array" ${array[@]}

function submitting {
	
	simulation=$1
    
	echo ' '
	echo ' '
	

	rundir=${emulatoroutputroot}/${simulation}

	LVL=$(grep level ${rundir}/NAMELIST)
	LVL=${LVL: -1} # last char
	echo 'Nyt suoritetaan funktiota submitting simulation: ' $simulation 'level' $LVL
    if [ $LVL -le 3 ]; then
        walltime=12:00:00
    else
        walltime=24:00:00
    fi
	## submit
	
    input=${rundir} outputname=${simulation} modifyoutput='false' COPY=false clean=false WT=$walltime exe=${exe} ${emulatoroutputroot}/submit_uclales-salsa.bash ${emulatorname}/${simulation} $nproc


}


for i in ${array[@]}
do

# sed -i "/runtype\s\{0,\}=\s\{0,\}/c\  runtype  = '"HISTORY"'"  ${outputroot}/${simulation}/NAMELIST
    
    
    status=$( tarkistastatus ${emulatorname}/emul${i} emul${i} )
    LS=${status:0:1}
    PPS=${status:1:2}
    

    start=$( date +%s )
    
    echo ' '
    echo ' '
    echo -n 'Toimitaan emulaattoriajossa ' $i ' '; date '+%T %d-%m-%Y'; echo -n 'status' $status
    echo ' '
    loopN=1 
    while [[ $LS -ne '1' ]] ; do
        if [[ $LS == '2' ]]; then
            sed -i "/runtype\s\{0,\}=\s\{0,\}/c\  runtype  = '"HISTORY"'"  ${outputroot}/${emulatorname}/emul${i}/NAMELIST
        fi
        submitting emul${i}
        odota LESemul${i}
        
        status=$( tarkistastatus ${emulatorname}/emul${i} emul${i} )
        LS=${status:0:1}
        PPS=${status:1:2}
	    echo "while-kierros loppuu: simulaatio" emul${i} "status" $status
        loopN=$((loopN+1))

        if [[ $loopN -gt 4 ]]; then
            printf "There's something wrong with simulation emul${i}, looping for ${loopN}:th time" | mail -s "Something wrong with EMULATOR" jaakko.ahola@fmi.fi,Muzaffer.Ege.Alper@fmi.fi
            echo 'breaking out of while'
            break
        fi
    done
    
    if [[ $loopN -gt 4 ]]; then
        echo 'continue with next item in the for loop'
        continue
    fi
     
    echo ' '
    echo "while paattynyt: simulaatio" emul${i} "status" $status
    end=$( date +%s )
    runtime=$(( end-start ))
    echo -n 'Simulaatio on valmis: ' LES_emul${i}' '; date '+%T %d-%m-%Y'
    echo -n 'Suoritusaika' LES_emul${i}' ' $runtime' '; printf '%02dh:%02dm:%02ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60))
    
    if [ $PPS -eq 0 ] || [ $PPS -eq 3 ] ; then # [ $PPS -eq 2 ] || [ $LS -eq 2 ]    
        scriptname=${scriptname} scriptfolder=${emulatoroutputroot} postprosessoi ${emulatorname}/emul${i} emul${i}
        odota nc_emul${i}
        echo -n 'Postprosessointi on valmis: ' nc_emul${i}' '; date '+%T %d-%m-%Y'
        odota ps_emul${i} 5s
        echo -n 'Postprosessointi on valmis: ' ps_emul${i}' '; date '+%T %d-%m-%Y'
        odota ts_emul${i} 5s
        echo -n 'Postprosessointi on valmis: ' ts_emul${i}' '; date '+%T %d-%m-%Y'
        echo ' '
        echo 'Kaikki on postprosessoitu'
    elif [[ $PPS == '1' ]]; then
        echo 'postprosessointi on jo valmis'
    else
        echo "WARNING Jotain on pielessä restartissa kun postprosessointi status on $PPS"
    fi 
    
    status=$( tarkistastatus ${emulatorname}/emul${i} emul${i} )

    echo "emulaattoriajon emul${i} status: ${status}"

    poistaturhat ${emulatorname}/emul${i} emul${i}
    echo -n 'Turhat poistettu'' '; date '+%T %d-%m-%Y'
    end=$( date +%s )
    runtime=$(( end-start))
    echo 'Valmis' emul$i 
    echo 'Valmis' emul$i >>  ${emulatoroutputroot}/emul${i}/valmis${i}
    echo $runtime >>   ${emulatoroutputroot}/emul${i}/valmis${i}

    
done

echo -n "Valmis threadNro" $threadNro ' '; date '+%T %d-%m-%Y'

echo "LOG Valmis threadNro" $threadNro >> ${emulatoroutputroot}/log
echo -n "LOG " >> ${emulatoroutputroot}/log; date '+%T %d-%m-%Y'   >> ${emulatoroutputroot}/log
date '+%s'            >> ${emulatoroutputroot}/log
echo ' '              >> ${emulatoroutputroot}/log
