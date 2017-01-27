#!/bin/bash

# Exit on error
set -e

root=/home/users/aholaj/UCLALES-SALSA
salsa=${root}/src/src_salsa
les=${root}/src/src_LES
bin=${root}/bin
script=${root}/script


#### user parameters
username=aholaj
folder=case_emulator

inputrootfolder=${bin}/${folder}
outputrootfolder=/lustre/tmp/${username}/UCLALES-SALSA/${folder}
script=${root}/script
ibrixrootfolder=/ibrix/arch/ClimRes/aholaj/${folder}

runNroBegin=${runNroBegin:-1}
simulNro=${simulNro:-1}
runNroEnd=${runNroEnd:-90}
threadNro=${threadNro:-0}
nproc=100    # number of processors
jobflag=${jobflag:-PBS}
scriptname=combine.py




compile=${compile:-false}

if [ $compile == 'true' ]; then
    cd ${root}
    make seq || exit 1
    make mpi || exit 1
    cd ${script}
    echo " "
    echo "Compiled"
elif [ $compile == 'false' ]; then
    echo 'use previously compiled les-binary'
    compile=$2
fi


function odota {
    
    nimi=$1
    aika=$2
    aika=${aika:-30s}
    
    echo ' '
    echo 'Nyt odotetaan' $nimi $aika
    while [[ ! -z $( qstat -u $username | grep $nimi ) ]]
    do
        qstat -u aholaj
#         sleep 15s
        sleep $aika
    done

}

function submitting {
	
	nimi=$1
    
	echo ' '
	echo ' '
	echo 'Nyt suoritetaan funktiota submitting '

	dir=${inputrootfolder}/${nimi}
	
	cp ${dir}/sound_in ${bin}/sound_in
	cp ${dir}/NAMELIST ${bin}/NAMELIST


	echo 'nimi ' $nimi
	

	
	## submit
	
    subsubfolder=${folder} modifyoutput='false' ${script}/submit_uclales-salsa.bash $nimi $nproc $jobflag
	
	
	sleep 3s

}


function postprosessoi {
    
    echo ' '
    nimi=$1
    echo 'Nyt postprosessoidaan .nc'
    ${script}/submit_postpros.bash ${outputrootfolder}/${nimi}/${nimi} $scriptname $jobflag $nimi
    
    echo ' '
    echo 'Nyt postprosessoidaan .ps'
    ${script}/submit_postpros.bash ${outputrootfolder}/${nimi}/${nimi}.ps $scriptname $jobflag $nimi
    
    echo ' '
    echo 'Nyt postprosessoidaan .ts'
    ${script}/submit_postpros.bash ${outputrootfolder}/${nimi}/${nimi}.ts $scriptname $jobflag $nimi
    echo ' '
    
    echo 'Kaikki submittoitu postprosessointiin'
    
    

}

function poistaturhat {
    
    echo ' '
    echo 'Poistetaan lustrelta turhat tiedostot'
    nimi=$1
    rm -rf ${outputrootfolder}/${nimi}/datafiles
    rm -rf ${outputrootfolder}/${nimi}/*.sh
    rm -rf ${outputrootfolder}/${nimi}/*.py
    rm -rf ${outputrootfolder}/${nimi}/les.*
    rm -rf ${outputrootfolder}/${nimi}/*.rst
    rm -rf ${outputrootfolder}/${nimi}/${nimi}.ts.0*0*.nc
    rm -rf ${outputrootfolder}/${nimi}/${nimi}.ps.0*0*.nc
    rm -rf ${outputrootfolder}/${nimi}/${nimi}.0*0*.nc
    rm -rf ${outputrootfolder}/${nimi}/0*_0*.${nimi}.*
    echo ' '

}


function kopioibrixille {
    
    echo ' '
    echo 'Kopioidaan ibrixille'
    nimi=$1
    mkdir -p ${ibrixrootfolder}/${nimi}/
    rsync -avz ${outputrootfolder}/${nimi}/ ${ibrixrootfolder}/${nimi}/
    echo ' '

}


function poistalustrelta {
    
    echo ' '
    echo 'Poistetaan lustrelta tiedostot'
    nimi=$1
    rm -rf ${outputrootfolder}/${nimi}
    echo ' '

}




for i in $(seq -f"%02g" ${runNroBegin} ${simulNro} ${runNroEnd}  )
do
    echo ' '
    echo ' '
    echo ' '
    echo ' '
    echo ' '
    echo 'Käynnistetään emulaattoriajo ' $i; date '+%T %d-%m-%Y'
    echo ' '
    submitting emul${i} 
    odota LES_emul${i}
    echo 'Simulaatio on valmis: ' LES_emul${i}; date '+%T %d-%m-%Y'
    
    
    postprosessoi emul${i}
    odota nc_emul${i}
    echo 'Postprosessointi on valmis: ' nc_emul${i}; date '+%T %d-%m-%Y'
    odota ps_emul${i} 5s
    echo 'Postprosessointi on valmis: ' ps_emul${i}; date '+%T %d-%m-%Y'
    odota ts_emul${i} 5s
    echo 'Postprosessointi on valmis: ' ts_emul${i}; date '+%T %d-%m-%Y'
    echo ' '
    echo 'Kaikki on postprosessoitu'
    
    poistaturhat emul${i}
    echo 'Turhat poistettu'; date '+%T %d-%m-%Y'
    kopioibrixille emul${i}
    echo 'kopioitu ibrixille'; date '+%T %d-%m-%Y'
    poistalustrelta emul${i}
    echo 'lustre tyhjennetty'; date '+%T %d-%m-%Y'
done

echo "Valmis threadNro" $threadNro 
