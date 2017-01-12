#!/bin/bash

# Exit on error
set -e

root=/home/users/aholaj/UCLALES-SALSA/
salsa=${root}/src/src_salsa
les=${root}/src/src_LES
bin=${root}/bin
script=${root}/script


#### user parameters
username=aholaj
folder=case_emulator

inputrootfolder=${bin}/${folder}
outputrootfolder=/lustre/tmp/${username}/UCLALES-SALSA/${folder}
ibrixrootfolder=/ibrix/arch/ClimRes/aholaj/${folder}

runNroBegin=1
runNroEnd=90
nproc=100    # number of processors
jobflag=${jobflag:-PBS}
scriptname=combine.py

if [ -z $postfix ]; then
  echo ' '
  echo "You didn't give any postfix name"
  echo ' '
else
  postfix=_${postfix}
fi


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
    
    echo ' '
    echo 'Nyt odotetaan'
    while [[ ! -z $( qstat -u $username ) ]]
    do
        qstat -u aholaj
#         sleep 15s
        sleep 3m
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
	
	
	nimi=${nimi}${postfix}

	echo 'nimi ' $nimi
	

	
	## submit
	
    subsubfolder=${folder} modifyoutput='false' ${script}/submit_uclales-salsa.bash $nimi $nproc $jobflag
	
	
	sleep 3s

}


function postprossoi {
    
    echo ' '
    nimi=$1
    echo 'Nyt postprosessoidaan .nc'
    ${script}/submit_postpros.bash ${outputrootfolder}/${nimi}/${nimi} $scriptname $jobflag
    
    echo ' '
    echo 'Nyt postprosessoidaan .ps'
    ${script}/submit_postpros.bash ${outputrootfolder}/${nimi}/${nimi}.ps $scriptname $jobflag
    
    echo ' '
    echo 'Nyt postprosessoidaan .nc'
    ${script}/submit_postpros.bash ${outputrootfolder}/${nimi}/${nimi}.ts $scriptname $jobflag
    echo ' '
    
    echo 'Kaikki submittoitu postprosessointiin'
    
    

}

function poistaturhat {
    
    echo ' '
    echo 'Ollaan postprosessoitu ja kopioitu ibrixille, nyt voidaan poistaa lustrelta tiedostot'
    nimi=$1
    rm -rf ${outputrootfolder}/${nimi}/datafiles
    rm -rf ${outputrootfolder}/${nimi}/*.sh
    rm -rf ${outputrootfolder}/${nimi}/*.py
    rm -rf ${outputrootfolder}/${nimi}/les.mpi
    rm -rf ${outputrootfolder}/${nimi}/les.seq
    rm -rf ${outputrootfolder}/${nimi}/*.rst
    rm -rf ${outputrootfolder}/${nimi}/${nimi}.ts.0*0*.nc
    rm -rf ${outputrootfolder}/${nimi}/${nimi}.ps.0*0*.nc
    rm -rf ${outputrootfolder}/${nimi}/${nimi}.0*0*.nc
    rm -rf ${outputrootfolder}/${nimi}/0*_0*.${nimi}.*
    echo ' '

}


function kopioibrixille {
    
    echo ' '
    echo 'Ollaan postprosessoitu, nyt voidaan kopioida ibrixille'
    nimi=$1
    rsync -avz ${outputrootfolder}/${nimi}/ ${ibrixrootfolder}/${nimi}/
    echo ' '

}


function poistalustrelta {
    
    echo ' '
    echo 'Ollaan postprosessoitu ja kopioitu ibrixille, nyt voidaan poistaa lustrelta tiedostot'
    nimi=$1
    rm -rf ${outputrootfolder}/${nimi}
    echo ' '

}




for i in $(seq -f"%02g" ${runNroBegin} ${runNroEnd}  )
do
    echo ' '
    echo ' '
    echo ' '
    echo ' '
    echo ' '
    echo $i
    echo ' '
    submitting emul${i} 
    odota
    
    postprossoi emul${i}
    odota
    
    
    poistaturhat emul${i}
    kopioibrixille emul${i}
    poistalustrelta emul${i}
done


