#!/bin/bash

# Exit on error
set -e


#### user parameters
username=aholaj
folder=case_emulator
subfolder=UCLALES-SALSA


root=/lustre/tmp/${username}/${subfolder}/${folder}



runNroBegin=${runNroBegin:-1}
simulNro=${simulNro:-1}
runNroEnd=${runNroEnd:-90}
threadNro=${threadNro:-0}
nproc=${nproc:-100}    # number of processors
jobflag=${jobflag:-PBS}
scriptname=${scriptname:-combine.py}







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

	dir=${root}/${nimi}

	echo 'nimi ' $nimi
	

	
	## submit
	
    input=${dir} subsubfolder=${folder} modifyoutput='false' COPY=false clean=false ${root}/submit_uclales-salsa.bash $nimi $nproc $jobflag


}


function postprosessoi {
    
    echo ' '
    nimi=$1
    echo 'Nyt postprosessoidaan .nc'
    scriptfolder=${root} ${root}/submit_postpros.bash ${root}/${nimi}/${nimi} $scriptname $jobflag $nimi
    
    echo ' '
    echo 'Nyt postprosessoidaan .ps'
    scriptfolder=${root} ${root}/submit_postpros.bash ${root}/${nimi}/${nimi}.ps $scriptname $jobflag $nimi
    
    echo ' '
    echo 'Nyt postprosessoidaan .ts'
    scriptfolder=${root} ${root}/submit_postpros.bash ${root}/${nimi}/${nimi}.ts $scriptname $jobflag $nimi
    echo ' '
    
    echo 'Kaikki submittoitu postprosessointiin'
    
    

}

function poistaturhat {
    
    echo ' '
    echo 'Poistetaan lustrelta turhat tiedostot'
    nimi=$1
    rm -rf ${root}/${nimi}/datafiles
    rm -rf ${root}/${nimi}/*.sh
    rm -rf ${root}/${nimi}/*.py
    rm -rf ${root}/${nimi}/les.*
    rm -rf ${root}/${nimi}/*.rst
    rm -rf ${root}/${nimi}/${nimi}.ts.0*0*.nc
    rm -rf ${root}/${nimi}/${nimi}.ps.0*0*.nc
    rm -rf ${root}/${nimi}/${nimi}.0*0*.nc
    rm -rf ${root}/${nimi}/0*_0*.${nimi}.*
    echo ' '

}







for i in $(seq -f"%02g" ${runNroBegin} ${simulNro} ${runNroEnd}  )
do
    echo ' '
    echo ' '
    echo ' '
    echo ' '
    echo ' '
    echo -n 'Käynnistetään emulaattoriajo ' $i; date '+%T %d-%m-%Y'
    echo ' '
    submitting emul${i} 
    odota LES_emul${i}
    echo -n 'Simulaatio on valmis: ' LES_emul${i}; date '+%T %d-%m-%Y'
    
    
    postprosessoi emul${i}
    odota nc_emul${i}
    echo -n 'Postprosessointi on valmis: ' nc_emul${i}; date '+%T %d-%m-%Y'
    odota ps_emul${i} 5s
    echo -n 'Postprosessointi on valmis: ' ps_emul${i}; date '+%T %d-%m-%Y'
    odota ts_emul${i} 5s
    echo -n 'Postprosessointi on valmis: ' ts_emul${i}; date '+%T %d-%m-%Y'
    echo ' '
    echo 'Kaikki on postprosessoitu'
    
    poistaturhat emul${i}
    echo -n 'Turhat poistettu'; date '+%T %d-%m-%Y'
    echo 'Valmis' 'emul'$i 
    echo 'Valmis' 'emul'$i >  ${root}/emul${i}/valmis${i}
    date '+%T %d-%m-%Y' >  ${root}/emul${i}/valmis${i}
done

echo "Valmis threadNro" $threadNro 
