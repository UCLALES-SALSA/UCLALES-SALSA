#!/bin/bash

# user parameters

# folders
model=UCLALES-SALSA



email='jaakko.ahola@fmi.fi'

jobflag=PBS
scriptname=combine.py

ibrixrootfolder=/ibrix/arch/ClimRes/${USER}/
# subroutines

if [ $jobflag == 'PBS' ]; then
    outputroot=/lustre/tmp/${USER}/${model}/
    root=/home/users/${USER}/${model}
    nodeNPU=20  # number of processing units in a node  

elif [ $jobflag == 'SBATCH' ]; then ## CSC's Sisu machine values
    nodeNPU=24 # number of processing units in a node 
    QUEUE=small_long # name of the queue of Sisu machine
    outputroot=/wrk/${USER}/${model}
    root=/homeappl/home/${USER}/appl_sisu/${model}
fi

salsa=${root}/src/src_salsa
les=${root}/src/src_LES
bin=${root}/bin
script=${root}/script


function odota {
    
    outputname=$1
    aika=$2
    aika=${aika:-30s}
    
    echo ' '
    echo 'Nyt odotetaan' $outputname $aika
    while [[ ! -z $( qstat -u $USER | grep $outputname ) ]]
    do
        date +%Y-%m-%d-%H-%M
        qstat -u $USER
        sleep $aika
    done

}

function kopioibrixille {
    
    echo ' '
    echo 'Kopioidaan ibrixille'
    simulation=$1
    mkdir -p ${ibrixrootfolder}/${simulation}/
    rsync -avz ${outputroot}/${simulation}/ ${ibrixrootfolder}/${simulation}/
    echo ' '

}

function postprosessoi {
    
    echo ' '
    simulation=$1
    outputname=$2
    jobnamepostfix=$3
    
    # if outputname doesn't exist set it to be same as simulation
    if [[ -z $outputname ]]; then
        outputname=${simulation} 
    fi
    
    # if scriptfolder doesn't exist set it to be the default
    if [[ -z $scriptfolder ]]; then
        scriptfolder=${script}
    fi
    
    # if jobnamepostfix doesn't exist set it to be the default
    if [[ -z $jobnamepostfix ]]; then
        jobnamepostfix=${outputname}
    fi
    
    
    echo 'Nyt postprosessoidaan nc'
    ${scriptfolder}/submit_postpros.bash ${outputroot}/${simulation}/${outputname}    $scriptname $jobflag $jobnamepostfix
    
    echo ' '
    echo 'Nyt postprosessoidaan ps'
    ${scriptfolder}/submit_postpros.bash ${outputroot}/${simulation}/${outputname}.ps $scriptname $jobflag $jobnamepostfix
    
    echo ' '
    echo 'Nyt postprosessoidaan ts'
    ${scriptfolder}/submit_postpros.bash ${outputroot}/${simulation}/${outputname}.ts $scriptname $jobflag $jobnamepostfix
    echo ' '
    
    echo 'Kaikki submittoitu postprosessointiin'

}

function poistaturhat {
    
    echo ' '
    simulation=$1
    outputname=$2
    restart=${restart:-false}
    if [[ -z $outputname ]]; then
        outputname=${simulation}
    fi
    
    echo 'Poistetaan lustrelta turhat tiedostot jos postprosessointi on tehty, kansio :' ${outputroot}/${simulation}/
    
    if [ $restart == 'true' ] ||  ([ $restart == 'false' ] && [ -f ${outputroot}/${simulation}/${outputname}.nc ] && [ -f ${outputroot}/${simulation}/${outputname}.ts.nc ] && [ -f ${outputroot}/${simulation}/${outputname}.ps.nc ])
    then
        
        if [[ $restart == 'false' ]]
        then
        
            echo "kaikki kolme postprosessoitua tiedostoa ovat olemassa"
            basename ${outputroot}/${simulation}/${outputname}.nc
            basename ${outputroot}/${simulation}/${outputname}.ts.nc
            basename ${outputroot}/${simulation}/${outputname}.ps.nc
        else
            echo "restartataan, poistetaan mahdolliset postprosessoidut setit"
            rm -rf ${outputroot}/${simulation}/*.nc
        fi
        
        echo 'Poistetaan'
        rm -rf ${outputroot}/${simulation}/datafiles
        rm -rf ${outputroot}/${simulation}/*.sh
        rm -rf ${outputroot}/${simulation}/*.py
        rm -rf ${outputroot}/${simulation}/*.rst
        rm -rf ${outputroot}/${simulation}/${outputname}.ts.0*0*.nc
        rm -rf ${outputroot}/${simulation}/${outputname}.ps.0*0*.nc
        rm -rf ${outputroot}/${simulation}/${outputname}.0*0*.nc
        rm -rf ${outputroot}/${simulation}/0*_0*.${outputname}.*
    else
        echo 'simulation' $simulation "ei ole valmis ei poisteta turhia"
    fi
}
