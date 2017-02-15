#!/bin/bash
start=$(date +%s)

#### user parameters
username=aholaj
email=jaakko.ahola@fmi.fi

folder=${folder:-case_emulator}
subfolder=UCLALES-SALSA

outputrootfolder=/lustre/tmp/${username}/${subfolder}/${folder}
ibrixrootfolder=/ibrix/arch/ClimRes/aholaj/${folder}

runNroBegin=${runNroBegin:-1}
simulNro=${simulNro:-1}
runNroEnd=${runNroEnd:-90}

# supercomputer related variable settings
WT=24:00:00 # walltime
WTs=$(echo $WT | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }')
poistavanhat=${poistavanhat:-false}

################################
###			                 ###
### output directory	     ###
###			                 ###
################################

if [ -z $1 ]; then
    echo "You didn't give any postfix"
else
    postfix=_$1
fi

################################
###			                 ###
### number of processors     ###
###			                 ###
################################
if [ -z $2 ]; then
  echo "You didn't give time limit"
  echo "Use default value" $WT
else
  WT=$2
fi

if [ $poistavanhat == 'true' ]; then
    echo 'Poistetaan vanhat'
    rm -rf ${ibrixrootfolder}${postfix}
fi

function kopioibrixille {
    
    echo ' '
    echo 'Kopioidaan ibrixille'
    nimi=$1
    mkdir -p ${ibrixrootfolder}${postfix}/${nimi}/
    rsync -avz ${outputrootfolder}/${nimi}/ ${ibrixrootfolder}${postfix}/${nimi}/
    echo ' '

}


function poistalustrelta {
    
    echo ' '
    echo 'Poistetaan lustrelta tiedostot'
    nimi=$1
    rm -rf ${outputrootfolder}/${nimi}
    echo ' '

}
timeleft=$WTs
usedtime=0
while [ $timeleft -ge 0 -a $(ls -l ${outputrootfolder} | grep "^d" | awk -F" " '{print $9}' | wc -l) -gt 0 ]
do
    echo -n 'loopissa ollaan '; date '+%T %d-%m-%Y'
    for i in $(seq -f"%02g" ${runNroBegin} ${runNroEnd}  )
    do
        if [ -e ${outputrootfolder}/emul${i}/valmis${i} ]; then
            echo ' '
            echo 'TADAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
            echo 'tiedosto' valmis$i 'on olemassa'
            echo ' '
            kopioibrixille emul${i}
            poistalustrelta emul${i}
        fi
    done
    sleep 5s
    now=$(date +%s)
    usedtime=$((now-start))
    timeleft=$((WTs-usedtime))
    
done

kopioibrixille
poistalustrelta