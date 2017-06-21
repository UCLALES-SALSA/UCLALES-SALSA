#!/bin/bash

# user parameters

# folders
model=UCLALES-SALSA



email='jaakko.ahola@fmi.fi'

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
    tulosta=${tulosta:-false}
    echo ' '
    echo 'Nyt odotetaan' $outputname $aika
    while [[ ! -z $( qstat -u $USER | grep $outputname ) ]]
    do
        if [[ $tulosta == 'true' ]]; then
            date +%Y-%m-%d-%H-%M
            qstat -u $USER | grep $outputname
        fi
        sleep $aika
    done

}

function kopioibrixille {
    
    echo ' '
    echo 'Kopioidaan ibrixille'
    simulation=$1
    outputname=$2
    
    # if outputname doesn't exist set it to be same as simulation
    if [[ -z $outputname ]]; then
        outputname=${simulation} 
    fi
    
    mkdir -p ${ibrixrootfolder}/${outputname}/
    rsync -avz ${outputroot}/${simulation}/ ${ibrixrootfolder}/${outputname}/
    echo ' '

}

function postprosessoi {
    
    echo ' '
    simulation=$1
    outputname=$2
    jobnamepostfix=$3
    
    ncsub=${ncsub:-true}
    pssub=${pssub:-true}
    tssub=${tssub:-true}
    
    
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
    
    if [[ $ncsub == "true" ]]; then
        echo ' '
        echo 'Nyt postprosessoidaan nc'
        ${scriptfolder}/submit_postpros.bash ${outputroot}/${simulation}/${outputname}    $scriptname $jobflag $jobnamepostfix
    fi
    
    if [[ $pssub == "true" ]]; then
        echo ' '
        echo 'Nyt postprosessoidaan ps'
        ${scriptfolder}/submit_postpros.bash ${outputroot}/${simulation}/${outputname}.ps $scriptname $jobflag $jobnamepostfix
    fi
    
    if [[ $tssub == "true" ]]; then    
        echo ' '
        echo 'Nyt postprosessoidaan ts'
        ${scriptfolder}/submit_postpros.bash ${outputroot}/${simulation}/${outputname}.ts $scriptname $jobflag $jobnamepostfix
    fi
    
    echo ' '
    echo 'Submittoitu postprosessointiin'

}

function poistaturhat {
    

    simulation=$1
    outputname=$2

    if [[ -z $outputname ]]; then
        outputname=${simulation}
    fi
    
    echo 'Tarkastellaan statusta, kansio :' ${simulation}
    
    status=$( tarkistastatus $simulation $outputname )
    LS=${status:0:1}
    PPS=${status:1:2}

    if [[ $status == '11' ]]
    then
        echo 'Simulaatio' $simulation  'on valmis poistetaan turhat'
#       rm -rf ${outputroot}/${simulation}/datafiles
        rm -rf ${outputroot}/${simulation}/*.sh
        rm -rf ${outputroot}/${simulation}/*.py
        rm -rf ${outputroot}/${simulation}/*.rst
        rm -rf ${outputroot}/${simulation}/${outputname}.ts.0*0*.nc
        rm -rf ${outputroot}/${simulation}/${outputname}.ps.0*0*.nc
        rm -rf ${outputroot}/${simulation}/${outputname}.0*0*.nc
        rm -rf ${outputroot}/${simulation}/0*_0*.${outputname}.*
    elif [ $PPS -eq 2 ] || [ $LS -eq 2 ]
    then
        echo "Restartataan, koska jokin simulaation vaiheista on kesken, poistetaan mahdolliset postprosessoidut setit"
        rm -rf ${outputroot}/${simulation}/${outputname}.nc ${outputroot}/${simulation}/${outputname}.ts.nc  ${outputroot}/${simulation}/${outputname}.ps.nc

    else
        echo 'Simulaatio' $simulation "on alkutilassa ei poisteta turhia"
    fi
    echo $outputname 'status' $status
}

# this function returns following values according to status of simulation and post-prosessing
#    00 (LES not ready  POSTPROS not ready )
#    11 (LES     ready  POSTPROS     ready )
#    22 (LES incomplete POSTPROS incomplete)
function tarkistastatus {
    simulation=$1
    outputname=$2
    
    if [[ -z $outputname ]]; then
        outputname=${simulation}
    fi
    #################
    #### LESin STATUS
    #################
    for f in ${outputroot}/${simulation}/NAMELIST; do
        [ -e "$f" ] && timmax=$( cat ${outputroot}/${simulation}/NAMELIST | grep timmax | cut -c11-30 | tr -d .) || timmax=100000000
        break
    done
    
    for f in ${outputroot}/${simulation}/LES*; do
        [ -e "$f" ] && last=$( cat "$(ls -rt ${outputroot}/${simulation}/LES* | tail -n1)" | grep --ignore-case "model time" | tail -1 | cut -c40-46 ) || last=0
        break
    done
    
    if [[ $last -ge $((timmax-1)) ]]; then
        for f in ${outputroot}/${simulation}/*.nc; do 
            [ -e "$f" ] && les=1 || les=0
            break
        done
    else
        for f in ${outputroot}/${simulation}/*.rst; do # les ei ole valmis, voidaan restartata jos tiedostot on olemassa
            [ -e "$f" ] && les=2 || les=0
            break
        done
    fi
    
    #touch ${outputroot}/${simulation}/debug${outputname}.log #debugkebab
    #echo $simulation $outputname 'ollaan tarkistastatuksessa LESSIN STATUKSEN tarkastelun jalkeen LES status' $les >> ${outputroot}/${simulation}/debug${outputname}.log   #debugkebab
    ############################
    #### postprosessointi STATUS
    ############################
    if [[ $les == 1 ]]; then
        if [ -f ${outputroot}/${simulation}/${outputname}.nc ] && [ -f ${outputroot}/${simulation}/${outputname}.ts.nc ] && [ -f ${outputroot}/${simulation}/${outputname}.ps.nc ]; then
            #lasttimeNC=$( ncdump -v time -f fortran ${outputroot}/${simulation}/${outputname}.nc    | tail -2 | head -1 | cut -f 1 --delimiter=/ |tr -dc '[:alnum:].' | cut -f 1 --delimiter=. ) 
	    lasttimeNC=$timmax
            lasttimePS=$( ncdump -v time -f fortran ${outputroot}/${simulation}/${outputname}.ps.nc | tail -2 | head -1 | cut -f 1 --delimiter=/ | tr -dc '[:alnum:].' | cut -f 1 --delimiter=. )
            lasttimeTS=$( ncdump -v time -f fortran ${outputroot}/${simulation}/${outputname}.ts.nc | tail -2 | head -1 | cut -f 1 --delimiter=/ | tr -dc '[:alnum:].' | cut -f 1 --delimiter=. )
	    PPtime=$( python -c "print min( $lasttimeNC, $lasttimePS, $lasttimeTS )" )

    	    if [[ $PPtime -ge $((timmax-1)) ]]; then
	        postpros=1
	    else
	        postpros=2
	    fi
	    #for f in ${outputroot}/${simulation}/${outputname}*.0*0*.nc; do # jos postprosessoituja filuja on, mutta prossufiluja on edelleen olemassa niin postprosessoidaan uusiks ja lessi on tosiaan valmis; muutoin ollaan valmiita
            #    [ -e "$f" ] && postpros=2 || postpros=1
            #    break
            #done

        else
            for f in ${outputroot}/${simulation}/${outputname}*.0*0*.nc; do # jos postprosessoituja filuja ei ole,  mutta prossufiluja on, postprosessoidaan uusiks; jos mitään ei ole kaikki menee urhomatti
                if [ -e "$f" ]; then
                    postpros=2 
                else
                    postpros=0
                    les=0
                fi
                break
            done
        fi
    else
        postpros=0
    fi
    
    #echo 'ollaan tarkistastatuksessa POSTPROS STATUKSEN tarkastelun jalkeen LES status' $les 'postpros status' $postpros  >> ${outputroot}/${simulation}/debug${outputname}.log #debugkebab
    echo ${les}${postpros}
    
}    

