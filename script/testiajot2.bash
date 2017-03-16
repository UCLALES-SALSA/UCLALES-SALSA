#!/bin/bash

# Exit on error
set -e


#### user parameters

username=aholaj
email=jaakko.ahola@fmi.fi
subfolder=UCLALES-SALSA
outputroot=/lustre/tmp/${username}/${subfolder}/${subsubfolder}
root=/home/users/${username}/${subfolder}

jobflag=${jobflag:-PBS}
scriptname=${scriptname:-combine.py}
bin=${root}/bin
script=${root}/script

copyOUT=${copyOUT:-false}

submit=${submit:-false}
postpros=${postpros:-false}
odotus=${odotus:-false}
poista=${poista:-false}

function odota {
    
    nimi=$1
    aika=$2
    aika=${aika:-30s}
    
    echo ' '
    echo 'Nyt odotetaan' $nimi $aika
    while [[ ! -z $( qstat -u $username | grep $nimi ) ]]
    do
        qstat -u $username
#         sleep 15s
        sleep $aika
    done

}

function copy {
	
	inputsubfolder=$1
	nimi=$2
	mode=$3



#     if [ $clean == true ]; then
    if [ -d ${outputroot}/${nimi} ] ; then
        rm -rf ${outputroot}/${nimi}/*
    else
        mkdir -p ${outputroot}/${nimi}/
    fi
#     fi    
    datadir=${outputroot}/${nimi}/datafiles
    mkdir -p ${datadir} 

    cp ${bin}/les.${mode} ${outputroot}/${nimi}/
    cp ${bin}/datafiles/* ${datadir}/
    
    cp ${bin}/${inputsubfolder}/sound_in ${outputroot}/${nimi}/
#    cp ${input}/NAMELIST ${rundir}/

}


    
function submitting {
	
	inputsubfolder=$1
	nimi=$2
	nproc=$3
	mode=$4
	
    if [ -z $mode ]; then
        $mode=mpi
    fi    

	echo ' '
	echo 'Nyt suoritetaan funktiota submitting '
	echo 'nimi ' $nimi
	
	
	## submit
	
    input=${bin}/${inputsubfolder} mode=${mode} modifyoutput='true' COPY=${copyOUT} clean=${copyOUT} ${script}/submit_uclales-salsa.bash $nimi $nproc $jobflag
	
	sleep 3s
    qstat -u $username

}


function postprosessoi {
    
    echo ' '
    nimi=$1
    echo 'Nyt postprosessoidaan nc'
    scriptfolder=${script} ${script}/submit_postpros.bash ${outputroot}/${nimi}/${nimi}    $scriptname $jobflag $nimi
    
    echo ' '
    echo 'Nyt postprosessoidaan ps'
    scriptfolder=${script} ${script}/submit_postpros.bash ${outputroot}/${nimi}/${nimi}.ps $scriptname $jobflag $nimi
    
    echo ' '
    echo 'Nyt postprosessoidaan ts'
    scriptfolder=${script} ${script}/submit_postpros.bash ${outputroot}/${nimi}/${nimi}.ts $scriptname $jobflag $nimi
    echo ' '
    
    echo 'Kaikki submittoitu postprosessointiin'

}

function poistaturhat {
    
    echo ' '
    echo 'Poistetaan lustrelta turhat tiedostot'
    nimi=$1
    if [ -f ${outputroot}/${nimi}/${nimi}.nc ] && [ -f ${outputroot}/${nimi}/${nimi}.ts.nc ] && [ -f ${outputroot}/${nimi}/${nimi}.ps.nc ]; then
        echo "kaikki kolme postprosessoitua tiedostoa ovat olemassa"
        ls ${outputroot}/${nimi}/${nimi}.nc
        ls ${outputroot}/${nimi}/${nimi}.ts.nc
        ls ${outputroot}/${nimi}/${nimi}.ps.nc
        echo 'poistetaan'
        rm -rf ${outputroot}/${nimi}/datafiles
        rm -rf ${outputroot}/${nimi}/*.sh
        rm -rf ${outputroot}/${nimi}/*.py
        rm -rf ${outputroot}/${nimi}/*.rst
        rm -rf ${outputroot}/${nimi}/${nimi}.ts.0*0*.nc
        rm -rf ${outputroot}/${nimi}/${nimi}.ps.0*0*.nc
        rm -rf ${outputroot}/${nimi}/${nimi}.0*0*.nc
        rm -rf ${outputroot}/${nimi}/0*_0*.${nimi}.*
    fi

    

    echo ' '

}

function main {

	inputsubfolder=$1
	nproc=$2
	mode=$3
	
	if [[ -n $postfix ]]; then
        postfix=_${postfix}
    fi
	
    nn=$(( ${#mode}- 4 ))
    if [[ $nn -gt 0 ]]; then
        nimi=${inputsubfolder}${postfix}_${mode:$((${#mode}-$nn)):$nn}
    else
        nimi=${inputsubfolder}${postfix}
    fi 
    
    length=$(( ${#nimi} < 6 ? ${#nimi} : 6))
    odotusLES=${nimi:$((${#nimi}-${length})):${length}}
    
    if [[ $copyOUT == 'false' ]] && [[ ${submit} == 'true' ]]; then
        copy $inputsubfolder $nimi $mode
        
        dir=${outputroot}/${nimi} nxp=${nxp} nyp=${nyp} nzp=${nzp} deltax=${deltax} deltay=${deltay} deltaz=${deltaz} nxpart=${nxpart} dzmax=${dzmax} dzrat=${dzrat} dtlong=${dtlong} distim=${distim} timmax=${timmax} Tspinup=${Tspinup} minispinup01=${minispinup01} minispinup02=${minispinup02} minispinupCase01=${minispinupCase01} minispinupCase02=${minispinupCase02} runtype=${runtype} level=${level} CCN=${CCN} prndtl=${prndtl} filprf=${filprf} hfilin=${hfilin} ssam_intvl=${ssam_intvl} savg_intvl=${savg_intvl} mcflg=${mcflg} frqhis=${frqhis} istpfl=${istpfl} lbinanl=${lbinanl} frqanl=${frqanl} corflg=${corflg} ipsflg=${ipsflg} itsflg=${itsflg} strtim=${strtim} sed_aero=${sed_aero} sed_cloud=${sed_cloud} sed_precp=${sed_precp} sed_ice=${sed_ice} sed_snow=${sed_snow} iradtyp=${iradtyp} case_name=${case_name} div=${div} sfc_albedo=${sfc_albedo} radsounding=${radsounding} cntlat=${cntlat} strtim=${strtim} isfctyp=${isfctyp} sst=${sst} dthcon=${dthcon} drtcon=${drtcon} ubmin=${ubmin} zrough=${zrough} th00=${th00} umean=${umean} vmean=${vmean} nlcoag=${nlcoag} nlcgcc=${nlcgcc} nlcgpp=${nlcgpp} nlcgaa=${nlcgaa} nlcgii=${nlcgii} nlcgss=${nlcgss} nlcgpc=${nlcgpc} nlcgca=${nlcgca} nlcgpa=${nlcgpa} nlcgia=${nlcgia} nlcgic=${nlcgic} nlcgip=${nlcgip} nlcgsa=${nlcgsa} nlcgsc=${nlcgsc} nlcgsi=${nlcgsi} nlcgsp=${nlcgsp} nlcnd=${nlcnd} nlcndgas=${nlcndgas} nlcndh2oae=${nlcndh2oae} nlcndh2ocl=${nlcndh2ocl} nlcndh2oic=${nlcndh2oic} nlauto=${nlauto} nlautosnow=${nlautosnow} nlactiv=${nlactiv} nlactbase=${nlactbase} nlactintst=${nlactintst} nlichom=${nlichom} nlichet=${nlichet} nlicimmers=${nlicimmers} nlicmelt=${nlicmelt} nlicbasic=${nlicbasic} nlfixinc=${nlfixinc} fixINC=${fixINC} rhlim=${rhlim} isdtyp=${isdtyp0} nspec=${nspec1} listspec=${listspec} volDistA=${volDistA} volDistB=${volDistB} nf2a=${nf2a} sigmag=${sigmag} dpg=${dpg} n=${n} notJJA=${notJJA} ${script}/generate_namelist_ISDAC.bash
        
    fi
    
    
	if [[ $submit == 'true' ]]; then
        
        echo ' '
        echo -n 'Käynnistetään Simulaatio ' $nimi ' '; date '+%T %d-%m-%Y'
        echo ' '
        submitting $inputsubfolder $nimi $nproc $mode
        qstat -u $username | grep $odotusLES
    fi

    if [[ $odotus == 'true' ]]; then
        echo 'odotetaan'
        odota $odotusLES
        echo -n 'Simulaatio on valmis: ' $nimi' '; date '+%T %d-%m-%Y'
    fi
    
    length=$(( ${#nimi} < 7 ? ${#nimi} : 7))
    odotusPP=${nimi:$((${#nimi}-${length})):${length}}
    
    if [[ $postpros == 'true' ]]; then
        echo 'postprosessoidaan'
        
        postprosessoi $nimi
        
        qstat -u $username | grep $odotusPP
        
    fi
    
    if [[ $odotus == 'true' ]]; then
        odota nc_$odotusPP
        echo -n 'Postprosessointi on valmis: ' nc_$odotusPP' '
        date '+%T %d-%m-%Y'

        odota ps_$odotusPP 5s
        echo -n 'Postprosessointi on valmis: ' ps_$odotusPP' '; date '+%T %d-%m-%Y'
    
        odota ts_$odotusPP 5s
        echo -n 'Postprosessointi on valmis: ' ts_$odotusPP' '; date '+%T %d-%m-%Y'
    
        echo ' '
        echo 'Kaikki on postprosessoitu'
    fi
        
    if [[ $poisto == 'true' ]]; then 
        poistaturhat $nimi
        echo -n 'Turhat poistettu'' '; date '+%T %d-%m-%Y'
        echo 'Valmis' $nimi 
        echo 'Valmis' $nimi >  ${outputroot}/${nimi}/valmis${i}
        date '+%T %d-%m-%Y' >>  ${outputroot}/${nimi}/valmis${i}
    fi
}


# main case_testaus_lvl3 100 mpi.cray.2017.01.13  
# main case_testaus_lvl3 100 mpi.mergattu        
# main case_testaus_lvl3 100 mpi.v1.0.4           
# main case_testaus_lvl3 100 mpi.hotfix.1.0.5     
# 
# main case_testaus_lvl4 100 mpi.cray.2017.01.13 
# main case_testaus_lvl4 100 mpi.mergattu         
# main case_testaus_lvl4 100 mpi.v1.0.4           
# main case_testaus_lvl4 100 mpi.hotfix.1.0.5    


# 2D
postfix=LVL4_2D nxp=5 level=4          main case_isdac 8 mpi.mergattu
postfix=LVL4_2D nxp=5 level=4 notJJA=! main case_isdac 8 mpi.v1.0.4
postfix=LVL4_2D nxp=5 level=4 notJJA=! main case_isdac 8 mpi.cray.2017.01.13

# 3D
postfix=LVL4_3D level=4 notJJA=! main case_isdac 64 mpi.mergattu
postfix=LVL4_3D level=4 notJJA=! main case_isdac 64 mpi.v1.0.4


qstat -u $username