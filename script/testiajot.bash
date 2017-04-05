#!/bin/bash

# Exit on error
set -e

# import subroutines & variables 
source ./subroutines_variables.bash


copyOUT=${copyOUT:-false}

submit=${submit:-false}
postpros=${postpros:-false}
odotus=${odotus:-false}
poista=${poista:-false}

compiler=${compiler:-intel.zero}
vers=${vers:-Jaakko}

RUNTYPE='"INITIAL"'

LVL=${LVL:-5}
isoT=${isoT:-28800.}
pikkuT=${pikkuT:-7200.}
ice=${ice:-1.0} # ice #/kg

dim1=${dim1:-false}
dim2=${dim2:-false}
dim3=${dim3:-false}

if [[ -n $hfilebase ]]; then
    RUNTYPE='"HISTORY"'
    
fi

jaakkoNL=${jaakkoNL:-false}

if [[ $jaakkoNL == 'false' ]]; then
    jaakkoNL='!'
else
    jaakkoNL=''
fi

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

#   copy hfilin file

	if [[ -n $hfilebase ]]; then
            echo $hfilebase
            cp ${bin}/0*_0*.${hfilebase}  ${outputroot}/${nimi}/
            modifyoutputHistory='false'
    else
            modifyoutputHistory='true'
	fi

}


    
function submitting {
	
	inputsubfolder=$1$
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
	
    input=${bin}/${inputsubfolder} mode=${mode} modifyoutput='true' modifyoutputHistory=${modifyoutputHistory} COPY=${copyOUT} clean=${copyOUT} ownjobname=$ownjobnameSUB ${script}/submit_uclales-salsa.bash $nimi $nproc $jobflag
	
	sleep 3s
    qstat -u $USER
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
    if [[ -n $ownjobnameMAIN ]]; then
        odotusLES=LES_${ownjobnameMAIN}
    else
        odotusLES=${nimi:$((${#nimi}-${length})):${length}}
    fi
    
    if [[ $copyOUT == 'false' ]] && [[ ${submit} == 'true' ]]; then
        copy $inputsubfolder $nimi $mode
        
        dir=${outputroot}/${nimi} nxp=${nxp} nyp=${nyp} nzp=${nzp} deltax=${deltax} deltay=${deltay} deltaz=${deltaz} nxpart=${nxpart} dzmax=${dzmax} dzrat=${dzrat} dtlong=${dtlong} distim=${distim} timmax=${timmax} Tspinup=${Tspinup} minispinup01=${minispinup01} minispinup02=${minispinup02} minispinupCase01=${minispinupCase01} minispinupCase02=${minispinupCase02} runtype=${runtype} level=${level} CCN=${CCN} prndtl=${prndtl} filprf=${filprf} hfilin=${hfilin} ssam_intvl=${ssam_intvl} savg_intvl=${savg_intvl} mcflg=${mcflg} frqhis=${frqhis} istpfl=${istpfl} lbinanl=${lbinanl} frqanl=${frqanl} corflg=${corflg} ipsflg=${ipsflg} itsflg=${itsflg} strtim=${strtim} sed_aero=${sed_aero} sed_cloud=${sed_cloud} sed_precp=${sed_precp} sed_ice=${sed_ice} sed_snow=${sed_snow} iradtyp=${iradtyp} case_name=${case_name} div=${div} sfc_albedo=${sfc_albedo} radsounding=${radsounding} cntlat=${cntlat} strtim=${strtim} isfctyp=${isfctyp} sst=${sst} dthcon=${dthcon} drtcon=${drtcon} ubmin=${ubmin} zrough=${zrough} th00=${th00} umean=${umean} vmean=${vmean} nlcoag=${nlcoag} nlcgcc=${nlcgcc} nlcgpp=${nlcgpp} nlcgaa=${nlcgaa} nlcgii=${nlcgii} nlcgss=${nlcgss} nlcgpc=${nlcgpc} nlcgca=${nlcgca} nlcgpa=${nlcgpa} nlcgia=${nlcgia} nlcgic=${nlcgic} nlcgip=${nlcgip} nlcgsa=${nlcgsa} nlcgsc=${nlcgsc} nlcgsi=${nlcgsi} nlcgsp=${nlcgsp} nlcnd=${nlcnd} nlcndgas=${nlcndgas} nlcndh2oae=${nlcndh2oae} nlcndh2ocl=${nlcndh2ocl} nlcndh2oic=${nlcndh2oic} nlauto=${nlauto} nlautosnow=${nlautosnow} nlactiv=${nlactiv} nlactbase=${nlactbase} nlactintst=${nlactintst} nlichom=${nlichom} nlichet=${nlichet} nlicimmers=${nlicimmers} nlicmelt=${nlicmelt} nlicbasic=${nlicbasic} nlfixinc=${nlfixinc} fixINC=${fixINC} rhlim=${rhlim} isdtyp=${isdtyp0} nspec=${nspec1} listspec=${listspec} volDistA=${volDistA} volDistB=${volDistB} nf2a=${nf2a} sigmag=${sigmag} dpg=${dpg} n=${n} notJJA=${notJJA} ${script}/generate_namelist_ISDAC.bash
        
    fi
    
    
	if [[ $submit == 'true' ]]; then
        
        echo ' '
        echo -n 'Käynnistetään Simulaatio ' $nimi ' '; date '+%T %d-%m-%Y'
        echo ' '
        ownjobnameSUB=$ownjobnameMAIN submitting $inputsubfolder $nimi $nproc $mode
        qstat -u $USER | tail -1
    fi

    if [[ $odotus == 'true' ]]; then
        echo 'odotetaan'
        odota $odotusLES
        echo -n 'Simulaatio on valmis: ' $nimi' '; date '+%T %d-%m-%Y'
    fi
    
    
    length=$(( ${#nimi} < 7 ? ${#nimi} : 7))
    if [[ -n $ownjobnameMAIN ]]; then
        odotusPP=${ownjobnameMAIN}
    else
        odotusPP=${nimi:$((${#nimi}-${length})):${length}}
    fi
    
    
    if [[ $postpros == 'true' ]]; then
        echo 'postprosessoidaan'
        
        postprosessoi $nimi $nimi $ownjobnameMAIN
        
        qstat -u $USER | tail -1
        
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




if [[ $LVL == 5 ]]; then
    icePF=_INC${ice}
fi

if [[ -n $testinumero ]]; then
    testinumero=_testi${testinumero}
fi    

# 1D
if [[ ${dim1} == 'true' ]]; then
    echo 1D
    fixINC=$ice timmax=$isoT Tspinup=$pikkuT runtype=$RUNTYPE hfilin="'"${hfilebase}"'" level=$LVL ownjobnameMAIN=${LVL}_1D${testinumero} jaakkoNL=${jaakkoNL}     nyp=5 nxp=5 postfix=LVL${LVL}_1D${icePF}${testinumero}  main case_isdac 1 seq.${vers}.${compiler}
fi

# 2D
if [[ ${dim2} == 'true' ]]; then
    echo 2D
    fixINC=$ice timmax=$isoT Tspinup=$pikkuT runtype=$RUNTYPE hfilin="'"${hfilebase}"'" level=$LVL ownjobnameMAIN=${LVL}_2D${testinumero} jaakkoNL=${jaakkoNL}           nxp=5 postfix=LVL${LVL}_2D${icePF}${testinumero}  main case_isdac 8 mpi.${vers}.${compiler}
fi

# 3D
if [[ ${dim3} == 'true' ]]; then
    echo 3D
    fixINC=$ice timmax=$isoT Tspinup=$pikkuT runtype=$RUNTYPE hfilin="'"${hfilebase}"'" level=$LVL ownjobnameMAIN=${LVL}_3D${testinumero} jaakkoNL=${jaakkoNL}                 postfix=LVL${LVL}_3D${icePF}${testinumero} main case_isdac 64 mpi.${vers}.${compiler}
fi
##################
qstat -u $USER