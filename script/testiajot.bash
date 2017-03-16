#!/bin/bash
# 
# input variables
# $1 = postfix
# compile true or false


# Exit on error
set -e

root="$(dirname $PWD)"
salsa=${root}/src/src_salsa
les=${root}/src/src_LES
bin=${root}/bin
script=${root}/script


if [ -z $1 ]; then
  echo "You didn't give any postfix name"
  exit 1
else
  postfix=$1
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

odotus=${odotus:-false}



function submitting {
	
	testi=$1
	nproc=$2
    namelistPF=$3
    historyrun=$4
	
	if [ $namelistPF == 'normal' ]; then
            namelistPF=''
        fi
 	
        
	echo ' '
	echo ' '
	echo 'Nyt suoritetaan funktiota submitting '

	nimi=${testi}
	dir=${bin}/${nimi}
	
	cp ${dir}/sound_in ${bin}/sound_in
	cp ${dir}/NAMELIST${namelistPF} ${bin}/NAMELIST
	
	
	nimi=${nimi}_${nproc}${namelistPF}_${postfix}

	echo 'nimi ' $nimi
	
	if [ $historyrun != 'initial' ]; then         
            cp ${dir}/$historyrun ${bin}/$historyrun
            
	fi
	
	
	## submit
	
    ${script}/submit_uclales-salsa.bash $nimi $nproc
	
	
	sleep 3s

}


if [ $odotus == 'true' ]; then
    while [[ ! -z $( qstat -u aholaj ) ]]
    do
        qstat -u aholaj
        sleep 3m

    done

fi


#submitting sheba
#submitting ascos

# function name  nproc namelistPF historyrun

##########
# isdac 
#########

# submitting case_isdac 64 _thrm5_all_on_fixINC_1_7800 initial
# 
submitting case_isdac 64 _thrm4                  initial
# submitting case_isdac 64 _thrm5_all_on_fixINC_1  initial
# submitting case_isdac 64 _thrm5_all_on_fixINC_4  initial

# submitting case_isdac 1 _1D_testing  initial
#  0000_0000.SPINUP7200.rst

####################
### speed tests  ###
###              ###
####################


# submitting case_speed 64    normal     initial
# 
# submitting case_speed 100   normal     initial

# submitting case_speed 400   normal     initial
# 
# submitting case_speed 100   _double     initial
# submitting case_speed 144   _double     initial




echo 'Simulaatioajojen tulostus: '
qstat -u aholaj
