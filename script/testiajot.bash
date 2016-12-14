#!/bin/bash

# Exit on error
set -e

if [ -z $1 ]; then
  echo "You didn't give any postfix name"
  exit 1
fi

# if [ -z $2 ]; then
#   echo "You didn't give waiting parameter ( true / false )"
#   exit 1
# fi

postfix=$1


if [ -z $2 ]; then
  odotus='false'
else
  odotus=$2
fi





root="$(dirname $PWD)"
salsa=${root}/src/src_salsa
les=${root}/src/src_LES
bin=${root}/bin
script=${root}/script


cd ${root}
make seq || exit 1
make mpi || exit 1
cd ${script}
echo " "
echo "Compiled"

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
	
	if [ $nproc -gt 1 ]; then
            ${script}/ajoskripti_MPI.bash $nimi $nproc
	else
            ${script}/ajoskripti.bash $nimi
	fi
	
	sleep 3s

}


if [ $odotus == 'true' ]; then
    while [[ ! -z $( qstat -u aholaj ) ]]
    do
        qstat -u aholaj
        sleep 15m

    done

fi


#submitting sheba
#submitting ascos

# function name  nproc namelistPF historyrun

##########
# isdac 
#########

submitting case_isdac 64 _thrm5_all_on_fixINC_1_7600 initial

submitting case_isdac 64 _thrm4                  initial
submitting case_isdac 64 _thrm5_all_on_fixINC_1  initial
submitting case_isdac 64 _thrm5_all_on_fixINC_4  initial


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
