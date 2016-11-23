#!/bin/bash

# Exit on error
set -e

postfix=$1

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
	echo 'nimi ' $nimi
	
	cp ${dir}/sound_in ${bin}/sound_in
	cp ${dir}/NAMELIST${namelistPF} ${bin}/NAMELIST
	
	
	nimi=${nimi}_${nproc}${namelistPF}_${postfix}
	
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


#submitting sheba
#submitting ascos



##########
# isdac 
#########
# submitting isdac 3 64 _thrm4
# submitting isdac 3 64 _thrm5

# submitting isdac 3 64 _all_on_fixINC_1 _thrm5
# submitting isdac 3 64 _all_on_fixINC_1 _thrm4
# 
# submitting isdac 3 64 _all_on_fixINC_4 _thrm5
# submitting isdac 3 64 _all_on_fixINC_4 _thrm4

# submitting isdac 1 1 _all_on_fixINC_4 

# submitting isdac 2 8 _all_on_fixINC_4 _thrm5

##################################################
# submitting isdac 3 64 _all_on_fixINC_4 _thrm5
# 
# submitting isdac 3 64 _all_on_fixINC_4 _thrm4

##################################
# submitting isdac 2 8 _init_ice_all_off

# submitting isdac 1 1 _init_ice_all_off

#submitting isdac 1 1 _init_ice_cond_activ_on
# submitting isdac 2 8 _init_iceliq_cond_on
# submitting isdac 2 8 _init_iceliq_immers_on
# submitting isdac 1 1 _init_iceliq_basic_on

# submitting isdac 1 1 _basic_on
#submitting isdac 1 1 _fixedinc_on
# submitting isdac 1 1 _immers_on


# submitting isdac 1 1 _init_iceliq_cond_on


#  0000_0000.SPINUP7200.rst

##
### speed tests


# function name  nproc namelistPF historyrun
submitting speed 64    normal     initial

submitting speed 100   normal     initial

submitting speed 400   normal     initial

echo 'Simulaatioajojen tulostus: '
qstat -u aholaj
