#!/bin/bash

# Exit on error
set -e




postfix=$1


salsa=../src/src_salsa
les=../src/src_LES
bin=../bin
scripts=../scripts


cd ${bin}
make seq || exit 1
make mpi || exit 1
cd ${scripts}
echo " "
echo "Compiled"

function submitting {
	
	testi=$1
	dim=$2
	nproc=$3
	NLsalsaPF=$4  
	NLpostfix=$5

	if [ $NLsalsaPF == normal ]; then
		$NLsalsaPF=''
	fi


	if [ $NLpostfix == normal ]; then
		$NLpostfix=''
	fi
	
 	pre=campaign_
        
	echo ' '
	echo ' '
	echo 'Nyt suoritetaan funktiota submitting '

	nimi=${pre}${testi}
	dir=${bin}/${nimi}
	echo 'nimi ' $nimi
	
	cp ${dir}/namelist.salsa${NLsalsaPF} ${bin}/namelist.salsa

	cp ${dir}/sound_in ${bin}/sound_in
	
	
	
	nimi=${nimi}_${dim}${NLsalsaPF}${NLpostfix}_${postfix}
	if [ $dim == 3 ]; then
	  cp ${dir}/NAMELIST_3D${NLpostfix} ${bin}/NAMELIST
	  ${scripts}/ajoskripti_MPI.bash $nimi $nproc
	  
	elif [ $dim == 2 ]; then
	  cp ${dir}/NAMELIST_2D ${bin}/NAMELIST
	  ${scripts}/ajoskripti_MPI.bash $nimi $nproc
	  
	else
	  cp ${dir}/NAMELIST_1D ${bin}/NAMELIST
	  cp ${dir}/0000_0000.SPINUP7200.rst ${bin}/0000_0000.SPINUP7200.rst
	  ${scripts}/ajoskripti.bash $nimi
	fi
	
	sleep 3s

}


#submitting sheba
#submitting ascos
# submitting isdac 3 64 _thrm4
# submitting isdac 3 64 _thrm5

# submitting isdac 3 64 _all_on_fixINC_1 _thrm5
# submitting isdac 3 64 _all_on_fixINC_1 _thrm4
# 
# submitting isdac 3 64 _all_on_fixINC_4 _thrm5
# submitting isdac 3 64 _all_on_fixINC_4 _thrm4

# submitting isdac 1 1 _all_on_fixINC_4 

# submitting isdac 2 8 _all_on_fixINC_4 _thrm5
submitting speed 3 64 normal normal

submitting speed 3 64 normal normal

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

echo 'Simulaatioajojen tulostus: '
qstat -u aholaj
