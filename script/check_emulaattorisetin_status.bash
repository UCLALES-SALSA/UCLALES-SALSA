#!/bin/bash

# Exit on error
set -e
shopt -s extglob
# import subroutines & variables 
if [[ -d ${SCRIPT} ]]; then
   scriptref=${SCRIPT}
else
   scriptref=.
fi
source ${scriptref}/subroutines_variables.bash
#######################

# input variables
#
# $1 = absolute path of emulator
#

emulatorname=$(basename $1)
folderROOT=$(dirname $1)



#######################
A=${A:-1}
B=${B:-90}


etunolla=3 # huomhuom
if [[ -z $list ]]
then
    array=($(seq -f"%0${etunolla}g" $A $B  ))
else
    u=0
    for kk in ${list[@]}
    do
        array[u]=$(printf %0${etunolla}d ${kk##+(0)})
        u=$((u+1))
    done    
fi
    
for i in ${array[@]}
do

    echo ' '
	echo ' '
    	status=$( tarkistastatus ${emulatorname}/emul${i} emul${i}  $folderROOT)

    	echo "statuksen tarkistus" emul${i} $status
    	LS=${status:0:1}
    	PPS=${status:1:2}
    	
    	if [[ $status == '11' ]]
    	then
        	echo emul${i} "on VALMIS"
        	statusValmiit=$((statusValmiit+1))
    	elif [ $LS -eq 1 ] && [ $PPS -ne 1 ]; then
            statusVainLESValmis=$((statusVainLESValmis+1))
    	elif [ $LS -ne 1 ] && [ $PPS -ne 1 ]; then
            statusKesken=$((statusKesken+1))
    	fi
    
    
done

echo 'valmiit' $statusValmiit
echo 'vain les valmis' $statusVainLESValmis
echo 'kesken' $statusKesken
echo 'tarkistussumma kaikki' $((statusValmiit+statusVainLESValmis+statusKesken))
echo ' '
echo ' '
########################################