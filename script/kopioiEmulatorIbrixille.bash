#!/bin/bash

### 
### usage: give input folder as $1 and optional output folder as $2
### copying is not executed until there's no EMULATOR submitted psub

# Exit on error
set -e

# import subroutines & variables 
if [[ -d ${SCRIPT} ]]; then
   scriptref=${SCRIPT}
else
   scriptref=.
fi

echo "scriptref" $scriptref
source ${scriptref}/subroutines_variables.bash

if [ -z $1 ]; then
  echo "You didn't give any name for output subfolder of lustre"
  exit 1
fi
simulation=$1

outputname=$2

# if outputname doesn't exist set it to be same as simulation
if [[ -z $outputname ]]; then
    outputname=${simulation}
    echo "You didn't give any name for name for destination subfolder of ibrix, use default" $outputname
     
fi

if [[ -d ${outputroot}/${simulation}/ ]]; then
    echo "Lustre folder ${outputroot}/${simulation} exists. We're good to go."
else
    echo "Lustre folder ${outputroot}/${simulation} does not exist. EXIT."
    exit 1
fi

aika=${aika:-10m}
while [[ ! -z $( qstat -u $USER | grep EMULATOR ) ]]
do
	echo ' '
	date "+%a %x %T"
	qstat -u $USER | grep EMULATOR
	sleep $aika
done

echo ' '

kopioibrixille $simulation $outputname

setfacl -Rm g:tut-kuo:rx /ibrix/arch/ClimRes/aholaj/
setfacl -Rm g:climres:rx /ibrix/arch/ClimRes/aholaj/



# printf "Training simulations are ready in ${ibrixrootfolder}${outputname} \n\n -Jaakko- \n\n PS This was an automated message" | mail -s "Training simulations ready $outputname" Jaakko.Ahola@fmi.fi,Muzaffer.Ege.Alper@fmi.fi -r jaakko.ahola@fmi.fi # Ege
printf "Training simulations are ready in ${ibrixrootfolder}${outputname} \n\n -Jaakko- \n\n PS This was an automated message" | mail -s "Training simulations ready $outputname" Jaakko.Ahola@fmi.fi,Muzaffer.Ege.Alper@fmi.fi,Tomi.Raatikainen@fmi.fi -r jaakko.ahola@fmi.fi