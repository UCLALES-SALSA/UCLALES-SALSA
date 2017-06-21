#!/bin/bash
if [[ -d ${SCRIPT} ]]; then
   scriptref=${SCRIPT}
else
   scriptref=.
fi
source ${scriptref}/subroutines_variables.bash

emulatorname=${emulatorname:-case_emulator_DESIGN_v1.4.0_LES_nudgetus.0.1_LVL4}
emulatoroutputroot=${outputroot}/${emulatorname}
A=${A:-1}
B=${B:-99}
shopt -s extglob
echo "paakansio " $emulatorname A $A B $B
lista=""
for k in $(seq -f"%03g" $A $B )
do
    i=$(printf %02d ${k##+(0)})
    #echo i $i k $k
    if [ ! -f ${emulatoroutputroot}/emul${k}/emul${k}.nc ] || [ ! -f ${emulatoroutputroot}/emul${k}/emul${k}.ts.nc ] || [ ! -f ${emulatoroutputroot}/emul${k}/emul${k}.ps.nc ]; then
        echo ' '
        echo ' '
	lista=${lista}" "${k}
	echo "alikansio" emul${k} "file" emul$i
        if [ ! -f ${emulatoroutputroot}/emul${k}/emul${k}.nc ]; then
  		echo emul${k} "NC PUUTTUU"
	fi
	
	if [ ! -f ${emulatoroutputroot}/emul${k}/emul${k}.ts.nc ]; then
		echo emul${k} "TS PUUTTUU"
	fi 
 
	if [ ! -f ${emulatoroutputroot}/emul${k}/emul${k}.ps.nc ]; then
		echo emul${k} "PS PUUTTUU"
	fi
	
	echo  "emul${k} tarkistastatus: " $(tarkistastatus  ${emulatorname}/emul${k} emul${k})
        
        #ncsub=true tssub=true pssub=true postprosessoi ${emulatorname}/emul${k} emul${k} FIX${k}
        #tulosta=true odota FIX${k}
	echo  "JALKEEN tarkistastatus: " $(tarkistastatus  ${emulatorname}/emul${k} emul${k})
	#sleep 5s
    fi

    if [  -f ${emulatoroutputroot}/emul${k}/emul${k}.nc ] && [  -f ${emulatoroutputroot}/emul${k}/emul${k}.ts.nc ] && [  -f ${emulatoroutputroot}/emul${k}/emul${k}.ps.nc ]; then
	echo ' '
	echo  "KAIKKI LOYTYY emul${k} tarkistastatus: " $(tarkistastatus  ${emulatorname}/emul${k} emul${k})
        #poistaturhat ${emulatorname}/emul${k} emul${i}
	#mv ${emulatoroutputroot}/emul${k}/emul${i}.nc ${emulatoroutputroot}/emul${k}/emul${k}.nc
	#mv ${emulatoroutputroot}/emul${k}/emul${i}.ps.nc ${emulatoroutputroot}/emul${k}/emul${k}.ps.nc
	#mv ${emulatoroutputroot}/emul${k}/emul${i}.ts.nc ${emulatoroutputroot}/emul${k}/emul${k}.ts.nc
        #echo  "KAIKKI LOYTYY emul${k} tarkistastatus: " $(tarkistastatus  ${emulatorname}/emul${k} emul${k})
 
        echo $k valmis
    fi
done
echo $lista
