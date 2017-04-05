#!/bin/bash

scr="python /home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/script/thermo4pathsemulator.py "
puuttuvat="puuttuvat keissit: "
for i in $(ls -d emul*)
do 
if [ -f ${i}/${i}.nc ] && [ -f ${i}/${i}.ts.nc ] && [ -f ${i}/${i}.ps.nc ]; then
    scr="${scr} ${i}/${i}"
else
    puuttuvat="${puuttuvat} ${i}"
fi
done
$scr
echo $puuttuvat
