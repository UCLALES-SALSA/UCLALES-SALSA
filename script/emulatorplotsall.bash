#!/bin/bash
shopt -s nullglob
scr="python ${LES}/script/plotLesOutput.py "
puuttuvat="puuttuvat keissit: "
pattern='emul([0-9])\w+'
for i in *
do 
    if [[ $i =~ $pattern ]]; then
        if  [ -f ${i}/${i}.ts.nc ] && [ -f ${i}/${i}.ps.nc ]; then #[ -f ${i}/${i}.nc ] &&
           scr="${scr} ${i}/${i}"
        else
            puuttuvat="${puuttuvat} ${i}"
        fi
    fi
done
$scr
echo $puuttuvat
