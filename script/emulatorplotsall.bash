#!/bin/bash
shopt -s nullglob
plotScriptName=plotLesOutput.py

plotScriptMount=${LES}/script/${plotScriptName}
plotScriptSync=${SOFTA}/${plotScriptName}

if  [ -f $plotScriptMount ]; then
    scr="python $plotScriptMount "

elif [ -f $plotScriptSync ]; then
    scr="python $plotScriptSync "
else
    echo 'skripti ei ole olemassa'
    exit
fi

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