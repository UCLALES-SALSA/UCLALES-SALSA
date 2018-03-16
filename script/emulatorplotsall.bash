#!/bin/bash
shopt -s nullglob
plotScriptName=plotLesOutput.py

plotScript=${LES}/script/${plotScriptName}

if  [ -f $plotScript ]; then
    scr="python $plotScript "


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