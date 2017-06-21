#!/bin/bash

# tarkista lapimenneet
shopt -s nullglob
mode=${mode:-last}
pattern='emul([0-9])\w+'
for i in *
do
    if [[ $i =~ $pattern ]]; then

        if [[ $mode == 'risen' ]]; then
            success=$(tail -15 ${i}/LES* | grep "Normal termination")
    
            if [[ -n $success ]];then 
                echo $i
            fi
        elif [[ $mode == 'fallen' ]]; then
    
    
            error=$(tail -15 ${i}/LES* | grep "exit codes")
    
            if [[ -n $error ]]; then 
                cat  ${i}/LES* | grep --ignore-case "model time" | tail -1
            fi
    
        elif [[ $mode == 'last' ]]; then
            timmax=$( cat ${i}NAMELIST | grep timmax | cut -c11-30 | tr -d .)
            last=$( cat  ${i}/LES* | grep --ignore-case "model time" | tail -1 | cut -c40-46 )
        
            if [[ $last -lt $((timmax-1)) ]]; then
                echo ' '
                echo "kaatui ennen aikojaan" $i
                echo $last
                echo ' '
            fi
            
    
        fi
    fi
done
