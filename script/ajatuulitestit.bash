#!/bin/bash

DESIGN='/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/designs/corr_design_15d.csv'


# A=28 B=30 ./submit_emulator_runs.bash 3
# ./kopio_ibrixille.bash 0.2K 25:00:00
# 
# sleep 5m
# 
# python emulator_inputs.py | tee emulator_output_50
# 
# A=28 B=30 ./submit_emulator_runs.bash 3
# ./kopio_ibrixille.bash 50m 25:00:00
# 


current=$(date +%s)
target=$(date -d '02/01/2017 22:00' +%s)
sleep_seconds=$(( $target - $current ))
sleep $sleep_seconds
date "+%a %x %T"
./kopio_ibrixille.bash zero

sleep 30m
date "+%a %x %T"
A=28 B=30 ./submit_emulator_runs.bash 3



current=$(date +%s)
target=$(date -d '02/02/2017 04:00' +%s)
sleep_seconds=$(( $target - $current ))
sleep $sleep_seconds
date "+%a %x %T"
./kopio_ibrixille.bash ascos

sleep 30m
date "+%a %x %T"

A=28 B=30 ./submit_emulator_runs.bash 3


current=$(date +%s)
target=$(date -d '02/02/2017 08:00' +%s)
sleep_seconds=$(( $target - $current ))
sleep $sleep_seconds
date "+%a %x %T"
./kopio_ibrixille.bash dycoms