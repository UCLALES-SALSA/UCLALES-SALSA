#!/bin/bash
vanha=$1
uusi=${1:0:$((${#1}-4))}_pilkku.csv
cat $vanha | tr ',' ';' | tr '.' ',' > $uusi
sed -i 's/I1/q_inv/g'    $uusi
sed -i 's/I2/tpot_inv/g' $uusi
sed -i 's/I3/q_pbl/g'    $uusi
sed -i 's/I4/tpot_pbl/g' $uusi
sed -i 's/I5/pblh/g'     $uusi
sed -i 's/I6/num_pbl/g'  $uusi