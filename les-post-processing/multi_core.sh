#!/bin/sh

python postp_UCLALES-SALSA.py -i wk82modseed -o wk82modseed_postp_binned_Na -v S_Naba,S_Nabb  &

python postp_UCLALES-SALSA.py -i wk82modseed -o wk82modseed_postp_binned_Nc -v S_Ncba,S_Ncbb  &

python postp_UCLALES-SALSA.py -i wk82modseed -o wk82modseed_postp_binned_NpNi -v S_Npba,S_Niba  &

python postp_UCLALES-SALSA.py -i wk82modseed -o wk82modseed_postp_binned_Rwc -v S_Rwcba,S_Rwcbb  &

python postp_UCLALES-SALSA.py -i wk82modseed -o wk82modseed_postp_binned_RwpRwi -v S_Rwpba,S_Rwiba  &

wait
