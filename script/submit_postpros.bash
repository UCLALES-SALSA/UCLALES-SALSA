#!/bin/bash

# Submit postprosessing to queue, Jaakko Ahola, FMI, 09/2016
#
#
#


# Exit on error
set -e


if [ -z $1 ]; then
  echo "You didn't give any name of netcdf input file"
  exit 1
fi

input=$1

dir=$(dirname ${input})

if [ $dir == '.' ]; then
  dir=$PWD
fi

length=$(( ${#1} < 5 ? ${#1} : 5))



# Create run script
# 
cd ${dir}

rm -rf post_* *pros.sh combine2.py

cp /home/users/aholaj/UCLALES-SALSA/bin/combine2.py ${dir}/

cat > runpostpros.sh <<FINAL
#!/bin/sh
#PBS -N post_${1:$((${#1}-${length})):${length}}
#PBS -l mppwidth=1
#PBS -l mppnppn=1
#PBS -l mppdepth=20
#PBS -l walltime=08:00:00
#PBS -j oe
#PBS -M jaakko.ahola@fmi.fi
#PBS -m ae

source /etc/profile
module load Python

cd ${dir}

aprun -n1 -N1 -d20 ./postpros.sh | tee ${PBS_JOBNAME:-post_interactive}.${PBS_JOBID:-help}



exit
FINAL

cat > postpros.sh <<EOF
#!/bin/bash

set -e

cd ${dir}

python combine2.py $input

EOF

cd ${dir}
# Make initial submit
chmod +x runpostpros.sh combine2.py postpros.sh
qsub runpostpros.sh


#
exit

