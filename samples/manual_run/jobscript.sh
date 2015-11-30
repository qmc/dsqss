#!/bin/sh
#PBS -l nodes=1:ppn=12
#PBS -q small
#PBS -N rtype_0
#PBS -j oe
#PBS -V
export PATH=$WORM_HOME:$PATH
cd $PBS_O_WORKDIR

ls -1 inp_H*.dat | xargs -n 1 -P 11 dla
