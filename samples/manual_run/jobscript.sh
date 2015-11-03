#!/bin/sh
#PBS -l nodes=1:ppn=12
#PBS -q small
#PBS -N rtype_0
#PBS -j oe
export PATH=/home2/u00300/Git/worm/worm/bin:$PATH
cd $PBS_O_WORKDIR

dla inp_H0000.dat
dla inp_H1000.dat
dla inp_H2000.dat
dla inp_H3000.dat
dla inp_H4000.dat
dla inp_H5000.dat
dla inp_H6000.dat
dla inp_H7000.dat
dla inp_H8000.dat
dla inp_H9000.dat
dla inp_H10000.dat
