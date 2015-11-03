#!/bin/bash


merpl=$WORM_HOME/tool/merge.pl

XPARAM=${1} # !=triv_p
#M, H, L, B



. ./ham.dat


P_s=0
P_e=$(( $np_triv_p - 1 ))

if [ ! -d "Arc" ]; then
mkdir Arc
fi

dP=$triv_p_step
J0=$J
H0=$H
L0=$L
B0=$B
G0=$G


for ii in `seq $P_s $P_e`; do

    PPN=`printf "%03d\n" $ii`
    if [ $P_s -eq $P_e ];then
	file="sample.log"
    else
	file="sample_${PPN}.log"
    fi

    echo "merge ${file}.XXX"


    $merpl $file

    if [ $np_triv_p -gt "1" ]; then

	echo "ii=$ii,dP= $dP, J0= $J0"

	if [ $triv_p = 'j' ]; then
	    J=`echo "scale=8; $J0+$dP*$ii"  | bc`
	elif [ $triv_p = 'f' ]; then
	    H=`echo "scale=8; $H0+$dP*$ii"  | bc`
	elif [ $triv_p = 'l' ]; then
	    L=`echo "scale=8; $L0+$dP*$ii"  | bc`
	elif [ $triv_p = 'b' ]; then
	    B=`echo "scale=8; $B0+$dP*$ii"  | bc`
	elif [ $triv_p = 'g' ]; then
	    G=`echo "scale=8; $G0+$dP*$ii"  | bc`
	fi
    fi

    if [ -z $G ]; then
	pfile=Arc/HM${M}D${D}_L${L}B${B}J${J}H${H}.dat
    else
	pfile=Arc/HM${M}D${D}_L${L}B${B}J${J}H${H}G${G}.dat
    fi
    echo "OUTPUT: $pfile"

    rm -f $pfile &>/dev/null
    echo "P M = $M" >  $pfile
    echo "P J = $J" >> $pfile
    echo "P H = $H" >> $pfile
    echo "P G = $G" >> $pfile
    cat $file       >> $pfile


done

cd Arc
tar zcvf sample.tar.gz HM*.dat
cd ..
