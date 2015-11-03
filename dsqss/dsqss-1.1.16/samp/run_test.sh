#!/bin/bash
mpipath=/usr/bin 
mpipath=/usr/bin 
#mpipath=$HOME/local/bin #nec
binpath=../src #../bin
inpfile=qmc.inp
outfile=dla.log

############################
set_parameter_m()
{
    echo "-----"
    echo "M ... the number of bosons on each site"
    echo "      (M = 1,2,3 .... for  S= 1/2, 1, 3/2 ...)"
    echo "J ... the coupling constat"
    echo "F ... the base magnetic field in the pair Hamiltonian"
    echo "      ( = H/z if the field H is shared equally by all pairs,
                    where z = 2 for a ring geometry.)"
    echo "-----"
    echo "set parameter for hamgen_H. M J F ?"
    echo "ex) 1 1.0 0.4"
    read MM JJ FF
}

set_parameter_b()
{
    echo "-----"
    echo "D ... the dimension."
    echo "L ... the liner size of the lattice"
    echo "      ( L,L0,L1, ... must be even number."
    echo "B ... the base inverse tempereture."
    echo "-----"
    echo "set parameter for lattgene. D L B ?"
    echo "ex) 1 4 0.1"
    read DD LL BB
    
}

set_parameter_m0()
{
    echo "-----"
    echo "M ... the number of bosons on each site"
    echo "      (M = 1,2,3 .... for  S= 1/2, 1, 3/2 ...)"
    echo "J ... the coupling constat"
    echo "F ... the base magnetic field in the pair Hamiltonian"
    echo "      ( = H/z if the field H is shared equally by all pairs,
                    where z = 2 for a ring geometry.)"
    echo "-----"
    echo "set parameter for hamgen_H. M J F ?"
    echo "ex) 1 -0.5 0.0"
    read MM JJ FF
}

set_parameter_b0()
{
    echo "-----"
    echo "D ... the dimension."
    echo "L ... the liner size of the lattice"
    echo "      ( L,L0,L1, ... must be even number."
    echo "B ... the base inverse tempereture."
    echo "-----"
    echo "set parameter for lattgene. D L B ?"
    echo "ex) 1 2 1.0"
    read DD LL BB
    
}
#NREP=$1
#NREP1=`expr $NREP-1`

#if [ $# -ne 1 ];then
#    echo "Command line: $0" 1>&2
#    echo "Usage: $./run_test.sh [n] " 1>&2
#    echo "       n : The parallel numbers or Replicas" 1>&2
#    exit 1
#fi

############################
#####  input parameter #####
############################
while true
do 
    echo "dou you calculate using the replica method? [y/n]"
    read ANS
    
    if [ "${ANS}" = 'y' -o "${ANS}" = 'yes' ]; then

        echo " parallel (replica) number? [>1]"
        echo " ex) 8 "
        read NREP
        NREP1=`expr $NREP-1`

	while true
	do 
	    echo " choose the magnetic field or inverse tempereture for replica. [m/t]"
	    read ANS
	    if [ "${ANS}" = 'm' ]; then
                RUNTYPE=1
		set_parameter_m;
		
		echo "the magnetic interval?"
		echo "ex) 0.02"
		read DF 
		
		set_parameter_b;
		break;
	    elif [ "${ANS}" = 't' ]; then
                RUNTYPE=2
		set_parameter_m;
		set_parameter_b;
		
		echo "the tempereture interval?"
		echo "ex) 0.12"
		read DB 
		break;
	    else
		continue 1
	    fi
	done
	break;
    elif [ "${ANS}" = 'n' -o "${ANS}" = 'no' ]; then
        echo " parallel (replica) number? [>0]"
        echo " ex) 1 "
        read NREP
        NREP1=`expr $NREP-1`
	set_parameter_m0;
	set_parameter_b0;
        RUNTYPE=0
	break;
    else
	continue 1
    fi
done

####===================== RUNTYPE = 1 ===========================###
# -- set parameter -- #
irep=0

if [ $RUNTYPE -eq 0 ];then
    $binpath/hamgen_H $MM $JJ $FF
    $binpath/lattgene $DD $LL $BB
    $binpath/dla_alg

    rm $inpfile &>/dev/null
    echo "runtype = 0" >$inpfile

    echo "nmcse  = 1000">>$inpfile
    echo "nmcsd  = 1000">>$inpfile
    echo "nmcs   = 10000">>$inpfile
    echo "nset   = 10">>$inpfile
    echo "seed   = 31415">>$inpfile
    echo "nvermax   = 10000">>$inpfile
    echo "nsegmax   = 10000">>$inpfile

    echo "nrep=$NREP">>$inpfile

    echo "algfile=algorithm.xml">>$inpfile
    echo "latfile=lattice.xml">>$inpfile
    echo "outfile=sample.log">>$inpfile

    $mpipath/mpirun -np $NREP $binpath/dla $inpfile &>dla.log
#    sed "6i P F        =             $FF \\" $oufile.* >

fi

if [ $RUNTYPE -eq 1 ];then
    
    while [ $irep -lt $NREP ];
    do 
#+++++++++++++++
	
	XF=`echo "scale=5; $FF+$irep*$DF "  | bc`
#	FF=`echo "scale=5; $XF/0.5 "  | bc`
	
	$binpath/hamgen_H $MM $JJ $XF
	$binpath/lattgene $DD $LL $BB
	$binpath/dla_alg
	cp algorithm.xml $irep.xml
	
#++++++ exact calc +++++++
#$binpath/exact_H $M $NSITE $J $FF $B > exact_${M}_${NSITE}_${J}_${XF}_${B}.out
#++++++ exact calc +++++++
	
	irep=`expr $irep + 1`
    done
    
    rm $inpfile &>/dev/null
    echo "runtype = 1" >$inpfile

    echo "nmcse  = 10">>$inpfile
    echo "nmcsd  = 100">>$inpfile
    echo "nmcs   = 500">>$inpfile
    echo "nset   = 100">>$inpfile

    echo "nrep=$NREP">>$inpfile
    echo "vf=$FF">>$inpfile
    echo "df=$DF">>$inpfile

    echo "outfile=qmc.log">>$inpfile

 # return 
    $mpipath/mpirun -np $NREP $binpath/dla $inpfile &>dla.log
    rm aaa &>/dev/null
    
fi
####===================== RUNTYPE = 1 ===========================###

####===================== RUNTYPE = 2 ===========================###
if [ $RUNTYPE -eq 2 ];then
    $binpath/hamgen_H $MM $JJ $FF
    $binpath/lattgene $DD $LL $BB
    $binpath/dla_alg

    while [ $irep -lt $NREP ];
    do
	cp algorithm.xml $irep.xml 
	irep=`expr $irep + 1`
    done

    rm $inpfile &>/dev/null
    echo "runtype = 2" >$inpfile

    echo "nmcse  = 10">>$inpfile
    echo "nmcsd  = 100">>$inpfile
    echo "nmcs   = 500">>$inpfile
    echo "nset   = 100">>$inpfile

    echo "nrep=$NREP">>$inpfile
    echo "vb=$BB">>$inpfile
    echo "db=$DB">>$inpfile

    echo "outfile=qmc.log">>$inpfile

 # return
    $mpipath/mpirun -np $NREP $binpath/dla $inpfile &>dla.log
    rm aaa &>/dev/null
     
fi
#----- R2X-----
