#!/bin/bash

#J=1.0  # ferromagnetic
J=-1.0  # antiferromagnetic
#G=2.0  # easy axis
G=-4.0  # easy plane
Hmin=0.0 # starting field
DH=1.0 # field step
Nstep=11 # number of steps
M=2
D=1
L=4
B=4.0


Jxy=`echo "sqrt(($J)*($J))" | bc` ;

HTEMP="hamiltonian.xml.templete"
JTEMP="jobscript.sh.templete"
ITEMP="qmc.inp.templete"
LTEMP="lattice-${L}.xml.templete"
JOBSCRIPT="jobscript.sh"
HAMDAT="ham.dat"
Lfile="lattice.xml"

if [ ! -f "$LTEMP" ] ; then
    echo; echo "Error. $LTEMP does not exist."; echo
    exit 
fi

cat $LTEMP | sed "s/###BETA/${B}/g" > $Lfile

cp $JTEMP $JOBSCRIPT

for (( i=0; i<$Nstep; i++ )) ; do
    Istr=`printf "%03d" $i` ;
    H=`echo "$Hmin + $DH * $i" | bc` ;
    intH=`echo "1000 * $H" | bc` ;
    Hstr=`printf "%04.0f" $intH` ;
    echo $i "  " $H ;
    h0000=`echo " ($J) + ($G) - ($H)" | bc` ;
    h2222=`echo " ($J) + ($G) + ($H)" | bc` ;
    h0011=`echo "0.5*($G) - 0.5*($H)" | bc` ;
    h1100=$h0011 ;
    h2211=`echo "0.5*($G) + 0.5*($H)" | bc` ;
    h1122=$h2211 ;
    h0022=`echo "(-1.0) * ($J) + ($G)" | bc` ;
    h2200=$h0022 ;
    h1001=`echo "$Jxy" | bc` ;
    h0110=$h1001 ;
    h2112=$h1001 ;
    h1221=$h1001 ;
    h1210=$h1001 ;
    h1012=$h1001 ;
    h2101=$h1001 ;
    h0121=$h1001 ;
    Hfile="ham_H${Hstr}.xml" ;
    Afile="alg_H${Hstr}.xml" ;
    Ofile="sample_${Istr}.log" ;
    Jfile="job_H${Hstr}.xml" ;
    Ifile="inp_H${Hstr}.dat" ;
    echo $Hfile ;
    h0000str=$( printf "%20.16f" $h0000 ) ;
    h2222str=$( printf "%20.16f" $h2222 ) ;
    h1100str=$( printf "%20.16f" $h1100 ) ;
    h0011str=$( printf "%20.16f" $h0011 ) ;
    h1122str=$( printf "%20.16f" $h1122 ) ;
    h2211str=$( printf "%20.16f" $h2211 ) ;
    h0022str=$( printf "%20.16f" $h0022 ) ;
    h2200str=$( printf "%20.16f" $h2200 ) ;
    h1001str=$( printf "%20.16f" $h1001 ) ;
    h0110str=$( printf "%20.16f" $h0110 ) ;
    h2112str=$( printf "%20.16f" $h2112 ) ;
    h1221str=$( printf "%20.16f" $h1221 ) ;
    h1210str=$( printf "%20.16f" $h1210 ) ;
    h1012str=$( printf "%20.16f" $h1012 ) ;
    h2101str=$( printf "%20.16f" $h2101 ) ;
    h0121str=$( printf "%20.16f" $h0121 ) ;
    cat $HTEMP | sed "s/###int0000/${h0000str}/g" \
               | sed "s/###int2222/${h2222str}/g" \
               | sed "s/###int1100/${h1100str}/g" \
               | sed "s/###int0011/${h0011str}/g" \
               | sed "s/###int1122/${h1122str}/g" \
               | sed "s/###int2211/${h2211str}/g" \
               | sed "s/###int0022/${h0022str}/g" \
               | sed "s/###int2200/${h2200str}/g" \
               | sed "s/###int1001/${h1001str}/g" \
               | sed "s/###int0110/${h0110str}/g" \
               | sed "s/###int2112/${h2112str}/g" \
               | sed "s/###int1221/${h1221str}/g" \
               | sed "s/###int1210/${h1210str}/g" \
               | sed "s/###int1012/${h1012str}/g" \
               | sed "s/###int2101/${h2101str}/g" \
               | sed "s/###int0121/${h0121str}/g" > $Hfile ;
    dla_alg $Hfile $Afile ;
#    echo "dla $Ifile" >> $JOBSCRIPT ;
    cat $ITEMP | sed "s/###ALGFILE/${Afile}/g" \
	       | sed "s/###OUTFILE/${Ofile}/g" > $Ifile ;
done
    echo "ls -1 inp_H*.dat | xargs -n 1 -P $Nstep dla" >> $JOBSCRIPT ;

echo "np_triv_p=$Nstep" > $HAMDAT
echo "triv_p=f" >> $HAMDAT
echo "triv_p_step=$DH" >> $HAMDAT
echo "M=$M" >> $HAMDAT
echo "J=$J" >> $HAMDAT
echo "H=$Hmin" >> $HAMDAT
echo "G=$G" >> $HAMDAT
echo "D=$D" >> $HAMDAT
echo "L=$L" >> $HAMDAT
echo "B=$B" >> $HAMDAT
