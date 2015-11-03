#!/bin/bash

input_data=qmc_test.inp
enerb=ene.rb
val="-3.480828706"

cd tests

#ln -sf ../src/exact_H hamgen_H
../src//hamgen_H 1 -1 1
../src/lattgene 2  4 1
../src/dla_alg
prog="../src/dla"

echo "runtype = 0 " >> $input_data
echo "outfile = dat " >> $input_data

#$HOME/local/bin/mpirun $prog $input_data
$prog $input_data

echo "while line=gets()" >> $enerb
echo "line.chomp!" >> $enerb
echo "line1=line.split()" >> $enerb
echo "print line1[3].to_f if /ene/=~line" >> $enerb
echo "end" >> $enerb

mv dat.000 dat &>/dev/null


sed -e "s/dat.000/dat/g" dat 1> aaa  2>/dev/null

ruby $enerb aaa > a1
echo $val > a2


diff -w aaa dat.org > /dev/null 2>&1
diff -w a1 a2 > /dev/null 2>&1


if [ $? -eq 0 ]; then
echo ""
echo "++++++++++++++++++++++++++++"
echo "|  dsqss-1.1 Test Passed   |"
echo "++++++++++++++++++++++++++++"
echo ""
else
echo ""
echo "++++++++++++++++++++++++++++"
echo "|  dsqss-1.1 Test Failed   |"
echo "++++++++++++++++++++++++++++"
echo ""
fi

rm -rf aaa dat a1 a2 $enerb $input_data >/dev/null
#rm hamgen_H
rm hamiltonian.xml
rm lattice.xml
rm algorithm.xml

