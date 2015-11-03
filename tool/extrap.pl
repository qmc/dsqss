#!/usr/bin/perl

$M = shift;

$pltfile = "__plotfile__.gp";

@filelist = glob("*.log");

foreach $filename (@filelist){


open (PLTFILE,"> $pltfile");


print PLTFILE "set print \"$filename.param\"\n";


print PLTFILE "set fit errorvariables\n";

#Paraboric fitting
print PLTFILE "f(x)=a*x**2+p0\n";
print PLTFILE "fit [0.01:0.04] f(x) \"$filename\" u 1:2:3 via a,p0\n";
#Linear fitting
print PLTFILE "g(x)=b*x+l0\n";
print PLTFILE "fit g(x) \"$filename\" u 1:2:3 via b,l0\n";


print PLTFILE "print $M,p0,p0_err,l0,l0_err,a,a_err,b,b_err\n";


print PLTFILE "set out\"$filename.eps\"\n";


print PLTFILE "set term post eps enh color\n";


print PLTFILE "p f(x)                             lc 1 title \"parabopric fit\",\\\n";
print PLTFILE "  g(x)                             lc 2 title \"linear fit\",\\\n";
print PLTFILE "  \"$filename\"       u 1:2:3 w e lc 3 title \"QMC\",\\\n";
print PLTFILE "  \"$filename.param\" u 1:2:3 w e lc 1 notitle,\\\n";
print PLTFILE "  \"$filename.param\" u 1:4:5 w e lc 2 notitle\n";
print PLTFILE "exit\n";


close (PLTFILE);

$result=`gnuplot $pltfile`;
}
