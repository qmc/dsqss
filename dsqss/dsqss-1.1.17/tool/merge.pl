#!/usr/bin/perl

$BASE = shift;

@FLIST = `ls $BASE.???`;
$NSMP=1;
$NLINE=0;
$NCOM=0;
$NFILE=0;
foreach $f (@FLIST) {
  printf(">  $f");
  open(FIN, "<$f");
  while(<FIN>) {
    if ( /^\s*(P|S)\s+(\w+)\s+=\s+(\d+)/ ) {
      $thiskey = $2;
      $thisval = $3;
      if ( $thiskey eq "NSMP" ) { $NSMP = $thisval; }
      if ( $parameter{$thiskey} eq "" ) {
        $parameter{$thiskey} = $thisval;
        $NLINE++;
        if ( $line[$NLINE] eq "" ) { $line[$NLINE] = $_; }
      }
      if (( $thiskey ne "SEED"  )&&($thiskey ne "NCYC" )) {
#      if ( $thiskey ne "SEED" ){
        if ( $parameter{$thiskey} ne $thisval ) {
          printf("Error. Parameter ( $thiskey ) mismatch.\n");
          exit(0);
        }
      }
    } elsif ( /^\s*R\s+(\w+)\s+=\s+(\S+)\s+(\S+)/ ) {
      $key = $1;
      if ( $QID{$key} eq "" ) {
        $NQ++;
        $qkey[$NQ] = $key;
        $QID{$key} = $NQ;
        $ndat[NQ] = 0;
      }
      $qid = $QID{$key};
      $ndat[$qid]++;
      $i = $ndat[$qid];
      $a[$qid][$i] = $2;
      $d[$qid][$i] = $3;
    } else {
      $NCOM++;
      if ( $com[$NCOM] eq "" ) { $com[$NCOM] = $_; }
    }
  }
}

open(FOUT,">$BASE");

for ($i=1; $i<=$NLINE; $i++) { print FOUT $line[$i]; }

printf FOUT ("S NPROC    = %16d\n", $ndat[1]);

for ($i=1; $i<=$NCOM; $i++) { print FOUT $com[$i]; }

for ($qid=1; $qid<=$NQ; $qid++) {
  $key=$qkey[$qid];
  $NFILE=$ndat[$qid];
  if ($NFILE == 1) {
    printf FOUT ("R %6s = %16.6f %16.6f\n", $key, $a[$qid][1], $d[$qid][1]);
  } else {
    $ave=0.0;
    $sgm=0.0;
    $err=0.0;
    for ($i=1; $i<=$NFILE; $i++) {
      $x = $a[$qid][$i];
      $e = $d[$qid][$i];
      $ave += $x;
      $sgm += $x * $x;
      $err += $NSMP * ( $x * $x + ($NSMP - 1) * $e * $e );
    }
    $ave /= $NFILE;
    $sgm /= $NFILE;
    if ( $NFILE == 1 ) {
	$sgm = 0.0;
    } else {
	$sgm = sqrt( ($sgm-$ave*$ave) / ($NFILE-1) );
    }
    $NTOT = $NSMP * $NFILE;
    $err /= $NTOT;
    if ( $NTOT == 1 ) {
	$err = 0.0;
    } else {
	$err = sqrt( ($err-$ave*$ave) / ($NTOT-1) );
    }
#   printf FOUT ("R %6s = %16.6e %16.6e %16.6e\n", $key, $ave, $sgm, $err);
     printf FOUT ("R %6s = %16.10e %16.10e\n", $key, $ave, $err);
#    printr ($key,$ave,$err);
  }
}

#system("rm -f $BASE.???");

