#!/usr/bin/perl

# this program

$PROG = "plot.pl";
$VERSION = "20060609";

# list of options with arguments

%OPTLIST = ( 
  -f => "",
);

# list of flags

%FLAGLIST = ( 
  -v => "",
  '-?' => "",
);

# gnuplot command

#$GNUPLOT="/usr/local/bin/wgnuplot"; 
$GNUPLOT="/usr/bin/gnuplot"; 

#-----------------------------------------------------------------------

#$X = "[a] [bb] <ccc> [dddd 1 2 3]";
#printf("$X\n");
#ExtractKeys( $X , \@KEYS );
#foreach $k ( @KEYS ) { printf("$k\n"); };

#$X="NBT";
#$X="NBT NSMP";
#$X=" ( BETA<3 && LX < 2 ) ";
#$X="\"NBT= \" NBT \", NSMP=\" NSMP";
#$X=" sgm(3) ";
#$Y=normalize($X);
#printf("$X --> $Y\n");
#exit(0);

#$X = "abc [s_k 1 2 3] / (<aaa>)**2 \$arg\$";
#$X = shift;
#$Y = Evaluate($X);
#printf("$X --> $Y\n");
#exit(0);

#$X = "[abc 1 \$second\$ 3]";
#ExtractRunningArguments( $X , @akey , @qkey , @place );
#exit(0);

########################################################################
#################    M A I N   F U N C T I O N    ######################
########################################################################

$WORKDIR="_PLOT_"; 

# コマンドラインの処理

ProcessCommandLine(@ARGV); 
if ( defined $FLAG{"-\?"} ) { Usage(); }
if ( defined $FLAG{"-v"}  ) { $verbose="1"; }

# コンフィグレーションファイルの名前

$CFILE=$OARG{"-f"}; 
if ( "$CFILE" eq "" ) { $CFILE = "plot.config"; }
if ( ! -f $CFILE ) { Abort("Configuration file \"$CFILE\" doesn't exist."); }

# コンフィグレーションファイルの読み取り

ReadConfigFile(); 

# アーカイブファイルがあるか？

if( -f $ARCHIVE ) { 
    $ARCHIVE =~ /([^\/]+)\.(tgz|tar.gz)$/;
    $TOCFILE = "$1.toc";
} else {
    Abort("The archive file \"$ARCHIVE\" does not exist."); 
};

# 出力ディレクトリを準備

if ( "$OUTDIR" eq "" ) { $OUTDIR = "Plot"; }
if ( -d "$OUTDIR" ) {
    system("rm -f $OUTDIR/*.dat");
} else { 
    system("mkdir $OUTDIR"); 
};
$GPFILE="$OUTDIR/plot.gp"; 

# 必要なデータファイルをアーカイブファイルから抜き出して作業ディレクトリに展開

CreateTableOfContents();

# 作業ディレクトリの準備

if ( -d $WORKDIR ) { system("rm -r $WORKDIR"); }
system("mkdir $WORKDIR"); 
ExtractDatFiles();

# データファイルの読み取り

$FILES=`cd $WORKDIR; ls -v`;
@FLIST=split(/\s+/,$FILES);
$PlotDatFiles="";
foreach $f (@FLIST) {
#    print "datafile=$f\n";

    chomp( $f ); 
    ReadOneFile( $f );
    WritePlottingData( $f );
};

# あとしまつ

#if ( -d $WORKDIR ) { system("rm -r $WORKDIR"); }

# ＧＰファイルの作成と表示

@PDFILE = split(/\s+/,$PlotDatFiles);
$NF = $#PDFILE + 1;
if ( $NF == 0 ) {
    printf("\n");
    printf("$PROG> Error.\n");
    printf("$PROG> No file contains the data to be plotted.\n");
    printf("$PROG> This happens, for example, when you try to\n");
    printf("$PROG> plot a quantity that doesn't exist in any data file.\n");
    printf("\n");
    exit(0);
}

CreateGPFile(); 
system("cd $OUTDIR; $GNUPLOT plot.gp; ps2pdf plot.ps");
system("cp $OUTDIR/plot.pdf .");

exit(0);

########################################################################
####################   S U B R O U T I N S  ############################
########################################################################

sub Usage() {
    print "\n";
    print "usage:\n";
    print "\n";
    print "  \$ plot.pl [-f CONFIG_FILE] [-v]\n";
    print "\n";
    print "      -v ... verbose mode\n";
    print "\n";
    exit(0);
};

#-----------------------------------------------------------------------

sub Abort() {
    $message = shift;
    print "\n";
    print "${PROG}> Error.\n";
    print "    ${message}\n";
    print "\n";
    exit(0);
};

#-----------------------------------------------------------------------

sub ProcessCommandLine() {

    foreach $opt (@OPTLIST) { $OARG{$opt} = ""; }

    @A = ( @_ , "-" );

    while ( @A ) {
        $a = shift @A;
        last if ("$a" eq "-");
#       print "\/$a\/ ";
        if ( $a =~ /^[-](\w|[?])/ ) {
            if ( defined $OPTLIST{$a} ) {
                $b = shift @A;
#               print "--> option with argument \/$b\/\n";
                $OARG{"$a"} = $b;
            } elsif ( defined $FLAGLIST{$a} ) {
#               print "--> flag\n";
                $FLAG{"$a"} = "";
            } else {
#               print "--> invalid option\n";
            }
        } else {
#           print "--> file\n";
        }
    }
};

#-----------------------------------------------------------------------

sub ConvertKey { # convert a group of a key and its arguments
    # into a single string ready to be evaluated
    # Example:  [abc 1 22] --> $Q{"abc:1:22"}

my $K0 = shift;
    $K0 =~ /^(.)(.*)(.)$/;
    my $left = $1;
    my $middle = $2;
    my $right = $3;
    my @W = ();
    @W = split(/\s+/,$middle);
    $K = $W[0];
    for $i (1 .. $#W) {
        $K .= ":".$W[$i];
    }
    if ( "$left" eq "\[" ) { $K = "\$Q\{\"".$K."\"\}"; }
    if ( "$left" eq "\<" ) { $K = "\$D\{\"".$K."\"\}"; }
    if ( "$left" eq "\$" ) { $K = "\$A\{\"".$K."\"\}"; }
    return $K;
    
};

#-----------------------------------------------------------------------

sub Replace { # replace keys in a string by their evaluatable form
    # Example: [abc]/[def] --> $Q{"abc"}/$Q{"def"}
    my $left, $key, $right;
    $right = shift;
    my $X = "";
    while ( "$right" ne "" ) {
        if ( $right =~ /(.*?)(\[.+?\]|\<.+?\>|\$.+?\$)(.*)/ ) {
            $left = $1;
            $key = $2;
            $right = $3;
        } else {
            $left = $right;
            $key = "";
            $right = "";
        }
        $X .= $left;
	if ( "$key" ne "" ) { $X .= ConvertKey( $key ); }
    }
    return $X;
};

#-----------------------------------------------------------------------

sub Evaluate { # convert a string into a value
    my $input = shift;
    my $str = Replace( $input );
    my $V = eval($str);
    return $V;
};

#-----------------------------------------------------------------------

sub normalize { # convert a key string into a normal form

    @ODEL = ( "\"" , "\$" , "\[" , "\<" );
    @CDEL = ( "\"" , "\$" , "\]" , "\>" );
    $nrepl = 0;
    @srepl = ();
    my $X = shift;

    # [...] や "..." を #0#, #1#, ... に置き換えていく
    for $i (0 .. $#ODEL) {
        $op = $ODEL[$i];
        $cl = $CDEL[$i];
        if ("$op" eq "\<") { # 大小記号と区別するために条件を厳しく
            $reg = "([$op][^&|]+?[$cl])";
        } else {
            $reg = "([$op].+?[$cl])";
        }
        while ( $X =~ /$reg/ ) {
            $srepl[$nrepl] = $1;
            $X =~ s/[$op].+?[$cl]/\#$nrepl\#/; 
            $nrepl++;
        }
    }

    # 残った裸のキーをブラケットで包む
    while ( $X =~ /(^|\W)(\p{IsAlpha}\w*?)($|\W)/ ) { 
        $srepl[$nrepl] = "\[$2\]"; 
        $X =~ s/(^|\W)(\p{IsAlpha}\w*?)($|\W)/$1\#$nrepl\#$3/;
        $nrepl++;
    }

    # #0#, #1#, ... をもとに戻す
    while ( $X =~ /\#[0-9]+\#/ ) {
        for $i ( 0 .. $#srepl ) {
            $X =~ s/\#$i\#/$srepl[$i]/; # 置き換えた文字列をもとにもどす
        }
    }

    return $X;

};

#-----------------------------------------------------------------------

sub GetKey { 

    $str = shift;
    if    ( $str =~ /\[(.+)\]/ ) { $bare_str = $1; }
    elsif ( $str =~ /\<(.+)\>/ ) { $bare_str = $1; }
    else { 
      printf("GetKey> Error. Does not seem a single bracketed string.\n"); 
      exit(0);
    }
    my @W = split(/\s+/,"$bare_str");
    $qkey = $W[0];
    return $qkey;

};

#-----------------------------------------------------------------------

sub ExtractKeys { # make a list of quantity keys contained in $str

    my ( $str , $ref_qkey ) = @_;
    my %CHK = ();
    my @QK = ();
    my $left, $middle, $right;
    $right = $str;
    while ( "$right" ne "" ) {
        if ( $right =~ /(.*?)(\[.+?\]|\<.+?\>)(.*)/ ) {
            $left = $1;
            $middle = $2;
            $right = $3;
        } else {
            $left = $right;
            $middle = "";
            $right = "";
        }
        $qkey = GetKey( $middle );
        if ( ! defined( $CHK{$qkey} ) ) {
            $CHK{$qkey} = "";
            @QK = ( @QK , $qkey );
        }
    }
    @$ref_qkey  = @QK;

};

#-----------------------------------------------------------------------

sub GetRA { 
# About $akey, $qkey, and $place
# Suppose an expression "[abc 0 1 $x$ 3]" appears in XPLOT field of 
# plot.config. Then, $akey = "x", $qkey = "abc", $place = "2".

    $str = shift;
#printf("GetRA> input= $str\n");

    if    ( $str =~ /\[(.+)\]/ ) { $bare_str = $1; }
    elsif ( $str =~ /\<(.+)\>/ ) { $bare_str = $1; }
    else { 
      printf("GetRA> Error. Does not seem a single bracketed string.\n"); 
      exit(0);
    }

#printf("GetRA> bare_string= $bare_str\n");
    my @W = split(/\s+/,"$bare_str");
    $akey = "";
    $qkey = "";
    $place = "";
    for $i (1 .. $#W) {
        if ( $W[$i] =~ /^\$(.+)\$$/ ) {
            $akey = $1;
            $qkey = $W[0];
            $place = $i - 1;
            break;
        }
    }
#printf(" \$akey=  $akey\n");
#printf(" \$qkey=  $qkey\n");
#printf(" \$place= $place\n");
    return ( $akey , $qkey , $place );

};

#-----------------------------------------------------------------------

sub ExtractRunningArguments {

    my ( $plots , $ref_akey, $ref_qkey , $ref_place ) = @_;

#printf("\nExtractRunningArguments> $plots\n");

    my %CHK = ();
    my @AK = ();
    my @QK = ();
    my @PL = ();

#printf("\$\#AK= $#AK\n");

    my $left, $middle, $right;
    $right = $plots;
    while ( "$right" ne "" ) {
        if ( $right =~ /(.*?)(\[.+?\]|\<.+?\>)(.*)/ ) {
            $left = $1;
            $middle = $2;
            $right = $3;
        } else {
            $left = $right;
            $middle = "";
            $right = "";
        }
#printf("$middle\n");
        if ( "$middle" ne "" ) {
          ( $akey , $qkey , $place ) = GetRA( $middle );
        } else {
            $akey = "";
            $qkey = "";
            $place = "";
        }
        if ( "$akey" ne "" ) {
            if ( ! defined( $CHK{$akey} ) ) {
                $CHK{$akey} = "";
                @AK = ( @AK , $akey );
                @QK = ( @QK , $qkey );
                @PL = ( @PL , $place );
            }
        }
    }

#printf("\$\#AK= $#AK\n");

    @$ref_akey  = @AK;
    @$ref_qkey  = @QK;
    @$ref_place = @PL;
};

#-----------------------------------------------------------------------

sub ReadConfigFile() { # read the configuration file
    $NBEFORE=0;
    $NAFTER=0;
    open(FIN,"<$CFILE");
    while(<FIN>) {
        if ( /^\s*\#/ ) { next; }
        if ( ! /\w+/ ) { next; }
        if ( /^\s*(\S+)\s*::\s*(\S|\S.*\S)\s*$/ ) {
            if ( $1 eq "ZDATA" || $1 eq "ARCHIVE") { $ARCHIVE = $2; };
            if ( $1 eq "BASE" || $1 eq "OUTDIR") { $OUTDIR = $2; };
            if ( $1 eq "XPLOT" ) { $XPLOT_org = $2; };
            if ( $1 eq "YPLOT" ) { $YPLOT_org = $2; };
            if ( $1 eq "DPLOT" ) { $DPLOT_org = $2; };
            if ( $1 eq "PLOTID" ) { $PLOTID_org = $2; };
            if ( $1 eq "LEGEND" ) { $LEGEND_org = $2; };
            if ( $1 eq "CONDITION" ) { $CONDITION_org = $2; };
            if ( $1 eq "BEFOREPLOT" ) { 
                $BEFOREPLOT[$NBEFORE] = $2; 
                $NBEFORE++;
            };
            if ( $1 eq "AFTERPLOT" ) { 
                $AFTERPLOT[$NAFTER] = $2; 
                $NAFTER++;
            };
        }
    };
    close(FIN);

    $XPLOT  = &normalize($XPLOT_org);
    $YPLOT  = &normalize($YPLOT_org);
    $DPLOT  = &normalize($DPLOT_org);
    $LEGEND = &normalize($LEGEND_org);
    $PLOTID = &normalize($PLOTID_org);

# 以下の小細工は大小記号が誤差を表す記号と同じことによる誤解を防ぐため．
# CONDITION 指定のなかでは <...> を誤差をあらわす記号として使うことはできない．
    $CONDITION = &normalize($CONDITION_org);
    $CONDITION =~ s/\</\@SMALLER\@/g;
    $CONDITION =~ s/\>/\@GREATER\@/g;
#printf("CONDITION= |$CONDITION|\n"); 
    $PLOTS = "$XPLOT $YPLOT $DPLOT";

#   $AKEY[.]  = the argument specifier key that appears in $PLOT
#   $QKEY[.]  = the quantity specifier key that the argument is associated with
#   $PLACE[.] = order in which the argument appears
    ExtractRunningArguments( $PLOTS , \@AKEY , \@QKEY , \@PLACE );

#   @RequiredKey = [the list of keys of the records needed for plotting]
    ExtractKeys( $PLOTID  , \@RequiredKey  );

};

#-----------------------------------------------------------------------

sub CreateGPFile() {
    open(FGP,">$GPFILE");
    printf FGP ("set title \"$CONDITION_org\"\n");
    printf FGP ("set size 1.0, 0.9\n");
    printf FGP ("set xlabel \"$XPLOT\"\n");
    printf FGP ("set ylabel \"$YPLOT\"\n");
    printf FGP ("set style data errorbars\n");
    printf FGP ("set key outside\n");
    printf FGP ("set key samplen 0\n");
    for ($i=0; $i<$NBEFORE; $i++) {
        printf FGP ("%s\n",$BEFOREPLOT[$i]);
    }
    $command="plot";
    @OL=split(" ",$PlotDatFiles);
    foreach $f (@OL) {
        printf FGP ("$command \"$f\" title \"$LEG{$f}\" with errorbars\n");
        if ($command eq "plot") { $command = "replot"; }
    };
    printf FGP ("pause -1\n");
    printf FGP ("set size 1.10,0.50\n");
    printf FGP ("set term postscript portrait color\n");
    printf FGP ("set output \"plot.ps\"\n");
    printf FGP ("replot\n");
    close(FGP);
};

#-----------------------------------------------------------------------

sub CreateTableOfContents() {
    
# 作業ディレクトリの準備
    if ( -d $WORKDIR ) { 
        system("rm -f $WORKDIR/*");
    } else {
        system("mkdir $WORKDIR"); 
    };
    
    $time_arc=(stat($ARCHIVE))[9];
    $time_toc=(stat($TOCFILE))[9];
    return if( $time_arc < $time_toc );
    
    print "\n";
    print "$PROG> The archive file is new or has been updated.\n";
    print "$PROG> Creating TOC file \"$TOCFILE\" ... \n";
    
    if ( -f $TOCFILE ) { system("rm -f $TOCFILE"); }
    system("cd $WORKDIR; gunzip -c ../$ARCHIVE | tar xf -");
    $FILES=`cd $WORKDIR; ls`;
    @FLIST=split(/\s+/,$FILES);
#   @FLIST=split(/\s+/,`cd $WORKDIR; ls`);
    foreach $f (@FLIST) {
        chomp($f); 
        ReadOneFileForTOC($f);
        $LINE = $f." ";
        foreach $k ( keys %P ) {
            $LINE .= " ".$k."=".$P{$k};
        };
        system("echo \'$LINE\' >> $TOCFILE");
    };
    
    print "$PROG> Done.\n\n";
    system("rm -r $WORKDIR");
}

#-----------------------------------------------------------------------

sub ExtractDatFiles() {
    $FLIST="";
    open(FIN,"<$TOCFILE");
    while(<FIN>) {
        @word = split;
        $f = $word[0];
        for ($i=1; $i<=$#word; $i++) {
            $w = $word[$i];
            $w =~ /(\S+)=(\S+)/;
            $k = $1;
            $v = $2;
            $Q{$k} = $v;
        }
        $C = Replace($CONDITION);
        $C =~ s/\@SMALLER\@/</g;
        $C =~ s/\@GREATER\@/>/g;
        if (eval($C)) {
            $FLIST .= " ".$f;
        }
    }
    close(FIN);

    if ( "$FLIST" eq "" ) {
        printf("\n");
        printf("$PROG> Error.\n");
        printf("$PROG> No files satisfies the condition.\n");
        printf("$PROG> Check the CONDITION field.\n");
        printf("$PROG> This error occurs, for example, when the field\n");
        printf("$PROG> contains a key that is not defined in any\n");
        printf("$PROG> of the data files.\n");
        printf("\n");
        exit(0);
    }
    system("cd $WORKDIR; gunzip -c ../$ARCHIVE | tar xf - $FLIST");

};

#-----------------------------------------------------------------------

sub ReadOneFile() {
    $FILE = shift;
# $P{$key} = [the value of the parameter labeled as $key]
# $Q{$key} = [the value of the quantity labeled as $key]
# $D{$key} = [the estimated error in $Q{$key}]
# $NARG{$key} = [ # of arguments the quantity $key takes]
# $SizeOfArgumentList{$key}[$i] = 
#            [ # of distinct values that the "$i"-th argument 
#             of the quantity "$key" takes ]
# $Argument{$key}[$i][$j] =
#            [ The "$j"-th value of the "$i"-th argument of the 
#             quantity "$key"]
    %P=();    
    %Q=();    
    %D=();    
    %NARG=(); 
    %SizeOfArgVal=(); 
    %ArgVal=();
    my %CHK=();

    $ID="";
    $PID="";
    open(FIN,"<$WORKDIR/$FILE");
    $fout="";
    while(<FIN>) {
#print;
        if (/^P\s+(\S+)\s*=\s*(\S+)/) {
            $key = $1;
            $val = $2;
            $ID="${ID}_${key}_${val}";
            $P{$key} = $val;
            $Q{$key} = $val;
            if ( "$PLOTID" =~ /\[$key\]/ ) {
                $PID="${PID}_${key}_${val}";
            }
        }
        if (/^R\s+(\S.*\S|\S)\s*=\s*(\S+)/) {
#print; #koko
            if (/^R\s+(\S.*\S|\S)\s*=\s*(\S+)\s+(\S+)/) {
                $err = $3;
            } else {
                $err = 0.0;
            }
            $ave = $2;
            @w = split(/\s+/,$1);
            $na = $#w;
            $key = $w[0];
#printf(" --> key= %s, na= %d, val= %f (%f)\n", $key, $na, $ave, $err); #koko
            if ( ! defined( $NARG{$key} ) ) {
#printf(" key= %s\n", $key);
                $NARG{$key} = $na;
                for $i ( 0 .. $na-1 ) {
                    $SizeOfArgVal{$key}[$i] = 0;
                }
            } else {
                if ( $na+0 != $NARG{$key}+0 ) {
                    printf("ReadOneFile> Error.");
                    printf(" The number of arguments does not match");
                    printf(" the previous record with the same key.\n");
                    exit(0);
                }
            }
            $K = $key;
            for $i ( 0 .. $na-1 ) {
                $arg = $w[$i+1];
                if ( ! defined( $CHK{$key}[$i]{$arg} ) ) {
                    $CHK{$key}[$i]{$arg}="";
                    $j = $SizeOfArgVal{$key}[$i];
                    $ArgVal{$key}[$i][$j] = $arg;
# $ArgVal{$key}[$i][$j] is the $j-th value of the $i-th
# argument of the quantity labeled as $key.
                    $SizeOfArgVal{$key}[$i] = $j + 1;
                }
                $K = $K.":".$arg;
            }
#printf(" --> $K\n"); #koko
            if ( $na == 0 ) {
              $Q{$K} = $ave;
              $D{$K} = $err;
            }
#$str = "\$Q\{\"$K\"\}";
#$val = eval($str);
#printf(" --> $str = $Q{$K} = $val\n");
        }
    }
    close(FIN);

    #foreach $k ( keys(%NARG) ) {
    #   printf(" key= %-8s", $k);
    #   printf(" narg= %1d\n", $NARG{$k});
    #   for $i ( 0 .. $NARG{$k}-1 ) {
    #       printf("    ... arg No. %1d", $i);
    #       printf(" nval= %2d\n", $SizeOfArgVal{$k}[$i]);
    #       for $j ( 0 .. $SizeOfArgVal{$k}[$i]-1 ) {
    #           printf("        ... val No. %1d -> %8d\n", 
    #                  $j, $ArgVal{$k}[$i][$j]);
    #       }
    #   }
    #}
    #    %$refP = %P
    #    %$refQ = %Q
    #    %$refD = %D
    #    %$refNARG = %NARG
    #    %$refSOAV = %SizeOfArgVal
    #    %$refAV = %ArgVal
};

#-----------------------------------------------------------------------

sub ReadOneFileForTOC() {
    $FILE = shift;
    %P=();    
    %Q=();    
    %D=();    
    %NARG=(); 
    %SizeOfArgVal=(); 
    %ArgVal=();
    my %CHK=();

    $ID="";
    $PID="";
    open(FIN,"<$WORKDIR/$FILE");
    $fout="";
    while(<FIN>) {
        if (/^P\s+(\S+)\s*=\s*(\S+)/) {
            $key = $1;
            $val = $2;
            $ID="${ID}_${key}_${val}";
            $P{$key} = $val;
            $Q{$key} = $val;
            if ("$PLOTID" =~ $key) {
                $PID="${PID}_${key}_${val}";
            }
        }
    }
    close(FIN);
};

#-----------------------------------------------------------------------

sub ValidityCheck { # check if the file contains required data
    $ans = "yes";
    foreach $k (@RequiredKey) {
        if ( ! defined( $P{$k} ) ) { 
#printf("    $k ... does not exist in the file\n");
            $ans = "no"; 
        } else {
#printf("    $k ... OK\n");
        }
    };
    return $ans;
};

#-----------------------------------------------------------------------

sub WritePlottingData {

    my $datfile = shift;
    my $anydata = "no";

    my $valid = ValidityCheck();
    if ( "$valid" eq "no" ) { return; }

    my $fout = "plot${PID}.dat"; 

    open(FOUT,">>${OUTDIR}/${fout}");
    if ( $#AKEY == -1 ) {
        my $X0 = Evaluate($XPLOT);
        my $Y0 = Evaluate($YPLOT);
        my $D0 = Evaluate($DPLOT);
        if ( "$X0" ne "" && "$Y0" ne "" ) {
#if ( "$Y0" ne "0" ) {
            $anydata = "yes";
            printf FOUT " $X0 $Y0 $D0\n";
#}
        }
    } else {
        my @NA = ();
        my $NREC = 1;
        for my $i (0 .. $#AKEY) {
            my $akey  = $AKEY[$i];
            my $qkey  = $QKEY[$i];
            my $place = $PLACE[$i];
            my $na = $SizeOfArgVal{$qkey}[$place];
            @NA = ( @NA , $na );
            $NREC *= $na;
        }
        for my $R (0 .. $NREC-1) {
            my $r = $R;
            for $i (0 .. $#AKEY) {
                my $j = $r % $NA[$i];
                my $r = $r / $NA[$i];
                $akey  = $AKEY[$i];
                $qkey  = $QKEY[$i];
                $place = $PLACE[$i];
                $aval  = $ArgVal{$qkey}[$place][$j];
                $A{$akey}  = $aval;
            }
            $X0 = $XPLOT;
            $Y0 = $YPLOT;
            $D0 = $DPLOT;
            for $i (0 .. $#AKEY) {
                $akey = $AKEY[$i];
                $aval = $A{$akey};
                $X0 =~ s/\$$akey\$/$aval/g;
                $Y0 =~ s/\$$akey\$/$aval/g;
                $D0 =~ s/\$$akey\$/$aval/g;
            }

            $X0 = Evaluate($X0);
            $Y0 = Evaluate($Y0);
            $D0 = Evaluate($D0);

            if ( "$X0" ne "" && "$Y0" ne "" ) {
                $anydata = "yes";
                printf FOUT " $X0 $Y0 $D0\n";
            }
        }
    }
    close(FOUT);

    if ( "$anydata" eq "yes" ) {
      my $LEG0 = Replace($LEGEND);
      my $L = eval("sprintf $LEG0");
      $LEG{$fout} = $L;
      if ( $PlotDatFiles !~ $fout ) {
#printf("%s\n", $fout); #koko
          $PlotDatFiles = $PlotDatFiles." ".$fout;
      }
    } else {
        printf("$PROG> Warning. Required data do not exist in $datfile\n");
    }

};
