#!/bin/bash

#====================
#====================
PARALLEL=YES
#COMPILER=GNU
COMPILER=INTEL
#COMPILER=FUJITSU

######################################################
# !!!!! Please check folloing compile commands !!!!! #
######################################################

case $COMPILER in 
    GNU)
    LIBS="-lblas -llapack "
    OPT="-O2"
    
    case $PARALLEL in
        YES)
        CXX=mpic++
        CC=mpic++
        MULTI=-DMULTI;;
        NO)
        CXX=g++
        CC=g++;;
        *)
        echo parameter error: PARALLEL = $PARALLEL
        echo PARALLEL should be YES or NO.
        exit -1;;
    esac;;
    FUJITSU)
    LIBS="-SSL2 "
    OPT="-Kfast -Ksimd=2 -O3 -Nlst "
    CRSS="--host --build "
    case $PARALLEL in
        YES)
        CXX=mpiFCCpx
        CC=mpifccpx
        MULTI=-DMULTI;;
        NO)
        CXX=FCCpx
        CC=fccpx;;
        *)
        echo parameter error: PARALLEL = $PARALLEL
        echo PARALLEL should be YES or NO.
        exit -1;;
    esac;;
    INTEL)
    LIBS=-mkl
    OPT="-O2 "
    case $PARALLEL in
        YES)
        CXX=mpicxx
        CC=mpicxx
        MULTI=-DMULTI;;
        NO)
        CXX=icpc
        CC=icpc;;
        *)
        echo parameter error: PARALLEL = $PARALLEL
        echo PARALLEL should be YES or NO.
        exit -1;;
    esac;;
    *)
    echo parameter error: COMPILER = $COMPILER
    echo COMPILER should be GCC, INTEL, or FUJITSU.
    exit -1;;
esac

######################################################
######################################################
#====================
ROOTDIR=$PWD
Dver=1.1.17
Pver=1.1.2
BINDIR=$PWD/bin
#====================
#====================

if [ ! -d $BINDIR ]; then
mkdir $BINDIR
fi

# ++++ installation directory ++++
DSQSS_INSTALL_DIR=$ROOTDIR/dsqss/dsqss-$Dver
PMWA_INSTALL_DIR=$ROOTDIR/pmwa/pmwa-$Pver
# ++++ installation directory ++++

cd $DSQSS_INSTALL_DIR
./runCtest "$ROOTDIR" "$CC" "$CXX" "$LIBS" "$OPT" "$MULTI" "$CRSS"
cd $ROOTDIR

cd $PMWA_INSTALL_DIR
./runCtest "$ROOTDIR" "$CC" "$CXX" "$LIBS" "$OPT" "$MULTI" "$CRSS"
cd $ROOTDIR

echo "WORM_HOME=$ROOTDIR"     >  $BINDIR/wormvars.sh
echo "export WORM_HOME"       >> $BINDIR/wormvars.sh
echo "PATH="'$WORM_HOME/bin:$PATH'>> $BINDIR/wormvars.sh
echo "PATH="'$WORM_HOME/bin/plot:$PATH'>> $BINDIR/wormvars.sh
echo "export PATH"            >> $BINDIR/wormvars.sh
