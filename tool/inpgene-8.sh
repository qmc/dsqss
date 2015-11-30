#!/bin/bash
mpipath=/usr/bin 

Machine=psi

hpath=$PWD
if [ -z "$WORM_HOME" ]; then
echo "Please set the environmental variables by bin/wormvars.sh."
echo "In worm directory, type "'"'"source  ./bin/wormvars.sh"'"'"."
exit 1 
fi
rootpath=${WORM_HOME}

inpfile=qmc.inp

dver=1.1.17
pver=1.1.3

dsqsspath=$rootpath/dsqss/dsqss-$dver/src
pmwapath=$rootpath/pmwa/pmwa-$pver/src

dbinpath=$rootpath/bin
pbinpath=$rootpath/bin

PFILE=param.dat
HFILE=ham.dat

PMODE="s" ## if submit with mpi "p", else "s" for parameter parallelization of dsqss.
if [ -z "$PMODE" ]; then
    PMODE="p"
elif [ $PMODE != "p" -a $PMODE != "s" ];then
    PMODE="p"
fi
echo "PMODE=$PMODE"

OPT=$1
if [[ "$@" =~ (^| )-D($| ) ]] ; then
  DEBUG_MODE="y"
else
  DEBUG_MODE="n"
fi
if [[ "$@" =~ (^| )-R($| ) ]] ; then
  QSUBMIT="on" ## if "on", automatically submit job.
else
  QSUBMIT="off" ## if "on", automatically submit job.
fi


MODEL="spin" #spin or boson
if [ -z "$MODEL" ]; then
    MODEL="spin"
elif [ $MODEL != "spin" -a $MODEL != "boson" ];then
    MODEL="spin"
fi
echo "MODEL=$MODEL"

############################
set_parameter_hs()
{

    echo "> defining a Heisenberg model"
    echo ">   M ... the number of bosons on each site"
    echo ">         (e.g., M = 1,2,3 .... for  S= 1/2, 1, 3/2 ...)"
    echo ">   J ... the coupling constant (e.g., -1.0 for antiferromagnets)"
    echo ">   H ... the magnetic field per site."
    echo ">"
    echo "> M J H ?"

    while true; do

    read MM JJ HH

    if [ $DEBUG_MODE = "y" ] ; then 
    echo ""
    echo "> ***********************************"
    echo "> INPUT:"
    echo "> M = $MM ," 
    echo "> J = $JJ ," 
    echo "> H = $HH ,"
    echo ">  have been set." 
    echo "> ***********************************"
    echo ""
    fi

    if [ -n "$MM" -a -n "$JJ" -a -n "$HH" ]; then 	
	break
    else
	echo "> THREE values are required !" 
	echo ">    ex) 1 4 2" 
    fi

    done

    echo ""

}

set_parameter_lat()
{
    echo "> defining a hyper cubic lattice"
    echo ">   D ... the dimension."
    echo ">   L ... the liner size of the lattice"
    echo ">         (must be an even integer.)"
    echo ">   B ... the inverse tempereture."
    echo ">"
    echo "> D L B ? "

    while true; do

    read DD LL BB

    if [ $DEBUG_MODE = "y" ] ; then 
    echo ""
    echo "> ***********************************"
    echo "> INPUT:"
    echo "> D = $DD," 
    echo "> L = $LL," 
    echo "> B = $BB,"
    echo "> have been set." 
    echo "> ***********************************"
    echo ""
    fi
    
    if [ -n "$DD" -a -n "$LL" -a -n "$BB" ]; then
	if [ $DD -gt "0" -a $DD -lt "4" ]; then

	    Lp=`expr $LL % 2`
	    if [ $LL -gt "1" -a $Lp -eq "0"  ]; then
      		X=`echo "$BB > 0.0" | bc`
		if [ $X -eq "1" ]; then 
		    break
		else
		    echo ">The inverse temperature must be B > 0.0 !" 
		    echo ">Input D L B." 
		fi
	    else
		echo ">L must be a positive even number."
		echo ">Input D L B."
	    fi
	else
	    echo ">The dimension must be  1<= D <=3."
	    echo ">Input D L B." 
	fi
    else
	echo "> THREE values are required!" 
	echo ">ex) 2 4 2" 
    fi

    done

    echo ""

}

set_parameter_bhp()
{
    echo "> Hardcore Boson model"
    echo "> J ... the hopping constat."
    echo "> V ... the nearest neighbor interaction."
    echo "> M ... the chemical potential in the Hamiltonian."
    echo ">"
    echo "> J V M ?"
 
    while true; do

    read JJ VV HH 

    if [ $DEBUG_MODE = "y" ] ; then 
    echo ""
    echo ">***********************************"
    echo "> INPUT:"
    echo "> J = $JJ," 
    echo "> V = $VV," 
    echo "> M = $HH,"
    echo ">  have been set."
    echo ">***********************************"
    echo ""
    fi

    if [ -n "$JJ" -a -n "$VV" -a -n "$HH" ]; then 	
	break
    else 
	echo ">Please input THREE values for J V M!" 
	echo ">ex) 1 4 2" 
    fi

    done

    echo ""
}

set_parameter_xxzp()
{
    echo "> defining a S=1/2 XXZ model"
    echo "> Jxy ... the coupling constat for x and y."
    echo "> Jz  ... the coupling constat for z."
#   echo "> Fxy ... the staggerd magnetic field in the pair Hamiltonian for easy plane."
    echo "> Hz  ... the magnetic field in the pair Hamiltonian for easy axis."
    echo ">"
    echo "> Jxy Jz Hz?"

 
    while true; do

    read JJ VV HH

   if [ $DEBUG_MODE = "y" ] ; then 
    echo ""
    echo ">***********************************"
    echo "> INPUT:"
    echo "> Jxy = $JJ," 
    echo "> Jz  = $VV," 
    echo "> Hz  = $HH,"
    echo ">  have been set."
    echo ">***********************************"
    echo ""
   fi

    if [ -n "$JJ" -a -n "$VV" -a -n "$HH" ]; then 	
	break
    else 
	echo ">Please input THREE values for Jxy Jz Hz!" 
	echo ">ex) 1 4 2" 
    fi

    done

    X=`echo "$JJ < 0.0" | bc`
    if [ $X -eq "1" ]; then
	JJ=`echo "scale=8; -1*$JJ"  | bc`
    fi

    echo ""
}

set_worm_source()
{

    ### Input worms ### 
    echo "> G ...The source fields for introducing multiple-worms."
    echo ">      Or staggerd transverse field when your model has the term."
    echo ">"
    echo "> G ?  [> 0.0]  (If no input, a default value is set.)"
    

    while true; do
	read GG;

	if [ -z "$GG" ]; then	    
	    GG=`echo "scale=5; 1.0/$LL^$DD"  | bc`
	    break
	else
	    X=`echo "$GG > 0.0" | bc`
	    if [ $X -eq "1" ]; then
		break
	    else
		echo "> The weight must be larger than 0.0."
		echo "> G ?" 
	    fi
        fi

    done

    echo ""

}


set_parameter_latp()
{
    echo "> defining a hyper cubic lattice"
    echo "> D ... the dimension."
    echo "> L ... the liner size of the lattice."
    echo "> B ... the inverse of tempereture."
    echo "> NL... the number of parallelization for spatial axes."
    echo ">       [ 1 <= NL < L ]"
    echo "> NB... the number of parallelization for the temporal axis."
    echo ">       [ 1 <= NB  ]"
#   echo "> NFIELD ..the number of external fields."
#   echo ">          for example, random potential etc..."
#   echo ">          if no, here is 0."

    echo "> D L B NL NB?"


    while true; do

	read DD LL BB NL NB

	NFIELD="0"
	if [ $DEBUG_MODE = "y" ] ; then 
	echo ""
	echo ">***********************************"
	echo "> INPUT:"
	echo "> D  = $DD ," 
	echo "> L  = $LL ," 
	echo "> B  = $BB ,"
	echo "> NL = $NL ,"
	echo "> NB = $NB ,"
	echo ">  have been set." 
	echo ">***********************************"
	fi

    	if [ -n "$NL" -a -n "$NB" -a -n "$DD" -a -n "$LL" -a -n "$BB" ]; then
	    if [ $DD -gt "0" -a $DD -lt "4" ]; then

		Lp=`expr $LL % 2`
		if [ $LL -gt "1" -a $Lp -eq "0"  ]; then
		    if [ $NL -gt "0" -a $NL -lt $LL ]; then 
			if [ $NB -gt "0" ]; then
			    
			    DL=`expr $LL % $NL`
			    if [ $DL -eq "0" ]; then
				
				X=`echo "$BB > 0.0" | bc`
				if [ $X -eq "1" ]; then 
				    break
				else
				    echo ">The inverse temperature must be B > 0.0 !" 
				    echo ">Input D L B NL NB."
				    echo ">ex) 1 4 8.0 2 2" 
				fi
			    else
				echo ">L/NL Must be integer!"
				echo ">Input D L B NL NB."
				echo ">ex) 1 4 8.0 2 2" 
			    fi
			    echo ">The number of parallelization for B-axis must be NB > 0!" 
			    echo ">(If no parallization for B, NB=1.)"
			    echo ">Input D L B NL NB."
			    echo ">ex) 1 4 8.0 2 2"
			fi
		    else
			echo ">The number of parallelization for L-axis must be 0 < NL < L!" 
			echo ">(If no parallization for L, NL=1.)"
			echo ">Input D L B NL NB."
			echo ">ex) 1 4 8.0 2 2"
		    fi
		else
		    echo ">L must be a positive even number."
		    echo ">Input D L B NL NB."
		    echo ">ex) 1 4 8.0 2 2"
		fi
	    else
		echo ">The dimension must be  1<= D <=3."
		echo ">Input D L B NL NB."
	    fi
	else
	    echo ">Please input FIVE values for D L B NL NB!" 
	    echo ">ex) 2 4 2 2 2" 
	fi

    done


    SS=`echo "scale=8; $LL/$NL"  | bc`
    if [ $DD -eq "2" ]; then SS=`echo "scale=8; $SS*$SS"  | bc`; fi
    if [ $DD -eq "3" ]; then SS=`echo "scale=8; $SS*$SS*$SS"  | bc`; fi
    dB=`echo "scale=8; $BB/$NB"  | bc`
    dV=`echo "scale=8; $SS*$dB"  | bc`

    if [ $DEBUG_MODE = "y" ] ; then 
	    echo "> Domain size for the B-axis is $dB."
	    echo "> Domain size for a L-axis is $SS."
	    echo "> Domain time-space volume is $dV."
    fi

    echo ""

}

set_parameter_parallelization(){

    if [ $NPTP -gt "1" ]; then
	echo "> vary  parameter ( default=f )"
	echo "> b ( inverse temperature )"
	echo "> l ( system size )"
	echo "> j ( coupling constat )"
	echo "> f ( magnetic field )"
	if [ $NPNT -gt "1" ]; then
	    echo "> g ( worm source field )"
	fi
	while true; do
	    read PARAM;
	    if [ -z "$PARAM" ]; then 
		PARAM='f'
		break;
	    elif [ $PARAM = "b" -o $PARAM = "l" -o $PARAM = "j" -o $PARAM = "f" -o $PARAM = "g" ]; then
		break;
	    else
		if [ $NPNT -gt "1" ]; then
		    echo "> Select from [ b, l, j, f, g ] !" ;
		else
		    echo "> Select from [ b, l, j, f ] !" ;
		fi 
	    fi
	done
        echo ""

	while true; do
	    if [ $PARAM = "l" ]; then
		echo "> step size? ( default=2 )"
		read dP
		if [ -z "$dP" ]; then 
		    dP=2;
		    break;
		    
		elif [ $dP -gt "1" -a `expr $dP % 2` -eq "0"  ]; then
		    break
		else
		    echo "> Must be a positive even number." 
		fi
	    else
		if [ $NPNT -gt "1" -a $PARAM = "g" ]; then
			echo "> step size? ( default=$GG )"
      		else
			echo "> step size? ( default=0.5 )"
		fi
		read dP
		if [ -z "$dP" ]; then 
		    if [ $NPNT -gt "1" -a $PARAM = "g" ]; then
			dP=$GG
		    else
			dP=0.5;
		    fi
		    break;
		elif [ `echo "$dP > 0.0" | bc` -eq "1" ]; then 
		    break
		else
		    echo "> Must be larger than 0.0." 
		fi
	    fi
	done
    fi

    echo ""

}

func_set_input_d()
{
    inpname=${1}
    LATNAME=${2}
    ALGNAME=${3}
    OUTNAME=${4}
    S_NP=${5}
    S_NPTP=${6}

    rm $inpname &>/dev/null
    echo "runtype   = $RUNTYPE" >$inpname
    echo "nmcse     = $NMCSE">>$inpname
    echo "nmcsd     = $NMCSD">>$inpname
    echo "nmcs      = $NMCS">>$inpname
    echo "nset      = $NSET">>$inpname
    echo "seed      = $SEED">>$inpname
    echo "nvermax   = $NVERMAX">>$inpname
    echo "nsegmax   = $NSEGMAX">>$inpname
    echo "nptp      = $S_NPTP">>$inpname  ##number of parallelization of a parameter
    echo "nrep      = $S_NP">>$inpname    ##Must be NREP=N_PROC in the source code.
    echo "algfile   = $ALGNAME">>$inpname
    echo "latfile   = $LATNAME">>$inpname
    echo "outfile   = $OUTNAME">>$inpname

}

func_set_input_p(){

    inpname=${1}
    LATNAME=${2}
    ALGNAME=${3}
    OUTNAME=${4}

    rm $inpname &>/dev/null

     echo  "RUNTYPE = $RUNTYPE" > $inpname
     echo  "NSET  = $NSET" >> $inpname    #the number of bins
     echo  "NMCSE = 100000" >> $inpname #thermalization
     echo  "NMCSD = $NMCSD" >> $inpname  #interval of bins
     echo  "NMCS  = $NMCS" >> $inpname  #the number of Monte Carlo samples.
     echo  "NTEST = 100000" >> $inpname #the number of Monte Carlo Steps for the preparation of N_cycle
     echo  "SEED  = $SEED" >> $inpname    #seed
     echo  "NC    = 0" >> $inpname     #0:all simulation, 1:preparation for Ncyc, 2:thermalization, 3:measurement

     echo  "NVERMAX = 10000000" >> $inpname #the max number of vertices. 
     echo  "NWORMAX = 1000" >> $inpname    # the max number of worms.

     echo  "algfile   = $ALGNAME" >> $inpname
     echo  "latfile   = $LATNAME" >> $inpname
     echo  "outfile   = $OUTNAME" >> $inpname

     echo  "CB      = 2" >> $inpname   #initial state of the configuration: 1-> checker board, 2->random, 0--> vacuum

     echo  "G       = $GG" >> $inpname  # the strenght of the transverse field for introducing worm.
     echo  "UBB     = 0" >> $inpname  #on-site interaction
     echo  "VB1     = $VV" >> $inpname  # nearest-neighbor interaction
     echo  "VB2     = 0" >> $inpname  #nextnearest neighbor
     echo  "tb      = $JJ" >> $inpname 
     echo  "MU      = $HH" >> $inpname  #chemical potential
     echo  "Wc      = 0" >> $inpname  #the external field
     echo  "PS      = 1" >> $inpname    #potential seed
     echo  "NMAX    = 1" >> $inpname    # the maximum number of bosons on a same site.

}

func_set_psi()
{

    S_NP=${1}
    
    jobscript="jobscript.sh"
    core="12"
    node=`expr $NP / $core`
    if [ `expr $NP % $core` -ne "0" ]; then
	node=`expr $node + 1`
    fi
    
    if [ $node -le "2" ]; then
	queue="small"
    elif [ $node -le "4" ]; then
	queue="middle"
    elif [ $node -le "8" ]; then
	queue="large"
    fi
    
    if [ $node -gt "8" ]; then
	echo ">psi has 96 cores."
	echo ">Your input core is ${NP}."
	echo ">Recommend to start over again!!"
        exit 1
    fi

    echo "#!/bin/sh"                   > $jobscript
    echo "#PBS -l nodes=$node:ppn=$core" >> $jobscript
    echo "#PBS -q $queue"                >> $jobscript
    echo "#PBS -N rtype_$RUNTYPE"        >> $jobscript
    echo "#PBS -j oe"                    >> $jobscript
    echo "export PATH=$WORM_HOME/bin:"'$PATH'     >> $jobscript
    echo 'cd $PBS_O_WORKDIR'             >> $jobscript

    if [ $NPTP -gt "1" -a $PMODE = "s" ];then
	XARG="ls -1 $inpfile.* | xargs -n 1 -P $NPTP"
	if [ $NPTS -eq "1" -a $NPNT -eq "1" ]; then
	    echo "$XARG $exefile" >> $jobscript
	else
	    echo "$XARG mpiexec -np $S_NP $exefile" >> $jobscript
	fi
    elif [ $NPTP -eq "1" -a $S_NP -eq "1" ]; then
	echo "$exefile $inpfile" >> $jobscript
    else
	echo "mpiexec -np $S_NP $exefile $inpfile" >> $jobscript
    fi
    
    if [ -z $QSUBMIT -o $QSUBMIT != "on" ]; then
	echo ">A job script "'"'"jobscript.sh"'"'" is generated in DATA/$OD."
	echo ">To execute your simulation,"
	echo "> qsub $jobscript"

    else
        qsub $jobscript
    fi


}

func_set_pc()
{

    S_NP=${1}

    if [ -z $QSUBMIT -o $QSUBMIT != "on" ]; then
	echo ">A job script "'"'"jobscript.sh"'"'" is generated in DATA/$OD."
	echo ">To execute your simulation, "
       
	if [ $NPTP -gt "1" -a $PMODE ="s" ];then
	    P_e=$(( $NPTP - 1 ))
	    for ii in `seq 0 $P_e`; do
		PPN=`printf "%03d\n" $ii`
		if [ $NPTS -eq "1" ]; then
		    echo "$exefile $inpfile.$PPN" 
		else
		    echo "mpirun -np $NPTS $exefile $inpfile.$PPN" 
		fi
	    done
	elif [ $NPTP -eq "1" -a $S_NP -eq "1" ]; then
	    echo "$exefile $inpfile" 
	else
	    echo "mpirun -np $S_NP $exefile $inpfile" 
	fi
	
    else
        if [ $NPTP -gt "1" -a $PMODE ="s" ];then
	    P_e=$(( $NPTP - 1 ))
	    for ii in `seq 0 $P_e`; do
		PPN=`printf "%03d\n" $ii`
		if [ $NPTS -eq "1" ]; then
		    nohup $exefile $inpfile.$PPN &
		else
		    nohup $mpipath/mpirun -np $NPTS $exefile $inpfile.$PPN & 
		fi
	    done
	elif [ $NPTP -eq "1" -a $S_NP -eq "1" ]; then
	    nohup $exefile $inpfile &
	else
	    nohup $mpipath/mpirun -np $S_NP $exefile $inpfile & 
	fi

    fi



}



####################################################################################
#####                                start here                                #####
####################################################################################
echo "*************************************************************"
echo "> Welcome to DSQSS/PMWA!"
echo "> DSQSS/PMWA package implements four kinds of parallelization for followings:"
echo ">"
echo "> 1. replica       (trivial)."
echo "> 2. domains of configuration space (nontrivial)."
echo "> 3. parameters     (trivial)."
echo "> 4. random number  (trivial)."
echo ">"
echo "> You can create input files and execute dsqss/pmwa"
echo "> w/o parallelization by mean of this navigater."
echo "> Please start from creating your project directory first."
echo "*************************************************************"
echo ""
echo "> project directory name? (default=test)"
read OD
if [ -z "$OD" ]; then 
    OD="test"
    PRODDIR=$OD
fi

ODIR=$PWD/DATA/$OD

if [ ! -d $ODIR ]; then
mkdir -p $ODIR
fi

cd $ODIR

echo ""
####################################
#####  select parallelization  #####
####################################
echo "> Input the degree of following parallelizations."
echo ">"
echo "> NREP  ... replica Monte Carlo"
echo "> NPNT ... nontrivial parallelization"
echo "> NPTS ... trivial parallelization for random seed"
echo "> NPTP ... trivial parallelization for parameter"
echo ">"
echo "> (NOTE : if NREP > 1, others must be 1.   )"
echo "> (NOTE : if NPTP = 1 no parallelization,  )"
echo "> (       NPTP > 1 the parallelization is available )"
echo "> (       but the number is dummy.         )"
echo ">"
echo "> NREP NPNT NPTS NPTP ? (default : 1 1 1 1)"
#-------- Parallelization 1 : reprica ---------
while true; do

    read NREP NPNT NPTS NPTP;


    if [ -z "$NREP" -a -z "$NPNT" -a -z "$NPTS" -a -z "$NPTP" ]; then
	NREP=1
	NPNT=1
	NPTS=1
	NPTP=1
	break
    fi

    if [ -z "$NREP" -o -z "$NPNT" -o -z "$NPTS" -o -z "$NPTP" ]; then
	echo "> FOUR values are required !" 
	echo ">    ex) 1 1 4 6" 
    else
	if [ $NREP -gt "0" -a $NPNT -gt "0" -a $NPTS -gt "0" -a $NPTP -gt "0" ]; then
	    NN=$(( $NPNT * $NPTS * $NPTP ))
	    if [ $NREP -gt "1" -a $NN -eq "1" ]; then
		break
	    elif [ $NREP -eq "1" -a $NN -ge "1" ]; then
		break
	    else
		echo "> When NREP > 1, others must be 1."
		echo "> NREP NPNT NPTS NPTP ?"
	    fi
	else
	    echo "> Must be larger than 0!" ;
	fi
    fi

    echo ""

done

if [ $DEBUG_MODE = "y" ] ; then 
    echo ""
    echo "> ***********************************"
    echo "> INPUT:"
    echo "> NREP = $NREP," 
    echo "> NPNT = $NPNT," 
    echo "> NPTS = $NPTS,"
    echo "> NPTP = $NPTP,"
    echo "> have been set." 
    echo "> ***********************************"
    echo ""
fi


echo ""

############################
#####  input parameter #####
############################

if [ $NPNT -eq "1" ]; then 

    if [ $NREP -eq "1" ]; then
    ############ conventional DLA #############
	RUNTYPE=0
	set_parameter_hs
	set_parameter_lat
	set_parameter_parallelization

    else
    ########## replica DLA ##########	
	echo ">  vary magnetic field (m) or temperature (t) [m/t] (default=m)."
	while true; do
	    read MORT;
	    if [ -z "$MORT" ]; then 
		MORT='m'
		break;
	    elif [ $MORT = "m" -o $MORT = "t" ]; then
		break;
	    else 
		echo "> Must be m or t!" ;
	    fi
	done

	if [ "$MORT" = 'm' ]; then  
         #### replica DLA with varying field ####
            RUNTYPE=1
	    set_parameter_hs;
	    echo "> magnetic-field interval? ( ex. 0.02 )."
	    while true; do
		read DH 
		if [ -n "$DH" ]; then
		    X=`echo "$DH > 0.0" | bc`
		    if [ $X -eq "1" ]; then 
			echo "> the magnetic interval is $DH."
			break
		    else
			echo "> choose a positive value!"
		    fi
		else
		    echo "> choose some value!"
		fi
	    done
	    set_parameter_lat
	else
        #### replica DLA with varying temperature ####
            RUNTYPE=2
	    set_parameter_hs;
	    set_parameter_lat;
	    echo "> inverse-temperature interval? ( ex. 0.02 )"
	    while true; do
		read DB
		if [ -n "$DB" ]; then
		    X=`echo "$DB > 0.0" | bc`
		    if [ $X -eq "1" ]; then 
			echo "> inverse-temperature interval is $DB."
			break
		    else
			echo "> choose a positive value!"
		    fi
		else
		    echo "> choose some value!"
		fi
	    done
	fi
    fi
else
############# PMWA ############ 
    RUNTYPE=3
    if [ $MODEL = "spin" ]; then
       set_parameter_xxzp;
    else
       set_parameter_bhp;
    fi

    set_parameter_latp;

    set_worm_source

    set_parameter_parallelization

fi

echo ""
################################
#####  create input files  #####
################################
#initialize
echo "> number of Monte Carlo step? (default=10000)"
read NMCS
if [ -z "$NMCS" ]; then 
    NMCS=10000
fi
NMCSE=$NMCS
NMCSD=$NMCS
NMCS=$NMCS
NSET=10
SEED=31415
NVERMAX=10000
NSEGMAX=10000

ALGFILE=algorithm.xml
LATFILE=lattice.xml
OUTFILE=sample.log

####===================== RUNTYPE = 0 (normal dla) ===========================###
# -- set parameter -- #
if [ $RUNTYPE -eq 0 ];then

    NP=$(( $NPTS * $NPTP * $NPNT ))

    JJ0=$JJ
    HH0=$HH
    LL0=$LL 
    BB0=$BB

    exefile=dla

############loop start

    if [ $NPTP -eq "1" ]; then

	FF=`echo "scale=8; 0.5*$HH/$DD" | bc`
	hamgen_H $MM $JJ $FF >> gen_xml.log
	lattgene $DD $LL $BB >> gen_xml.log
	dla_alg              >> gen_xml.log
	func_set_input_d $inpfile $LATFILE $ALGFILE $OUTFILE $NP $NPTP

    else

	P_e=$(( $NPTP - 1 ))

	for ii in `seq 0 $P_e`; do

	    if [ $PARAM = 'j' ]; then
		JJ=`echo "scale=8; $JJ0+$dP*$ii"  | bc`
	    elif [ $PARAM = 'f' ]; then
		HH=`echo "scale=8; $HH0+$dP*$ii"  | bc`
	    elif [ $PARAM = 'l' ]; then
		LL=`echo "scale=8; $LL0+$dP*$ii"  | bc`
	    elif [ $PARAM = 'b' ]; then
		BB=`echo "scale=8; $BB0+$dP*$ii"  | bc`
	    fi

            FF=`echo "scale=8; 0.5*$HH/$DD" | bc`
	    hamgen_H $MM $JJ $FF >> gen_xml.log
	    lattgene $DD $LL $BB >> gen_xml.log
	    dla_alg              >> gen_xml.log
	    
	    PPN=`printf "%03d\n" $ii`
	    mv $LATFILE lattice_${PPN}.xml
	    mv $ALGFILE algorithm_${PPN}.xml

	    if [ $PMODE = "s" ]; then
 		func_set_input_d $inpfile.$PPN lattice_${PPN}.xml algorithm_${PPN}.xml sample_${PPN}.log $NPTS 1
	    fi

	done

	if [ $PMODE = "p" ]; then
	    func_set_input_d $inpfile $LATFILE $ALGFILE $OUTFILE $NP $NPTP
	fi

    fi

    JJ=$JJ0
    HH=$HH0
    LL=$LL0 
    BB=$BB0

#########loop end


    rm -f $PFILE &>/dev/null
    echo "runtype      = $RUNTYPE" >$PFILE
    echo "project_dir  = $PRODDIR">>$PFILE
    echo "np_triv_s    = $NPTS">>$PFILE  ##for seed
    echo "np_triv_p    = $NPTP">>$PFILE  ##for parameter
    if [ $NPTP -gt "1" ]; then
	echo "triv_p       = $PARAM">>$PFILE 
	echo "triv_p_step  = $dP">>$PFILE 
    fi
    echo "np_nontriv   = $NPNT">>$PFILE
    echo "nrep         = $NREP">>$PFILE
    echo "M            = $MM">>$PFILE
    echo "J            = $JJ">>$PFILE 
    echo "H            = $HH">>$PFILE
    echo "D            = $DD">>$PFILE
    echo "L            = $LL">>$PFILE
    echo "B            = $BB">>$PFILE
    echo "nmcse        = $NMCSE">>$PFILE
    echo "nmcsd        = $NMCSD">>$PFILE
    echo "nmcs         = $NMCS">>$PFILE
    echo "nset         = $NSET">>$PFILE
    echo "seed         = $SEED">>$PFILE
    echo "nvermax      = $NVERMAX">>$PFILE
    echo "nsegmax      = $NSEGMAX">>$PFILE
    echo "algfile      = $ALGFILE">>$PFILE
    echo "latfile      = $LATFILE">>$PFILE
    echo "outfile      = $OUTFILE">>$PFILE

    rm -f $HFILE &>/dev/null
    echo "np_triv_p=$NPTP">>$HFILE  ##for parameter
    if [ $NPTP -gt "1" ]; then
	echo "triv_p=$PARAM">>$HFILE 
	echo "triv_p_step=$dP">>$HFILE 
    fi
    echo "M=$MM">>$HFILE
    echo "J=$JJ">>$HFILE 
    echo "H=$HH">>$HFILE
    echo "D=$DD">>$HFILE
    echo "L=$LL">>$HFILE
    echo "B=$BB">>$HFILE


fi
####===================== RUNTYPE = 1 (replica m) ===========================###
####===================== RUNTYPE = 2 (replica t) ===========================###
if [ $RUNTYPE -eq 1 -o $RUNTYPE -eq 2 ];then

    PMODE="s"    
    exefile=dla
    lattgene $DD $LL $BB >> gen_xml.log

    irep=0

    if [ $RUNTYPE -eq 1 ];then

	while [ $irep -lt $NREP ];do
 
	XH=`echo "scale=5; $HH+$irep*$DH "  | bc`

	echo "XH=$XH"
	
	XF=`echo "scale=8; 0.5*$XH/$DD" | bc`
	hamgen_H $MM $JJ $XF >> gen_xml.log
	dla_alg              >> gen_xml.log
	cp algorithm.xml $irep.xml

	irep=`expr $irep + 1`
	done
    else

	FF=`echo "scale=8; 0.5*$HH/$DD" | bc`
	hamgen_H $MM $JJ $FF >> gen_xml.log
	dla_alg              >> gen_xml.log
	
	while [ $irep -lt $NREP ]; do
	    cp algorithm.xml $irep.xml 
	    irep=`expr $irep + 1`
	done
	
    fi


    NP=$NREP

    rm $inpfile &>/dev/null
    echo "runtype = 1" >$inpfile
    echo "nmcse  = 10">>$inpfile
    echo "nmcsd  = 100">>$inpfile
    echo "nmcs   = 500">>$inpfile
    echo "nset   = 100">>$inpfile
    echo "nrep=$NREP">>$inpfile
    if [ $RUNTYPE -eq 1 ]; then
	echo "vh=$HH">>$inpfile
	echo "dh=$DH">>$inpfile
    else
	echo "vb=$BB">>$inpfile
	echo "db=$DB">>$inpfile
    fi
    echo "outfile=$OUTFILE">>$inpfile
    
    rm aaa &>/dev/null
    
fi
####===================== RUNTYPE = 3 (pmwa) ===========================###
# -- set parameter -- #
if [ $RUNTYPE -eq 3 ];then

    PMODE="s"
    Ntdiv=$NB
    Nxdiv=$NL
    
    if [ $DD -ge "2" ]; then
	 Nydiv=$Nxdiv 
    else
	Nydiv=1
    fi
    if [ $DD -eq "3" ]; then
	Nzdiv=$Nxdiv
    else
	Nzdiv=1
    fi
    
    NPNT=`expr $Ntdiv \* $Nxdiv \* $Nydiv \* $Nzdiv`
    NP=`expr $NPNT \* $NPTS \* $NPTP`
    if [ $NFIELD -gt 0 ]; then
	NP=`expr $NP \* $NFIELD`
    fi


    JJ0=$JJ
    HH0=$HH
    LL0=$LL 
    BB0=$BB
    GG0=$GG

    exefile=dla_P


############loop start

    if [ $NPTP -eq "1" ]; then

            #hamgen_B $JJ $VV $HH
	lattgene_P $DD $LL $BB $NL $NB $NFIELD >> gen_xml.log
            #pmwa_alg
 	func_set_input_p $inpfile lattice.xml algorithm.xml sample.log

    else

	P_e=$(( $NPTP - 1 ))

	for ii in `seq 0 $P_e`; do
	    
	    if [ $PARAM = 'j' ]; then
		JJ=`echo "scale=8; $JJ0+$dP*$ii"  | bc`
	    elif [ $PARAM = 'f' ]; then
		HH=`echo "scale=8; $HH0+$dP*$ii"  | bc`
	    elif [ $PARAM = 'l' ]; then
		LL=`echo "scale=8; $LL0+$dP*$ii"  | bc`
	    elif [ $PARAM = 'b' ]; then
		BB=`echo "scale=8; $BB0+$dP*$ii"  | bc`
	    elif [ $PARAM = 'g' ]; then
		GG=`echo "scale=8; $GG0+$dP*$ii"  | bc`
	    fi
		
	        #hamgen_B $JJ $VV $HH
	    lattgene_P $DD $LL $BB $NL $NB $NFIELD >> gen_xml.log
                #pmwa_alg
	    PPN=`printf "%03d\n" $ii`
	    mv $LATFILE lattice_${PPN}.xml
	    
#	    if [ $PMODE = "s" ]; then
 	    func_set_input_p $inpfile.$PPN lattice_${PPN}.xml algorithm_${PPN}.xml sample_${PPN}.log
#	    fi

	done
       	
    fi

    JJ=$JJ0
    HH=$HH0
    LL=$LL0 
    BB=$BB0
    GG=$GG0

#########loop end


    rm -f $HFILE &>/dev/null
    echo "np_triv_p=$NPTP">>$HFILE  ##for parameter
    if [ $NPTP -gt "1" ]; then
	echo "triv_p=$PARAM">>$HFILE 
	echo "triv_p_step=$dP">>$HFILE 
    fi
    echo "V=$VV">>$HFILE
    echo "J=$JJ">>$HFILE 
    echo "H=$HH">>$HFILE
    echo "D=$DD">>$HFILE
    echo "L=$LL">>$HFILE
    echo "B=$BB">>$HFILE
    echo "NL=$NL">>$HFILE
    echo "NB=$NB">>$HFILE
    echo "G=$GG">>$HFILE
  

#    NPTP=$(( $NPTP * 3 ))

fi

#  $mpipath/mpirun -np $NP $exefile $inpfile &>dla.log

################################
#####        common        #####
################################

echo ""
echo ">*******************************************************************"
echo ">NP = $NP" 
echo ">Input files and the exective are successflly generated in the directory "'"'DATA/$OD.'"'
echo ">The following files are created in "'"'$ODIR.'"'
echo ">--------------------------"
#echo ">Executive file   : $exefile"
echo ">Input file       : $inpfile"
echo ">Associated files : *.xml"
echo ">--------------------------"

cp ${WORM_HOME}/bin/plot/plot.config .

if [ $Machine = "psi" ]; then

    echo ">Your machine is ${Machine}."

    if [ $NREP -gt 1 -o $PMODE = "p" ]; then
	func_set_psi $NP
    elif [ $PMODE = "s" ]; then
	func_set_psi $(( $NPTS * $NPNT))
    else
	func_set_psi 1
    fi

else ##################other machines

    if [ $NREP -gt 1 -o $PMODE = "p" ]; then
	func_set_pc $NP
    elif [ $PMODE = "s" ]; then
	func_set_pc $(( $NPTS * $NPNT))
    else
	func_set_pc 1
    fi



fi

echo ">*******************************************************************"
echo ""

cd $hpath
