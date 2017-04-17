//######################################################################
//####  MPI Header Files
//######################################################################

#ifdef MULTI
#include <mpi.h>
int N_PROC; // the total number of processors
int I_PROC; // the processor id of the current processor
#endif
//#include <iostream.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//######################################################################
//####  Global Variables
//######################################################################

FILE* FERR = stdout;


//######################################################################
//####  Graphics Header Files
//######################################################################

#ifdef GRAPHIC
#include <eggx.h>
#endif

//######################################################################
//####  Random Number Generator Header Files
//######################################################################

#include "random.h"
Random RND;

//######################################################################
//####  Local Header Files
//######################################################################

inline void abort(const char* s="") { printf("%s\n",s); exit(0); }

#ifdef DEB
bool DEBUG=true;
#else
bool DEBUG=false;
#endif

bool ALERT=false;

#include "name.h"
#include "io.h"
#include "array.h"
#include "link.hpp"
#include "parameter.hpp"
#include "algorithm.hpp"
#include "objects.hpp"
#include "lattice.hpp"
#include "measure.hpp"

//######################################################################

class Simulation {

private:

  int ISET;
  int IMCSE;
  int IMCSD;
  int IMCS;
  int ICYC;

  int ISETstart;
  int IMCSDstart;
  int IMCSstart;

  double AMP; // the average mean path of the worm per cycle
  bool EndOfCycle; // set to be true when the cycle ends

public:
  int np;
  Simulation(Parameter& P0);
  ~Simulation();

  void calc();

  Parameter& P;
  Algorithm ALG;

  //Algorithm alg[8];
  //Algorithm* alg;
  Lattice LAT;
  Measurement MSR;
  Worm W;
  
  Vertex V4REF;

  void reset_counters();
  // === edit rep sakakura ===
  //void Set();
  void Set(int index);
  // === edit rep sakakura ===
  void Set();
  void set_NCYC();

  // == edit sakakura ==
  //void CalcRepEx();
  void RepChangeF(Algorithm& A0);
  void RepChangeB(double h0);
  // == edit sakakura ==

  void Sweep();
  double Cycle();
  bool PlaceWorm();
  double MoveHead ();

  double UP_ONESTEP();
  double DOWN_ONESTEP();


#if defined(CF) || defined(CK)
  void   Sweep_CF();
  double Cycle_CF();
  double MoveHead_CF();
  double UP_ONESTEP_CF();
  double DOWN_ONESTEP_CF();

  
  double tail_tau;
  int    tail_site;
#endif

#ifdef CF
  
  int**    counter4CF;
  double** cfsamp;

  int count4CF(int icf, double tT ,double bT ){
    //#ifdef NORM
    //    double dtau = 1.0 /((double) MSR.Ntau1);
    //#else
    //    double dtau = LAT.BETA /((double) MSR.Ntau1);
    //#endif

    double dtau=MSR.dtau;
    
    double bTr;
    double tTr;
    
    bTr= tail_tau - tT;
    tTr= tail_tau - bT;
#ifdef NORM
    bTr +=(double)((bTr < -1.0e-20)); // [0:1)
    tTr +=(double)((tTr < -1.0e-20)); // (0:1]
#else
    bTr +=(double)((bTr < -1.0e-20)) * LAT.BETA; // [0:BETA)
    tTr +=(double)((tTr < -1.0e-20)) * LAT.BETA; // (0:BETA]
#endif   
    int bTri = (int)(bTr / dtau)  ;
    int tTri = (int)(tTr / dtau)  ; //Relative bottom (top) time integer

    if(  tTr < bTr ){
      for(int ktau = 0         ; ktau<= tTri   ; ktau++){ counter4CF[icf][ktau*( ktau < MSR.Ntau ) ] += 1;}
      for(int ktau = MSR.Ntau-1; ktau >= bTri+1; ktau--){ counter4CF[icf][ktau*( ktau < MSR.Ntau ) ] += 1;}
    }else{
      for(int ktau = bTri+1;     ktau <= tTri  ; ktau++){ counter4CF[icf][ktau*( ktau < MSR.Ntau ) ] += 1;}
    }
    
    return 0;
  };
  
  void reset4CF(){

    for(int icf=0; icf<MSR.NCF; icf++){
      for(int it=0; it<MSR.Ntau; it++){
	counter4CF[icf][it] = 0;
	cfsamp[icf][it]     = 0.0;
      }
    }

  };


#endif

#ifdef CK

  double** counter4CKC;
  double** counter4CKS;
  double** cksamp;


  int count4CK(int s, double tT ,double bT ){
    //#ifdef NORM
    //    double dtau = 1.0 /((double) MSR.Ntau1);
    //#else
    //    double dtau = LAT.BETA /((double) MSR.Ntau1);
    //#endif

    double dtau=MSR.dtau;
    
    double bTr;
    double tTr;
    
    bTr= tail_tau - tT;
    tTr= tail_tau - bT;
#ifdef NORM
    bTr +=(double)((bTr < -1.0e-20)); // [0:1)
    tTr +=(double)((tTr < -1.0e-20)); // (0:1]
#else
    bTr +=(double)((bTr < -1.0e-20)) * LAT.BETA; // [0:BETA)
    tTr +=(double)((tTr < -1.0e-20)) * LAT.BETA; // (0:BETA]
#endif   
    int bTri = (int)(bTr / dtau)  ;
    int tTri = (int)(tTr / dtau)  ; //Relative bottom (top) time integer

    if(  tTr < bTr ){
      for(int ick=0;ick<MSR.NKMAX;ick++){      
	for(int ktau = 0         ; ktau<= tTri   ; ktau++){
	  counter4CKC[ick][ktau*( ktau < MSR.Ntau ) ] += MSR.COSrk[s][ick];
	}
	for(int ktau = MSR.Ntau-1; ktau >= bTri+1; ktau--){ 
	  counter4CKC[ick][ktau*( ktau < MSR.Ntau ) ] += MSR.COSrk[s][ick];
	}
      }
    }else{
      for(int ick=0;ick<MSR.NKMAX;ick++){
	for(int ktau = bTri+1;     ktau <= tTri  ; ktau++){
	  counter4CKC[ick][ktau*( ktau < MSR.Ntau ) ] += MSR.COSrk[s][ick];
	}
      }
    }
    
    return 0;
  };
  

  void reset4CK(){

    for(int ick=0; ick<MSR.NCK; ick++){
      for(int it=0; it<MSR.Ntau1; it++){
	counter4CKC[ick][it] = 0.0;
	//counter4CKS[ick][it] = 0.0;
      }
    }

  };


#endif

  void Check();
  void dump();
  void dump(const char*);

#ifdef GRAPHIC
  int g_win;
  void g_init();
  void g_draw();
  void g_clear();
#endif

#ifdef CHAINJOB
private:
  bool isChainjob;
  bool isEnd; //true->exit(0); false-> read cjob.bin 
  char CJOBFILE[32];
  ofstream* cjobout;
  ifstream* cjobin;
  void BinaryIO();
  void read_cjob();
  void write_cjob();
  void end_cjob();
  void end_job();
#endif

};
