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
