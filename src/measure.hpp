
#ifndef MEASURE_H
#define MEASURE_H

//######################################################################

#include <stdio.h>
#include "name.h"
#include "parameter.hpp"
#include "measure_specific.h"

//######################################################################

class Estimate {

public:

  const char* key;
  double value  ,s1;
  double error;

  void set_key(string s) { key = s.c_str(); };

  void reset() {
    key = "undefined"; 
    value = 0.0;
    error = 0.0;
  };

  Estimate() { reset(); };

  Estimate(char* k) { set_key(k); reset(); };

  void operator =( Estimate E ) {
    value = E.value;
    error = E.error;
  };

  void dump() {
    printf(" %6s = %16.6e %16.6e\n", key, value, error);
  };

};

//######################################################################

class Accumulator : public Estimate {

private:

  int n;
  //kota  double s1;
  //  double s1;
  double s2;

public:

  void reset() {
    n = 0;
    s1 = 0.0;
    s2 = 0.0;
  };

  Accumulator() : Estimate::Estimate() {
    reset();
  };

  Accumulator(char* k) : Estimate::Estimate(k) {
    reset();
  };

  ~Accumulator() {};

  void accumulate(double x) {
    n++;
    s1 += x;
    s2 += (x*x);
  };

  void accumulatex(double x) {
    n++;
    s1 += x;
    s2 += (x*x);
  };
  // ++++ edit rep sakakura ++++
  double acc_get_value(){return s1;}
  double acc_put_value(double x){value=x;}
  double acc_get_error(){return s2;}
  double acc_put_error(double x, int nx){
    error=x;
    //    cout <<"nx = " <<nx <<endl;
    //  value=s1;
    if (nx>2) {
      error = error - value * value;
      error /= (double)(nx-1);
      error = sqrt(error);
    } else{
      error = 0.0;}
  }

  // ++++ edit rep sakakura ++++

  void average() {
    value = 0.0;
    error = 0.0;
    if ( n > 0 ) {
      value = s1 / (double)n;
      error = s2 / (double)n;
      if ( n > 1 ) {
        error = error - value * value;
        error /= (double)(n-1);
        error = sqrt(error);
      } else {
        error = 0.0;
      }
    }
  };

};

//######################################################################

class Measurement {

private:

  int NACC;
  int NPHY;
  Parameter& P;
  Lattice& LAT;
  Algorithm& ALG; 
  double* Q;    //

public:

  Accumulator* ACC;  // accumurator of snapshot values
  Accumulator* PHY;  // accumurator of set averages
  // == edit rep sakakura ===
  Accumulator** PHYX;  // accumurator of set averages
  //  Accumulator* PHYE;  // accumurator of set averages


  double ediag;
  int  nkink;
  // == edit rep sakakura ===

  double EBASE;

  Measurement(Parameter& P0, Lattice& L, Algorithm& A);
  ~Measurement();
  void measure();
  void summary();
  void setinit();

  // == edit rep sakakura ===
  void summary_REP(int index);
  //  void dumpX(int i);
  void setsummary(int index);
  double get_data();
  double get_ediag();
  int get_nkink();
  //  void setsummary();
  // == edit rep sakakura ===

  void accumulate_length(double len);
  void dump();
  void show(FILE*);
  int imin(int n, double* t);


  //  double get_data();
};

//######################################################################
inline int Measurement::imin(int n, double* t) {
  if ( n==0 ) { printf("imin> Error. n = 0.\n"); exit(0); }
  int i = 0;
  double tmin = t[0];
  for ( int j=1; j<n; j++) {
    if ( t[j] < tmin ) {
      i = j;
      tmin = t[j];
    }
  }
  return i;
}
//######################################################################
inline Measurement::Measurement(Parameter& P0, Lattice& L, Algorithm& A) : 
  P(P0), LAT(L), ALG(A) 
{

  if (DEBUG) printf("\nMeasurement::Measurement> Start.\n");

  NACC = Specific::NACC; 
  NPHY = Specific::NPHY;
  Q   = new double[NPHY];
  ACC = new Accumulator[NACC];
  PHY = new Accumulator[NPHY];

  // == edit rep sakakura == 
  //  PHYE = new Accumulator[P.NREP];
  PHYX = new Accumulator*[NPHY];

  for (int i=0;i<NPHY;i++){
    PHYX[i]=new Accumulator[P.NREP];
  }

  for (int i=0; i<P.NREP; i++) { PHYX[1][i].reset(); }

  //  for (int i=0; i<P.NREP; i++) { PHYE[i].reset(); }
 
  // == edit rep sakakura == 


  for (int i=0; i<NACC; i++) { ACC[i].reset(); }

  for (int i=0; i<NPHY; i++) { PHY[i].reset(); }

  for (int i=0; i<NACC; i++) { ACC[i].set_key( Specific::ANAME[i] ); }

  for (int i=0; i<NPHY; i++) { PHY[i].set_key( Specific::PNAME[i] ); }

  // setting EBASE

  EBASE = 0.0;
  for (int i=0; i<LAT.NINT; i++) {
    InteractionProperty& IP = LAT.I(i).property();
    double eb = IP.EBASE;
    EBASE += eb;
  }

  if (DEBUG) printf("Measurement::Measurement> End.\n");

}

//======================================================================

inline Measurement::~Measurement() {
  //  printf("*** Destroying Measurement\n");
  delete [] Q;
  delete [] ACC;
  delete [] PHY;
  //======== edit rep sakakura =======//
  for (int i=0;i<NPHY;i++){
    delete [] PHYX[i];
  }
  delete [] PHYX;
  //  delete [] PHYE;
  //======== edit rep sakakura =======//
}

//======================================================================

inline void Measurement::setinit() {

  for (int i=0; i<NACC; i++) ACC[i].reset();
}

//======================================================================

inline void Measurement::summary() {

  for (int i=0; i<NPHY; i++) {
    PHY[i].average();
  }

}

//======================================================================

inline void Measurement::show(FILE* F) {
  for (int i=0; i<NPHY; i++) {
    fprintf(F,"R %-6s = %16.10e %16.10e\n", 
	    PHY[i].key, PHY[i].value, PHY[i].error);
  }
}

//======================================================================

inline void Measurement::dump() {
  printf("\n");
  printf(">>> Dumping accumurators :\n");
  for (int i=0; i<NACC; i++) {
    printf(" ACC[%2d] = %16.6f\n", i, ACC[i].value);
  }
  printf("\n");
  for (int i=0; i<=NPHY; i++) {
    printf(" %s = phy[%2d] = %16.6f %16.6f\n", PHY[i].key, i, PHY[i].value, PHY[i].error);
  }
}

//======================================================================

inline void Measurement::accumulate_length(double len) {
  ACC[Specific::LE1].accumulate(len);
}

//inline double Measurement::get_data(){
//  return PHY[1].value;
//};

#include "measure_specific.cc"

#endif

// ----------------------------------------------------------------------
// === edit rep sakakura ===

//inline void Measurement::dumpX(int i) {
//  printf ("index=%d,%12.5f %12.5f\n",i,PHYX[1][i].value,PHYX[1][i].error);

//    cout << PHYE[0].value <<" dumpx "<<i<< endl;
    
//      printf ("index=%d,%12.5f %12.5f\n",i,PHYE[i].value,PHYE[i].error);
//}
inline double Measurement::get_data(){
  return PHY[1].value;
}

inline double Measurement::get_ediag(){
  return ediag;
}

inline int Measurement::get_nkink(){
  return nkink;
}

inline void Measurement::summary_REP(int index) {

  double buf0[NPHY][P.NREP],buf1[NPHY][P.NREP];

  for (int i=0; i<P.NREP; i++){
    for (int j=0; j<NPHY; j++){
      buf0[j][i] = PHYX[j][i].acc_get_value();
    }
  }

  double size = P.NREP * NPHY;

#ifdef MULTI
  // gather value
#ifdef KASHIWA
  MPI_Allreduce(&buf0,&buf1,size,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Allreduce(&buf0,&buf1,size,MPI::DOUBLE,MPI::SUM);
#endif

#endif

  for (int i=0; i<NPHY; i++){
    PHY[i].acc_put_value(buf1[i][index]/P.NSET);
  }

  for (int i=0; i<P.NREP; i++){
    for (int j=0; j<NPHY; j++){
      buf0[j][i] = PHYX[j][i].acc_get_error();
      buf1[j][i] = 0;
    }     
  }

#ifdef MULTI
  // gather error
#ifdef KASHIWA
  MPI_Allreduce(&buf0,&buf1,size,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Allreduce(&buf0,&buf1,size,MPI::DOUBLE,MPI::SUM);
#endif
#endif

  for (int i=0; i<NPHY; i++){
    PHY[i].acc_put_error(buf1[i][index]/P.NSET,P.NSET);
  }

}
