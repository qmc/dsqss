
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

  std::string key;
  double value  ,s1;
  double error;

  void set_key(std::string const &s) { key = s; };

  void reset() {
    key = "undefined"; 
    value = 0.0;
    error = 0.0;
  };

  Estimate() { reset(); };

  Estimate(const char* k) {reset(); set_key(k); };
  Estimate(std::string const& k) {reset(); set_key(k); };

  void operator =( Estimate E ) {
    value = E.value;
    error = E.error;
  };

  void dump() {
    printf(" %6s = %16.6e %16.6e\n", key.c_str(), value, error);
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

  Accumulator(const char* k) : Estimate::Estimate(k) {
    reset();
  };
  Accumulator(std::string const& k) : Estimate::Estimate(k) {
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
  void acc_put_value(double x){value=x;}
  double acc_get_error(){return s2;}
  void acc_put_error(double x, int nx){
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
#if defined(CF) || defined(CK) || defined(SF)
  XML::Block X;
#endif

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

#if defined(CF) ||  defined(CK) || defined(SF)
  int  Ntau;
  int Ntau1;
  int NKMAX;
  double dtau;
#endif
#if defined(CK) || defined(SF)
  void read_sf();
  double** COSrk;
  double** SINrk;
#endif
#ifdef CF

  int  NCF;
  int* DCF;
  int** ICF;//ICF[NSPIN][NSPIN];
  
  Accumulator** ACC4CF;  // accumurator for correlation functions
  Accumulator** PHY4CF;  // accumurator of set averages

  void show4CF(FILE*);

  void accumulate4CF(double** vx);

  void read_cf();

#endif
#ifdef SF
  
  int  NSF;
  
  Accumulator** ACC4SF;  // accumurator for correlation functions
  Accumulator** PHY4SF;  // accumurator of set averages

  double** counter4SFC;
  double** counter4SFS;
  double** sfsamp;
  
  void show4SF(FILE*);

  void accumulate4SF(double** vx);

  void reset4SF(){

    for(int isf=0; isf<NSF; isf++){
      for(int it=0; it<Ntau1; it++){
	counter4SFC[isf][it] = 0.0;
	counter4SFS[isf][it] = 0.0;
      }
    }

  };

#endif
#ifdef CK
  
  int  NCK;
  
  Accumulator** ACC4CK;  // accumurator for correlation functions
  Accumulator** PHY4CK;  // accumurator of set averages
  
  void show4CK(FILE*);

  void accumulate4CK(double** vx);

#endif



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


#ifdef CF
  if (DEBUG) printf("Measurement::Measurement> Reading %s\n",P.CFINPFILE );
  X.initialize( P.CFINPFILE );
  read_cf();
#endif

#ifdef SF
  if (DEBUG) printf("Measurement::Measurement> Reading %s\n",P.SFINPFILE );
  X.initialize( P.SFINPFILE );
#endif

#ifdef CK
  if (DEBUG) printf("Measurement::Measurement> Reading %s\n",P.CKINPFILE );
  X.initialize( P.CKINPFILE );
#endif

#if defined(CK) || defined(SF)
  read_sf();
#endif

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

#ifdef CF
  for(int i=0;i<NCF;i++){
    delete [] ACC4CF[i];
    delete [] PHY4CF[i];
  }

  delete [] ACC4CF;
  delete [] PHY4CF;
  
  for(int i=0;i<LAT.NSITE;i++)
    delete [] ICF[i];
  delete [] ICF;

  delete [] DCF;
#endif

#if defined(CK) || defined(SF)
  for(int i=0;i<LAT.NSITE;i++)
    delete [] COSrk[i];
  delete [] COSrk;

  for(int i=0;i<LAT.NSITE;i++)
    delete [] SINrk[i];
  delete [] SINrk;
#endif

#ifdef SF
  for(int i=0;i<NSF;i++){
    delete [] ACC4SF[i];
    delete [] PHY4SF[i];
  }

  delete [] ACC4SF;
  delete [] PHY4SF;

  for(int isf=0; isf<NSF+1; isf++){
    delete [] counter4SFC[isf] ; 
    delete [] counter4SFS[isf] ;
    delete [] sfsamp[isf] ;
    
  }
  delete [] counter4SFC; 
  delete [] counter4SFS;
  delete [] sfsamp; 
#endif

#ifdef CK
  for(int i=0;i<NCK;i++){
    delete [] ACC4CK[i];
    delete [] PHY4CK[i];
  }

  delete [] ACC4CK;
  delete [] PHY4CK;
#endif
  
}

//======================================================================

inline void Measurement::setinit() {

  for (int i=0; i<NACC; i++) ACC[i].reset();
#ifdef CF
  for (int i=0; i<NCF;  i++)
    for (int itau=0; itau<Ntau;  itau++)
      ACC4CF[i][itau].reset();
#endif
#ifdef SF
  for (int i=0; i<NSF;  i++)
    for (int itau=0; itau<Ntau;  itau++)
      ACC4SF[i][itau].reset();
#endif
#ifdef CK
  for (int i=0; i<NCK;  i++)
    for (int itau=0; itau<Ntau;  itau++)
      ACC4CK[i][itau].reset();
#endif

}

//======================================================================

inline void Measurement::summary() {

  for (int i=0; i<NPHY; i++) {
    PHY[i].average();
  }

#ifdef CF
  for (int i=0; i<NCF;  i++)
    for (int itau=0; itau<Ntau;  itau++)
      PHY4CF[i][itau].average();
#endif
#ifdef SF
  for (int i=0; i<NSF;  i++)
    for (int itau=0; itau<Ntau;  itau++)
      PHY4SF[i][itau].average();
#endif
#ifdef CK
  for (int i=0; i<NCK;  i++)
    for (int itau=0; itau<Ntau;  itau++)
      PHY4CK[i][itau].average();
#endif

}

//======================================================================

inline void Measurement::show(FILE* F) {
  for (int i=0; i<NPHY; i++) {
    fprintf(F,"R %-6s = %16.10e %16.10e\n", 
	    PHY[i].key.c_str(), PHY[i].value, PHY[i].error);
  }
}

#ifdef CF

inline void Measurement::show4CF(FILE* F) {
  for (int i=0; i<NCF; i++) 
    for (int it=0; it<Ntau; it++) 
      fprintf(F,"R C%dt%d = %16.10e %16.10e\n", i , it , PHY4CF[i][it].value , PHY4CF[i][it].error );
  
};

void Measurement::accumulate4CF(double** vt){
  for (int icf=0; icf<NCF; icf++) 
    for(int it=0;it<Ntau;it++)
      ACC4CF[icf][it].accumulate(vt[icf][it]);
};

#endif
#ifdef SF

inline void Measurement::show4SF(FILE* F) {
  
  for (int i=0; i<NSF; i++){
    for (int it=0; it<Ntau; it++){
      fprintf(F,"R C%dt%d = %16.10e %16.10e\n", i , it , PHY4SF[i][it].value , PHY4SF[i][it].error );
    }
    fprintf(F,"\n") ;
  }
  
};

void Measurement::accumulate4SF(double** vt){
  for (int isf=0; isf<NSF; isf++) 
    for(int it=0;it<Ntau;it++)
      ACC4SF[isf][it].accumulate(vt[isf][it]);
};

#endif
#ifdef CK

inline void Measurement::show4CK(FILE* F) {
  
  for (int i=0; i<NCK; i++){
    for (int it=0; it<Ntau; it++){
      fprintf(F,"R C%dt%d = %16.10e %16.10e\n", i , it , PHY4CK[i][it].value , PHY4CK[i][it].error );
    }
    fprintf(F,"\n") ;
  }
  
};

void Measurement::accumulate4CK(double** vt){
  for (int isf=0; isf<NCK; isf++) 
    for(int it=0;it<Ntau;it++)
      ACC4CK[isf][it].accumulate(vt[isf][it]);
};

#endif

//======================================================================

inline void Measurement::dump() {
  printf("\n");
  printf(">>> Dumping accumurators :\n");
  for (int i=0; i<NACC; i++) {
    printf(" ACC[%2d] = %16.6f\n", i, ACC[i].value);
  }
  printf("\n");
  for (int i=0; i<=NPHY; i++) {
    printf(" %s = phy[%2d] = %16.6f %16.6f\n", PHY[i].key.c_str(), i, PHY[i].value, PHY[i].error);
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

#ifdef CF
void Measurement::read_cf(){
  
  Ntau     = X["Ntau"].getInteger();
  int Nline = X["NumberOfElements"].getInteger();
#ifdef NORM
  dtau = 1.0 /((double) Ntau);
#else
  dtau = LAT.BETA /((double) Ntau);
#endif

  NCF       = X["NumberOfKinds"].getInteger();
  ACC4CF    = new Accumulator* [NCF];
  PHY4CF    = new Accumulator* [NCF];
  for(int i=0;i<NCF;i++){
    ACC4CF[i]    = new Accumulator [Ntau];
    PHY4CF[i]    = new Accumulator [Ntau];
    for(int it =0 ;it<Ntau;it++){
      ACC4CF[i][it].reset();
      PHY4CF[i][it].reset();
    }
  }
  
  DCF = new int [NCF];
  for(int i = 0; i <NCF; i++)
    DCF[i]  = 0;

  ICF = new int *[LAT.NSITE];
  for(int i=0;i<LAT.NSITE;i++){
    ICF[i] = new int[LAT.NSITE];
    for(int j=0;j<LAT.NSITE;j++)
      ICF[i][j]=NCF;
  }

  int count = 0;
  for (int i=0; i<X.NumberOfBlocks(); i++) {
    XML::Block& B = X[i];
    if ( B.getName() == "CF" ) {
      int icf_  = B.getInteger(0);
      int isite = B.getInteger(1);
      int jsite = B.getInteger(2);
      //cout<<isite<<"-"<<jsite<<endl;
      ICF[isite][jsite]=icf_;
      DCF[icf_]++;
      count++;
    }
  }

  if(count != Nline){
    cout<<"ERROR Nline( "<<Nline<<" ) != count( "<<count<<" )"<<endl;
    exit(0);
  }
};
#endif

#if defined(SF) || defined(CK)
void Measurement::read_sf(){

  Ntau1     = X["Ntau"].getInteger();
  int Nline = X["NumberOfElements"].getInteger();
#ifdef NORM
  dtau = 1.0 /((double) Ntau1);
#else
  dtau = LAT.BETA /((double) Ntau1);
#endif

  Ntau      = X["CutoffOfNtau"].getInteger();
  NKMAX     = X["NumberOfInverseLattice"].getInteger();

  COSrk = new double*[LAT.NSITE];
  for (int i=0;i<LAT.NSITE;i++){
    COSrk[i]=new double[NKMAX];
  }
  SINrk = new double*[LAT.NSITE];
  for (int i=0;i<LAT.NSITE;i++){
    SINrk[i]=new double[NKMAX];
  }

  int count = 0;
  for (int i=0; i<X.NumberOfBlocks(); i++) {
    XML::Block& B = X[i];
    if ( B.getName() == "SF" ) {
      double COSrk_  = B.getDouble(0);
      double SINrk_  = B.getDouble(1);
      int isite = B.getInteger(2);
      int ksite = B.getInteger(3);
      //cout<<isite<<"-"<<jsite<<endl;
      COSrk[isite][ksite]=COSrk_;
      SINrk[isite][ksite]=SINrk_;
      count++;
    }
  }


#ifdef SF
  NSF       = NKMAX;
  ACC4SF    = new Accumulator* [NSF];
  PHY4SF    = new Accumulator* [NSF];
  for(int i=0;i<NSF;i++){
    ACC4SF[i]    = new Accumulator [Ntau];
    PHY4SF[i]    = new Accumulator [Ntau];
    for(int it =0 ;it<Ntau;it++){
      ACC4SF[i][it].reset();
      PHY4SF[i][it].reset();
    }
  }
  
  counter4SFC  = new double* [NSF+1];
  counter4SFS  = new double* [NSF+1];
  sfsamp  = new double* [NSF+1];
  for(int isf=0;isf<NSF+1;isf++){
    counter4SFC[isf] = new double [Ntau1];
    counter4SFS[isf] = new double [Ntau1];
    sfsamp[isf] = new double [Ntau];
  }

#endif
  

#ifdef CK
  NCK       = NKMAX;
  ACC4CK    = new Accumulator* [NCK];
  PHY4CK    = new Accumulator* [NCK];
  for(int i=0;i<NCK;i++){
    ACC4CK[i]    = new Accumulator [Ntau];
    PHY4CK[i]    = new Accumulator [Ntau];
    for(int it =0 ;it<Ntau;it++){
      ACC4CK[i][it].reset();
      PHY4CK[i][it].reset();
    }
  }
#endif
  
  if(count != Nline){
    cout<<"ERROR Nline( "<<Nline<<" ) != count( "<<count<<" )"<<endl;
    exit(0);
  }


};
#endif
