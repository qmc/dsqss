
#ifndef LATTICE_H
#define LATTICE_H

//######################################################################

#include <stdio.h>
#include "io.h"
#include "objects.hpp"

//######################################################################

class Lattice {

public:

  class Edge{
    
  public:
    void init(int _d,int _a, int _b){
      
      bd=_d;
      A=_a;
      B=_b;
      
    };
    
    int A;
    int B;
    int bd;
    
  };


private:

  Site* site; // the list of the poiters to sites
  Interaction* interaction; // the list of the pointers to interactions
  Edge* edge; // the list of the poiters to sites

  XML::Block X;
  Algorithm& ALG;

public:

  int D; // dimension
  int BD; // bond dimension
  int* L; // linear size
  double BETA; // inverse temperature
  int NCELL; // total number of cells
  int NSITE; // total number of sites
  int NINT; // total number of interactions
  int NSTYPE; // number of site types
  int NITYPE; // number of interaction types
  int NEDGE;

  double **vec;

  Lattice(const char* FNAME, Algorithm& A);

  ~Lattice();

  void read();

  void initialize();

  Site& S(int i) { return site[i]; };

  Interaction& I(int i) { return interaction[i]; };

  Edge& EDGE(int i) { return edge[i]; }

  int countVertices();

  void show_param(FILE* F) {
    fprintf(F,"P D       = %12d\n",D);
    fprintf(F,"P L       =      ");
    for (int i=0; i<D; i++) {
      fprintf(F," %6d", L[i]);
    }
    fprintf(F,"\n");
    fprintf(F,"P BETA    = %24.16f\n", BETA);
  };

  void dump();

};

//######################################################################

inline Lattice::Lattice(const char* FNAME, Algorithm& A) : ALG(A) {

  if (DEBUG) printf("\nLattice::Lattice> Start.\n"); // koko

  X.initialize( FNAME , "LATTICE" );
  read();
  initialize();

  if (DEBUG) printf("Lattice::Lattice> End.\n"); // koko

}

//======================================================================

void Lattice::read() {
  
  if (DEBUG) printf("Lattice::read> Start.\n");

  D = X["Dimension"].getInteger();  
  //#ifdef WINDING
  BD = X["BondDimension"].getInteger();
  //#endif
  L = new int [D];
  for (int i=0; i<D; i++) {
    L[i] = X["LinearSize"].getInteger(i);
  }
  BETA = X["Beta"].getDouble();

  NCELL  = X["NumberOfCells"].getInteger();
  NSITE  = X["NumberOfSites"].getInteger();
  NINT   = X["NumberOfInteractions"].getInteger();
  NSTYPE = X["NumberOfSiteTypes"].getInteger();
  NITYPE = X["NumberOfInteractionTypes"].getInteger();
  //#ifdef WINDING
  NEDGE = X["NumberOfEdgeInteractions"].getInteger();
  //#endif
  
  site = new Site[NSITE];
  interaction = new Interaction[NINT];
  edge = new Edge [NEDGE];

  vec = new double* [D];
  for(int di=0;di<D;di++){
    vec[di] = new double [BD];
  }

  //  cout<<"Lattice::NumberOfBlocks"<<X.NumberOfBlocks()<<endl;

  for (int i=0; i<X.NumberOfBlocks(); i++) {
    XML::Block& B = X[i];

    if ( B.getName() == "S" ) {
      int id = B.getInteger(0);
      int st = B.getInteger(1);
      int mt = B.getInteger(2);
      
      if( B.NumberOfValues() == 3 ){
	//S(id).init( ALG.getSiteProperty(st) , mt );
	S(id).init( ALG.getSiteProperty(st) , mt, 0.0, 0.0 );
      }
      else{
	double cq = B.getDouble(3);
	double sq = B.getDouble(4);
	S(id).init( ALG.getSiteProperty(st) , mt, cq, sq );
      }

      
#ifdef NORM
      S(id).setBeta( 1.0 );
#else
      S(id).setBeta( BETA );
#endif
    }
    if ( B.getName() == "I" ) {
      int id = B.getInteger(0);
      int it = B.getInteger(1);
      int nb = B.getInteger(2);
      I(id).init(ALG.getInteractionProperty(it));
      for (int ii=0; ii<nb; ii++) {
        int sid = B.getInteger(3+ii);
        I(id).setSite(ii,S(sid));
      }
      //#ifdef WINDING
      int eid = B.getInteger(3+nb);
      int edim= B.getInteger(4+nb);
      if(eid >= 0 && nb==2){
	EDGE(eid).init(edim,I(id).site(0).id()-1,I(id).site(1).id()-1);
	//cout<<"EDGE:: eid="<<eid<<"  Asite="<<EDGE(eid).A<<"  Bsite="<<EDGE(eid).B<<"  D="<<EDGE(eid).bd<<endl;
      }
      //#endif
    }
    if ( B.getName() == "V" ) {

      int bd = B.getInteger(0);
      for(int di=0;di<D;di++){
	double length = B.getDouble(1+di);
	vec[di][bd]=length;
      }

    }

  }

  if(DEBUG) dump();
  if (DEBUG) printf("Lattice::read> End.\n");
}

//======================================================================
#include <set>
void Lattice::initialize() {

  if (DEBUG) printf("Lattice::initialize> Start.\n");

  //---// サイトに相互作用が働くサイトを登録する．//
  //set<Interaction*> InteractionOnEachSite[NSITE];
  set<Interaction*>* InteractionOnEachSite;
  InteractionOnEachSite = new set<Interaction*> [NSITE];
  
  int NRVIMAX = 0; // maximum number of registered vertex informations

  for(int iid=0; iid<NINT; iid++){
    for (int i=0; i<I(iid).NBODY(); i++) {

      InteractionOnEachSite[I(iid).site(i).id()-1].insert( &I(iid) );
    }
  }
 

  for(int sid=0; sid<NSITE; sid++){
    S( sid ).setNCI( InteractionOnEachSite[sid].size() );

    int counter = 0 ;
    set<Interaction*>::iterator it = InteractionOnEachSite[sid].begin();

    while( it != InteractionOnEachSite[sid].end() ){
      S( sid ).setCI( counter, (*it) );
      ++it;
      ++counter;
    }

    InteractionOnEachSite[sid].clear();
    if(NRVIMAX < S(sid).getNCI()) NRVIMAX = S(sid).getNCI();
  }

  TheRVIPool.init(NRVIMAX*NRVIMAX);
  delete [] InteractionOnEachSite;
  if (DEBUG) printf("Lattice::initialize> End.\n");
}

//======================================================================

inline Lattice::~Lattice() {
  // printf("*** Destroying Lattice\n");
  if ( L != 0 ) { delete [] L; }
  if ( site != 0 ) { delete [] site; }
  if ( interaction != 0 ) { delete [] interaction; }
  if ( edge != 0 ) { delete [] edge; }

  for(int i=0;i<D;i++){

    delete [] vec[i];
  }
  delete [] vec;

}

//======================================================================

inline int Lattice::countVertices() {
  int NV = 0;
  for (int b=0; b<NINT; b++) {
    NV += I(b).count();
  }
  return NV;
}

//======================================================================

inline void Lattice::dump() {
  printf("\n");
  printf("Lattice Information:\n");
  printf("  D      = %d\n",D);
  printf("  BD     = %d\n",BD);
  printf("  L      =");
  for (int i=0; i<D; i++) {
    printf(" %d", L[i]);
  }
  printf("\n");
  printf("  BETA   = %24.16f\n", BETA);
  printf("  NCELL  = %d\n",NCELL);
  printf("  NSITE  = %d\n",NSITE);
  printf("  NINT   = %d\n",NINT);
  printf("  NSTYPE = %d\n",NSTYPE);
  printf("  NITYPE = %d\n",NITYPE);
  printf("  NEDGE = %d\n",NEDGE);
  printf("\n");
  for (int i=0; i<NSITE; i++) { S(i).dump(); }
  printf("\n");
  for (int i=0; i<NINT; i++) { I(i).dump(); }
}

#endif
