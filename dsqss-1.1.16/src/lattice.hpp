
#ifndef LATTICE_H
#define LATTICE_H

//######################################################################

#include <stdio.h>
#include "io.h"
#include "objects.hpp"

//######################################################################

class Lattice {

private:

  Site* site; // the list of the poiters to sites
  Interaction* interaction; // the list of the pointers to interactions
  XML::Block X;
  Algorithm& ALG;

public:

  int D; // dimension
  int* L; // linear size
  double BETA; // inverse temperature
  int NCELL; // total number of cells
  int NSITE; // total number of sites
  int NINT; // total number of interactions
  int NSTYPE; // number of site types
  int NITYPE; // number of interaction types

  Lattice(const char* FNAME, Algorithm& A);

  ~Lattice();

  void read();

  void initialize();

  Site& S(int i) { return site[i]; };

  Interaction& I(int i) { return interaction[i]; };

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
  
  //  if (DEBUG) printf("Lattice::read> Start.\n");

  D = X["Dimension"].getInteger();  
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
  site = new Site[NSITE];
  interaction = new Interaction[NINT];
  for (int i=0; i<X.NumberOfBlocks(); i++) {
    XML::Block& B = X[i];

    if ( B.getName() == "S" ) {
      int id = B.getInteger(0);
      int st = B.getInteger(1);
      int mt = B.getInteger(2);
      S(id).init( ALG.getSiteProperty(st) , mt );

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
    }
  }
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
  printf("\n");
  for (int i=0; i<NSITE; i++) { S(i).dump(); }
  printf("\n");
  for (int i=0; i<NINT; i++) { I(i).dump(); }
}

#endif
