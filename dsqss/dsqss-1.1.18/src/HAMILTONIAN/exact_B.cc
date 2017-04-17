#include <stdio.h>
#include <iostream>
#include <string>
using namespace std;
#include "matrix.h"
#include "boson_B.h"
#include "canonical.h"

//----------------------------------------------------------------------

// Bose-Hubbard Model
//H= - J\sum_{<i,j>} b_i b_j + V\sum_{<i,j>} n_i n_j
//    + u/2z \sum_{<i,j>}( n_i (n_i + 1) + n_j (n_j + 1) )
//     - mu/z \sum_{<i,j>}( n_i + n_j ).
//F=mu/z, 
//U=u/z

class BoseHubbardModel {
public:
  int M; // Number of bosons ( [M]-representaion )
  int NSITE;
  double J;    // the bilinear coupling for XY
  double V;    // the bilinear coupling for Z
  double U;    // the on-site coupling for Z
  double F;    // the field
  int DIM;     // dimension of the Hilbert space
  dgematrix H;
  dgematrix MXU;
  dgematrix MZU;
  dgematrix MXS;
  dgematrix MZS;
  dgematrix I;

  BoseHubbardModel( int M0, int NSITE0, double J0, double V0, double U0, double F0) {
    M = M0;
    NSITE = NSITE0;
    J = J0;
    V = V0;
    U = U0;
    F = F0;
    //double Vh=V*0.5;
    
    BosonOperatorSet S( M , NSITE );

    DIM = S.DIM;
    cmatrix h(DIM);
    h.zero();
    printf("  defining the hamiltonian...\n");
    for (int k=0; k<NSITE; k++) {
      int l =  (k+1) % NSITE;
      printf("    k= %d\n", k);
      //h += (-J) * ( S.X[k] * S.X[l] + S.Y[k] * S.Y[l] );
      h += (-J) * ( S.UP[k] * S.DN[l] );
      h += (+V) * ( S.Z[k] * S.Z[l] );

    }
    for (int k=0; k<NSITE; k++) {
      h += (+U) * ( S.Z[k] * S.Z[k] ); 
      h += (-F-U) * ( S.Z[k] ); 
    }
    printf("  ... done.\n");

    printf("  defining the magnetizations...\n");
    cmatrix mxu(DIM);
    cmatrix mzu(DIM);
    cmatrix mxs(DIM);
    cmatrix mzs(DIM);
    mxu.zero();
    mzu.zero();
    mxs.zero();
    mzs.zero();
    for (int k=0; k<NSITE; k++) {
      printf("   k= %d\n", k);
      double sgn = 1.0;
      if ( k % 2 == 1 ) sgn = -1.0;
      mxu +=  ( 1.0 * S.X[k] );
      mzu +=  ( 1.0 * S.Z[k] );
      mxs +=  ( sgn * S.X[k] );
      mzs +=  ( sgn * S.Z[k] );
    }
    printf("  ... done.\n");

    H   = h.re;
    MXU = mxu.re;
    MZU = mzu.re;
    MXS = mxs.re;
    MZS = mzs.re;
    I   = S.I.re;

  };
};

//============================================================================

void Average( 
	     int DIM , dgematrix& R , dgematrix& Q , double& Ave , double& Var ) {
  dgematrix W(DIM,DIM);
  W = Q * Q;
  Ave = CanonicalAverage( R , Q );
  Var = CanonicalAverage( R , W );
  Var = Var - Ave * Ave;
}

//============================================================================

void WriteXML(int M, dgematrix& Q, dgematrix& H) {

  FILE* FOUT = fopen("hamiltonian.xml","w");
  int D = M + 1;
  int DD = D * D ;
  fprintf(FOUT,"<Hamiltonian>\n");
  fprintf(FOUT,"  <General>\n");
  fprintf(FOUT,"    <Comment> Extended Bose-Hubbard model with NMAX=%d </Comment>\n", M);
  fprintf(FOUT,"    <NSTYPE> 1 </NSTYPE>\n");
  fprintf(FOUT,"    <NITYPE> 1 </NITYPE>\n");
  fprintf(FOUT,"    <NXMAX>  %d </NXMAX>\n", D );
  fprintf(FOUT,"  </General>\n");
  fprintf(FOUT,"\n");
  fprintf(FOUT,"  <Site>\n");
  fprintf(FOUT,"    <STYPE> 0 </STYPE>\n");
  fprintf(FOUT,"    <TTYPE> 0 </TTYPE>\n");
  fprintf(FOUT,"    <NX>   %d </NX>\n", D);
  fprintf(FOUT,"  </Site>\n");
  fprintf(FOUT,"\n");
  fprintf(FOUT,"  <Source>\n");
  fprintf(FOUT,"    <TTYPE> 0 </TTYPE>\n");
  fprintf(FOUT,"    <STYPE> 0 </STYPE>\n");
  for (int i=0; i<D; i++) {
    for (int j=0; j<D; j++) {
      double x = Q(i,j);
      if ( abs( x ) > 1.0e-8 ) {
        fprintf(FOUT, "    <Weight> %d %d %24.16f </Weight>\n", i, j, x );
      }
    }
  }
  fprintf(FOUT,"  </Source>\n");
  fprintf(FOUT,"\n");
  fprintf(FOUT,"  <Interaction>\n");
  fprintf(FOUT,"    <ITYPE> 0 </ITYPE>\n");
  fprintf(FOUT,"    <NBODY> 2 </NBODY>\n");
  fprintf(FOUT,"    <STYPE> 0 0 </STYPE>\n");
  for (int i=0; i<DD; i++) {
    int i0 = i % D;
    int i1 = i / D;
    for (int j=0; j<DD; j++) {
      int j0 = j % D;
      int j1 = j / D;
      double x = -H(i,j);
      if ( abs( x ) > 1.0e-8 ) {
        if ( i != j ) x = abs(x);
        fprintf(FOUT, "    <Weight> %d %d %d %d %24.16f </Weight>\n", 
		i0, j0, i1, j1, x );
      }
    }
  }
  fprintf(FOUT,"  </Interaction>\n");
  fprintf(FOUT,"</Hamiltonian>\n");
}

//============================================================================
//    Main
//============================================================================

int main(int argc, char** argv) {

  //  using namespace CPPL;

  int NARG;

  string name = argv[0];

//  if ( name == "./hamgen_H" ) {

    NARG = 5;
    if ( argc != NARG+1 ) {
      printf("usage: $ %s [ M , J , F ]\n", argv[0]);
      printf("    M ... the maxmum number of bosons on each site \n");
      printf("    J ... the hopping constant (positive)\n");
      printf("    V ... the nearest neighbor interaction (positive for replusive)\n");
      printf("    U ... the on-site interaction (positive for replusive)\n");
      printf("          ( = u/z if the field u is the field per site.)\n");
      printf("    F ... the chemical potential in the pair Hamiltonian\n");
      printf("          ( = H/z if the field H is the field per site\n");
      printf("            and H is shared equally by all pairs,\n");
      printf("            e.g., F = H/4 for a square lattice. )\n");
      exit(0);
    }
    int    M  = atoi(argv[1]);
    double J  = (double)atof(argv[2]);
    double V  = (double)atof(argv[3]);
    double U  = (double)atof(argv[4]);
    double F  = (double)atof(argv[5]);
    printf(" M     = %4d\n", M );
    printf(" J     = %8.3f\n", J );
    printf(" V     = %8.3f\n", V );
    printf(" U     = %8.3f\n", U );
    printf(" F     = %8.3f\n", F );
    BosonOperator S( M );
    BoseHubbardModel MDL( M , 2 , J, 0.5*V, 0.5*U, F );
    WriteXML( M, S.X.re, MDL.H);
    exit(0);

//  }

// -- comment out -- // edit sakakura
//NARG = 5;
//if ( argc != NARG+1 ) {
//  printf("usage: \n");
//  printf("  $ %s [ M, NSITE , J , F , B ]\n", argv[0]);
//  printf("    M ... the number of bosons on each site \n");
//  printf("            ( M= 1, 2, 3, ... for S= 1/2, 1, 3/2, ... ) \n");
//  printf("    NSITE ... the number of sites ( the length of the ring ) \n");
//  printf("    J ... the coupling constant (positive for ferromagnets)\n");
//  printf("    F ... the magnetic field in the pair Hamiltonian\n");
//  printf("          ( = H/z if the field H is shared equally by all\n");
//  printf("            pairs, where z = 2 for a ring geometry. )\n");
//  printf("    B ... the inverse temperature\n");
//  printf("\n");
//  printf("    CAUTION: NSITE=2 means a ring with length 2,\n");
//  printf("             which is a 2-site system of \n");
//  printf("             *DOUBLY* coupled spins. \n");
//  exit(0);
//}

//int    M  = atoi(argv[1]);
//int    NSITE = atoi(argv[2]);
//double J  = (double)atof(argv[3]);
//double F  = (double)atof(argv[4]);
//double B  = (double)atof(argv[5]);

//printf(" M     = %4d\n", M );
//printf(" NSITE = %4d\n", NSITE );
//printf(" J     = %8.3f\n", J );
//printf(" F     = %8.3f\n", F );
//printf(" B     = %8.3f\n", B  );

//BoseHubbardModel Model( M , NSITE , J , F );

//int DIM = Model.DIM;

//vector<double> V(DIM);
//dgematrix U(DIM,DIM);

//printf("  diagonalizing the hamiltonian...\n");
//diagonalize( Model.H, V, U);

//dump(Model.H);

//double emin = V[0];
//for (int i=0; i<DIM; i++) if ( emin > V[i] ) emin = V[i];
//emin /= (double)NSITE;
//
//printf("  ... done. (emin = %16.6f)\n", emin);

//dgematrix R(DIM,DIM);

//R = DensityMatrix( B, V, U);

////  dump("R=",R);

//double AVE, VAR;
//double aen, spe;
//double amxu, smxu, xmxu;
//double amzu, smzu, xmzu;
//double amxs, smxs, xmxs;
//double amzs, smzs, xmzs;

//Average( DIM , R , Model.H , AVE , VAR );
//aen = AVE / (double)Model.NSITE;
//spe = VAR * B * B / (double)Model.NSITE;

//Average( DIM , R , Model.MXU , AVE , VAR );
//amxu = AVE / (double)Model.NSITE;
//smxu = VAR / (double)Model.NSITE;
//xmxu = Susceptibility( B , V , U , Model.MXU ) / (double)Model.NSITE;

//Average( DIM , R , Model.MZU , AVE , VAR );
//amzu = AVE / (double)Model.NSITE;
//smzu = VAR / (double)Model.NSITE;
//xmzu = Susceptibility( B , V , U , Model.MZU ) / (double)Model.NSITE;

//Average( DIM , R , Model.MXS , AVE , VAR );
//amxs = AVE / (double)Model.NSITE;
//smxs = VAR / (double)Model.NSITE;
//xmxs = Susceptibility( B , V , U , Model.MXS ) / (double)Model.NSITE;

//Average( DIM , R , Model.MZS , AVE , VAR );
//amzs = AVE / (double)Model.NSITE;
//smzs = VAR / (double)Model.NSITE;
//xmzs = Susceptibility( B , V , U , Model.MZS ) / (double)Model.NSITE;

//printf("\n");
//printf("================ RESULT ===============\n");
//printf("\n");
//printf("M     = %4d\n", M );
//printf("NSITE = %4d\n", NSITE );
//printf("J     = %8.3f\n", J );
//printf("F     = %8.3f\n", F );
//printf("B     = %8.3f\n", B  );
//printf("\n");
//printf("emin  = %16.6f\n",emin);
//printf("aen   = %16.6f\n",aen);
//printf("spe   = %16.6f\n",spe);
//printf("\n");
//printf("amxu  = %16.6f\n",amxu);
//printf("smxu  = %16.6f\n",smxu);
//printf("xmxu  = %16.6f\n",xmxu);
//printf("\n");
//printf("amzu  = %16.6f\n",amzu);
//printf("smzu  = %16.6f\n",smzu);
//printf("xmzu  = %16.6f\n",xmzu);
//printf("\n");
//printf("amxs  = %16.6f\n",amxs);
//printf("smxs  = %16.6f\n",smxs);
//printf("xmxs  = %16.6f\n",xmxs);
//printf("\n");
//printf("amzs  = %16.6f\n",amzs);
//printf("smzs  = %16.6f\n",smzs);
//printf("xmzs  = %16.6f\n",xmzs);

//FILE* FOUT = fopen("res_exact.dat","w");

//fprintf(FOUT, " M     = %4d\n", M );
//fprintf(FOUT, " NSITE = %4d\n", NSITE );
//fprintf(FOUT, " J     = %8.3f\n", J );
//fprintf(FOUT, " F     = %8.3f\n", F );
//fprintf(FOUT, " B     = %8.3f\n", B  );
//fprintf(FOUT, "emin   = %16.6f\n",emin);
//fprintf(FOUT, "aen    = %16.6f\n",aen);
//fprintf(FOUT, "spe    = %16.6f\n",spe);
//fprintf(FOUT,"\n");
//fprintf(FOUT, "amxu  = %16.6f\n",amxu);
//fprintf(FOUT, "smxu  = %16.6f\n",smxu);
//fprintf(FOUT, "xmxu  = %16.6f\n",xmxu);
//fprintf(FOUT,"\n");
//fprintf(FOUT, "amzu  = %16.6f\n",amzu);
//fprintf(FOUT, "smzu  = %16.6f\n",smzu);
//fprintf(FOUT, "xmzu  = %16.6f\n",xmzu);
//fprintf(FOUT,"\n");
//fprintf(FOUT, "amxs  = %16.6f\n",amxs);
//fprintf(FOUT, "smxs  = %16.6f\n",smxs);
//fprintf(FOUT, "xmxs  = %16.6f\n",xmxs);
//fprintf(FOUT,"\n");
//fprintf(FOUT, "amzs  = %16.6f\n",amzs);
//fprintf(FOUT, "smzs  = %16.6f\n",smzs);
//fprintf(FOUT, "xmzs  = %16.6f\n",xmzs);

// -- comment out -- // edit sakakura
  return 0;
}

