#ifndef QUANTITIES_H
#define QUANTITIES_H

#include<complex>

class Quantities{

  //  public:
 private:


  void WindingNumber( vector<GraphSpace::Vertex> &ev, int mcs);
  void WindingNumber( vector<GraphSpace::Vertex *> WORM, int mcs);
  void Density( GraphSpace::Vertex *world, GraphSpace::Vertex *worldB);
  void Energy(int countv, int mcs);
  void NumberOfVertices( int countv, int mcs  );
 void NumberOfKinks( int Nv, int mcs  );
  void NumberOfWorms( int x, int mcs );
  void CondensateFraction(int mcs, GraphSpace::Vertex *world);
  void CorrelationFunction1( );
  void CorrelationFunction2( GraphSpace::Vertex *world);
  void StructureFactor();
  void SpecificHeat();
  void Susceptibility();
  void NoiseCorrelation();
  void Compressibility();
  void WindingNumber2();

  int Nc, Nend, Lsize;
  int Cknum;
  int NVMAX, NWMAX;

  enum{ nver, nwor, nkin, ene, spe, xmx, amzu, bmzu, magx, magp, magm, Nq };
  //  enum{ nver, nwor, nkin, ene, spe, xmx, wnd2, wndx, wndy,wndz, amzu, bmzu, magx, magp, magm, Nq };

  inline int f_ld( int i ){ return i; }
  inline int f_gf( int i ){ return i+Lsum[0]; }
  inline int f_sk( int i ){ return i+Lsum[1]; }
  inline int f_nk( int i ){ return i+Lsum[2]; }
  inline int f_gf2( int i ){ return i+Lsum[3]; }
  inline int f_gk2( int i ){ return i+Lsum[4]; }
  inline int f_noise( int i ){ return i+Lsum[5];}


  inline int f_nkr( int k ){ 
    return k + Lsum[2]; 
  }
  inline int f_nki( int k ){ 
    return k + Lsum[2] + Nkmax; 
  }
  inline int f_gf2( int r1, int r2 ){ 
    return r1 + r2*V + Lsum[3]; 
  }
  inline int f_gk2r( int k, int k_ ){
    return k + k_*Nkxmax +Lsum[4]; 
  }
  inline int f_gk2i( int k, int k_ ){ 
    return k + k_*Nkxmax + Lsum[4] + Nkkmax*Nkxmax; 
  }

  inline int f_noise( int k, int k_ ){
    return k + k_*Nkxmax +Lsum[5]; 
  }
  
  inline int f_ck(int k, int num, int a){  return 5 + k + num*Nkmax + a*Nkmax*Cknum; };
  inline int f_ck(int k, int num){  return 5 + k + num*Nkmax; };

  //  inline int theta(int i, int k){ return i+(k+Nkxmax)*V; };
  inline int theta(int i, int k, int a){ return i+( k + (2*a + 1)*Nkmax)*V; }
 
  int *Lmax, *Lsum;

  void MCsum_S();
  double *values_S;
  double *MCmean_S;
  double *BINmean_S;
  double *RNDmean_S;

  void MCsum_L();
  double *values_k;
  double *values_L;
  double *MCmean_L;
  double *BINmean_L;
  double *RNDmean_L;


  double *m_val;
  double *AC;
  double *MCmeanACT;
  double *RNDmeanACT;

  char **Qname;
  char **file;
  char **acfile;

  void BINsum(double *MCmean, double *BINmean, int Lmax, int bin);
  void Average(double *g, int Nval, int S, double *MCmean, int kstep);

  void show_S(ofstream &F);
  void show_L();

  void RNDsum(double *local, double *global,int Ndiv, int my_rank, int rmax2, int Npara, int);

  void TreeSum( int num, double *child, double Normalization, int cnum, int mnum, int Nsum );
  void TreeSum( int num, double *child, double *mother, double Normalization, int cnum, int mnum, int Nsum );

  int V, Nx,  NUM;
  System *sp;
  Size *N;
  Parallel *PR;
  Lattice *LT;
  MC_p *MC;


  void autocorrelation( int mcs, double q ){  AC[mcs]=q; };
  void AutoCorrelationAverage( char *fname );


  //  private:

  int Npara;
  int my_rank;
  int Kcut, Kcut2;
  int Nxmax, Nkxmax, Nkmax, Nkkmax, Nkstep,  Nk_set;
  double *COSrk, *SINrk;
  complex <double> *EXPrk, *Ck, *Ck_m;

  double SFD_Norm;
  double cr0, an0;

  public:

  double *an, *cr, *Q, *ca, *ac;

  char parainfo[64];

  Quantities( Size *m_N, MC_p *m_MC, System *m_sp, Lattice *m_LT, Parallel *m_PR, char outfile[128] );
  ~Quantities();

  void Init();

  void Output(char fname[128], double g);

  void Measure( int Nv, int Nk, vector<GraphSpace::Vertex> &ev, vector<GraphSpace::Vertex *> WORM, GraphSpace::Vertex *world, GraphSpace::Vertex *worldB, double length, int m_Wnum, int mcs );

  void Measure();

  void BINsum(int);
  void BINaverage();

  void show(ofstream &F);


};


#endif
