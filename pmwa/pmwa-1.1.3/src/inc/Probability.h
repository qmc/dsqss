#ifndef PROB_H
#define PROB_H

#include <simtype.h>
#include <systemparameter.h>
#include <stdma.h>

class Probability{

  int nmax, XMAX;
  double Ubb, V1, V2, z;


  public:

  double tb, local_Et, rtmax, rv2max, dim, local_Ev2, *rumax, *local_Eu;
  double rh_odd, rh_even;

  Probability( Size *N, System *sp, Parallel *m_PR );
  ~Probability();
  void LookUpTable( Size *N, System *sp );
  
  double ******t;
  double ***u;
  double **v2;
  double **w;
  double ***at;
  double **av2;
  double **ru;

//###################################################
  double ***FracW;
//###################################################

  class Omega{

    public:
    double val; 
    int num;
    bool off;
    int i;
  };


  private:

  Parallel *PR;
  int Nx;

  double *ex_Wall;
  double *ex_Penki;
  double **Wall;

  int *Tr;
  Omega *Om;


  double mu1[2], mu2[2];
  double *ep;
  double *weight;

  void LocalPotential( double m_Vc, double mu, long seed );

  double at_make(int p, int q, int x);
  double au_make(int p, int x);
  double av2_make(int p, int q);

  void look( Size *N, System *sp );


  void Color(int);
  void SolveWeightEquation(int cmax);

  double Tuab(int i, int j, int x);
  double Tv2ab(int i, int j);

};



#endif
