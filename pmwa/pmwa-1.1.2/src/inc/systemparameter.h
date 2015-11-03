#ifndef SYSTEMPARAM_H
#define SYSTEMPARAM_H

#define VAL 10
#define INF 1.0e+14
#define EPS 1.0e-14

struct Size{

  int x; //the number of lattice for x.
  int y; //the number of lattice for y.
  int z; //the number of lattice for z.
  int V; //the total number of lattice.
  int d;  //dimension.
  double B; //the inverse of temperature.
  double oldB; //for annealing.

};


struct System{

  int nmax; // the maximum number of bosons on a same site.
  double Vc; // external field.
  double Ubb; //on-site interaction.
  double Vb1; // nearest-neighbor interaction.
  double Vb2; // next nearest
  double tb; // hopping
  double mu; // chemical potential.
  int w_num; //
  double Htr; // the strenght of the transverse field for introducing worm.
  double Eu; // arbitrary energy shift for Ubb. 
  double Et; // arbitrary energy shift for t-Vb1. 
  double Ev2; // arbitrary energy shift for Vb2. 

  long seed;

};

struct MC_p{

  int seed;
  int Nstep;
  int Nthermal;
  int Nsample;
  int Nbin;
  int Ntest;
  int nc;
  int runtype;

};

struct Parallel{

  int p_num; // the total number of processors. 
  int my_rank; // the rank of each proseccor.
  int Ntdiv; // the number of the temporal decomposition.
  int Nsdiv, Nxdiv, Nydiv, Nzdiv; // the number of the spatial decomposition. (total, x, y, z)
  int NtNs; // the number of domain decomposition. 
  int Npara; //the number of trivial parallelization.
  int Rpara; //
  double B; // beta of "a domain".
  double oldB; // for annealing.

  int x,y,z,V; // the coordinate number.

  //please refer to "DLA_multi.cpp"
  int np;
  int nq, nr;
  int nt, ns, nx, ny, nz;
  int nst, nst0, nt0, ns0, nx0;

  int upper;
  int lower;
  int right[3];
  int left[3];

};

#endif
