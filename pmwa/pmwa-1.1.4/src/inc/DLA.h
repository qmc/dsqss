#ifndef DLA_H
#define DLA_H

#include <Configuration.h>
#include <Probability.h>

#define NUMPARA 25
//////////////////////////////////////////////////////////////////////////////////////////

class Dla{


  private:

  Size N;
  System sp;
  MC_p MC;

  Parallel PR;

  long long IMAX, WMAX;
  double allstart, allend;
  double pstart, pend, tstart, tend, ostart, oend;
  int cb;

  char confdir[128];
  char outfile[128];
  char sfinpfile[128];
  char sfoutfile[128];
  char algfile[128];
  char latfile[128];
  char Eventfile[128];
  char Eventfile_old[128];

  std::ofstream ftest;
  FILE* FOUT4SF; // file handler for the output file
  
  int num_para ;
  bool wcount[NUMPARA];
  string para[NUMPARA];

  void init();
  void set_value(int,int,double,string);

  void set_Parallel();

  int  numcheck(const string& );
  string trim(const string& );

  struct ToLower{char operator()(char c) { return tolower(c);}};


  void ReadBuffer( int p_num, int my_rank, int NP, char **PLIST );
  void NameOutfiles();


  public:

  double PMWA();

  Dla( int NP, char **PLIST );
  ~Dla();

 private:

  
  inline void show_MPIP(){

    cout<<"MPI check ---> rank= "<<PR.my_rank<<", V= "<<PR.V
	<<" B=  "<<PR.B<<", Nx= "<<PR.x<<", Ny= "<<PR.y<<", Nz= "<<PR.z
	<<",, nt= "<<PR.nt<<", nx= "<<PR.nx<<", ny= "<<PR.ny<<", nz= "<<PR.nz
	<<", up= "<<PR.upper<<", lower= "<<PR.lower
	<<", right(x)= "<<PR.right[0]<<", left(x)= "<<PR.left[0]
	<<", right(y)= "<<PR.right[1]<<", left(y)= "<<PR.left[1]<<endl;
  };

  inline void show_SP(){
    std::cout << "step =" << MC.Nstep <<" thermal="<< MC.Nthermal<< "  bin=" << MC.Nbin 
	      <<"  B="<< N.B <<",  Nx="<< N.x <<",  Ny="<< N.y <<"  Nz="<< N.z 
	      <<"  IMAX="<< IMAX <<"  WMAX="<< WMAX <<"  Ubb="<< sp.Ubb <<"  Vb1="<< sp.Vb1 
	      <<" Ntest ="<< MC.Ntest <<"  Nsample="<< MC.Nsample <<"  nc="<< MC.nc
	      <<"  tb="<< sp.tb <<"  mu="<< sp.mu <<" seed ="<< MC.seed <<" Vc="<< sp.Vc
	      <<"  nmax="<< sp.nmax << "  Nd=" << N.d << " Vb2=" << sp.Vb2 << "  Htr=" << sp.Htr 
	      << "  Ntvid=" << PR.Ntdiv << "  Npara=" << PR.Npara <<  "  Nsdiv=" << PR.Nsdiv 
	      << " outfile="<< outfile <<  std::endl;
  };

  inline void show_TIME(){

    ftest<<"Testing Time for Ncyc = "<< pend-pstart <<endl;
    ftest<<"Thermalization Time   = "<< tend-tstart <<endl;
    ftest<<"Measuring Time        = "<< oend-ostart <<endl;
    ftest<<"Total Simulation Time = "<< allend-allstart <<endl;

  }

};


#endif
