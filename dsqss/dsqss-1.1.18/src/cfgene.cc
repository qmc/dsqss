/*---------------------------------------------

   Generating cf.xml for a hypercubic lattice

----------------------------------------------*/

#include<iostream>
#include<stdlib.h>
#include<fstream>
#include<math.h>

using namespace std;

//--------------------------------------------------------------
void ShowUsage(int argc,char** argv){
  cout<<"usage: $ "<<argv[0]<<" [ D,  L0 , L1 , Ntau]     \n";
  cout<<"    Hypercubic lattice L0 x L1 x L2...          \n";
  cout<<"                                         \n";
  cout<<"    L0, L1 ... the liner size of the lattice. \n";
  cout<<"          ( L0, L1, ..., must be even number. )\n";
}
//-------------------------------------------------------------
void WriteXML(int D, int L[], int Ntau) {
  
  ofstream fout("cf.xml");
  fout.precision(15);
  int N = 1; //number of sites.
  for(int i = 0 ; i<D ; i++) { N *= L[i] ; }

  int NumberOfElements = N*N;
  int NumberOfKinds    = N;
  
  fout<<"<CorrelationFunction>"<<endl<<endl;
  fout<<"<Comment>"<<endl;
  fout<<"  "<<D<<"-dimension hypercubic lattice"<<endl;
  fout<<"</Comment>"<<endl<<endl;
  
  fout<<"<Ntau>             " <<Ntau            <<" </Ntau>"<<endl;
  fout<<"<NumberOfElements> " <<NumberOfElements<<" </NumberOfElements>"<<endl;
  fout<<"<NumberOfKinds>    " <<NumberOfKinds   <<" </NumberOfKinds>"<<endl;
  fout<<endl;
  
  fout<<"<!-- <CF> [kind] [isite] [jsite] </CF> -->"<<endl<<endl;
 
  int  NB = 0; //3 * N ;   // number of bonds
  int* x  = new int[D];
  int* dx = new int[D];
  
  int kind = 0;
  
  for(int di=0; di<N; di++){
    int dk = di;
    for (int q=0; q<D; q++) {
      dx[q] =  dk % L[q];
      dk   /= L[q];
    }
    
    for (int i=0; i<N; i++) {
      int k = i;
      for (int q=0; q<D; q++) {
	x[q] = k % L[q];
	k /= L[q];
      }
      for (int q=0; q<D; q++) {
	x[q] = (x[q]+dx[q])%L[q];
      }
      int j = 0;
      for (int q=D-1; q>=0; q--) {
	j *= L[q];
	j += x[q];
      }
      
      fout<<"<CF> "<<kind<<" "<<i<<" "<<j<<" </CF>"<<endl;
    }
    kind++;
  }
  cout<<"Nkind = "<<kind<<endl;
  fout<<endl;
  fout<<"</CorrelationFunction>"<<endl;
  
  delete [] x;
  
}
//--------------------------------------------------------------

int main(int argc,char** argv){

  int NARG = 3;
  if ( argc < NARG + 1 ) {
    ShowUsage(argc,argv);
    exit(0);
  }
  int iarg = 1;
  const int  D    = atoi(argv[iarg]); iarg++;
  if ( argc != D + 3 ) {
    ShowUsage(argc,argv);
    exit(0);
  }

  int*       L =  new int [D] ;

  for(int i=0;i<D;i++){
    L[i] = atoi(argv[iarg]) ; iarg++;
  }
  int Ntau = atoi(argv[iarg]);
  
  int EvenOrOdd = 0;
  cout.precision(15);
  cout<<"D     = "<<D<<endl;
  for(int i = 0 ; i < D ; i++){
    cout<<"L"<<i<<"    = "<<L[i]<<endl;
    EvenOrOdd += L[i]%2 ;
  }
  
  if( EvenOrOdd ) { cout<<"Warnig: L should be an even number."<<endl;}

  WriteXML( D, L, Ntau);
  cout<<"... done."<<endl;
  delete [] L;
  return 0 ;
  
}
