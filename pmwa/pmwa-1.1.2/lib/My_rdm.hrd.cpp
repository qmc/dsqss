#include <My_rdm.hrd.h>

My_rdm::My_rdm(long seed){

  RND.setSeed(seed,32);

}

My_rdm::~My_rdm(){


}

void My_rdm::outgen(char *fname){

  ofstream fout;
  fout.open(fname,ios::out|ios::binary);

  fout.write((char *) this,sizeof(My_rdm));
  fout.close();

}

void My_rdm::ingen(char *fname){

  ifstream fin;
  fin.open(fname,ios::in|ios::binary);

  if(fin) fin.read((char *) this,sizeof(My_rdm));
  else{cout<<"no file!"<<endl;};
  fin.close();

}
