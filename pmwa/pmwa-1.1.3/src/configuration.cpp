#include <Configuration.h>
#include <Probability.h>

//************************************Configurationクラスここから***********************************************************
//#define Bcheck
#include "graphspace.cpp"

Configuration::Configuration( MC_p *m_MC,Size *m_N, int m_nmax,Lattice *m_LT, Probability *m_P,Parallel *PR,  long long  m_IMAX,  long long  m_WMAX,char m_Eventfile[128], My_rdm *MR, bool cb, char m_outfile[128])
:GraphSpace(m_N,m_nmax,m_LT,m_P,PR,m_IMAX,m_WMAX,m_Eventfile,m_MC->nc,MR,cb,m_outfile){

  N=m_N;
  LT=m_LT;
  BV=PR->B*PR->V;
  p_num=PR->p_num;
  MC=m_MC;

  sprintf(ALGFILE,"algoorithm.xml");
  sprintf(LATFILE,"lattice.xml");


}
Configuration::~Configuration(){

}


void Configuration::DeterminationNworm( int MCS, My_rdm *MR, Quantities *QNT  ){

  //######## Ncyc is obtained here ########
  QNT->Init();
  QNT2=QNT;

  updateAnner(MCS,MR,QNT);

}

void Configuration::updateAnner( int MCS, My_rdm *MR, Quantities *QNT ){

  MPI_Status status;

  double ave_Wlen;
  int count,step=100;
  long pcyc;
  
  For_Nworm = true;
  MCS/=step;

  for(int i=0; i<step; i++){

    Wlen=0.0;
    for(int k=0; k<MCS; k++) if(Diagonal_Update(MR,QNT)) Worm_Update(MR);

    if(Wlen!=0.0){
      ave_Wlen = Wlen/(double)(Ncyc*MCS);
      Ncyc = (long)(BV/ave_Wlen);
      if(Ncyc==0) Ncyc=1;
    }
    else Ncyc=0;

    ////////////////////////////////////////////////////////////
    if( PR->my_rank==0 ){//Sum over t and s
      
      cout<<"rank"<<PR->my_rank<<"::Wnum="<<W[0]<<endl;
      count=(Ncyc)? 1:0;
      for(int tag=1; tag< PR->p_num; tag++ ){
	MPI_Recv(&pcyc, 1, MPI_LONG, tag, 0, MPI_COMM_WORLD,&status);
	if(pcyc>0){
	  count++;
	  cout<<"rank"<<tag<<":: Ncyc="<<pcyc<<endl;
	  Ncyc += pcyc;
	}
	else cout<<"rank::"<<tag<<" !!!!!No worm!!!!!"<<endl;
      }

      if( count > PR->p_num/2 ) Ncyc /= count;
      else{ cout<<"Larger eta is recomended."<<endl; exit(1);  }

    }
    else{
      cout<<"rank"<<PR->my_rank<<"::Wnum="<<W[0]<<endl;
      MPI_Send(&Ncyc, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD );
    }
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    MPI_Bcast(&Ncyc, 1, MPI_LONG, 0, MPI_COMM_WORLD ); 
    ////////////////////////////////////////////////////////////
  
    cout<< i <<": Wlen="<<Wlen/(double)MCS<<",  Ncyc="<<Ncyc<< ",  Wnum=" << W[0] << ",  EV=" << EVmax <<endl;
    
  }


}



void Configuration::sweep( int MCS, My_rdm *MR, Quantities *QNT ){

  //

  For_Nworm = false;
  Wlen=0.0;
  while(MCS--){

    if(Diagonal_Update(MR,QNT)) Worm_Update(MR);

  }
 
}

void Configuration::measurement( Quantities *QNT, My_rdm *MR ){

  QNT->Init();

  for(int bin=0; bin<MC->Nbin; bin++){

    sweep(MC->Nstep,MR,QNT);

    QNT->Init();

    for(int i=0; i<MC->Nsample; i++){
      
      For_Nworm = true;
      Wlen=0.0;
      if(Diagonal_Update(MR,QNT)) Worm_Update(MR);
      

      QNT->Measure( Nv, Nk, ev, WORM, world, worldB, Wlen, W[0], i);
    }
    

    QNT->Measure();

    if( PR->nst==0 )  QNT->BINsum(bin);

  }

    
  if( PR->nst==0 )   QNT->BINaverage();


}


void Configuration::dump(ofstream &F)
{
  F << "P NSET    = " << MC->Nbin <<endl;
  F << "P NTEST   = " << MC->Ntest <<endl;
  F << "P NMCSE   = " << MC->Nthermal <<endl;
  F << "P NMCSD   = " << MC->Nstep <<endl;
  F << "P NMCS    = " << MC->Nsample <<endl;
  F << "P SEED    = " << MC->seed <<endl;
  F << "P NVERMAX = " << IMAX <<endl;
  F << "P NWORMAX = " << WMAX <<endl;
  F << "P NCYC    = " << Ncyc <<endl;
  F << "P ALGFILE = " << ALGFILE <<endl;
  F << "P LATFILE = " << LATFILE <<endl;
  F << "P OUTFILE = " << outfile <<endl;
}
