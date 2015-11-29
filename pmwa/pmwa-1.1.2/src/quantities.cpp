//#include<stdma.h>
#include<Configuration.h>
//¤³¤³¤Ç¤Îpnum¤ÏNpara¤Î¤³¤È
//#define ATtime
#define TYPE_A
//#define TYPE_B
//#define STEPALL
//#define CFOUT 
//#include <cmath>


Quantities::Quantities( Size *m_N, MC_p *m_MC, System *m_sp, Lattice *m_LT, Parallel *m_PR, char outfile[128] ){

  Npara=m_PR->Npara;
  my_rank=m_PR->my_rank;

  MC=m_MC;
  N=m_N;
  LT=m_LT;
  PR=m_PR;

  sp=m_sp;
  V=PR->V;
  Nx=PR->x;

  SFD_Norm=(double)N->V*N->B*sp->tb*2.0*N->d;

  NVMAX=0;
  NWMAX=0;


#ifdef STEPALL
  Nkstep=1;
#else
  Nkstep=Nx/2;
#endif


  Nkxmax=N->x/Nkstep;
  Nxmax=Nkxmax;
  Nkmax=Nkxmax*2; //i.e. kx=ky=kz, kx!=ky=kz=0
  Nkkmax=Nkxmax;

  Nc=7; //the number of multi-point functions.
  Nq=12; //the number of one-point quantities +2.
  Nend=Nq-3; //Mx's number 


  newcall(Lmax,Nc);
  newcall(Lsum,Nc);

  Lsum[0]=Lmax[0]=V;
  Lmax[1]=Nxmax;
  Lmax[2]=Nkmax*2;//i.e. real and imaginary part
  Lmax[3]=Nkmax*2;
#ifdef GF22Gk2
  Lmax[4]=N->x*N->x;
#else 
  Lmax[4]=Nxmax*PR->Nxdiv;
#endif
  Lmax[5]=Nkxmax*Nkkmax*2;
  Lmax[6]=Nkxmax*Nkkmax;

  Lsize=0; //sin,cos,Nk,Sk

  for(int i=1; i<Nc; i++)Lsum[i]=Lsum[i-1]+Lmax[i];
  for(int i=0; i<Nc; i++)Lsize+=Lmax[i];

  newcall(file,Nc+Nq,128);
  newcall(Qname,Nc+Nq,128);


  nver=0;
  nwor=1;
  nkin=2;
  ene=3;
  spe=4;
  susx=5;
  wind=6;
  comp=7;
  magz=8;
  magx=Nend;
  magp=Nend+1;
  magm=Nend+2;

  int NVMAX, NWMAX; 

  //  sprintf(parainfo,"B%.1lf_Nx%d_Ny%d_Vbb%.1lf",N->B,N->x,N->y,sp->Vb1);
  sprintf(Qname[wind],"wind "); //Winding Number
  sprintf(Qname[magz],"amzu "); //Density or z-Magnetization
  sprintf(Qname[ene] ,"ene "); //Energy
  sprintf(Qname[nver],"nver "); //Number of vertices
  sprintf(Qname[spe] ,"spe "); //Specific heat
  sprintf(Qname[nwor],"nworm "); //Number of worms
  sprintf(Qname[susx],"xmx "); //Susceptibility
  sprintf(Qname[nkin],"nkink ");  //Number of kinks
  sprintf(Qname[comp],"comp ");  //Compressibility
  sprintf(Qname[magx],"amxu ");  //BEC order parameter or x-Magnet
  sprintf(Qname[magp],"ampu ");  //BEC order parameter or x-Magnet
  sprintf(Qname[magm],"ammu ");  //BEC order parameter or x-Magnet

  sprintf(Qname[Nq],"amz");  //Local density or local mz 
  sprintf(Qname[Nq+1],"cr2"); //2-points correlation (mx-mx)
  sprintf(Qname[Nq+2],"szk");  //structure factor
  sprintf(Qname[Nq+3],"ck2"); //k-space of cr2
  sprintf(Qname[Nq+4],"cr4"); //4-points correlation 
  sprintf(Qname[Nq+5],"ck4"); //k-space of cr4
  sprintf(Qname[Nq+6],"dkk"); //Noise correlation

  int S=Nc+Nq;
  /*
  if(sp->M_Anneal==1){
    cout<<"M_Anneal ON"<<endl;
    sprintf(outfile,"%s.m%.4lf",outfile,sp->mu);
  }
  else cout<<"M_Anneal OFF"<<endl;
  */
  for(int name=0; name<S;name++){
    sprintf(file[name],"%s_%s",Qname[name],outfile);
  }

  
  newcall(values_S,Nq);
  newcall(MCmean_S,Nq*2);
  newcall(BINmean_S,Nq*2*MC->Nbin);
  newcall(RNDmean_S,Nq*2*PR->Npara);

  newcall(values_L,Lsize);
  newcall(MCmean_L,Lsize*2);
  newcall(BINmean_L,Lsize*MC->Nbin);
  newcall(RNDmean_L,Lsize*PR->Npara);

  
  newcall(m_val,max(max(Lsize,Nq),max(Lmax[1]*3,V*2)));

  newcalls(EXPrk,4*Nkmax*V);

  Cknum=16;
  Nk_set=5+2*Nkmax*Cknum;

  newcalls(Ck,Nk_set);
  newcalls(Ck_m,Nk_set);


  double PI=M_PI;
  for(int i=0;i<V;i++){
    int x=(i%Nx) + PR->nx*Nx;
    int y=(int)(i/Nx)%PR->y + PR->ny*PR->y;
    int z=(int)(i/(Nx*PR->y)) + PR->nz*PR->z;

    for(int k=-Nkmax+1;k<Nkmax;k++){ //max(kx+kk)=Nkkmax-1 + Nkxmax-1
   
      double phase = 2.0*PI*x*k/(double)Nkxmax;//kx!=0,ky=kz=0
      complex<double> phase_c(0.0,phase);
      EXPrk[theta(i,k,0)]=exp(phase_c);

      phase = 2.0*PI*(x + y + z)*k/(double)Nkxmax;//kx=ky
      phase_c=complex< double >(0.0,phase);
      EXPrk[theta(i,k,1)]=exp(phase_c);
    }

  }

#ifdef ATtime
  newcall(acfile,Nq,128);
  sprintf(acfile,fname);
  newcall(AC,MC->Nsample);
  newcall(MCmeanACT,2);
  newcall(RNDmeanACT,2*PR->Npara);
#endif

  newcall(an,V);
  newcall(cr,V);
  newcall(ca,V);
  newcall(ac,V);
  newcall(Q,Lmax[1]*3);

}


Quantities::~Quantities(){

  delcall(Lmax);  
  delcall(Lsum);  

  delcall(file,Nc+Nq);
  delcall(Qname,Nc+Nq);

  delcall(values_S);  
  delcall(MCmean_S);
  delcall(BINmean_S);
  delcall(RNDmean_S);

  delcall(values_L);
  delcall(MCmean_L);
  delcall(BINmean_L);
  delcall(RNDmean_L);

  delcall(m_val);

  delcall(EXPrk);

  delcall(Ck);
  delcall(Ck_m);

#ifdef ATtime
  delcall(acfile,Nq);
  delcall(AC);
  delcall(MCmeanACT);
  delcall(RNDmeanACT);
#endif

  delcall(an);
  delcall(cr);
  delcall(ca);
  delcall(ac);
  delcall(Q);

}


void Quantities::Init(){

  for( int i=0; i<Nq; i++ ) values_S[i]=0;
  for( int i=0; i<Lsize; i++ ) values_L[i]=0;
  for( int i=0; i<Nq*2; i++ ) MCmean_S[i]=0;
  for( int i=0; i<Lsize*2; i++ ) MCmean_L[i]=0;
  //  for( int i=0; i<Nq*2*PR->Npara; i++ ) RNDmean_S[i]=0;
  //  for( int i=0; i<Lsize*PR->Npara; i++ ) RNDmean_L[i]=0;

}


void Quantities::Measure( int Nv, int Nk, vector<GraphSpace::Vertex> &ev, vector<GraphSpace::Vertex *> WORM, GraphSpace::Vertex *world, GraphSpace::Vertex *worldB, double length, int m_Wnum, int mcs ){

  NVMAX=max(NVMAX, m_Wnum+Nv);
  NWMAX=max(NWMAX, m_Wnum);

#ifdef TYPE_A
    WindingNumber(ev,mcs);
#elif TYPE_B
    WindingNumber(WORM,mcs);
#endif
    //    cout<<"Eng"<<endl;
    Energy(m_Wnum+Nv,mcs);
    // cout<<"Nv"<<endl;
    NumberOfVertices(m_Wnum+Nv,mcs);
    //  cout<<"Nw"<<endl;
    NumberOfWorms(m_Wnum,mcs);
   // cout<<"Nk"<<endl;
    NumberOfKinks(Nk,mcs);
    // cout<<"BEC"<<endl;
    CondensateFraction(mcs,world);
    //  cout<<"Dns"<<endl;
    Density(world,worldB);
    // cout<<"Sk"<<endl;
#ifdef CFOUT
    StructureFactor();
#endif
    ////////////////////////
    MCsum_S();
#ifdef CFOUT
    MCsum_L();
#endif
    ////////////////////////

}

void Quantities::Measure(){

  SpecificHeat();
  Susceptibility();
  Compressibility();
#ifdef CFOUT
  NoiseCorrelation( );
#endif

}


void Quantities::RNDsum(double *MCmean, double *RNDmean, int NtNs, int my_rank, int Lmax, int Npara, int np){

  MPI_Status status;

  if (my_rank==0){

      for(int i=0; i< Lmax; i++ ) RNDmean[i]=MCmean[i];

      for(int tag=1; tag< Npara; tag++ ){

        MPI_Recv(&(RNDmean[tag*Lmax]), Lmax, MPI_DOUBLE, tag*NtNs, 0, MPI_COMM_WORLD,&status);

     }
  }
  else{
      MPI_Send(MCmean, Lmax, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
  }


}

void Quantities::BINsum(double *MCmean, double *BINmean, int Lmax, int bin){


  for(int i=0; i< Lmax; i++ ) BINmean[i+bin*Lmax]=MCmean[i];

}

void Quantities::Average(double *g, int Nval, int S, double *MCmean, int kstep){

  double num =(Nval==1)? 1.0:  Nval-1.0;

  for(int k=0; k<S; k+=2){
    double mean=0.0, error=0.0;

    for(int i=0; i<Nval; i++){
   
      mean  += g[(k + i*S)/kstep];
      error += g[(k + i*S)/kstep]*g[(k + i*S)/kstep];
    }

    //    if(isnan(mean)){cout<<"k="<<k<<"   mean="<<mean<<endl; exit(1);}

    mean  /= (double)Nval;
    error /= (double)Nval;

    MCmean[k] = mean;
    MCmean[k+1] = sqrt( ( error - mean*mean )/num );

  }
  
}

void Quantities::show_S(ofstream &F){//ÄÌ¾ï

  int S=Nq*2;

  for(int k=0; k<S; k+=2) F <<"R " <<Qname[k/2]<<" = "<< MCmean_S[k] <<" "<< MCmean_S[k+1]<<endl;

  //  F <<"R "<<Qname[Nend]<<" = ";
  //  for(int k=0; k<6; k++) F <<" "<< MCmean_S[k+S];   
  //  F << endl;

}

void Quantities::show_L(){//ÄÌ¾ï

  int k, S=Nq;
  long R, R_, Ry;
  string ll, Cname;

  std::ofstream F; 

  //Local density
  F.open( file[S] , std::ios::app );
  F <<"# Density at each site." <<endl;
  F <<"# R x-SITE  <VAL>  <ERROR>" <<endl;
  F << setprecision(16);
  for(int l=0; l<Nx; l++){

    R=l;   
    k=l*2;
    //    F<<R<<" "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 
    Cname = "real_" + tostr(R);
    F<<"R "<< Cname <<" = "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 
  }
  F.close();

  //GF
  F.open( file[S+1] , std::ios::app );
  F <<"#2-point Correlation in the real space." <<endl;
  F <<"# R x-distance  <VAL>  <ERROR>" <<endl;
  F << setprecision(16);
  for(int l=0; l<Lmax[1]; l++){
      
    R =l*Nkstep;   
    k=(l+Lsum[0])*2;
    //    F<<R<<" "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 
    Cname = "real_" + tostr(R);
    F<<"R "<< Cname <<" = "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 
  }
  F.close();

  //structure factor
  F.open( file[S+2] , std::ios::app );
  F <<"# Structure factor." <<endl;
  F <<"# real/imag  kx  ky  <VAL>  <ERROR>" <<endl;
  F << setprecision(16);
  for(int l=0; l<Lmax[2]; l++){
      
    R=l%Nkxmax;
    Ry=((int)(l/Nkxmax)%2)*R;
    ll =((bool)(l/Nkmax))? "imag" : "real";
    k=(l+Lsum[1])*2;
    //    F<<ll<<" "<<R*Nkstep<<" "<<Ry*Nkstep<<" "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 
    Cname = ll + "_" + tostr(R*Nkstep) + "_" + tostr(Ry*Nkstep);
    F<<"R "<< Cname <<" = "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 

  }
  F.close();

  //nk
  F.open( file[S+3] , std::ios::app );
  F <<"# 2-point Correlation of the wave number space." <<endl;
  F <<"# real/imag  kx  ky  <VAL>  <ERROR>" <<endl;
  F << setprecision(16);
  for(int l=0; l<Lmax[3]; l++){
      
    R=l%Nkxmax;
    Ry=((int)(l/Nkxmax)%2)*R;
    ll =((bool)(l/Nkmax))? "imag" : "real";
    k=(l+Lsum[2])*2;
    //    F<<ll<<" "<<R*Nkstep<<" "<<Ry*Nkstep<<" "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 
    Cname = ll + "_" + tostr(R*Nkstep) + "_" + tostr(Ry*Nkstep);
    F<<"R "<< Cname <<" = "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 

  }
  F.close();
  /*  // 2body G
  F.open( file[S+4] , std::ios::app );
  F <<"# 4-point Correlation of the real space." <<endl;
  F <<"# Rx  Ry  <VAL>  <ERROR>" <<endl;
  F << setprecision(16);
  for(int l=0; l<Lmax[4]; l++){
    
    //    F<<sp->Htr<<" "<<sp->mu<<" "<<sp->Vb1<<" "<<sp->tb<<" "<<N->B<<" ";
    k=(l+Lsum[3])*2;
#ifdef GF22Gk2
    R=l%Nxmax;
    Ry=(int)(l/Nxmax);
    F<<R<<" "<<Ry<<" "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 
#else
    R =l*Nkstep;   
    F<<R<<" "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 
#endif
  }
  F.close();*/
  //4-point correlation
  F.open( file[S+5] , std::ios::app );
  F <<"# 4-point Correlation of the wave number space." <<endl;
  F <<"# kx  kx'  <VAL>  <ERROR>" <<endl;
  F << setprecision(16);
  for(int l=0; l<Lmax[5]; l++){
    
    //    F<<sp->Htr<<" "<<sp->mu<<" "<<sp->Vb1<<" "<<sp->tb<<" "<<N->B<<" ";
    k=(l+Lsum[4])*2;
    R=l%Nkxmax;
    R_=(int)(l/Nkxmax)%Nkkmax;
    bool kl=l/(Nkkmax*Nkxmax);
    ll=(kl)? "imag" : "real";// real or imag
    //    F<<ll<<" "<<R*Nkstep<<" "<<R_*Nkstep<<" "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; ;
    Cname = ll + "_" + tostr(R*Nkstep) + "_" + tostr(R_*Nkstep);
    F<<"R "<< Cname <<" = "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 

    if(R==Nkxmax-1) F<<endl;
  }
  F.close();

  //noise correlation
  F.open( file[S+6] , std::ios::app );
  F <<"# Noise Correlation with wave numbers." <<endl;
  F <<"# kx  kx'  <VAL>  <ERROR>" <<endl;
  F << setprecision(16);
  for(int l=0; l<Lmax[5]; l++){
    
    k=(l+Lsum[5])*2;
    R=l%Nkxmax;
    R_=(int)(l/Nkxmax);
    ll ="real";
    Cname = ll + "_" + tostr(R*Nkstep) + "_" + tostr(R_*Nkstep);
    F<<"R "<< Cname <<" = "<<MCmean_L[k]<<" "<<MCmean_L[k+1]<<endl; 

    if(R==Nkxmax-1) F<<endl;
  }
  F.close();

    
}

///////
void Quantities::BINsum(int bin){

  BINsum(MCmean_S, BINmean_S, Nq*2, bin);
#ifdef CFOUT 
  BINsum(MCmean_L, BINmean_L, Lsize, bin);
#endif
}

void Quantities::BINaverage(){

  Average(BINmean_S, MC->Nbin, Nq*2, MCmean_S, 1);
#ifdef CFOUT 
  Average(BINmean_L, MC->Nbin, Lsize*2, MCmean_L, 2 );
#endif
}

void Quantities::show(ofstream &F){//ÄÌ¾ï

  F<<"I [the maximum number of vertices]          ="<<NVMAX<<endl;
  F<<"I [the maximum number of worms]             ="<<NWMAX<<endl;
  show_S(F);
#ifdef CFOUT 
  show_L();
#endif

}
//////
void Quantities::Output(char fname[128], double g){//¥¨¥é¡¼¥Ð¡¼¤Ê¤·

  std::ofstream file( fname, std::ios::app );
  file<<sp->Htr<<" "<<sp->mu<<" "<<sp->Vb1<<" "<<sp->tb<<" "<<N->B<<" "<<g<<endl;
  file.close();

}


//##########################################################################################

//*******************¥ï¥¤¥ó¥Ç¥£¥ó¥°¥Ê¥ó¥Ð¡¼************************

void Quantities::WindingNumber(vector<GraphSpace::Vertex> &ev, int mcs){

  MPI_Status status;
  int Wi[3]={0,0,0};

  vector<GraphSpace::Vertex>::iterator it = ev.begin();
  while( it != ev.end() ){

      int l = it->i/V;

      if(( it->type==2 && /*LT->bd[it->i%V][l*2] == it->nleg->i%V &&*/ it->i%V < it->nleg->i%V  ) || it->type==-1){	

	int crr = ( it->p < it->next[0]->p )?  1 : -1;
	for(int d=0; d<N->d; d++) if(LT->bond_vec[l][d]!=0.0) Wi[d] += crr;
	
      }

      ++it;
  }


  if( PR->nst==0 ){//Sum over t and s

    int m_Wi[3];

    for(int tag=1; tag< PR->NtNs; tag++ ){
 
      MPI_Recv(m_Wi, 3, MPI_INT, tag + PR->nst0, 0, MPI_COMM_WORLD,&status);
      for(int i=0; i<3; i++) Wi[i] += m_Wi[i];
    }

    values_S[wind] = ( Wi[0]*Wi[0] + Wi[1]*Wi[1] + Wi[2]*Wi[2]); ///SFD_Norm;

  }
  else{

    MPI_Send(Wi, 3, MPI_INT, PR->nst0, 0, MPI_COMM_WORLD );
    values_S[wind] = 0.0; //Winding
  }
 

#ifdef ATtime
  autocorrelation(mcs,values_S[0]);
#endif

}

void Quantities::WindingNumber(vector<GraphSpace::Vertex *> WORM, int mcs){

  MPI_Status status;
  int Wi[3]={0,0,0};
  bool DIR;
  int PAIR;

  vector<GraphSpace::Vertex *>::iterator it = WORM.begin();

  while( it != WORM.end() ){ (*it)->dir=1; ++it;}

  it=WORM.begin();

  while( it != WORM.end() ){

    if((*it)->dir){ //creation worm
      
      DIR=(*it)->type-4;
      PAIR=(DIR)? 4:5;
      GraphSpace::Vertex *wl=(*it)->next[DIR];
      while(1){
	
	if(wl->type==2){
	    
	  int l = wl->i/V;
	  int crr = ( wl->p < wl->next[0]->p )?  1 : -1;
	  
	  for(int d=0; d<N->d; d++) if(LT->bond_vec[l][d]!=0.0) Wi[d] += crr;
	  wl=wl->nleg;
	    
	}
	else if(wl->type==PAIR){(*it)->dir=0; wl->dir=0; break;}
	else if(wl->type==-3){
	  int l = wl->i/V;
	  int crr = ( wl->p < wl->next[0]->p )?  1 : -1;
	  for(int d=0; d<N->d; d++) if(LT->bond_vec[l][d]!=0.0) Wi[d] += crr;
	  (*it)->dir=0;
	  break;
	}
	else if(wl->type==-1||wl->type==-5){(*it)->dir=0; break;}
	
	wl=wl->next[DIR];
      }
    }
    
    ++it;
  }
 

  
  if( PR->nst==0 ){//Sum over t and s

    int m_Wi[3];
    
    for(int tag=1; tag< PR->NtNs; tag++ ){
      
      MPI_Recv(m_Wi, 3, MPI_INT, tag + PR->nst0, 0, MPI_COMM_WORLD,&status);
      for(int i=0; i<3; i++) Wi[i] += m_Wi[i];
    }
    
    values_S[wind] = ( Wi[0]*Wi[0] + Wi[1]*Wi[1] + Wi[2]*Wi[2]); ///SFD_Norm;
    
  }
  else{

    MPI_Send(Wi, 3, MPI_INT, PR->nst0, 0, MPI_COMM_WORLD );
    values_S[wind] = 0.0; //Winding
  }
 

#ifdef ATtime
  autocorrelation(mcs,values_S[wind]);
#endif

}


//*******************Ì©ÅÙ************************

void Quantities::Density(GraphSpace::Vertex *world, GraphSpace::Vertex *worldB ){

  MPI_Status status;

  for( int i=0; i<V; i++ ){

    double n0=0.0;
    for( GraphSpace::Vertex *wl=&(world[i]); wl!=&(worldB[i]); wl=wl->next[1] ){
      n0 += (wl->next[1]->t - wl->t)*(wl->p-0.5);
    }

    values_L[i] = n0/N->B;
  }


//***MPI sum for density
  double ntot=0.0;
 
  if( PR->nt==0 ){//Sum over t and s

    for(int tag=1; tag< PR->Ntdiv; tag++ ){
        MPI_Recv(m_val, V, MPI_DOUBLE, tag + PR->nt0, 0, MPI_COMM_WORLD,&status);
        for(int i=0; i<V; i++) values_L[i] += m_val[i];
      }
      for(int i=0; i<V; i++)   ntot += values_L[i];
      values_S[magz] = ntot/(double)N->V; //Density
  }
  else{
    MPI_Send(&(values_L[0]), V, MPI_DOUBLE, PR->nt0, 0, MPI_COMM_WORLD );
    values_S[magz] = 0.0; //Density
  }
  
//*******

#ifdef ATtime
  autocorrelation(mcs,values_S[magz]);
#endif

}

void Quantities::StructureFactor(){

  MPI_Status status;

  for(int i=0;i<Lmax[2];i++) Ck[i] = complex<double>(0.0,0.0);

  
  for(int k=0; k<Nkxmax; k++){    
    for(int r=0; r<V; r++){

      Ck[k] += values_L[r]*EXPrk[theta(r,k,0)];//kx=0
      Ck[k+Nkmax] += values_L[r]*EXPrk[theta(r,-k,0)]; // h.c.

      Ck[k+Nkxmax] += values_L[r]*EXPrk[theta(r,k,1)];//ky=kx
      Ck[k+Nkxmax+Nkmax] += values_L[r]*EXPrk[theta(r,-k,1)]; // h.c.
    }
  }


  int Skmax=Nkmax*2;
  if( PR->nst==0 ){//Sum over site at same tau

    for(int tag=1; tag< PR->Nsdiv; tag++ ){

        MPI_Recv(Ck_m, Skmax, MPI_DOUBLE_COMPLEX, tag*PR->Ntdiv + PR->nst0, 0, MPI_COMM_WORLD,&status);        
	for(int k=0;k<Skmax;k++){
	  Ck[k] += Ck_m[k];
	}
    }
    for(int k=0;k<Nkmax;k++){
      
      Ck[k]*=Ck[k+Nkmax]; //Ck*Ck^{\dagger}=|Ck|^2
      values_L[f_sk(k)] = real(Ck[k])/(double)N->V;
      values_L[f_sk(k+Nkmax)] = imag(Ck[k])/(double)N->V;
    }
  }
  else if( PR->nt==0 ){
    MPI_Send(Ck, Skmax, MPI_DOUBLE_COMPLEX, PR->nst0, 0, MPI_COMM_WORLD );
  }

}

void Quantities::Compressibility(){

  MCmean_S[comp*2] = N->V * (MCmean_S[magz*2+1] / (MCmean_S[magz*2]*MCmean_S[magz*2]) - 1.0);
  MCmean_S[comp*2+1] = MCmean_S[comp*2]*MCmean_S[comp*2];

}


//*******************¥¨¥Í¥ë¥®¡¼************************

void Quantities::Energy( int countv, int mcs ){


   values_S[ene] = (sp->Eu + sp->Et - (double)countv/PR->B)/(double)(N->V*PR->Ntdiv);

#ifdef ATtime
  autocorrelation(mcs,values_S[ene]);
#endif

}


//*******************´¶¼õÎ¨************************

void Quantities::NumberOfVertices( int countv, int mcs  ){

  values_S[nver] = countv;

#ifdef ATtime
  autocorrelation(mcs,values_S[nver]);
#endif

}

void Quantities::SpecificHeat(){

  MCmean_S[spe*2]   = (MCmean_S[nver*2+1] - MCmean_S[nver*2]*MCmean_S[nver*2] - MCmean_S[nver*2])/(double)N->V;
  MCmean_S[spe*2+1] = MCmean_S[spe*2]*MCmean_S[spe*2];

}


//*******************¥ï¡¼¥à¿ô************************

void Quantities::NumberOfWorms( int m_Wnum, int mcs  ){


  values_S[nwor] = m_Wnum;

#ifdef ATtime
  autocorrelation(mcs,values_S[nwor]);
#endif

}

//*******************¥ï¡¼¥à¿ô************************

void Quantities::NumberOfKinks( int Nk, int mcs  ){


  values_S[nkin] = Nk;

#ifdef ATtime
  autocorrelation(mcs,values_S[nkin]);
#endif

}


//*******************´¶¼õÎ¨************************

void Quantities::Susceptibility(){

  MCmean_S[susx*2]   = (MCmean_S[nwor*2+1] - MCmean_S[nwor*2]*MCmean_S[nwor*2] - MCmean_S[nwor*2])/(N->V*N->B*sp->Htr*sp->Htr);
  MCmean_S[susx*2+1] = MCmean_S[susx*2]*MCmean_S[susx*2];

}


//*******************Ä¶Î®Æ°OP¡¢¥°¥ê¡¼¥ó´Ø¿ô************************

void Quantities::CondensateFraction(int mcs, GraphSpace::Vertex *world){


  double atot=0.0, ctot=0.0;

  for( int i=0; i<V; i++ ){

    atot+=an[i];
    ctot+=cr[i];
   }

  atot/=N->V;
  ctot/=N->V;


  values_S[magx] = (atot+ctot)/2.0;
  values_S[magx+1] = atot;
  values_S[magx+2] = ctot;

#ifdef CFOUT
  CorrelationFunction1();
  CorrelationFunction2(world);
#endif
#ifdef ATtime
  autocorrelation(mcs,values_S[magx]);
#endif

}

void Quantities::CorrelationFunction1( ){//real space

  MPI_Status status;

//***MPI sum for G
  double ntot=0.0;
  int ixmax=Lmax[1]/PR->Nxdiv;
  double Norm;


  for(int i=0; i<Lmax[1]; i++) values_L[f_gf(i)]=0.0;


  if( PR->nx==0 ){//Sum over x at same tau

    for(int i=0;i<ixmax;i++){
      int ixx=i*Nkstep;
      this->Q[i] = an[ixx];
      this->Q[i+Lmax[1]] = cr[ixx];
      this->Q[i+2*Lmax[1]] = ca[ixx];
    }

    for(int tag=1; tag< PR->Nxdiv; tag++ ){

      MPI_Recv(m_val, ixmax*3, MPI_DOUBLE, tag*PR->Ntdiv + PR->nx0, 0, MPI_COMM_WORLD,&status);
      for(int i=0;i<ixmax;i++)  this->Q[i+tag*ixmax]=m_val[i];
      for(int i=0;i<ixmax;i++)  this->Q[i+tag*ixmax+Lmax[1]]=m_val[i+ixmax];
      for(int i=0;i<ixmax;i++)  this->Q[i+tag*ixmax+2*Lmax[1]]=m_val[i+2*ixmax];
    }

    for(int ix=0;ix<Lmax[1];ix++)
      for(int ixx=0;ixx<Lmax[1];ixx++){

	int rx= (ix - ixx + Lmax[1])%Lmax[1];
	if(ix==ixx)  values_L[f_gf(0)] += this->Q[ix+2*Lmax[1]]*2.0;
	else values_L[f_gf(rx)] += this->Q[ix]*this->Q[ixx+Lmax[1]] + this->Q[ix+Lmax[1]]*this->Q[ixx]; 
      }

    Norm = 2.0*Lmax[1]*PR->Ntdiv*PR->Nydiv*PR->Nzdiv;
    for(int i=0; i<Lmax[1]; i++) values_L[f_gf(i)]/=Norm;

  }
  else{

    for(int i=0;i<ixmax;i++){
      int ixx=i*Nkstep;
      this->Q[i] = an[ixx];
      this->Q[i+ixmax] = cr[ixx];  
      this->Q[i+2*ixmax] = cr[ixx];  
    }
    MPI_Send(this->Q, ixmax*3, MPI_DOUBLE, PR->nx0, 0, MPI_COMM_WORLD );
  }

  
  if( PR->nst==0 ){//Sum over x=0

    for(int tag=1; tag< PR->NtNs; tag++ ){

       if((tag/PR->Ntdiv)%PR->Nxdiv!=0) continue;

        MPI_Recv(m_val, Lmax[1], MPI_DOUBLE, tag + PR->nst0, 0, MPI_COMM_WORLD,&status);
        for(int i=0; i<Lmax[1]; i++) values_L[f_gf(i)] += m_val[i];
    }


  }
  else if( PR->nx==0 ){
    MPI_Send(&(values_L[f_gf(0)]), Lmax[1], MPI_DOUBLE, PR->nst0, 0, MPI_COMM_WORLD );
    }

}

void Quantities::CorrelationFunction2( GraphSpace::Vertex *world){//momentum space

  MPI_Status status;
  complex<double>Nk2;
   ////G2///

  for(int i=0;i<Nk_set;i++) Ck[i] = complex<double>(0.0,0.0);

  double rootV = sqrt((double)N->V);
  double rootV2 = N->V;
  double rootV3 = rootV2*rootV;
  double rootV4 = rootV2*rootV2;

#ifdef GF22Gk2

  for(int i=0;i<Lmax[4];i++) values_L[f_gf2(i)] = 0.0;

  double GF2;
  
  for(int r1=0; r1<Nx; r1++){
    for(int r2=0; r2<Nx; r2++){
      for(int r3=0; r3<Nx; r3++){
	for(int r4=0; r4<Nx; r4++){

	  if(r1==r2){
	    
	    if(r3==r4){
	      if(r1==r3) GF2 = ca[r1];//(r1=r2=r3=r4)
	      else  GF2 = ca[r1]*ca[r3];//(r1=r2,r3=r4)
	    }
	    else{
	      if(r1==r3)  GF2 = cr[r3]*an[r4];//(r1=r2=r3)
	      else if(r1==r4)  GF2 = 0.0;//(r1=r2=r4)
	      else  GF2 = ca[r1]*cr[r3]*an[r4];//(r1=r2)
	    }

	  }
	  else if(r1==r3) GF2 = 0.0;//(r1=r3, r2=r4), (r1=r3), (r1=r3=r4)
	  else if(r1==r4){
	    if(r2==r3) GF2 = ca[r1]*ac[r2];//(r1=r4,r2=r3)
	    else GF2 = ca[r1]*an[r2]*cr[r3];//(r1=r4)
	  }
	  else if(r2==r3) GF2 = cr[r1]*ac[r2]*an[r4];//(r2=r3)
	  else if(r2==r4) GF2 = 0.0;//(r2=r4)
	  else if(r3==r4) GF2 = cr[r1]*an[r2]*ca[r3];//(r3=r4)
	  else GF2 = cr[r1]*an[r2]*cr[r3]*an[r4];

	  values_L[f_gf2((r1-r2+Nx)%Nx, (r3-r4+Nx)%Nx)] += GF2;

	}
      }
    }
  }
  
#endif

  for(int r=0; r<V; r++){

    double ancr=an[r]*cr[r];
    double n=world[r].p;
    double ck0=ancr - n;
    double nn=1.0 - n;

    Ck[0] += ck0/rootV2;
    Ck[1] += ancr/rootV2;
    Ck[2] += (ancr*ancr - n)/rootV4;
    Ck[3] += n/rootV2;
    Ck[4] += (5.0*ancr*ancr -6.0*n*ancr -2.0*nn*ancr +2.0*ancr + n*n + n*nn)/rootV4;


    for(int a=0; a<2; a++){

      for(int k=0; k<Nkmax; k++){

	Ck[f_ck(k,0,a)] += an[r]*EXPrk[theta(r,-k,a)]/rootV;
	Ck[f_ck(k,1,a)] += cr[r]*EXPrk[theta(r,k,a)]/rootV;
   
	Ck[f_ck(k,2,a)] += an[r]*an[r]*EXPrk[theta(r,-k,a)]/rootV2;
	Ck[f_ck(k,3,a)] += cr[r]*cr[r]*EXPrk[theta(r,k,a)]/rootV2;

	Ck[f_ck(k,5,a)] += ancr*EXPrk[theta(r,k,a)]/rootV2;
	Ck[f_ck(k,7,a)] += n*EXPrk[theta(r,k,a)]/rootV2;
	Ck[f_ck(k,9,a)] += nn*EXPrk[theta(r,k,a)]/rootV2;

	Ck[f_ck(k,10,a)] += ck0*an[r]*EXPrk[theta(r,-k,a)]/rootV3;
	Ck[f_ck(k,11,a)] += ck0*cr[r]*EXPrk[theta(r,k,a)]/rootV3;

	Ck[f_ck(k,12,a)] += (ancr - nn)*an[r]*EXPrk[theta(r,-k,a)]/rootV3;
	Ck[f_ck(k,13,a)] += (ancr - nn)*cr[r]*EXPrk[theta(r,k,a)]/rootV3;

	Ck[f_ck(k,14,a)] += ancr*an[r]*EXPrk[theta(r,-k,a)]/rootV3;
	Ck[f_ck(k,15,a)] += ancr*cr[r]*EXPrk[theta(r,k,a)]/rootV3;

      }
      for(int k=1; k<Nkmax; k++){
	Ck[f_ck(-k,5,a)] += ancr*EXPrk[theta(r,-k,a)]/rootV2;//ck(4)<=ck(-k)<=ck(5)
	Ck[f_ck(-k,7,a)] += n*EXPrk[theta(r,-k,a)]/rootV2;
	Ck[f_ck(-k,9,a)] += nn*EXPrk[theta(r,-k,a)]/rootV2;
      }
    }
  }


  if( PR->ns==0 ){//Sum over site at same tau

    for(int tag=1; tag< PR->Nsdiv; tag++ ){

        MPI_Recv(Ck_m, Nk_set, MPI_DOUBLE_COMPLEX, tag*PR->Ntdiv + PR->ns0, 0, MPI_COMM_WORLD,&status);        
	for(int i=0;i<Nk_set;i++) Ck[i] += Ck_m[i];
    }

    ///////////////// Gk //////////////////////
    for(int a=0; a<2; a++){
      for(int k=0; k<Nkxmax; k++){
        Nk2 = (Ck[f_ck(k,0,a)]*Ck[f_ck(k,1,a)] - Ck[0])/(double)N->V;
        values_L[f_nkr(k+a*Nkxmax)] = real(Nk2);
        values_L[f_nki(k+a*Nkxmax)] = imag(Nk2);
     }
    }
    //////////////////////////////////////////


  }
  else{

    MPI_Send(Ck, Nk_set, MPI_DOUBLE_COMPLEX, PR->ns0, 0, MPI_COMM_WORLD );
  }


  //
  complex<double> Gk2, Dkkk;
  if( PR->ns==0 ){


      for(int kx=0; kx<Nkxmax; kx++){ //(kx,ky,kz)=(0~Nkmax,0,0); (0~Nkmax,kx,kx)

	for(int kk=0; kk<Nkkmax; kk++){  //(k'x,k'y,k'z)=(0,0,0);...(depends on CKK)...; (kx,ky,kz) 

#ifdef GF22Gk2
	  Gk2 = complex<double>(0.0,0.0);	  
	  for(int r=0;r<Nx;r++) for(int rr=0; rr<Nx; rr++)
				 Gk2 += values_L[f_gf2(r,rr)]*EXPrk[theta(r,kx,0)]*EXPrk[theta(rr,kk,0)];
	  Dkkk = Gk2/(double)(Nx*Nx);
#else

	  Gk2 = Ck[f_ck(kx,0)]*Ck[f_ck(kx,1)]*Ck[f_ck(kk,0)]*Ck[f_ck(kk,1)]

	    - ( Ck[0]*(Ck[f_ck(kk,1)]*Ck[f_ck(kk,0)]-Ck[1]) - Ck[f_ck(kk,11)]*Ck[f_ck(kk,0)] - Ck[f_ck(kk,10)]*Ck[f_ck(kk,1)] ) //i=j
	    - ( Ck[f_ck(kk+kx,3)]*(Ck[f_ck(kx,0)]*Ck[f_ck(kk,0)] - Ck[f_ck(kk+kx,2)]) - Ck[f_ck(kk,15)]*Ck[f_ck(kk,0)] - Ck[f_ck(kx,15)]*Ck[f_ck(kx,0)] ) //i=l
	    - ( (Ck[f_ck(kx-kk,5)]-Ck[f_ck(kx-kk,7)])*(Ck[f_ck(kx,0)]*Ck[f_ck(kk,1)] - Ck[f_ck(kk-kx,5)]) - Ck[f_ck(kk,10)]*Ck[f_ck(kk,1)] - Ck[f_ck(kx,11)]*Ck[f_ck(kx,0)] ) //i=m
	    - ( Ck[f_ck(kx+kk,2)]*(Ck[f_ck(kx,1)]*Ck[f_ck(kk,1)] - Ck[f_ck(kk+kx,3)]) - Ck[f_ck(kx,14)]*Ck[f_ck(kx,1)] - Ck[f_ck(kk,14)]*Ck[f_ck(kk,1)] ) //j=m
	    - ( (Ck[f_ck(kk-kx,5)]-Ck[f_ck(kk-kx,9)])*(Ck[f_ck(kx,1)]*Ck[f_ck(kk,0)] - Ck[f_ck(kx-kk,5)]) - Ck[f_ck(kx,12)]*Ck[f_ck(kx,1)] - Ck[f_ck(kk,13)]*Ck[f_ck(kk,0)] ) //j=l
	    - ( Ck[0]*(Ck[f_ck(kx,1)]*Ck[f_ck(kx,0)]-Ck[1]) - Ck[f_ck(kx,11)]*Ck[f_ck(kx,0)] - Ck[f_ck(kx,10)]*Ck[f_ck(kx,1)] ) //l=m

	    -  (Ck[f_ck(kk,15)] - Ck[f_ck(kk,1)]/rootV2)*Ck[f_ck(kk,0)] //i=j=l
	    -  Ck[f_ck(kk,14)]*Ck[f_ck(kk,1)] //i=j=m
	    -  Ck[f_ck(kx,15)]*Ck[f_ck(kx,0)] //i=m=l
	    -  (Ck[f_ck(kx,14)] - Ck[f_ck(kx,0)]/rootV2)*Ck[f_ck(kx,1)] //j=l=m

	    -  (Ck[1]*Ck[1] - Ck[3]*Ck[3]) //i=j,l=m
	    -  (Ck[f_ck(kx-kk,5)]*Ck[f_ck(kk-kx,5)] - Ck[f_ck(kx-kk,7)]*Ck[f_ck(kk-kx,9)])  //i=m,j=l
	    -  Ck[f_ck(kk+kx,3)]*Ck[f_ck(kk+kx,2)] //i=l,j=m
	    -  Ck[2] //i=j=l=m
	    -  Ck[4];


	  Dkkk = Gk2; //vol2;
#endif	  

	  //	  for(int kky=0; kky<1; kky++){ //ky=kky= 0 or kx(kkx)
	    values_L[f_gk2r(kx,kk)]=real(Dkkk);
	    values_L[f_gk2i(kx,kk)]=imag(Dkkk);
	    // }
	    
	}
      }
      
  }


//*******
  int L3=Lmax[3]+Lmax[4]+Lmax[5]; //for nk, gf2, gk2

  if( PR->nst==0 ){//Sum over t=0

    for(int tag=1; tag< PR->Ntdiv; tag++ ){

        MPI_Recv(m_val, L3, MPI_DOUBLE, tag + PR->nst0, 0, MPI_COMM_WORLD,&status);
        for(int i=0; i<L3; i++) values_L[f_nk(i)] += m_val[i];
    }

    for(int i=0; i<L3; i++) values_L[f_nk(i)]/=(double)PR->Ntdiv;

  }
  else if( PR->ns==0 ){//Ntdiv prosessors for nk, gf2, gk2
    MPI_Send(&(values_L[f_nk(0)]), L3, MPI_DOUBLE, PR->nst0, 0, MPI_COMM_WORLD );
  }


}


/*
void Quantities::StructureFactor(){

  MPI_Status status;

  for(int i=0;i<Lmax[2];i++) Ck[i] = complex<double>(0.0,0.0);

  for(int k=0; k<Nkxmax; k++){    
    for(int r=0; r<V; r++){

      Ck[k] += values_L[r]*EXPrk[theta(r,k,0)];//kx=0
      Ck[k+Nkmax] += values_L[r]*EXPrk[theta(r,-k,0)]; // h.c.

      Ck[k+Nkxmax] += values_L[r]*EXPrk[theta(r,k,1)];//ky=kx
      Ck[k+Nkxmax+Nkmax] += values_L[r]*EXPrk[theta(r,-k,1)]; // h.c.
    }
  }


  int Skmax=Nkmax*2;
  if( PR->nst==0 ){//Sum over site at same tau

    for(int tag=1; tag< PR->NtNs; tag++ ){

        MPI_Recv(Ck_m, Skmax, MPI_DOUBLE_COMPLEX, tag + PR->nst0, 0, MPI_COMM_WORLD,&status);        
	for(int k=0;k<Skmax;k++){
	  Ck[k] += Ck_m[k];
	}
    }
    for(int k=0;k<Nkmax;k++){
      Ck[k]*=Ck[k+Nkmax]; //Ck*Ck^{\dagger}=|Ck|^2
      values_L[f_sk(k)] = real(Ck[k])/(double)(N->V*PR->Ntdiv);
      values_L[f_sk(k+Nkmax)] = imag(Ck[k])/(double)(N->V*PR->Ntdiv);
    }
  }
  else{
    MPI_Send(Ck, Skmax, MPI_DOUBLE_COMPLEX, PR->nst0, 0, MPI_COMM_WORLD );
  }

}
*/

void Quantities::NoiseCorrelation( ){

  for(int kx=0;kx<Nkxmax;kx++){
    for(int kk=0;kk<Nkkmax;kk++){

      MCmean_L[f_noise(kx,kk)] =  MCmean_L[f_gk2r(kx,kk)] -  MCmean_L[f_nkr(kx)]*MCmean_L[f_nkr(kk)];

    }
  }


}



//****MC sum
void Quantities::MCsum_S(){

//  TreeSum( Nq, values_S, MC->Nsample , PR->nst, PR->nst0, PR->NtNs );
 
  MPI_Status status;
 
  if ( PR->nst==0 ){//Sum over t and s
      
    for(int tag=1; tag< PR->NtNs; tag++ ){

      MPI_Recv(m_val, Nq, MPI_DOUBLE, tag + PR->nst0, 0, MPI_COMM_WORLD,&status);
      for(int i=0; i<Nq; i++){
	if( tag%PR->Ntdiv==0 || i<Nend ) values_S[i] += m_val[i];
      }
    }

  }
  else{
    MPI_Send(values_S, Nq, MPI_DOUBLE, PR->nst0, 0, MPI_COMM_WORLD );
  }


  if( PR->nst==0 ){

    for(int i=0; i<Nq; i++){
      MCmean_S[2*i] += values_S[i]/(double)MC->Nsample;
      MCmean_S[2*i+1] += values_S[i]*values_S[i]/(double)MC->Nsample;
    }
    
  }


}



void Quantities::MCsum_L(){


  if( PR->nst==0 ){
    for(int i=0; i<Lsize; i++){
      MCmean_L[i] += values_L[i]/(double)MC->Nsample;
    }
  }


}

void Quantities::TreeSum( int num, double *child, double *mother, double Normalization, int cnum, int mnum, int Nsum ){

  MPI_Status status;

  double *m_val;
  newcall(m_val,num);
 
  for(int i=0; i<num; i++) child[i]/=Nsum;

  int tag=1, global_i;
  int logNsum=log2((double)Nsum);

  for(int i=0; i< logNsum; i++ ){

    if( ( cnum%(2*tag) ) == 0 ){

        MPI_Recv(m_val, num, MPI_DOUBLE, cnum + tag  + mnum, 0, MPI_COMM_WORLD,&status);

        for(int j=0; j<num; j++) mother[global_i] += m_val[j];

        tag=tag*2;

    }
    else{

      MPI_Send(child, num, MPI_DOUBLE,  cnum - tag  + mnum, 0, MPI_COMM_WORLD );
      break;
    }


  }


  delcall(m_val);


}

void Quantities::TreeSum( int num, double *child, double Normalization, int cnum, int mnum, int Nsum ){

  MPI_Status status;

  double *m_val;
  newcall(m_val,num);
 
  for(int i=0; i<num; i++) child[i]/=Nsum;

  int tag=1;
  int logNsum=log2((double)Nsum);

  for(int i=0; i< logNsum; i++ ){

    if( ( cnum%(2*tag) ) == 0 ){

        MPI_Recv(m_val, num, MPI_DOUBLE, cnum + tag  + mnum, 0, MPI_COMM_WORLD,&status);

        for(int j=0; j<num; j++) child[j] += m_val[j];

        tag=tag*2;

    }
    else{

      MPI_Send(child, num, MPI_DOUBLE,  cnum - tag  + mnum, 0, MPI_COMM_WORLD );
      break;
    }


  }


  delcall(m_val);


}


/*
void Quantities::AutoCorrelationAverage( char *fname ){


  MPI_Status status;
  FILE *file;
  double *i_AC;
  double *m_AC,*n_AC,*RNDmeanm_AC,*ACT,A=0.0,s;
  int Nmcs=MC->Nsample;


  newcall(i_AC,Nmcs);
  newcall(m_AC,Nmcs);
  newcall(n_AC,Nmcs);
  newcall(RNDmeanm_AC,Nmcs);
  newcall(ACT,2);


  if ( PR->nst==0 ){//Sum over t and s

     for(int tag=1; tag< PR->NtNs; tag++ ){

        MPI_Recv(i_AC, Nmcs, MPI_DOUBLE, tag + PR->nst0, 0, MPI_COMM_WORLD,&status);
        for(int i=0; i< Nmcs; i++ ) AC[i] += i_AC[i];
     }

    //############################################
    for(int i=0; i< Nmcs; i++ ){

       s=0.0;
       for(int j=0; j< Nmcs; j++ ){ 
              s += AC[j]*AC[(i+j)%Nmcs];  
        }
        n_AC[i]=s;
        A += AC[i];
    }


    for(int i=0; i< Nmcs; i++ ) m_AC[i] = ( n_AC[i]  - A*A/(double)Nmcs )/( n_AC[0] - A*A/(double)Nmcs ) ;

    ACT[0]=0.0;
    for(int i=0; i< Nmcs; i++ ) {
      if(m_AC[i]<0) break;
      ACT[0] += m_AC[i];
    }
    ACT[1]=ACT[0]*ACT[0];

  //########################################

  }
  else{
      MPI_Send(AC, Nmcs, MPI_DOUBLE, PR->nst0, 0, MPI_COMM_WORLD );
   
  }


  MPI_Reduce(ACT, RNDmeanACT, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(m_AC, RNDmeanm_AC, Nmcs, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(PR->my_rank==0){

      file = fopen(fname,"a+");

      for(int i=0; i< Nmcs; i++ ){

        if(i%100==0) fprintf(file,"%d %le\n", i, RNDmeanm_AC[i]/(double)(PR->Npara*PR->NtNs) );

     }


      fclose(file);

  }

 
  delcall(i_AC);
  delcall(n_AC);
  delcall(m_AC);
  delcall(RNDmeanm_AC);
  delcall(ACT);


}
*/


