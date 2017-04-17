//######################################################################
//####
//####  World-Line Monte Carlo simulation
//####                       by the Directed-Loop Algorithm 
//####
//####                                 Mar.03 / 2005, Naoki Kawashima
//####
//######################################################################

//######################################################################
//####
//####  World-Line Monte Carlo simulation
//####                       by the non-Vertex Directed-Loop Algorithm 
//####
//####                                 Nov.11 / 2007, Yasuyuki Kato
//####
//######################################################################

/*! \mainpage DSQSS-1.1 プログラムドキュメント
 *
 * \section intro Introduction
 *
 * DSQSSは、連続虚数時間向き付きループアルゴリズムに基づく量子モンテカルロ法によって離散空間上で定義された量子多体問題を解くためのプログラムです。。単位格子の情報、２体相互作用ハミルトニアンの行列要素などを入力ファイルとして広い範囲のモデルに対応しています。例えば、次元、格子の１辺の長さ、相互作用の異方性、スピンの大きさ、磁場の大きさ、温度などのパラメータを任意にとって、ＸＸＺハイゼンベルクモデルの計算を行うことができます。また、ボーズ系のシミュレーションも可能です。

 * \section install Program Document
 *
 * DSQSSは、C++でプログラミングされています。ここでは、クラス階層、構成メンバなどのプログラム構造について説明しています。  
 *  
 *  
 * \section release Release 
 *  
 * 2011.3.31 DSQSSバージョン1.0 リリース
 *  
 * 2011.9.30 DSQSSバージョン1.1 リリース
 *
 */


#define Version "1.1.19"

//#define CF //Computing G(tau)

//#include <mpi.h>
#include "dla.hpp"

#ifdef GRAPHIC
#include "graphic.cc"
#endif

#ifdef CHAINJOB
#define killtime 3600 //(sec.)
//#define killtime 32400
#include "chainjob.hpp"
#endif

#ifdef PAPI
#include <papi.h>
#endif

//######################################################################

//int main(int argc, char** argv) {
int main(int argc, char* argv[]) {

  // MPI  Initialization & Getting Size and Rank
#ifdef MULTI
#ifdef KASHIWA
    MPI_Init( &argc , &argv );
    MPI_Comm_size( MPI_COMM_WORLD , &N_PROC );
    MPI_Comm_rank( MPI_COMM_WORLD , &I_PROC );
#else
    MPI::Init(argc,argv);
    I_PROC = MPI::COMM_WORLD.Get_rank();
    N_PROC = MPI::COMM_WORLD.Get_size();
#endif

  //  int itag = 0;
  printf(
	 ">>> The program is being run with MPI mode.( IP= %3d / NP= %3d ) \n\n", 
	 I_PROC, N_PROC
	 );
#endif

#ifdef DEB
  printf("\n\n>>> The program is being run in DEBUG mode.\n\n\n");
#endif

#ifdef GRAPHIC
  printf(">>> EGGX library is used.\n");
#endif

#ifdef CHAINJOB
  printf(">>CHAIN JOB\n");
  gettimeofday(&s_time, NULL);
#endif

  Parameter P(argc,argv);

  TheSegmentPool.init(P.NSEGMAX);
  TheVertexPool.init(P.NVERMAX);

  RND.setSeed(P.SEED,32);

#ifdef PAPI
  int ierr;
  long_long nflops;
  float real_time,proc_time,mflops;
#endif
  
#ifdef PAPI
  ierr = PAPI_flips(&real_time,&proc_time,&nflops,&mflops);
#endif

  Simulation Sim(P);

#ifdef PAPI
  ierr = PAPI_flips(&real_time,&proc_time,&nflops,&mflops);

  cout << "real_time = "<< real_time << endl;
  cout << "proc_time = "<< proc_time << endl;
  cout << "MFLOPS    = "<< mflops << endl;
#endif 

#ifdef MULTI
#ifdef KASHIWA
  MPI_Finalize();
#else
  MPI::Finalize();
#endif
#endif

}
//
//######################################################################
 

inline Simulation::Simulation(Parameter& P0) :
  P(P0), ALG(P.ALGFILE), LAT(P.LATFILE, ALG), MSR(P,LAT,ALG){

printf("[%2d] Pass 1\n", I_PROC);  //koko

#ifdef GRAPHIC
  g_init();
#endif
  reset_counters();

#ifdef CF
  counter4CF  = new int*    [MSR.NCF+1];
  cfsamp      = new double* [MSR.NCF+1];
  for(int icf=0;icf<MSR.NCF+1;icf++){
    counter4CF[icf] = new int    [MSR.Ntau];
    cfsamp[icf]     = new double [MSR.Ntau];
  }
#endif
#ifdef CK
  counter4CKC  = new double* [MSR.NCK+1];
  counter4CKS  = new double* [MSR.NCK+1];
  cksamp  = new double* [MSR.NCK+1];
  for(int ick=0;ick<MSR.NCK+1;ick++){
    counter4CKC[ick] = new double [MSR.Ntau1];
    counter4CKS[ick] = new double [MSR.Ntau1];
    cksamp[ick] = new double [MSR.Ntau];
  }
#endif

#ifdef CHAINJOB
  BinaryIO();
  if(!isChainjob){
    set_NCYC();
  }else{
    gettimeofday(&e_time, NULL);
    cout<<"Reading bin file takes "<<dtime(s_time,e_time)<<" sec."<<endl;
  }
  isEnd = false;
#else
  if (P.RUNTYPE==0)set_NCYC();
#endif

  int index;

  if (P.RUNTYPE > 0 ){
#ifdef MULTI
 
    // =========================  edit  sakakura  =========================
    // =========================  edit  sakakura  =========================

    int nkink,nkink0,nkink1;
    double ediag,ediag0,ediag1;

    double dee;
//    bool judge;
    int judge;
    int index0,index1,indexx;
    double hh[P.NREP];

    index  = I_PROC;

printf("[%2d] Pass 2\n", I_PROC);  //koko

    ///  +++++ replica exchange for F ++++ ///
    Algorithm* alg = new Algorithm[P.NREP];

printf("[%2d] Pass 3\n", I_PROC);  //koko

    for (int i=0; i<P.NREP; i++ ){
      alg[i].set_i(i);
      alg[i].set_alg(i);
    }

    if (P.RUNTYPE == 1){
      for (int i=0; i<P.NREP; i++) {
	hh[i]=P.VF + P.DF*i;
      }
      ///  +++++ replica exchange for F ++++ ///
    }

    if (P.RUNTYPE == 2){
      ///  +++++ replica exchange for B ++++ ///
      for (int i=0; i<P.NREP; i++) {
	hh[i]=P.VB + P.DB*i;
      }
      ///  +++++ replica exchange for B ++++ ///
    }  

    if (P.RUNTYPE == 3){
      ///  +++++ nontrivial parallel ++++ ///

printf("[%2d] Pass 4\n", I_PROC);  //masaki
    }  

    if (P.RUNTYPE == 4){
      ///  +++++ simple parameter parallel ++++ ///

printf("[%2d] Pass 4\n", I_PROC);  //koko
    } 

    //set_NCYC();
    int is,ie,ikey,i;

    ikey = 1;
    is = 1;
    ie = P.NREP-2;
    //  judge = false;
    judge = 1;//true;
    index = I_PROC;

    MPI_Status recv_status;

    for (ISET=ISETstart; ISET<P.NSET; ISET++) {
      for (i=is; i<ie; i+=2) {
///     if (judge == true && ISET !=ISETstart){
        if (judge == 1 && ISET !=ISETstart){

#ifdef KASHIWA
  	  if (I_PROC == i+1) MPI_Send(&index , 1, MPI_INT, i  , 0, MPI_COMM_WORLD);
  	  if (I_PROC == i  ) MPI_Recv(&index0 , 1, MPI_INT, i+1  , 0, MPI_COMM_WORLD, &recv_status);
  	  if (I_PROC == i  ) MPI_Send(&index , 1, MPI_INT, i+1, 0, MPI_COMM_WORLD);
  	  if (I_PROC == i+1) MPI_Recv(&index1, 1, MPI_INT, i  , 0, MPI_COMM_WORLD, &recv_status);
#else
          if (I_PROC == i+1) MPI::COMM_WORLD.Send(&index , 1, MPI::INT, i  , 0);
          if (I_PROC == i  ) MPI::COMM_WORLD.Recv(&index0, 1, MPI::INT, i+1, 0);
          if (I_PROC == i  ) MPI::COMM_WORLD.Send(&index , 1, MPI::INT, i+1, 0);
          if (I_PROC == i+1) MPI::COMM_WORLD.Recv(&index1, 1, MPI::INT, i  , 0);
#endif

	  if (I_PROC == i   ) index=index0;
	  if (I_PROC == i+1 ) index=index1;
	}
      }

printf("[%2d] Pass 5\n", I_PROC);  //koko
    //if (judge == true) {
      if (judge == 1) {
	if (P.RUNTYPE==1) RepChangeF(alg[index]);
        if (P.RUNTYPE==2) RepChangeB(hh[index]);
        if (P.RUNTYPE==3) {}
        if (P.RUNTYPE==4) {}

	set_NCYC();
      }

      if (ikey == 0) {
	ikey = 1;
	is = 1; ie=P.NREP-2;
      }
      else{
	ikey = 0;
	is = 0; ie=P.NREP-1;
      }

      Set(index);
      MSR.summary();
      //  ene = MSR.get_data();

      ediag = MSR.get_ediag();

      if (P.RUNTYPE==1){
	ediag = MSR.get_ediag();
      }   

      if (P.RUNTYPE==2){
	ediag = MSR.get_ediag();
	nkink = MSR.get_nkink();
      }   

      if (P.RUNTYPE==3){
      }

      if (P.RUNTYPE==4){
      }

//    judge = false;
      judge = 0;

printf("[%2d] Pass 6\n", I_PROC);  //koko

      if (ISET == P.NSET-1) break;

      for (i=is; i<ie; i+=2) {

	if (P.RUNTYPE==2){
#ifdef KASHIWA
	  if (I_PROC == i   ) MPI_Send(&nkink , 1, MPI_INT, i+1,0, MPI_COMM_WORLD);
	  if (I_PROC == i+1 ) MPI_Recv(&nkink0, 1, MPI_INT, i  ,0, MPI_COMM_WORLD, &recv_status);
#else
          if (I_PROC == i   ) MPI::COMM_WORLD.Send(&nkink , 1, MPI::INT, i+1,0);
          if (I_PROC == i+1 ) MPI::COMM_WORLD.Recv(&nkink0, 1, MPI::INT, i  ,0);
#endif
	}

#ifdef KASHIWA
  	if (I_PROC == i  ) MPI_Send(&ediag , 1, MPI_DOUBLE, i+1,0, MPI_COMM_WORLD);
  	if (I_PROC == i+1) MPI_Recv(&ediag0, 1, MPI_DOUBLE, i  ,0, MPI_COMM_WORLD, &recv_status);
  	if (I_PROC == i  ) MPI_Send(&index , 1, MPI_INT, i+1, 0, MPI_COMM_WORLD);
  	if (I_PROC == i+1) MPI_Recv(&index0, 1, MPI_INT, i  , 0, MPI_COMM_WORLD, &recv_status);
#else
        if (I_PROC == i  ) MPI::COMM_WORLD.Send(&ediag , 1, MPI::DOUBLE, i+1,0);
        if (I_PROC == i+1) MPI::COMM_WORLD.Recv(&ediag0, 1, MPI::DOUBLE, i  ,0);
        if (I_PROC == i  ) MPI::COMM_WORLD.Send(&index , 1, MPI::INT, i+1, 0);
        if (I_PROC == i+1) MPI::COMM_WORLD.Recv(&index0, 1, MPI::INT, i  , 0);
#endif

	if (I_PROC == i+1 ) {

	  if (P.RUNTYPE == 1)dee = exp(  LAT.BETA * (ediag - ediag0) \
					 * ( hh[index] - hh[index0]) ) ;

	  if (P.RUNTYPE == 2)dee = exp( (hh[index]-hh[index0])*(ediag-ediag0) )\
			       * ( pow( hh[index0]/hh[index] ,nkink-nkink0 ) );
 
	  if (P.RUNTYPE == 3) {}

	  if (P.RUNTYPE == 4) {}

////////  if (dee > RND.Uniform()) judge = true;
	  if (dee > RND.Uniform()) judge = 1;
	  // for debug//       judge=false;

	  /// ++++for debug ++++
	  //  printf("rank=%d index(i,i+1)=%d %d judge=%d\n",\
	  //  I_PROC,index0,index,judge);
	  /// ++++for debug ++++

	}

#ifdef KASHIWA
  	if (I_PROC == i+1 ) MPI_Send(&judge , 1, MPI_INT, i,   0, MPI_COMM_WORLD);
  	if (I_PROC == i   ) MPI_Recv(&judge,  1, MPI_INT, i+1, 0, MPI_COMM_WORLD, &recv_status);
#else 
        if (I_PROC == i+1 ) MPI::COMM_WORLD.Send(&judge , 1, MPI::INT, i,   0);
        if (I_PROC == i   ) MPI::COMM_WORLD.Recv(&judge,  1, MPI::INT, i+1, 0);
//      if (I_PROC == i+1 ) MPI::COMM_WORLD.Send(&judge , 1, MPI::BOOL, i,   0);
//      if (I_PROC == i   ) MPI::COMM_WORLD.Recv(&judge,  1, MPI::BOOL, i+1, 0);
#endif
      }

    }

printf("[%2d] Pass 7\n", I_PROC);  //koko
    MSR.summary_REP(index);

    ISETstart = 0;

    // =========================  edit  sakakura  =========================
 
#endif

  }else{
    //// =========================  add sakakura  =========================
    int idum=0;
    for (ISET=ISETstart; ISET<P.NSET; ISET++) Set(idum);
    ISETstart=0;
    MSR.summary();
  }
  //// =========================

  P.openfile();
  fprintf(P.FOUT,"C This is DSQSS ver.%s\n\n",Version);
  LAT.show_param(P.FOUT);
  P.dump(P.FOUT);
  MSR.show(P.FOUT);

  int ns_used = TheSegmentPool.number_of_used_elements();
  int nv_used = TheVertexPool.number_of_used_elements();
  int nrvi_used = TheRVIPool.number_of_used_elements();
  fprintf(P.FOUT,"I [the maximum number of segments]          = %d\n", ns_used);
  fprintf(P.FOUT,"I [the maximum number of vertices]          = %d\n", nv_used);
  fprintf(P.FOUT,"I [the maximum number of reg. vertex info.] = %d\n", nrvi_used);
  // === edit sakakura ===
  if (P.RUNTYPE > 0) {
    fprintf(P.FOUT,"I [the processor ID]                        = %d\n", index);
  } 

#ifdef CF
  fprintf(P.FOUT4CF,"C This is dla.cc ver.%d\n\n",Version); 
  LAT.show_param(P.FOUT4CF);
  P.dump(P.FOUT4CF);
  MSR.show4CF(P.FOUT4CF);
#endif
#ifdef SF
  fprintf(P.FOUT4SF,"C This is dla.cc ver.%d\n\n",Version); 
  LAT.show_param(P.FOUT4SF);
  P.dump(P.FOUT4SF);
  MSR.show4SF(P.FOUT4SF);
#endif
#ifdef CK
  fprintf(P.FOUT4CK,"C This is dla.cc ver.%d\n\n",Version); 
  LAT.show_param(P.FOUT4CK);
  P.dump(P.FOUT4CK);
  MSR.show4CK(P.FOUT4CK);
#endif

  P.closefile();

#ifdef GRAPHIC
  g_clear();
#endif

#ifdef CHAINJOB
  isEnd = true;
  cjobout = new ofstream();
  (*cjobout).open(CJOBFILE, ios::out | ios::binary);//reset
  end_cjob();
#endif

}

//======================================================================

void Simulation::reset_counters() {
  ISET = -1;
  IMCSE = -1;
  IMCSD = -1;
  IMCS = -1;
  ICYC = -1;

  ISETstart  = 0;
  IMCSDstart = 0;
  IMCSstart  = 0;

}

//======================================================================

inline Simulation::~Simulation() {
  //  printf("*** Destroying Simulation\n");
#ifdef CF
    for(int icf=0; icf<MSR.NCF+1; icf++){
      delete [] counter4CF[icf] ; 
      delete [] cfsamp[icf] ;
  }
  delete [] counter4CF; 
  delete [] cfsamp ;
#endif
#ifdef CK
  for(int ick=0; ick<MSR.NCK+1; ick++){
    delete [] counter4CKC[ick] ; 
    delete [] counter4CKS[ick] ;
    delete [] cksamp[ick] ;
    
  }
  delete [] counter4CKC; 
  delete [] counter4CKS;
  delete [] cksamp; 
#endif

}

//======================================================================

void Simulation::set_NCYC() {
  // edit Suzuki
  // 逆温度に対する規格化
#ifdef NORM
  double vol = (double)LAT.NSITE;
#else
  double vol = LAT.BETA * (double)LAT.NSITE;
#endif
  // end edit Suzuki

  double path;
  int ncyc=1;

  int NSAMP = P.NMCSE/10;
  int*  ncycSAMP = new int[NSAMP]; //edit sakakura
  //  int ncycSAMP[NSAMP];
  
  for (IMCSE=0; IMCSE<P.NMCSE; IMCSE++) {

    path = 0.0;
    ncyc = 0;
    // 次の行の LAT.NSITE == 2 && ncyc%2 == 0 の指定は
    // S=1/2 反強磁性ハイゼンベルクモデルの特殊な非エルゴード性
    // を排除するため
    while ( path < vol || ( LAT.NSITE==2 && ncyc%2==0 ) ) {
      double p = Cycle();
      path += p;
      ncyc++;
    }

    ncycSAMP[IMCSE%NSAMP]= ncyc;
  };

  int ncycSUM = 0;
  for(int isamp = 0; isamp < NSAMP; isamp++){
    ncycSUM += ncycSAMP[isamp];
  }

  P.NCYC = ncycSUM / NSAMP ;
#ifdef GRAPHIC
  g_draw();
#endif
}

//======================================================================

// === edit rep sakakrua ===
//void Simulation::Set() {
void Simulation::Set(int ix) {
  // === edit rep sakakrua ===

  ICYC=-1;

  for (IMCSD=IMCSDstart; IMCSD<P.NMCSD; IMCSD++) {
#ifdef CHAINJOB
    gettimeofday(&e_time, NULL);
    if(dtime(s_time,e_time) > killtime ){IMCS=0;write_cjob();end_job();}
    isChainjob=false;
#endif
    Sweep();
  }
  IMCSDstart = 0;
#ifdef CHAINJOB
  if(!isChainjob){MSR.setinit();}else{isChainjob=false;}
#else
  MSR.setinit();//kota for (int i=0; i<NACC; i++) ACC[i].reset();
#endif

  for (IMCS=IMCSstart; IMCS<P.NMCS;  IMCS++) {
#ifdef CHAINJOB
    gettimeofday(&e_time, NULL);
    if(dtime(s_time,e_time) > killtime ){write_cjob();end_job();}
#endif
#if defined(CF) || defined(CK)
    Sweep_CF();
#else
    Sweep();
#endif
    MSR.measure();
    MSR.accumulate_length(AMP);
#ifdef CF
    MSR.accumulate4CF( cfsamp );
#endif
#ifdef CK
    MSR.accumulate4CK( cksamp );
#endif

  }
  IMCSstart = 0;

  // === edit rep sakakrua ===
  //MSR.setsummary();
  MSR.setsummary(ix);
  // === edit rep sakakrua ===

}

//======================================================================

void Simulation::Sweep() {

  double len = 0.0;
  for (ICYC=0; ICYC<P.NCYC; ICYC++) {
    len += Cycle();
  }
  AMP = len/(double)P.NCYC;
#ifdef GRAPHIC
  g_draw();
#endif

}

//======================================================================

double Simulation::Cycle() {
  if ( ! PlaceWorm() ) return 0.0;
  return MoveHead(); 

}

//======================================================================

bool Simulation::PlaceWorm() {

#ifdef MONITOR
  printf("\n#####################################################\n");
  printf("Simulation::PlaceWorm> Start.\n");
#endif

  // ワーム初期バーテックスの設定

  Vertex& Vorg = W.origin();
  W.setCurrentVertex(Vorg);
    
  // ワーム初期位置の決定

  int s = RND.Int(LAT.NSITE);
  //  edit Suzuki
  // 逆温度に対する補正(beta -> 1で規格化)
#ifdef NORM
  double t =  RND.Uniform();
#else
  double t = (double)LAT.BETA * RND.Uniform();
#endif
 
  // end edit Suzuki
  //kota//  乱数を用いて初期位置（サイト番号、位置）を決定

#ifdef MONITOR
  printf("  The worm is placed at (s= %d, t= %8.3f)\n", s+1, t);
#endif

  // ワームテールの初期化
  Site& ST = LAT.S(s);

  SiteProperty& SP = ST.Property();

  VertexProperty& VP = SP.getVertexProperty();
  //サイト、サイトプロパティ、バーテックスプロパティの各オブジェクトの設定
  Vorg.BareVertex::init( t , VP );


  //バーテックス初期化
  Vorg.setTime(t);
  //tをセット

  // セグメントを探す
  Segment& S0 = ST.findS(t);
  // 初期方向
  int x = S0.X();

  SiteInitialConfiguration& IC = SP.getInitialConfiguration(x);

  int NCH = IC.NCH; //Mに依存 when M=1,NCH=2
  int c;
  double r = RND.Uniform();
  for (c=0; c<NCH; c++) {
    if ( r < IC.CH[c].PROB ) break;
  } //初期状態の選択 alg.xmlの<Site>から

  ScatteringChannel& CH = IC.CH[c];

  int out  = CH.OUT;
  int xout = CH.XOUT;

  if ( out == DIR::UNDEF ) return false; // ワームは生成されない
  // セグメントをワームテールで切断 "UNDEF=-1"
  Segment& S1 = S0.cut( Vorg , 0 );

  // 方向に応じて初期セグメントを選ぶ

  if ( out == UORD::DOWN ) { 
    W.setCurrentSegment(S0);
  } else {       // moving up initially
    W.setCurrentSegment(S1);
  }

  // 初期タイプの設定

  W.setXBEHIND(xout);

#if defined(CF) || defined(CK)
  tail_tau  =t;
  tail_site =s;
#endif

#ifdef MONITOR
  printf("\nSimulation::PlaceWorm> End.\n");
#endif
  return true; // ワームが生成された
}

//======================================================================

double Simulation::MoveHead() {

#ifdef MONITOR
  printf("\n#####################################################\n");
  printf("Simulation::MoveHead> Start.\n");
  dump();
  printf("--------------------------------\n");
#endif

  double len=0.0;
  double len0=0.0;
  double len1=0.0;
  EndOfCycle = false;

  while(true) {
    //     W.getV();
    //    cout <<"getuord() = " << W.getUORD() << endl;
    if(W.getUORD()){

      len0=DOWN_ONESTEP();
      len += len0;
      len1-= len0;
      //             printf("-- DOWN -- len =%.4f len1=%.4f\n",len,len1);
      //     cout << "-- DOWN -- "<<len <<" "<<len1<<endl;
      //     len += DOWN_ONESTEP();
    }else{
      len0=UP_ONESTEP();
      len += len0;
      len1+= len0;
      //            printf("--  UP  -- len =%.4f len1=%.4f\n",len,len1);
      //     cout << "--  UP  -- "<<len <<" "<<len1<<endl;
      //     len += UP_ONESTEP();
    }
    if ( EndOfCycle ) break;
  }
#ifdef DEB
  if ( ! W.atOrigin() ) {  
    fprintf(FERR,"\nMoveHead> ### Error.\n");
    fprintf(FERR," ... Hasn't come back to the origin.\n");  
    dump();
    exit(0);
  }
#endif

  W.remove();

#ifdef MONITOR
  printf("\nSimulation::MoveHead> End.\n");
#endif
  return len;

}

//======================================================================
double Simulation::UP_ONESTEP(){
  double len = 0.0;
  double RHO = 0.0;

  Segment& c_S   = W.getCurrentSegment();
  Site&   c_Site = (* c_S.getONSITE() );
 
  Interaction** CI = c_Site.getCI();
  int NCI = c_Site.getNCI(); // the number of Interaction:

  Vertex& c_V = W.getCurrentVertex();
  Vertex& n_V = W.getNextVertex();
 
  double c_time;
  double n_time;
  int xinc = W.getXBEHIND();
  UniformInterval* UI = new UniformInterval[NCI];
  
  Ring<RegVInfo> RingRVI ;
  RingRVI.ROOT.V_x=&V4REF;
  // edit Suzuki
  // 逆温度に対する規格化
#ifdef NORM
  double n_Vtime =(n_V.isTerminal())? 1.0 : n_V.time();
#else
  double n_Vtime =(n_V.isTerminal())? LAT.BETA : n_V.time();
#endif
  // end edit Suzuki
  if(c_V.isTerminal()){
    c_time  = 0.0;
    n_time  = n_Vtime ;
    for(int iCI = 0; iCI < NCI; iCI++){ 
      UI[iCI].init(CI[iCI],xinc);
      for(int i_body = 0; i_body <UI[iCI].nbody ; i_body++){
	if( &( (*CI[iCI]).site( i_body ) ) == &c_Site ) {
	  UI[iCI].inc = 2 * i_body ;
	  UI[iCI].n_S[i_body] = &c_S ;
	}else{
	  UI[iCI].n_S[i_body] = &( (*CI[iCI]).site( i_body ).head() );
	  Vertex& near_V = (*UI[iCI].n_S[i_body]).top();
	  double V_time = near_V.time();
	  if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	    RegVInfo& RVI = TheRVIPool.pop();
	    RVI.setRVI( &near_V, iCI, i_body, V_time);
	    RingRVI.add_tail(RVI);
	  }
	}
      }
      UI[iCI].setx();
      UI[iCI].setVIC();

      //  edit Suzuki
      // 逆温度に対する規格化
#ifdef NORM
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO * LAT.BETA ; }
#else
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO; }
#endif
      //  end edit Suzuki
    }
  }else{
    c_time  = c_V.time();
    n_time  = n_Vtime-c_time;
    for(int iCI = 0; iCI < NCI; iCI++){
      UI[iCI].init(CI[iCI],xinc);
      for(int i_body = 0; i_body <UI[iCI].nbody ; i_body++){
	if( &( (*CI[iCI]).site( i_body ) ) == &c_Site ) {
	  UI[iCI].inc = 2 * i_body ;
	  UI[iCI].n_S[i_body] = &c_S ;
	}else{
	  if( CI[iCI] ==  c_V.getONINTERACTION()  ){
	    UI[iCI].n_S[i_body] = &( c_V.S( 2*i_body + 1 ) );
	  }else{
	    UI[iCI].n_S[i_body] = &( (*CI[iCI]).site( i_body ).findS( c_time ) );
	  }
	  Vertex& near_V = (*UI[iCI].n_S[i_body]).top();
	  double V_time = near_V.time() - c_time;
	  if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	    RegVInfo& RVI = TheRVIPool.pop();
	    RVI.setRVI( &near_V, iCI, i_body, V_time);
	    RingRVI.add_tail(RVI);
	  }
	}
      }
      UI[iCI].setx();
      UI[iCI].setVIC(); //どのVICかセットする

      //  edit Suzuki
      // 逆温度に対する規格化
#ifdef NORM
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO * LAT.BETA ; }
#else
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO; }
#endif
      //  end edit Suzuki
    } 
  }
  double near_tau = 0.0;
  double try_tau  = 0.0;
  int out, xout;

  while(true){

    Ring<RegVInfo>::iterator it = RingRVI.sort_min();
  
    near_tau =(RingRVI.empty())? n_time - len : (*it).V_time - len ;

    if( RHO > EPS ){
      try_tau = RND.Exp()/ RHO ;
      if( try_tau < near_tau ){
	len += try_tau;
	double RHORND = RHO * RND.Uniform();
		
	/*
	  double dRHO;
	  while(true){
	  dRHO = (UI[i_UI].DefinedVIC)?(*UI[i_UI].VIC).dRHO : 0.0 ;
	  if( RHORND > dRHO ){
	  RHORND -= dRHO;
	  ++i_UI;
	  }else{
	  break;
	  }
	  }
	*/
	int i_UI=0;
	int s_UI=0;
	double iRHO = 0.0;
	double sRHO = 0.0;
	do{
	  sRHO = iRHO;
	  if(UI[i_UI].DefinedVIC ){
	    //  edit Suzuki
	    // 逆温度に対する規格化
#ifdef NORM
	    iRHO += (*UI[i_UI].VIC).dRHO * LAT.BETA ;
#else
	    iRHO += (*UI[i_UI].VIC).dRHO;
#endif
	    // end edit Suzuki
	    s_UI = i_UI;
	  }
	  ++i_UI;
	}while( RHORND >= iRHO );
	RHORND -= sRHO;
	//  edit Suzuki
	// 逆温度に対する規格化
#ifdef NORM
	ScatteringChannel& CH = (*UI[s_UI].VIC).getScatteringChannel(RHORND / LAT.BETA);

#else
	ScatteringChannel& CH = (*UI[s_UI].VIC).getScatteringChannel(RHORND);
#endif
	// これは必要？
	// end edit Suzuki
	
	double setVtime =  c_time + len ;
	Vertex& setV = TheVertexPool.pop();
	setV.BareVertex::init(UI[s_UI].I_n, setVtime, (*UI[s_UI].VP));
	for(int ibody = 0; ibody<UI[s_UI].nbody; ++ibody){ (* UI[s_UI].n_S[ibody]).cut(setV, ibody);}
	(* UI[s_UI].I_n).add_tail(setV);

	setV.S(UI[s_UI].inc).setX(xinc);
	if(! ( c_V.isTerminal() || c_V.isKink() ) ){ c_V.erase(); }
	W.setCurrentVertex( setV );
	W.setCurrentSegment( setV.S( CH.OUT ) );
	W.setXBEHIND( CH.XOUT );
	
	break;
      }
    }
    if(RingRVI.empty()){
      W.setCurrentVertex(n_V);
      int inc  = n_V.which(c_S);
      if ( n_V.isTerminal() ) {
	out = 1-inc;
	xout = xinc;
      } else{
	VertexInitialConfiguration& IC  = n_V.getInitialConfiguration( inc , xinc );
	ScatteringChannel& CH = IC.getScatteringChannel();
	out = CH.OUT;
	xout = CH.XOUT;
      }

      c_S.setX( xinc );
      if(! ( c_V.isTerminal() || c_V.isKink() ) ){ c_V.erase(); }
      if ( out == DIR::UNDEF ) {
	EndOfCycle = true;
      } else {
	W.setCurrentSegment( n_V.S( out ) );
	W.setXBEHIND( xout );
      }
      len += near_tau;
      break;
    }
    
    //UPDATE RegVI and UI
    Vertex* changeV = (*it).V_x;
    do{
      RegVInfo& RVImin = (*it);
      // +++ edit sakakura +++ //
#ifdef NORM
      if(UI[RVImin.i_UI].DefinedVIC){RHO -= (*(UI[RVImin.i_UI].VIC)).dRHO * LAT.BETA; }
#else
      if(UI[RVImin.i_UI].DefinedVIC){RHO -= (*(UI[RVImin.i_UI].VIC)).dRHO; }
#endif
      // +++ edit sakakura +++ //
      //
      //
      // +++ edit sakakura +++ //
      //    if(UI[RVImin.i_UI].DefinedVIC){RHO -= (*(UI[RVImin.i_UI].VIC)).dRHO; }
      UI[RVImin.i_UI].n_S[RVImin.i_body] = &( (* UI[RVImin.i_UI].n_S[RVImin.i_body] ).next());//
      UI[RVImin.i_UI].setx(RVImin.i_body);
      UI[RVImin.i_UI].setVIC();
      // +++ edit sakakura +++ //
#ifdef NORM
      if(UI[RVImin.i_UI].DefinedVIC){RHO += (*(UI[RVImin.i_UI].VIC)).dRHO * LAT.BETA; }
#else 
      if(UI[RVImin.i_UI].DefinedVIC){RHO += (*(UI[RVImin.i_UI].VIC)).dRHO; }
#endif 
      // +++ edit sakakura +++ //
      //kota RHOの処理で逆温度規格化しなくてよい？

      //add RegVI
      Vertex& near_V = (*UI[RVImin.i_UI].n_S[RVImin.i_body]).top();//
      double V_time = near_V.time() - c_time;//
      if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	RegVInfo& RVI = TheRVIPool.pop();
	RVI.setRVI( &near_V, RVImin.i_UI, RVImin.i_body, V_time);
	RingRVI.add_tail(RVI);
      }
      ++it;
      RVImin.erase();
    }while(changeV == (*it).V_x);
    
    it = RingRVI.head();
    len += near_tau;
    
  }

  while(!RingRVI.empty()){
    RegVInfo& DRVI = RingRVI.first();
    DRVI.remove();
    TheRVIPool.push(DRVI);
  };
  
  delete [] UI;
  return len;
}
//======================================================================

//======================================================================
double Simulation::DOWN_ONESTEP(){
  double len = 0.0;
  double RHO = 0.0;

  Segment& c_S   = W.getCurrentSegment();
  Site&   c_Site = (* c_S.getONSITE() );
  
  Interaction** CI = c_Site.getCI(); //  Interaction* CI [n_int]; which belong to curSite
  int NCI = c_Site.getNCI(); // the number of Interaction:

  Vertex& c_V = W.getCurrentVertex();
  Vertex& n_V = W.getNextVertex();
  
  double c_time ;
  double n_time ;
  int xinc = W.getXBEHIND();
  UniformInterval* UI = new UniformInterval[NCI];

  Ring<RegVInfo> RingRVI ;
  RingRVI.ROOT.V_x=&V4REF;
  double n_Vtime =(n_V.isTerminal())? 0.0 : n_V.time();

  if(c_V.isTerminal()){
    //  +++ add sakakura +++
#ifdef NORM
    c_time = 1.0; //?
    n_time = 1.0 - n_Vtime ;
#else
    c_time  = LAT.BETA;
    n_time  = LAT.BETA - n_Vtime ;
#endif
    //  +++ add sakakura +++
    //    n_time  = LAT.BETA - n_Vtime ;
    for(int iCI = 0; iCI < NCI; iCI++){ 
      UI[iCI].init(CI[iCI],xinc);
      for(int i_body = 0; i_body <UI[iCI].nbody ; i_body++){
	if( &( (*CI[iCI]).site( i_body ) ) == &c_Site ) {
	  UI[iCI].inc = 2 * i_body + 1;
	  UI[iCI].n_S[i_body] = &c_S ;
	}else{
	  UI[iCI].n_S[i_body] = &( (*CI[iCI]).site( i_body ).tail() );//
	  Vertex& near_V = (*UI[iCI].n_S[i_body]).bottom();
	  double V_time = c_time - near_V.time();
	  if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	    RegVInfo& RVI = TheRVIPool.pop();
	    RVI.setRVI( &near_V, iCI, i_body, V_time);
	    RingRVI.add_tail(RVI);
	  }
	}
      }
      UI[iCI].setx();
      UI[iCI].setVIC();
      //  edit Suzuki
      // 逆温度に対する規格化
#ifdef NORM
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO * LAT.BETA ; }
#else
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO; }
#endif
      //  end edit Suzuki
    }
    
  }else{
    c_time  = c_V.time();
    n_time  = c_time - n_Vtime;
    for(int iCI = 0; iCI < NCI; iCI++){
      UI[iCI].init(CI[iCI],xinc);
      for(int i_body = 0; i_body <UI[iCI].nbody ; i_body++){
	if( &( (*CI[iCI]).site( i_body ) ) == &c_Site ) {
	  UI[iCI].inc = 2 * i_body + 1 ;//
	  UI[iCI].n_S[i_body] = &c_S ;
	}else{
	  if( CI[iCI] ==  c_V.getONINTERACTION()  ){
	    UI[iCI].n_S[i_body] = &( c_V.S( 2*i_body ) );//
	  }else{
	    UI[iCI].n_S[i_body] = &( (*CI[iCI]).site( i_body ).findS( c_time ) );
	  }
	  Vertex& near_V = (*UI[iCI].n_S[i_body]).bottom();//
	  double V_time = c_time - near_V.time();//
	  if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	    RegVInfo& RVI = TheRVIPool.pop();
	    RVI.setRVI( &near_V, iCI, i_body, V_time);
	    RingRVI.add_tail(RVI);
	  }
	}
      }
      UI[iCI].setx();
      UI[iCI].setVIC();
      //  edit Suzuki
      // 逆温度に対する規格化
#ifdef NORM
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO * LAT.BETA ; }
#else
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO; }
#endif 
      //  end edit Suzuki

    } 
  }
  double near_tau = 0.0;
  double try_tau  = 0.0;
  int out, xout;

  while(true){

    Ring<RegVInfo>::iterator it = RingRVI.sort_min();

    near_tau =(RingRVI.empty())? n_time - len : (*it).V_time - len ;

    if( RHO > EPS ){
      try_tau = RND.Exp()/ RHO ;
      if( try_tau < near_tau ){
	len += try_tau;
	double RHORND = RHO * RND.Uniform();
	/*
	  int i_UI=0;
	  double dRHO;
	  while(true){
	  dRHO = (UI[i_UI].DefinedVIC)?(*UI[i_UI].VIC).dRHO : 0.0 ;
	  if( RHORND > dRHO ){
	  RHORND -= dRHO;
	  ++i_UI;
	  }else{
	  break;
	  }
	  }
	*/

	int i_UI=0;
	int s_UI=0;
	double iRHO = 0.0;
	double sRHO = 0.0;
	do{
	  sRHO = iRHO;
	  if(UI[i_UI].DefinedVIC ){
	    //  edit Suzuki
	    // 逆温度に対する規格化
#ifdef NORM
	    iRHO += (*UI[i_UI].VIC).dRHO * LAT.BETA ;
#else
	    iRHO += (*UI[i_UI].VIC).dRHO;
#endif
	    // end edit Suzuki
	    s_UI = i_UI;
	  }
	  ++i_UI;
	}while( RHORND >= iRHO );
	RHORND -= sRHO;
	//  edit Suzuki
	// 逆温度に対する規格化
#ifdef NORM
	ScatteringChannel& CH = (*UI[s_UI].VIC).getScatteringChannel(RHORND / LAT.BETA);
#else
	ScatteringChannel& CH = (*UI[s_UI].VIC).getScatteringChannel(RHORND);
#endif
	//
	// end edit Suzuki

	double setVtime =  c_time - len ;
	Vertex& setV = TheVertexPool.pop();
	setV.BareVertex::init(UI[s_UI].I_n, setVtime, (*UI[s_UI].VP));
	for(int ibody = 0; ibody<UI[s_UI].nbody; ++ibody){ (* UI[s_UI].n_S[ibody]).cut(setV, ibody);}
	(* UI[s_UI].I_n).add_tail(setV);

	setV.S(UI[s_UI].inc).setX(xinc);
	if(! ( c_V.isTerminal() || c_V.isKink() ) ){ c_V.erase(); }
	W.setCurrentVertex( setV );
	W.setCurrentSegment( setV.S( CH.OUT ) );
	W.setXBEHIND( CH.XOUT );
	
	break;
      }
    }
    if(RingRVI.empty()){
      W.setCurrentVertex(n_V);
      int inc  = n_V.which(c_S);
      if ( n_V.isTerminal() ) {
	out = 1-inc;
	xout = xinc;
      } else {

	VertexInitialConfiguration& IC  = n_V.getInitialConfiguration( inc , xinc );
	ScatteringChannel& CH = IC.getScatteringChannel();
	out = CH.OUT;
	xout = CH.XOUT;
      }
      
      c_S.setX( xinc );
      if(! ( c_V.isTerminal() || c_V.isKink() ) ){ c_V.erase(); }
      if ( out == DIR::UNDEF ) {
	EndOfCycle = true;
      } else {
	W.setCurrentSegment( n_V.S( out ) );
	W.setXBEHIND( xout );
      }
      len += near_tau;
      break;
    }
    
    //UPDATE RegVI and UI
    Vertex* changeV = (*it).V_x;
    do{
      RegVInfo& RVImin = (*it);
      // +++ edit sakakura +++ //
#ifdef NORM 
      if(UI[(*it).i_UI].DefinedVIC){RHO -= (*(UI[(*it).i_UI].VIC)).dRHO * LAT.BETA; }
#else 
      if(UI[(*it).i_UI].DefinedVIC){RHO -= (*(UI[(*it).i_UI].VIC)).dRHO; }
#endif 
      // +++ edit sakakura +++ //
      UI[(*it).i_UI].n_S[(*it).i_body] = &( (* UI[(*it).i_UI].n_S[(*it).i_body] ).prev());//
      UI[(*it).i_UI].setx((*it).i_body);
      UI[(*it).i_UI].setVIC();
      // +++ edit sakakura +++ //
#ifdef NORM 
      if(UI[(*it).i_UI].DefinedVIC){RHO += (*(UI[(*it).i_UI].VIC)).dRHO * LAT.BETA; }
#else 
      if(UI[(*it).i_UI].DefinedVIC){RHO += (*(UI[(*it).i_UI].VIC)).dRHO; }
#endif 
      // +++ edit sakakura +++ //
      
      //add RegVI
      Vertex& near_V = (*UI[(*it).i_UI].n_S[(*it).i_body]).bottom();//
      double V_time = c_time - near_V.time();//
      if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	RegVInfo& RVI = TheRVIPool.pop();	
	RVI.setRVI( &near_V, (*it).i_UI, (*it).i_body, V_time);
	RingRVI.add_tail(RVI);
      }
      ++it;
      RVImin.erase();
    }while(changeV == (*it).V_x);
        
    it = RingRVI.head();
    len += near_tau;
        
  }

  while(!RingRVI.empty()){
    RegVInfo& DRVI = RingRVI.first();
    DRVI.remove();
    TheRVIPool.push(DRVI);
  };

  delete [] UI;
  return len;
}

//CF sweep here======================================================================
#if defined(CF) || defined(CK)

void Simulation::Sweep_CF() {
  //---CF
#ifdef CF
  reset4CF();
#endif
#ifdef CK
  reset4CK();
#endif
  //---
  double len = 0.0;
  for (ICYC=0; ICYC<P.NCYC; ICYC++) { 
    len += Cycle_CF(); 
  }
  AMP = len/(double)P.NCYC;

  //----CF
#ifdef CF
  for(int icf=0; icf< MSR.NCF; icf++)
    for(int it=0 ; it < MSR.Ntau; it++)
      cfsamp[icf][it] = ((double)counter4CF[icf][it])/((double)P.NCYC);
#endif
#ifdef CK
  for(int ick=0; ick< MSR.NCK; ick++)
    for(int it=0 ; it < MSR.Ntau; it++)
      cksamp[ick][it] = ((double)counter4CKC[ick][it])/((double)P.NCYC);
#endif
  //----

#ifdef GRAPHIC
  g_draw();
#endif

};

//======================================================================

double Simulation::Cycle_CF() {

  if ( ! PlaceWorm() ) return 0.0;
  return MoveHead_CF(); 

};

//======================================================================

double Simulation::MoveHead_CF() {

#ifdef MONITOR
  printf("\n#####################################################\n");
  printf("Simulation::MoveHead> Start.\n");
  dump();
  printf("--------------------------------\n");
#endif

  double len=0.0;
  double len0=0.0;
  double len1=0.0;
  EndOfCycle = false;
  ////koko1

  while(true) {
    //     W.getV();
    //    cout <<"getuord() = " << W.getUORD() << endl;
    if(W.getUORD()){

      len0=DOWN_ONESTEP_CF();
      len += len0;
      len1-= len0;
      //             printf("-- DOWN -- len =%.4f len1=%.4f\n",len,len1);
      //     cout << "-- DOWN -- "<<len <<" "<<len1<<endl;
      //     len += DOWN_ONESTEP();
    }else{
      len0=UP_ONESTEP_CF();
      len += len0;
      len1+= len0;
      //            printf("--  UP  -- len =%.4f len1=%.4f\n",len,len1);
      //     cout << "--  UP  -- "<<len <<" "<<len1<<endl;
      //     len += UP_ONESTEP();
    }
    if ( EndOfCycle ) break;
  }


  ////koko2
  //cout<<"len = "<<len<<endl;
#ifdef DEB
  if ( ! W.atOrigin() ) {  
    fprintf(FERR,"\nMoveHead> ### Error.\n");
    fprintf(FERR," ... Hasn't come back to the origin.\n");  
    dump();
    exit(0);
  }
#endif

  W.remove();

#ifdef MONITOR
printf("\nSimulation::MoveHead> End.\n");
#endif

  return len;

};

//======================================================================
double Simulation::UP_ONESTEP_CF(){
  double len = 0.0;
  double RHO = 0.0;

  Segment& c_S   = W.getCurrentSegment();
  Site&   c_Site = (* c_S.getONSITE() );
 
  Interaction** CI = c_Site.getCI();
  int NCI = c_Site.getNCI(); // the number of Interaction:

  Vertex& c_V = W.getCurrentVertex();
  Vertex& n_V = W.getNextVertex();
 
  double c_time;
  double n_time;
  int xinc = W.getXBEHIND();
  UniformInterval* UI = new UniformInterval[NCI];
  
  Ring<RegVInfo> RingRVI ;
  RingRVI.ROOT.V_x=&V4REF;
  // edit Suzuki
  // 逆温度に対する規格化
#ifdef NORM
  double n_Vtime =(n_V.isTerminal())? 1.0 : n_V.time();
#else
  double n_Vtime =(n_V.isTerminal())? LAT.BETA : n_V.time();
#endif
  // end edit Suzuki
  if(c_V.isTerminal()){
    c_time  = 0.0;
    n_time  = n_Vtime ;
    for(int iCI = 0; iCI < NCI; iCI++){ 
      UI[iCI].init(CI[iCI],xinc);
      for(int i_body = 0; i_body <UI[iCI].nbody ; i_body++){
	if( &( (*CI[iCI]).site( i_body ) ) == &c_Site ) {
	  UI[iCI].inc = 2 * i_body ;
	  UI[iCI].n_S[i_body] = &c_S ;
	}else{
	  UI[iCI].n_S[i_body] = &( (*CI[iCI]).site( i_body ).head() );
	  Vertex& near_V = (*UI[iCI].n_S[i_body]).top();
	  double V_time = near_V.time();
	  if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	    RegVInfo& RVI = TheRVIPool.pop();
	    RVI.setRVI( &near_V, iCI, i_body, V_time);
	    RingRVI.add_tail(RVI);
	  }
	}
      }
      UI[iCI].setx();
      UI[iCI].setVIC();

      //  edit Suzuki
      // 逆温度に対する規格化
#ifdef NORM
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO * LAT.BETA ; }
#else
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO; }
#endif
      //  end edit Suzuki
    }
  }else{
    c_time  = c_V.time();
    n_time  = n_Vtime-c_time;
    for(int iCI = 0; iCI < NCI; iCI++){
      UI[iCI].init(CI[iCI],xinc);
      for(int i_body = 0; i_body <UI[iCI].nbody ; i_body++){
	if( &( (*CI[iCI]).site( i_body ) ) == &c_Site ) {
	  UI[iCI].inc = 2 * i_body ;
	  UI[iCI].n_S[i_body] = &c_S ;
	}else{
	  if( CI[iCI] ==  c_V.getONINTERACTION()  ){
	    UI[iCI].n_S[i_body] = &( c_V.S( 2*i_body + 1 ) );
	  }else{
	    UI[iCI].n_S[i_body] = &( (*CI[iCI]).site( i_body ).findS( c_time ) );
	  }
	  Vertex& near_V = (*UI[iCI].n_S[i_body]).top();
	  double V_time = near_V.time() - c_time;
	  if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	    RegVInfo& RVI = TheRVIPool.pop();
	    RVI.setRVI( &near_V, iCI, i_body, V_time);
	    RingRVI.add_tail(RVI);
	  }
	}
      }
      UI[iCI].setx();
      UI[iCI].setVIC(); //どのVICかセットする

      //  edit Suzuki
      // 逆温度に対する規格化
#ifdef NORM
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO * LAT.BETA ; }
#else
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO; }
#endif
      //  end edit Suzuki
    } 
  }
  double near_tau = 0.0;
  double try_tau  = 0.0;
  int out, xout;

  while(true){

    Ring<RegVInfo>::iterator it = RingRVI.sort_min();
  
    near_tau =(RingRVI.empty())? n_time - len : (*it).V_time - len ;

    if( RHO > EPS ){
      try_tau = RND.Exp()/ RHO ;
      if( try_tau < near_tau ){
	len += try_tau;
	double RHORND = RHO * RND.Uniform();
		
	/*
	  double dRHO;
	  while(true){
	  dRHO = (UI[i_UI].DefinedVIC)?(*UI[i_UI].VIC).dRHO : 0.0 ;
	  if( RHORND > dRHO ){
	  RHORND -= dRHO;
	  ++i_UI;
	  }else{
	  break;
	  }
	  }
	*/
	int i_UI=0;
	int s_UI=0;
	double iRHO = 0.0;
	double sRHO = 0.0;
	do{
	  sRHO = iRHO;
	  if(UI[i_UI].DefinedVIC ){
	    //  edit Suzuki
	    // 逆温度に対する規格化
#ifdef NORM
	    iRHO += (*UI[i_UI].VIC).dRHO * LAT.BETA ;
#else
	    iRHO += (*UI[i_UI].VIC).dRHO;
#endif
	    // end edit Suzuki
	    s_UI = i_UI;
	  }
	  ++i_UI;
	}while( RHORND >= iRHO );
	RHORND -= sRHO;
	//  edit Suzuki
	// 逆温度に対する規格化
#ifdef NORM
	ScatteringChannel& CH = (*UI[s_UI].VIC).getScatteringChannel(RHORND / LAT.BETA);

#else
	ScatteringChannel& CH = (*UI[s_UI].VIC).getScatteringChannel(RHORND);
#endif
	// これは必要？
	// end edit Suzuki
	
	double setVtime =  c_time + len ;
	Vertex& setV = TheVertexPool.pop();
	setV.BareVertex::init(UI[s_UI].I_n, setVtime, (*UI[s_UI].VP));
	for(int ibody = 0; ibody<UI[s_UI].nbody; ++ibody){ (* UI[s_UI].n_S[ibody]).cut(setV, ibody);}
	(* UI[s_UI].I_n).add_tail(setV);

	setV.S(UI[s_UI].inc).setX(xinc);
	if(! ( c_V.isTerminal() || c_V.isKink() ) ){ c_V.erase(); }
	W.setCurrentVertex( setV );
	W.setCurrentSegment( setV.S( CH.OUT ) );
	W.setXBEHIND( CH.XOUT );
	
	break;
      }
    }
    if(RingRVI.empty()){
      W.setCurrentVertex(n_V);
      int inc  = n_V.which(c_S);
      if ( n_V.isTerminal() ) {
	out = 1-inc;
	xout = xinc;
      } else{
	VertexInitialConfiguration& IC  = n_V.getInitialConfiguration( inc , xinc );
	ScatteringChannel& CH = IC.getScatteringChannel();
	out = CH.OUT;
	xout = CH.XOUT;
      }

      c_S.setX( xinc );
      if(! ( c_V.isTerminal() || c_V.isKink() ) ){ c_V.erase(); }
      if ( out == DIR::UNDEF ) {
	EndOfCycle = true;
      } else {
	W.setCurrentSegment( n_V.S( out ) );
	W.setXBEHIND( xout );
      }
      len += near_tau;
      break;
    }
    
    //UPDATE RegVI and UI
    Vertex* changeV = (*it).V_x;
    do{
      RegVInfo& RVImin = (*it);
      // +++ edit sakakura +++ //
#ifdef NORM
      if(UI[RVImin.i_UI].DefinedVIC){RHO -= (*(UI[RVImin.i_UI].VIC)).dRHO * LAT.BETA; }
#else
      if(UI[RVImin.i_UI].DefinedVIC){RHO -= (*(UI[RVImin.i_UI].VIC)).dRHO; }
#endif
      // +++ edit sakakura +++ //
      //
      //
      // +++ edit sakakura +++ //
      //    if(UI[RVImin.i_UI].DefinedVIC){RHO -= (*(UI[RVImin.i_UI].VIC)).dRHO; }
      UI[RVImin.i_UI].n_S[RVImin.i_body] = &( (* UI[RVImin.i_UI].n_S[RVImin.i_body] ).next());//
      UI[RVImin.i_UI].setx(RVImin.i_body);
      UI[RVImin.i_UI].setVIC();
      // +++ edit sakakura +++ //
#ifdef NORM
      if(UI[RVImin.i_UI].DefinedVIC){RHO += (*(UI[RVImin.i_UI].VIC)).dRHO * LAT.BETA; }
#else 
      if(UI[RVImin.i_UI].DefinedVIC){RHO += (*(UI[RVImin.i_UI].VIC)).dRHO; }
#endif 
      // +++ edit sakakura +++ //
      //kota RHOの処理で逆温度規格化しなくてよい？

      //add RegVI
      Vertex& near_V = (*UI[RVImin.i_UI].n_S[RVImin.i_body]).top();//
      double V_time = near_V.time() - c_time;//
      if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	RegVInfo& RVI = TheRVIPool.pop();
	RVI.setRVI( &near_V, RVImin.i_UI, RVImin.i_body, V_time);
	RingRVI.add_tail(RVI);
      }
      ++it;
      RVImin.erase();
    }while(changeV == (*it).V_x);
    
    it = RingRVI.head();
    len += near_tau;
    
  }

  while(!RingRVI.empty()){
    RegVInfo& DRVI = RingRVI.first();
    DRVI.remove();
    TheRVIPool.push(DRVI);
  };
  
  delete [] UI;

#ifdef CF
  count4CF( MSR.ICF[ tail_site ][ c_Site.id() -1 ] , c_time+len, c_time );
#endif
#ifdef CK
  count4CK( (c_Site.id() -1 - tail_site + LAT.NSITE)%LAT.NSITE , c_time+len, c_time );
#endif
  return len;
}
//======================================================================

//======================================================================
double Simulation::DOWN_ONESTEP_CF(){
  double len = 0.0;
  double RHO = 0.0;

  Segment& c_S   = W.getCurrentSegment();
  Site&   c_Site = (* c_S.getONSITE() );
  
  Interaction** CI = c_Site.getCI(); //  Interaction* CI [n_int]; which belong to curSite
  int NCI = c_Site.getNCI(); // the number of Interaction:

  Vertex& c_V = W.getCurrentVertex();
  Vertex& n_V = W.getNextVertex();
  
  double c_time ;
  double n_time ;
  int xinc = W.getXBEHIND();
  UniformInterval* UI = new UniformInterval[NCI];

  Ring<RegVInfo> RingRVI ;
  RingRVI.ROOT.V_x=&V4REF;
  double n_Vtime =(n_V.isTerminal())? 0.0 : n_V.time();

  if(c_V.isTerminal()){
    //  +++ add sakakura +++
#ifdef NORM
    c_time = 1.0; //?
    n_time = 1.0 - n_Vtime ;
#else
    c_time  = LAT.BETA;
    n_time  = LAT.BETA - n_Vtime ;
#endif
    //  +++ add sakakura +++
    //    n_time  = LAT.BETA - n_Vtime ;
    for(int iCI = 0; iCI < NCI; iCI++){ 
      UI[iCI].init(CI[iCI],xinc);
      for(int i_body = 0; i_body <UI[iCI].nbody ; i_body++){
	if( &( (*CI[iCI]).site( i_body ) ) == &c_Site ) {
	  UI[iCI].inc = 2 * i_body + 1;
	  UI[iCI].n_S[i_body] = &c_S ;
	}else{
	  UI[iCI].n_S[i_body] = &( (*CI[iCI]).site( i_body ).tail() );//
	  Vertex& near_V = (*UI[iCI].n_S[i_body]).bottom();
	  double V_time = c_time - near_V.time();
	  if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	    RegVInfo& RVI = TheRVIPool.pop();
	    RVI.setRVI( &near_V, iCI, i_body, V_time);
	    RingRVI.add_tail(RVI);
	  }
	}
      }
      UI[iCI].setx();
      UI[iCI].setVIC();
      //  edit Suzuki
      // 逆温度に対する規格化
#ifdef NORM
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO * LAT.BETA ; }
#else
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO; }
#endif
      //  end edit Suzuki
    }
    
  }else{
    c_time  = c_V.time();
    n_time  = c_time - n_Vtime;
    for(int iCI = 0; iCI < NCI; iCI++){
      UI[iCI].init(CI[iCI],xinc);
      for(int i_body = 0; i_body <UI[iCI].nbody ; i_body++){
	if( &( (*CI[iCI]).site( i_body ) ) == &c_Site ) {
	  UI[iCI].inc = 2 * i_body + 1 ;//
	  UI[iCI].n_S[i_body] = &c_S ;
	}else{
	  if( CI[iCI] ==  c_V.getONINTERACTION()  ){
	    UI[iCI].n_S[i_body] = &( c_V.S( 2*i_body ) );//
	  }else{
	    UI[iCI].n_S[i_body] = &( (*CI[iCI]).site( i_body ).findS( c_time ) );
	  }
	  Vertex& near_V = (*UI[iCI].n_S[i_body]).bottom();//
	  double V_time = c_time - near_V.time();//
	  if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	    RegVInfo& RVI = TheRVIPool.pop();
	    RVI.setRVI( &near_V, iCI, i_body, V_time);
	    RingRVI.add_tail(RVI);
	  }
	}
      }
      UI[iCI].setx();
      UI[iCI].setVIC();
      //  edit Suzuki
      // 逆温度に対する規格化
#ifdef NORM
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO * LAT.BETA ; }
#else
      if(UI[iCI].DefinedVIC){RHO += (*(UI[iCI].VIC)).dRHO; }
#endif 
      //  end edit Suzuki

    } 
  }
  double near_tau = 0.0;
  double try_tau  = 0.0;
  int out, xout;

  while(true){

    Ring<RegVInfo>::iterator it = RingRVI.sort_min();

    near_tau =(RingRVI.empty())? n_time - len : (*it).V_time - len ;

    if( RHO > EPS ){
      try_tau = RND.Exp()/ RHO ;
      if( try_tau < near_tau ){
	len += try_tau;
	double RHORND = RHO * RND.Uniform();
	/*
	  int i_UI=0;
	  double dRHO;
	  while(true){
	  dRHO = (UI[i_UI].DefinedVIC)?(*UI[i_UI].VIC).dRHO : 0.0 ;
	  if( RHORND > dRHO ){
	  RHORND -= dRHO;
	  ++i_UI;
	  }else{
	  break;
	  }
	  }
	*/

	int i_UI=0;
	int s_UI=0;
	double iRHO = 0.0;
	double sRHO = 0.0;
	do{
	  sRHO = iRHO;
	  if(UI[i_UI].DefinedVIC ){
	    //  edit Suzuki
	    // 逆温度に対する規格化
#ifdef NORM
	    iRHO += (*UI[i_UI].VIC).dRHO * LAT.BETA ;
#else
	    iRHO += (*UI[i_UI].VIC).dRHO;
#endif
	    // end edit Suzuki
	    s_UI = i_UI;
	  }
	  ++i_UI;
	}while( RHORND >= iRHO );
	RHORND -= sRHO;
	//  edit Suzuki
	// 逆温度に対する規格化
#ifdef NORM
	ScatteringChannel& CH = (*UI[s_UI].VIC).getScatteringChannel(RHORND / LAT.BETA);
#else
	ScatteringChannel& CH = (*UI[s_UI].VIC).getScatteringChannel(RHORND);
#endif
	//
	// end edit Suzuki

	double setVtime =  c_time - len ;
	Vertex& setV = TheVertexPool.pop();
	setV.BareVertex::init(UI[s_UI].I_n, setVtime, (*UI[s_UI].VP));
	for(int ibody = 0; ibody<UI[s_UI].nbody; ++ibody){ (* UI[s_UI].n_S[ibody]).cut(setV, ibody);}
	(* UI[s_UI].I_n).add_tail(setV);

	setV.S(UI[s_UI].inc).setX(xinc);
	if(! ( c_V.isTerminal() || c_V.isKink() ) ){ c_V.erase(); }
	W.setCurrentVertex( setV );
	W.setCurrentSegment( setV.S( CH.OUT ) );
	W.setXBEHIND( CH.XOUT );
	
	break;
      }
    }
    if(RingRVI.empty()){
      W.setCurrentVertex(n_V);
      int inc  = n_V.which(c_S);
      if ( n_V.isTerminal() ) {
	out = 1-inc;
	xout = xinc;
      } else {

	VertexInitialConfiguration& IC  = n_V.getInitialConfiguration( inc , xinc );
	ScatteringChannel& CH = IC.getScatteringChannel();
	out = CH.OUT;
	xout = CH.XOUT;
      }
      
      c_S.setX( xinc );
      if(! ( c_V.isTerminal() || c_V.isKink() ) ){ c_V.erase(); }
      if ( out == DIR::UNDEF ) {
	EndOfCycle = true;
      } else {
	W.setCurrentSegment( n_V.S( out ) );
	W.setXBEHIND( xout );
      }
      len += near_tau;
      break;
    }
    
    //UPDATE RegVI and UI
    Vertex* changeV = (*it).V_x;
    do{
      RegVInfo& RVImin = (*it);
      // +++ edit sakakura +++ //
#ifdef NORM 
      if(UI[(*it).i_UI].DefinedVIC){RHO -= (*(UI[(*it).i_UI].VIC)).dRHO * LAT.BETA; }
#else 
      if(UI[(*it).i_UI].DefinedVIC){RHO -= (*(UI[(*it).i_UI].VIC)).dRHO; }
#endif 
      // +++ edit sakakura +++ //
      UI[(*it).i_UI].n_S[(*it).i_body] = &( (* UI[(*it).i_UI].n_S[(*it).i_body] ).prev());//
      UI[(*it).i_UI].setx((*it).i_body);
      UI[(*it).i_UI].setVIC();
      // +++ edit sakakura +++ //
#ifdef NORM 
      if(UI[(*it).i_UI].DefinedVIC){RHO += (*(UI[(*it).i_UI].VIC)).dRHO * LAT.BETA; }
#else 
      if(UI[(*it).i_UI].DefinedVIC){RHO += (*(UI[(*it).i_UI].VIC)).dRHO; }
#endif 
      // +++ edit sakakura +++ //
      
      //add RegVI
      Vertex& near_V = (*UI[(*it).i_UI].n_S[(*it).i_body]).bottom();//
      double V_time = c_time - near_V.time();//
      if( (! near_V.isTerminal()) && ( V_time < n_time )&&( &near_V != &n_V) ){
	RegVInfo& RVI = TheRVIPool.pop();	
	RVI.setRVI( &near_V, (*it).i_UI, (*it).i_body, V_time);
	RingRVI.add_tail(RVI);
      }
      ++it;
      RVImin.erase();
    }while(changeV == (*it).V_x);
        
    it = RingRVI.head();
    len += near_tau;
        
  }

  while(!RingRVI.empty()){
    RegVInfo& DRVI = RingRVI.first();
    DRVI.remove();
    TheRVIPool.push(DRVI);
  };

  delete [] UI;

#ifdef CF
  count4CF( MSR.ICF[ tail_site ][ c_Site.id() -1 ] , c_time, c_time - len );
#endif
#ifdef CK
  count4CK( (c_Site.id() -1 - tail_site + LAT.NSITE)%LAT.NSITE , c_time, c_time - len );
#endif

  return len;
}
#endif 
//end CF
//======================================================================

inline void Simulation::dump(const char* s) {

  printf("\n>>> %s (iset=%d imcse=%d imcsd=%d imcs=%d icyc= %d)",
	 s, ISET, IMCSE, IMCSD, IMCS, ICYC
	 );

  LAT.dump();

  W.dump();

}

//======================================================================

inline void Simulation::dump() { dump(""); }

//======================================================================

inline void Simulation::Check() {

#ifndef DEBUG
  return;
#endif

  bool ERROR = false ;

  for (int b=0; b<LAT.NINT; b++) {
    Interaction& I = LAT.I(b);
    if ( I.empty() ) continue;
    Interaction::iterator p(I);
    while ( ! (++p).atOrigin() ) {
      ERROR = ERROR || p->check();
    }
  }

  for (int s=0; s<LAT.NSITE; s++) {
    Site& S = LAT.S(s);
    Site::iterator p(S);
    while ( ! (++p).atOrigin() ) {
      ERROR = ERROR || p->check();
    }
  }

  if ( ERROR ) {
    dump();
    exit(0);
  }

}
//======================================================================
// == edit sakakura ==
//void Simulation::RepChange(int ix){
void Simulation::RepChangeF(Algorithm& algx){
  //  InteractionProperty& IP4=alg[ix].getInteractionProperty(alg[ix].NSTYPE-1);
  InteractionProperty& IP4=algx.getInteractionProperty(algx.NSTYPE-1);
  InteractionProperty& IP2 = ALG.getInteractionProperty(ALG.NSTYPE-1);

  for (int i=0; i<ALG.NXMAX; i++){
    int* x = new int[2];
    for (int j=0; j<ALG.NXMAX; j++){

      x[0]=i;x[1]=j;
      double d= IP4.VertexDensity(x);
      IP2.VertexDensity(x)=d;
      IP2.AverageInterval(x)=1.0/d;

    }
  }

  IP2.EBASE=  IP4.EBASE;
  MSR.EBASE = IP4.EBASE * LAT.NINT;
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VertexProperty& VP2 = ALG.getVertexProperty(1);
  //VertexProperty& VP4 = alg[ix].getVertexProperty(1);
  VertexProperty& VP4 = algx.getVertexProperty(1);

  for (int i=0; i<VP4.NIC; i++){

    VP2.IC[i].NLEG =VP4.IC[i].NLEG;
    VP2.IC[i].INC  =VP4.IC[i].INC ;
    VP2.IC[i].XINC =VP4.IC[i].XINC;
    VP2.IC[i].NCH  =VP4.IC[i].NCH ;
    VP2.IC[i].dRHO =VP4.IC[i].dRHO ;
    VP2.IC[i].State[0] = VP4.IC[i].State[0] ; 
    VP2.IC[i].State[1] = VP4.IC[i].State[1] ;
    VP2.IC[i].State[2] = VP4.IC[i].State[2] ;
    VP2.IC[i].State[3] = VP4.IC[i].State[3] ;

    for (int j = 0; j<VP2.IC[i].NCH; j++){
      VP2.IC[i].CH[j].OUT =VP4.IC[i].CH[j].OUT; 
      VP2.IC[i].CH[j].XOUT =VP4.IC[i].CH[j].XOUT; 
      VP2.IC[i].CH[j].PROB =VP4.IC[i].CH[j].PROB;
    }
  }

  ALG.initializer(VP4.NIC);
}
// == edit sakakura ==

void Simulation::RepChangeB(double hx){

  LAT.BETA=hx;
  //for (int ix=0; ix<pow(LAT.NSITE,LAT.D);ix++){
  for (int ix=0; ix<LAT.NSITE;ix++){
    Site& SX = LAT.S(ix);
    SX.setBeta(1.0);
    }
}
