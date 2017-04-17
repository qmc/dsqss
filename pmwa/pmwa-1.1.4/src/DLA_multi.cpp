#include <DLA.h>
#include <time.h>

//#define check

int main(int argc, char **argv){


  Dla Sim(argc, argv);



  ////////////////////////////////////////////////////////////
  Sim.PMWA();
  ////////////////////////////////////////////////////////////


  return 0;

}



Dla::Dla( int NP, char **PLIST ){


  int p_num, my_rank;

  MPI_Init(&NP, &PLIST);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p_num);


  MPI_Barrier(MPI_COMM_WORLD); 
  allstart = MPI_Wtime();

  pstart=pend=tstart=tend=ostart=oend=0.0; 

  ReadBuffer( p_num, my_rank, NP, PLIST);


}

Dla::~Dla(){
 
  MPI_Barrier(MPI_COMM_WORLD);
  allend = MPI_Wtime();

  if(PR.nst==0){
    show_TIME();
    ftest.close();
    fclose(FOUT4SF);
  }


  MPI_Finalize();  

}


//*************************************************************************************************

double Dla::PMWA(){


  
  Lattice LT( latfile, &N, &PR );

  NameOutfiles();

  My_rdm MR(MC.seed);
  Probability P( &N, &sp, &PR );

  P.LookUpTable( &N, &sp );

  Configuration CS(&MC,&N,sp.nmax,&LT,&P,&PR,IMAX,WMAX,Eventfile_old,&MR,(bool)cb,outfile);
  Quantities QNT( &N, &MC, &sp, &LT, &PR, sfinpfile ); 

  //##########################################################################
  ////////////////// Determination of Ncyc ////////////////////////
  if(MC.nc<=1){

    MPI_Barrier(MPI_COMM_WORLD); 
    pstart = MPI_Wtime(); 
    
    CS.DeterminationNworm(MC.Ntest,&MR,&QNT);
 
    MPI_Barrier(MPI_COMM_WORLD);
    pend = MPI_Wtime();

  }
  if(PR.nst==0){
    ftest << "C This is PMWA ver.0" << endl << endl;
    LT.show_param(ftest);
    CS.dump(ftest);
#ifdef SF
    fprintf(FOUT4SF,"C This is PMWA ver.0\n\n");
    QNT.dump( FOUT4SF );
    fclose(FOUT4SF);
    FOUT4SF = fopen(sfoutfile,"a"); 
#endif
  }

  if(MC.nc==1){ 
    CS.Output(Eventfile,&MR); 
    return 1; 
  }

  cout<<"start thermal..."<<endl;
  ////////////////// Thermarization ////////////////////////
  if(MC.nc < 3){

    MPI_Barrier(MPI_COMM_WORLD); 
    tstart = MPI_Wtime(); 

    CS.sweep(MC.Nthermal,&MR, &QNT);

    MPI_Barrier(MPI_COMM_WORLD);
    tend = MPI_Wtime();

  }

  if(MC.nc==2){ 
    CS.Output(Eventfile,&MR); 
    return 2;
  }

  /////////////////////////////////////////////////////////
  ////////////////// Observation ////////////////////////
  /////////////////////////////////////////////////////////
  MPI_Barrier(MPI_COMM_WORLD); 
  ostart = MPI_Wtime(); 

  CS.measurement( &QNT, &MR );

  if( PR.nst==0 )  QNT.show( ftest, FOUT4SF );

  MPI_Barrier(MPI_COMM_WORLD);
  oend = MPI_Wtime();
  
  ////////////////// Output the Configuration ////////////////////////
  CS.Output(Eventfile,&MR);


  return 0;

  //////////////////////////////////////////////////////////////////////

}

void Dla::NameOutfiles(){

  sprintf(outfile,"%s.%d",outfile,PR.np);
  sprintf(sfoutfile,"%s.%d",sfoutfile,PR.np);

  
  if(PR.nst==0){
    ftest.open( outfile, std::ios::app );
    ftest << setprecision(16);
#ifdef SF
    FOUT4SF = fopen(sfoutfile,"a"); 
    //   cout<<"sf start"<<endl;
 #endif

  }
  

  #ifdef check
  show_MPIP();
  show_SP();
  #endif
  #ifdef SERIAL_S
  if(PR.Nxdiv!=1||PR.Nydiv!=1||PR.Nzdiv!=1){
    cout<< "Must be no spacially-decomposition" <<endl;
    exit(1);
  }
  #endif
  //##########################################################################


  sprintf(Eventfile,"evout_%s_rank%d.dat",outfile,PR.my_rank);
#ifdef Annealing
  sprintf(Eventfile_old,"evout_%s_rank%d.dat",outfile,PR.my_rank);
#else
  sprintf(Eventfile_old,"evout_%s_rank%d.dat",outfile,PR.my_rank);
#endif

  //##########################################################################

}

void Dla::ReadBuffer(int m_pnum, int m_myrank, int NP, char **PLIST){

  PR.p_num = m_pnum;
  PR.my_rank = m_myrank;

  init(); // == set default values ==

  ifstream ifs(PLIST[1]); 

  if (NP < 2)
    {
      printf("error on Open: inputfile\n");
      exit(1) ;                //　inputfileの指定なし
    }

  if (ifs.fail())
    {
      printf("no such File: %s\n",PLIST[1]);
      exit (1);                //ファイルIOエラー
    }

  double val;
  string buf0;
  string buf1;
  string buf2;
  string buf2a,buf2b;

  while(ifs && getline(ifs, buf0))
    {
      buf1 = trim(buf0);

      //===== コメントの処理 ====
      string::size_type index=buf1.find("#");

      if ( index == string::npos )
	{
          buf2=buf1;
	}
      else
	buf2=buf1.substr(0,index);

      if (buf2.empty()) continue;
      //===== コメントの処理 ====

      //===== 構文解析 ====
      string :: size_type index1=buf2.find("=");

      if ( index1 == string::npos )
	{
          printf("input struct miss 1: %s\n",buf0.c_str());
          exit(1);            //"="が未記述
	}
                               
      //    buf2a=ToLower(trim(buf2.substr(0,index1)));
      buf2a=trim(buf2.substr(0,index1));//"="の左側
      transform(buf2a.begin(), buf2a.end(), buf2a.begin(),ToLower());
                               
      buf2b=trim(buf2.substr(index1+1,buf2.length()));//"="の右側

      int i,nval;

      int ikey = -1;
      for (i=0; i<num_para; i++)
        {
	  if (buf2a == para[i]) ikey = i; //どのパラメータか？
	  continue;
        }
      if (ikey==-1)
        {
	  //printf("input struct miss 2:\n %s\n",buf0.c_str());
	  //exit(1);            //パラメータ名がただしくない
	  printf("Unknown parameter:\n %s\n",buf0.c_str());
	  continue;	//知らないパラメータ名はスキップ

        }

      index1=buf2b.find(" ");
      if ( index1 != string::npos )
        {
	  printf("input struct miss 3:\n %s\n",buf0.c_str());
	  exit(1);            //右端に複数の文字列が存在

        }

      if (ikey < 10 || ( ikey > 21 && ikey < NUMPARA) )
        {
	  if (numcheck(buf2b))
            {
	      printf("must be integer value  :\n %s\n",buf0.c_str());
	      exit(1); //整数指定必須
            }
          nval = atoi(buf2b.c_str());
	  //        cout << ikey<<nval <<endl;
        }


      if (ikey < 22 && ikey > 14) {
	val = atof(buf2b.c_str());
	//     cout <<ikey<<" "<<val << endl;
      }


      if (wcount[ikey]){
	printf("dual definition error   :\n %s\n",buf0.c_str());
	exit(1);//2重定義
      }

      wcount[ikey] = true;

      //    cout << ikey<<nval <<endl;
      set_value(ikey,nval,val,buf2b);

    }
  if (!wcount[0])  cout << "error: must set runtype" << endl;
        


  /////////////////////
  N.d = (N.y==1)? 1: (( N.z==1 )? 2: 3);

  MC.seed += PR.my_rank;

  sp.Eu = PR.V*sp.Ubb/2.0*sp.nmax*( sp.nmax - 1.0 );

  N.V = N.x*N.y*N.z;

  //  set_Parallel();


}

string Dla::trim(const string& s)
{
  string::const_iterator left = s.begin();
  string::const_iterator right = s.end();
  while(isspace(*left))
    {
      ++left;
    }
  if(left != right)
    {                            // 右端から非空白文字を検索
      do { --right; }
      while(isspace(*right));
      ++right;
    }
  return string(left, right);
}

void Dla::init(){

  
  ///// Parameter //////
  
  MC.runtype  = 3 ;   //0
  MC.Nbin  = 10;      //1
  MC.Nthermal = 1000 ;//2
  MC.Nstep  = 1000 ;  //3
  MC.Nsample  = 1000 ;//4
  MC.Ntest = 1000 ;   //5
  MC.seed  = 13 ;     //6
  MC.nc = 0 ;      //7

  IMAX = 100000 ; //8
  WMAX = 5000 ;   //9
  ///////////////////////////
  //10-12 .... input file
  ////////
  cb = 0 ;        //13
  ///// Hamiltonian //////
  sp.Htr = 0.2 ; //14
  sp.Ubb = 0.0 ; //15
  sp.Vb1 = 3.0 ; //16
  sp.Vb2 = 0.0 ; //17
  sp.tb = 1.0 ;  //18
  sp.mu =1.0 ;   //19
  sp.Vc = 0.0 ;  //20
  sp.seed = 1 ;  //21
  sp.nmax = 1 ;  //22

  strcpy(algfile,"algorithm.xml");
  strcpy(latfile,"lattice.xml");
  strcpy(sfinpfile,"sf.xml"); 
  strcpy(outfile,"res.dat");
  strcpy(sfoutfile,"sfout.dat");
  
  num_para = NUMPARA;

  for  (int i=0;i<num_para; i++) wcount[i]=false;
  
  para[0]  = "runtype";
  para[1]  = "nset";
  para[2]  = "nmcse";
  para[3]  = "nmcsd";
  para[4]  = "nmcs";
  para[5]  = "ntest";
  para[6]  = "seed";
  para[7]  = "nc";

  para[8]  = "nvermax";
  para[9]  = "nwormax";

  para[10] = "algfile";
  para[11] = "latfile";
  para[12] = "sfinpfile";
  para[13] = "outfile";
  para[14] = "sfoutfile";

  para[15] = "g";
  para[16] = "ubb";
  para[17] = "vb1";
  para[18] = "vb2";
  para[19] = "tb";
  para[20] = "mu";
  para[21] = "wc";
  
  para[22] = "cb";
  para[23] = "ps";
  para[24] = "nmax";


    
}

void Dla::set_value(int i,int nx,double val, string buf){
  switch(i)
    {
    case 0:
      MC.runtype=nx;
      break;
    case 1:
      MC.Nbin=nx;
      break;
    case 2:
      MC.Nthermal=nx;
      break;
    case 3:
      MC.Nstep=nx;
      break;
    case 4:
      MC.Nsample=nx;
      break;
    case 5:
      MC.Ntest=nx;
      break;
    case 6:
      MC.seed=nx;
      break;
    case 7:
      MC.nc=nx;
      break;
    case 8:
      IMAX=nx;
      break;
    case 9:
      WMAX=nx;
      break;
    case 10:
      strcpy(algfile,buf.c_str());
      break;
    case 11:
      strcpy(latfile,buf.c_str());
      break;
    case 12:
      strcpy(sfinpfile,buf.c_str());
      break;
    case 13:
      strcpy(outfile,buf.c_str());
      break;
    case 14:
      strcpy(sfoutfile,buf.c_str());
      break;
    case 15:
      sp.Htr=val;
      break;
    case 16:
      sp.Ubb=val;
      break;
    case 17:
      sp.Vb1=val;
      break;
    case 18:
      sp.Vb2=val;
      break;
    case 19:
      sp.tb=val;
      break;
    case 20:
      sp.mu=val;
      break;
    case 21:
      sp.Vc=val;
      break;
    case 22:
      cb=nx;
      break;
    case 23:
      sp.seed=nx;
      break;
    case 24:
      sp.nmax=nx;
      break;
    }
}

int  Dla::numcheck(const string& s)
{
  int i;
  for (i=0; i<s.length(); i++)
    {
      if (isdigit(s[i])==0)
        {
	  return 1;
        }
    }
  return 0;
}

void Dla::set_Parallel(){


  PR.B = N.B/(double)PR.Ntdiv; // beta for a domain.
  PR.oldB = N.oldB/(double)PR.Ntdiv; // for annealing.

  PR.Nsdiv = PR.Nxdiv*PR.Nydiv*PR.Nzdiv;// the number of spatial decompositions.

  PR.x = N.x/PR.Nxdiv;
  PR.y = N.y/PR.Nydiv;
  PR.z = N.z/PR.Nzdiv;

  PR.V = PR.x*PR.y*PR.z;

  PR.NtNs = PR.Ntdiv*PR.Nsdiv;//the number of decompositions.
  PR.Npara=PR.p_num/PR.NtNs;//the number of seeds for trivial parallelization.


  PR.nt = PR.my_rank%PR.Ntdiv; // the temporal domain number for the processor.
  PR.ns = (int)(PR.my_rank/PR.Ntdiv)%PR.Nsdiv;  // the spatial domain number for the processor.

  PR.nx = PR.ns%PR.Nxdiv;  // the x-directional domain number for the processor.
  PR.ny = (PR.ns/PR.Nxdiv)%PR.Nydiv;  // the y-directional domain number for the processor.
  PR.nz = PR.ns/(PR.Nxdiv*PR.Nydiv); // the z-directional domain number for the processor.

  PR.nst = PR.my_rank%PR.NtNs;  // the domain number for the processor.
  PR.np = PR.my_rank/(PR.NtNs);  // the seed number of the trivial parallelization for the processor.

  if(PR.Rpara>0){
    PR.nr = PR.np%PR.Rpara; // a random potential number (one of trivial parallelization)
    PR.nq = PR.np/PR.Rpara; // a seed number of the trivial parallelization for the random potential (one of trivial parallelization)
  }
  else{
    PR.nr = 0;
    PR.nq = PR.np;
  }

  //the coordinate is (nt,nx,ny,nz,np)

  PR.nst0 = PR.np*PR.NtNs;// nst=0 process number for the processor.
  PR.nt0 = PR.ns*PR.Ntdiv + PR.np*PR.NtNs; //nt=0 process number for the processor.
  PR.ns0 = PR.nt + PR.np*PR.NtNs; //ns=0 process number for the processor.
  PR.nx0 = PR.nt + ( PR.ny*PR.Nxdiv + PR.nz*PR.Nxdiv*PR.Nydiv )*PR.Ntdiv + PR.np*PR.NtNs; //nx=0 process number for the processor.

  PR.upper= (PR.nt + 1)%PR.Ntdiv + PR.ns*PR.Ntdiv + PR.np*PR.Ntdiv*PR.Nsdiv; //the upper process number for the temporal direction.
  PR.lower= (PR.nt - 1 + PR.Ntdiv)%PR.Ntdiv + PR.ns*PR.Ntdiv + PR.np*PR.Ntdiv*PR.Nsdiv; //the lower process number the temporal direction.

  PR.right[0] = PR.nt + ( (PR.nx+1)%PR.Nxdiv           + PR.ny*PR.Nxdiv + PR.nz*PR.Nxdiv*PR.Nydiv )*PR.Ntdiv + PR.np*PR.Ntdiv*PR.Nsdiv; //the right side process number for the x direction.
  PR.left[0]  = PR.nt + ( (PR.nx-1+PR.Nxdiv)%PR.Nxdiv + PR.ny*PR.Nxdiv + PR.nz*PR.Nxdiv*PR.Nydiv )*PR.Ntdiv + PR.np*PR.Ntdiv*PR.Nsdiv; //the left side process number for the x direction.

  PR.right[1] = PR.nt + ( PR.nx + ((PR.ny+1)%PR.Nydiv)*PR.Nxdiv           + PR.nz*PR.Nxdiv*PR.Nydiv )*PR.Ntdiv + PR.np*PR.Ntdiv*PR.Nsdiv; //the right side process number for the y direction.
  PR.left[1]  = PR.nt + ( PR.nx + ((PR.ny-1+PR.Nydiv)%PR.Nydiv)*PR.Nxdiv + PR.nz*PR.Nxdiv*PR.Nydiv )*PR.Ntdiv + PR.np*PR.Ntdiv*PR.Nsdiv; //the left side process number for the y direction.

  PR.right[2] = PR.nt + ( PR.nx + PR.ny*PR.Nxdiv + ((PR.nz+1)%PR.Nzdiv)*PR.Nxdiv*PR.Nydiv )*PR.Ntdiv + PR.np*PR.Ntdiv*PR.Nsdiv; //the right side process number for the z direction.
  PR.left[2]  = PR.nt + ( PR.nx + PR.ny*PR.Nxdiv + ((PR.nz-1+PR.Nzdiv)%PR.Nzdiv)*PR.Nxdiv*PR.Nydiv )*PR.Ntdiv + PR.np*PR.Ntdiv*PR.Nsdiv; //the left side process number for the z direction.


}
