#ifndef PARAMETER_H
#define PARAMETER_H

#include <stdio.h>
#include <string>

#include <string.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <cctype>
#include <cstdio>

#include <dsqss_read.hpp>

//######################################################################

class Parameter
{
private:
  int num_para ;//= 17;
  bool wcount[17];
  string para[17];
public:
  Parameter(int, char**);
  ~Parameter();
  void dump(FILE*);
  void dump();
  void openfile(){ FOUT = fopen(OUTFILE,"w"); };
  void closefile(){ fclose(FOUT); };

  void init();
  void set_value(int,int,double,string);

  int  numcheck(const string& );
  string trim(const string& );

  struct ToLower{char operator()(char c) { return tolower(c);}};
  //
  // total number of sweeps executed in the whole simulation is
  //    NMCSE + NSET * ( NMCSD + NMCS )

  // ++++ edit sakakura ++++
  int NREP;                // number of replica
  //  int MXVIC;               // max number of Vertexinitialconfiguration
  // ++++ edit sakakura ++++
  int RUNTYPE;             // runtype  number
  int NSET;                // number of sets
  int NMCSE;               // number of sweeps for determining NCYC
  int NMCSD;               // number of sweeps for decorrelating two subsequent sets
  int NMCS;                // number of sweeps for measurement in a set
  int SEED;                // seed for the random number generator
  int NSEGMAX;             // maximum allowed number of segments
  int NVERMAX;             // maximum allowed number of vertices

  //kota edit //
  double VF; //value of VF
  double DF; //delta F
  double VB; //delta B
  double DB; //delta B
  //kota edit //

  //  int NRVIMAX; // maximum allowed number of registered vertex informations (used at ONESTEP())
  char ALGFILE[24];        // algorithm file name
  char LATFILE[24];        // lattice file name
  char OUTFILE[24];        // output file name
  int NCYC;                // number of cycles in a sweep (not provided from the file)
  //FILE* FALG; // file handler for the algorithm file
  //FILE* FLAT; // file handler for the lattice file
  FILE* FOUT;              // file handler for the output file
  //  int num_para = 11;
  Input  inp;
};

//######################################################################

inline Parameter::Parameter(int NP, char** PLIST)
{
  if (DEBUG) printf("\nParameter::Parameter> Start.\n");
  double val;
  //  init();
  //  read();
  //  set();
  inp.init();

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
	  printf("input struct miss 2:\n %s\n",buf0.c_str());
	  exit(1);            //パラメータ名がただしくない

        }

      index1=buf2b.find(" ");
      if ( index1 != string::npos )
        {
	  printf("input struct miss 3:\n %s\n",buf0.c_str());
	  exit(1);            //右端に複数の文字列が存在

        }

      if (ikey < 9)
        {
	  if (numcheck(buf2b))
            {
	      printf("must be integer value  :\n %s\n",buf0.c_str());
	      exit(1); //整数指定必須
            }
          nval = atoi(buf2b.c_str());
#ifdef MULTI
          if(ikey == 8 ){
	    if (nval != N_PROC){
	      printf("NREP must be equal with NPROC  :\n %s\n",buf0.c_str());
	      exit(1); //レプリカ数=NPROC必須（要検討）
	    }
          }
#endif

	  //        cout << ikey<<nval <<endl;
        }


      if (ikey < 16 && ikey > 11) {
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
        
  dump();


#ifdef MULTI
  if (RUNTYPE==0)SEED += I_PROC;
  char* IPSTR = new char[8];
  sprintf( IPSTR, ".%03d", I_PROC);
  strcat ( OUTFILE , IPSTR);
  delete [] IPSTR;
#endif

  if (DEBUG) printf("Parameter::Parameter> End.\n");
  NCYC = 0;
}

//======================================================================

inline Parameter::~Parameter()
{
  //  printf("*** Destroying Parameter\n");
}

int  Parameter::numcheck(const string& s)
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

string Parameter::trim(const string& s)
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

void Parameter::init(){
  NREP  = 1 ;
  VF    = 0.0;
  DF    = 0.1;
  VB    = 0.0;
  DB    = 0.1;
  RUNTYPE  = 0 ;
  NSET  = 10 ;
  NMCSE = 1000 ;
  NMCSD = 1000 ;
  NMCS  = 1000 ;
  SEED  = 198212240 ;
  NSEGMAX = 10000 ;
  NVERMAX = 10000 ;
#ifdef MULTI
  // +++ edit sakakura +++
  stringstream ss;
  ss << I_PROC << ".xml";
  strcpy(ALGFILE,ss.str().c_str());
  //
  // strcpy(ALGFILE,"algorithm.xml");
  // +++ edit sakakura +++
#endif
 
  strcpy(ALGFILE,"algorithm.xml");
  strcpy(LATFILE,"lattice.xml");
  strcpy(OUTFILE,"res.dat");
  num_para = 17;

  for  (int i=0;i<num_para; i++) wcount[i]=false;
  
  para[0]  = "runtype";
  para[1]  = "nset";
  para[2]  = "nmcse";
  para[3]  = "nmcsd";
  para[4]  = "nmcs";
  para[5]  = "seed";
  para[6]  = "nsegmax";
  para[7]  = "nvermax";
  para[8]  = "nrep";
  para[9]  = "algfile";
  para[10]  = "latfile";
  para[11] = "outfile";
  para[12] = "vf";
  para[13] = "df";
  para[14] = "vb";
  para[15] = "db";
  //  para[15] = "db";
  //  para[16] = "db";
    
}

void Parameter::set_value(int i,int nx,double val, string buf){
  switch(i)
    {
    case 0:
      RUNTYPE=nx;
      break;
    case 1:
      NSET=nx;
      break;
    case 2:
      NMCSE=nx;;
      break;
    case 3:
      NMCSD=nx;
      break;
    case 4:
      NMCS=nx;
      break;
    case 5:
      SEED=nx;
      break;
    case 6:
      NSEGMAX=nx;
      break;
    case 7:
      NVERMAX=nx;
      break;
    case 8:
      NREP=nx;
      break;
    case 9:
      strcpy(ALGFILE,buf.c_str());
      break;
    case 10:
      strcpy(LATFILE,buf.c_str());
      break;
    case 11:
      strcpy(OUTFILE,buf.c_str());
      break;
    case 12:
      VF=val;
      break;
    case 13:
      DF=val;
      break;
    case 14:
      VB=val;
      break;
    case 15:
      DB=val;
      break;
    }
}

//======================================================================

inline void Parameter::dump()
{
  cout << endl<<"+++++++++ input data +++++++++" << endl;
  cout << "RUNTYPE = " << RUNTYPE << endl<<endl;
  cout << "NREP    = " << NREP    << endl;
  cout << "NSET    = " << NSET    << endl;
  cout << "NMCSD   = " << NMCSD   << endl;
  cout << "NMCS    = " << NMCS    << endl;
  cout << "SEED    = " << SEED    << endl;
  cout << "NVERMAX = " << NVERMAX << endl;
  cout << "NSEGMAX = " << NSEGMAX << endl<<endl;
  cout << "ALGFILE = " << ALGFILE << endl;
  cout << "LATFILE = " << LATFILE << endl;
  cout << "OUTFILE = " << OUTFILE << endl;
  cout << "+++++++++ input data +++++++++" << endl;
}
inline void Parameter::dump(FILE* F)
{
  fprintf(F,"P NSET    = %12d\n",NSET);
  fprintf(F,"P NMCSE   = %12d\n",NMCSE);
  fprintf(F,"P NMCSD   = %12d\n",NMCSD);
  fprintf(F,"P NMCS    = %12d\n",NMCS);
  fprintf(F,"P SEED    = %12d\n",SEED);
  fprintf(F,"P NSEGMAX = %12d\n",NSEGMAX);
  fprintf(F,"P NVERMAX = %12d\n",NVERMAX);
  fprintf(F,"P NCYC    = %12d\n",NCYC);
  fprintf(F,"P ALGFILE = %s\n",    ALGFILE);
  fprintf(F,"P LATFILE = %s\n",    LATFILE);
  fprintf(F,"P OUTFILE = %s\n",    OUTFILE);
}
#endif

//inline void  Parameter::numcheck(int s){cout << s << endl;;}
//int Parameter::numcheck(const string& s)
//{
//    int i;
//    for (i=0; i<s.length(); i++)
//    {
//        if (isdigit(s[i])==0)
//        {
//            return 1;
//        }
//    }
//    return 0;
//}
