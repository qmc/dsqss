
#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <iostream>
#include <stdio.h>
#include "name.h"
#include "io.h"
#include "array.h"
#include "link.hpp"
#include "xml.h"

//#include "objects.hpp"
//#######################################################################

class Algorithm;
class SiteProperty;
class SiteInitialConfiguration;
class InteractionProperty;
class VertexProperty;
class VertexInitialConfiguration;
class ScatteringChannel;

//######################################################################
//####  Object Property Declarations
//######################################################################

class SiteInitialConfiguration {
public:
  int State;
  int NCH;
  Array<ScatteringChannel> CH;
  SiteInitialConfiguration()  { 
    CH.setLabel("SiteInitialConfiguration::CH"); 
  };
  int id() { return State; };
  void initialize(XML::Block& X);
  void dump();
};

//#######################################################################

class SiteProperty {
public:
  int STYPE;
  int VTYPE;
  VertexProperty* _VP;
  int NX;
  int NIC;
  Array<SiteInitialConfiguration> IC;
  SiteProperty() {
    IC.setLabel("SiteProperty::IC");
  };

  int id() { return STYPE; };
  void initialize(XML::Block& X);
  void setVertexProperty(VertexProperty& VP) { _VP = &VP; };
  VertexProperty& getVertexProperty() { return *_VP; };
  SiteInitialConfiguration& getInitialConfiguration(int x) { 
    return IC[x];
  }
  void dump();
};

//#######################################################################

class ScatteringChannel {
public:
  int OUT;
  int XOUT;
  double PROB;

  void initialize(XML::Block& X);
  void dump();
};

//#######################################################################

class InteractionProperty {
public:
  int ITYPE;
  int VTYPE;
  int NBODY;
  int NXMAX;
  double EBASE;
  VertexProperty* _VP;
  Array<int> STYPE;
  Array<double> VertexDensity;
  Array<double> AverageInterval;

  InteractionProperty() {
    STYPE.setLabel("InteractionProperty::STYPE");
    VertexDensity.setLabel("InteractionProperty::VertexDensity");
    AverageInterval.setLabel("InteractionProperty::AverageInterval");
  };

  void setVertexProperty(VertexProperty& vp) { _VP = &vp; };
  VertexProperty& getVertexProperty() { return *_VP; };

  int id() { return ITYPE; };
  void initialize(XML::Block& X);
  void dump();
};

//#######################################################################

class VertexProperty {
private:
public:
  int VTYPE;
  int VCAT;
  int NBODY;
  int NLEG;
  int NST;
  int NIC;
  int NHTYPE;
  int NXMAX;
  // === edit sakakura  ===
  int MIC;
  // === edit sakakura  ===
  Array<int> STYPE;
  Array<VertexInitialConfiguration> IC;
  Array<VertexInitialConfiguration*> _IC;
  Array<int> StateCode;
  Array<int> SCNK; // State Code for NON-KINK

  VertexProperty() {
    STYPE.setLabel("VertexProperty::STYPE");
    IC.setLabel("VertexProperty::IC");
    _IC.setLabel("VertexProperty::_IC");
    StateCode.setLabel("VertexProperty::StateCode");
    SCNK.setLabel("VertexProperty::StateCode4NON-KINK");
  };

  int id() { return VTYPE; };

  // == edit sakakura ==
  void initialize(XML::Block& X,int i);
  //void initialize(XML::Block& X);
  // == edit sakakura ==

  int getSiteType(int out);
  bool isTerminal() { return VCAT == VCAT::TERM; };
  VertexInitialConfiguration& getIC(int st, int inc, int xinc);

  void dump();
};

//#######################################################################

class VertexInitialConfiguration {
private:


public:
  int ID;

  Array<ScatteringChannel> CH;

  VertexInitialConfiguration() {
    ID = -1;
    CH.setLabel("VertexInitialConfiguration::CH");
  };

  ~VertexInitialConfiguration() { 
    if ( State != 0 ) delete [] State;
  };

  int NLEG;
  int* State;
  int INC;
  int XINC;
  int NCH;

  // -- edit sakakura --
  int MCH; //NCHの最大値
  // -- edit sakakura --
  //
  double dRHO;// [vertex density] * [1-p1]
  void setID(int i) { ID = i; };
  int id() { return ID; };
  void initialize(XML::Block& X);

  // edit sakakura
  void initialize();
  // edit sakakura
  //
  ScatteringChannel& getScatteringChannel();
  ScatteringChannel& getScatteringChannel(double drho);
  void dump();

};

//int VertexInitialConfiguration::LastID = 0;

//#######################################################################

class Algorithm {

private:

  XML::Block X;
  Array<SiteProperty> SPROP;
  Array<InteractionProperty> IPROP;
  Array<VertexProperty> VPROP;

  // --edit sakakrua --
  int ix;
  // --edit sakakrua --
public:

  Algorithm(char* FNAME);
  
  // --edit sakakrua --
  Algorithm();//no
  // --edit sakakrua --

  ~Algorithm();

  int NSTYPE;   // Number of site types
  int NITYPE;   // Number of interaction types
  int NVTYPE;   // Number of bond types
  int NXMAX;    // Maximum number of segment states
  double WDIAG; // Artifitial weight attached to the diagonal state
  // Used in relating the worm-head mean-path to the susceptibility

  void read();
  void initialize();

  // --edit sakakrua --
  int MXNIC1;    // Maximum number of Vertex InitialConfigureations
  void initializer(int i );
  void set_i(int n) { ix=n; }
  void set_alg(int n) ;
  // --edit sakakrua --

  SiteProperty&        getSiteProperty(int i)        { return SPROP[i]; };//return val[i]
  InteractionProperty& getInteractionProperty(int i) { return IPROP[i]; };
  VertexProperty&      getVertexProperty(int i)      { return VPROP[i]; };

  void dump() { X.dump(); };
  //
  //  void get_spid(){return id();};

};

//######################################################################
//###########  Member Functions  #######################################
//######################################################################
//
// ============= edit sakakura ============

inline Algorithm::Algorithm(char* FNAME) {
  SPROP.setLabel("Algorithm::SPROP");
  IPROP.setLabel("Algorithm::IPROP");
  VPROP.setLabel("Algorithm::VPROP");

  if (DEBUG) printf("\nAlgorithm::Algorithm> Start.\n");
  X.initialize( FNAME ); //pass2 //xml::Block X

  // == edit sakakura ==
  int ixx =  X["General"]["NXMAX" ].getInteger(); 
  MXNIC1=24*(ixx-1)*(ixx-1);
  // == edit sakakura ==

  read();
  initialize();

  if (DEBUG) printf("Algorithm::Algorithm> End.\n");
}

//======================================================================
// ============= edit sakakura ============
inline Algorithm::Algorithm() {
}
// ============= edit sakakura ============
void Algorithm::set_alg(int ix) {
  SPROP.setLabel("Algorithm::SPROP");
  IPROP.setLabel("Algorithm::IPROP");
  VPROP.setLabel("Algorithm::VPROP");

  // == edit sakakura ==
  MXNIC1=0;
  // == edit sakakura ==
 
  //char* FNAME="algorithm.xml";

  stringstream ss;
  ss << ix << ".xml";
  //  char* xxx ;

  //  string zzz = ss.str();

  //
  //  if (DEBUG) printf("\nAlgorithm::Algorithm> Start.\n");
  //  X.initialize( FNAME ); 
  X.initialize(ss.str().c_str() ); 
  read();
  initialize();

  //  if (DEBUG) printf("Algorithm::Algorithm> End.\n");
}
// ============= edit sakakura ============


//======================================================================

void Algorithm::read() {

  //  if (DEBUG) printf("Algorithm::read> Start.\n");

  NSTYPE = X["General"]["NSTYPE"].getInteger(); //1 []operator
  NITYPE = X["General"]["NITYPE"].getInteger(); //1
  NVTYPE = X["General"]["NVTYPE"].getInteger(); //2
  NXMAX  = X["General"]["NXMAX" ].getInteger(); //2
  WDIAG  = X["General"]["WDIAG" ].getDouble(); // 0.25

  SPROP.init(1,NSTYPE); //1,SiteProperty
  IPROP.init(1,NITYPE); //1,InteractionProperty
  VPROP.init(1,NVTYPE); //2,VertexProperty

  for (int i=0; i<NITYPE; i++) { //1 IPROP(i) : return val[i]
    IPROP(i).NXMAX = NXMAX; 
  }

  for (int i=0; i<NVTYPE; i++) { //2
    VPROP(i).NXMAX  = NXMAX;
  }

  for (int i=0; i<X.NumberOfBlocks(); i++) { //X.numberofblocks=6
    XML::Block& B = X[i];
    const string& name = B.getName();

    if ( name == "Site" ) {
      int id = B["STYPE"].getInteger(); //id=0
      SPROP(id).initialize(B); //463行目 この時点でSitePropertyClass set完了
      //SPROP(0) は、SiteProperty& val[0]を返す。
    }

    if ( name == "Interaction" ) {
      int id = B["ITYPE"].getInteger();//id=0
      IPROP(id).initialize(B);//son
    }

    if ( name == "Vertex" ) { //riyo
      int id = B["VTYPE"].getInteger(); //id=0,1
      // -- edit sakakura --
      VPROP(id).initialize(B,MXNIC1);//mayu
      //    VPROP(id).initialize(B);//mayu
      // -- edit sakakura --

    }
  }
}


//= edit sakakura =======================================================
void Algorithm::initializer(int ix) {

  for (int i=1; i<2; i++) {
    VertexProperty& VP = VPROP(i);

    VP.NIC = ix;

    //   VP.initialize.StateCode.set_all(STATE::UNDEF);
    VP.StateCode.set_all(-1);
    VP._IC.set_all(0);
    VP.SCNK.set_all(STATE::UNDEF); //-1 //size=2,4
    //    cout <<"VP.NIC = " <<VP.NIC << endl;
    //    VP.dump();
    //    exit(31);

    int nst = 0;
 
    // --- edit sakakura ---
    for (int ic=0; ic<ix; ic++) {
      //    int ic = ix;
      //    cout <<"initializer ="<< ic << endl;
      // --- edit sakakura ---
      VertexInitialConfiguration& IC = VP.IC[ic]; //tomo

      int nl = VP.NLEG;
      int* x = IC.State;
      int inc = IC.INC;
      int xinc = IC.XINC;
      int st;
      bool isKink = false;
      int* xx = new int[nl/2];

      for(int il = 0; il < nl ; il+=2){
	if(x[il]!=x[il+1])isKink = true;
	else xx[il/2] = x[il];
      }
      //    printf ("ic=%d  ( %d %d %d %d ) vps=%d\n",ic,x[0],x[1],x[2],x[3],VP.StateCode(x));
      if ( VP.StateCode(x) == STATE::UNDEF ) { //operator return val(ID(x))
        VP.StateCode(x) = nst;
        st = nst;
        if(!isKink) VP.SCNK(xx)=nst;
        nst++;
      } else {
	st = VP.StateCode(x);
      }
      //      printf("ic=%d st=%d inc=%d xinc=%d\n",ic,st,inc,xinc) ;

      //      st = VP.StateCode(x);
      VP._IC( st, inc, xinc) = &IC; //st,inc,xincを利用してINDEX生成//rei
      delete [] xx;
    } //end VP.NICloop
    //   VP.dump();
  }
  // InteractionProperty の _VP を初期化
  for (int i=0; i<NITYPE; i++) {
    InteractionProperty& I = IPROP(i);
    I.setVertexProperty(VPROP(I.VTYPE));
  }

  return;

  // 散乱確率の積算と ScatteringChannel::PROB の再々定義
  // NOT KINK のバーテックスに対してはPROB = [vertexdensity]*[probablity]とする
  for (int i=1; i<2; i++) {
    VertexProperty& VP = VPROP(i);
    // --- edit sakakura ---

    for (int j=VP.NIC; j<VP.MIC; j++) {

      VertexInitialConfiguration& IC = VP.IC[j];

      bool isKink = false;
      int* x = new int[VP.NBODY]; //edit sakakura
      //      int x[VP.NBODY];
      for(int ileg = 0; ileg<IC.NLEG; ileg+=2){
	if(IC.State[ileg]!=IC.State[ileg+1]){isKink=true;}
	x[ileg/2]=IC.State[ileg];
      }
      double p = 0.0;

      bool isVertex = false;// false: kink or term , true: means interaction
      int i_type =0;
      for(i_type = 0; i_type<NITYPE; i_type++){
	if(IPROP[i_type].VTYPE==i){isVertex = true; break;}
      }// If VTYPE is for Vertex (not for term or tail), there should be InteractionProperty with density of vertex.

      if( isKink||(!isVertex )){

	IC.dRHO=1.0;
	for (int k=0; k<IC.NCH; k++) {
	  ScatteringChannel& SC = IC.CH[k];
	  p += SC.PROB;
	  SC.PROB = p;
	}
	
      }else{
	double Vdensity = IPROP[i_type].VertexDensity(x);
	int count_ch = 0;
	double dR = Vdensity;

	for (int k=0; k<IC.NCH; k++) {

	  ScatteringChannel& SC = IC.CH[k];
	  if(SC.OUT == IC.INC + 1 - 2*(IC.INC%2)){
	    dR = (1.0-SC.PROB)*Vdensity; 
	  }else{
	    p += (SC.PROB * Vdensity);
	    IC.CH[count_ch] = SC;
	    IC.CH[count_ch].PROB = p; 
	    count_ch++;
	  }

	}
	IC.dRHO = dR;
	IC.NCH  = count_ch;
      }
      if ( fabs(p-IC.dRHO) > EPS * 10.0 ) {
	printf("Algorithm::initialize> Error.\n");
	printf("  The sum of probabilities is not equal to dRHO.\n");
	VP.dump();
	IC.dump();
	exit(0);
      }
      
      if(IC.NCH != 0){ IC.CH[ IC.NCH - 1 ].PROB = IC.dRHO + EPS; 
      }else{ IC.dRHO = 0.0;}
      delete []x; 
    }
  }
  
  //  if (DEBUG) printf("Algorithm::initialize> End.\n");
}
  
//= edit sakakura =======================================================

//======================================================================

void Algorithm::initialize() {

  // if (DEBUG) printf("Algorithm::initialize> Start.\n");
  // Site::_VP の初期化

  for (int i=0; i<NSTYPE; i++) { //1
    int vt = SPROP(i).VTYPE;//2
    SPROP(i).setVertexProperty( VPROP(vt) ); //VPROP(vt) {return val}

    //kota! setVertexProperty(VertexProperty& VP) { _VP = &VP; };//
  }
  // 散乱確率の積算と ScatteringChannel::PROB の再定義 //1.0にする

  for (int i=0; i<NSTYPE; i++) {//1
    SiteProperty& SP = SPROP(i);
    //    SP.dump(); 定義前
    for (int j=0; j<SP.NIC; j++) { //NIC=NX=Numberofstates 2
      SiteInitialConfiguration& IC = SP.IC[j]; //[]operator return val[i]
      double p = 0.0;

      for (int k=0; k<IC.NCH; k++) { //numberofchannnels
        ScatteringChannel& SC = IC.CH[k];
        p += SC.PROB;
        SC.PROB = p;
      }
      if ( fabs(p-1.0) > EPS ) { //EPS=1.d-14
        printf("Algorithm::initialize> Error.\n");
        printf("  The sum of probabilities is not equal to 1.\n");
        SP.dump();
        IC.dump();
        exit(0);
      }
      IC.CH[ IC.NCH - 1 ].PROB = 1.0;
    }
  }

  // バーテックスのサイトタイプ情報の初期化

  for (int i=0; i<NSTYPE; i++) {
    SiteProperty& S = SPROP(i);
    VertexProperty& V = VPROP(S.VTYPE);

    int NLEG = 2;
    for (int l=0; l<NLEG; l++) {
      V.STYPE[l] = S.STYPE;
    }
  }
  for (int i=0; i<NITYPE; i++) {
    InteractionProperty& I = IPROP(i);
    VertexProperty& V = VPROP(I.VTYPE);
    int NLEG = 2 * I.NBODY;
    for (int l=0; l<NLEG; l++) {
      V.STYPE[l] = I.STYPE[l/2];
    }
  }

  // バーテックスプロパティーの _IC の初期化
  for (int i=0; i<NVTYPE; i++) {
    VertexProperty& VP = VPROP(i);

    int nst = 0;

    for (int ic=0; ic<VP.NIC; ic++) {

      VertexInitialConfiguration& IC = VP.IC[ic]; //tomo

      int nl = VP.NLEG;
      int* x = IC.State;
      int inc = IC.INC;
      int xinc = IC.XINC;
      int st;
      bool isKink = false;
      int* xx = new int[nl/2];

      for(int il = 0; il < nl ; il+=2){
	if(x[il]!=x[il+1])isKink = true;
	else xx[il/2] = x[il];
      }

      if ( VP.StateCode(x) == STATE::UNDEF ) { //operator return val(ID(x))
        VP.StateCode(x) = nst;
        st = nst;
	if(!isKink) VP.SCNK(xx)=nst;
        nst++;
      } else {
        st = VP.StateCode(x);
      }

      //    printf("ic=%d st=%d inc=%d xinc=%d\n",ic,st,inc,xinc) ;
      VP._IC( st, inc, xinc) = &IC; //st,inc,xincを利用してINDEX生成//rei
      delete [] xx;
    } //end VP.NICloop

  }
  
  // InteractionProperty の _VP を初期化
  for (int i=0; i<NITYPE; i++) {
    InteractionProperty& I = IPROP(i);
    I.setVertexProperty(VPROP(I.VTYPE));
  }
 
  // 散乱確率の積算と ScatteringChannel::PROB の再々定義
  // NOT KINK のバーテックスに対してはPROB = [vertexdensity]*[probablity]とする
  for (int i=0; i<NVTYPE; i++) {
    VertexProperty& VP = VPROP(i);
    // --- edit sakakura ---

    for (int j=0; j<VP.NIC; j++) {

      VertexInitialConfiguration& IC = VP.IC[j];

      bool isKink = false;
      int* x = new int[VP.NBODY]; //edit sakakura
      //      int x[VP.NBODY];
      for(int ileg = 0; ileg<IC.NLEG; ileg+=2){
	if(IC.State[ileg]!=IC.State[ileg+1]){isKink=true;}
	x[ileg/2]=IC.State[ileg];
      }
      double p = 0.0;

      bool isVertex = false;// false: kink or term , true: means interaction
      int i_type =0;
      for(i_type = 0; i_type<NITYPE; i_type++){
	if(IPROP[i_type].VTYPE==i){isVertex = true; break;}
      }// If VTYPE is for Vertex (not for term or tail), there should be InteractionProperty with density of vertex.

      if( isKink||(!isVertex )){
	//      printf("i=%d j=%d iskink \n",i,j);
	IC.dRHO=1.0;
	for (int k=0; k<IC.NCH; k++) {
	  ScatteringChannel& SC = IC.CH[k];
	  p += SC.PROB;
	  SC.PROB = p;
	}
	
      }else{
	double Vdensity = IPROP[i_type].VertexDensity(x);
	int count_ch = 0;
	double dR = Vdensity;
	//       printf("i=%d j=%d dR=%f \n",i,j,dR);

	for (int k=0; k<IC.NCH; k++) {

	  ScatteringChannel& SC = IC.CH[k];
	  if(SC.OUT == IC.INC + 1 - 2*(IC.INC%2)){
	    dR = (1.0-SC.PROB)*Vdensity; 
	  }else{
	    p += (SC.PROB * Vdensity);
	    IC.CH[count_ch] = SC;
	    IC.CH[count_ch].PROB = p; 
	    count_ch++;
	  }

	}
	IC.dRHO = dR;
	IC.NCH  = count_ch;
      }
      if ( fabs(p-IC.dRHO) > EPS * 10.0 ) {
	printf("Algorithm::initialize> Error.\n");
	printf("  The sum of probabilities is not equal to dRHO.\n");
	VP.dump();
	IC.dump();
	exit(0);
      }
      
      if(IC.NCH != 0){ IC.CH[ IC.NCH - 1 ].PROB = IC.dRHO + EPS; 
      }else{ IC.dRHO = 0.0;}
      delete []x;
    }
  }
  
  //////=============================================///////

  //  if (DEBUG) printf("Algorithm::initialize> End.\n");
}

//======================================================================

inline Algorithm::~Algorithm() {
  //    printf("*** Destroying Algorithm\n");
}

//######################################################################
//######################################################################
//######################################################################

void SiteProperty::initialize(XML::Block& X) { //<Site> 

  STYPE = X["STYPE"].getInteger();
  NX    = X["NumberOfStates"].getInteger();
  VTYPE = X["VertexTypeOfSource"].getInteger();
  NIC   = NX;

  IC.init(1,NIC);

  for (int i=0; i<X.NumberOfBlocks(); i++) {
    XML::Block& B = X[i];
    if ( B.getName() == "InitialConfiguration" ) {
      int st = B["State"].getInteger();
      IC[st].initialize(B);
    }
  }
}

//======================================================================

void SiteProperty::dump() {
  printf("\n");
  printf("<SiteProperty>\n");
  printf("  STYPE= %d\n", STYPE);
  printf("  NX=    %d\n", NX);
  printf("  VTYPE= %d\n", VTYPE);
  printf("  NIC=   %d\n", NIC);
  for (int i=0; i<NIC; i++) IC[i].dump(); 
  printf("</SiteProperty>\n");
}

//#######################################################################

void SiteInitialConfiguration::initialize(XML::Block& X) {

  //printf("SiteInitialConfiguration::initialize> Start.\n");

  State = X["State"].getInteger();
  NCH   = X["NumberOfChannels"].getInteger();
  CH.init("SCH",1,NCH);

  int ch = 0;
  for (int i=0; i<X.NumberOfBlocks(); i++) {
    XML::Block& B = X[i];
    if ( B.getName() == "Channel" ) CH[ch++].initialize(B);
  }

  if ( ch != NCH ) {
    printf("SiteInitialConfiguration::initialize> Error.\n");
    printf("  The actual number of channels (=%d)\n", ch);
    printf("  does not agree with NCH (=%d)\n", NCH);
  }

}

//======================================================================

void SiteInitialConfiguration::dump() {
  printf("\n");
  printf("<SiteInitialConfiguration>\n");
  printf("  State= %d\n", State);
  printf("  NCH=   %d\n", NCH);
  for (int i=0; i<NCH; i++) CH[i].dump(); 
  printf("</SiteInitialConfiguration>\n");
}

//#######################################################################

void ScatteringChannel::initialize(XML::Block& X) {
  OUT  = X.getInteger(0);
  XOUT = X.getInteger(1);
  PROB = X.getDouble(2);
}

//======================================================================

void ScatteringChannel::dump() {
  //  printf("    CH= %2d", Code);
  printf("  OUT= %2d", OUT);
  printf(", XOUT= %2d", XOUT);
  printf(", P= %8.3f", PROB);
  printf("\n");
}

//#######################################################################

void InteractionProperty::initialize(XML::Block& X) { //<interaction> son
  ITYPE = X["ITYPE"].getInteger();
  VTYPE = X["VTYPE"].getInteger();
  NBODY = X["NBODY"].getInteger();
  EBASE = X["EBASE"].getDouble();
  STYPE.init( 1, NBODY);
  STYPE.set_all(STYPE::UNDEF);
  VertexDensity.init( NBODY, NXMAX, ARRAY::EOL);
  VertexDensity.set_all(0.0);
  AverageInterval.init( NBODY, NXMAX, ARRAY::EOL);
  AverageInterval.set_all(-1.0);
  
  for (int i=0; i<X.NumberOfBlocks(); i++) {
    XML::Block& B = X[i];
    if (B.getName() == "VertexDensity") {
      int*x = new int[NBODY]; //edit sakakura
      //    int x[NBODY];
      for (int ii=0; ii<NBODY; ii++) {
        x[ii] = B.getInteger(ii);
      }
      double d = B.getDouble(NBODY);
      // VertexDensity   val[n] = d    を格納(n:1-n)
      // AverageInterval val[n] = 1/d  を格納(n:1-n)
      VertexDensity(x) = d; //x[0],x[1]分が処理される return C* val[id]
      AverageInterval(x) = 1.0/d;
      delete []x;
    }
  }
  //  Numberofblocksは/Interaction の行数
  
}

//======================================================================

inline void InteractionProperty::dump() {
  printf("\n");
  printf("<InteractionProperty>\n");
  printf("  ITYPE= %2d", ITYPE);
  printf(", VTYPE= %2d", VTYPE);
  printf(", NBODY= %2d", NBODY);
  printf(", NXMAX= %2d", NXMAX);
  printf(", EBASE= %24.16f\n", EBASE);
  IndexSystem& I = VertexDensity.index_system();
  int* x = new int[I.dimension()];

  for (int i=0; i<I.size(); i++) {
    I.coord(i,x);
    printf("     (");
    for (int j=0; j<I.dimension(); j++) {
      printf(" %1d", x[j]);
    }
    printf(") --> %8.3f, %8.3f\n", 
           VertexDensity(x), AverageInterval(x)  ); //kota averagei=1/vd
  }
  printf("</InteractionProperty>\n");
  delete [] x;
}

//#######################################################################

void VertexProperty::initialize( XML::Block& X,int mx) {//<vertex> //mayu
  VTYPE  = X["VTYPE"].getInteger(); 
  VCAT   = X["VCATEGORY"].getInteger();
  NBODY  = X["NBODY"].getInteger();
  NIC    = X["NumberOfInitialConfigurations"].getInteger();
  NLEG = 2 * NBODY;
 
  STYPE.init(1,NLEG);//array(int) STYPE //size 2,4
  STYPE.set_all(STYPE::UNDEF); //val[i]に(size()個分)-1をセット

  StateCode.init(NLEG,NXMAX,ARRAY::EOL); //array(int) StateCode -1をセット
  StateCode.set_all(STATE::UNDEF); //2^2,4^2
  
  SCNK.init(NBODY,NXMAX,ARRAY::EOL); //array(int) scnk
  SCNK.set_all(STATE::UNDEF); //-1 //size=2,4

  NST = StateCode.size(); //4,16
  _IC.init(3,NST,NLEG,NXMAX); //index用の数値用意[4,2,2],[16,4,2]
  _IC.set_all(0); //size 16,128

  // === edit sakakura ===
  //  if ( RUNTYPE == 1 ) {
  int id = X["VTYPE"].getInteger(); //id=0,1

  if (id == 0){
    IC.init(1,NIC);
  }else{
    MIC = NIC;

    if (mx > NIC)MIC=mx;
    IC.init(1,MIC);
  }
  //   }
  // === edit sakakura ===

  int ic = 0;
  for (int i=0; i<X.NumberOfBlocks(); i++) { //2番目のtag
    XML::Block& B = X[i];
    if ( B.getName() == "InitialConfiguration" ) {
      IC[ic].setID(ic);
      IC[ic].NLEG = NLEG;
      IC[ic].initialize( B );
      ic++;
    }
  }

  // == edit sakakura ==
  if (id==1) {
    for ( int i=NIC; i<MIC; i++){
      IC[i].setID(i);
      IC[i].NLEG = NLEG;//NLEG;
      IC[i].initialize(); //mai
    }
  }
  // == edit sakakura ==

}

//======================================================================

inline int VertexProperty::getSiteType(int out) {
#ifdef DEB
  if ( out < 0 || out >= NLEG ) {
    printf("VertexProperty::getSiteType> Error.\n");
    printf("  The argument (= %d) is out of the bounds.\n", out);
    dump();
    exit(0);
  }
#endif
  return STYPE[out]; 
}

//======================================================================

inline VertexInitialConfiguration& VertexProperty::getIC( int st, int inc, int xinc) {
  return *_IC(st, inc, xinc);
}

//======================================================================

void VertexProperty::dump() {
  printf("\n");
  printf("<VertexProperty>\n");
  printf("  VTYPE  = %d\n", VTYPE);
  printf("  NBODY  = %d\n", NBODY);
  printf("  NLEG   = %d\n", NLEG);
  printf("  NST    = %d\n", NST);
  printf("  NIC    = %d\n", NIC);
  printf("  NXMAX  = %d\n", NXMAX);
  printf("  STYPE  = ");
  for (int i=0; i<NLEG; i++) printf(" %d", STYPE[i]);
  printf("\n");
  int* x = new int[NLEG];
  IndexSystem& I = StateCode.index_system();
  for (int i=0; i<I.size(); i++) {
    I.coord(i,x);
    int sc = StateCode[i];
    printf("  StateCode");
    for (int j=0; j<I.dimension(); j++) {
      printf("[%1d]", x[j]);
    }
    printf(" = %d\n", sc);
  }
  IndexSystem& I2 = SCNK.index_system();
  for (int i=0; i<I2.size(); i++) {
    I2.coord(i,x);
    int sc = SCNK[i];
    printf("  StateCode4NON-KINK");
    for (int j=0; j<I2.dimension(); j++) {
      printf("[%1d]", x[j]);
    }
    printf(" = %d\n", sc);
  }
  int y[3];
  IndexSystem& J = _IC.index_system();
  cout <<"J.size() = "<<J.size() << endl;
  for (int i=0; i<J.size(); i++) {
    J.coord(i,y);
    if ( _IC[i] != 0 ) {
      int icc = (*_IC[i]).id();
      printf("  (_IC[%1d][%1d][%1d])->id() = %d\n",
             y[0],y[1],y[2],icc);
    }
  }
  for (int i=0; i<NIC; i++) IC[i].dump();
  printf("</VertexProperty>\n");
  delete [] x;
}

//#######################################################################

void VertexInitialConfiguration::initialize(XML::Block& X) {//<initialconfig>

  //printf("VertexInitialConfiguration::initialize> Start.\n");

  if (NLEG == 0) {
    printf("VertexInitialConfiguration::read> Error.\n");
    printf(" ...  NLEG has not been set.\n");
    exit(0);
  }

  State = new int [NLEG];
  for (int i=0; i<NLEG; i++) {
    State[i] = X["State"].getInteger(i);
  }

  INC  = X["IncomingDirection"].getInteger();
  XINC = X["NewState"].getInteger();
  NCH  = X["NumberOfChannels"].getInteger();

  // -- edit sakakura --
  MCH  = 4;
  CH.init("VCH",1,MCH);
  //CH.init("VCH",1,NCH);//org
  // -- edit sakakura --

  int ch = 0;
  for (int i=0; i<X.NumberOfBlocks(); i++) {
    XML::Block& B = X[i];
    if (B.getName() == "Channel") CH[ch++].initialize(B);
  }

  if ( ch != NCH ) {
    printf("SiteInitialConfiguration::initialize> Error.\n");
    printf("  The actual number of channels (=%d)\n", ch);
    printf("  does not agree with NCH (=%d)\n", NCH);
  }
}

//======================================================================
//###edit sakakura######################################################mai

void VertexInitialConfiguration::initialize() {//<initialconfig>

  //printf("VertexInitialConfiguration::initialize> Start.\n");

  if (NLEG == 0) {
    printf("VertexInitialConfiguration::read> Error.\n");
    printf(" ...  NLEG has not been set.\n");
    exit(0);
  }

  State = new int [NLEG];

  NCH  = 4;  //X["NumberOfChannels"].getInteger();
  CH.init("VCH",1,NCH);

}

//====edit sakakura========================================================


inline ScatteringChannel& VertexInitialConfiguration::getScatteringChannel() {
  double p;
  if ( NCH == 1 ) return CH[0];
  int ch;
  double r = RND.Uniform();
  for ( ch=0; ch<NCH; ch++) {
    p = CH[ch].PROB;
    if (r < p) break;
  }
  if ( ch == NCH ) {
    printf("ERROR! at Usual Vertex (KINK) \n");
    dump();
    printf("VertexInitialConfiguration::getScatteringChannel> Error.\n");
    printf("  ... ch=%d exceeds the limit.\n",ch);
    exit(0);
  }
  return CH[ch];
}

//======================================================================
inline ScatteringChannel& VertexInitialConfiguration::getScatteringChannel(double drho) {
  if ( NCH == 1 ) return CH[0];
  int ch;
  double p;
  for ( ch=0; ch<NCH; ch++) {
    p = CH[ch].PROB;
    if (drho < p) break;
  }
  if ( ch == NCH ) {
    printf("ERROR! at place Vertex (non-KINK)\n");
    dump();
    printf("VertexInitialConfiguration::getScatteringChannel> Error.\n");
    printf("  ... ch=%d exceeds the limit.\n",ch);
    exit(0);
  }
  return CH[ch];
}

//======================================================================

void VertexInitialConfiguration::dump() {
  printf("\n");
  printf("<VertexInitialConfiguration> ID= %d\n", ID);
  printf("  NLEG= %d",NLEG);
  printf(", INC= %d",INC);
  printf(", XINC= %d",XINC);
  printf(", NCH= %d",NCH);
  printf(", dRHO= %f",dRHO);
  
  printf(", X= (");
  for (int i=0; i<NLEG; i++) printf(" %1d", State[i]);
  printf(")\n");
  for (int i=0; i<NCH; i++) CH[i].dump();
  printf("</VertexInitialConfiguration>\n");
}

#endif
