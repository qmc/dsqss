//======================================================================

void Measurement::measure() {

  using namespace Specific;
  int NV = LAT.countVertices();

  ACC[NV1].accumulate( ((double)NV) );

  //  double *phase= new double [Dim]; 

  double phase[2] = { 1.0 , -1.0 };

  double MZU  = 0.0;
  double MZUA = 0.0;
  double MZS  = 0.0;
  double MZSA = 0.0;
  double MZQC = 0.0;
  double MZQS = 0.0;

#ifdef SF
    reset4SF();
#endif

  
  for (int s=0; s<LAT.NSITE; s++) {
    Site& SITE = LAT.S(s);
    int mt = SITE.getMTYPE();
    double ph = phase[mt];
    Segment& S0 = SITE.first();
    double mz0 = S0.X();
    Site::iterator p(SITE);
    double mza0 = 0.0;

#ifdef SF
    double bTri=0.0;
    double tau=0.0;
    int it=0;
#endif

    
    while ( ! (++p).atOrigin() ) {
      Segment& S = *p;
      double mz = S.X();

#ifdef SF
      double tTri=bTri+S.length();
      
      while(tau<tTri){

	for(int k=0;k<NKMAX;k++){
	  counter4SFC[k][it] += mz * COSrk[s][k];
	  counter4SFS[k][it] += mz * SINrk[s][k];
	}

	it++;
	tau = it*dtau;
	
      }

      bTri=tTri;
#endif
      
      mza0 += mz * S.length();
    }

    MZU  +=      mz0;
    MZUA +=      mza0;
    MZS  += ph * mz0;
    MZSA += ph * mza0;

    //#ifdef SPECIFICQ  
    double CQ=SITE.getCOSQ();
    double SQ=SITE.getSINQ();
    MZQC  += CQ * mz0;
    MZQS  += SQ * mz0;
    //#endif


  }

#ifdef SF
  
  for(int k=0;k<NSF;k++){
    for(int it=0; it<Ntau;it++){
      double SZKT=0.0;
      for(int tt=0; tt<Ntau1;tt++){

	SZKT += counter4SFC[k][tt]*counter4SFC[k][(tt+it)%Ntau1]
	  + counter4SFS[k][tt]*counter4SFS[k][(tt+it)%Ntau1];
      }
      sfsamp[k][it] = SZKT/(double)Ntau1;
    }
  }
  
  accumulate4SF( sfsamp );
  
#endif

  
  ACC[MZU1].accumulate(MZU);
  ACC[MZU2].accumulate(MZU * MZU);
  ACC[MZUA1].accumulate(MZUA);
  ACC[MZUA2].accumulate(MZUA * MZUA);
  ACC[MZS1].accumulate(MZS);
  ACC[MZS2].accumulate(MZS * MZS);
  ACC[MZSA1].accumulate(MZSA);
  ACC[MZSA2].accumulate(MZSA * MZSA);
  ACC[MZQ2].accumulate(MZQC * MZQC + MZQS * MZQS);

  double EBSAMP = - (double) NV ;

  for (int b=0; b<LAT.NINT; b++) {
    Interaction& I = LAT.I(b);
    InteractionProperty& IP = I.property();
    //VertexProperty& VP = IP.getVertexProperty();
    int NBODY = IP.NBODY;
    double* tau = new double[NBODY];
    int* x = new int[NBODY];
    Site::iterator* p = new Site::iterator[NBODY];
    
    for (int i=0; i<NBODY; i++) {
      Site& S = I.site(i);
      p[i].init(S);
      p[i]++;
      tau[i] = p[i]->topTime();
      x[i] = p[i]->X();
    }

    double t  = 0.0;
    int it;

    // edit Suzuki
    // 逆温度に対する規格化
#ifdef NORM
    while ( t < 1.0 ) {
#else
    while ( t < LAT.BETA ) {
#endif
	// end edit Suzuki
	it = imin( NBODY , tau ); //nbody=2 it=0
	// edit Suzuki
	// 逆温度に対する規格化
#ifdef NORM
	EBSAMP -= (tau[it]-t) * IP.VertexDensity( x ) * LAT.BETA;
#else
	EBSAMP -= (tau[it]-t) * IP.VertexDensity( x );
#endif
	// end edit Suzuki

	if ( p[it]->top().isTerminal() ) break;
	t = tau[it];
	p[it]++;
	tau[it] = p[it]->topTime();
	x[it] = p[it]->X();
      }

      delete [] tau;
      delete [] x;
      delete [] p;
    
    }

    ACC[EB1].accumulate( EBSAMP );//kota n++; s1+=x; s2+=(x*x)
    ACC[EB2].accumulate( EBSAMP * EBSAMP );
    ACC[NH1].accumulate( MZUA * EBSAMP);

    //+++++ edit sakakura +++++//

    if (P.RUNTYPE == 1){
      ediag = MZUA;
    }
    if (P.RUNTYPE == 2){
      ediag = EBSAMP + double(NV);
      nkink   = NV; //number of kink
    }

    //+++++ edit sakakura +++++//

// WINDING NUMBER
//#ifdef WINDING
  int Dim = LAT.BD;
  double* wind = new double[Dim];


  for(int xi = 0; xi<Dim;xi++){wind[xi]=0;}


  for(int xi = 0; xi < LAT.NEDGE; xi++ ){


    int Asite=LAT.EDGE(xi).A;
    int Bsite=LAT.EDGE(xi).B;
    int dim=LAT.EDGE(xi).bd;

    Site& SITE = LAT.S(Asite);
    Segment& S0 = SITE.first();
    Site::iterator p(SITE);
    int xlast = S0.X() ;

    //count winding

    while ( ! (++p).atOrigin() ) {
      Segment& S = *p;
      if(xlast != S.X()){
	if((* S.bottom().S(0).getONSITE()).id()-1 == Bsite){
	  wind[dim] += (double) (S.bottom().S(1).X() - S.bottom().S(0).X());
	}else if((* S.bottom().S(2).getONSITE()).id()-1 == Bsite){
	  wind[dim] += (double) (S.bottom().S(3).X() - S.bottom().S(2).X());
	}
	xlast = S.X();
      }
    }

  }

  double wind2=0.0;
  for (int di=0; di<LAT.D; di++){
    double W=0;
    for (int db=0; db<Dim; db++){
      W+=wind[db]*LAT.vec[di][db]*LAT.L[di];
    }
    wind2 += W*W;
  }
   

  ACC[Wxy2].accumulate(wind2);
// WINDING NUMBER

  delete [] wind;
  //#endif

}

//======================================================================
// === edit rep sakakura ===
//void Measurement::setsummary() {
void Measurement::setsummary(int index) {
    // === edit rep sakakura ===

    using namespace Specific;

    double* X = new double[NACC]; //edit sakakura
    //  double X[NACC];

    for (int i=0; i<NACC; i++) ACC[i].average();
    for (int i=0; i<NACC; i++) X[i] = ACC[i].value;

    double B = LAT.BETA;
    double V = (double)LAT.NSITE;
    double D = (double)LAT.D;

    Q[ANV]  = X[NV1]/V;
    Q[ENE]  = (EBASE + X[EB1]/B)/V;

    Q[SPE]  = (X[EB2]-X[EB1]*X[EB1]-X[NV1])/V;

    // edit Suzuki
    // 逆温度に対する規格化
#ifdef NORM
    Q[LEN]  = X[LE1] * B ;
    Q[XMX]  = ALG.WDIAG * X[LE1] ;
#else
    Q[LEN]  = X[LE1];
    Q[XMX]  = ALG.WDIAG * X[LE1] / B;
#endif
    // end edit Suzuki
    //
    Q[AMZU] = X[MZU1]/V;

    // edit Suzuki
    // 逆温度に対する規格化
#ifdef NORM
    Q[BMZU] = X[MZUA1]/V;
    Q[SMZU] = X[MZU2]/V;
    Q[XMZU] = X[MZUA2]/V;

    Q[AMZS] = X[MZS1]/V;
    Q[BMZS] = X[MZSA1]/V;
    Q[SMZS] = X[MZS2]/V;
    Q[XMZS] = X[MZSA2]/V;
    Q[SMZQ] = X[MZQ2]/V;
#else
    Q[BMZU] = X[MZUA1]/V/B;
    Q[SMZU] = X[MZU2]/V;
    Q[XMZU] = X[MZUA2]/V/B/B;

    Q[AMZS] = X[MZS1]/V;
    Q[BMZS] = X[MZSA1]/V/B;
    Q[SMZS] = X[MZS2]/V;
    Q[XMZS] = X[MZSA2]/V/B/B;
    Q[SMZQ] = X[MZQ2]/V;
#endif


    Q[DS1]  = B*(X[NH1]-X[MZUA1]*X[EB1])/V;
    Q[W2]   = X[Wxy2] ;
    //    Q[RHOS] = ((double)X[Wxy2])*pow((double)LAT.L[0],(double)(2-LAT.D) ) / B;  // When Lx=Ly=Lz
    Q[RHOS] = X[Wxy2]*0.5/D/V/B;
    double XRHOXRHO=X[MZUA1]*X[MZUA1];
    Q[COMP] = B * V * ( X[MZUA2] - XRHOXRHO )/XRHOXRHO;


    // end edit Suzuki

    // === edit rep sakakura ===
    if (P.RUNTYPE == 0)for (int i=0; i<NPHY; i++) PHY[i].accumulate(Q[i]);
    if (P.RUNTYPE >> 0)for (int i=0; i<NPHY; i++) PHYX[i][index].accumulate(Q[i]);
    // === edit rep sakakura ===


#ifdef CF
  for(int icf=0; icf<NCF; icf++){
    for(int it=0; it<Ntau; it++){
      ACC4CF[icf][it].average();
      double factor = ALG.WDIAG * V /((double) DCF[icf]);
      PHY4CF[icf][it].accumulate(ACC4CF[icf][it].value * factor);
    }
  }
#endif
#ifdef SF
  for(int isf=0; isf<NSF; isf++){
    for(int it=0; it<Ntau; it++){
      ACC4SF[isf][it].average();
      PHY4SF[isf][it].accumulate(ACC4SF[isf][it].value/V);
    }
  }
#endif

  delete []X;
}
