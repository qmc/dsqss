//======================================================================

void Measurement::measure() {

  using namespace Specific;
  int NV = LAT.countVertices();

  ACC[NV1].accumulate( ((double)NV) );

  double* val = new double[ALG.NXMAX]; //edit sakakura
  //  double val[ALG.NXMAX];

  for (int i=0; i<ALG.NXMAX; i++) {
    val[i] = 0.5 * (double)(2 * i - ALG.NXMAX + 1);
  }
  double phase[2] = { 1.0 , -1.0 };

  double MZU  = 0.0;
  double MZUA = 0.0;
  double MZS  = 0.0;
  double MZSA = 0.0;

  for (int s=0; s<LAT.NSITE; s++) {
    Site& SITE = LAT.S(s);
    int mt = SITE.getMTYPE();
    double ph = phase[mt];
    Segment& S0 = SITE.first();
    double mz0 = val[S0.X()];
    Site::iterator p(SITE);
    double mza0 = 0.0;
    while ( ! (++p).atOrigin() ) {
      Segment& S = *p;
      double mz = val[S.X()];
      mza0 += mz * S.length();
    }

    MZU  +=      mz0;
    MZUA +=      mza0;
    MZS  += ph * mz0;
    MZSA += ph * mza0;
  }

  ACC[MZU1].accumulate(MZU);
  ACC[MZU2].accumulate(MZU * MZU);
  ACC[MZUA1].accumulate(MZUA);
  ACC[MZUA2].accumulate(MZUA * MZUA);
  ACC[MZS1].accumulate(MZS);
  ACC[MZS2].accumulate(MZS * MZS);
  ACC[MZSA1].accumulate(MZSA);
  ACC[MZSA2].accumulate(MZSA * MZSA);

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

    //+++++ edit sakakura +++++//

    if (P.RUNTYPE == 1){
      ediag = MZUA;
    }
    if (P.RUNTYPE == 2){
      ediag = EBSAMP + double(NV);
      nkink   = NV; //number of kink
    }

    //+++++ edit sakakura +++++//
    delete []val;

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
#else
    Q[BMZU] = X[MZUA1]/V/B;
    Q[SMZU] = X[MZU2]/V;
    Q[XMZU] = X[MZUA2]/V/B/B;

    Q[AMZS] = X[MZS1]/V;
    Q[BMZS] = X[MZSA1]/V/B;
    Q[SMZS] = X[MZS2]/V;
    Q[XMZS] = X[MZSA2]/V/B/B;
#endif
    // end edit Suzuki

    // === edit rep sakakura ===
    if (P.RUNTYPE == 0)for (int i=0; i<NPHY; i++) PHY[i].accumulate(Q[i]);
    if (P.RUNTYPE >> 0)for (int i=0; i<NPHY; i++) PHYX[i][index].accumulate(Q[i]);
    // === edit rep sakakura ===

    delete []X;
  }
