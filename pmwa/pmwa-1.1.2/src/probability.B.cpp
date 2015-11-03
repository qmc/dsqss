#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include <iostream>
using namespace std;


#include <Probability.h>

//#define Pdebug
//#define checkED
//#define checkR



int Pcomp(const void* p, const void* q){

  double x, y;
  x = ((Probability::Omega*)p) -> val;
  y = ((Probability::Omega*)q) -> val;

  if(x > y)
    return 1;
  else if(x < y)
    return -1;
  else
    return 0;
}

void Probability::look( Size *N, System *sp ){


  int oh, type; 
  double sql;
  int h;

  local_Et = 0.0; //local energy shift
  double rmin=0.0, Ri;

  // search the minimum value of the vertex density.
  for(int i=0;i<=nmax;i++){
   for(int j=0;j<=nmax;j++){
     for(int x=0;x<XMAX;x++){
       Ri = at_make(i, j, x);
       if(Ri<=rmin){ rmin=Ri; }
     }
   }
  }

  local_Et = tb  - rmin;
  sp->Et = z/2.0*PR->V*local_Et; 

  ////////

  for(int x=0; x<PR->V; x++) {
    double rumin=0.0;
    for(int i=0;i<=nmax;i++){
      Ri = au_make(i,x);
      if(Ri<=rumin) rumin=Ri;
    }

    cout << "rumin=" << rumin <<endl;
    local_Eu[x] = - rumin;
  }

  sp->Eu = 0.0;
  for(int x=0; x<PR->V; x++) sp->Eu += local_Eu[x];
  ////////

  rh_odd=sp->Htr;
  rh_even=sp->Htr;

  //************at,av2$BBeF~(B**************
  for(int i=0;i<=nmax;i++)for(int j=0;j<=nmax;j++)for(int x=0;x<XMAX;x++) at[i][j][x] = at_make(i, j, x);
  for(int i=0;i<=nmax;i++)for(int j=0;j<=nmax;j++)av2[i][j] = av2_make(i, j);
  for(int i=0; i<=nmax; i++ )for(int x=0; x<PR->V; x++ ) ru[i][x] = au_make( i, x );
  //*******************romax$B7hDj(B******************
  rtmax=rmin;
  for(int i=0; i<=nmax; i++ )for(int j=0; j<=nmax; j++ )for(int x=0;x<XMAX;x++)  if( rtmax < at[i][j][x] ) rtmax=at[i][j][x];

  rv2max=0.0;
  if(N->d==2) for(int i=0; i<=nmax; i++ )for(int j=0; j<=nmax; j++ ) if( rv2max < av2[i][j] ) rv2max=av2[i][j];


  for(int i=0; i<=nmax; i++ )for(int x=0; x<PR->V; x++ ) if( rumax[x] < ru[i][x] ) rumax[x]=ru[i][x];
  //******************* on-site scattering probability ******************
  for(int i=0; i<=nmax; i++ )for(int x=0; x<PR->V; x++ ) ru[i][x] /= rumax[x];

  for(int x=0; x<PR->V; x++ ) 
    for(int b=0;b<2;b++)//b=0 corresconds "oh=dir" : b =(oh==dir)? 0: 1;
      for(int i=0;i<=nmax;i++){//state before scattering
	int j=(b)? i-1: i+1; //state after scattering
	if(j<0||j>nmax) u[ b ][ i ][ x ] = 0.0;
	else if(ru[i][x]<=ru[j][x]) u[ b ][ i ][ x ] = 1.0; 
	else u[ b ][ i ][ x ] = Tuab( i, j, x );//when i is larger(L), if L->S then P = S/L, if L->L then 1.0 - S/L.

	if(PR->my_rank==0)cout <<"x="<<x<<",  b="<<b<<", ni="<<i<<", nj="<<j<<", u="<<u[b][i][x]<<endl;
      }


  //*****************************************************
  //*********************** two-site scattering probability *********************************

  int flaver=4;
  for(int x=0;x<XMAX;x++){
    for( h=0;h<2;h++){//$B%X%C%I$N1i;;;R(B


      oh = h*2 -1;
      
      for(int i=0;i<=nmax;i++){//$B:8$N%5%$%H$NN3;R?t(B
     
	sql = sqrt( i - h + 1.0 );
 
	for(int j=0;j<=nmax;j++){//$B1&$N%5%$%H$NN3;R?t(B
	  
	  if( h == i ) type = 5; //$B7h$^$j(B
	  else{
	    type=0;
	    Om[0].val = at[i][j][x];
	    Om[1].val = at[i+oh][j][x];
	    Om[2].val =( j != h )? tb : 0.0;
	    Om[3].val =( j == h )? tb : 0.0;
	  }


	for(int k=0; k<4; k++ )Om[k].num = k; //$B$[$s$H$NHV9f$r5-O?(B
	qsort( Om, 4, sizeof (Omega), Pcomp );
	for(int k=0; k<4; k++ )Tr[ Om[k].num ] = k; //$B$[$s$H$NHV9f$rF~$l$k$HBg$-$$=g$G$NHV9f$rJV$9(B


	if(type!=5) SolveWeightEquation(flaver); //
        //if(type!=5) Color(flaver);


//	for(int x=0;x<N->V;x++){//$B1&$N%5%$%HHV9f(B($B<!85$GDL$7(B)
	  for(int b=0;b<4;b++)//$B99?7A0$NG[CV(B
	    for(int a=0;a<4;a++){//$B99?78e$NG[CV(B 

    //******************* table ******************
	      if( type == 5 ) t[ h ][ a ][ b ][ i ][ j ][ x ] = 0.0;
	      else if( Om[Tr[b]].val != 0.0 ) t[ h ][ a ][ b ][ i ][ j ][ x ] = Wall[Tr[b]][Tr[a]]/Om[Tr[b]].val;
    //*****************************************************
#ifdef Pdebug
	      if(x==0)std::cout <<"h="<< h <<" a="<< a <<" b="<< b <<" i="<< i <<" j="<< j << "  type=" << type << "  Om[0]=" << Om[0].val<< "  Om[1]=" << Om[1].val<< "  Om[2]=" << Om[2].val<< "  Om[3]=" << Om[3].val << "  sql=" << sql << "  a_t=" << at[i][j][x] << "   t ="<< t[ h ][ a ][ b ][ i ][ j ][ x ]<< std::endl;
#endif
	    }
      }
    }
    }
  }

  



#ifdef Pdebug
double P;
for(int h=0;h<2;h++)
for(int i=0;i<=nmax;i++)
for(int j=0;j<=nmax;j++)
for( int b=0; b<4; b++){
 
  P=0.0;

  for( int a=0; a<4; a++){
    P += t[ h ][ a ][ b ][ i ][ j ][ x ];
  }

  cout<< "h="<< h << "  i="<<i << "  j="<< j << "  b="<< b <<"  P="<<P<<endl;
}

exit(1);
#endif


}

void Probability::Color( int cmax ){


  for( int i=0; i<cmax; i++){ ex_Wall[i] = ex_Penki[i] = Om[i].val; } //$BJI$NEI$j;D$7!"%Z%s%-$N;D$j!"=E$_(B
  for( int i=0; i<cmax; i++) for( int p=0; p<cmax; p++) Wall[i][p] = 0.0;

  double total_Penki, paint;
  
  for( int penki=0; penki<cmax-1; penki++){

    if(ex_Penki[penki]==0.0)continue;
    total_Penki=0.0;
    for( int kabe=penki+1; kabe<cmax; kabe++) total_Penki += Om[kabe].val;
    for( int kabe=penki+1; kabe<cmax; kabe++){

        paint = ex_Wall[penki]*(Om[kabe].val/total_Penki);
        Wall[kabe][penki] = paint;//penki$BHVL\$N%Z%s%-$G(Bkabe$BHVL\$NJI$rEI$k!#(B//$B$3$l$G(Bpenki$BHVL\$N%Z%s%-$r;H$$2L$?$7$?!#(B
        Wall[penki][kabe] = paint;//kabe$BHV$a$N%Z%s%-$G(Bpenki$BHVL\$NJI$r<+J,$,EI$i$l$?J,NL$HF1$8$@$1EI$jBX$($9!#(B//$B$3$l$G(Bpenki$BHVL\$NJI$rEI$j=*$($?!#(B
        ex_Wall[kabe]-=paint; //kabe$BHVL\$NJI$NEI$j;D$7(B
        
    }

  }

  Wall[cmax-1][cmax-1]=ex_Wall[cmax-1];


#ifdef Pdebug
    for( int i=0; i<4; i++) for( int p=0; p<4; p++) cout<<"wall=" << i << "  penki=" << p <<  "  wij="<< Wall[i][p] <<endl;
 #endif


}

void  Probability::SolveWeightEquation(int cmax){

  for( int i=0; i<cmax; i++){ ex_Wall[cmax-1-i] = ex_Penki[cmax-1-i] = Om[i].val; } //$BJI$NEI$j;D$7!"%Z%s%-$N;D$j!"=E$_(B
  for( int i=0; i<cmax; i++) for( int p=0; p<cmax; p++) Wall[i][p] = 0.0;
    
  int p,q;
  double V_first;
  double V_second;
  double V_third;
  double V_first_new;
  double V_second_new;
  int N_first;   // the number of the largest elements
  int N_second;  // the number of the second largest

// $B=EJ#$r$O$V$$$F:G=i$+$i#3$D$N=E$_$H$=$NHV9f$r<hF@!%(B
  while ( ex_Wall[1] > EPS ) {
    V_first = ex_Wall[0];
    for (p=0; p<cmax; p++) if ( ex_Wall[p] != V_first ) break;
    N_first = p;
    if ( p == cmax ) {//all are same
      V_second = 0.0;
      V_third = 0.0;
      N_second = 0;
    } else {
      V_second = ex_Wall[p];
      for (q=p; q<cmax; q++) if ( ex_Wall[q] != V_second ) break;
      if ( q == cmax ) {
        V_third = 0.0;
      } else {
        V_third = ex_Wall[q];
      }
      N_second = q-p;
    }

//$B0J2<$G$O!$(Bex_Wall[i] <= w_i $B$+$i=PH/$7$F!$(BWall(i,j) = w_{ij} $B$r=g<!A}2C$5$;$J$,$i!$(B
//$B$=$NJ,!$(Bex_Wall[i] $B$r2<$2$F$$$/!%$9$Y$F$N(B ex_Wall[i] $B$,#0$K$J$C$?$i=*N;!%(B

    double dw1; // decrement of the first largest element
    double dw2; // decrement of the second largest element
    if ( N_first == 1 ) {  
//$B:GBg%&%'%$%H$,C1FH$N>uBV$N$H$-!$:GBg%&%'%$%H$HBh#2%&%'%$%H$N4V$NA+0\$r(B
//$BF3F~$7$F$=$l$i$r2<$2$k!%(B
      double x = V_first - V_second;
      double y = (double)(N_second - 1) * (V_second - V_third);
      if ( x < y ) { 
//$B:GBg%&%'%$%H$,Bg$-$/$J$$$H$-!$:GBg%&%'%$%H$,Bh#2%&%'%$%H$KEy$7$/$J$k$^$G2<$2$k!%(B
        dw1 = (V_first-V_second) / ( 1.0 - 1.0 / (double)(N_second));
        dw2 = dw1 / (double)N_second;
        V_second_new = V_second - dw2;
        V_first_new  = V_second_new;
      } else { 
// $B:GBg%&%'%$%H$,Bg$-$$$H$-!$Bh#2%&%'%$%H$,Bh#3%&%'%$%H$HEy$7$/$J$k$^$G2<$2$k!%(B
        dw2 = V_second - V_third;
        dw1 = dw2 * (double)N_second;
        V_second_new = V_third;
        V_first_new  = V_first - dw1;
      }
      ex_Wall[0] = V_first_new;
      for (int i=1; i<1+N_second; i++) {
        Wall[cmax-1][cmax-1-i] += dw2;
        Wall[cmax-1-i][cmax-1] += dw2;
        ex_Wall[i] = V_second_new;
      }
    } else { 
//$BJ#?t$N>uBV$,:GBg%&%'%$%H$r$H$k$H$-!$Aj8_$K0\$jJQ$o$kA+0\3NN($rF3F~$7$F!$(B
//$B$=$l$i$rBh#2%&%'%$%H$HF1$8$K$J$k$^$G2<$2$k!%(B
      dw1 = (V_first - V_second) / (double)( N_first-1 );
      for (int i=0; i<N_first; i++) {
        ex_Wall[i] = V_second;
        for (int j=0; j<N_first; j++) {
          if ( i==j ) continue;
          Wall[cmax-1-i][cmax-1-j] += dw1;
        }
      }
    }
  }
  if ( ex_Wall[0] > 0.0 ) {
    Wall[cmax-1][cmax-1] += ex_Wall[0];
    ex_Wall[0] = 0.0;
  }

}

double Probability::Tuab(int p, int q, int x){ // p(L)->q(S)

  return au_make( q, x )/au_make( p, x );

}


double Probability::Tv2ab(int p, int q){ //v2$B%P!<%F%C%/%9$N3NN(9TNs@.J,!!(B($BA}$($k$H$-$NF)2a3NN((B)


  return (av2_make( 1, q )/av2_make( 0, q ));

}

double Probability::at_make(int p, int q, int  x){ //t$B%P!<%F%C%/%9L)EY$rJV$9(B

  //double Ht = dim*(double)( p + q ) - V1*(double)(p*q); 
  //  double Ht = (double)( mu1[x]*p + mu2[x]*q )/z - V1*(double)(p*q); 
  double Ht = - V1*(double)(p*q); 
  // if(Ht + local_Et < 0.0){cout<<"negative sign! (Hv)"<<endl; exit(1);}

  return (Ht + local_Et); 

}

double Probability::av2_make(int p, int q){ //t$B%P!<%F%C%/%9L)EY$rJV$9(B

  double Ht = V2*(double)(2.0 - p*q); 

  return Ht; 

}


double Probability::au_make(int p, int x){ //u$A%P$B!<(B$B$A%F%C%/%9C\6H$r75$9(B

  double Hu = - Ubb/2.0*p*( p - 1 ) -  ep[x]*p;

  //  if(Hu + local_Eu[x] < 0.0){cout<<"negative sign! (Hu)"<<endl; exit(1);}

  return Hu + local_Eu[x];

}


Probability::Probability( Size *N, System *sp, Parallel *m_PR ){
  
  PR=m_PR;
  Ubb = sp->Ubb;
  V1 = sp->Vb1;
  V2 =(N->d==2)? sp->Vb2/2.0: 0.0;//2D$B$N$H$-$@$1$K$7$F$*$/(B
  nmax = sp->nmax;

  tb = sp->tb;
  z = 2.0*N->d;
  XMAX = 2; 

  newcall( ep, PR->V ); 
  newcall( local_Eu, PR->V );
  newcall(rumax,PR->V);
  newcall( ru, nmax+1, PR->V );
  

  newcall( t, 2, 4, 4, nmax+1, nmax+1, XMAX );
  newcall( u, 2, nmax+1, PR->V );
  newcall( v2, nmax+1, nmax+1 );
  newcall( w, 5, nmax+1 );
  newcall( at, nmax+1, nmax+1, XMAX );
  newcall( av2, nmax+1, nmax+1 );
 

  newcalls( Om, 4 );
  newcall( Tr, 4 );
  newcall(ex_Wall,4);
  newcall(ex_Penki,4);
  
  newcall(Wall,4,4);

  newcall(FracW,6,6,XMAX );
  newcall(weight,6);
  
 
}


void Probability::LocalPotential( double _Vc, double mu, long seed ){

  mu1[0]= - mu + _Vc;
  mu1[1]= - mu - _Vc;
  
  seed = PR->ns + seed*PR->Nsdiv;
  My_rdm LP(seed);
  double A0=1.0;//staggerd lattice 

  for(int z=0; z<PR->z; z++)
    for(int y=0; y<PR->y; y++)
      for(int x=0; x<PR->x; x++){
	int i = x + y*PR->x + z*PR->x*PR->y;
#ifdef checkED
	ep[i] = - mu + _Vc*(double)( x - Nx/2 )*( x - Nx/2 );
#else
	//ep[i]=mu1[(x+y+z)%2]; //sp->Vc*cos(A0*x*M_PI);//super lattice
	ep[i] = - mu + _Vc*(1.0 - 2.0*LP.rdm());
#endif
      }
#ifdef checkR
  char fname[256];
  sprintf(fname,"potential_S%02dT%02dR%03d",PR->ns,PR->nt,PR->nr);
  std::ofstream fout(fname);
  for(int i=0;i<PR->V;i++) fout<<i<<" "<<ep[i]<<endl;
  fout.close();
#endif

}


void Probability::LookUpTable( Size *N, System *sp ){

  dim = sp->mu/z;
  sp->Ev2 = 4.0*N->d*PR->V*V2;
  Nx=N->x;

  LocalPotential(sp->Vc, sp->mu, sp->seed);

  look( N, sp );

  for(int x=0;x<XMAX;x++){

    weight[0]=at[0][0][x];
    weight[1]=at[0][1][x];
    weight[2]=at[1][0][x];
    weight[3]=at[1][1][x];
    weight[4]=weight[5]=tb;

    for(int i=0;i<6;i++)for(int j=0;j<6;j++)FracW[j][i][x]=weight[i]/weight[j];
  }

}

Probability::~Probability(){

 delcall(ep);
  
  delcall( t, 2, 4, 4, nmax+1, nmax+1 );
  delcall( u, 2, nmax+1 );
  delcall( v2, nmax+1 );

  delcall( w, 5 );

  delcall( at, nmax+1, nmax+1 );
  delcall( av2, nmax+1 );

  delcall( ru, nmax+1 );
  delcall( rumax );
  delcall( local_Eu );
  

  delcall( Om );
  delcall( Tr );
  delcall(ex_Wall);
  delcall(ex_Penki);
  
  delcall(Wall,4);

  delcall(FracW,6,6);
  delcall(weight);
  

}



