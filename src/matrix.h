
//============================================================================
//    Exact Calculation of Finite Size Spin Systems
//============================================================================

#include <cpplapack.h>
#include <blaswrap.h>
#include <cstdlib>
#include <complex>
//#include <ctime>
//#include <blaswrap.h>
//#include <cpplapack.h>
using namespace std;
using namespace CPPL;

//============================================================================
//
//============================================================================

complex<double> IUNIT(0.0,1.0);

//============================================================================
//    Display
//============================================================================

void dump(const vector<double>& V) {
  printf("\n");
  for (int i=0; i<V.size(); i++) {
    printf(" %8.3f", V[i]);
  }
  printf("\n");
}

//----------------------------------------------------------------------------

void dump(char* s, const vector<double>& V) {
  printf("\n");
  printf("%s\n",s);
  dump(V);
}

//----------------------------------------------------------------------------

void dump(const dgematrix& A, int Mmax = 10) {
  int M = A.m;
  if ( M > Mmax ) M = Mmax;
  int N = A.n;
  if ( N > Mmax ) N = Mmax;
  printf("\n");
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      printf(" %8.3f", A(i,j));
//    printf(" %2d", (int)(A(i,j)+0.1));
    }
    printf("\n");
  }
}

//----------------------------------------------------------------------------

void dump01(const dgematrix& A, int Mmax = 64) {
  int M = A.m;
  if ( M > Mmax ) M = Mmax;
  int N = A.n;
  if ( N > Mmax ) N = Mmax;
  printf("\n");
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      char x = '.';
      if ( abs(A(i,j)) > 1.0e-8 ) x = 'X';
      printf("%1c", x);
    }
    printf("\n");
  }
}

//----------------------------------------------------------------------------

void dump(char* s, const dgematrix& A) {
  printf("\n");
  printf("%s\n",s);
  dump(A);
}

//============================================================================
//    Diagonalization
//============================================================================

void diagonalize(dsymatrix& A, vector<double>& E, dgematrix& U) {
  dsymatrix A0 = A;
  vector<dcovector> V;
  A0.dsyev(E,V);
  for (int i=0; i<A.n; i++) {
    for (int j=0; j<A.n; j++) {
      U(i,j) = V[j](i);
    }
  }
  return;
}

//----------------------------------------------------------------------------

void diagonalize(dgematrix& A, vector<double>& E, dgematrix& U) {
  dsymatrix A0(A.n);
  for (int i=0; i<A.n; i++) {
    for (int j=0; j<=i; j++) {
      A0(i,j) = A(i,j);
    }
  }
  diagonalize( A0, E, U);
}

//============================================================================
//    Tensor Product
//============================================================================

dgematrix operator^(const dgematrix& A, const dgematrix& B) { 
  int m = A.m * B.m;
  int n = A.n * B.n;
  dgematrix C(m,n);
  for (int i0=0; i0<A.m; i0++) {
    for (int i1=0; i1<B.m; i1++) {
      int i = i0 + A.m * i1;
      for (int j0=0; j0<A.n; j0++) {
        for (int j1=0; j1<B.n; j1++) {
          int j = j0 + A.n * j1;
          C(i,j) = A(i0,j0) * B(i1,j1);
        }
      }
    }
  }
  return C;
}

//============================================================================
//    Complex Matrix
//============================================================================

class cmatrix {
public:
  long m;
  long n;
  CPPL::dgematrix re;
  CPPL::dgematrix im;

  void resize(const long m0, const long n0) {
    m = m0;
    n = n0;
    re.resize(m,n);
    im.resize(m,n);
  };

  void resize(const long n0) { resize(n0,n0); }

  cmatrix() {};
  cmatrix(const long n0) { resize(n0,n0); };
  cmatrix(const long m0, const long n0) { resize(m0,n0); };
  cmatrix(const cmatrix& X) { resize(X.m,X.n); };

  cmatrix& operator+=(const cmatrix& A) {
    if ( m != A.m ) { printf("cmatrix::operator+= >> ERROR;"); exit(0); }
    if ( n != A.n ) { printf("cmatrix::operator+= >> ERROR;"); exit(0); }
    re += A.re;
    im += A.im;
  };

  cmatrix& operator-=(const cmatrix& A) {
    if ( m != A.m ) { printf("cmatrix::operator+= >> ERROR;"); exit(0); }
    if ( n != A.n ) { printf("cmatrix::operator+= >> ERROR;"); exit(0); }
    re -= A.re;
    im -= A.im;
  };

  cmatrix& operator=(const cmatrix& A) {
    m = A.m;
    n = A.n;
    re = A.re;
    im = A.im;
  };

  void clear() { re.clear(); im.clear(); };

  void zero();

  void unity() { 
    resize(1,1);
    re(0,0) = 1.0;
    im(0,0) = 0.0;
  };

  void identity() {
    re.identity();
    im.identity();
  }

  const void dump(int Mmax);
  const void dump(char*);
  const void dump01(int Mmax);

};

//----------------------------------------------------------------------------

void cmatrix::zero() {
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      re(i,j) = 0.0;
      im(i,j) = 0.0;
    }
  }
}

//----------------------------------------------------------------------------

const void cmatrix::dump(int Mmax = 10) {
  int M = m;
  if ( M > Mmax ) M = Mmax;
  int N = n;
  if ( N > Mmax ) N = Mmax;
  printf("\n");
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
//    if ( j!= 0) printf(" |");
      printf(" %6.3f", re(i,j) );
//    printf(" %4.1f %4.1f", re(i,j), im(i,j));
    }
    printf("\n");
  }
}

//----------------------------------------------------------------------------

const void cmatrix::dump01(int Mmax = 64) {
  int M = m;
  if ( M > Mmax ) M = Mmax;
  int N = n;
  if ( N > Mmax ) N = Mmax;
  printf("\n");
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      char x = '.';
      if ( re(i,j)*re(i,j)+im(i,j)*im(i,j) > 1.0e-16 ) x = 'X';
      printf("%1c", x);
    }
    printf("\n");
  }
}

//----------------------------------------------------------------------------

const void cmatrix::dump(char* s) {
  printf("%s\n",s);
  dump();
}

//----------------------------------------------------------------------------

cmatrix operator+(const cmatrix& A, const cmatrix& B) {
  cmatrix C(A.m,B.n);
  C.re = A.re + B.re;
  C.im = A.im + B.im;
  return C;
}

//----------------------------------------------------------------------------

cmatrix operator-(const cmatrix& A, const cmatrix& B) {
  cmatrix C(A.m,B.n);
  C.re = A.re - B.re;
  C.im = A.im - B.im;
  return C;
}

//----------------------------------------------------------------------------

cmatrix operator*(const cmatrix& A, const cmatrix& B) {
  cmatrix C(A.m,B.n);
  C.re = A.re * B.re - A.im * B.im;
  C.im = A.re * B.im + A.im * B.re;
  return C;
}

//----------------------------------------------------------------------------

cmatrix operator*(const double a, const cmatrix& A) {
  cmatrix C(A.m,A.n);
  C.re = a * A.re;
  C.im = a * A.im;
  return C;
}

//----------------------------------------------------------------------------

cmatrix operator*(const complex<double> c, const cmatrix& A) {
  cmatrix C(A.m,A.n);
  C.re = c.real() * A.re - c.imag() * A.im;
  C.im = c.real() * A.im + c.imag() * A.re;
  return C;
}

//----------------------------------------------------------------------------

cmatrix t(const cmatrix& A) {
  cmatrix C(A.n,A.m);
  C.re = t(A.re);
  C.im = t(A.im);
  return C;
}

//----------------------------------------------------------------------------

cmatrix operator^(cmatrix& A, cmatrix& B) { 
  int m = A.m * B.m;
  int n = A.n * B.n;
  cmatrix C(m,n);
  C.re = ((A.re)^(B.re)) - ((A.im)^(B.im));
  C.im = ((A.re)^(B.im)) + ((A.im)^(B.re));
  return C;
}
