#include "matrix.h"

#define MAX(x,y) ( ((x) > (y)) ? (x) : (y) )
#define MIN(x,y) ( ((x) < (y)) ? (x) : (y) )

void canorth(double * s, double * b, double * d, unsigned int n, double jeps, int jrot, FILE * f){
  /* S*C = D*C
   * D := diag(S)
   * B := (C+)/sqrt(s)
   * S := waste
   */
  mx_id(n, b);
  jacobi(s, b, d, n, jeps, jrot, f);
  for(int i=0; i<n; i++){
    vecscal(n, b+i*n, 1.0/sqrt(d[i]));
  }
  return;
}

double s_invsqrt_canorth(size_t n, double * S, double * S1, double * X, double jeps, int jrot, FILE * f){
  /* S1 := 1/sqrt(S)
   * X  := (C+)/sqrt(S)
   */
  double * d = malloc(sizeof(double)*n);
  veccp(symsize(n), S1, S);
  mx_id(n, X);
  jacobi(S1, X, d, n, jeps, jrot, f);
  double dmax = d[0];
  double dmin = d[0];
  for(int i=1; i<n; i++){
    dmax = MAX(dmax, d[i]);
    dmin = MIN(dmin, d[i]);
  }
  for(int i=0; i<n; i++){
    d[i] = 1.0/sqrt(d[i]);
  }
  mx_BDiagB(S1, X, d, n);
  for(int i=0; i<n; i++){
    vecscal(n, X+i*n, d[i]);
  }
  free(d);
  return dmin/dmax;
}

void mx_BDiagB(double * s, double * b, double * d, unsigned int n){
  /* S = (B+)DB */
  for(int i=0; i<n; i++){
    for(int j=i; j<n; j++){
      double t = 0.0;
      for(int k=0; k<n; k++){
        t += b[k*n+i] * b[k*n+j] * d[k];
      }
      s[mpos(i,j)] = t;
    }
  }
  return;
}

static void mx_BDiagB_odd(double * s, double * b, double * d, unsigned int n, unsigned int odd){
  for(int i=0; i<n; i++){
    for(int j=i; j<n; j++){
      double t = 0.0;
      for(int k=odd; k<n; k++){
        t += b[k*n+i] * d[k] * b[k*n+j];
      }
      s[mpos(i,j)] = t;
    }
  }
  return;
}

int mx_invsqrt1(double * s, double * b, double * d,
                 unsigned int n, double oddeps,
                 double jeps, int jrot, FILE * f){
  /* finds first 'odd' eigenvalues < 'oddeps'
   * and returns their number
   */
  mx_id(n, b);
  jacobi(s, b, d, n, jeps, jrot, f);
  eigensort(n, d, b);
  int odd = 0;
  for(unsigned int i=0; i<n; i++){
    if(d[i]<oddeps){
      odd++;
    }
    else{
      d[i] = 1.0/sqrt(d[i]);
    }
  }
  mx_BDiagB_odd(s, b, d, n, odd);
  return odd;
}

static int mx_inv_eigen_inner(
    double * s, double * b, double * d,
    unsigned int n, double oddeps,
    double jeps, int jrot, FILE * f){
  /* finds first 'odd' eigenvalues < 'oddeps'
   * and returns their number
   */
  mx_id(n, b);
  jacobi(s, b, d, n, jeps, jrot, f);
  eigensort(n, d, b);
  int odd = 0;
  for(size_t i=0; i<n; i++){
    if(d[i]<oddeps){
      odd++;
    }
    else{
      d[i] = 1.0/d[i];
    }
  }
  mx_BDiagB_odd(s, b, d, n, odd);
  return odd;
}

int mx_inv_eigen(int n, double * Ai, double * A, double oddeps){
  double * d = malloc(sizeof(double)*(n*n+n));
  double * b = d + n;
  veccp(symsize(n), Ai, A);
  int odd = mx_inv_eigen_inner(Ai, b, d, n, oddeps, 1e-15, 20, NULL);
  free(d);
  return odd;
}

