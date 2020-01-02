#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vecn.h"

#define symsize(M) (((M)*(M)+(M))/2)

#ifndef MPOS_IS
#define MPOS_IS
static inline unsigned int mpos(unsigned int i, unsigned int j){
/* A[i+j*(j+1)/2], i <= j, 0 <= j < N */
  return (i)+(((j)*((j)+1))>>1);
}
#define MPOSIF(i,j)   ((i)<=(j)? mpos((i),(j)):mpos((j),(i)))
#endif

void     mx_id            (unsigned int n, double * a);
void     mx_print         (unsigned int n, double * a, FILE * f);
void     mx_nosym_print   (unsigned int n, double * a, FILE * f);
void     mx_transp        (unsigned int n, double * a);
void     mx_symmultrectmx (unsigned int n, unsigned int m, double * p, double * a, double * b);
void     mx_multsymmx     (unsigned int n, double * p, double * a, double * b);
void     mx_BHBt_sym      (unsigned int n, double * h, double * b);
void     mx_BHBt_sym2     (unsigned int n, unsigned int m, double * r, double * h, double * b);
void     mx_BdiagBt       (unsigned int n, double * h, double * v, double * b);
void     mx_sqr2          (unsigned int m, unsigned int n, double * p,  double * a);
void     mx_multvec       (unsigned int n, unsigned int m, double * r,  double * a, double * v);
void     mx_symmultvec    (unsigned int n, double * r, double * a, double * v);
double   mx_vecdot        (unsigned int n, double * u, double * a, double * v);
double   mx_vecdot_nosym  (unsigned int n, double * u, double * a, double * v);

void     jacobi    (double * a, double * b, double * d, unsigned int n, double eps, unsigned int rot, FILE * f);
void     eigensort (int n, double * val, double * vec);
double * smalldiag (int n, double * A);

void   canorth           (double * s, double * b, double * d, unsigned int n, double jeps, int jrot, FILE * f);
double s_invsqrt_canorth (size_t n, double * S, double * S1, double * X, double jeps, int jrot, FILE * f);
int    mx_invsqrt1       (double * s, double * b, double * d, unsigned int n, double oddeps, double jeps, int jrot, FILE * f);
void   mx_BDiagB         (double * s, double * b, double * d, unsigned int n);
int    mx_inv_eigen      (int n, double * Ai, double * A, double oddeps);

