#ifndef INTEGRALS_H
#define INTEGRALS_H

#include "boys.h"

double def0(int n1, int n2, double PAx, double PBx, double p21);

void filldef1(int n1, int n2, double d[L_MAX*2+1],
              double PAx, double PBx, double p21);
void filldef1_self(int n1, int n2, double d[L_MAX*2+1], double p21);
void filldef_all(int l1, int l2,
    double D1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    double PAx, double PBx, double p21);

double RNLMj(int N, int L, int M, double alpha, double r[3], double * F);

/*-------------------------------------------------------------------------------------------------*/

static inline void filldef(int n1, int l1, int m1,
                           int n2, int l2, int m2,
                           double d[L_MAX*2+1],
                           double e[L_MAX*2+1],
                           double f[L_MAX*2+1],
                           double PA[3], double PB[3], double p21){
  filldef1(n1, n2, d, PA[0], PB[0], p21);
  filldef1(l1, l2, e, PA[1], PB[1], p21);
  filldef1(m1, m2, f, PA[2], PB[2], p21);
  return;
}

static inline void filldef_self(int n1, int l1, int m1,
                                int n2, int l2, int m2,
                                double d[L_MAX*2+1],
                                double e[L_MAX*2+1],
                                double f[L_MAX*2+1],
                                double p21){
  filldef1_self(n1, n2, d, p21);
  filldef1_self(l1, l2, e, p21);
  filldef1_self(m1, m2, f, p21);
  return;
}

static inline double txyz(int n1, int n2, double PAx, double PBx, double a, double b, double p2_1){
  double t = + 4.0*a*b  * def0(n1+1, n2+1, PAx, PBx, p2_1);
  if(n1){
    t +=     - 2.0*n1*b * def0(n1-1, n2+1, PAx, PBx, p2_1);
      if(n2){
        t += + n1*n2    * def0(n1-1, n2-1, PAx, PBx, p2_1);
      }
  }
  if(n2){
    t +=     - 2.0*n2*a * def0(n1+1, n2-1, PAx, PBx, p2_1);
  }
  return t;
}

static inline double txyz1(
    int n1, int n2, double a, double b,
    double d0[L_MAX+1+1][L_MAX+1+1]){

  double t = + 4.0*a*b  * d0[n1+1][n2+1];
  if(n1){
    t +=     - 2.0*n1*b * d0[n1-1][n2+1];
      if(n2){
        t += + n1*n2    * d0[n1-1][n2-1];
      }
  }
  if(n2){
    t +=     - 2.0*n2*a * d0[n1+1][n2-1];
  }
  return t;
}

double add_s00(
    int n, int l, int m,
    int n1, int l1, int m1,
    double a, double w);
double add_ij_overlap_self(
    int n, int l, int m, int n1, int l1, int m1,
    double p21);
void add_ij(
    int n, int l, int m, int n1, int l1, int m1,
    double a, double b, double Rij2, double w,
    double ri[3], double rj[3], mol * ml,
    int finite, atomo * u_rel,
    double ans[6], double * boys_array);

double add_ijs(
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    double d1[L_MAX*2+1],
    double e1[L_MAX*2+1],
    double f1[L_MAX*2+1],
    double PQ[3], double alpha, double * F);
double add_ijs_grad(
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    double d1[L_MAX*2+1],
    double e1[L_MAX*2+1],
    double f1[L_MAX*2+1],
    double PQ[3], double alpha, double * F, double a);
double add_ijs1(
    double rnlm[BOS_NMAX+1][BOS_NMAX+1][BOS_NMAX+1],
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    double d1[L_MAX*2+1],
    double e1[L_MAX*2+1],
    double f1[L_MAX*2+1]);

static inline void gprod(double ai, double aj, double p1, double ri[3], double rj[3], double P[3], double PA[3], double PB[3]){
  r3sums(P, ri, ai, rj, aj);
  r3scal(P, p1);
  r3diff(PA, P, ri);
  r3diff(PB, P, rj);
  return;
}

double nucattr_point(
    int nmax, int lmax, int mmax,
    double d[L_MAX*2+1],
    double e[L_MAX*2+1],
    double f[L_MAX*2+1],
    double P[3], double ap,
    mol * m, double * boys_array);
double nucattr_finite(
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    double d[L_MAX*2+1],
    double e[L_MAX*2+1],
    double f[L_MAX*2+1],
    double P[3], double ap,
    mol * ml, double * boys_array);
double kinetic_scale(
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    double a, double b, double u,
    double ra[3], double rb[3], double ru[3],
    double rab[3]);
double esp_at_point(
    int nmax, int lmax, int mmax,
    double d[L_MAX*2+1],
    double e[L_MAX*2+1],
    double f[L_MAX*2+1],
    double P[3], double ap,
    double R[3], double * boys_array);

#endif
