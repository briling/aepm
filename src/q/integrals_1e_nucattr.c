#include "q.h"
#include "integrals.h"
#include "gc.h"

void rnlm_fill_nucattr_point(int lsum,
            rnlm_t rnlm, double ap, double P[3],
            mol * m, double * boys_array){

  for(int k=0; k<m->n; k++){
    preijkl_t prem;
    r3diff(prem.PQ, P, m->r+3*k);
    double r2 = r3dot(prem.PQ,prem.PQ);
    double T = ap * r2;
    boys_os(prem.F, lsum, T, boys_array);
    prem.alpha = ap;
    rnlm_fill_adds(lsum, m->q[k], &prem, rnlm);
  }
  return;
}

void rnlm_fill_nucattr_finite(int lsum,
    rnlm_t rnlm, double ap, double P[3],
    mol * m, double * boys_array){

  for(int k=0; k<m->n; k++){
    int qk = m->q[k];
    double aq, wq;
    finite_nuc_par(qk, &aq, &wq);
    double * Q  = m->r+3*k;
    preijkl_t prem;
    r3diff(prem.PQ, P, Q);
    double r2 = r3dot(prem.PQ,prem.PQ);
    double t1 = 1.0/(ap+aq);
    double t2 = ap*aq;
    prem.alpha = t1*t2;
    double T = prem.alpha * r2;
    double lambda = sqrt(t1) / t2;
    boys_os(prem.F, lsum, T, boys_array);
    rnlm_fill_adds(lsum, lambda * wq, &prem, rnlm);
  }
  return;
}

double add_prim_nucattr(
    int N, int L, int M,
    double d[L_MAX*2+1],
    double e[L_MAX*2+1],
    double f[L_MAX*2+1],
    rnlm_t rnlm){
  // finite: * 2.0*M_PI*M_PI*SQRTPI;
  // point:  * 2.0*M_PI / ap;
  double nucattr = 0.0;
  for(int N1=0; N1<=N; N1++){
    for(int L1=0; L1<=L; L1++){
      for(int M1=0; M1<=M; M1++){
        nucattr += d[N1]*e[L1]*f[M1] * rnlm[N1][L1][M1];
      }
    }
  }
  return nucattr;
}

/*-------------------------------------------------------*/

double nucattr_point(
    int nmax, int lmax, int mmax,
    double d[L_MAX*2+1],
    double e[L_MAX*2+1],
    double f[L_MAX*2+1],
    double P[3], double ap,
    mol * m, double * boys_array){

  int lsum = nmax + lmax + mmax;
  double nucattr = 0.0;
  for(int k=0; k<m->n; k++){
    double PC[3];
    r3diff(PC, P, m->r+k*3);
    double T = ap * r3dot(PC,PC);
    double F[BOS_NMAX+1];
    boys_os(F, lsum, T, boys_array);
    double t = 0.0;
    for(int N=0; N<=nmax; N++){
      for(int L=0; L<=lmax; L++){
        for(int M=0; M<=mmax; M++){
          t += d[N]*e[L]*f[M] * RNLMj(N, L, M, ap, PC, F);
        }
      }
    }
    nucattr += t * m->q[k];
  }
  return nucattr * 2.0*M_PI / ap;
}

double nucattr_finite(
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    double d[L_MAX*2+1],
    double e[L_MAX*2+1],
    double f[L_MAX*2+1],
    double P[3], double ap,
    mol * ml, double * boys_array){

  double nucattr = 0.0;
  int lsum   = n1+l1+m1 + n2+l2+m2;
  for(int k=0; k<ml->n; k++){
    int qk = ml->q[k];
    double aq, wq;
    finite_nuc_par(qk, &aq, &wq);
    double * Q  = ml->r+3*k;
    double PQ[3];
    r3diff(PQ, P, Q);
    double r2 = r3dot(PQ,PQ);
    double t1 = 1.0/(ap+aq);
    double t2 = ap*aq;
    double alpha = t1*t2;
    double T = alpha * r2;
    double F[BOS_NMAX+1];
    boys_os(F, lsum, T, boys_array);
    double lambda = 2.0*M_PI*M_PI*SQRTPI * sqrt(t1) / t2;
    nucattr += lambda * wq * add_ijs(
        n1,l1,m1, n2,l2,m2,
        d, e, f, PQ, alpha, F);
  }
  return nucattr;
}

