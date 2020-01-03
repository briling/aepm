#include "q.h"
#include "integrals.h"
#include "gc.h"

double add_ijs(
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    double d1[L_MAX*2+1],
    double e1[L_MAX*2+1],
    double f1[L_MAX*2+1],
    double PQ[3], double alpha, double * F){
  double d = 0.0;
  for(int N1=0; N1<=n1+n2; N1++){
    for(int L1=0; L1<=l1+l2; L1++){
      for(int M1=0; M1<=m1+m2; M1++){
        d += d1[N1]*e1[L1]*f1[M1] * RNLMj(N1, L1, M1, alpha, PQ, F);
      }
    }
  }
  return d;
}

double add_ijs_grad(
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    double d1[L_MAX*2+1],
    double e1[L_MAX*2+1],
    double f1[L_MAX*2+1],
    double PQ[3], double alpha, double * F, double a){
  double d  = 0.0;
  for(int N1=0; N1<=n1+n2; N1++){
    for(int L1=0; L1<=l1+l2; L1++){
      for(int M1=0; M1<=m1+m2; M1++){
        d += d1[N1]*e1[L1]*f1[M1] * (
          RNLMj(N1+2, L1, M1, alpha, PQ, F)+
          RNLMj(N1, L1+2, M1, alpha, PQ, F)+
          RNLMj(N1, L1, M1+2, alpha, PQ, F) );
      }
    }
  }
  return -0.25/(a*a) * d;
}

double add_ijs1(
    rnlm_t rnlm,
    int n1, int l1, int m1, int n2, int l2, int m2,
    double d1[L_MAX*2+1], double e1[L_MAX*2+1], double f1[L_MAX*2+1]){
  double d = 0.0;
  for(int N1=0; N1<=n1+n2; N1++){
    for(int L1=0; L1<=l1+l2; L1++){
      for(int M1=0; M1<=m1+m2; M1++){
        d += d1[N1]*e1[L1]*f1[M1] * rnlm[N1][L1][M1];
      }
    }
  }
  return d;
}

