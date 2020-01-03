#include "q.h"
#include "qinit.h"

int measure_ortho(int N, int M, int core,
    double * S, double * C0, double * C){
  double * Asym = malloc(symsize(N/2)*sizeof(double));
  int core2 = C_maxoverlap(N/2, core, M, Asym, C, C0, S);
  free(Asym);
  return core2-core;
}

double measure_minoverlap(int N, int M, int core,
    double * S, double * C0, double * C){
  double * Asym = malloc(symsize(N/2)*sizeof(double));
#if 1
  C_maxoverlap(N/2, core, M, Asym, C, C0, S);
#else
  mx_BHBt_sym2(N/2, M, Asym, P0, C);
#endif
  double * v4 = smalldiag(N/2, Asym);
  double d = 1.0 - v4[core];
  free(v4);
  free(Asym);
  return d;
}

double measure_E(int N, int M, int core,
    double * F0, double * V0, double * C){

  double * CFC = malloc(symsize(N/2)*sizeof(double));
  mx_BHBt_sym2(N/2, M, CFC, F0, C);
  double * v1 = smalldiag(N/2, CFC);
  double d1 = 0.0;
  for(int i=core; i<N/2; i++){
    d1 += v1[i] - V0[i];
  }
  free(v1);
  free(CFC);
  return d1;
}

double measure_E0(int N, int M, int core,
    double * F0, double * V0, double * C){

  double * CFC = malloc(symsize(N/2)*sizeof(double));
  mx_BHBt_sym2(N/2, M, CFC, F0, C);
  double d1 = 0.0;
  for(int i=core; i<N/2; i++){
    d1 += CFC[mpos(i,i)] - V0[i];
  }
  free(CFC);
  return d1;
}

double measure_S(int N, int M, int core, double * P0, double * C){
  double * CPC = malloc(symsize(N/2)*sizeof(double));
  mx_BHBt_sym2(N/2, M, CPC, P0, C);
  double * v1 = smalldiag(N/2, CPC);
  double d6 = 0.0;
  for(int i=core; i<N/2; i++){
    d6 += 1.0 - v1[i];
  }
  free(v1);
  free(CPC);
  return d6;
}

double measure_S0(int N, int M, int core, double * P0, double * C){

  double * CPC = malloc(symsize(N/2)*sizeof(double));
  mx_BHBt_sym2(N/2, M, CPC, P0, C);
  double d1 = 0.0;
  for(int i=core; i<N/2; i++){
    d1 += 1.0 - CPC[mpos(i,i)];
  }
  free(CPC);
  return d1;
}

double measure_E0_grad(
    int npar, int N, int M, int core, int * mpl,
    double * F0, double * V0,
    double * C,  double * V,
    double * Gs, double * g1){

  double * C_F0_C = malloc(2*symsize(M)*sizeof(double));
  double * C_dHdA_C = C_F0_C + symsize(M);
  mx_BHBt_sym2(M, M, C_F0_C, F0, C);

  for(int p=0; p<npar; p++){
    if(!mpl[p]) continue;
    mx_BHBt_sym2(M, M, C_dHdA_C, Gs+symsize(M)*p, C);
    for(int i=core; i<N/2; i++){
      for(int j=0; j<M; j++){
        if(j>=core && j<N/2) continue;
        g1[p] += 2.0 * C_F0_C[MPOSIF(i,j)] * C_dHdA_C[MPOSIF(i,j)] / (V[i]-V[j]);
      }
    }
  }

  double d1 = 0.0;
  for(int i=core; i<N/2; i++){
    d1 += C_F0_C[mpos(i,i)] - V0[i];
  }
  free(C_F0_C);
  return d1;
}

double measure_S0_grad(
    int npar, int N, int M, int core, int * mpl,
    double * P0,
    double * C,  double * V,
    double * Gs, double * g1){

  double * C_P0_C = malloc(2*symsize(M)*sizeof(double));
  double * C_dHdA_C = C_P0_C + symsize(M);
  mx_BHBt_sym2(M, M, C_P0_C, P0, C);

  for(int p=0; p<npar; p++){
    if(!mpl[p]) continue;
    mx_BHBt_sym2(M, M, C_dHdA_C, Gs+symsize(M)*p, C);
    for(int i=core; i<N/2; i++){
      for(int j=0; j<M; j++){
        if(j>=core && j<N/2) continue;
        g1[p] -= 2.0 * C_P0_C[MPOSIF(i,j)] * C_dHdA_C[MPOSIF(i,j)] / (V[i]-V[j]);
      }
    }
  }

  double d1 = 0.0;
  for(int i=core; i<N/2; i++){
    d1 += 1.0 - C_P0_C[mpos(i,i)];
  }
  free(C_P0_C);
  return d1;
}

double measure_E_grad(
    int npar, int N, int M, int core, int * mpl,
    double * F0, double * V0,
    double * C,  double * V,
    double * Gs, double * g1){

  size_t size = 2*symsize(M) + (N/2)*(N/2);
  double * Q        = malloc(size*sizeof(double));
  double * C_dHdA_C = Q        + symsize(M);
  double * dQdA     = C_dHdA_C + symsize(M);

  // Q = Ct*F0*C
  mx_BHBt_sym2(M, M, Q, F0, C);
  // Q*B = B*D
  double * D = smalldiag(N/2, Q);
  double * B = D + N/2;

  for(int p=0; p<npar; p++){
    if(!mpl[p]) continue;
    mx_BHBt_sym2(M, M, C_dHdA_C, Gs+symsize(M)*p, C);
    for(int i=0; i<N/2; i++){
      for(int j=0; j<N/2; j++){
        double s = 0.0;
        for(int k=N/2; k<M; k++){
          s += Q[MPOSIF(i,k)] * C_dHdA_C[MPOSIF(k,j)] / (V[j]-V[k]);
        }
        dQdA[i*(N/2)+j] = 2.0*s;
      }
    }
    double g = 0.0;
    for(int i=core; i<N/2; i++){
      double * Bi = B+(N/2)*i;
      g += mx_vecdot_nosym(N/2, Bi, dQdA, Bi);
    }
    g1[p] += g;
  }
  free(Q);

  double d1 = 0.0;
  for(int i=core; i<N/2; i++){
    d1 += D[i] - V0[i];
  }
  free(D);
  return d1;
}

double measure_S_grad(
    int npar, int N, int M, int core, int * mpl,
    double * P0,
    double * C,  double * V,
    double * Gs, double * g1){

  size_t size = 2*symsize(M) + (N/2)*(N/2);
  double * Q        = malloc(size*sizeof(double));
  double * C_dHdA_C = Q        + symsize(M);
  double * dQdA     = C_dHdA_C + symsize(M);

  // Q = Ct*P0*C
  mx_BHBt_sym2(M, M, Q, P0, C);
  // Q*B = B*D
  double * D = smalldiag(N/2, Q);
  double * B = D + N/2;

  for(int p=0; p<npar; p++){
    if(!mpl[p]) continue;
    mx_BHBt_sym2(M, M, C_dHdA_C, Gs+symsize(M)*p, C);
    for(int i=0; i<N/2; i++){
      for(int j=0; j<N/2; j++){
        double s = 0.0;
        for(int k=N/2; k<M; k++){
          s += Q[MPOSIF(i,k)] * C_dHdA_C[MPOSIF(k,j)] / (V[j]-V[k]);
        }
        dQdA[i*(N/2)+j] = 2.0*s;
      }
    }
    double g = 0.0;
    for(int i=core; i<N/2; i++){
      double * Bi = B+(N/2)*i;
      g += mx_vecdot_nosym(N/2, Bi, dQdA, Bi);
    }
    g1[p] -= g;
  }
  free(Q);

  double d1 = 0.0;
  for(int i=core; i<N/2; i++){
    d1 += 1.0 - D[i];
  }
  free(D);
  return d1;
}

