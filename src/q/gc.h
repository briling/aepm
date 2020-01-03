#include "boys.h"
#include "integrals.h"

typedef double rnlm_t[BOS_NMAX+1][BOS_NMAX+1][BOS_NMAX+1];

typedef struct{
  double ap;
  double p21;
  double P[3], PA[3], PB[3];
  double E;
} preij_t;

typedef struct{
  double PQ[3];
  double alpha;
  double F[BOS_NMAX+1];
  double fact;
} preijkl_t;

typedef struct {
  double s;
  double t;
  double v;
  double d[3];
} shd_t;

preij_t preij(double * r, double a1, double a2, int k1, int k2);
preijkl_t preij0 (double * Q, double aq, int lsum, preij_t * pP, double * boys_array);

double int_ij0_prim(int l[2], int m[2],
    double D1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    double E1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    double F1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    rnlm_t rnlm);

void rnlm_fill_adds(int nlm, double s, preijkl_t * prem, rnlm_t rnlm);
void rnlm_fill(int nlm, preijkl_t * prem, rnlm_t rnlm);

void rnlm_fill_nucattr_point(int lsum,
    rnlm_t rnlm, double ap, double P[3],
    mol * m, double * boys_array);
void rnlm_fill_nucattr_finite(int lsum,
    rnlm_t rnlm, double ap, double P[3],
    mol * m, double * boys_array);
double add_prim_nucattr(
    int N, int L, int M,
    double d[L_MAX*2+1],
    double e[L_MAX*2+1],
    double f[L_MAX*2+1],
    rnlm_t rnlm);
void kinetic_scale_prim(
    int l1, int l2, double a, double b,
    preij_t pre, shd_t * pint,
    double r1[3], double r2[3], mol * m, atomo * u_rel);

static inline size_t ij_pint_memory(mol * m, basis_gc * bas){
  size_t np_max = 0;
  for(int k=0; k<m->n; k++){
    int q = m->q[k];
    for(int ll=bas->ll[q-1]; ll<bas->ll[q]; ll++){
      int np = bas->np[ll];
      np_max = MAX(np, np_max);
    }
  }
  return np_max * np_max;
}

static inline size_t ij_pint_memory1(mol * m, basis_gc * bas){
  size_t nl_np_max = 0;
  for(int k1=0; k1<m->n; k1++){
    int q1 = m->q[k1];
    for(int ll1=bas->ll[q1-1]; ll1<bas->ll[q1]; ll1++){
      int l1  = bas->l[ll1];
      int np1 = bas->np[ll1];
      size_t t = np1 * (2*l1+1);
      nl_np_max = MAX(t, nl_np_max);
    }
  }
  return nl_np_max * nl_np_max;
}

static inline void progress(int i, int M, int n, FILE * fo){
  if(fo){
    fprintf(fo, "%3d ", i+1);
    if( (!((i+1)%n)) || i==M-1){
      fprintf(fo, "\n");
    }
    fflush(fo);
  }
  return;
}

