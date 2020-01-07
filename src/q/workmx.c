#include "q.h"
#include "matrix.h"

void D_fill(int N, int M, double * C, double * D, unsigned int w){
  for(int i=0; i<M; i++){
    for(int j=i; j<M; j++){
      double s = 0.0;
      for(int k=0; k<N/w; k++){
        s += C[k*M+i]*C[k*M+j];
      }
      D[mpos(i,j)] = s*w;
    }
  }
  return;
}

void F_restore(int M, double * F, double * S, double * V, double * C){
  double * SC = calloc(sizeof(double)*M*M,1);
  mx_multsymmx(M, SC, C, S);  // C^T * S
  mx_transp(M, SC);
  mx_BdiagBt(M, F, V, SC);
  free(SC);
  return;
}

void MO_proj(int m, int M, double * P, double * S, double * C){
  /* all occ. states: MO_proj(N/2, M, P, S, C);
   * valence  states: MO_proj(N/2-core, M, P, S, C+core*M);
   */
  double * SC = malloc(M*m*sizeof(double));
  for(int i=0; i<M; i++){
    for(int j=0; j<m; j++){
      double t = 0.0;
      for(int k=0; k<M; k++){
        t += S[MPOSIF(i,k)] * C[j*M+k];
      }
      SC[i*m+j] = t;
    }
  }
  mx_sqr2(M, m, P, SC);
  free(SC);
  return;
}

static int mx_invsqrt_wrap(int n, double * Aisqrt, double * A){
  double * d = malloc(sizeof(double)*(n*n+n));
  double * b = d + n;
  veccp(symsize(n), Aisqrt, A);
  int odd = mx_invsqrt1(Aisqrt, b, d, n, 1e-12, 1e-15, 20, NULL);
  free(d);
  return odd;
}

static void mx_BHBt_sym_aug0(
    unsigned int n, unsigned int m,
    double * h, double * b, double * t){
  /* H := BHB^T
   * H: n*n
   * B: m*n, m<=n -> augment with zeros
   */
  unsigned int i,j,k;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      double s = 0.0;
      for(k=0; k<=j; k++){
        s += b[i*n+k] * h[mpos(k,j)];
      }
      for(k=j+1; k<n; k++){
        s += b[i*n+k] * h[mpos(j,k)];
      }
      t[i*n+j] = s;
    }
  }
  for(i=0; i<m; i++){
    for(j=i; j<m; j++){
      double s = 0.0;
      for(k=0; k<n; k++){
        s += b[i*n+k] * t[j*n+k];
      }
      h[mpos(i,j)] = s;
    }
    for(j=m; j<n; j++){
      h[mpos(i,j)] = 0.0;
    }
  }
  for(i=m; i<n; i++){
    for(j=i; j<n; j++){
      h[mpos(i,j)] = 0.0;
    }
  }
  return;
}

int C_maxoverlap(int n, int n1, int M, double * Asym,
                 double * C, double * C0, double * S){
  /* See Int. J. Quantum Chem. 111 (2011), 2851 (Sec. 2.3):
   *
   * <C0|S|C'> = <C'|S|C0> (maximum overlap), when
   * C' = C * U,
   * U  = (At*A)^(-1/2) * A,
   * A  = C0^T * S * C.
   * Thus <C0|S|C'> = At * (At*A)^(-1/2) * A,
   * and this function returns it.
   *
   * M  is the dimension of the AO space
   * n  is the dimension of the occ. MO subspace (C, C0)
   * n1 is the number of the first C0 vector
   *    we take into account, thus real dim(C0) is m=n-n1
   */

  int m = n-n1;

  double * A   = malloc( (2*n*m+symsize(n)) * sizeof(double));
  double * At  = A  + n*m;
  double * AtA = At + n*m;

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      double csc = mx_vecdot(M, C+M*i, S, C0+M*(j+n1));
      A [j*n+i] = csc;
      At[i*m+j] = csc;
    }
  }

  mx_sqr2(n, m, AtA, At);
  int n2 = mx_invsqrt_wrap(n, Asym, AtA);   // (At*A)^(-1/2)
  // n2 has to be == n1
  mx_BHBt_sym_aug0(n, m, Asym, A, At);

  free(A);
  return n2;
}

