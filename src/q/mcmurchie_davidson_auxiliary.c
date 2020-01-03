#include "q.h"
#include "integrals.h"
#include "boys.h"
#include "gc.h"

static double RNLMj_r(int N, int L, int M, unsigned int j, double alpha, double r[3], double * F) __attribute__ ((unused));
static double RNLMj_r(int N, int L, int M, unsigned int j, double alpha, double r[3], double * F){
  if((N<0)||(L<0)||(M<0)){
    return 0.0;
  }
  else if((!N)&&(!L)&&(!M)){
    if(j>BOS_NMAX){
      GOTOHELL;
    }
    return intpow(-2.0*alpha, j) * F[j];
  }
  else if((!N)&&(!L)){
    return r[2]  * RNLMj_r(N, L, M-1, j+1, alpha, r, F) +
           (M-1) * RNLMj_r(N, L, M-2, j+1, alpha, r, F);
  }
  else if(!N){
    return r[1]  * RNLMj_r(N, L-1, M, j+1, alpha, r, F) +
           (L-1) * RNLMj_r(N, L-2, M, j+1, alpha, r, F);
  }
  else{
    return r[0]  * RNLMj_r(N-1, L, M, j+1, alpha, r, F) +
           (N-1) * RNLMj_r(N-2, L, M, j+1, alpha, r, F);
  }
}

double RNLMj(int N, int L, int M, double alpha, double r[3], double * F){

  if((!N)&&(!L)&&(!M)){
    return F[0];
  }

  double arr[BOS_NMAX+1][BOS_NMAX+1]={};
  int j0 = (N/2)+(M/2)+(L/2);
  double aj = intpow(-2.0*alpha, j0);
  for(int j=j0; j<=N+L+M; j++){
    arr[j][0] = aj * F[j];
    aj *= -2.0*alpha;
  }

  double x = r[0];
  double y = r[1];
  double z = r[2];

  for(int j=N+L+M-1; j>=(N/2)+(L/2); j--){
    int dj = j-(N/2)-(L/2);
    int M1 = MAX(2, M-2*dj);
    int M2 = j>N+L ? M+(N+L-j) : M;
    if( M <= 2*dj+1 ){
      arr[j][ 1] = z * arr[j+1][0];
    }
    for(int Mk=M1; Mk<=M2; Mk++){
      arr[j][Mk] = z * arr[j+1][Mk-1] + (Mk-1) * arr[j+1][Mk-2];
    }
  }

  for(int j=N+L-1; j>=N/2; j--){
    int dj = j-N/2;
    int L1 = MAX(2, L-2*dj);
    int L2 = j>N ? L+N-j : L;
    if( L <= 2*dj+1 ){
      arr[j][ 1+M] = y * arr[j+1][M];
    }
    for(int Lk=L1; Lk<=L2; Lk++){
      arr[j][Lk+M] = y * arr[j+1][Lk+M-1] + (Lk-1) * arr[j+1][Lk+M-2] ;
    }
  }

  for(int j=N-1; j>=0; j--){
    int N1 = MAX(N-2*j, 2);
    int N2 = N-j;
    if(N <= 2*j+1 ){
      arr[j][ 1+L+M] = x * arr[j+1][L+M];
    }
    for(int Nk=N1; Nk<=N2; Nk++){
      arr[j][Nk+L+M] = x * arr[j+1][Nk+L+M-1] + (Nk-1) * arr[j+1][Nk+L+M-2];
    }
  }

  return arr[0][N+L+M];
}

//////////////////////////////////////////////////////////////////////

typedef double arr_t[BOS_NMAX+1][BOS_NMAX+1];

static void RNLMj_0(int nlm, double alpha, double * F, arr_t arr){
  double aj = 1.0;
  for(int j=0; j<=nlm; j++){
    arr[j][0] = aj * F[j];
    aj *= -2.0*alpha;
  }
  return;
}

static void RNLMj_1(int nlm, int M, double z, arr_t arr){

  for(int j=nlm-1; j>=0; j--){
    int M0 = M-2*j;
    int M1 = MAX(M0, 2);
    int M2 = j>nlm-M ? nlm-j : M;
    if( M0 <= 1 ){
      arr[j][ 1] = z * arr[j+1][0];
    }
    for(int Mk=M1; Mk<=M2; Mk++){
      arr[j][Mk] = z * arr[j+1][Mk-1] + (Mk-1) * arr[j+1][Mk-2];
    }
  }
  return;
}

static void RNLMj_2(int nlm, int L, int M, double y, arr_t arr){

  for(int j=(nlm-M)-1; j>=0; j--){
    int L0 = L-2*j;
    int L1 = MAX(L0, 2);
    int L2 = j>(nlm-M)-L ? (nlm-M)-j : L;
    if( L0 <= 1 ){
      arr[j][ 1+M] = y * arr[j+1][M];
    }
    for(int Lk=L1; Lk<=L2; Lk++){
      arr[j][Lk+M] = y * arr[j+1][Lk+M-1] + (Lk-1) * arr[j+1][Lk+M-2] ;
    }
  }
  return;
}

static double RNLMj_3(int N, int LM, double x, arr_t arr){

  for(int j=N-1; j>=0; j--){
    int N0 = N-2*j;
    int N1 = MAX(N0, 2);
    int N2 = N-j;
    if( N0 <=1 ){
      arr[j][ 1+LM] = x * arr[j+1][LM];
    }
    for(int Nk=N1; Nk<=N2; Nk++){
      arr[j][Nk+LM] = x * arr[j+1][Nk+LM-1] + (Nk-1) * arr[j+1][Nk+LM-2];
    }
  }
  return arr[0][N+LM];
}

void rnlm_fill_adds(int nlm, double s, preijkl_t * prem, rnlm_t rnlm){
  if(nlm > BOS_NMAX) GOTOHELL;
  double x = prem->PQ[0];
  double y = prem->PQ[1];
  double z = prem->PQ[2];
  arr_t arr;
  RNLMj_0(nlm, prem->alpha, prem->F, arr);
  for(int M=0; M<=nlm; M++){
    RNLMj_1(nlm, M, z, arr);
    for(int L=0; L<=nlm-M; L++){
      RNLMj_2(nlm, L, M, y, arr);
      for(int N=0; N<=nlm-M-L; N++){
        rnlm[N][L][M] += s * RNLMj_3(N, L+M, x, arr);
      }
    }
  }
  return;
}

void rnlm_fill(int nlm, preijkl_t * prem, rnlm_t rnlm){
  if(nlm > BOS_NMAX) GOTOHELL;
  double x = prem->PQ[0];
  double y = prem->PQ[1];
  double z = prem->PQ[2];
  arr_t arr;
  RNLMj_0(nlm, prem->alpha, prem->F, arr);
  for(int M=0; M<=nlm; M++){
    RNLMj_1(nlm, M, z, arr);
    for(int L=0; L<=nlm-M; L++){
      RNLMj_2(nlm, L, M, y, arr);
      for(int N=0; N<=nlm-M-L; N++){
        rnlm[N][L][M] = RNLMj_3(N, L+M, x, arr);
      }
    }
  }
  return;
}

