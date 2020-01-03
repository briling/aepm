#include "q.h"

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

