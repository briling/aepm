#include "q.h"

#define EPSPRINT  1e-10

int * oneint_fill(int M, mol * m, void * bas, basis_type btype, atomo * ao,
                  double * S, double * H, double * Dxyz,
                  int finite_nuclei, atomo * u_rel,
                  double * boys_array, FILE * fo){
  int * al = NULL;
  if(btype == BASIS_TYPE_GC){
    al = al_gc(m, bas);
    oneint_fill_gc(m, bas, ao, al, S, H, Dxyz, finite_nuclei, u_rel, boys_array, fo);
  }
  else{
    oneint_fill_gen(m, M, ao, S, H, Dxyz, finite_nuclei, u_rel, boys_array);
  }
  return al;
}

void oneint_fill_gen(mol * m, int M, atomo * ao,
                     double * S, double * H, double * Dxyz,
                     int finite_nuclei, atomo * u_rel,
                     double * boys_array){
  for(int i=0; i<M; i++){
    for(int j=i; j<M; j++){
      int_ij(m, ao+i, ao+j, S+mpos(i,j), H+mpos(i,j), Dxyz?(Dxyz+3*mpos(i,j)):NULL, finite_nuclei, u_rel, boys_array);
    }
  }
  return;
}

void oneint_print(int M, double * S, double * H){
  for(int i=0; i<M; i++){
    for(int j=i; j<M; j++){
      double s = S[mpos(i,j)];
      double h = H[mpos(i,j)];
      if(fabs(s)>EPSPRINT || fabs(h)>EPSPRINT){
        printf("%02d %02d % .10lf % 16.10lf\n", i+1, j+1, S[mpos(i,j)], H[mpos(i,j)]);
      }
    }
  }
  return;
}

void atcv_ij_add(mol * m, int M, double * H, atomo * ao, atomo * aox, double * boys_array){
  for(int i=0; i<M; i++){
    for(int j=i; j<M; j++){
      H[mpos(i,j)] += int_ij0_sum(m, ao+i, ao+j, aox, boys_array);
    }
  }
  return;
}

void atcv_ij_add2(mol * m, int M, double * Hs,
    int npar, int atpar[],
    atomo * ao, atomo * aox, double * boys_array){
  for(int i=0; i<M; i++){
    for(int j=i; j<M; j++){
      for(int p=0; p<npar; p++){
        // TODO: use md->mpl
        double * H = Hs + symsize(M)*p;
        H[mpos(i,j)] += int_ij0_sum2(m, p, atpar, ao+i, ao+j, aox, boys_array);
      }
    }
  }
  return;
}

void atcv_ij_add_grad(mol * m, int M,
    double * Gs, int npar, int atpar[],
    atomo * ao, atomo * aox, double * boys_array){
  for(int i=0; i<M; i++){
    for(int j=i; j<M; j++){
      for(int p=0; p<npar; p++){
        // TODO: use md->mpl
        double * G = Gs + symsize(M)*p;
        G[mpos(i,j)] += int_ij0_sum_grad(m, p, atpar, ao+i, ao+j, aox, boys_array);
      }
    }
  }
  return;
}

