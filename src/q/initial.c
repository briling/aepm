#include "q.h"

static void mol_atcv_ij_add(int M, double * F,
    mol * m, atomo * ao, int * al, atomo * aox,
    basis_type btype, void * bas, double * boys_array){

  if(btype == BASIS_TYPE_GC){
    atcv_ij_add_gc(m, M, F, NULL, bas, ao, al, aox, boys_array);
  }
  else{
    atcv_ij_add(m, M, F, ao, aox, boys_array);
  }
  return;
}

void init_lb20_heff(int M, double * F,
    double * H, mol * m, atomo * ao, int * al,
    basis_type btype, void * bas, double * boys_array){

  /* Adds to matrix F the matrix of the effective potential minus bare-nuclear potential
   * Theor Chem Acc 139, 17 (2020)
   */

  atomo * aox = malloc(sizeof(atomo)*m->n);
  veccp(symsize(M), F, H);

  atcv_prep();
  atcv_fill(aox, m);
  mol_atcv_ij_add(M, F, m, ao, al, aox, btype, bas, boys_array);

  static double a[NELEMENTS+1] = {
    [  1 ...   2 ] = 1.0 /  3.0,
    [  3 ...   4 ] = 1.0 / 16.0,
    [  5 ...  10 ] = 1.0 /  3.0,
    [ 11 ...  12 ] = 1.0 / 32.0,
    [ 13 ...  18 ] = 1.0 /  8.0,
    [ 19 ...  20 ] = 1.0 / 32.0,
    [ 21 ...  30 ] = 1.0 /  6.0,
    [ 31 ...  36 ] = 1.0 / 12.0,
    [ 37 ...  38 ] = 1.0 / 32.0,
    [ 39 ...  48 ] = 1.0 /  8.0,
    [ 49 ...  54 ] = 1.0 / 12.0,
    [ 55 ...  70 ] = 1.0 / 32.0,
    [ 71 ...  86 ] = 1.0 / 12.0,
    [ 87 ... 102 ] = 1.0 / 32.0
  };
  double w[NELEMENTS+1];
  for(int i=1; i<=NELEMENTS; i++){
    charge1_norm(a+i, w+i);
  }
  atcv_nozzle_fill_all1(aox, a, w, m);
  mol_atcv_ij_add(M, F, m, ao, al, aox, btype, bas, boys_array);

  free(aox);
  return;
}

