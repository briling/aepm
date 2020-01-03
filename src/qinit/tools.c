#include "q.h"
#include "qinit.h"

void Hmod_only_solve(mol_data * md, double * V, double * C, double * F){
  int M = md->M;
#if 0
  veccp(M*M, C, md->X0);
#else
  veccp(M*M, C, md->C0);
#endif
  mx_BHBt_sym(M, F, C);
  jacobi(F, C, V, M, 1e-15, 20, NULL);
  eigensort(M, V, C);
  return;
}
