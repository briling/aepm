#include "q.h"
#include "qinit.h"
#include "task_q.h"

void sol_save(mol_data * md, const char fname[],
    basis_type btype, void * bas, double * boys_array){

  int M = md->M;
  double * V = malloc((M*M+M+symsize(M))*sizeof(double));
  double * C = V + M;
  double * F = C + M*M;
  veccp(symsize(M), F, md->H);
  nozzle(md, F, btype, bas, boys_array);
  veccp(M*M, C, md->X0);
  mx_BHBt_sym(M, F, C);
  jacobi(F, C, V, M, 1e-15, 20, NULL);
  eigensort(M, V, C);

  filename vname;
  change_suffix(vname, fname, ".new.vec", sizeof(vname));
  pvec_write(M, md->ao, V, C, V, C, vname);

  free(V);
  return;
}

