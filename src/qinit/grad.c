#include "q.h"
#include "qinit.h"

double calc_gradient_mol(
    int measure, int npar, double * g,
    mol_data * md, int atpar[],
    basis_type btype, basis_gc * bas, double * boys_array){

  int M = md->M;
  double * C = malloc((M*M+symsize(M)+M)*sizeof(double));
  double * V = C + M*M;
  double * F = V + M;
  double * Gs = calloc(npar*symsize(M), sizeof(double));
  veccp(symsize(M), F, md->H);
  nozzle_gradient(npar, md, F, Gs, atpar, btype, bas, boys_array);
  Hmod_only_solve(md, V, C, F);
  double E = calc_grad_measure(measure, npar, g, C, V, Gs, md);
  free(C);
  free(Gs);
  return E;
}

