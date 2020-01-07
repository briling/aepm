#include "q.h"
#include "qap.h"

static void calc_gradient_mol_num(int measure, int npar,
    double * g, mol_data * md, double * a, double * w,
    basis_type btype, basis_gc * bas, double * boys_array){

  int M = md->M;
  double * V = malloc((M*M+M+symsize(M))*sizeof(double));
  double * C = V + M;
  double * F = C + M*M;

  double d = 1e-5;
  for(int k=0; k<npar; k++){
    double a0 = a[k];
    for(int i=-1; i<=1; i+=2){
      a[k] = a0 + i*d;
      charge1_norm(a+k, w+k);
      veccp(symsize(M), F, md->H);
      cap(md, F, btype, bas, boys_array);
      Hmod_only_solve(md, V, C, F);
      g[k] += calc_only_measure(measure, C, md) * i;
    }
    a[k] = a0;
    charge1_norm(a+k, w+k);
    g[k] *= 0.5 / d;
  }
  free(V);
  return;
}

void test_grad_a(int nmol, int npar, mol_data * md,
    int atpar[NATOMS+1], double * a, double * w,
    basis_type btype, basis_gc * bas, double * boys_array){

  int measures[] = {MEASURE_E0, MEASURE_S0, MEASURE_S, MEASURE_E};
  for(int j=0; j<nmol; j++){
    printf("mol: %4d\n", j);
    for(int z=0; z<sizeof(measures)/sizeof(measures[0]); z++){
      printf("f = %s\n", measure_names[measures[z]]);
      double g1[NATOMS]={}, g2[NATOMS]={};
      calc_gradient_mol_num(measures[z], npar, g1, md+j, a, w,  btype, bas, boys_array);
      calc_gradient_mol    (measures[z], npar, g2, md+j, atpar, btype, bas, boys_array);
      for(int k=0; k<npar; k++){
        printf("%3d % 20.15lf % 20.15lf % 11.2e (%.0e)\n", k, g1[k], g2[k], g1[k]-g2[k], fabs((g1[k]-g2[k])/g1[k]));
      }
    }
  }
  return;
}
