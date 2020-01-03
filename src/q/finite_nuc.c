#include "q.h"

static const int masses[] = {
  #include <masses_int.h>
};

void finite_nuc_par(int q, double * alpha, double * omega){
  /* At. Data Nucl. Data Tables 67, 207--224 (1997) */
  int mass = masses[q];
  double a = 5700.0 + 8360.0 * pow(mass, 1./3.);
  a = 529177249.0 / a;
  a = 1.5 * a * a;
  double w = 0.0;
  charge1_norm(&a, &w);
  w *= q;
  *alpha = a;
  *omega = w;
  return;
}

