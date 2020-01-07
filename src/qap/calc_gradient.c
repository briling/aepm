#include "q.h"
#include "qap.h"
#include <pthread.h>

extern int nproc;
extern pthread_t * threads;

double calc_gradient(int npar, double * g, thread_arg * ta){

  double E = 0.0;
  vecset(npar, g, 0.0);

  double * gs = calloc(sizeof(double)*npar*nproc,1);
  double * es = calloc(sizeof(double)*nproc,1);
  for(int i=0; i<nproc; i++){
    ta[i].g = gs+i*npar;
    ta[i].e = es+i;
    pthread_create(threads+i, NULL, calc_gradient_mol_ll, ta+i);
  }
  for(int i=0; i<nproc; i++){
    pthread_join(threads[i], NULL);
    E += es[i];
    vecadd(npar, g, gs+i*npar);
  }
  free(gs);
  free(es);
  return E;
}
