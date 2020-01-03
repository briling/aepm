#include <pthread.h>
#include "q.h"
#include "qinit.h"

int mol2proc(int nproc, int nmol, pthread_t ** threads, thread_arg ** args, mol_data * md){
  if(nproc > nmol){
    nproc = nmol;
  }
  *args    = malloc(sizeof(thread_arg)*nproc);
  *threads = malloc(sizeof(pthread_t)*nproc);
  int mpp  = nmol/nproc;
  int rmpp = nmol%nproc;
  for(int i=0; i<nproc; i++){
    (*args)[i].md    = ( i<rmpp ? i*(mpp+1) : rmpp*(mpp+1) + (i-rmpp)*mpp ) + md;
    (*args)[i].nmol  = i<rmpp ? mpp+1 : mpp;
  }
  return nproc;
}

void * mol_atcv_ll(void * arg){
  thread_arg * ta = arg;
  for(int i=0; i<ta->nmol; i++){
    nozzle(ta->md+i, ta->md[i].H, ta->btype, ta->bas, ta->boys_array);
  }
  return NULL;
}

void * mol_1el_ll(void * arg){
  thread_arg * ta = arg;
  for(int i=0; i<ta->nmol; i++){
    mol_1el(ta->md+i,
            ta->finite_nuclei, ta->urelconst,
            ta->btype, ta->bas, ta->boys_array);
  }
  return NULL;
}

void * calc_gradient_mol_ll(void * arg){
  thread_arg * ta = arg;
  for(int i=0; i<ta->nmol; i++){
    *(ta->e) += calc_gradient_mol(
        ta->measure, ta->npar, ta->g, ta->md+i, ta->atpar,
        ta->btype, ta->bas, ta->boys_array);
  }
  return NULL;
}

