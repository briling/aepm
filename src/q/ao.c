#include "q.h"

static int ao_fill_gen(mol * m, basis * bas, atomo * ao){
  int M = 0;
  for(int k=0; k<m->n; k++){
    int q = m->q[k];
    int n = bas->lsto[q]-bas->lsto[q-1];
    if(!n){
      PRINT_ERR(" no basis functions for %d !\n", q);
      GOTOHELL;
    }
    for(int i=bas->lsto[q-1]; i<bas->lsto[q]; i++){
      int l = bas->l[i];
      if(ao){
        for(int m=-l; m<=l; m++){
          ao[M+l+m].m  = m;
          ao[M+l+m].l  = l;
          ao[M+l+m].k  = k;
          ao[M+l+m].a  = bas->a + bas->lgto[i] ;
          ao[M+l+m].w  = bas->w + bas->lgto[i] ;
          ao[M+l+m].ng = bas->lgto[i+1] - bas->lgto[i];
        }
      }
      M += 2*l+1;
    }
  }
  return M;
}

static int ao_fill_gc(mol * m, basis_gc * bas, atomo * ao){
  int M = 0;
  for(int k=0; k<m->n; k++){
    int q = m->q[k];
    int n = bas->ll[q]-bas->ll[q-1];
    if(!n){
      PRINT_ERR(" no basis functions for %d !\n", q);
      GOTOHELL;
    }
    for(int ll=bas->ll[q-1]; ll<bas->ll[q]; ll++){
      int l = bas->l[ll];
      int nc = bas->nc[ll];
      for(int j=0; j<nc; j++){
        if(ao){
          double * a = bas->a + bas->lp[ll];
          double * w = bas->c + bas->lc[ll] + bas->np[ll]*j;
          int ng = count_primitives_gc(bas, ll, j);
          for(int m=-l; m<=l; m++){
            ao[M+l+m].m  = m;
            ao[M+l+m].l  = l;
            ao[M+l+m].k  = k;
            ao[M+l+m].a  = a;
            ao[M+l+m].w  = w;
            ao[M+l+m].ng = ng;
          }
        }
        M += 2*l+1;
      }
    }
  }
  return M;
}

void ao_print(atomo * ao, int M, FILE * f){
  for(int i=0; i<M; i++){
    fprintf(f, "%4d   l =%2d   m =%2d   na =%3d   f =%p\n", i,  ao[i].l, ao[i].m, ao[i].k, ao[i].a);
    for(int j=0; j<ao[i].ng; j++){
      fprintf(f, "          % lf   % lf\n", ao[i].a[j], ao[i].w[j]);
    }
  }
  fprintf(f, "\n");
}

atomo * ao_fill(mol * m, void * bas, basis_type btype, int * M){
  atomo * ao;
  if(btype == BASIS_TYPE_GENERAL){
    *M  = ao_fill_gen(m, bas, NULL);
    ao = malloc(*M*sizeof(atomo));
    ao_fill_gen(m, bas, ao);
  }
  else{
    *M  = ao_fill_gc(m, bas, NULL);
    ao = malloc(*M*sizeof(atomo));
    ao_fill_gc(m, bas, ao);
  }
  return ao;
}

