#include "q.h"

#define STRLEN 256

atomo * u_rel_fill(mol * m, urelconst_t * urelconst){
  atomo * u_rel = malloc(sizeof(atomo)*m->n);
  for(int k=0; k<m->n; k++){
    int q  = m->q[k];
    int i0 = urelconst->list[q-1];
    int ng = urelconst->list[q] - i0;
    u_rel[k].l = u_rel[k].m = 0;
    u_rel[k].k = k;
    u_rel[k].ng = ng;
    u_rel[k].a = urelconst->a + i0;
    u_rel[k].w = urelconst->w + i0;
  }
  return u_rel;
}

static long findkinetic(FILE * fb){
  char s[STRLEN];
  while (1){
    if (!fgets(s, sizeof(s), fb)) {
      return -1;
    }
    if(!strcmp(s, "$kinetic\n")){
      break;
    }
  }
  return ftell(fb);
}

static int kinetic_read(urelconst_t * urel, FILE * fb){

  char   s[STRLEN];
  double a, w;
  int q, nq, n = 0;

  if(urel){
    for(int i=0; i<NELEMENTS+1; i++){
      urel->list[i] = 0;
    }
  }

  while(fscanf(fb, "%d%d", &q, &nq) == 2){
    if(urel && q<=NELEMENTS){
      urel->list[q] = nq;
    }
    for(int i=0; i<nq; i++){
      if (fscanf(fb, "%lf%lf", &w, &a) != 2){
        goto hell;
      }
      if(q<=NELEMENTS){
        if(urel){
          urel->a[n] = a;
          urel->w[n] = w;
        }
        n++;
      }
    }
  }
  if ( (!fgets(s, sizeof(s), fb)) || strcmp(s, "$end\n")){
    return -1;
  }

  if(urel){
    for(int i=1; i<NELEMENTS+1; i++){
      urel->list[i] += urel->list[i-1];
    }
  }

  return n;
hell:
  return -1;
}

void urelconst_print(urelconst_t * urel, FILE * fo){
  if(!urel){
    return;
  }
  fprintf(fo, "$kinetic\n");
  for(int q=1; q<=NELEMENTS; q++){
    int nq = urel->list[q]-urel->list[q-1];
    fprintf(fo, "%03d %d\n", q, nq);
    for(int i=urel->list[q-1]; i<urel->list[q]; i++){
      fprintf(fo, "%+23.15e %+23.15e\n", urel->w[i], urel->a[i]);
    }
  }
  fprintf(fo, "$end\n");
  return;
}

urelconst_t * urelconst_read(FILE * fb){

  long kstart = findkinetic(fb);
  if(kstart == -1){
    return NULL;
  }

  int n = kinetic_read(NULL, fb);
  if(n<0){
    return NULL;
  }

  urelconst_t * urel = malloc(sizeof(urelconst_t) + sizeof(double)*n*2);
  urel->a = (double *)(urel + 1);
  urel->w = urel->a + n;

  rewind(fb);
  fseek(fb, kstart, SEEK_SET);
  kinetic_read(urel, fb);

  return urel;
}

