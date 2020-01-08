#include "q.h"
#include "boys.h"

#define EPS 1e-15
#define STRLEN 256

  /*  assumptions:
   *  1.  there is only one set in the file
   *  2.  elements are ordered
   *  3.  angular momenta are ordered
   *  4.  no empty lines
   */

static const char shellname[] = { 's', 'p', 'd', 'f', 'g', 'h' };

static void g_norm(void * bas_p, basis_type bt){
  if(bt == BASIS_TYPE_GENERAL){
    basis * bas = bas_p;
    for(int n=1; n<=NELEMENTS; n++){
      for(int i=bas->lsto[n-1]; i<bas->lsto[n]; i++){
        int l = bas->l[i];
        for(int j=bas->lgto[i]; j<bas->lgto[i+1]; j++){
          double a = bas->a[j];
          double s = int_s00(l, a);
          bas->w[j] /= sqrt(s);
        }
      }
    }
  }
  else if(bt == BASIS_TYPE_GC){
    basis_gc * bas = bas_p;
    for(int n=1; n<=NELEMENTS; n++){
      for(int ll=bas->ll[n-1]; ll<bas->ll[n]; ll++){
        int l = bas->l[ll];
        for(int i=bas->lp[ll]; i<bas->lp[ll+1]; i++){
          double a  = bas->a[i];
          double sc = 1.0/sqrt(int_s00(l, a));
          for(int j=0; j<bas->nc[ll]; j++){
            bas->c[
              +bas->lc[ll]
              +bas->np[ll]*j
              +(i-bas->lp[ll])
              ] *= sc;
          }
        }
      }
    }
  }
  return;
}

static void ao_norm(void * bas_p, basis_type bt){
  if(bt == BASIS_TYPE_GENERAL){
    basis * bas = bas_p;
    for(int q=1; q<=NELEMENTS; q++){
      for(int i=bas->lsto[q-1]; i<bas->lsto[q]; i++){
        int l = bas->l[i];
        atomo ao = {
          .k = 0,
          .m = 0,
          .l = l,
          .a = bas->a + bas->lgto[i],
          .w = bas->w + bas->lgto[i],
          .ng = bas->lgto[i+1] - bas->lgto[i],
        };
        double sc = 1.0/sqrt(int_ij_overlap_self(&ao));
        vecscal(ao.ng, ao.w, sc);
      }
    }
  }
  else if(bt == BASIS_TYPE_GC){
    basis_gc * bas = bas_p;

    for(int q=1; q<=NELEMENTS; q++){
      for(int ll=bas->ll[q-1]; ll<bas->ll[q]; ll++){
        int l = bas->l[ll];
        int nc = bas->nc[ll];
        for(int j=0; j<nc; j++){
          atomo ao = {
            .k = 0,
            .m = 0,
            .l = l,
            .a = bas->a + bas->lp[ll],
            .w = bas->c + bas->lc[ll] + bas->np[ll]*j,
            .ng = count_primitives_gc(bas, ll, j),
          };
          double sc = 1.0/sqrt(int_ij_overlap_self(&ao));
          vecscal(ao.ng, ao.w, sc);
        }
      }
    }
  }
  return;
}

void basel_print_gen(int n, basis * bas, int c, const char s[], FILE * f){
  for(int i=bas->lsto[n-1]; i<bas->lsto[n]; i++){
    int l = bas->l[i];
    if(c){
      fprintf(f, "%s %c %2d\n",  s, shellname[l], bas->lgto[i+1] - bas->lgto[i]);
    }
    else{
      fprintf(f, "%s%1d %2d\n", s, l,     bas->lgto[i+1] - bas->lgto[i]);
    }
    for(int j=bas->lgto[i]; j<bas->lgto[i+1]; j++){
      double w  = bas->w[j];
      double a  = bas->a[j];
      double s2 = int_s00(l, a);
      fprintf(f, "%s %.10e % .10e\n", s, a, w*sqrt(s2));
    }
  }
  return;
}

int count_primitives_gc(basis_gc * bas, int ll, int j __attribute__ ((unused))){
  int n  = bas->lp[ll+1]-bas->lp[ll];
#if 0
  int ng = 0;
  for(int i=0; i<n; i++){
    double c = bas->c[
      +bas->lc[ll]
        +bas->np[ll]*j
        +i];
    if(fabs(c)>EPS){
      ng++;
    }
  }
  return ng;
#else
  return n;
#endif
}

void basel_print_gc2gen(int n, basis_gc * bas, int c, const char s[], FILE * f){
  for(int ll=bas->ll[n-1]; ll<bas->ll[n]; ll++){
    int l = bas->l[ll];

    for(int j=0; j<bas->nc[ll]; j++){

      int np = 0;
      for(int i=bas->lp[ll]; i<bas->lp[ll+1]; i++){
        double c = bas->c[ +bas->lc[ll] +bas->np[ll]*j +(i-bas->lp[ll]) ];
        if(fabs(c)>EPS) np++;
      }

      if(c){
        fprintf(f, "%s %c %2d\n",  s, shellname[l], np);
      }
      else{
        fprintf(f, "%s%1d %2d\n", s, l, np);
      }

      for(int i=bas->lp[ll]; i<bas->lp[ll+1]; i++){
        double c = bas->c[ +bas->lc[ll] +bas->np[ll]*j +(i-bas->lp[ll]) ];
        if(fabs(c)>EPS){
          double a  = bas->a[i];
          double s2 = int_s00(l, a);
          fprintf(f, "%s %.10e % .10e\n", s, a, c*sqrt(s2));
        }
      }


    }
  }
  return;
}

static void bas_print_gen(basis * bas, const char s[], FILE * f){
  for(int n=1; n<=NELEMENTS; n++){
    int nm = bas->lsto[n]-bas->lsto[n-1];
    if(nm){
      fprintf(f, "%s\t%3d %2d\n", s, n, nm);
    }
    basel_print_gen(n, bas, 0, s, f);
  }
  return;
}

static void bas_print_gc(basis_gc * bas, const char s[], FILE * f){
  for(int n=1; n<=NELEMENTS; n++){
    int ll0 = bas->ll[n-1];
    int ll1 = bas->ll[n];
    if(ll1-ll0){
      fprintf(f, "%s%9d %d\n", s, n, ll1-ll0-1);
      for(int ll=ll0; ll<ll1; ll++){
        int l  = bas->l[ll];
        int p0 = bas->lp[ll];
        int np = bas->lp[ll+1]-bas->lp[ll];
        fprintf(f, "%s%2d %d\n", s, np, bas->nc[ll]);
        for(int i=0; i<np; i++){
          double a  = bas->a[i+p0];
          double sc = sqrt(int_s00(l, a));
          fprintf(f, "%s%+.15e  ", s, a);
          for(int j=0; j<bas->nc[ll]; j++){
            double c = bas->c[
                +bas->lc[ll]
                +bas->nc[ll]*i
                +j];
            fprintf(f, "%+.14e ", c*sc);
          }
          fprintf(f, "\n");
        }
      }
    }
  }
  return;
}

void bas_print(void * bas, basis_type btype, const char s[], FILE * f){
  fprintf(f, "%s$basis\n", s);
  fprintf(f, "%stype=%s\n", s, btype == BASIS_TYPE_GENERAL ? "general" : "gc");
  char * name = (btype == BASIS_TYPE_GENERAL ? ((basis * )bas)->name : ((basis_gc * )bas)->name);
  if(name[0]){
    fprintf(f, "%sdefault=%s\n", s, name);
    fprintf(f, "%sset=%s\n", s, name);
  }
  btype == BASIS_TYPE_GENERAL ? bas_print_gen(bas, s, f) : bas_print_gc (bas, s, f) ;
  fprintf(f, "%s$end\n", s);
  fprintf(f, "\n");
  return;
}

static long findbasis(char type[STRLEN], FILE * fb){
  char s[STRLEN];
  while (1){
    if (!fgets(s, sizeof(s), fb)) {
      return -1;
    }
    if(!strcmp(s, "$basis\n")){
      break;
    }
  }
  if(fscanf(fb, " type=%255s", type) != 1){
    return -1;
  };
  if(fscanf(fb, " default=%255s", s) == 1){
    /**/
  };
  return ftell(fb);
}

static int bas_read_gen(FILE * fb, basis * bas, int * Nsto, int * Ng){

  char   s[STRLEN];
  int    n,nm,l,lm;
  double a,w;

  int nsto = 0;
  int ng   = 0;

  if(fscanf(fb, " set=%255s", s) == 1){
    if(bas){
      strncpy(bas->name, s, sizeof(bas->name));
    }
  }

  while(fscanf(fb, "%d%d", &n, &nm) == 2){ /* element , number of functions */
    if(n > NELEMENTS){
      PRINT_WARN(" N=%d > %d=Nmax!\n", n, NELEMENTS);
      goto hell;
    }
    if(bas){
      if(bas->lsto[n]){
        goto hell;
      }
      bas->lsto[n] = nsto + nm;
    }
    for(int i=0; i<nm; i++){
      if (fscanf(fb, "%d%d", &l, &lm) != 2){ /* angular momentum , number of gto's */
        goto hell;
      }
      if(l>L_MAX){
        PRINT_WARN(" L=%d > %d=Lmax!\n", l, L_MAX);
        goto hell;
      }
      if(bas){
        bas->l[i+nsto] = l;
        bas->lgto[i+nsto  ] = ng;
        bas->lgto[i+nsto+1] = ng + lm;
      }
      for(int j=ng; j<ng+lm; j++){
        if (fscanf(fb, "%lf%lf", bas ? bas->a+j : &a, bas ? bas->w+j : &w) != 2){
          goto hell;
        }
      }
      ng += lm;
    }
    nsto += nm;
  }
  if ( (!fgets(s, sizeof(s), fb)) || strcmp(s, "$end\n")){
    return -1;
  }

  if(bas){
    for(int i=1; i<=NELEMENTS; i++){
      if(!bas->lsto[i]){
        bas->lsto[i] = bas->lsto[i-1];
      }
    }
  }

  *Ng   = ng;
  *Nsto = nsto;
  return 0;

hell:
  return -1;
}

static int bas_read_gc(FILE * fb, basis_gc * bas, int * Sum_l, int * Sum_p, int * Sum_c){

  char   s[STRLEN];

  int sum_l = 0;
  int sum_p = 0;
  int sum_c = 0;

  int n, lmax, np, nc;
  double a, c;

  if(fscanf(fb, " set=%255s", s) == 1){
    if(bas){
      strncpy(bas->name, s, sizeof(bas->name));
    }
  }

  while(fscanf(fb, "%d%d", &n, &lmax) == 2){ /* element , max angular momentum */
    if(n > NELEMENTS){
      PRINT_WARN(" N=%d > %d=Nmax!\n", n, NELEMENTS);
      goto hell;
    }
    if(lmax>L_MAX){
      PRINT_WARN(" L=%d > %d=Lmax!\n", lmax, L_MAX);
      goto hell;
    }

    if(bas){
      if(bas->ll[n]){
        goto hell;
      }
      bas->ll[n] = sum_l+lmax+1;
    }

    for(int l=0; l<=lmax; l++){

      if(fscanf(fb, "%d%d", &np, &nc) != 2){ /* number of primitives , number of contracted functions */
        goto hell;
      }
      if(bas){
        bas->l [sum_l+l] = l;
        bas->np[sum_l+l] = np;
        bas->nc[sum_l+l] = nc;
        bas->lc[sum_l+l] = sum_c;
        bas->lp[sum_l+l] = sum_p;
        bas->lp[sum_l+l+1] = sum_p + np;
      }

      for(int i=0; i<np; i++){
        if(fscanf(fb, "%lf", bas ? bas->a+sum_p+i : &a) != 1){
          goto hell;
        }
        for(int j=0; j<nc; j++){
          if(fscanf(fb, "%lf", bas ? bas->c+sum_c+np*j+i : &c) != 1){
            goto hell;
          }
        }
      }
      sum_p += np;
      sum_c += np * nc;
    }
    sum_l += lmax+1;
  }
  if ( (!fgets(s, sizeof(s), fb)) || strcmp(s, "$end\n")){
    return -1;
  }

  if(bas){
    for(int i=1; i<=NELEMENTS; i++){
      if(!bas->ll[i]){
        bas->ll[i] = bas->ll[i-1];
      }
    }
  }

  *Sum_l = sum_l;
  *Sum_p = sum_p;
  *Sum_c = sum_c;
  return 0;

hell:
  return -1;
}

void * bas_read(FILE * fb, basis_type * btype){

  char type[STRLEN] = {};
  rewind(fb);
  long bstart = findbasis(type, fb);

  if(bstart == -1){
    return NULL;
  }

  if(!strcmp(type, "general")){

    *btype = BASIS_TYPE_GENERAL;

    int nsto, ng;

    if( bas_read_gen(fb, NULL, &nsto, &ng)){
      return NULL;
    }

    size_t bsize = sizeof(basis)
      + sizeof(double)*ng     // double * a
      + sizeof(double)*ng     // double * w
      + sizeof(int)*nsto      // int * l
      + sizeof(int)*(nsto+1); // int * lgto
    basis * bas = malloc(bsize);
    if(!bas){
      return NULL;
    }
    memset(bas->lsto, 0, (NELEMENTS+1)*sizeof(int));
    memset(bas->name, 0, sizeof(bas->name));
    bas->a = (double *)(bas+1);
    bas->w = bas->a + ng;
    bas->l = (int *)(bas->w + ng);
    bas->lgto = bas->l+nsto;

    fseek(fb, bstart, SEEK_SET);
    if(bas_read_gen(fb, bas, &nsto, &ng)){
      free(bas);
      return NULL;
    }

    g_norm(bas, *btype);
    ao_norm(bas, *btype);

#if 0
    for(int i=0; i<=NELEMENTS; i++) printf("%d %d\n", i, bas->lsto[i]);
    for(int i=0; i<nsto; i++) printf("%d %d %d\n", i, bas->l[i], bas->lgto[i+1] - bas->lgto[i]);
    for(int i=0; i<ng; i++) printf("%d %lf %lf\n", i, bas->a[i], bas->w[i]);
#endif
    return bas;
  }

  else if(!strcmp(type, "gc")){

    *btype = BASIS_TYPE_GC;

    int nl, np, nc;

    if( bas_read_gc(fb, NULL, &nl, &np, &nc)){
      return NULL;
    }

    size_t bsize = sizeof(basis_gc)
      + sizeof(int)*nl     // int * l
      + sizeof(int)*nl     // int * nc
      + sizeof(int)*nl     // int * np
      + sizeof(int)*nl     // int * lc
      + sizeof(int)*(nl+1) // int * lp
      + sizeof(double)*np  // double * a
      + sizeof(double)*nc; // double * c

    basis_gc * bas = malloc(bsize);
    if(!bas){
      return NULL;
    }
    memset(bas->ll, 0, (NELEMENTS+1)*sizeof(int));
    memset(bas->name, 0, sizeof(bas->name));
    bas->a  = (double *)(bas+1);
    bas->c  = bas->a + np;
    bas->l  = (int *)(bas->c + nc);
    bas->np = bas->l  + nl;
    bas->nc = bas->np + nl;
    bas->lc = bas->nc + nl;
    bas->lp = bas->lc + nl;

    fseek(fb, bstart, SEEK_SET);
    if( bas_read_gc(fb, bas, &nl, &np, &nc)){
      free(bas);
      return NULL;
    }
    g_norm(bas, *btype);
    ao_norm(bas, *btype);

    return bas;
  }

  else{
    return NULL;
  }
}

