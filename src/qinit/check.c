#include "q.h"
#include "qinit.h"

static void regression(int nmol, mol_data * md,
    double * f1, double * f6, FILE * fo){

  int npar, atpar[NATOMS+1];

  {
    double a0[NATOMS],  a[NATOMS], fa[NATOMS];
    int bp[NATOMS], fixflag[NATOMS]={0}, fatpar[NATOMS+1];
    int nfpar, nparams = pars_default(2, a0, bp);
    npar = atpar_fill(nparams, nmol, fa, a, a0, bp, fixflag, &nfpar, fatpar, atpar, md, NULL);
  }

  double * C = calloc(nmol*npar*sizeof(double), 1);
  for(int i=0; i<nmol; i++){
    mol * m = md[i].m;
    for(int k=0; k<m->n; k++){
      int q = m->q[k];
      int j = atpar[q];
      C[j*nmol+i]++;
    }
  }

  double * D = malloc(symsize(npar)*sizeof(double));
  mx_sqr2(npar, nmol, D, C);

  double * Di = malloc(symsize(npar)*sizeof(double));
  int p = mx_inv_eigen(npar, Di, D, 1e-15);
  if(p){
    goto ret;
  }

  double * DiC = malloc(npar*nmol*sizeof(double));
  mx_symmultrectmx(npar, nmol, DiC, Di, C);

  double * w1 = malloc(npar*sizeof(double));
  double * w6 = malloc(npar*sizeof(double));
  mx_multvec(npar, nmol, w1, DiC, f1);
  mx_multvec(npar, nmol, w6, DiC, f6);

  fprintf(fo, "\nlinear regression\n\n");
  for(int i=1; i<=NATOMS; i++){
    int j = atpar[i];
    if(j>-1){
      fprintf(fo, "atom %4d : % 20.15lf % 20.15lf\n", i, w1[j], w6[j]);
    }
  }

  fprintf(fo, "\n");
  double f21 = 0.0;
  double f26 = 0.0;
  double err21 = 0.0;
  double err26 = 0.0;
  for(int i=0; i<nmol; i++){
    mol * m = md[i].m;
    double err1 = f1[i];
    double err6 = f6[i];
    for(int k=0; k<m->n; k++){
      int q = m->q[k];
      int j = atpar[q];
      err1 -= w1[j];
      err6 -= w6[j];
    }
    f21 += f1[i]*f1[i];
    f26 += f6[i]*f6[i];
    err21 += err1*err1;
    err26 += err6*err6;
  }
  fprintf(fo, "Σe^2/Σf^2 : % 20.15lf % 20.15lf\n", err21/f21, err26/f26);

  free(w1);
  free(w6);
  free(DiC);
ret:
  free(Di);
  free(D);
  free(C);
  return;
}

static void mols_check_print(int nmol,
    double * f1, double * f6, double * fs, int * dc,
    int pr, char ** fnames, FILE * fo){

  size_t l = 0;
  for(int i=0; i<nmol; i++){
    size_t tl = strlen(fnames[i]);
    l = MAX(l, tl);
  }
  fprintf(fo, "\n%*s       f(energy)            f(overlap)           max. overlap. d.\n", l, "");
  for(int i=0; i<nmol; i++){
    fprintf(fo, "%-*s  :    %-17.*lf    %-17.*lf    %-17.*lf   (%d)   #%d\n",
        l, fnames[i], pr, f1[i], pr, f6[i], pr, fs[i], dc[i], i);
    if(!(i%16)) fflush(fo);
  }
  return;
}

void mols_check(int nmol, mol_data * md,
    basis_type btype, basis_gc * bas, double * boys_array,
    int check, char ** fnames, FILE * fo){

  double * f1 = malloc(sizeof(double)*nmol);
  double * f6 = malloc(sizeof(double)*nmol);
  double * fs = malloc(sizeof(double)*nmol);
  int    * dc = malloc(sizeof(int)*nmol);

  for(int i=0; i<nmol; i++){
    int M = md[i].M;
    double * V = malloc((M*M+M+symsize(M))*sizeof(double));
    double * C = V + M;
    double * F = C + M*M;
    veccp(symsize(M), F, md[i].H);

    if(check==1){
      nozzle(md+i, F, btype, bas, boys_array);
    }

#if 0
    veccp(M*M, C, md[i].X0);
#else
    veccp(M*M, C, md[i].C0);
#endif
    mx_BHBt_sym(M, F, C);
    jacobi(F, C, V, M, 1e-15, 20, NULL);
    eigensort(M, V, C);
    f1[i] = calc_only_measure(1, C, md+i);
    f6[i] = calc_only_measure(6, C, md+i);
    fs[i] = measure_minoverlap(md[i].N, md[i].M, md[i].core, md[i].S, md[i].C0, C);
    dc[i] = measure_ortho     (md[i].N, md[i].M, md[i].core, md[i].S, md[i].C0, C);
    free(V);
  }

  mols_check_print(nmol, f1, f6, fs, dc, 15, fnames, fo);
  regression(nmol, md, f1, f6, fo);

  free(dc);
  free(fs);
  free(f1);
  free(f6);
  return;
}

