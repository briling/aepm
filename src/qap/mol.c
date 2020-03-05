#include "q.h"
#include "qap.h"
#include "task_q.h"

int mol_init(mol_data * md, basis_type btype, void * bas, char * fname){

  mol * m;
  FILE * fm = fopen(fname, "r");

  if( !fm || !(m = mol_read(fm)) ){
    PRINT_WARN("\tmol?!\n");
    return 1;
  }
  fclose(fm);

  int M; atomo * ao = ao_fill(m, bas, btype, &M);
  int N = elnumber(m);
  int core = corepairs_simple(m);
  if( m->mult-1 || N%2 ) {
    PRINT_WARN("\tN = %d, mult = %d !\n", N, m->mult);
    return 1;
  }
  if( N/2 > M) {
    PRINT_WARN("\tN = %d, M = %d !\n", N, M);
    return 1;
  }
  if(core>N/2){
    PRINT_WARN("\tN = %d, core = %d !\n", N, core);
    return 1;
  }

  md->M = M;
  md->N = N;
  md->core = core;
  md->m = m;
  md->ao = ao;
  md->al  = NULL;
  md->H   = NULL;
  md->aox = NULL;
  md->V0  = NULL;
  return 0;
}

int remove_oddmols(int nmol, mol_data * md, int * atpar, char ** fnames){
  int * tmp = calloc(sizeof(int)*nmol,1);
  for(int i=0; i<nmol; i++){
    for(int k=0; k<md[i].m->n; k++){
      int q = md[i].m->q[k];
      if(atpar[q]>=0){
        tmp[i] = 1;
        break;
      }
    }
  }
  int newnmol = nmol;
  for(int i=nmol-1; i>=0; i--){
    if(!tmp[i]){
      free(md[i].m);
      free(md[i].ao);
      for(int j=i+1; j<newnmol; j++){
        md[j-1] = md[j];
        fnames[j-1] = fnames[j];
      }
      newnmol--;
    }
  }
  free(tmp);
  return newnmol;
}

int * molparlist(int nmol, int npar, mol_data * md, int * atpar){
  int * mpl = calloc(sizeof(int)*nmol*npar,1);
  for(int i=0; i<nmol; i++){
    for(int k=0; k<md[i].m->n; k++){
      int q = md[i].m->q[k];
      int aq = atpar[q];
      if(aq>=0){
        mpl[i*npar+aq] = 1;
      }
    }
    md[i].mpl = mpl + i*npar;
  }
#if 0
  for(int i=0; i<nmol; i++){
    printf("mol# %d:  ", i);
    for(int j=0; j<npar; j++){
      if(md[i].mpl[j]){
        printf("%d ", j);
      }
    }
    printf("\n");
  }
#endif
  return mpl;
}

void mol_1el(mol_data * md,
    int finite_nuclei, urelconst_t * urelconst,
    basis_type btype, void * bas, double * boys_array){
  int M = md->M;
  md->H  = calloc((2*symsize(M)+M*M)*sizeof(double), 1);
  md->S  = md->H + symsize(M);
  md->X0 = md->S + symsize(M);

  atomo * u_rel = NULL;
  if(urelconst){
    u_rel = u_rel_fill(md->m, urelconst);
  }
  md->al = oneint_fill(M, md->m, bas, btype, md->ao, md->S, md->H, NULL, finite_nuclei, u_rel, boys_array, NULL);
  free(u_rel);

  double * wS = malloc((symsize(M)+M)*sizeof(double));
  veccp(symsize(M), wS, md->S);
  canorth(wS, md->X0, wS+symsize(M), M, 1e-15, 20, NULL);
  free(wS);
  return;
}

void mol_readvec(mol_data * md, char * fname, FILE * fo){

  filename vname;
  change_suffix(vname, fname, ".vec", sizeof(vname));

  int M = md->M;
  int n = md->N/2;
  int core = md->core;

  md->V0 = malloc(sizeof(double)*(M+M*M+2*symsize(M)));
  md->C0 = md->V0 + M;
  md->F0 = md->C0 + M*M;
  md->P0 = md->F0 + symsize(M);
  if(!pvec_read(M, md->ao, md->V0, md->C0, NULL, NULL, vname)){
    PRINT_ERR("Cannot read vectors from file '%s'\n", vname);
    GOTOHELL;
  }

  F_restore(M, md->F0, md->S, md->V0, md->C0);
  MO_proj(n-core, M, md->P0, md->S, md->C0+core*M);

  if(core){
    fprintf(fo, " | %+.1e", md->V0[core-1]);
  }
  else{
    fprintf(fo, " |         ");
  }
  fprintf(fo, " | %+.1e  %+.1e", md->V0[core], md->V0[n-1]);
  if(n<M){
    fprintf(fo, " | %+.1e", md->V0[n]);
  }
  else{
    fprintf(fo, " |         ");
  }
  fprintf(fo, " |  %s\n", fname);
  return;
}

