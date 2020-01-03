#include "q.h"
#include "task_q.h"
#include "qinit.h"
#include <stddef.h>
#include <pthread.h>
#include "mytime.h"

#define MEASURE MEASURE_S
#define NPROC 1
int nproc = NPROC;
pthread_t * threads;

static void add_fixed(int nmol, int nproc,
    mol_data * md, int fatpar[], double * fa, double * fw,
    thread_arg * args){
  for(int i=0; i<nmol; i++){
    atcv_nozzle_fill_all(md[i].aox, fatpar, fa, fw, md[i].m);
  }
  for(int i=0; i<nproc; i++){
    pthread_create(threads+i, NULL, mol_atcv_ll, args+i);
  }
  for(int i=0; i<nproc; i++){
    pthread_join(threads[i], NULL);
  }
  return;
}

static void read_solution(int nproc, char ** fnames, thread_arg * args, FILE * fo){
  for(int i=0; i<nproc; i++){
    ptrdiff_t j0 = args[i].md-args[0].md;
    for(int j=0; j<args[i].nmol; j++){
      fprintf(fo, "(%d:%d) ", i, j);
      mol_readvec(args[i].md+j, fnames[j0+j], fo);
    }
  }
  fprintf(fo, "\n");
  return;
}

static void add_atcv(int nmol, int nproc, int HFS,
                     mol_data * md, thread_arg * args){

  for(int i=0; i<nmol; i++){
    md[i].aox = malloc(sizeof(atomo)*md[i].m->n);
    if(!HFS){
      atcv_fill(md[i].aox, md[i].m);
    }
    else{
      at0cv_fill(md[i].aox, md[i].m);
    }
  }
  for(int i=0; i<nproc; i++){
    pthread_create(threads+i, NULL, mol_atcv_ll, args+i);
  }
  for(int i=0; i<nproc; i++){
    pthread_join(threads[i], NULL);
  }
  return;
}

int main(int argc, char * argv[]){

  if(argc<3){
    printf("  usage:\n\
        %s basis.in list.in [out] [arguments]\n", argv[0]);
    printf("  command-line arguments:\n\
        np:%%d       -- number of threads (if compiled)\n\
        f:%%s        -- measure\n\
        save:%%d     -- save guesses after optimization\n\
        check:%%d    -- check the measure on molecules (1 -- atcv+cap, 2 -- at0cv)\n\
        clust:%%d    -- group molecules wrt parameters\n\
        o:%%d,%%lf    -- optimization parameters (K,MG)\n\
        o1:%%d,%%lf,%%lf,%%lf,%%lf -- 1d optimization parameters (K,G,h,D,H)\n\
        \n");
    return 1;
  }

  FILE * fo = stdout;
  int check   = 0;
  int clust   = 0;
  int save    = 0;
  int measure = MEASURE;
  char measure_s[256] = {};
  ostr op = {
    .MG = 1e-5,  .K = 128,
    .o1dpars = {
      .G = 1e-6, .K = 32,
      .h = 3.0,  .D = 0.125, .H = 1e-3
    }
  };
  for(int i=3; i<argc; i++){
    if( sscanf (argv[i], "np:%d",      &nproc) ) { continue; }
    if( sscanf (argv[i], "f:%s",    measure_s) ) { continue; }
    if( sscanf (argv[i], "save:%d",     &save) ) { continue; }
    if( sscanf (argv[i], "check:%d",   &check) ) { continue; }
    if( sscanf (argv[i], "clust:%d",   &clust) ) { continue; }
    if( sscanf (argv[i], "o:%d,%lf", &op.K, &op.MG) ) { continue; }
    if( sscanf (argv[i], "o1:%d,%lf,%lf,%lf,%lf", &op.o1dpars.K, &op.o1dpars.G, &op.o1dpars.h, &op.o1dpars.D, &op.o1dpars.H ) ) { continue; }
    if(! (fo = fopen(argv[i], "w"))){
      fo = stdout;
    }
  }
  str_toupper(measure_s);
  for(int i=0; i<sizeof(measure_names)/sizeof(measure_names[0]); i++){
    if(!strcmp(measure_s, measure_names[i])){
      measure = i;
      break;
    }
  }

  double * boys_array = boys_fill();
  atcv_prep();
  at0cv_prep();

  FILE * fb;
  fb = fopen(argv[1], "r");
  if(!fb){
    free(boys_array);
    fprintf(stderr, "\tbasis?\n");
    return 1;
  }

  urelconst_t * urelconst = urelconst_read(fb);
  int finite_nuclei = !!urelconst;
#if 0
  urelconst_print(urelconst, fo);
#endif

  basis_type btype;
  void * bas = bas_read(fb, &btype);
  fclose(fb);
  if(!bas){
    fprintf(stderr, "\tbasis!\n");
    free(urelconst);
    free(boys_array);
    fclose(fo);
    return 1;
  }
#if 0
  bas_print(bas, btype, ">", fo);
#endif

  FILE * fl;
  fl = fopen(argv[2], "r");
  if(!fl){
    free(urelconst);
    free(bas);
    free(boys_array);
    fprintf(stderr, "\tlist?\n");
    return 1;
  }

  double a0[NATOMS];
  int bp[NATOMS], fixflag[NATOMS]={0};
  int nparams = pars_read(a0, bp, fixflag, fl);
  if(!nparams){
    nparams = pars_default(1, a0, bp);
  }

  int nmol;
  char ** fnames = flist_read(&nmol, fl);
  char ** fnames0  = fnames;
  fclose(fl);

//------------------------------------------------------

  mol_data * md = malloc(sizeof(mol_data)*nmol);

  for(int i=0; i<nmol; i++){
    if(mol_init(md+i, btype, bas, fnames[i])){
      fprintf(fo, "(%d) %s\n", i, fnames[i]);
      GOTOHELL;
    }
  }

  double  a[NATOMS],  w[NATOMS];
  double fa[NATOMS], fw[NATOMS];
  int  atpar[NATOMS+1];
  int fatpar[NATOMS+1];
  int nfpar, npar = atpar_fill(nparams, nmol, fa, a, a0, bp, fixflag, &nfpar, fatpar, atpar, md, fo);

  for(int i=0; i<npar; i++){
    charge1_norm(a+i, w+i);
  }
  for(int i=0; i<nfpar; i++){
    charge1_norm(fa+i, fw+i);
  }

  nmol = remove_oddmols(nmol, md, atpar, fnames);
  if(!nmol || !npar){
    GOTOHELL;
  }

  int * mpl = molparlist(nmol, npar, md, atpar);

  for(int i=0; i<nmol; i++){
    fprintf(fo, "(%d) %s\n", i, fnames[i]);
    fprintf(fo, "N = %d, M = %d, core = %d\n", md[i].N, md[i].M, md[i].core);
  }
  fprintf(fo, "\n");
  fflush(fo);

  if(clust){
    mols_clusterize(npar, nmol, md, fnames, fo);
    goto hell1;
  }

  thread_arg * args;
  nproc = mol2proc(nproc, nmol, &threads, &args, md);
  for(int i=0; i<nproc; i++){
    args[i].finite_nuclei = finite_nuclei;
    args[i].urelconst = urelconst;
    args[i].btype = btype;
    args[i].bas   = bas;
    args[i].boys_array = boys_array;
    args[i].npar  = npar;
    args[i].atpar = atpar;
    args[i].measure = measure;
  }

  fprintf(fo, "\n<<< one-electron integrals\n");
  double time_sec, rtime_sec;
  time_sec = myutime();
  rtime_sec = myrealtime();
  for(int i=0; i<nproc; i++){
    pthread_create(threads+i, NULL, mol_1el_ll, args+i);
  }
  for(int i=0; i<nproc; i++){
    pthread_join(threads[i], NULL);
  }
  printf("time (cpu)    : %.2f s\n", myutime()-time_sec);
  printf("time (real)   : %.2f s\n", myrealtime()-rtime_sec);

  fprintf(fo, "\n<<< add atcv\n");
  time_sec = myutime();
  rtime_sec = myrealtime();
  add_atcv(nmol, nproc, check==2, md, args);
  printf("time (cpu)    : %.2f s\n", myutime()-time_sec);
  printf("time (real)   : %.2f s\n", myrealtime()-rtime_sec);

  if(nfpar && check!=2){
    fprintf(fo, "\n<<< add fixed\n");
    add_fixed(nmol, nproc, md, fatpar, fa, fw, args);
  }

  fprintf(fo, "\n<<< read solution\n");
  read_solution(nproc, fnames, args, fo);

  for(int i=0; i<nmol; i++){
    atcv_nozzle_fill_all(md[i].aox, atpar, a, w, md[i].m);
  }

  if(check){
    fprintf(fo, "check=%d\n", check);
    mols_check(nmol, md, btype, bas, boys_array, check, fnames, fo);
    goto hell;
  }

#if 0
  test_grad_a(nmol, npar, md, atpar, a, w, btype, bas, boys_array);
  goto hell;
#endif

  fprintf(fo, "o:%d,%e\n", op.K, op.MG);
  fprintf(fo, "o1:%d,%e,%e,%e,%e\n", op.o1dpars.K, op.o1dpars.G, op.o1dpars.h, op.o1dpars.D, op.o1dpars.H );
  fprintf(fo, "\n");
  opt_grad_conj (npar, a, w, args, op, fo);

  if(save){
    fprintf(fo, "\n<<< save guess\n");
    for(int i=0; i<npar; i++){
      charge1_norm(a+i, w+i);
    }
    for(int i=0; i<nmol; i++){
      sol_save(md, fnames[i], btype, bas, boys_array);
    }
  }
  fprintf(fo, "\nf:%s\n", measure_names[measure]);
  fprintf(fo, "rel=%d, fn=%d\n", !!urelconst, finite_nuclei);
  conv_par_print(atpar, fatpar, a, fa, fo);

hell:
  free(threads);
  free(args);
hell1:
  for(int i=0; i<nmol; i++){
    free(md[i].m);
    free(md[i].al);
    free(md[i].ao);
    free(md[i].H);
    free(md[i].aox);
    free(md[i].V0);
  }
  free(md);
  free(mpl);
  free(fnames0);
  free(boys_array);
  free(bas);
  free(urelconst);
  fclose(fo);
  return 0;
}

