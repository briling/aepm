#define NATOMS NELEMENTS

enum measure_enum { MEASURE_E, MEASURE_S, MEASURE_E0, MEASURE_S0 };

static const char measure_names[][4] = {
  [MEASURE_E ] = "E",
  [MEASURE_S ] = "S",
  [MEASURE_E0] = "E0",
  [MEASURE_S0] = "S0",
};

typedef struct{
  int M, N, core;
  mol * m;
  atomo * ao;
  atomo * aox;
  int * al;
  double * S;
  double * H;
  double * X0;
  double * V0;
  double * C0;
  double * F0;
  double * P0;
  int * mpl;
} mol_data;

typedef struct{
  int nmol;
  mol_data * md;
  int finite_nuclei;
  urelconst_t * urelconst;
  basis_type btype;
  void * bas;
  double * boys_array;

  int measure;
  int npar;
  double * g;
  double * e;
  int * atpar;

} thread_arg;

typedef struct{
  double h;
  double D;
  double H;
  double G;
  int    K;
} o1dstr;

typedef struct{
  double MG;
  int    K;
  o1dstr o1dpars;
} ostr;

int pars_read(double a0[NATOMS], int bp[NATOMS], int fixflag[NATOMS], FILE * fl);
int pars_default(int mode, double a0[NATOMS], int bp[NATOMS]);
int atpar_fill(int nparams, int nmol, double fa[NATOMS], double a[NATOMS], double a0[NATOMS], int bp[NATOMS], int fixflag[NATOMS], int * nfpar_r, int fatpar[NATOMS+1], int  atpar[NATOMS+1], mol_data * md, FILE * fo);
void conv_par_print(int atpar[NATOMS+1], int fatpar[NATOMS+1], double * a, double * fa, FILE * fo);

char ** flist_read(int * nlines, FILE * fl);

int mol_init(mol_data * md, basis_type btype, void * bas, char * fname);
int remove_oddmols(int nmol, mol_data * md, int * atpar, char ** fnames);
int * molparlist(int nmol, int npar, mol_data * md, int * atpar);
void mol_1el(mol_data * md, int finite_nuclei, urelconst_t * urelconst, basis_type btype, void * bas, double * boys_array);
void mol_readvec(mol_data * md, char * fname, FILE * fo);

void mols_clusterize(int npar, int nmol, mol_data * md, char ** fnames, FILE * fo);

int mol2proc(int nproc, int nmol, pthread_t ** threads, thread_arg ** args, mol_data * md);
void * mol_1el_ll(void * arg);
void * mol_atcv_ll(void * arg);
void * calc_gradient_mol_ll(void * arg);

int measure_ortho(int N, int M, int core, double * S, double * C0, double * C);
double measure_minoverlap(int N, int M, int core, double * S, double * C0, double * C);
double measure_E(int N, int M, int core, double * F0, double * V0, double * C);
double measure_E_grad( int npar, int N, int M, int core, int * mpl, double * F0, double * V0, double * C,  double * V, double * Gs, double * g1);
double measure_S(int N, int M, int core, double * P0, double * C);
double measure_S_grad( int npar, int N, int M, int core, int * mpl, double * P0, double * C,  double * V, double * Gs, double * g1);
double measure_S0(int N, int M, int core, double * P0, double * C);
double measure_S0_grad( int npar, int N, int M, int core, int * mpl, double * P0, double * C,  double * V, double * Gs, double * g1);
double measure_E0(int N, int M, int core, double * F0, double * V0, double * C);
double measure_E0_grad( int npar, int N, int M, int core, int * mpl, double * F0, double * V0, double * C,  double * V, double * Gs, double * g1);

double calc_gradient_mol(int measure, int npar, double * g, mol_data * md, int atpar[], basis_type btype, basis_gc * bas, double * boys_array);
void test_grad_a(int nmol, int npar, mol_data * md, int atpar[NATOMS+1], double * a, double * w, basis_type btype, basis_gc * bas, double * boys_array);
void opt_grad_conj(int npar, double a[], double weigth[], thread_arg * ta, ostr pars, FILE * f);
double calc_gradient(int npar, double * g, thread_arg * ta);

void mols_check(int nmol, mol_data * md, basis_type btype, basis_gc * bas, double * boys_array, int check, char ** fnames, FILE * fo);
void Hmod_only_solve(mol_data * md, double * V, double * C, double * F);
void sol_save(mol_data * md, const char fname[], basis_type btype, void * bas, double * boys_array);

static inline double calc_only_measure(int measure, double * C, mol_data * md){
  switch(measure){
    case MEASURE_E:
      return measure_E(md->N, md->M, md->core, md->F0, md->V0, C);
    case MEASURE_S:
      return measure_S(md->N, md->M, md->core, md->P0, C);
    case MEASURE_E0:
      return measure_E0(md->N, md->M, md->core, md->F0, md->V0, C);
    case MEASURE_S0:
      return measure_S0(md->N, md->M, md->core, md->P0, C);
    default:
      GOTOHELL;
  }
}

static inline double calc_grad_measure(int measure, int npar,
    double * g, double * C, double * V, double * Gs, mol_data * md){
  switch(measure){
    case MEASURE_E:
      return measure_E_grad(npar, md->N, md->M, md->core, md->mpl, md->F0, md->V0, C, V, Gs, g);
    case MEASURE_S:
      return measure_S_grad(npar, md->N, md->M, md->core, md->mpl, md->P0, C, V, Gs, g);
    case MEASURE_E0:
      return measure_E0_grad(npar, md->N, md->M, md->core, md->mpl, md->F0, md->V0, C, V, Gs, g);
    case MEASURE_S0:
      return measure_S0_grad(npar, md->N, md->M, md->core, md->mpl, md->P0, C, V, Gs, g);
    default:
      GOTOHELL;
  }
}

static inline void nozzle(mol_data * md, double * F,
    basis_type btype, basis_gc * bas, double * boys_array){
  if(btype == BASIS_TYPE_GC){
    atcv_ij_add_gc(md->m, md->M, F, NULL, bas, md->ao, md->al, md->aox, boys_array);
  }
  else{
    atcv_ij_add(md->m, md->M, F, md->ao, md->aox, boys_array);
  }
  return;
}

static inline void nozzle_gradient(
    int npar, mol_data * md, double * F, double * Gs, int atpar[],
    basis_type btype, basis_gc * bas, double * boys_array){
  if(btype == BASIS_TYPE_GC){
    atcv_ij_add_gc_grad(md->m, md->M, F, Gs, atpar, bas, md->ao, md->al, md->aox, boys_array);
  }
  else{
    atcv_ij_add(md->m, md->M, F, md->ao, md->aox, boys_array);
    atcv_ij_add_grad(md->m, md->M, Gs, npar, atpar, md->ao, md->aox, boys_array);
  }
  return;
}

