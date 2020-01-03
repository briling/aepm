#include "mol.h"
#include "matrix.h"
#include "vec3.h"

#define NELEMENTS     102

typedef enum {BASIS_TYPE_GENERAL, BASIS_TYPE_GC} basis_type;

void * bas_read(FILE * fb, basis_type * btype);
void bas_print(void * bas, basis_type btype, const char s[], FILE *f);

typedef struct{
  char name[64];
  int  lsto[NELEMENTS+1]; /* list: element -> basis functions   */
  int    *  lgto;         /* list: basis function -> primitive  */
  int    *  l;            /* angular momenta                    */
  double *  a;            /* exponents                          */
  double *  w;            /* coefficients                       */
} basis;

void basel_print_gen(int n, basis * bas, int c, const char s[], FILE *f);

typedef struct{
  char name[64];
  int  ll[NELEMENTS+1];   /* list: element -> angular momenta        */
  int * l;                /* angular momenta                         */
  int * np;               /* number of primitives                    */
  int * nc;               /* number of coefficient sets              */
  int * lp;               /* list: angular momentum -> primitives    */
  int * lc;               /* list: angular momentum -> coefficients  */
  double * a;             /* exponents                               */
  double * c;             /* coefficients                            */
} basis_gc;

int count_primitives_gc(basis_gc * bas, int k, int ll);
void bas_print_gc2gen(basis_gc * bas, const char s[], FILE * f);
void basel_print_gc2gen(int n, basis_gc * bas, int c, const char s[], FILE * f);

typedef struct{
  int list[NELEMENTS+1];
  double * a;
  double * w;
} urelconst_t;

typedef struct{
  int      l; /* angular momentum            */
  int      m; /* angular momentum projection */
  int      k; /* index of nucleus            */
  int     ng; /* number of primitives        */
  double * a; /* exponents                   */
  double * w; /* coefficients                */
} atomo;

atomo * ao_fill(mol * m, void * bas, basis_type btype, int * M);
void ao_print(atomo * ao, int M, FILE * f);

/*---------------------------------------------------------------------------*/

int * oneint_fill(int M, mol * m, void * bas, basis_type btype, atomo * ao,
                  double * S, double * H, double * Dxyz,
                  int finite_nuclei, atomo * u_rel,
                  double * boys_array, FILE * fo);
void oneint_print(int M, double * S, double * H);
void oneint_fill_gc(mol * m, basis_gc * bas, atomo * ao, int * al, double * S, double * H, double * Dxyz, int finite_nuclei, atomo * u_rel, double * boys_array, FILE * fo);
void oneint_fill_gen(mol * m, int M, atomo * ao, double * S, double * H, double * Dxyz, int finite_nuclei, atomo * u_rel, double * boys_array);
int * al_gc(mol * m, basis_gc * bas);

void finite_nuc_par(int q, double * alpha, double * omega);
atomo * u_rel_fill(mol * m, urelconst_t * urelconst);
urelconst_t * urelconst_read(FILE *fb);
void urelconst_print(urelconst_t * urel, FILE * fo);

double * boys_fill(void);
void   boys_os(double f[], unsigned int n, double T, double * boys_array);

double int_s00(int l, double a);
void int_ij(mol * m, atomo * ao1, atomo * ao2, double * S, double * H, double * Dxyz, int finite_nuclei, atomo * u_rel, double * boys_array);
double int_ij_overlap_self(atomo * ao);
double int_ij0_sum(mol * m, atomo * aoi, atomo * aoj, atomo * aoqs, double * boys_array);
double int_ij0_sum2(mol * m, int mypar, int atpar[], atomo * aoi, atomo * aoj, atomo * aoqs, double * boys_array);
double int_ij0_sum_grad(mol * m, int mypar, int atpar[], atomo * aoi, atomo * aoj, atomo * aoqs, double * boys_array);

int elnumber(mol * m);
int corepairs_simple(mol * m);

void D_fill  (int N, int M, double * C, double * D, unsigned int w);
void F_restore(int M, double * F, double * S, double * V, double * C);
void MO_proj(int m, int M, double * P, double * S, double * C);
int  C_maxoverlap(int n, int n1, int M, double * Asym, double * C, double * C0, double * S);

void init_lb20_heff(int M, double * F, double * H, mol * m, atomo * ao, int * al, basis_type btype, void * bas, double * boys_array);

int pvec_write(int M, atomo * ao, double * Va, double * Ca, double * Vb, double * Cb, const char s[]);
int pvec_read (int M, atomo * ao, double * Va, double * Ca, double * Vb, double * Cb, const char s[]);

void at0cv_prep();
void at0cv_fill(atomo * ao, mol * m);
void atcv_prep();
void atcv_fill(atomo * ao, mol * m);
void atcv_nozzle_fill_all(atomo * ao, int * atpar, double * a, double * w, mol * m);
void atcv_nozzle_fill_all1(atomo * ao, double * a, double * w, mol * m);
void atcv_nozzle_fill_all2(atomo * ao, double * a, double * w, mol * m);
void charge1_norm(double * a, double * w);
void atcv_ij_add(mol * m, int M, double * H, atomo * ao, atomo * aox, double * boys_array);
void atcv_ij_add2(mol * m, int M, double * Hs, int npar, int atpar[], atomo * ao, atomo * aox, double * boys_array);
void atcv_ij_add_grad(mol * m, int M, double * Gs, int npar, int atpar[], atomo * ao, atomo * aox, double * boys_array);
void atcv_ij_add_gc(mol * m, int M, double * Hs, int atpar[], basis_gc * bas, atomo * ao, int * al, atomo * aox, double * boys_array);
void atcv_ij_add_gc_grad(mol * m, int M, double * Hs, double * Gs, int atpar[], basis_gc * bas, atomo * ao, int * al, atomo * aox, double * boys_array);
void atcv_ij_add_gc_peratom(mol * m, int M, double * Hs, basis_gc * bas, atomo * ao, int * al, atomo * aox, double * boys_array);

