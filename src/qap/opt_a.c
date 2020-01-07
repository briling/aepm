#include "q.h"
#include "qap.h"

static double opt1d(int npar, double * p,
    double * a, double * weigth, thread_arg * ta,
    double * energy, double * gr, double * dir,
    o1dstr pars, FILE * f){
#define ADDS(N,L,D)         vecadds(N, p, L, D);\
                            for(int i=0; i<npar; i++){\
                              a[i] = exp(p[i]);\
                              charge1_norm(a+i, weigth+i);\
                            }
#define GETN                npar;
#define WALK(ENERGY,GRAD)   *(ENERGY) = calc_gradient(npar, GRAD, ta);\
                            for(int i=0; i<npar; i++){\
                              GRAD[i] *= a[i];\
                            }
#include "opt1d.h"
}
#define  GETMG(N,G)            vecabsmax(N,G);
#define  PRINTEGK(F,E,G,K)     fprintf(F, "k=%3d   E=%20.15lf   mg=%8.6e\n\n", K, E, G);
#define  PRINTENT(F)           fprintf(F, "a:  "); vecprint(npar, a, "", F);
#define  OPT1D(ENERGY, GRAD, DIR, PARS, F)  opt1d(npar, p, a, weigth, ta, ENERGY, GRAD, DIR, PARS, F)

void opt_grad_conj(int npar, double a[], double weigth[],
                          thread_arg * ta, ostr pars, FILE * f){
  double p[NATOMS];
  for(int i=0; i<npar; i++){
    p[i] = log(a[i]);
  }
#include "grad_conj.h"
  return;
}

