#include "q.h"
#include "boys.h"
#include "integrals.h"
#include "nlmc.h"
#include "gc.h"

double int_s00(int l, double a){
  int K;
  nlmc_str * g = pxyz(l, 0, &K);
  double s = 0.0;
  for(int k1=0; k1<K; k1++){
    for(int k2=0; k2<K; k2++){
      double w = g[k1].c * g[k2].c;
      s += add_s00(
          g[k1].n, g[k1].l, g[k1].m,
          g[k2].n, g[k2].l, g[k2].m,
          a, w);
    }
  }
  return s;
}

void int_ij(mol * m, atomo * ao1, atomo * ao2,
    double * S, double * H, double * Dxyz,
    int finite_nuclei, atomo * u_rel,
    double * boys_array){

  double * r1 = m->r+3*ao1->k;
  double * r2 = m->r+3*ao2->k;
  double R2   = r3d2(r1, r2);
  double ints[6] = {0.0};

  int K1,K2;
  nlmc_str * g1 = pxyz(ao1->l, ao1->m, &K1);
  nlmc_str * g2 = pxyz(ao2->l, ao2->m, &K2);

  for(int i1=0; i1<ao1->ng ; i1++){
    for(int i2=0; i2<ao2->ng ; i2++){
      double a1 = ao1->a[i1];
      double a2 = ao2->a[i2];
      double w  = ao1->w[i1]*ao2->w[i2];
      for(int k1=0; k1<K1; k1++){
        for(int k2=0; k2<K2; k2++){
          double ww = w * g1[k1].c * g2[k2].c;
          add_ij(
              g1[k1].n, g1[k1].l, g1[k1].m,
              g2[k2].n, g2[k2].l, g2[k2].m,
              a1,a2, R2, ww, r1, r2, m,
              finite_nuclei, u_rel,
              ints, boys_array);
        }
      }
    }
  }
  *S = ints[0];
  *H = ints[1] - ints[2];
  if(Dxyz){
    Dxyz[0] = ints[3];
    Dxyz[1] = ints[4];
    Dxyz[2] = ints[5];
  }
  return;
}

double int_ij_overlap_self(atomo * ao){

  int K;
  nlmc_str * g1 = pxyz(ao->l, ao->m, &K);

  double S = 0.0;

  for(int i1=0; i1<ao->ng ; i1++){
    for(int i2=0; i2<ao->ng ; i2++){
      double a1 = ao->a[i1];
      double a2 = ao->a[i2];
      double w  = ao->w[i1]*ao->w[i2];

      double p1  = 1.0/(a1+a2);
      double t1  = SQRTPI*sqrt(p1);
      double t2  = t1*t1*t1;
      double sc = w * t2;

      double s = 0.0;
      for(int k1=0; k1<K; k1++){
        for(int k2=0; k2<K; k2++){
          s += g1[k1].c*g1[k2].c*add_ij_overlap_self(
              g1[k1].n, g1[k1].l, g1[k1].m,
              g1[k2].n, g1[k2].l, g1[k2].m,
              0.5*p1);
        }
      }
      S += sc * s;
    }
  }
  return S;
}

/*-------------------------------------------------------------------------------------------------*/

double int_ij0_sum(mol * m, atomo * aoi, atomo * aoj, atomo * aoqs, double * boys_array){
  /* \sum_k \int aoi(1)*aoj(1)*aoq_k(2)/r_{12}
   * l(aoq) = 0
   */
  double * ri = m->r+3*(aoi->k);
  double * rj = m->r+3*(aoj->k);
  int K1,K2;
  nlmc_str * g1 = pxyz(aoi->l, aoi->m, &K1);
  nlmc_str * g2 = pxyz(aoj->l, aoj->m, &K2);
  double Rij2 = r3d2(ri, rj);
  int lsum = aoi->l + aoj->l;

  double sum = 0.0;

  for(int i=0; i<aoi->ng ; i++){
    double ai = aoi->a[i];
    double w1 = aoi->w[i];
    for(int j=0; j<aoj->ng ; j++){
      double aj  = aoj->a[j];
      double ap  = ai + aj;
      double p1  = 1.0/ap;
      double p21 = 0.5*p1;
      double P[3], PA[3], PB[3];
      gprod(ai, aj, p1, ri, rj, P, PA, PB);
      double E1 = exp(-ai*aj*p1*Rij2);
      double w12 = w1*aoj->w[j];

      for(int k1=0; k1<K1; k1++){
        for(int k2=0; k2<K2; k2++){
          double d1[L_MAX*2+1];
          double e1[L_MAX*2+1];
          double f1[L_MAX*2+1];
          filldef(
              g1[k1].n,g1[k1].l,g1[k1].m,
              g2[k2].n,g2[k2].l,g2[k2].m,
              d1,e1,f1, PA,PB, p21);

          double d = 0.0;
          for(int k=0; k<m->n; k++){
            double * Q  = m->r+3*k;
            atomo * aoq = aoqs+k;
            for(int q=0; q<aoq->ng ; q++){
              double aq = aoq->a[q];
              double wq = aoq->w[q];
              double PQ[3];
              r3diff(PQ, P, Q);
              double r2 = r3dot(PQ,PQ);
              double t1 = 1.0/(ap+aq);
              double t2 = ap*aq;
              double alpha = t1*t2;
              double T = alpha * r2;
              double F[BOS_NMAX+1];
              boys_os(F, lsum, T, boys_array);
              double lambda = 2.0*M_PI*M_PI*SQRTPI * sqrt(t1) / t2;
              d += lambda * wq *
                g1[k1].c * g2[k2].c * add_ijs(
                    g1[k1].n,g1[k1].l,g1[k1].m,
                    g2[k2].n,g2[k2].l,g2[k2].m,
                    d1, e1, f1, PQ, alpha, F);
            }
          }
          sum += d * E1 * w12;
        }
      }
    }
  }
  return sum;
}

double int_ij0_sum_grad(mol * m, int mypar, int atpar[],
    atomo * aoi, atomo * aoj, atomo * aoqs, double * boys_array){
  /* gradient of \sum_k \int aoi(1)*aoj(1)*aoq_k(2)/r_{12}
   * if k uses parameter mypar, i.e. atpar[m->q[k]] == mypar
   */
  double * ri = m->r+3*(aoi->k);
  double * rj = m->r+3*(aoj->k);
  int K1,K2;
  nlmc_str * g1 = pxyz(aoi->l, aoi->m, &K1);
  nlmc_str * g2 = pxyz(aoj->l, aoj->m, &K2);
  double Rij2 = r3d2(ri, rj);
  int lsum = aoi->l + aoj->l;

  double sum = 0.0;

  for(int i=0; i<aoi->ng ; i++){
    double ai = aoi->a[i];
    double w1 = aoi->w[i];
    for(int j=0; j<aoj->ng ; j++){
      double aj  = aoj->a[j];
      double ap  = ai + aj;
      double p1  = 1.0/ap;
      double p21 = 0.5*p1;
      double P[3], PA[3], PB[3];
      gprod(ai, aj, p1, ri, rj, P, PA, PB);
      double E1 = exp(-ai*aj*p1*Rij2);
      double w12 = w1*aoj->w[j];

      for(int k1=0; k1<K1; k1++){
        for(int k2=0; k2<K2; k2++){
          double d1[L_MAX*2+1];
          double e1[L_MAX*2+1];
          double f1[L_MAX*2+1];
          filldef(
              g1[k1].n,g1[k1].l,g1[k1].m,
              g2[k2].n,g2[k2].l,g2[k2].m,
              d1,e1,f1, PA,PB, p21);

          double d = 0.0;
          for(int k=0; k<m->n; k++){
            if(atpar[m->q[k]]!=mypar){
              continue;
            }
            double * Q  = m->r+3*k;
            atomo * aoq = aoqs+k;
            for(int q=0; q<aoq->ng ; q++){
              double aq = aoq->a[q];
              double wq = aoq->w[q];
              double PQ[3];
              r3diff(PQ, P, Q);
              double r2 = r3dot(PQ,PQ);
              double t1 = 1.0/(ap+aq);
              double t2 = ap*aq;
              double alpha = t1*t2;
              double T = alpha * r2;
              double F[BOS_NMAX+1];
              boys_os(F, lsum+2, T, boys_array);
              double lambda = 2.0*M_PI*M_PI*SQRTPI * sqrt(t1) / t2;
              d += lambda * wq *
                g1[k1].c * g2[k2].c * add_ijs_grad(
                    g1[k1].n,g1[k1].l,g1[k1].m,
                    g2[k2].n,g2[k2].l,g2[k2].m,
                    d1, e1, f1, PQ, alpha, F, aq);
            }
          }
          sum += d * E1 * w12;
        }
      }
    }
  }
  return sum;
}

double int_ij0_prim(int l[2], int m[2],
    double D1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    double E1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    double F1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    rnlm_t rnlm){

  int K1,K2;
  nlmc_str * g1 = pxyz(l[0], m[0], &K1);
  nlmc_str * g2 = pxyz(l[1], m[1], &K2);

  double sum = 0.0;

  for(int k1=0; k1<K1; k1++){
    for(int k2=0; k2<K2; k2++){
      double * d1 = D1[ g1[k1].n ][ g2[k2].n ];
      double * e1 = E1[ g1[k1].l ][ g2[k2].l ];
      double * f1 = F1[ g1[k1].m ][ g2[k2].m ];
      sum += g1[k1].c*g2[k2].c * add_ijs1(
          rnlm,
          g1[k1].n,g1[k1].l,g1[k1].m,
          g2[k2].n,g2[k2].l,g2[k2].m,
          d1,e1,f1);
    }
  }
  return sum;
}

/*-------------------------------------------------------------------------------------------------*/

preij_t preij(double * r, double a1, double a2, int k1, int k2){
  preij_t ret;
  double * ri = r+3*k1;
  double * rj = r+3*k2;
  double Rij2 = r3d2(ri, rj);
  ret.ap = a1 + a2;
  double p1 = 1.0/ret.ap;
  gprod(a1, a2, p1, ri, rj, ret.P, ret.PA, ret.PB);
  ret.E = exp(-a1*a2*p1*Rij2);
  ret.p21 = 0.5*p1;
  return ret;
}

preijkl_t preij0(
    double * Q, double aq, int lsum,
    preij_t * pP, double * boys_array){
  preijkl_t ret;
  double ap = pP->ap;
  r3diff(ret.PQ, pP->P, Q);
  double r2 = r3dot(ret.PQ,ret.PQ);
  double t1 = 1.0/(ap+aq);
  double t2 = ap*aq;
  ret.alpha = t1*t2;
  double T = ret.alpha * r2;
  boys_os(ret.F, lsum, T, boys_array);
  double lambda = 2.0*M_PI*M_PI*SQRTPI * sqrt(t1) / t2;
  ret.fact = lambda * pP->E;
  return ret;
}

