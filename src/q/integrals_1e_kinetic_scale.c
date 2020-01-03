#include "q.h"
#include "integrals.h"
#include "gc.h"
#include "nlmc.h"

double kinetic_scale(
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    double a, double b, double u,
    double ra[3], double rb[3], double ru[3],
    double rab[3]){
  /* Computes ``scaled'' kinetic energy  (up to a constant)
   *     1/2 < \phi_b | \nabla u(r) \nabla | \phi_a > =
   *  = -1/2 < \nabla \phi_b | u(r) \nabla \phi_a >,
   * where u(r) is an atom-centered gaussian with l=0.
   *
   * J. Chem. Phys. 150, 061103 (2019)
   */
  double ab    = a+b;
  double q     = ab+u;
  double q1    = 1.0/q;
  double q21   = 0.5*q1;
  double t1    = SQRTPI*sqrt(q1);
  double rabu2 = r3d2(rab,ru);
  double E     = exp(-ab*u*q1*rabu2);

  double Q[3], QA[3], QB[3];
  r3sums(Q, rab, ab, ru, u);
  r3scal(Q, q1);
  r3diff(QA, Q, ra);
  r3diff(QB, Q, rb);

  double d0 = def0(n1, n2, QA[0], QB[0], q21);
  double e0 = def0(l1, l2, QA[1], QB[1], q21);
  double f0 = def0(m1, m2, QA[2], QB[2], q21);
  double txx = e0*f0 * txyz(n1, n2, QA[0], QB[0], a, b, q21);
  double tyy = d0*f0 * txyz(l1, l2, QA[1], QB[1], a, b, q21);
  double tzz = e0*d0 * txyz(m1, m2, QA[2], QB[2], a, b, q21);
  return (txx+tyy+tzz) * E * (t1*t1*t1);
}

/*-------------------------------------------------------*/

static double kinetic_scale_prim_fill(
    int l1, int l2,
    double a, double b, double u,
    double d0[L_MAX+1+1][L_MAX+1+1],
    double e0[L_MAX+1+1][L_MAX+1+1],
    double f0[L_MAX+1+1][L_MAX+1+1],
    double tx[L_MAX+1][L_MAX+1],
    double ty[L_MAX+1][L_MAX+1],
    double tz[L_MAX+1][L_MAX+1],
    double ra[3], double rb[3], double ru[3],
    double rab[3]){

  double ab    = a+b;
  double q     = ab+u;
  double q1    = 1.0/q;
  double q21   = 0.5*q1;
  double t1    = SQRTPI*sqrt(q1);
  double rabu2 = r3d2(rab,ru);
  double E     = exp(-ab*u*q1*rabu2);

  double Q[3], QA[3], QB[3];
  r3sums(Q, rab, ab, ru, u);
  r3scal(Q, q1);
  r3diff(QA, Q, ra);
  r3diff(QB, Q, rb);

  for(int i=0; i<=l1+1; i++){
    for(int j=0; j<=l2+1; j++){
      d0[i][j] = def0(i, j, QA[0], QB[0], q21);
      e0[i][j] = def0(i, j, QA[1], QB[1], q21);
      f0[i][j] = def0(i, j, QA[2], QB[2], q21);
    }
  }
  for(int i=0; i<=l1; i++){
    for(int j=0; j<=l2; j++){
      tx[i][j] = txyz1(i, j, a, b, d0);
      ty[i][j] = txyz1(i, j, a, b, e0);
      tz[i][j] = txyz1(i, j, a, b, f0);
    }
  }

  return  E * (t1*t1*t1);
}

static double kinetic_scale_prim_sum(
    int l1, int l2, int m1, int m2,
    double ds0[L_MAX+1+1][L_MAX+1+1],
    double es0[L_MAX+1+1][L_MAX+1+1],
    double fs0[L_MAX+1+1][L_MAX+1+1],
    double tx[L_MAX+1][L_MAX+1],
    double ty[L_MAX+1][L_MAX+1],
    double tz[L_MAX+1][L_MAX+1]){

  double t = 0.0;
  int K1,K2;
  nlmc_str * g1 = pxyz(l1, m1, &K1);
  nlmc_str * g2 = pxyz(l2, m2, &K2);
  for(int k1=0; k1<K1; k1++){
    int n1 = g1[k1].n;
    int l1 = g1[k1].l;
    int m1 = g1[k1].m;
    for(int k2=0; k2<K2; k2++){
      int n2 = g2[k2].n;
      int l2 = g2[k2].l;
      int m2 = g2[k2].m;
      double d0 = ds0[n1][n2];
      double e0 = es0[l1][l2];
      double f0 = fs0[m1][m2];
      double txx = e0*f0 * tx[n1][n2];
      double tyy = d0*f0 * ty[l1][l2];
      double tzz = e0*d0 * tz[m1][m2];
      t += (txx+tyy+tzz) * g1[k1].c * g2[k2].c;
    }
  }
  return t;
}

void kinetic_scale_prim(
    int l1, int l2, double a, double b,
    preij_t pre, shd_t * pint,
    double r1[3], double r2[3], mol * m, atomo * u_rel){

  for(int k=0; k<m->n; k++){
    for(int i=0; i<u_rel[k].ng; i++){
      double * rk = m->r+k*3;
      double wu = u_rel[k].w[i];
      double au = u_rel[k].a[i];
      double d0[L_MAX+1+1][L_MAX+1+1];
      double e0[L_MAX+1+1][L_MAX+1+1];
      double f0[L_MAX+1+1][L_MAX+1+1];
      double tx[L_MAX+1][L_MAX+1];
      double ty[L_MAX+1][L_MAX+1];
      double tz[L_MAX+1][L_MAX+1];
      double Et3 = kinetic_scale_prim_fill(l1,l2,a,b,au,
          d0,e0,f0,tx,ty,tz, r1,r2,rk,pre.P);
      for(int m1=-l1; m1<=l1; m1++){
        for(int m2=-l2; m2<=l2; m2++){
          int ind1 = (l1+m1)*(2*l2+1) + (l2+m2);
          int ind  = ind1;
          pint[ind].t += 0.5 * wu * Et3 * pre.E *
            kinetic_scale_prim_sum(l1,l2,m1,m2,d0,e0,f0,tx,ty,tz);
        }
      }
    }
  }
  return;
}

