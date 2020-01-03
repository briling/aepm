#include "q.h"
#include "integrals.h"

void add_ij(int n, int l, int m, int n1, int l1, int m1,
            double a, double b, double Rij2, double w,
            double ri[3], double rj[3], mol * ml,
            int finite, atomo * u_rel,
            double ans[6], double * boys_array){

  double p   = a+b;
  double p1  = 1.0/p;
  double p21 = 0.5*p1;
  double t1  = SQRTPI*sqrt(p1);
  double t2  = t1*t1*t1;
  double Eij = exp(-a*b*p1*Rij2);
  double sc  = w * Eij;
  double sc2 = sc*t2;

  double P[3], PA[3], PB[3];
  gprod(a, b, p1, ri, rj, P, PA, PB);

  double d[L_MAX*2+1];
  double e[L_MAX*2+1];
  double f[L_MAX*2+1];
  filldef(n,l,m, n1,l1,m1, d,e,f, PA,PB, p21);

  double overlap = d[0]*e[0]*f[0] * sc2;
  double ef0 = e[0]*f[0] * sc2;
  double df0 = d[0]*f[0] * sc2;
  double ed0 = e[0]*d[0] * sc2;

  double txx = ef0 * txyz(n, n1, PA[0], PB[0], a, b, p21);
  double tyy = df0 * txyz(l, l1, PA[1], PB[1], a, b, p21);
  double tzz = ed0 * txyz(m, m1, PA[2], PB[2], a, b, p21);
  double kinetic = txx+tyy+tzz;

  if(u_rel){
    double rel = 0.0;
    for(int k=0; k<ml->n; k++){
      double * rk = ml->r+k*3;
      for(int i=0; i<u_rel[k].ng; i++){
        double wu = u_rel[k].w[i];
        double au = u_rel[k].a[i];
        rel += wu * kinetic_scale(n,l,m, n1,l1,m1, a,b,au, ri,rj,rk,P);
      }
    }
    kinetic += rel * sc;
  }

  double dipole[3];
  r3cpsc(dipole, P, overlap);
  if(n+n1)  { dipole[0] += d[1]*ef0; }
  if(l+l1)  { dipole[1] += e[1]*df0; }
  if(m+m1)  { dipole[2] += f[1]*ed0; }

  double nucattr;
  if(finite){
    nucattr  = nucattr_finite(n,l,m,n1,l1,m1, d, e, f, P, p, ml, boys_array);
  }
  else{
    nucattr  = nucattr_point(n+n1, l+l1, m+m1, d, e, f, P, p, ml, boys_array);
  }

  ans[0] += overlap;
  ans[1] += kinetic*0.5;
  ans[2] += nucattr*sc ;
  ans[3] += dipole[0];
  ans[4] += dipole[1];
  ans[5] += dipole[2];
  return;
}

double add_ij_overlap_self(
    int n, int l, int m, int n1, int l1, int m1,
    double p21){

  double d[L_MAX*2+1];
  double e[L_MAX*2+1];
  double f[L_MAX*2+1];

  filldef_self(n,l,m, n1,l1,m1, d,e,f, p21);

  return d[0]*e[0]*f[0];
}

double add_s00(int n, int l, int m, int n1, int l1, int m1, double a, double w){
  double p   = 2.0*a;
  double p1  = 1.0/p;
  double p21 = 0.5*p1;
  double t   = SQRTPI*sqrt(p1);
  double d0  = def0(n, n1, 0.0, 0.0, p21);
  double e0  = def0(l, l1, 0.0, 0.0, p21);
  double f0  = def0(m, m1, 0.0, 0.0, p21);
  return w * t*t*t * d0*e0*f0;
}

