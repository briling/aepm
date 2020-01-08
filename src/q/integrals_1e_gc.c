#include "q.h"
#include "gc.h"
#include "nlmc.h"

static inline void ivecset(size_t n, int * r, int s){
  for(size_t i=0; i<n; i++){
    r[i] = s;
  }
  return;
}

int * al_gc(mol * m, basis_gc * bas){
  /* where functions with given k and l begin
   *
   * to find all functions with given k,l,m:
   * U = al[k*(L_MAX+1)+l] + (l+m)  // 1st
   * u = U + (2*l+1)*j, 0<=j<nc     // all
   */
  int * al = malloc(sizeof(int)*(m->n)*(L_MAX+1));
  ivecset(m->n*(L_MAX+1), al, -1);
  int M = 0;
  for(int k=0; k<m->n; k++){
    int q = m->q[k];
    for(int ll=bas->ll[q-1]; ll<bas->ll[q]; ll++){
      int l = bas->l[ll];
      int nc = bas->nc[ll];
      al[k*(L_MAX+1)+l] = M;
      M += (2*l+1)*nc;
    }
  }
#if 0
  for(int k=0; k<m->n; k++){
    int q = m->q[k];
    int lmax = bas->ll[q] - bas->ll[q-1] - 1;
    for(int l=0; l<=lmax; l++){
      printf("k=%d l=%d start=%d\n", k, l, al[k*(L_MAX+1)+l]);
    }
  }
#endif
  return al;
}

static void fill_txyz(int l1, int l2,
    double a1, double a2, preij_t pre,
    double D1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    double E1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    double F1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    double tx[L_MAX+1][L_MAX+1],
    double ty[L_MAX+1][L_MAX+1],
    double tz[L_MAX+1][L_MAX+1]
    ){

  double d0[L_MAX+1+1][L_MAX+1+1];
  double e0[L_MAX+1+1][L_MAX+1+1];
  double f0[L_MAX+1+1][L_MAX+1+1];
  for(int j=0; j<=l2+1; j++){
    for(int i=0; i<=l1+1; i++){
      if(i<=l1 && j<=l2){
        d0[i][j] = D1[i][j][0];
        e0[i][j] = E1[i][j][0];
        f0[i][j] = F1[i][j][0];
      }
      else{
        d0[i][j] = def0(i, j, pre.PA[0], pre.PB[0], pre.p21);
        e0[i][j] = def0(i, j, pre.PA[1], pre.PB[1], pre.p21);
        f0[i][j] = def0(i, j, pre.PA[2], pre.PB[2], pre.p21);
      }
    }
  }

  for(int i=0; i<=l1; i++){
    for(int j=0; j<=l2; j++){
      tx[i][j] = txyz1(i, j, a1, a2, d0);
      ty[i][j] = txyz1(i, j, a1, a2, e0);
      tz[i][j] = txyz1(i, j, a1, a2, f0);
    }
  }
  return;
}

static void std_prim(
    int n1, int l1, int m1,
    int n2, int l2, int m2,
    double sc2,
    double d[L_MAX*2+1],
    double e[L_MAX*2+1],
    double f[L_MAX*2+1],
    double tx[L_MAX+1][L_MAX+1],
    double ty[L_MAX+1][L_MAX+1],
    double tz[L_MAX+1][L_MAX+1],
    double P[3], double ans[5]){

  double overlap = d[0]*e[0]*f[0] * sc2;
  double ef0 = e[0]*f[0] * sc2;
  double df0 = d[0]*f[0] * sc2;
  double ed0 = e[0]*d[0] * sc2;

  double txx = ef0 * tx[n1][n2];
  double tyy = df0 * ty[l1][l2];
  double tzz = ed0 * tz[m1][m2];
  double kinetic = txx+tyy+tzz;

  double dipole[3];
  r3cpsc(dipole, P, overlap);
  if(n1+n2)  { dipole[0] += d[1]*ef0; }
  if(l1+l2)  { dipole[1] += e[1]*df0; }
  if(m1+m2)  { dipole[2] += f[1]*ed0; }

  ans[0] += overlap;
  ans[1] += kinetic*0.5;
  ans[2] += dipole[0];
  ans[3] += dipole[1];
  ans[4] += dipole[2];
  return;
}

void oneint_fill_gc(mol * m, basis_gc * bas, atomo * ao, int * al,
                    double * S, double * H, double * Dxyz,
                    int finite_nuclei, atomo * u_rel,
                    double * boys_array, FILE * fo){

  size_t pint_size = ij_pint_memory1(m, bas);
  shd_t * pint = malloc(sizeof(shd_t)*pint_size);

  for(int k1=0; k1<m->n; k1++){
    for(int k2=k1; k2<m->n; k2++){
      if(fo){
        progress(k2, m->n, 12, fo);
      }
      int q1 = m->q[k1];
      int q2 = m->q[k2];
      double * r1 = m->r+3*k1;
      double * r2 = m->r+3*k2;

      for(int ll1=bas->ll[q1-1]; ll1<bas->ll[q1]; ll1++){
        for(int ll2=(k2==k1? ll1 : bas->ll[q2-1]) ; ll2<bas->ll[q2]; ll2++){
          int l1 = bas->l[ll1];
          int l2 = bas->l[ll2];

          int nl1 = 2*l1+1;
          int nl2 = 2*l2+1;

          double * a1 = bas->a + bas->lp[ll1];
          double * a2 = bas->a + bas->lp[ll2];

          int np1 = bas->np[ll1];
          int np2 = bas->np[ll2];

          int nc1 = bas->nc[ll1];
          int nc2 = bas->nc[ll2];

          for(int p1=0; p1<np1; p1++){
            for(int p2=0; p2<np2; p2++){

              int ind0 = (p1*np2+p2) * (nl1*nl2);
              preij_t pre = preij(m->r, a1[p1], a2[p2], k1, k2);
              double t2 = SQRTPI*M_PI*2.0*M_SQRT2 * pre.p21*sqrt(pre.p21);

              double D1[L_MAX+1][L_MAX+1][L_MAX*2+1];
              double E1[L_MAX+1][L_MAX+1][L_MAX*2+1];
              double F1[L_MAX+1][L_MAX+1][L_MAX*2+1];
              filldef_all(l1, l2, D1, pre.PA[0], pre.PB[0], pre.p21);
              filldef_all(l1, l2, E1, pre.PA[1], pre.PB[1], pre.p21);
              filldef_all(l1, l2, F1, pre.PA[2], pre.PB[2], pre.p21);

              double tx[L_MAX+1][L_MAX+1] = {0};
              double ty[L_MAX+1][L_MAX+1] = {0};
              double tz[L_MAX+1][L_MAX+1] = {0};
              fill_txyz(l1, l2, a1[p1], a2[p2], pre, D1,E1,F1, tx,ty,tz);

              rnlm_t rnlm_nuc={0};
              if(finite_nuclei){
                rnlm_fill_nucattr_finite (l1+l2, rnlm_nuc, a1[p1]+a2[p2], pre.P, m, boys_array);
              }
              else{
                rnlm_fill_nucattr_point  (l1+l2, rnlm_nuc, a1[p1]+a2[p2], pre.P, m, boys_array);
              }

              for(int m1=-l1; m1<=l1; m1++){
                for(int m2=-l2; m2<=l2; m2++){
                  int ind = ind0 + (l1+m1)*nl2 + (l2+m2);

                  int K1,K2;
                  nlmc_str * g1 = pxyz(l1, m1, &K1);
                  nlmc_str * g2 = pxyz(l2, m2, &K2);
                  double ints[5] = {0.0};
                  double pot = 0.0;

                  for(int k1=0; k1<K1; k1++){
                    for(int k2=0; k2<K2; k2++){
                      double ww = g1[k1].c * g2[k2].c * pre.E;

                      std_prim(
                          g1[k1].n, g1[k1].l, g1[k1].m,
                          g2[k2].n, g2[k2].l, g2[k2].m,
                          ww*t2,
                          D1[ g1[k1].n ][ g2[k2].n ],
                          E1[ g1[k1].l ][ g2[k2].l ],
                          F1[ g1[k1].m ][ g2[k2].m ],
                          tx, ty, tz,
                          pre.P, ints);

                      pot += ww * add_prim_nucattr(
                          g1[k1].n+g2[k2].n,
                          g1[k1].l+g2[k2].l,
                          g1[k1].m+g2[k2].m,
                          D1[ g1[k1].n ][ g2[k2].n ],
                          E1[ g1[k1].l ][ g2[k2].l ],
                          F1[ g1[k1].m ][ g2[k2].m ],
                          rnlm_nuc);

                    }
                  }

                  if(finite_nuclei){
                    pot *= 2.0*M_PI*M_PI*SQRTPI;
                  }
                  else{
                    pot *= 4.0*M_PI*pre.p21;
                  }

                  pint[ind].v    = -pot;
                  pint[ind].s    =  ints[0];
                  pint[ind].t    =  ints[1];
                  pint[ind].d[0] =  ints[2];
                  pint[ind].d[1] =  ints[3];
                  pint[ind].d[2] =  ints[4];


                }
              }

              if(u_rel){
                kinetic_scale_prim(l1, l2, a1[p1], a2[p2], pre, pint+ind0, r1, r2, m, u_rel);
              }

            }
          }

          for(int m1=-l1; m1<=l1; m1++){
            for(int m2=-l2; m2<=l2; m2++){
              int U = al[k1*(L_MAX+1)+l1] + (l1+m1);
              int V = al[k2*(L_MAX+1)+l2] + (l2+m2);

              for(int j1=0; j1<nc1; j1++){
                for(int j2=(U==V?j1:0); j2<nc2; j2++){
                  int u = U + nl1*j1;
                  int v = V + nl2*j2;
                  if(u>v){
                    continue;
                  }
                  shd_t tmp={};
                  for(int p1=0; p1<np1; p1++){
                    for(int p2=0; p2<np2; p2++){
                      int p1p2 = p1*np2+p2;
                      int ind1 = (l1+m1)*nl2 + (l2+m2);
                      int ind  = p1p2*(nl1*nl2) + ind1;
                      double c1 = ao[u].w[p1];
                      double c2 = ao[v].w[p2];
                      tmp.s += c1*c2 * pint[ind].s;
                      tmp.t += c1*c2 * pint[ind].t;
                      tmp.v += c1*c2 * pint[ind].v;
                      r3adds(tmp.d, pint[ind].d, c1*c2);
                    }
                  }
                  int uv = mpos(u,v);
                  S[uv] = tmp.s;
                  H[uv] = tmp.v + tmp.t;
                  if(Dxyz){
                    r3cp(Dxyz+3*uv, tmp.d);
                  }
                }
              }
            }
          }

        }
      }
    }
  }
  free(pint);
  return;
}

