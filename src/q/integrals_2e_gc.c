#include "q.h"
#include "gc.h"
#include "nlmc.h"

#define LLL  ( (L_MAX+1) * (L_MAX+1) * (L_MAX*2+1) )

static void rnlm_ij0_g_cpy(
    rnlm_t rnlm_tmp, rnlm_t rnlm, rnlm_t rnlm_g,
    int lsum, double a, double c){

  double c2 = -0.25/(a*a) * c;
  for(int N1=0; N1<=lsum; N1++){
    for(int L1=0; L1<=lsum; L1++){
      for(int M1=0; M1<=lsum; M1++){
        rnlm_g[N1][L1][M1] += c2 * (
            rnlm_tmp[N1+2][L1][M1] +
            rnlm_tmp[N1][L1+2][M1] +
            rnlm_tmp[N1][L1][M1+2] );
        rnlm[N1][L1][M1] += c * rnlm_tmp[N1][L1][M1];
      }
    }
  }
  return;
}

static inline void filldef_allprim(int np1, int np2,
                                   int k1,  int k2,
                                   int l1,  int l2,
                                   double a1[], double a2[],
                                   mol * m, preij_t * preP,
                                   void ** Ds, void ** Es, void ** Fs){
  for(int p1=0; p1<np1; p1++){
    for(int p2=0; p2<np2; p2++){
      int p1p2 = p1*np2+p2;
      preP[p1p2] = preij(m->r, a1[p1], a2[p2], k1, k2);
      filldef_all(l1, l2, Ds[p1p2], preP[p1p2].PA[0], preP[p1p2].PB[0], preP[p1p2].p21);
      filldef_all(l1, l2, Es[p1p2], preP[p1p2].PA[1], preP[p1p2].PB[1], preP[p1p2].p21);
      filldef_all(l1, l2, Fs[p1p2], preP[p1p2].PA[2], preP[p1p2].PB[2], preP[p1p2].p21);
    }
  }
  return;
}

/*-------------------------------------------------------*/

void atcv_ij_add_gc(mol * m, int M, double * Hs,
                     int atpar[],
                     basis_gc * bas, atomo * ao, int * al,
                     atomo * aox, double * boys_array){

  size_t ij_size = ij_pint_memory1(m, bas);
  double * pint = malloc(ij_size * sizeof(double));

  size_t np2 = ij_pint_memory(m, bas);
  preij_t * preP = malloc(sizeof(preij_t)*np2);
  void ** Ds = malloc( (sizeof(void *) + sizeof(double)*LLL ) * np2*3 );
  void ** Es = Ds + np2;
  void ** Fs = Es + np2;
  for(int i=0; i<np2; i++){
    Ds[i] = Fs + np2 + LLL * (3*i+0);
    Es[i] = Fs + np2 + LLL * (3*i+1);
    Fs[i] = Fs + np2 + LLL * (3*i+2);
  }

  double * H = Hs;

  for(int k1=0; k1<m->n; k1++){
    for(int k2=k1; k2<m->n; k2++){

      int q1 = m->q[k1];
      int q2 = m->q[k2];

      for(int ll1=bas->ll[q1-1]; ll1<bas->ll[q1]; ll1++){
        for(int ll2=(k2==k1? ll1 : bas->ll[q2-1]) ; ll2<bas->ll[q2]; ll2++){
          int l1 = bas->l[ll1];
          int l2 = bas->l[ll2];
          int ls[2] = {l1, l2 };

          double * a1 = bas->a + bas->lp[ll1];
          double * a2 = bas->a + bas->lp[ll2];

          int np1 = bas->np[ll1];
          int np2 = bas->np[ll2];

          int nc1 = bas->nc[ll1];
          int nc2 = bas->nc[ll2];

          int nl1 = 2*l1+1;
          int nl2 = 2*l2+1;

          filldef_allprim(np1, np2, k1, k2, l1, l2, a1, a2, m, preP, Ds, Es, Fs);

          if (atpar ){

            for(int k3=0; k3<m->n; k3++){
              int q3 = m->q[k3];
              int ap = atpar[q3];
              if(ap<0){
                continue;
              }
              H = Hs + symsize(M)*ap;
              double * a3 = aox[k3].a;
              int np3 = aox[k3].ng;

              for(int p1=0; p1<np1; p1++){
                for(int p2=0; p2<np2; p2++){

                  rnlm_t rnlm = {0};
                  int p1p2 = p1*np2+p2;
                  for(int p3=0; p3<np3; p3++){
                    preijkl_t prem = preij0(m->r+3*k3, a3[p3], l1+l2, preP+p1p2, boys_array);
                    double c3 = aox[k3].w[p3] * prem.fact;
                    rnlm_fill_adds(l1+l2, c3, &prem, rnlm);
                  }

                  for(int m1=-l1; m1<=l1; m1++){
                    for(int m2=-l2; m2<=l2; m2++){
                      int ind1 = (l1+m1)*nl2 + (l2+m2);
                      int ind  = p1p2*(nl1*nl2) + ind1;
                      int ms[2] = {m1, m2 };
                      pint[ind] = int_ij0_prim(ls, ms, Ds[p1p2], Es[p1p2], Fs[p1p2], rnlm);
                    }
                  }

                }
              }

              for(int m1=-l1; m1<=l1; m1++){
                for(int m2=-l2; m2<=l2; m2++){
                  int s1 = al[k1*(L_MAX+1)+l1] + (l1+m1);
                  int s2 = al[k2*(L_MAX+1)+l2] + (l2+m2);
                  int ind1 = (l1+m1)*nl2 + (l2+m2);
                  for(int j1=0; j1<nc1; j1++){
                    for(int j2=(s1==s2?j1:0); j2<nc2; j2++){
                      int i1 = s1 + (2*l1+1)*j1;
                      int i2 = s2 + (2*l2+1)*j2;
                      if(i1>i2){
                        continue;
                      }
                      double tmp = 0.0;
                      for(int p1=0; p1<np1; p1++){
                        for(int p2=0; p2<np2; p2++){
                          int p1p2 = p1*np2 + p2;
                          int ind  = p1p2*(nl1*nl2) + ind1;
                          double c1 = ao[i1].w[p1];
                          double c2 = ao[i2].w[p2];
                          tmp += c1*c2 * pint[ind];
                        }
                      }
                      H[mpos(i1,i2)] += tmp;
                    }
                  }
                }
              }

            }
          }

          else{

            for(int p1=0; p1<np1; p1++){
              for(int p2=0; p2<np2; p2++){
                int p1p2 = p1*np2+p2;
                rnlm_t rnlm = {0};

                for(int k3=0; k3<m->n; k3++){
                  double * a3 = aox[k3].a;
                  int np3 = aox[k3].ng;
                  for(int p3=0; p3<np3; p3++){
                    //printf("%lf\n", a3[p3]);
                    preijkl_t prem = preij0(m->r+3*k3, a3[p3], l1+l2, preP+p1p2, boys_array);
                    double c3 = aox[k3].w[p3] * prem.fact;
                    rnlm_fill_adds(l1+l2, c3, &prem, rnlm);
                  }
                }

                for(int m1=-l1; m1<=l1; m1++){
                  for(int m2=-l2; m2<=l2; m2++){
                    int ind1 = (l1+m1)*nl2 + (l2+m2);
                    int ind  = p1p2*(nl1*nl2) + ind1;
                    int ms[2] = {m1, m2 };
                    pint[ind] = int_ij0_prim(ls, ms, Ds[p1p2], Es[p1p2], Fs[p1p2], rnlm);
                  }
                }

              }
            }

            for(int m1=-l1; m1<=l1; m1++){
              for(int m2=-l2; m2<=l2; m2++){
                int s1 = al[k1*(L_MAX+1)+l1] + (l1+m1);
                int s2 = al[k2*(L_MAX+1)+l2] + (l2+m2);
                int ind1 = (l1+m1)*nl2 + (l2+m2);
                for(int j1=0; j1<nc1; j1++){
                  for(int j2=(s1==s2?j1:0); j2<nc2; j2++){
                    int i1 = s1 + (2*l1+1)*j1;
                    int i2 = s2 + (2*l2+1)*j2;
                    if(i1>i2){
                      continue;
                    }
                    double tmp = 0.0;
                    for(int p1=0; p1<np1; p1++){
                      for(int p2=0; p2<np2; p2++){
                        int p1p2 = p1*np2 + p2;
                        int ind  = p1p2*(nl1*nl2) + ind1;
                        double c1 = ao[i1].w[p1];
                        double c2 = ao[i2].w[p2];
                        tmp += c1*c2 * pint[ind];
                      }
                    }
                    H[mpos(i1,i2)] += tmp;
                  }
                }
              }
            }

          }

        }
      }
    }
  }
  free(preP);
  free(pint);
  free(Ds);
  return;
}

void atcv_ij_add_gc_grad(mol * m, int M,
    double * H, double * Gs, int atpar[],
    basis_gc * bas, atomo * ao, int * al,
    atomo * aox, double * boys_array){

  size_t ij_size  = ij_pint_memory1(m, bas);
  double * pint   = malloc(2*ij_size * sizeof(double));
  double * pint_g = pint+ij_size;

  size_t np2 = ij_pint_memory(m, bas);
  preij_t * preP = malloc(sizeof(preij_t)*np2);
  void ** Ds = malloc( (sizeof(void *) + sizeof(double)*LLL ) * np2*3 );
  void ** Es = Ds + np2;
  void ** Fs = Es + np2;
  for(int i=0; i<np2; i++){
    Ds[i] = Fs + np2 + LLL * (3*i+0);
    Es[i] = Fs + np2 + LLL * (3*i+1);
    Fs[i] = Fs + np2 + LLL * (3*i+2);
  }

  for(int k1=0; k1<m->n; k1++){
    for(int k2=k1; k2<m->n; k2++){

      int q1 = m->q[k1];
      int q2 = m->q[k2];

      for(int ll1=bas->ll[q1-1]; ll1<bas->ll[q1]; ll1++){
        for(int ll2=(k2==k1? ll1 : bas->ll[q2-1]) ; ll2<bas->ll[q2]; ll2++){
          int l1 = bas->l[ll1];
          int l2 = bas->l[ll2];
          int ls[2] = {l1, l2 };
          int lsum = l1+l2;

          double * a1 = bas->a + bas->lp[ll1];
          double * a2 = bas->a + bas->lp[ll2];

          int np1 = bas->np[ll1];
          int np2 = bas->np[ll2];

          int nc1 = bas->nc[ll1];
          int nc2 = bas->nc[ll2];

          int nl1 = 2*l1+1;
          int nl2 = 2*l2+1;

          filldef_allprim(np1, np2, k1, k2, l1, l2, a1, a2, m, preP, Ds, Es, Fs);

          for(int k3=0; k3<m->n; k3++){
            int q3 = m->q[k3];
            int ap = atpar[q3];
            if(ap<0){
              continue;
            }
            double * G  = Gs + symsize(M)*ap;
            double * a3 = aox[k3].a;
            int np3 = aox[k3].ng;

            for(int p1=0; p1<np1; p1++){
              for(int p2=0; p2<np2; p2++){

                rnlm_t rnlm   = {};
                rnlm_t rnlm_g = {};
                int p1p2 = p1*np2+p2;
                for(int p3=0; p3<np3; p3++){
                  preijkl_t prem = preij0(m->r+3*k3, a3[p3], lsum+2, preP+p1p2, boys_array);
                  double c3 = aox[k3].w[p3] * prem.fact;
                  rnlm_t rnlm_tmp;
                  rnlm_fill(lsum+2, &prem, rnlm_tmp);
                  rnlm_ij0_g_cpy(rnlm_tmp, rnlm, rnlm_g, lsum, a3[p3], c3);
                }

                for(int m1=-l1; m1<=l1; m1++){
                  for(int m2=-l2; m2<=l2; m2++){
                    int ind1 = (l1+m1)*nl2 + (l2+m2);
                    int ind  = p1p2*(nl1*nl2) + ind1;
                    int ms[2] = {m1, m2 };
                    pint[ind]   = int_ij0_prim(ls, ms, Ds[p1p2], Es[p1p2], Fs[p1p2], rnlm);
                    pint_g[ind] = int_ij0_prim(ls, ms, Ds[p1p2], Es[p1p2], Fs[p1p2], rnlm_g);
                  }
                }

              }
            }

            for(int m1=-l1; m1<=l1; m1++){
              for(int m2=-l2; m2<=l2; m2++){
                int s1 = al[k1*(L_MAX+1)+l1] + (l1+m1);
                int s2 = al[k2*(L_MAX+1)+l2] + (l2+m2);
                int ind1 = (l1+m1)*nl2 + (l2+m2);
                for(int j1=0; j1<nc1; j1++){
                  for(int j2=(s1==s2?j1:0); j2<nc2; j2++){
                    int i1 = s1 + (2*l1+1)*j1;
                    int i2 = s2 + (2*l2+1)*j2;
                    if(i1>i2){
                      continue;
                    }
                    double tmp = 0.0;
                    double gmp = 0.0;
                    for(int p1=0; p1<np1; p1++){
                      for(int p2=0; p2<np2; p2++){
                        int p1p2 = p1*np2 + p2;
                        int ind  = p1p2*(nl1*nl2) + ind1;
                        double c1 = ao[i1].w[p1];
                        double c2 = ao[i2].w[p2];
                        tmp += c1*c2 * pint[ind];
                        gmp += c1*c2 * pint_g[ind];
                      }
                    }
                    H[mpos(i1,i2)] += tmp;
                    G[mpos(i1,i2)] += gmp;
                  }
                }
              }
            }

          }

        }
      }
    }
  }
  free(preP);
  free(pint);
  free(Ds);
  return;
}

void atcv_ij_add_gc_peratom(mol * m, int M, double * Hs,
                     basis_gc * bas, atomo * ao, int * al,
                     atomo * aox, double * boys_array){

  size_t ij_size = ij_pint_memory1(m, bas);
  double * pint = malloc(ij_size * sizeof(double));

  size_t np2 = ij_pint_memory(m, bas);
  preij_t * preP = malloc(sizeof(preij_t)*np2);
  void ** Ds = malloc( (sizeof(void *) + sizeof(double)*LLL ) * np2*3 );
  void ** Es = Ds + np2;
  void ** Fs = Es + np2;
  for(int i=0; i<np2; i++){
    Ds[i] = Fs + np2 + LLL * (3*i+0);
    Es[i] = Fs + np2 + LLL * (3*i+1);
    Fs[i] = Fs + np2 + LLL * (3*i+2);
  }

  double * H = Hs;

  for(int k1=0; k1<m->n; k1++){
    for(int k2=k1; k2<m->n; k2++){

      int q1 = m->q[k1];
      int q2 = m->q[k2];

      for(int ll1=bas->ll[q1-1]; ll1<bas->ll[q1]; ll1++){
        for(int ll2=(k2==k1? ll1 : bas->ll[q2-1]) ; ll2<bas->ll[q2]; ll2++){
          int l1 = bas->l[ll1];
          int l2 = bas->l[ll2];
          int ls[2] = {l1, l2 };

          double * a1 = bas->a + bas->lp[ll1];
          double * a2 = bas->a + bas->lp[ll2];

          int np1 = bas->np[ll1];
          int np2 = bas->np[ll2];

          int nc1 = bas->nc[ll1];
          int nc2 = bas->nc[ll2];

          int nl1 = 2*l1+1;
          int nl2 = 2*l2+1;

          filldef_allprim(np1, np2, k1, k2, l1, l2, a1, a2, m, preP, Ds, Es, Fs);

          for(int k3=0; k3<m->n; k3++){
            H = Hs + symsize(M)*k3;
            double * a3 = aox[k3].a;
            int np3 = aox[k3].ng;

            for(int p1=0; p1<np1; p1++){
              for(int p2=0; p2<np2; p2++){

                rnlm_t rnlm = {0};
                int p1p2 = p1*np2+p2;
                for(int p3=0; p3<np3; p3++){
                  preijkl_t prem = preij0(m->r+3*k3, a3[p3], l1+l2, preP+p1p2, boys_array);
                  double c3 = aox[k3].w[p3] * prem.fact;
                  rnlm_fill_adds(l1+l2, c3, &prem, rnlm);
                }

                for(int m1=-l1; m1<=l1; m1++){
                  for(int m2=-l2; m2<=l2; m2++){
                    int ind1 = (l1+m1)*nl2 + (l2+m2);
                    int ind  = p1p2*(nl1*nl2) + ind1;
                    int ms[2] = {m1, m2 };
                    pint[ind] = int_ij0_prim(ls, ms, Ds[p1p2], Es[p1p2], Fs[p1p2], rnlm);
                  }
                }

              }
            }

            for(int m1=-l1; m1<=l1; m1++){
              for(int m2=-l2; m2<=l2; m2++){
                int s1 = al[k1*(L_MAX+1)+l1] + (l1+m1);
                int s2 = al[k2*(L_MAX+1)+l2] + (l2+m2);
                int ind1 = (l1+m1)*nl2 + (l2+m2);
                for(int j1=0; j1<nc1; j1++){
                  for(int j2=(s1==s2?j1:0); j2<nc2; j2++){
                    int i1 = s1 + (2*l1+1)*j1;
                    int i2 = s2 + (2*l2+1)*j2;
                    if(i1>i2){
                      continue;
                    }
                    double tmp = 0.0;
                    for(int p1=0; p1<np1; p1++){
                      for(int p2=0; p2<np2; p2++){
                        int p1p2 = p1*np2 + p2;
                        int ind  = p1p2*(nl1*nl2) + ind1;
                        double c1 = ao[i1].w[p1];
                        double c2 = ao[i2].w[p2];
                        tmp += c1*c2 * pint[ind];
                      }
                    }
                    H[mpos(i1,i2)] += tmp;
                  }
                }
              }
            }

          }

        }
      }
    }
  }
  free(preP);
  free(pint);
  free(Ds);
  return;
}

