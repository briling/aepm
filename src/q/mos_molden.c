#include "q.h"

#define FORMAT5 "% .8lf"

static void abmolden_print(mol * m, void * bas, basis_type btype, const char ps[], FILE * fo){
  const char aname[][3]={
    #include "elements.h"
  };
  fprintf(fo, "%s[Molden Format]\n", ps);
  fprintf(fo, "%s[Atoms] AU\n", ps);
  for(int i=0; i<m->n; i++){
    fprintf(fo, "%s  %s %3d %3d % 13.8lf % 13.8lf % 13.8lf\n",
        ps, aname[m->q[i]], i+1, m->q[i], m->r[3*i], m->r[3*i+1], m->r[3*i+2]);
  }
  fprintf(fo, "%s[5D]\n", ps);
  fprintf(fo, "%s[GTO]\n", ps);
  for(int i=0; i<m->n; i++){
    fprintf(fo, "%s%4d\n", ps, i+1);
    if(btype == BASIS_TYPE_GENERAL){
      basel_print_gen(m->q[i], bas, 1, ps, fo);
    }
    else{
      basel_print_gc2gen(m->q[i], bas, 1, ps, fo);
    }
    fprintf(fo, "%s\n", ps);
  }
  fprintf(fo, "%s[MO]\n", ps);
  return;
}

static inline int M1(int n){
  return (n%2) ? (n+1)/2 : -(n+1)/2;
}

static void momolden_print(int N, int M, double * C, double * V, atomo * ao, int oc, const char spin[], const char ps[], FILE * fo){
  for(int i=0; i<M; i++){
    fprintf(fo, "%sSym=C1\n", ps);
    fprintf(fo, "%sEne=% 14.7lf\n", ps, V[i]);
    fprintf(fo, "%sSpin=%s\n", ps, spin);
    fprintf(fo, "%sOccup=%d\n", ps, ((i < N) ? oc:0) );
    for(int j=0; j<M; j++){
      int l  = ao[j].l;
      int m  = ao[j].m;
      int j1 = j;
      /* see pvec_io.c
       *  P: x, y, z => (-1,0,1)->(0,1,-1)->(1,-1,0)
       * 5D: D 0, D+1, D-1, D+2, D-2
       * 7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
       * 9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
       */
      if(l==1){
        j1 += M1(M1(m+l)+l)-m;
      }
      else{
        j1 += M1(m+l)-m;
      }
      fprintf(fo, "%s %3d  "FORMAT5"\n", ps, j+1, C[M*i+j1]);
    }
  }
  fprintf(fo, "\n");
  fflush (fo);
  return;
}

void molden_print(int N, int M, int occ,
    mol * m, atomo * ao, void * bas, basis_type btype,
    double * C, double * V, const char spin[], const char ps[], FILE * fo){
  abmolden_print(m, bas, btype, ps, fo);
  momolden_print(N, M, C, V, ao, occ, spin, ps, fo);
  return;
}

