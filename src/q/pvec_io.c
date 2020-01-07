#include "q.h"
#include <stdint.h>

/* In sake of compatibility with Priroda,
 * we change the order of basis functions
 * when write/read MO coefficient files.
 *
 * order in memory (m)  | order in file (m1)
 * ----------------------------------------
 *          -1          |        0
 *           0          |       +1
 *          +1          |       -1
 * ----------------------------------------
 *          -2          |        0
 *          -1          |       +1
 *           0          |       -1
 *          +1          |       +2
 *          +2          |       -2
 * -----------------------------------------
 *
 * convert from file:
 * > veccp( n, C + n*(i+(m1-m)), Ct + n*i );
 * convert to file:
 * > veccp( n, Ct + n*i, C + n*(i+(m1-m)) );
 *
 * we can create an array of m1's:
 * > int M1[L_MAX*2+1] = { 0, +1, -1, +2, -2};
 * and take values
 * > m1 = M1[l+m];
 * but
 * they can be calculated:
 * n     : 0, 1, 2, 3, 4, 5, 6, ...
 * M1[n] : 0,+1,-1,+2,-2,+3,-3, ...
 * thus
 * M1[n] = (n%2) ? (n+1)/2 : -(n+1)/2
 */

static inline int M1(int n){
  return (n%2) ? (n+1)/2 : -(n+1)/2;
}

void vec_to_p(unsigned int n, atomo * ao, double * Ct, double * C){
  mx_transp(n, C);
  for(unsigned int i=0; i<n; i++){
    int l  = ao[i].l;
    int m  = ao[i].m;
    int m1 = M1(m+l);
    veccp( n, Ct + n*i, C + n*(i+(m1-m)) );
  }
  mx_transp(n, Ct);
  mx_transp(n, C);
  return;
}

static void vec_from_p(unsigned int n, atomo * ao, double * C, double * Ct){
  mx_transp(n, Ct);
  for(unsigned int i=0; i<n; i++){
    int l  = ao[i].l;
    int m  = ao[i].m;
    int m1 = M1(m+l);
    veccp( n, C + n*(i+(m1-m)), Ct + n*i );
  }
  mx_transp(n, C);
  return;
}

int pvec_read(int M, atomo * ao, double * Va, double * Ca, double * Vb, double * Cb, const char s[]){

  FILE * f;
  if( !(f = fopen(s, "r"))){
    return 0;
  }

  uint32_t n;
  if( !fread(&n, sizeof(n), 1, f) || n!=M ){
    fclose(f);
    return 0;
  }

  size_t vsize = sizeof(double)*n;
  size_t csize = sizeof(double)*n*n;
  double * Ct = malloc(csize);

  if( !fread(Va, vsize, 1, f) || !fread(Ct, csize, 1, f) ){
    fclose(f);
    free(Ct);
    return 0;
  }
  vec_from_p(n, ao, Ca, Ct);

  if( Vb && Cb ){
    if( !fread(Vb, vsize, 1, f) || !fread(Ct, csize, 1, f) ){
      fclose(f);
      free(Ct);
      return 0;
    }
    vec_from_p(n, ao, Cb, Ct);
  }

  free(Ct);
  fclose(f);
  return 1;
}

int pvec_write(int M, atomo * ao, double * Va, double * Ca, double * Vb, double * Cb, const char s[]){

  FILE * f;
  if( !(f = fopen(s, "w"))){
    return 0;
  }

  int32_t n = M;
  size_t vsize = sizeof(double)*n;
  size_t csize = sizeof(double)*n*n;

  double * Cta = malloc(csize);
  double * Ctb = malloc(csize);
  vec_to_p(n, ao, Cta, Ca);
  vec_to_p(n, ao, Ctb, Cb);

  if( !fwrite(&n, sizeof(n), 1, f) ||
      !fwrite(Va, vsize, 1, f) || !fwrite(Cta, csize, 1, f) ||
      !fwrite(Vb, vsize, 1, f) || !fwrite(Ctb, csize, 1, f) ){
    free(Cta);
    free(Ctb);
    fclose(f);
    return 0;
  }

  free(Cta);
  free(Ctb);

  fclose(f);
  return 1;
}

int qvec_write(int M, double * Va, double * Ca, double * Vb, double * Cb, const char s[]){

  FILE * f;
  if( !(f = fopen(s, "w"))){
    return 0;
  }

  int32_t n = M;
  size_t vsize = sizeof(double)*n;
  size_t csize = sizeof(double)*n*n;

  if( !fwrite(&n, sizeof(n), 1, f) ||
      !fwrite(Va, vsize, 1, f) || !fwrite(Ca, csize, 1, f) ||
      !fwrite(Vb, vsize, 1, f) || !fwrite(Cb, csize, 1, f) ){
    fclose(f);
    return 0;
  }

  fclose(f);
  return 1;
}

