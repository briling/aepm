#include "matrix.h"

void mx_id(unsigned int n, double * a){
  unsigned int i,j;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      a[i*n+j] = (i==j ? 1.0 : 0.0);
    }
  }
  return;
}

void mx_print(unsigned int n, double * a, FILE   * f){
  unsigned int i,j;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      fprintf(f, "%25.15lf", a[i*n+j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");
  fflush(f);
  return;
}

void mx_nosym_print(unsigned int n, double * a, FILE   * f){
  unsigned int i,j;
  for(i=0; i<n; i++){
    for(j=0; j<=i; j++){
      fprintf(f, "%25.15lf", a[mpos(j,i)]);
    }
    for(j=i+1; j<n; j++){
      fprintf(f, "%25.15lf", a[mpos(i,j)]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");
  fflush(f);
  return;
}

void mx_transp(unsigned int n, double * a){
  unsigned int  i,j;
  double t;
  for(i=0; i<n; i++){
    for(j=i+1; j<n; j++){
      t = a[i*n+j];
      a[i*n+j] = a[j*n+i];
      a[j*n+i] = t;
    }
  }
  return;
}

void mx_symmultrectmx(unsigned int n, unsigned int m, double * p, double * a, double * b){
  // P(n*m) = A(n*n)*B(n*m)
  unsigned int i,j,k;
  double t;
  for(i=0; i<n; i++){
    for(j=0; j<m; j++){
      t=0.0;
      for(k=0; k<i; k++){
        t += a[mpos(k,i)] * b[k*m+j];
      }
      for(k=i; k<n; k++){
        t += a[mpos(i,k)] * b[k*m+j];
      }
      p[i*m+j] = t;
    }
  }
  return;
}

void mx_multsymmx(unsigned int n, double * p, double * a, double * b){
  unsigned int i,j,k;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      double t = 0.0;
      for(k=0; k<n; k++){
        t += a[i*n+k] * b[MPOSIF(k,j)];
      }
      p[i*n+j] = t;
    }
  }
  return;
}

void mx_BHBt_sym(unsigned int n, double * h, double * b){
  /* H := BHB^T */
  unsigned int i,j,k;
  double   s;
  double * t = malloc(n*n*sizeof(double));
  if(!t){
    abort();
  }
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      s = 0.0;
      for(k=0; k<=j; k++){
        s += b[i*n+k] * h[mpos(k,j)];
      }
      for(k=j+1; k<n; k++){
        s += b[i*n+k] * h[mpos(j,k)];
      }
      t[i*n+j] = s;
    }
  }
  for(i=0; i<n; i++){
    for(j=i; j<n; j++){
      s = 0.0;
      for(k=0; k<n; k++){
        s += b[i*n+k] * t[j*n+k];
      }
      h[mpos(i,j)] = s;
    }
  }
  free(t);
  return;
}

void mx_BHBt_sym2(unsigned int n, unsigned int m,
                  double * r, double * h, double * b){
  /* R := BHB^T
   * R: n*n
   * H: m*m
   * B: n*m
   */

  double * bh = malloc(n*m*sizeof(double));
  if(!bh){
    abort();
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      double t = 0.0;
      for(int k=0; k<j; k++){
        t += b[m*i+k] * h[mpos(k,j)];
      }
      for(int k=j; k<m; k++){
        t += b[m*i+k] * h[mpos(j,k)];
      }
      bh[i*m+j] = t;
    }
  }

  for(int i=0; i<n; i++){
    for(int j=i; j<n; j++){
      double t = 0.0;
      for(int k=0; k<m; k++){
        t += b[m*j+k] * bh[i*m+k];
      }
      r[mpos(i,j)] = t;
    }
  }

  free(bh);
  return;
}

void mx_BdiagBt(unsigned int n, double * h, double * v, double * b){
  /* H := BVB^T, V is diagonal */
  unsigned int i,j,k;
  for(i=0; i<n; i++){
    for(j=i; j<n; j++){
      double s = 0;
      for(k=0; k<n; k++){
        s += b[i*n+k] * v[k] * b[j*n+k];
      }
      h[mpos(i,j)] = s;
    }
  }
  return;
}

void mx_sqr2(unsigned int m, unsigned int n,
             double * p, double * a){
  /* a: m*n
   * p = a*(a^T): m*m
   */
  unsigned int i,j;
  for(i=0; i<m; i++){
    for(j=i; j<m; j++){
      p[mpos(i,j)] = vecdot(n, a+i*n, a+j*n);
    }
  }
  return;
}

void mx_multvec(unsigned int n, unsigned int m,
                double * r, double * a, double * v){
  /*  r(n*1) = A(n*m) * v(m*1) */
  double t;
  unsigned int i, j;
  for(i=0; i<n; i++){
    t = 0.0;
    for(j=0; j<m; j++){
      t += v[j]*a[m*i+j];
    }
    r[i] = t;
  }
  return;
}

double mx_vecdot(unsigned int n, double * u, double * a, double * v){
  /* (u+)*A*v,  A=A^T */
  unsigned int i, j;
  double s = 0.0;
  for(i=0; i<n; i++){
    double tv = 0.0;
    double tu = 0.0;
    for(j=0; j<i; j++){
      tv += u[j] * a[mpos(j,i)];
      tu += v[j] * a[mpos(j,i)];
    }
    s += v[i]*u[i]*a[mpos(i,i)] + v[i]*tv + u[i]*tu;
  }
  return s;
}

double mx_vecdot_nosym(unsigned int n, double * u, double * a, double * v){
  /* (u+)*A*v */
  unsigned int i, j;
  double s = 0.0;
  for(i=0; i<n; i++){
    double av = 0.0;
    for(j=0; j<n; j++){
      av += a[i*n+j]*v[j];
    }
    s += u[i]*av;
  }
  return s;
}

