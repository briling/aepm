#include "vecn.h"

void vecset(int n, double * r, double s){
  int i;
  for (i=0; i<n; i++){
    r[i] = s;
  }
  return;
}

void vecprint(int n, double * r, const char * s, FILE * f){
  int i;
  for (i=0; i<n; i++){
    fprintf(f, "%25.15lf%s", r[i], s);
  }
  fprintf(f, "\n");
  fflush(f);
  return;
}

void veccp(int n, double * u, double * v){
  int i;
  for (i=0; i<n; i++){
    u[i] = v[i];
  }
  return;
}

double vecdot(int n, double * u, double * v){
  double s = 0;
  int i;
  for (i=0; i<n; i++){
    s += u[i] * v[i];
  }
  return s;
}

void vecadd(int n, double * u, double * v){
  int i;
  for (i=0; i<n; i++){
    u[i] += v[i];
  }
  return;
}

void vecadds(int n, double * u, double * v, double s){
  int i;
  for (i=0; i<n; i++){
    u[i] += v[i]*s;
  }
  return;
}

void vecscal(size_t n, double * u, double s){
  size_t i;
  for (i=0; i<n; i++){
    u[i] *= s;
  }
  return;
}

void vecsums(int n, double * w, double * u, double * v, double s){
  int i;
  for (i=0; i<n; i++){
    w[i] = u[i]+v[i]*s;
  }
  return;
}

double vecabsmax(int n, double * u){
  int    i;
  double t;
  double s = 0.0;
  for (i=0; i<n; i++){
    t = fabs(u[i]);
    if (s<t){
      s = t;
    }
  }
  return s;
}

