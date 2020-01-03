#include "common.h"
#include "boys.h"

/* J. Chem. Phys. 84, 3963 (1986) */

static inline double boys0(double z){
  if(z<EPS1){
    return 1.0;
  }
  double t = sqrt(z);
  return 0.5 * SQRTPI / t * erf(t);
}

static inline double boys_next(double z, int n, double prev, double ez){
  if(z<EPS1){
    return 1.0/(2.0*n+1.0);
  }
  double z1 = 1.0/z;
  return ((2.0*n-1.0)*prev-ez)*z1*0.5;
}

double * boys_fill(void){
  double f, ex, x;
  int size = BOS_N*(6+BOS_NMAX+1);
  double * boys_array = malloc(sizeof(double)*size);
  for(int i=0; i<BOS_N; i++){
    x = i*BOS_D;
    f = boys0(x);
    boys_array[i*(6+BOS_NMAX+1)] = f;
    ex = exp(-x);
    for(int k=1; k<=6+BOS_NMAX; k++){
      f = boys_next(x, k, f, ex);
      boys_array[i*(6+BOS_NMAX+1)+k] = f;
    }
  }
  return boys_array;
}

/*-------------------------------------------------------------------------------------------------*/

static const double Tf[] = { 33.0, 37.0, 41.0, 43.0, 46.0,
                             49.0, 51.0, 54.0, 56.0, 58.0,
                             61.0, 63.0, 66.0, 68.0, 70.0,
                             72.0, 74.0 };

static void boys_os_small(double f[], unsigned int m, unsigned int n, double T, double * boys_array){

  if(n > BOS_NMAX){
    GOTOHELL;
  }

  int i = (int)(round(T*BOS_D1));
  if (i>=BOS_N){
    return;
  }

  double dT = BOS_D*i-T;
  int j = i*(6+BOS_NMAX+1);

  double p[6];
  p[0] = dT;
  p[1] = dT*p[0]*0.5;
  p[2] = dT*p[1]*0.33333333333333333333;
  p[3] = dT*p[2]*0.25;
  p[4] = dT*p[3]*0.2;
  p[5] = dT*p[4]*0.16666666666666666666;

  for(int k=m; k<=n; k++){
    int l = j+k;
    f[k]  = boys_array[l];
    f[k] += boys_array[l+1] * p[0];
    f[k] += boys_array[l+2] * p[1];
    f[k] += boys_array[l+3] * p[2];
    f[k] += boys_array[l+4] * p[3];
    f[k] += boys_array[l+5] * p[4];
    f[k] += boys_array[l+6] * p[5];
  }

  return;
}

static void boys_os_big(double f[], unsigned int k, double T){
  double T1 = 1.0/T;
  double p  = SQRTPI*0.5*sqrt(T1);
  f[0] = p;
  for(unsigned int i=1; i<k; i++){
    p *= 0.5*T1 * (2*i-1);
    f[i] = p;
  }
  return;
}

void boys_os(double f[], unsigned int n, double T, double * boys_array){
  unsigned int k;
  for(k=0; k<=n; k++){
    if (T <= Tf[k]){
      break;
    }
  }
  if(k){
    boys_os_big(f, k, T);
  }
  boys_os_small(f, k, n, T, boys_array);
  return;
}

