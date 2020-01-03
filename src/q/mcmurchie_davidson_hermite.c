#include "q.h"
#include "integrals.h"

static void def_all_inner(
    int n1, int n3,
    double PAx, double p2_1,
    double arr0[(L_MAX+1)*2+1]
    ){
  double arr1[(L_MAX+1)*2+1]={};
  for(int j=0; j<n1; j++){
    for(int n=1; n<n3; n++){
      arr1[n]    = PAx    * arr0[n]
                 + p2_1   * arr0[n-1]
                 + (n+1)  * arr0[n+1];
    }
    arr1[0]      = PAx    * arr0[0]
                 +          arr0[1];
    arr1[n3]     = PAx    * arr0[n3]
                 + p2_1   * arr0[n3-1];
    for(int n=0; n<=n3; n++){
      arr0[n] = arr1[n];
    }
  }
  return;
}

static void def_all(int n1, int n2,
    double PAx, double PBx, double p2_1,
    double ret[]){

  double arr0[(L_MAX+1)*2+1]={};
  arr0[0] = 1.0;
  def_all_inner(n2, n2,    PBx, p2_1, arr0);
  def_all_inner(n1, n1+n2, PAx, p2_1, arr0);

  for(int n=0; n<n1+n2; n++){
    ret[n] = arr0[n];
  }
  return;
}

static double def(int n1, int n2, int N, double PAx, double PBx, double p2_1) __attribute__ ((unused));
static double def(int n1, int n2, int N, double PAx, double PBx, double p2_1){

  if(N<0 || N>n1+n2 || n1<0 || n2<0){
    return 0.0;
  }

  double arr0[(L_MAX+1)*2+1]={};
  double arr1[(L_MAX+1)*2+1]={};
  arr0[0] = 1.0;

  int N1 = MAX(N-n1, 0    );
  int N2 = MIN(N+n1, n1+n2);
  for(int j=1; j<=n2; j++){

    int N3 = MAX(N1-(n2-j), 1   );
    int N4 = MIN(N2+(n2-j), n2-1);

    for(int Nk=N3; Nk<=N4; Nk++){
      arr1[Nk] = PBx    * arr0[Nk]
               + p2_1   * arr0[Nk-1]
               + (Nk+1) * arr0[Nk+1] ;
    }
    if(N1 <= n2-j){
      arr1[0] = PBx     * arr0[0]
              +           arr0[1];
    }
    if(N2 >= j){
      arr0[n2] = PBx    * arr0[n2]
               + p2_1   * arr0[n2-1];
    }
    arr0[0] = arr1[0];
    for(int i=N3; i<=N4; i++){
       arr0[i] = arr1[i];
    }

  }

  for(int j=1; j<=n1; j++){

    int N1 = MAX(N-(n1-j), 1     );
    int N2 = MIN(N+(n1-j), n1+n2-1);

    for(int Ni=N1; Ni<=N2; Ni++){
      arr1[Ni]    = PAx    * arr0[Ni]
                  + p2_1   * arr0[Ni-1]
                  + (Ni+1) * arr0[Ni+1];
    }
    if(N <= n1-j){
      arr1[0]     = PAx    * arr0[0]
                  +          arr0[1];
    }
    if(N >= n2+j){
      arr0[n1+n2] = PAx    * arr0[n1+n2]
                  + p2_1   * arr0[n1+n2-1];
    }
    arr0[0] = arr1[0];
    for(int i=N1; i<=N2; i++){
      arr0[i] = arr1[i];
    }

  }
  return arr1[N];
}

static double def0_inner(int n1, int n2, double PAx, double PBx, double p2_1){

  // n1 >= n2

  double arr0[(L_MAX+1)*2+1]={};
  double arr1[(L_MAX+1)*2+1]={};
  arr0[0] = 1.0;

  for(int j=0; j<n2; j++){
    for(int Nk=1; Nk<n2; Nk++){
      arr1[Nk] = PBx    * arr0[Nk]
               + p2_1   * arr0[Nk-1]
               + (Nk+1) * arr0[Nk+1] ;
    }
    arr1[0]    = PBx    * arr0[0]
               +          arr0[1];
    arr1[n2]   = PBx    * arr0[n2]
               + p2_1   * arr0[n2-1];
    for(int n=0; n<=n2; n++){
       arr0[n] = arr1[n];
    }
  }

  for(int j=0; j<n1; j++){
    arr1[0]       = PAx    * arr0[0]
                  +          arr0[1];
    int n11 = n1-j;
    for(int n=1; n<n11; n++){
      arr1[n]     = PAx    * arr0[n]
                  + p2_1   * arr0[n-1]
                  + (n+1)  * arr0[n+1];
    }
    for(int n=0; n<n11; n++){
      arr0[n] = arr1[n];
    }
  }
  return arr1[0];
}

static double def_r(int n, int n1, int N, double PAx, double PBx, double p2_1) __attribute__ ((unused));
static double def_r(int n, int n1, int N, double PAx, double PBx, double p2_1){
  if(N<0 || N>n+n1){
    return 0.0;
  }
  if( (n==0) && (n1==0) && (N==0) ){
    return 1.0;
  }
  if( n1==0 ){
    return p2_1  * def_r(n-1, n1, N-1, PAx, PBx, p2_1) +
           PAx   * def_r(n-1, n1, N,   PAx, PBx, p2_1) +
           (N+1) * def_r(n-1, n1, N+1, PAx, PBx, p2_1);
  }
  else {
    return p2_1  * def_r(n, n1-1, N-1, PAx, PBx, p2_1) +
           PBx   * def_r(n, n1-1, N,   PAx, PBx, p2_1) +
           (N+1) * def_r(n, n1-1, N+1, PAx, PBx, p2_1);
  }
}

double def0(int n1, int n2, double PAx, double PBx, double p21){

  int n = n1+n2;
  if(n==0){
    return 1.0;
  }

  if(n1<n2){
    int    ti = n1;  n1  = n2;  n2  = ti;
    double td = PAx; PAx = PBx; PBx = td;
  }

  switch(n){
    case 1:
      return PAx;
    case 2:
      if(n1==1){
        return PBx * PAx  + p21;
      }
      else if(n1==2){
        return PAx * PAx + p21;
      }
    case 3:
      if(n1==2){
        return PAx * PAx * PBx + p21 * (PBx + 2.0 * PAx);
      }
      else if(n1==3){
        return PAx * (PAx * PAx + 3.0 * p21);
      }
    case 4:
      if(n1==2){
        return (PAx * PAx * PBx * PBx + p21 * PAx * PAx + p21 * PBx * PBx + 4.0 * p21 * PAx * PBx + 3.0 * p21 * p21);
      }
      break;
  }

  return def0_inner(n1, n2, PAx, PBx, p21);
}

void filldef1(int n1, int n2, double d[L_MAX*2+1],
              double PAx, double PBx, double p21){
  int n = n1+n2;
  if(n==0){
    d[0] = 1.0;
    return;
  }
  d[n] = intpow(p21, n);

  if(n1<n2){
    int    ti = n1;  n1  = n2;  n2  = ti;
    double td = PAx; PAx = PBx; PBx = td;
  }

  switch(n){
    case 1:
      d[0] = PAx;
      return;
    case 2:
      if(n1==1){
        d[0] =  PBx * PAx  + p21;
        d[1] = (PAx + PBx) * p21;
      }
      else if(n1==2){
        d[0] = PAx * PAx + p21;
        d[1] = 2.0 * PAx * p21;
      }
      return;
    case 3:
      if(n1==2){
        d[0] =  PAx * PAx * PBx + p21 * (PBx + 2.0 * PAx);
        d[1] =  p21 * (PAx * (PAx + 2.0 * PBx) + 3.0 * p21);
        d[2] =  (PBx + 2.0*PAx) * p21 * p21;
      }
      else if(n1==3){
        d[0] = PAx * (PAx * PAx + 3.0 * p21);
        d[1] = 3.0 * p21 * (PAx * PAx + p21);
        d[2] = 3.0 * PAx * p21 * p21;
      }
      return;
    case 4:
      if(n1==2){
        d[0] = (PAx * PAx * PBx * PBx + p21 * PAx * PAx + p21 * PBx * PBx + 4.0 * p21 * PAx * PBx + 3.0 * p21 * p21);
        d[1] = 2.0 * PAx * PAx * PBx * p21 + 2.0 * PAx * PBx * PBx * p21 + 6.0 * PBx * p21 * p21 + 6.0 * PAx * p21 * p21;
        d[2] = p21 * p21 * (PAx * PAx + PBx * PBx + 4.0 * PAx * PBx + 6.0 * p21);
        d[3] = 2.0 * p21 * p21 * p21 * (PAx + PBx);
        return;
      }
      break;
  }
#if 0
  for(int N=0; N<n; N++){
    d[N] = def(n1, n2, N, PAx, PBx, p21);
  }
#else
  def_all(n1, n2, PAx, PBx, p21, d);
#endif
  return;
}

void filldef11(int n1, int n2, double d[L_MAX*2+1],
              double PAx, double PBx, double p21){
  if(n1<n2){
    int    ti = n1;  n1  = n2;  n2  = ti;
    double td = PAx; PAx = PBx; PBx = td;
  }
  def_all(n1, n2, PAx, PBx, p21, d);
  return;
}

void filldef_all(int l1, int l2,
    double D1[L_MAX+1][L_MAX+1][L_MAX*2+1],
    double PAx, double PBx, double p21){

  double pn[2*L_MAX+1]={1.0,};
  for(int i=0; i<l1+l2; i++){
    pn[i+1] = pn[i] * p21;
  }

  for(int n2=0; n2<=l2; n2++){
    double arr1[(L_MAX+1)*2+1]={};
    double arr0[(L_MAX+1)*2+1]={};
    arr0[0] = 1.0;
    def_all_inner(n2, n2,    PBx, p21, arr0);
    for(int n1=0; n1<=l1; n1++){
      int n12 = n1+n2;
      D1[n1][n2][n12] = pn[n12];
      memcpy(arr1, arr0, sizeof(arr0));
      def_all_inner(n1, n12, PAx, p21, arr1);
      for(int n=0; n<n12; n++){
        D1[n1][n2][n] = arr1[n];
      }
    }
  }
  return;
}

static void def_all_inner_self(
    int n1, int n3, double p2_1,
    double arr0[(L_MAX+1)*2+1]
    ){
  double arr1[(L_MAX+1)*2+1]={};
  for(int j=0; j<n1; j++){
    for(int n=1; n<n3; n++){
      arr1[n]    = p2_1   * arr0[n-1]
                 + (n+1)  * arr0[n+1];
    }
    arr1[0]      =          arr0[1];
    arr1[n3]     = p2_1   * arr0[n3-1];
    for(int n=0; n<=n3; n++){
      arr0[n] = arr1[n];
    }
  }
  return;
}

void filldef1_self(int n1, int n2, double d[L_MAX*2+1], double p21){
  int n = n1+n2;
  if(n==0){
    d[0] = 1.0;
    return;
  }
  d[n] = intpow(p21, n);
  if(n1<n2){
    int ti = n1; n1 = n2; n2 = ti;
  }

  switch(n){
    case 1:
      d[0] = 0.0;
      return;
    case 2:
      d[0] = p21;
      d[1] = 0.0;
      return;
    case 3:
      d[0] = 0.0;
      d[1] = 3.0 * p21 * p21;
      d[2] = 0.0;
      return;
    case 4:
      if(n1==2){
        d[0] = 3.0 * p21 * p21;
        d[1] = 0.0;
        d[2] = 6.0 * p21 * p21 * p21;
        d[3] = 0.0;
        return;
      }
      break;
  }

  double arr0[(L_MAX+1)*2+1]={};
  arr0[0] = 1.0;
  def_all_inner_self(n2, n2,    p21, arr0);
  def_all_inner_self(n1, n1+n2, p21, arr0);

  for(int n=0; n<n1+n2; n++){
    d[n] = arr0[n];
  }
  return;
}

