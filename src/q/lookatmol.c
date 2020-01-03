#include "q.h"

int elnumber(mol * m){
  int N = -m->z;
  int i;
  for(i=0; i<m->n; i++){
    N += m->q[i];
  }
  if(N<=0){
    fprintf(stderr, "\tN = %d ?\n", N);
    exit(1);
  }
  return N;
}

int corepairs_simple(mol * m){
  int core = 0;
  for(int i=0; i<m->n; i++){
    int q = m->q[i];
    if(q<=2){
      continue;
    }
    else if(q<=10){
      core += 2;
    }
    else if(q<=18){
      core += 10;
    }
    else if(q<=30){
      core += 18;
    }
    else if(q<=36){
      core += 28;
    }
    else if(q<=48){
      core += 36;
    }
    else if(q<=54){
      core += 46;
    }
    else if(q<=70){
      core += 54;
    }
    else if(q<=80){
      core += 68;
    }
    else if(q<=86){
      core += 78;
    }
    else if(q<=103){
      core += 86;
    }
  }
  return core/2;
}

