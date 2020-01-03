#include <unistd.h>
#include "q.h"
#include "qinit.h"

int pars_read(double a0[NATOMS], int bp[NATOMS], int fixflag[NATOMS], FILE * fl){
  double ta;
  unsigned int i1, i2;
  int npar = 0;
  while(fscanf(fl, " fix %lf %u - %u", &ta, &i1, &i2) || fscanf(fl, "%lf %u - %u", &ta, &i1, &i2)){
    npar++;
  }
  rewind(fl);
  int npar1 = 0;

  for(int j=0; j<npar; j++){
    int i = j+npar1;
    int ff = 0;
    double a;
    const double a_def = 0.125;

    int c = fscanf(fl, "%lf %u - %u", &a, &i1, &i2);
    if(!c){
      c = fscanf(fl, " fix %lf %u - %u", &a, &i1, &i2);
      ff = 1;
    }

    int i10 = j ? bp[i-1] : 1;
    if(i1 != i10){
      if(i1<i10){
        GOTOHELL;
      }
      a0[i] = a_def;
      fixflag[i] = 0;
      bp[i] = i1;
      npar1++;
      i++;
    }

    a0[i] = a;
    fixflag[i] = ff;
    bp[i] = (c==3 ? i2:i1) + 1;

    if( (j==npar-1) && (c==3?i2:i1)!=NATOMS ){
      a0[i+1] = a_def;
      fixflag[i+1] = 0;
      bp[i+1] = NATOMS+1;
      npar1++;
    }

  }
  fscanf(fl, "\n");
  return npar+npar1;
}

int pars_default(int mode, double a0[NATOMS], int bp[NATOMS]){
  if(mode==1){
    bp[ 0] =  3;  a0[ 0] = 1.0 /  3.0;
    bp[ 1] =  5;  a0[ 1] = 1.0 / 16.0;
    bp[ 2] = 11;  a0[ 2] = 1.0 /  3.0;
    bp[ 3] = 13;  a0[ 3] = 1.0 / 32.0;
    bp[ 4] = 19;  a0[ 4] = 1.0 /  8.0;
    bp[ 5] = 21;  a0[ 5] = 1.0 / 32.0;
    bp[ 6] = 31;  a0[ 6] = 1.0 /  6.0;
    bp[ 7] = 37;  a0[ 7] = 1.0 / 12.0;
    bp[ 8] = 39;  a0[ 8] = 1.0 / 32.0;
    bp[ 9] = 49;  a0[ 9] = 1.0 /  8.0;
    bp[10] = 55;  a0[10] = 1.0 / 12.0;
    bp[11] = 71;  a0[11] = 1.0 / 32.0;
    bp[12] = 87;  a0[12] = 1.0 / 12.0;
    bp[13] =103;  a0[13] = 1.0 / 32.0;
    return 14;
  }
  else if(mode==2){
    for(int i=0; i<NATOMS; i++){
      bp[i] = i+2;
      a0[i] = 0.125;
    }
    return NATOMS;
  }
  else{
    bp[0] = NATOMS+1; a0[0] = 0.125;
    return 1;
  }
}

static int find_aim(int q, int nmol, mol_data * md){
  for(int i=0; i<nmol; i++){
    mol * m = md[i].m;
    for(int k=0; k<m->n; k++){
      if(m->q[k] == q){
        return 1;
      }
    }
  }
  return 0;
}

static void used_pars_print(int pn, int nmol, mol_data * md,
    const char aname[][3], int atpar[NATOMS+1],
    int * up, double * a, FILE * fo){

  for(int i=0; i<pn; i++){
    fprintf(fo, "p%d(%d)=%10.9lf: ", i, up[i], a[i]);
    for(int j=1; j<=NATOMS; j++){
      if(atpar[j] == up[i]){
        if(find_aim(j, nmol, md)){
          fprintf(fo, "%s ", aname[j]);
        }
      }
    }
    fprintf(fo, "\n");
  }
  return;
}

static void used_pars(int nparams,
    int nmol, mol_data * md,
    const char aname[][3],
    int fixflag[NATOMS], int atpar[NATOMS+1],
    int * npar, int * nfpar,
    int * up, int * fup,
    double * a, FILE * fo){
  int istty = fo ? isatty(fileno(fo)) : 0;
  const char green[] = "\e[1;32m";
  const char red[]   = "\e[1;31m";
  int np = 0;
  int nf = 0;
  for(int i=0; i<nparams; i++){
    fprintf_if(fo, "p%d=%10.9lf: ", i, a[i]);
    int pused = 0;
    for(int j=1; j<=NATOMS; j++){
      if(atpar[j] == i){
        if(find_aim(j, nmol, md)){
          if(!fixflag[i]){
            pused = 1;
          }
          else{
            pused = -1;
          }
          const char * color   = fixflag[i] ? red:green;
          const char * nocolor = fixflag[i] ? "/":"'";
          fprintf_if(fo, "%s%s%s ", istty?color:nocolor, aname[j], istty?"\e[0m":nocolor);
        }
        else{
          fprintf_if(fo, "%s ", aname[j]);
        }
      }
    }
    if(pused==1){
      up[np] = i;
      np++;
    }
    else if(pused==-1){
      fup[nf] = i;
      nf++;
    }
    fprintf_if(fo, "\n");
  }
  fprintf_if(fo, "\n");
  *npar  = np;
  *nfpar = nf;
  return;
}

int atpar_fill(int nparams, int nmol,
    double fa[NATOMS],
    double a[NATOMS],
    double a0[NATOMS],
    int bp[NATOMS], int fixflag[NATOMS],
    int * nfpar_r,
    int fatpar[NATOMS+1],
    int  atpar[NATOMS+1],
    mol_data * md, FILE * fo){

  int atpar0[NATOMS+1];
  for(int j=0; j<nparams; j++){
    for(int i=(j ? bp[j-1] : 1); i<bp[j]; i++){
      atpar0[i] = j;
    }
  }

  const char aname[][3]={
    #include "elements.h"
  };

  int up[NATOMS], fup[NATOMS];
  int npar, nfpar;
  used_pars(nparams, nmol, md, aname, fixflag, atpar0, &npar, &nfpar, up, fup, a0, fo);
  for(int i=0; i<npar; i++){
    a[i] = a0[up[i]];
  }
  for(int i=0; i<nfpar; i++){
    fa[i] = a0[fup[i]];
  }
  if(fo){
    used_pars_print(npar,  nmol, md, aname, atpar0,  up,  a, fo);
    if(nfpar){
      fprintf(fo, "fixed:\n");
      used_pars_print(nfpar, nmol, md, aname, atpar0, fup, fa, fo);
    }
    fprintf(fo, "\n");
    fflush(fo);
  }

  for(int i=1; i<=NATOMS; i++){
    atpar[i] = -1;
    fatpar[i] = -1;
    for(int j=0; j<npar; j++){
      if(atpar0[i] == up[j]){
        atpar[i] = j;
        continue;
      }
    }
    if(atpar[i]>=0){
      continue;
    }
    for(int j=0; j<nfpar; j++){
      if(atpar0[i] == fup[j]){
        fatpar[i] = j;
        continue;
      }
    }
  }
#if 0
  for(int i=1; i<=NATOMS; i++){
    printf("%3d: %d -> % d\n", i, atpar0[i], atpar[i]);
  }
  for(int i=1; i<=NATOMS; i++){
    printf("%3d: %d -> % d\n", i, atpar0[i], fatpar[i]);
  }
#endif
  *nfpar_r = nfpar;
  return npar;
}

void conv_par_print(int atpar[NATOMS+1], int fatpar[NATOMS+1],
                    double * a, double * fa, FILE * fo){
  int fullpar[NATOMS+1];
  for(int i=1; i<=NATOMS; i++){
    int i1 =  atpar[i]+1;
    int i2 = fatpar[i]+1;
    fullpar[i] = i1>i2?i1:-i2;
  }
  for(int i=1; i<=NATOMS; i++){
    int j = fullpar[i];
    if(!j) continue;
    int first = ( i==1      || j!=fullpar[i-1] );
    int last  = ( i==NATOMS || j!=fullpar[i+1] );
    if(first){
      fprintf(fo, "   %s%17.15lf  ", j>0?"    ":"fix ", j>0?a[j-1]:fa[-j-1]);
    }
    if(first&&!last){
      fprintf(fo, "%d-", i);
    }
    else if(last){
      fprintf(fo, "%d\n", i);
    }
  }
  fprintf(fo, "\n");
  return;
}
