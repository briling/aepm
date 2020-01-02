#include "mol.h"

void mol_print_m(mol * m, int bohr, const char s[], FILE * f){
  fprintf(f, "%s$molecule\n", s);
  if(m->z){
    fprintf(f, "%s charge=%d\n", s, m->z);
  }
  if(m->mult != 1){
    fprintf(f, "%s mult=%d\n", s, m->mult);
  }
  if(bohr){
    fprintf(f, "%s unit=b\n", s);
  }
  fprintf(f, "%s cartesian\n", s);

  for(int i=0; i<m->n; i++){
    if(bohr){
      fprintf(f, "%s%4d%15.8lf%15.8lf%15.8lf",
          s, m->q[i], m->r[3*i], m->r[3*i+1], m->r[3*i+2]);
    }
    else{
      fprintf(f, "%s%4d%15.8lf%15.8lf%15.8lf",
          s, m->q[i], m->r[3*i]*BA, m->r[3*i+1]*BA, m->r[3*i+2]*BA);
    }
    if(m->m[i] > 0){
      fprintf(f, "   mass=%lf", m->m[i]);
    }

    if (m->s[i][0] != 0){
      fprintf(f, "   type=%s", m->s[i]);
    }

    int a = m->l[i];
    int b = m->l[i+1];
    if(b-a > 0){
      fprintf(f, "   k=");
      for(int j=a; j<b; j++){
        fprintf(f, "%d", m->k[j]+1);
        if (m->b[j]!=1){
          fprintf(f, "(%d)", m->b[j]);
        }
        if(j!=b-1){
          fprintf(f, ",");
        }
      }
    }
    fprintf(f, "\n");
  }
  fprintf(f, "%s$end\n\n", s);
  fflush (f);
  return;
}

void mol_print_xyz(mol * m, const char comment[], FILE * f){
  const char aname[][3]={
    #include "elements.h"
  } ;
  fprintf(f, "%d\n", m->n);
  fprintf(f, "%s\n", comment);
  for(int i=0; i<m->n; i++){
    fprintf(f, "%3s % 10.5lf % 10.5lf % 10.5lf\n",
        aname[m->q[i]], (m->r[3*i])*BA, (m->r[3*i+1])*BA, (m->r[3*i+2])*BA);
  }
  fflush (f);
  return;
}

