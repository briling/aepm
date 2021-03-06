#ifndef MOL_H
#define MOL_H

#include "common.h"

typedef struct {
  styp   * s;
  double * r;
  int    * q;
  double * m;
  int    * k;
  int    * b;
  int    * l;
  int      n;
  int      c;
  int      z;
  int      mult;
} mol;

/*==========================================================================================*/

mol  * mol_read     (FILE   * f);
void   mol_print_m  (mol    * m,  int bohr, const char s[], FILE   * f);
void   mol_print_xyz(mol * m, const char comment[], FILE * f);

#endif

