#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void vecset(int n, double * r, double s);
void vecprint(int n, double * r, const char * s, FILE * f);
void veccp(int n, double * u, double * v);
double vecdot(int n, double * u, double * v);
void vecadd(int n, double * u, double * v);
void vecadds(int n, double * u, double * v, double s);
void vecscal(size_t n, double * u, double s);
void vecsums(int n, double * w, double * u, double * v, double s);
double vecabsmax(int n, double * u);

