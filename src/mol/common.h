#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define   BA    0.5291772
#define   AB    1.88972616356109068947

#define MAX(x,y) ( ((x) > (y)) ? (x) : (y) )
#define MIN(x,y) ( ((x) < (y)) ? (x) : (y) )

#define GOTOHELL { \
      fprintf(stderr, "%s:%d %s() -- ", \
      __FILE__, __LINE__, __FUNCTION__); \
      abort(); }

#define printalive printf("alive @ %s:%d\n", __FILE__, __LINE__);

#define fprintf_if(F, ...) if(F) fprintf(F, __VA_ARGS__);

#define PRINT_ERR(...) fprintf(stderr, "\e[1;31merror:\e[0m " __VA_ARGS__ );
#define PRINT_WARN(...) fprintf(stderr, "\e[1;35mwarning:\e[0m " __VA_ARGS__ );

typedef char styp[8];

static inline double intpow(double x, unsigned int n){
  switch(n){
    case 0:   return 1.0;
    case 1:   return x;
    default: {
               double y = x;
               for(unsigned int i=1; i<n; i++){
                 y *= x;
               }
               return y;
             }
  }
}

