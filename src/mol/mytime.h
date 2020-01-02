#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

static inline double myutime(){
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*1e-6;
}

static inline double myrealtime(){
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec*1e-9;
}

