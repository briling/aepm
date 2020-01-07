#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

typedef char filename[2048];

typedef struct{
  filename basis;
  filename vectors;
  int finite_nuclei;
  int aaar;
  int print;
} ctrlstr;

typedef struct{
  ctrlstr   control;
} taskstr;

taskstr task_q_init  (char * fbname, char * fname);
void    task_q_print (FILE * f, taskstr * task);
void    task_q_proc  (taskstr * task, void * urelconst);

void change_suffix(char * newname, const char * name, const char * suf, size_t size);

static inline void str_toupper(char * s){
  while(*s){
    *s = toupper(*s);
    s++;
  }
  return;
}

