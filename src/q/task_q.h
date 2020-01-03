#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

typedef enum {NIL, INIT} theory_t;

typedef char filename[256];

static const filename theory_s[] = {
  [INIT]    = "INIT",
};

typedef struct{
  filename theory;
  filename basis;
  filename vectors;
  char finite_nuclei;
  unsigned char aaar;
} ctrlstr;

typedef struct{
  ctrlstr   control;
} taskstr;

taskstr  task_q_init  (theory_t theory, void * urelconst, char * fbname, char * fname);
void     task_q_read  (FILE * f, taskstr * task);
void     task_q_print (FILE * f, taskstr * task);
theory_t task_q_proc  (taskstr * task);

void change_suffix(char * newname, const char * name, const char * suf, size_t size);

static inline void str_toupper(char * s){
  while(*s){
    *s = toupper(*s);
    s++;
  }
  return;
}

