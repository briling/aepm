#include "q.h"
#include "task_q.h"
#include "task_io.h"
#include <ctype.h>

void task_q_print(FILE * f, taskstr * task){
  fprintf(f, "#control\n");
  printkey(f,  control, theory);
  printkey(f,  control, vectors);
  if(task->control.aaar){
    printkey(f,  control, aaar);
  }
  if(task->control.aaar || task->control.finite_nuclei){
    printkey(f,  control, finite_nuclei);
  }
  printkey(f,  control, basis);
  fprintf(f, "#\n");
}

static void tread(FILE * f, taskstr * task, char * s){
  filename str;
  char * key;
  char * val;
  if (! strcmp(s, "#control")){
    while(fgets(str, sizeof(str), f)) {
      tokens(key, val, str);
      iskey(key, val, control, theory);
      iskey(key, val, control, finite_nuclei);
      iskey(key, val, control, aaar);
      iskey(key, val, control, vectors);
    }
  }
  return;
}

void task_q_read(FILE * f, taskstr * task){

  filename s;

  while(fgets(s, sizeof(s), f)) {
    if(s[0]=='\n'){
      continue;
    }
    s[strnlen(s, sizeof(s))-1] = '\0';
    char * p = strchr(s, '#');
    if(p && p[1]){
      tread(f, task, p);
    }
  }
  return;
}

taskstr task_q_init(theory_t theory, void * urelconst, char * fbname, char * fname){

  taskstr task = {};

  strcpy(task.control.theory, theory_s[theory]);
  task.control.finite_nuclei = -1;
  task.control.aaar = !!urelconst;
  strcpy(task.control.basis, fbname);

  change_suffix(task.control.vectors, fname, ".vec", sizeof(task.control.vectors));

  return task;
}

static inline void str_toupper(char * s){
  while(*s){
    *s = toupper(*s);
    s++;
  }
  return;
}

theory_t task_q_proc(taskstr * task){
  str_toupper(task->control.theory);
  int nth = sizeof(theory_s)/sizeof(theory_s[0]);
  int f   = 0;
  for(int i=0; i<nth; i++){
    if(!strcmp(task->control.theory, theory_s[i])){
      f = i;
      break;
    }
  }
  if(task->control.finite_nuclei<0){ // unspecified
    task->control.finite_nuclei = !!(task->control.aaar);
  }
  return f;
}

static inline void remove_in(char * fname){
  int t = strlen(fname) - 1;
  if(fname[t]=='n' && t-1>=0 && fname[t-1]=='i' && t-2>=0 && fname[t-2]=='.'){
    fname[t] = fname[t-1] = fname[t-2] = '\0';
  }
  return;
}

static inline void add_suffix(char * fname, const char * suf, size_t size){
  if(strlen(fname) + strlen(suf) >= size){
    GOTOHELL;
  }
  strcat(fname, suf);
  return;
}

void change_suffix(char * newname, const char * name, const char * suf, size_t size){
  strncpy(newname, name, size-1);
  remove_in(newname);
  add_suffix(newname, suf, size);
  return;
}

