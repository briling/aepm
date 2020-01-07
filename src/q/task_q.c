#include "q.h"
#include "task_q.h"
#include "task_io.h"

void task_q_print(FILE * f, taskstr * task){
  fprintf(f, "#control\n");
  printkey(f,  control, vectors);
  printkey(f,  control, aaar);
  printkey(f,  control, finite_nuclei);
  printkey(f,  control, print);
  printkey(f,  control, basis);
  fprintf(f, "#\n");
  return;
}

taskstr task_q_init(char * fbname, char * fname){
  taskstr task = {};
  task.control.finite_nuclei = -1;
  task.control.aaar          = -1;
  strcpy(task.control.basis, fbname);
  change_suffix(task.control.vectors, fname, ".vec", sizeof(task.control.vectors));
  return task;
}

void task_q_proc(taskstr * task, void * urelconst){
  if(task->control.aaar==1 && !urelconst){
    PRINT_WARN("Cannot use aaar without parameters\n");
    task->control.aaar = 0;
  }
  else if(task->control.aaar<0){
    task->control.aaar = !!urelconst;
  }
  if(task->control.finite_nuclei<0){
    task->control.finite_nuclei = !!(task->control.aaar);
  }
  return;
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

