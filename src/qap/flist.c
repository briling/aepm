#include "q.h"
#include "qap.h"

char ** flist_read(int * nlines, FILE * fl){

  long flstart = ftell(fl);
  long flsize  = 0;
  while(fgetc(fl)!=EOF) flsize++;
  fseek(fl, flstart, SEEK_SET);
  char * flcont = malloc(flsize);
  fread(flcont, flsize, 1, fl);

  int n = 0;
  if(flcont[flsize-1] != '\n'){
    PRINT_ERR("\tlist?\n");
    GOTOHELL;
  }
  for(int i=0; i<flsize; i++){
    if( (flcont[i] == '#') && (!i || flcont[i-1] == '\n') ){
      n--;
    }
    else if(flcont[i] == '\n'){
      n++;
    }
  }
  if(!n){
    PRINT_ERR("\tlist?\n");
    GOTOHELL;
  }

  char ** fnames = realloc(flcont, flsize + sizeof(char *)*n);
  if(!fnames) GOTOHELL;
  flcont = (char *)(fnames + n);
  memmove(flcont, fnames, flsize);

  int t = 0;
  int comflag = 0;
  for(int i=0; i<flsize; i++){
    if(flcont[i] == '#'){
      comflag = 1;
    }
    if(flcont[i] == '\n'){
      comflag = 0;
    }
    if( (!comflag) && ( (!i) || (flcont[i-1] == '\n') ) ){
      fnames[t++] = flcont + i;
    }
    if(i && flcont[i-1] == '\n'){
      flcont[i-1] = '\0';
    }
    if(flcont[i] == '#'){
      flcont[i] = '\0';
    }
  }
  flcont[flsize-1] = '\0';

  for(int i=0; i<n; i++){
    int l = strlen(fnames[i]);
    for(int j=l-1; j>=0; j--){
      if(fnames[i][j] == ' '){
        fnames[i][j] = '\0';
      }
      else{
        break;
      }
    }
  }
  *nlines = n;
  return fnames;
}
