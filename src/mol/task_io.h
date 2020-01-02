
/* input */

static inline void printikey(FILE * f, const char s[], int x){
  fprintf(f, " %s = %d\n", s, x) ;
}
static inline void printukey(FILE * f, const char s[], unsigned int x){
  fprintf(f, " %s = %u\n", s, x) ;
}
static inline void printfkey(FILE * f, const char s[], double x){
  fprintf(f, " %s = %20.15e\n", s, x) ;
}
static inline void printskey(FILE * f, const char s[], char * x){
  fprintf(f, " %s = %s\n", s, x) ;
}

#define printxkey(F,S,X) _Generic((X),  int:             printikey, \
                                        unsigned int:    printukey, \
                                        char:            printikey, \
                                        unsigned char:   printukey, \
                                        char*:           printskey, \
                                        filename:        printskey, \
                                        double:          printfkey  ) (F,S,X)

#define printkey_old(F,T,W,X)    fprintf(F, " " #X " = " T "\n", (task->W.X))
#define printkey(F,W,X)          printxkey(F,#X,(task->W.X))
#define printtask(F,X)           if(task->X ) {fprintf(F, " " #X " = %d\n", task->X);}


/* output */

static inline void setikey(int    * key, char * val){
  *key = atoi(val);
}
static inline void setckey(char   * key, char * val){
  *key = atoi(val);
}
static inline void setfkey(double * key, char * val){
  *key = atof(val);
}
static inline void setskey(filename * key, char * val){
  filename tmp;
  if(sscanf(val, "%255s", tmp)){
    strncpy((char *)key, tmp, sizeof(*key)-1);
  }
}

#define setxkey(KEY,VAL) _Generic((KEY),  int:             setikey, \
                                          unsigned  int:   setikey, \
                                          unsigned char:   setckey, \
                                          char:            setckey, \
                                          char*:           setskey, \
                                          double:          setfkey  ) (&KEY,VAL)

#define iskey_old(key,val,T,W,KEY)  if(!strcmp(key, #KEY)){ task->W.KEY = ato##T(val);   }
#define iskey(key,val,W,KEY)        if(!strcmp(key, #KEY)){ setxkey((task->W.KEY), val); }
#define istask(key,val,KEY)         if(!strcmp(key, #KEY)){ task->KEY = atoi(val);       }

#define tokens(key,val,str)                  \
  if(strchr(str, '#'))            break;     \
  if(!(key = strtok(str," =")))   continue;  \
  if(!(val = strtok(NULL," =")))  continue;
