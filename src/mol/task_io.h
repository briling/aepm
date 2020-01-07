
static inline void printikey(FILE * f, const char s[], int x){
  fprintf(f, " %s = %d\n", s, x) ;
}
static inline void printskey(FILE * f, const char s[], char * x){
  fprintf(f, " %s = %s\n", s, x) ;
}

#define printxkey(F,S,X) _Generic((X),  int:             printikey, \
                                        char:            printikey, \
                                        char*:           printskey, \
                                        filename:        printskey  ) (F,S,X)

#define printkey(F,W,X)          printxkey(F,#X,(task->W.X))


