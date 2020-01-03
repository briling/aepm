#include "q.h"
#include "task_q.h"

#define  THEORY_DEFAULT INIT

int main(int argc, char * argv[]){

  if(argc<3){
    return 1;
  }

  FILE * fo;
  if(argc<4 || (fo = fopen(argv[3], "w")) == NULL){
    fo = stdout;
  }

  FILE * fb;
  fb = fopen(argv[1], "r");
  if(!fb){
    fprintf(stderr, "\tbasis?\n");
    fclose(fo);
    return 1;
  }

  urelconst_t * urelconst = urelconst_read(fb);
#if 0
  urelconst_print(urelconst, fo);
#endif

  basis_type btype;
  void * bas = bas_read(fb, &btype);
  fclose(fb);
  if(!bas){
    fprintf(stderr, "\tbasis!\n");
    free(urelconst);
    fclose(fo);
    return 1;
  }
#if 0
  bas_print(bas, btype, ">", fo);
#endif

  FILE * fm;
  fm = fopen(argv[2], "r");
  if(!fm){
    free(urelconst);
    free(bas);
    fclose(fo);
    fprintf(stderr, "\tmol?\n");
    return 1;
  }
  mol * m  = mol_read(fm);
  if(!m){
    free(urelconst);
    free(bas);
    fclose(fm);
    fclose(fo);
    fprintf(stderr, "\tmol!\n");
    return 1;
  }
  mol_print_m(m, 0, ">", fo);

  taskstr task = task_q_init(THEORY_DEFAULT, urelconst, argv[1], argv[2]);
  rewind(fm);
  task_q_read (fm, &task);
  fclose(fm);
  theory_t theory = task_q_proc(&task);
  if(!urelconst && task.control.aaar){
    GOTOHELL;
  }
#if 1
  task_q_print(fo, &task);
#endif

  int M; atomo * ao = ao_fill(m, bas, btype, &M);
#if 0
  ao_print(ao, M, fo);
#endif

  int N  = elnumber(m);
  int Nu = m->mult-1;
  int Nb = (N-Nu)/2;
  int Na = N-Nb;
  if( (N-Nu)%2 ) {
    fprintf(stderr, "\tN = %d, mult = %d !\n", N, m->mult);
    return 1;
  }
  if( Na > M) {
    fprintf(stderr, "\tN = %d, M = %d !\n", Na, M);
    return 1;
  }
  fprintf(fo, "N = %d (%d + %d)\n", N, Na, Nb);
  fprintf(fo, "M = %d\n", M);
  fflush(fo);

  if(theory == NIL){
    goto home;
  }

  /********************************************************/

  double * boys_array = boys_fill();

  fprintf(fo, "\n<<< one-electron integrals\n");
  double * H = calloc(5*symsize(M)*sizeof(double), 1);
  double * S = H + symsize(M);
  double * Dxyz = S + symsize(M);

  atomo * u_rel = NULL;
  if(task.control.aaar){
    u_rel = u_rel_fill(m, urelconst);
  }
  int * al = oneint_fill(M, m, bas, btype, ao, S, H, Dxyz, task.control.finite_nuclei, u_rel, boys_array, fo);
  free(u_rel);

  double * X = malloc((symsize(M)+M*M)*sizeof(double));
  double * S_isqrt = X + M*M;
  s_invsqrt_canorth(M, S, S_isqrt, X, 1e-15, 20, NULL);

  if(theory==INIT){

    double * V = malloc((M*M+M+symsize(M))*sizeof(double));
    double * C = V + M;
    double * F = C + M*M;

    init_lb20_heff(M, F, H, m, ao, al, btype, bas, boys_array);
#if 0
    mx_nosym_print(M, F, stdout);
#endif
#if 0
    oneint_print(M, S, F);
#endif
    veccp(M*M, C, X);
    mx_BHBt_sym(M, F, C);
    jacobi(F, C, V, M, 1e-15, 20, NULL);
    eigensort(M, V, C);
    if(pvec_write(M, ao, V, C, V, C, task.control.vectors)){
      fprintf(fo, "\nwrite coefficients to '%s'\n\n", task.control.vectors);
    }
    else{
      fprintf(fo, "\nfailed to write coefficients to '%s'\n", task.control.vectors);
    }

    free(V);
  }

  free(al);
  free(H);
  free(X);
  free(boys_array);
home:
  free(m);
  free(bas);
  free(urelconst);
  free(ao);
  fclose(fo);
  return 0;
}

