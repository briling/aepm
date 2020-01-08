#include "q.h"
#include "task_q.h"

int main(int argc, char * argv[]){

  if(argc<3){
    return 1;
  }

  FILE * fb;
  fb = fopen(argv[1], "r");
  if(!fb){
    PRINT_ERR("\tbasis?\n");
    return 1;
  }

  urelconst_t * urelconst = urelconst_read(fb);
  basis_type btype;
  void * bas = bas_read(fb, &btype);
  fclose(fb);
  if(!bas){
    free(urelconst);
    PRINT_ERR("\tbasis!\n");
    return 1;
  }


  FILE * fm;
  fm = fopen(argv[2], "r");
  if(!fm){
    free(urelconst);
    free(bas);
    PRINT_ERR("\tmolecule?\n");
    return 1;
  }

  mol * m  = mol_read(fm);
  fclose(fm);
  if(!m){
    free(urelconst);
    free(bas);
    PRINT_ERR("\tmolecule!\n");
    return 1;
  }


  FILE * fo = stdout;
  taskstr task = task_q_init(argv[1], argv[2]);
  for(int i=3; i<argc; i++){
    if( sscanf (argv[i], "vectors:%d",        task.control.vectors      )) { continue; }
    if( sscanf (argv[i], "aaar:%d",          &task.control.aaar         )) { continue; }
    if( sscanf (argv[i], "finite_nuclei:%d", &task.control.finite_nuclei)) { continue; }
    if( sscanf (argv[i], "print:%d",         &task.control.print        )) { continue; }
    if(! (fo = fopen(argv[i], "w"))){
      fo = stdout;
    }
  }
  task_q_proc(&task, urelconst);


#if 0
  urelconst_print(urelconst, fo);
#endif
#if 0
  bas_print(bas, btype, ">", fo);
#endif
#if 1
  mol_print_m(m, 0, ">", fo);
#endif
#if 1
  task_q_print(fo, &task);
#endif

  /********************************************************/

  int M; atomo * ao = ao_fill(m, bas, btype, &M);
#if 0
  ao_print(ao, M, fo);
#endif

  int N  = elnumber(m);
  int Nu = m->mult-1;
  int Nb = (N-Nu)/2;
  int Na = N-Nb;
  if( (N-Nu)%2 ) {
    PRINT_ERR("\tN = %d, mult = %d !\n", N, m->mult);
    return 1;
  }
  if( Na > M) {
    PRINT_ERR("\tN = %d, M = %d !\n", Na, M);
    return 1;
  }
  fprintf(fo, "N = %d (%d + %d)\n", N, Na, Nb);
  fprintf(fo, "M = %d\n", M);
  fflush(fo);

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
  int * al = oneint_fill(M, m, bas, btype, ao, S, H, Dxyz, task.control.finite_nuclei, u_rel, boys_array, NULL);
  free(u_rel);

  double * X = malloc((symsize(M)+M*M)*sizeof(double));
  double * S_isqrt = X + M*M;
  s_invsqrt_canorth(M, S, S_isqrt, X, 1e-15, 20, NULL);

  /********************************************************/

  double * V = malloc((M*M+M+symsize(M))*sizeof(double));
  double * C = V + M;
  double * F = C + M*M;
  fprintf(fo, "\n<<< two-electron integrals\n");
  init_lb20_heff(M, F, H, m, ao, al, btype, bas, boys_array);
#if 0
  oneint_print(M, S, F);
#endif
  veccp(M*M, C, X);
  mx_BHBt_sym(M, F, C);
  jacobi(F, C, V, M, 1e-15, 20, NULL);
  eigensort(M, V, C);

  double * C1 = malloc(M*M*sizeof(double));
  vec_to_p(M, ao, C1, C);
  if(qvec_write(M, V, C1, V, C1, task.control.vectors)){
    fprintf(fo, "\nwrite coefficients to '%s'\n\n", task.control.vectors);
  }
  else{
    fprintf(fo, "\nfailed to write coefficients to '%s'\n", task.control.vectors);
  }

  if(task.control.print){
    double * D1 = malloc(symsize(M)*sizeof(double));
    D_fill(N, M, C1, D1, 2);
    switch(task.control.print){
      case 1: mx_nosym_print(M, D1, fo); break;
      case 2: mx_sym_print  (M, D1, fo); break;
      case 3: mx_print      (M, C1, fo); break;
      case 4: molden_print(N/2, M, 2, m, ao, bas, btype, C, V, "Alpha", "mos>", fo); break;
    }
    free(D1);
  }
  free(C1);

  free(V);

  /********************************************************/

  free(al);
  free(H);
  free(X);
  free(boys_array);

  free(m);
  free(bas);
  free(urelconst);
  free(ao);
  fclose(fo);
  return 0;
}

