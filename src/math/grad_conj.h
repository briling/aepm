{
  double MG = pars.MG;
  int    K  = pars.K;
  int N = GETN;                             /*=============*/
  double h_def = pars.o1dpars.h;

  double E, s1, s2, mg, w;

  double * grad = malloc(3*N * sizeof(double));
  if(grad == NULL){
    abort();
  }
  double * dir = grad + N;
  double * dir_prev = dir + N;

  WALK(&E, grad);                           /*=============*/

  int k = 0;
  do{
    vecset(N, dir_prev, 0.0);
    s1 = 1.0;
    for(int n=0; n<N; n++){
      mg = GETMG(N,grad);                   /*=============*/
      PRINTENT(f);                          /*=============*/
      PRINTEGK(f, E, mg, k);                /*=============*/
      fflush(f);
      if (mg < MG){
        fprintf(f, "converged (k  =%4d)\n\n", k);
        goto solovki;
      }
      if (++k > K){
        goto solovki;
      }
      s2 = vecdot(N, grad, grad);
      w  = s2/s1;
      vecsums(N, dir, grad, dir_prev, w);

      pars.o1dpars.h = OPT1D(&E, grad, dir, pars.o1dpars, f); /*======*/
      if(pars.o1dpars.h < 0){
        // OPT1D did not converge
        pars.o1dpars.h = h_def;
      }

      veccp(N, dir_prev, dir);
      s1 = s2;
    }
  } while(1);

  solovki:
    free(grad);
}

