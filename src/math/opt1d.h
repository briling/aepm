{
  double h = pars.h;
  double D = pars.D;
  double H = pars.H;
  double G = pars.G;
  int    K = pars.K;

  int N = GETN;                             /*=============*/
  double * l = malloc(N*sizeof(double));

  veccp (N, l, dir);
  double s = vecdot(N, l, l);
  vecscal(N, l, 1.0/sqrt(s));
  double g1 = vecdot(N, gr, l);
  int it = 1;
  while (fabs(g1) > G){

    double d = -g1/h;
    if (fabs(d) > D){
      d = d>=0.0 ? D : -D;
    }

    fprintf(f, "  it:%3d   g1:%+.6e   h:%+.6e   d:%+.6e   e:%+.6e\n", it, g1, h, d, *energy);
    fflush(f);
    ADDS(N,l,d);                            /*=============*/
    double g2 = g1;
    WALK(energy, gr);                           /*=============*/
    g1 = vecdot(N, gr, l);
    h  = (g1-g2)/d;
    if (h<H){
      h = H;
    }
    if (it++ == K){
      h = -1.0;
      fprintf(f, "  no convergence\n");
      break;
    }
  }
  free(l);
  return h;
}

