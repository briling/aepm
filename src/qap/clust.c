#include "q.h"
#include "qap.h"

void mols_clusterize(int npar, int nmol,
    mol_data * md, char ** fnames, FILE * fo){

  int cl[NATOMS] = {};
  int c = 1;
  for(int j=0; j<npar; j++){
    if(cl[j]==0){
      cl[j] = c++;
    }
    for(int i=0; i<nmol; i++){
      if(md[i].mpl[j]){
        for(int k=0; k<npar; k++){
          if(cl[k] && md[i].mpl[k]){
            cl[j] = cl[k];
          }
        }
        for(int k=0; k<npar; k++){
          if(!cl[k] && md[i].mpl[k]){
            cl[k] = cl[j];
          }
        }

      }
    }
  }

  for(int k=1; k<c; k++){
    for(int i=0; i<nmol; i++){
      for(int j=0; j<npar; j++){
        if(cl[j]==k && md[i].mpl[j]){
          fprintf(fo, "%3d : %s\n", k, fnames[i]);
          break;
        }
      }
    }
    fprintf(fo, "\n");
  }
  return;
}
