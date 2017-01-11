
#include "../defMacro.h"

void export_gnulot_semilog(char *NAME, double *TE, double *ENERGY, int nT)
{ 
  FILE *f;
  int j;

  f = fopen(NAME,"w");
  if (!f) 
  {
      printf("couldn't read %s\n",NAME);
      exit(1);
  }
  
  for(j=0;j<nT;j++)
	  fprintf(f,"%lf\t%lf\n",ENERGY[j],log10(TE[j])); 

  
  fclose(f);
}

void export_gnulot(char *NAME, double *TE, double *ENERGY, int nT)
{ 
  FILE *f;
  int j;

  f = fopen(NAME,"w");
  if (!f) 
  {
      printf("couldn't read %s\n",NAME);
      exit(1);
  }
  
  for(j=0;j<nT;j++)
	  fprintf(f,"%lf\t%lf\n",ENERGY[j],(TE[j])); 

  
  fclose(f);
}

