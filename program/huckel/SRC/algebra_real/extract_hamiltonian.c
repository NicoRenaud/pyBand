#include "../defMacro.h"
#include "../manag_output/print_screen.h"

////////////////////////////////////////////////////////////////////////////
//  EXTRACT A MATRIX FROM A LARGER ONE
////////////////////////////////////////////////////////////////////////////

void extract(double *Hout, double *Hin, int *ind1, int n1, int *ind2, int n2, int nin,char *file_name)
{
  int ii,jj,i,j;
	FILE *fout;
	fout = fopen(file_name,"w");
	

  for(ii=0;ii<n1;ii++)
  {
    i = ind1[ii];
    for(jj=0;jj<n2;jj++)
    {
      j = ind2[jj];
			
      Hout[ii*n2+jj] = Hin[i*nin+j];
			
			if(abs(Hout[ii*n2+jj])<10)
        fprintf(fout," %+5.4lf   ",Hout[ii*n2+jj]);
      else
        fprintf(fout,"%+5.4lf   ",Hout[ii*n2+jj]);
			
    }
		fprintf(fout,"\n");
  }
}



