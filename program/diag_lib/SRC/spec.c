#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>
#include <complex.h>
#include <cblas.h>
//#include <clapack.h>
#include "spec.h"

/////////////////////////////////////////////////////////
//               Diagonalisation                       //
/////////////////////////////////////////////////////////
void spec_pencil(double *VECT_PRP, double *VAL_PRP, double *Hreal, double *Himag, double *Sreal, double *Simag, int n) 
{
  int i,j;
  int nn = n, one = 1, lwork, info;
  double rwork;
  complex double tmp, *work; 
  complex double *VR, *VL, *H, *S, alpha[n], beta[n]; 

  // memory alloc
  H = calloc(n*n,sizeof(complex double));
  S = calloc(n*n,sizeof(complex double));

  VR = calloc(n*n,sizeof(complex double));
  VL = calloc(n*n,sizeof(complex double));

  // form the complex matrix
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      H[i*n+j] = Hreal[i*n+j] + I*Himag[i*n+j];
      S[i*n+j] = Sreal[i*n+j] + I*Simag[i*n+j];
    }
  }
  
  
  // precomppute the optimal lwork
  lwork = -1;
  zggev_("N", "N", &nn, H, &nn, S, &nn, 
	  alpha, beta, VL, &one, VR, &one, 
	  &tmp, &lwork, &rwork, &info); 
  
  //real computation
  lwork = (int) tmp; 
  work = (complex double *) malloc(sizeof(complex double)*lwork); 

  
  zggev_("V", "V", &nn, H, &nn, S, &nn, 
	alpha, beta, VL, &nn, VR, &nn, 
	work, &lwork, &rwork, &info); 

  if(info == 3)
    printf("warning : the lapack routine zggev_ failed \n"); 
  
  // store the eigenvalues and eigenvector in doubles
  for(i=0;i<n;i++)
    VAL_PRP[i] = 0; //(double) creal((alpha[i])/beta[i]);
    for (j=0;j<n;j++)
      VECT_PRP[i*n+j] = 0; //(double) creal(VR[i*n+j]);
  
  
  // free memory
  free(H);
  free(S);
  free(work);
  free(VL);
  free(VR);
 
}
