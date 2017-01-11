#include "../defMacro.h"



/////////////////////////////////////////////////////////
//  Transforme une matrice double  en matrice complex    //
/////////////////////////////////////////////////////////
void dbl2cplx(complex double  *CPLX, double  *DBL, int n)
{
    int i;
    for(i=0;i<n;i++)
    {
	CPLX[i] = (complex double ) DBL[i];
    }
}

//////////////////////////////////////////////////////////
//     Compute the transpose of a matrix		//
//////////////////////////////////////////////////////////
void transMat(complex double *M_TR, complex double *M_OR, int size)
{ 
  int i,j;
  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
      M_TR[i+j*size] = M_OR[i*size+j];
  }
}

void transMatrect(complex double *M_TR, complex double *M_OR, int size1, int size2)
{ 
  int i,j;
  for(i=0;i<size1;i++)
  {
    for(j=0;j<size2;j++)
      M_TR[i+j*size1] = M_OR[i*size2+j];
  }
}


//////////////////////////////////////////////////////////
//     Compute the inverse of a matrix			//
//////////////////////////////////////////////////////////
void invMatComp(complex double *M_INV, complex double *M_OR, int size)
{
  
    int i,j;
    int n=size;
    int lwork, info;
    complex double wkopt;
    complex double *work;
    int pivot[size];
    
    // copie de la matrice originale    
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
	M_INV[i*size+j] = M_OR[i*size+j];
    }
    
    // LU factorisation of the matrix
    zgetrf_(&n, &n, M_INV, &n, pivot, &info);
    
    // optimal work computation 
    lwork = -1;
    zgetri_(&n, M_INV, &n, pivot, &wkopt, &lwork, &info);
    
    
    // memory allocation 
    lwork = (int) creal(wkopt);
    work = (complex double *) calloc(lwork,sizeof(complex double));
    
    // invert the matrix   
    zgetri_(&n, M_INV, &n, pivot, work, &lwork, &info);
    
    free(work);
  
}


//////////////////////////////////////////////////////////
//     Compute the trace_cplx of a matrix		//
//////////////////////////////////////////////////////////
void trace_cplx(complex double *tr, complex double *MAT, int size)
{
  
 int i;
 *tr = 0;
 for(i=0;i<size;i++)
   *tr +=  MAT[i*size+i];

}
