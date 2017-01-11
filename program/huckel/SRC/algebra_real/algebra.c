#include "../defMacro.h"



/////////////////////////////////////////////////////////
//               Diagonalisation                       //
/////////////////////////////////////////////////////////

void spec(double *VECT_PRP, double *VAL_PRP, double *MAT, int n) {


	double  *work=NULL;
        int lwork,info,i;

	// allocation memoire
	work = (double *)malloc(1*sizeof(double));

	// copie de la matrice dans a forme vectorielle
	for(i=0;i<n*n;i++)
		VECT_PRP[i] = MAT[i];


	// precalcul du lwork optimal
	lwork=-1;
        dsyev_( "N", "U", &n, VECT_PRP, &n, VAL_PRP, work, &lwork, &info );
	lwork= work[0];

	// realocation de work
        free(work);
        work = (double *)malloc(lwork*sizeof(double));

	//calcul des vecteurs/valeurs propres;
         dsyev_( "V", "U", &n, VECT_PRP, &n, VAL_PRP, work, &lwork, &info );
	 
	 // free memory
	 free(work);

}

/////////////////////////////////////////////////////////
//               Diagonalisation                       //
/////////////////////////////////////////////////////////

void spec_pencil(double *VECT_PRP, double *VAL_PRP, double *H, double *S, int n) 
{
  int i;
  int nn = n, one = 1, lwork, info;
  double tmp, *work; 
  double *VR, *VL, alphar[n], alphai[n], beta[n]; 

  // memory alloc
  VR = calloc(n*n,sizeof(double));
  VL = calloc(n*n,sizeof(double));
  
  
  // precomppute the optimal lwork
  lwork = -1;
  dggev_("N", "N", &nn, H, &nn, S, &nn, 
	  alphar, alphai, beta, VL, &one, VR, &one, 
	  &tmp, &lwork, &info); 
  
  //real computation
  lwork = (int) tmp; 
  work = (double *) malloc(sizeof(double)*lwork); 

  
  dggev_("V", "V", &nn, H, &nn, S, &nn, 
	alphar, alphai, beta, VL, &nn, VR, &nn, 
	work, &lwork, &info); 

  if(info == 3)
    printf("warning : the lapack routine dggev_ failed \n"); 
  
  // store the eigenvalues
  for(i=0;i<n;i++)
   VAL_PRP[i] = (alphar[i])/beta[i];
  
  //store the eigenvalues
  cblas_dcopy(n*n,VR,1,VECT_PRP,1);
  
  // free memory
   free(work);
   free(VL);
   free(VR);
 
}


///////////////////////////////////////////////////////////////////
//  Transforme une matrice de la base  locale dans la base propre //
///////////////////////////////////////////////////////////////////
void real_local2eigen(double *M, double *VP, int n)
{
  double *M_temp, *VP_temp;
  
  double a = 1.0;
  double b = 0.0;
  
  // temporary matrix for multiplication
  M_temp = calloc(n*n,sizeof(double));
  VP_temp = calloc(n*n,sizeof(double));
  
  // copy the eigenvector matrix
  cblas_dcopy(n*n,VP,1,VP_temp,1);
	
  // multiplication : M_temp = VP*M
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,
              a,VP_temp,n,M,n,b,M_temp,n);
  
  // deuxieme multiplication: M = M_temp*VP'
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,n,n,n,
              a,M_temp,n,VP_temp,n,b,M,n);
  
  free(M_temp);
  free(VP_temp);
}




///////////////////////////////////////////////////////////////////
//  Transforme une matrice de la base propre dans la base locale //
///////////////////////////////////////////////////////////////////
void real_eigen2local(double *M, double *VP, int n)
{
  
  double *M_temp, *VP_temp;
  
  double a = 1.0;
  double b = 0.0;
  
  // temporary matrix for multiplication
  M_temp = calloc(n*n,sizeof(double));
  VP_temp = calloc(n*n,sizeof(double));
  
  // copy the eigenvector matrix
  cblas_dcopy(n*n,VP,1,VP_temp,1);
	
  
  // multiplication : M_temp = VP*M
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,n,n,n,
              a,VP_temp,n,M,n,b,M_temp,n);
  

  
  // deuxieme multiplication: M = M_temp*VP'
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,
              a,M_temp,n,VP_temp,n,b,M,n);
  
  free(M_temp);
  free(VP_temp);
}


/////////////////////////////////////////////////////////
//               create linspace vector                //
/////////////////////////////////////////////////////////
void linspace(double *V, double min, double max, int n)
{
  
 int i;
 double dV = (max-min)/(n-1);
 for(i=0;i<n;i++)
   V[i] = min+i*dV;
  
}


////////////////////////////////////////////////////////
//              transpose a real matrix               //
////////////////////////////////////////////////////////

void transMat_real(double *M_TR, double *M_OR, int size)
{ 
  int i,j;
  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
      M_TR[i+j*size] = M_OR[i*size+j];
  }
}


////////////////////////////////////////////////////////
//             Sorting function qsort                 //
////////////////////////////////////////////////////////
int compare_doubles (const void *X, const void *Y)
{
       double x = *((double *)X);
       double y = *((double *)Y);

       if (x > y)
       {
               return 1;
       }
       else
       {
               if (x < y)
               {
                       return -1;
               }
               else
               {
                       return 0;
               }
       }
}


////////////////////////////////////////////////////////
//            Reorder Eigenvector                     //
////////////////////////////////////////////////////////
void reorder_vect_prp(double *VECT_PRP,double *VAL_PRP,double *val_prp,int nb_orb)
{
   int i,j,k;
   double *vect_prp_cpy;
   vect_prp_cpy = calloc(nb_orb*nb_orb,sizeof(double));
   cblas_dcopy(nb_orb*nb_orb,VECT_PRP,1,vect_prp_cpy,1);
   
   
   // for all the eigenavalues
   for(i=0;i<nb_orb;i++)
   {
     // we find the good one
     for(j=0;j<nb_orb;j++)
     {
	if(VAL_PRP[i]==val_prp[j])
	{
	  for(k=0;k<nb_orb;k++)
	    VECT_PRP[i*nb_orb+k] = vect_prp_cpy[j*nb_orb+k];
	}
       
     }
  }
  
  free(vect_prp_cpy);
}