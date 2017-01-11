
#include "../defMacro.h"



#include "../defMacro.h"

////////////////////////////////////
// Print the header of the program //
////////////////////////////////////



void print_header(SYS_XYZ sys) {
 
  int i;
  printf("\n====  ===============================  ====\n");    
  printf("====                                   ====\n");   
  printf("====  	        Husky.1	               ====\n"); 
  printf("====                                   ====\n");   
  printf("====  ===============================  ====\n\n");
  
    printf(" == Coordinates of the molecule    \t: %s\n",sys.pos);
  printf(" == Index of MO for frag1    \t");
  for(i=0;i<sys.nbr_mo_1;i++)
    printf("%d \t",sys.index_mo_1[i]);
  printf("\n");
  printf(" == Index of MO for frag2    \t");
  for(i=0;i<sys.nbr_mo_2;i++)
    printf("%d \t",sys.index_mo_2[i]);
  printf("\n");
  printf(" == Index of MO for frag3    \t");
  for(i=0;i<sys.nbr_mo_3;i++)
    printf("%d \t",sys.index_mo_3[i]);
  printf("\n");
  printf(" \n====  ===============================  ====\n\n");
  
}









////////////////////////////////////
//      Print a Matrix            //
////////////////////////////////////

void print_mat(double * MAT,int n1,int n2, char *name){
      
     int i,j,k=0;
     printf("\n %s= [\n",name);
     for(i=0;i<n1;i++)
     {
//                       printf("\t\t");
                      for(j=0;j<n2;j++)
                      {
                                       printf("%e ",MAT[k]);
					k++;
                      }
                      printf("\n");
      }  
      printf("\t     ];\n"); 
}

////////////////////////////////////
//      Print a ComplexMatrix       //
////////////////////////////////////
void print_mat_C(complex double  *MAT, int n1, int n2, char *name)
{
      
     int i,j,k=0;
     printf("\n %s= [\n",name);
     for(i=0;i<n1;i++)
     {
                      printf("\t");
                      for(j=0;j<n2;j++)
                      {
                                       printf("%1.6f+%%i*%1.6f ",crealf(MAT[k]),cimagf(MAT[k]));
					k++;
                      }
                      printf("\n");
      }  
      printf("\t ];\n"); 
}


////////////////////////////////////
//      Print a Vector            //
///////////////////////////////////
void print_vect(double  *VECT, int n,char *name){
     int i;
     printf("\n %s= [ ",name);
     for(i=0;i<n;i++)
                     printf("%2.5g ",VECT[i]);
     printf(" ]\n");
}


void print_vect_real(int  *VECT, int n,char *name){
     int i;
     printf("\n %s= [ ",name);
     for(i=0;i<n;i++)
                     printf("%d ",VECT[i]);
     printf(" ]\n");
}

////////////////////////////////////
//      Print a Vector            //
///////////////////////////////////
void print_vect_C(complex double  *VECT, int n,char *name){
     int i;
     printf("\n %s= [ ",name);
     for(i=0;i<n;i++)
                     printf("%g +i %g ",creal(VECT[i]),cimag(VECT[i]));
     printf(" ]\n");
}
