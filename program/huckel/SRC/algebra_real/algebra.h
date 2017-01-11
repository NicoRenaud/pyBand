///////////////////////////////////////
//
//    Proto des fonctions algebra.c
//
///////////////////////////////////////

#ifndef _algebra_h

#define _algebra_h


void spec(double *VECT_PRP, double *VAL_PRP, double *MAT, int size);
void spec_pencil(double *VECT_PRP, double *VAL_PRP, double *H, double *S, int n) ;
void linspace(double *V, double min, double max, int n);
void transMat_real(double *M_TR, double *M_OR, int size);
int compare_doubles (const void *X, const void *Y);
void reorder_vect_prp(double *VECT_PRP,double *VAL_PRP,double *val_prp,int nb_orb);
void real_local2eigen(double *M, double *VP, int n);
void real_eigen2local(double *M, double *VP, int n);

#endif