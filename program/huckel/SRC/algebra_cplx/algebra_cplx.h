///////////////////////////////////////
//
//    Proto des fonctions algebra_cplx.c
//
///////////////////////////////////////

#ifndef _algebraCplx_h

#define _algebraCplx_h

void dbl2cplx(complex double *CPLX, double *DBL, int n);
void transMat(complex double *M_TR, complex double *M_OR, int size);
void transMatrect(complex double *M_TR, complex double *M_OR, int size1, int size2);
void invMatComp(complex double *M_INV, complex double *M_OR, int size);
void trace_cplx(complex double *tr, complex double *MAT, int size);
#endif
