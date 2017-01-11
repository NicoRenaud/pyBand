
#ifndef _spec_
#define _spec_


//extern int zggev_(char *jobvl, char *jobvr, int *n, 
//				 double *a, int *lda, 
//				 double *b, int *ldb, 
//				 double *alphar, double *alphai, double *beta, 
//				 double *vl, int *ldvl,
//                  double *vr, int *ldvr, double *work, int *lwork,
//                  int *info);


extern int zggev_(char *jobvl, char *jobvr, int *n,
        complex double *a, int *lda,
        complex double *b, int *ldb,
        complex double *alpha, complex double *beta,
        complex double *vl, int *ldvl,
        complex double *vr, int *ldvr,
        complex double *work, int *lwork,
        double *rwork,
        int *info);

void spec_pencil(double *VECT_PRP, double *VAL_PRP, double *Hreal, double *Himag, double *Sreal, double *Simag, int n); 


#endif