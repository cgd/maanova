#include <R_ext/Lapack.h>
#include <Rdefines.h> /* Rinternals.h + GET_SLOT etc */
#include <R.h>  /* includes Rconfig.h */
#include <R_ext/Complex.h>
#include <R_ext/RS.h>



SEXP masvd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v, SEXP method);

/* DGESVD - compute the singular value decomposition (SVD); of a   */
/* real M-by-N matrix A, optionally computing the left and/or      */
/* right singular vectors                                          */
void F77_NAME(dgesvd)(const char *jobu, const char *jobvt,
		      const int *m, const int *n,
		      double *a, const int *lda, double *s,
		      double *u, const int *ldu,
		      double *vt, const int *ldvt,
		      double *work, const int *lwork, int *info);
