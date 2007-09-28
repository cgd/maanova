/**********************************************************************
 *
 *
 * copyright (c) 2001-2007, Lei Wu and Gary A. Churchill, The Jackson Lab.
 * written Feb, 2007
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * Part of the R/maanova package
 *
 *
 **********************************************************************/
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
//#include <Defn.h>

#include "masvd.h"

SEXP masvd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v, SEXP method)
{
    int *xdims, n, p, lwork, info = 0;
    double *work, *xvals, tmp;
    SEXP val, nm;
    const char *meth;
    //char *meth;

/*    if (!(isString(jobu) && isString(jobv)))
	error(_("'jobu' and 'jobv' must be character strings"));
    if (!isString(method))
	error(_("'method' must be a character string"));*/
    meth = CHAR(STRING_ELT(method, 0));
/*#ifndef IEEE_754
    if (strcmp(meth, "dgesdd") == 0)
	error("method = \"dgesdd\" requires IEEE 754 arithmetic");
#endif*/
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0]; p = xdims[1];
    xvals = (double *) R_alloc(n * p, sizeof(double));
    /* work on a copy of x */
    Memcpy(xvals, REAL(x), (size_t) (n * p));

    if(strcmp(meth, "dgesdd")) {
	/* ask for optimal size of work array */
	lwork = -1;
	F77_CALL(dgesvd)(CHAR(STRING_ELT(jobu, 0)), CHAR(STRING_ELT(jobv, 0)),
			 &n, &p, xvals, &n, REAL(s),
			 REAL(u), INTEGER(getAttrib(u, R_DimSymbol)),
			 REAL(v), INTEGER(getAttrib(v, R_DimSymbol)),
			 &tmp, &lwork, &info);
//	if (info != 0)
//	    error(_("error code %d from Lapack routine '%s'"), info, "dgesvd");
	lwork = (int) tmp;

	work = (double *) R_alloc(lwork, sizeof(double));
	F77_CALL(dgesvd)(CHAR(STRING_ELT(jobu, 0)), CHAR(STRING_ELT(jobv, 0)),
			 &n, &p, xvals, &n, REAL(s),
			 REAL(u), INTEGER(getAttrib(u, R_DimSymbol)),
			 REAL(v), INTEGER(getAttrib(v, R_DimSymbol)),
			 work, &lwork, &info);
//	if (info != 0)
//	    error(_("error code %d from Lapack routine '%s'"), info, "dgesvd");
    } else {
	int ldu = INTEGER(getAttrib(u, R_DimSymbol))[0],
	    ldvt = INTEGER(getAttrib(v, R_DimSymbol))[0];
	int *iwork= (int *) R_alloc(8*(n<p ? n : p), sizeof(int));

	/* ask for optimal size of work array */
	lwork = -1;
	F77_CALL(dgesdd)(CHAR(STRING_ELT(jobu, 0)),
			 &n, &p, xvals, &n, REAL(s),
			 REAL(u), &ldu,
			 REAL(v), &ldvt,
			 &tmp, &lwork, iwork, &info);
//	if (info != 0)
//	    error(_("error code %d from Lapack routine '%s'"), info, "dgesdd");
	lwork = (int) tmp;
	work = (double *) R_alloc(lwork, sizeof(double));
	F77_CALL(dgesdd)(CHAR(STRING_ELT(jobu, 0)),
			 &n, &p, xvals, &n, REAL(s),
			 REAL(u), &ldu,
			 REAL(v), &ldvt,
			 work, &lwork, iwork, &info);
//	if (info != 0)
//	    error(_("error code %d from Lapack routine '%s'"), info, "dgesdd");
    }

    val = PROTECT(allocVector(VECSXP, 3));
    nm = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nm, 0, mkChar("d"));
    SET_STRING_ELT(nm, 1, mkChar("u"));
    SET_STRING_ELT(nm, 2, mkChar("vt"));
    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, s);
    SET_VECTOR_ELT(val, 1, u);
    SET_VECTOR_ELT(val, 2, v);
    UNPROTECT(2);
    return val;
}


