/**********************************************************************
 *
 *
 * copyright (c) 2001-2004, Hao Wu and Gary A. Churchill, The Jackson Lab.
 * written Nov, 2001
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * Part of the R/maanova package
 *
 *
 **********************************************************************/
/*
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Random.h>
*/

#include "util.h"

/* Equal probability sampling; with-replacement case */

void SampleReplace(int k, int n, int *y)
{
    int i;
    GetRNGstate();
    for (i = 0; i < k; i++)
        y[i] = n * unif_rand() + 1;
    PutRNGstate();
}

/* Equal probability sampling; without-replacement case */

void SampleNoReplace(int k, int n, int *y, int *x)
{
    int i, j;
    GetRNGstate();
    for (i = 0; i < n; i++)
        x[i] = i;
    for (i = 0; i < k; i++) {
        j = n * unif_rand();
        y[i] = x[j] + 1;
        x[j] = x[--n];
    }
    PutRNGstate();
}


/* multiply two matrices */
void matmult(double *result, double *a, int nrowa,
	     int ncola, double *b, int ncolb)

{
  int i, j, k;

  for(i=0; i<nrowa; i++)
    for(j=0; j<ncolb; j++) {
      /* clear the content of result */
      result[i*ncolb+j] = 0;
      for(k=0; k<ncola; k++)
	result[i*ncolb+j] += a[i*ncola+k]*b[k*ncolb+j];
    }

}

/* calculate the mean of an array */
double mean(double *x, int n)
{
  int i;
  double sum, result;

  for(i=0, sum=0.0; i<n; i++)
    sum += x[i];

  result = sum/n;
  return(result);
}

/* calculate the median of an array */
double median(double *x,  int n)
{
  int half;
  double result;

  /* sort x */
  R_rsort(x, n);

  /* find half */
  if( (n/2.0) == (n/2) ) { /* n is even number */
    half = n/2;
    result = (x[half]+x[half-1])/2;
  }
  else { /* n is odd number */
    half = (n+1)/2;
    result = x[half];
  }

  return(result);

}

/* James-Stein shrinkage estimator
 * This is working for a vector only */
void JS(double *X, int nX, double meanlog, double varlog, double *result)
{

  int i;
  double m, dtmp;

  /* take log on X then minus meanlog */
  for(i=0; i<nX; i++)
    X[i] = log(X[i]) - meanlog;

  /* calculate mean of log(X) */
  m = mean(X, nX);
  /* calculate (X-m)*(X-m) */
  for(i=0, dtmp=0.0; i<nX; i++)
    dtmp += (X[i]-m)*(X[i]-m);

  dtmp = 1 - (nX-3)*varlog/dtmp;
  if(dtmp < 0) dtmp = 0.0;

  /* calculate result */
  for(i=0; i<nX; i++)
    result[i] = exp(m + dtmp*(X[i]-m));
}
