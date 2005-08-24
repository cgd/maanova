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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Random.h>

void shift(double *r, double *g, int *nrows, double *c, 
	   double *sad, double *offset)
{
  /* local variables */
  int i, j, k;
  int idx[100];
  double *r_shift, *g_shift, *x, med;

  x = (double *)R_alloc(3*(*nrows), sizeof(double));
  r_shift = x+(*nrows);
  g_shift = r_shift+(*nrows);

  for(i=0, sad[i]=0.0; i<100; i++) {
    idx[i] = i;

    for(j=0; j<*nrows; j++) {
      r_shift[j] = r[j] - c[i];
      if(r_shift[j]<1) r_shift[j] = 1;
      g_shift[j] = g[j] - c[i];
      if(g_shift[j]<1) g_shift[j] = 1;
      x[j] = log(r_shift[j]/g_shift[j]);
    }
    /* calculate the median of x */
    med = median(x, *nrows);
    for(j=0; j<*nrows; j++)
      sad[i] += abs(x[j]-med);
  }

  /* sort sad with index */
  rsort_with_index(sad, idx, 100);
 
  *offset = c[idx[0]];

}
      
