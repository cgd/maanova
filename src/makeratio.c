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

void makeratio(double *adjdata, double *colmeans, double *std, 
	       int *norm_std, int *nrow, int *ncol, double *result)
{

  /* local variables */
  int i, j;
  double *tmp;

  /* allocate memory for tmp */
  tmp = (double *)R_alloc((*nrow)*(*ncol), sizeof(double));

  /* subtract column means from adjdata */
  for(i=0; i<*ncol; i++) {
    for(j=0; j<*nrow; j++) {
      tmp[i*(*nrow)+j] = adjdata[i*(*nrow)+j] - colmeans[i];
      if(*norm_std) /* if need to normalize the data by std */
	tmp[i*(*nrow)+j]  /= std[i];
    }
  }
   
  /* calculate the ratio */
  for(i=0; i<(*ncol)/2; i++)
    for(j=0; j<*nrow; j++) 
      result[i*(*nrow)+j] = (tmp[(2*i)*(*nrow)+j] - 
			     tmp[(2*i+1)*(*nrow)+j]);
  
}
