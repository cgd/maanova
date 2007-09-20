######################################################################
#
# dyeswapfilter.R
#
# copyright (c) 2001-2003, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written May, 2003
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
#
######################################################################

dyeswapfilter <- function(dataobj, r=4)
{

  if(dataobj$n.dye != 2)
    stop("dyeswapfilter only works for 2-dye arrays")

  # calculate the logratios without normalization
  lr <- make.ratio(dataobj, FALSE)
  
  # initialize flag
  flag <- dataobj$flag
  if(is.null(flag))
    flag <- matrix(0, dim(lr)[1], dataobj$n.array)
  
  # loop thru all arrays
  for (i in 1:(dataobj$n.array/2)) {
    # different in logratios for two dye-swapped array
    d <- lr[,2*i-1] + lr[,2*i]
    # IQR of d
    iqrd <- IQR(d)
    sdd <- iqrd / 1.35
    # flag for this pair
    idx <- which(abs(d)>r*sdd)
    flag[idx,2*i-1] <- 1
    flag[idx,2*i] <- 1
  }

  # result
  result <- dataobj
  result$flag <- flag
  result
}
 
