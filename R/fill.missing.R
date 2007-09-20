######################################################################
#
# fill.missing.R
#
# copyright (c) 2001-2005, Hao Wu and Gary A. Churchill, The Jackson Lab.
#
# written July 2005
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# This is the function to fill in missing data in an object of madata
#
######################################################################
fill.missing <-
  function(data, method="knn", k=20, dist.method="euclidean")
{
  if (class(data) != "madata")
    stop("You have to provide an object of class madata.")

  method <- match.arg(method)
  
  # calculate the distance matrix for all data
  cat("Calculating pairwise distances ...\n")
  d <- as.matrix(dist(data$data, method=dist.method))

  # find index for genes with NA values
  idx.missing <- which(apply(data$data, 1, function(x) any(is.na(x))))
  
  # loop thru genes with missing
  cat("Missing data imputation ...\n")
  
  for(i in idx.missing) {
    # this distance between this gene and all other genes
    dist.thisgene <- d[i,]
    # find the K nearest neighbour
    dist.sort <- sort(dist.thisgene, index.return=TRUE)
    idx.neighbour <- dist.sort$ix[2:(k+1)]
    
    # find the loci for missing data
    missing.loci <- which(is.na(data$data[i,]))
    # replace missing by a weighted average
    # the weight will be based on distance
    wt <- 1/dist.sort$x[2:(k+1)]
    wt <- wt/sum(wt)
    if(length(missing.loci) > 1)
      data$data[i,missing.loci] <-
        apply(data$data[idx.neighbour,missing.loci], 2,
              function(x) weighted.mean(x, wt, na.rm=TRUE))
    else
      data$data[i,missing.loci] <-
        weighted.mean(data$data[idx.neighbour,missing.loci],
                      wt, na.rm=TRUE)
  }
  data
}
