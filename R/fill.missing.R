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
# This is the function to fill in missing data in an object of rawdata
#
######################################################################

fill.missing <-
  function(rawdata, method="knn", k=20, dist.method="euclidean")
{
  if (class(rawdata) != "rawdata")
    stop("You have to provide an object of class rawdata.")

  method <- match.arg(method)
  
  # calculate the distance matrix for all data
  cat("Calculating pairwise distances ...\n")
  d <- as.matrix(dist(rawdata$data, method=dist.method))

  # find index for genes with NA values
  idx.missing <- which(apply(rawdata$data, 1, function(x) any(is.na(x))))
  
  # loop thru genes with missing
  cat("Missing data imputation ...\n")
  
  for(i in idx.missing) {
    # this distance between this gene and all other genes
    dist.thisgene <- d[i,]
    # find the K nearest neighbour
    dist.sort <- sort(dist.thisgene, index.return=TRUE)
    idx.neighbour <- dist.sort$ix[2:(k+1)]
    
    # find the loci for missing data
    missing.loci <- which(is.na(rawdata$data[i,]))
    # replace missing by a weighted average
    # the weight will be based on distance
    wt <- 1/dist.sort$x[2:(k+1)]
    wt <- wt/sum(wt)
    if(length(missing.loci) > 1)
      rawdata$data[i,missing.loci] <-
        apply(rawdata$data[idx.neighbour,missing.loci], 2,
              function(x) weighted.mean(x, wt, na.rm=TRUE))
    else
      rawdata$data[i,missing.loci] <-
        weighted.mean(rawdata$data[idx.neighbour,missing.loci],
                      wt, na.rm=TRUE)
  }

  rawdata
}

