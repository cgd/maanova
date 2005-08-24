######################################################################
#
# subset.madata.R
#
# copyright (c) 2001-2002, Hao Wu and Gary A. Churchill, The Jackson Lab.
#
# written Nov, 2001
# Modified Dec, 2002
# Modified Apr, 2004 for N-dye systems
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
# This is the function to subset an object of "madata"
#
######################################################################

subset.madata <-
  function(x, arrays, genes, ...)
{
  madata <- x
  if (class(madata) != "madata")
    stop("You have to provide an object of class madata.")

  # local variables
  narrays <- madata$n.array
  ngenes <- madata$n.gene
  ndyes <- madata$n.dye
  nreps <- madata$n.rep
  
  # create output variables
  result <- madata

  idx.array <- matrix(1:(narrays*ndyes), ndyes, narrays)
  idx.gene <- matrix(1:(ngenes*nreps), nreps, ngenes)
  
  if(!missing(arrays)) { # if arrays is given either in list or logical
    idx.array <- as.vector(idx.array[,arrays])
    result$n.array <- length(idx.array)/ndyes
  }
  else {
    arrays <- 1:narrays
  }
  
  if(!missing(genes)) { # if genes is given either in list or logical
    idx.gene <- as.vector(idx.gene[,genes])
    result$metarow <- madata$metarow[idx.gene]
    result$metacol <- madata$metacol[idx.gene]
    result$row <- madata$row[idx.gene]
    result$col <- madata$col[idx.gene]
    result$n.gene <- length(idx.gene)/nreps
  }
  else {
    genes <- 1:ngenes
  }

  # take out data
  result$data <- result$data[idx.gene,idx.array]
  # recalculate colmeans
  result$colmeans <- apply(result$data, 2, mean)
  # recalculate n.spot
  result$n.spot <- result$n.array * result$n.rep
  # take out flag if any
  if("flag" %in% names(madata))
    result$flag <- as.matrix(madata$flag[idx.gene, arrays])
  # take out ArrayName if any
#  if("ArrayName" %in% names(madata))
#    result$ArrayName <- madata$ArrayName[idx.array]
  # subset design (if any)
  if("design" %in% names(madata)) {
    tmp <- as.data.frame(madata$design[idx.array,])
    names(tmp) <- names(madata$design)
    result$design <- tmp
  }
  
  # take out the rest of fields
  fields.fixed <- c("n.dye","n.array","data","flag","design",
                    "n.gene","n.spot","n.rep","colmeans","offset","collapse",
                    "metarow", "metacol", "row", "col", "ArrayName",
                    "TransformMethod")
  fields.all <- names(madata)
  # find the user specified fields in madata
  fields.user <- setdiff(fields.all, fields.fixed)
  # Note that all these fields are supposed to be gene-based,
  # that is, they have dimension (n.gene*n.rep) x 1.
  # When averaging the replicates, these fields will be subset
  if(length(fields.user) >= 1) {
    for(i in 1:length(fields.user)) {
      idx <- match(fields.user[[i]], fields.all)
      result[[idx]] <- result[[idx]][genes]
    }
  }
  
  result

}
  
                                  
    
