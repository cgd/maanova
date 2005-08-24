######################################################################
#
# createData.R
#
# copyright (c) 2001, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written Nov, 2001
#
# Modified in Mar 2004 for N-dye system
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# This is the function to make a data object from raw data object
#
######################################################################

createData <-
  function(rawdata, n.rep=1, avgreps=0, log.trans=TRUE)
{
  if (class(rawdata) != "rawdata")
    stop("You have to provide an object of class rawdata.")

  n.dye <- rawdata$n.dye
  n.array <- rawdata$n.array
  n.gene <- dim(rawdata$data)[1]/n.rep
  if(round(n.gene) != n.gene)
    stop("Number of rows in data file doesnot match the number of replicates")
  
  # do log2 transformation on one dye array
  if( (n.dye==1) & (log.trans) )
    warning("You are doing log2 transformation on one-dye data")
  
  # output object
  data <- NULL
  data <- rawdata
  data$n.gene <- n.gene
  data$n.rep <- n.rep
  data$collapse <- FALSE
  #data$log.trans <- log.trans
  
  # take single replicates for cloneid
  idx.gene <- seq(1,n.gene*n.rep,n.rep)
  #data$flag <- rawdata$flag[idx.gene,]
  # subset all other gene-based fields in rawdata
  fields.fixed <- c("n.dye","n.array","data", "flag", 
                    "metarow", "metacol", "row", "col", "design",
                    "TransformMethod")
  fields.all <- names(rawdata)
  # find the user specified fields in rawdata
  fields.user <- setdiff(fields.all, fields.fixed)
  # Note that all these fields are supposed to be gene-based,
  # that is, they have dimension (n.gene*n.rep) x 1.
  # When averaging the replicates, these fields will be subset
  if(length(fields.user) >= 1) {
    for(i in 1:length(fields.user)) {
      idx <- match(fields.user[[i]], fields.all)
      data[[idx]] <- data[[idx]][idx.gene]
    }
  }

  # reconstruct data$data if avgreps is not zero and rep is not 1
  # need to collapse the replicates
  if( (avgreps != 0) & (n.rep != 1) ) {
    # initialize data and flag
    data$data <- matrix(0, n.gene, n.array*n.dye)
    # n.rep became one (because of collapsing)
    data$n.rep <- 1
    data$collapse <- TRUE
    # drop grid location fields to save some memory,
    # they are no longer useful anyway
    data$metarow <- NULL
    data$metacol <- NULL
    data$row <- NULL
    data$col <- NULL
    # row indices
    idx <- seq(1,dim(rawdata$data)[[1]],n.rep)

    if(avgreps == 1) { # take mean of the replicates
      for(i in 0:(n.rep-1)) 
        data$data <- data$data + rawdata$data[idx+i,]
      data$data <- data$data/n.rep
    }
    if(avgreps == 2) { # take median of the replicates
      for(i in 1:(n.array*n.dye)) {
        for(j in 1:n.gene) {
          tmp <- rawdata$data[(n.rep*(j-1)+1):(n.rep*j),i]
          data$data[j,i] <- median(tmp)
        }
      }
    }
                                  
    # take flag
    if( !is.null(rawdata$flag) ) {
      data$flag <- matrix(0, n.gene, n.array)
      for(i in 0:(n.rep-1))
        data$flag <- data$flag + rawdata$flag[idx+i,]
    }
  }

  # log2 transform data if requested
  if(log.trans) {
    data$data <- log2(data$data)
    data$TransformMethod <- "log2"
  }
  # number of spots per array
  data$n.spot <- data$n.rep*data$n.array
  # column means
  data$colmeans <- apply(data$data, 2, mean)
  
  class(data) <- "madata"

  invisible(data)
  
}


print.madata <- function(x, ...)
{
  print.summary.madata(x, ...)
}


# summerize the MAdata object
summary.madata <-
  function(object, ...)
{
  if( is.na(match("madata",class(object))) )
    stop("The input variable is not an object of class madata!")

  TransformMethod <- object$TransformMethod
  
  result <- list(n.dye=object$n.dye, n.array=object$n.array,
                 n.gene=object$n.gene,
                 n.spot=object$n.spot, n.rep=object$n.rep,
                 TransformMethod=TransformMethod, collapse=object$collapse )
  class(result) <- "summary.madata"

  result
}


print.summary.madata <-
  function(x, ...)
{
  
  cat("\n\t\tSummary for this experiment\n\n")
  cat("Number of dyes:\t", x$n.dye, "\n")
  cat("Number of arrays:\t", x$n.array, "\n")
  cat("Number of genes:\t", x$n.gene, "\n")
  cat("Number of replicates:\t", x$n.rep, "\n")
  cat("Transformation method:\t",x$TransformMethod,"\n")
  cat("Replicate collapsed:\t",x$collapse, "\n")
  cat("\n\n")
  
}
