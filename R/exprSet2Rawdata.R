######################################################################
#
# exprSet2Rawdata.R
#
# copyright (c) 2001-2005, Hao Wu and Gary A. Churchill, The Jackson Lab.
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# This is the function to convert an object of exprSet class
# in BioConductor to Rawdata in R/maanova
#
######################################################################

exprSet2Rawdata <- function(exprdata, ndye=1, trans.method="None") 
{
  if(class(exprdata) != "exprSet")
    stop("The input must be an object of class exprSet")

  # make an Rawdata object from exprdata
  data <- NULL
  data$n.dye <- as.integer(ndye)
  data$data <- exprs(exprdata)
  
  # number of columns in data
  n.col <- ncol(data$data)
  # calculate number of arrays
  data$n.array <- n.col/ndye

  # number of rows
  n.row <- nrow(data$data)
  # make default metarow and metacol
  data$metarow <- rep(1, n.row)
  data$metacol <- rep(1, n.row)
  # transformation method
  data$TransformMethod <- trans.method
  
  # other gene dependent info
  data$GeneName <- geneNames(exprdata)

  # experimental design
  design <- pData(exprdata)
  # the rows of the design file must be in correct positions
  cellfile <- colnames(data$data)
  data$design <- design[cellfile,]

  
  class(data) <- "rawdata"
  
  data
}
  
    
