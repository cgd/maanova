######################################################################
#
# write.madata.R
#
# copyright (c) 2001, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written Dec, 2001
# Modified Apr, 2004 for N-dye system
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
# This is the function to write an object of madata to files
#
######################################################################

write.madata <- function(madata, datafile="madata.txt", designfile="design.txt")
{
  if( (class(madata) != "madata") && (class(madata) != "rawdata") )
    stop("The first input variable must be an object of madata or rawdata!")

  # local variables
  nreps <- madata$n.rep
  ndyes <- madata$n.dye
  
  if(is.null(nreps))
    # for rawdata, there's no nreps, make it 1
    nreps <- 1
  
  ##############################
  # write the data file
  ##############################
  # find the fields in data shouldn't be in the text file
  field.fixed <- c("n.dye", "n.array", "n.gene", "n.rep", "n.spot",
                   "flag", "design", "collapse",
                   "colmeans", "offset", "TransformMethod",
                   "ArrayName", "data")
  field.all <- names(madata)
  field.col <- setdiff(field.all, field.fixed)

  # find the gene-specific field and spot specific field
  # note that the spot-specific field has nothing to do with replicates
  field.spot <- c("metarow", "metacol", "row","col", "Flag", "data")
  field.spot <- intersect(field.spot, field.col)
  field.gene <- setdiff(field.col, field.spot)
  
  # create the column header and the data frame to be written
  data <- NULL
  for(i in 1:length(field.col)) {
    if(field.col[i] %in% field.spot)
      data <- cbind(data, madata[[field.col[i]]])
    else if(field.col[i] %in% field.gene)
      data <- cbind(data, rep(madata[[field.col[i]]], each=nreps))
  }
  colnames(data) <- field.col

  dyename <- NULL
  for(i in 1:ndyes)
    dyename <- c(dyename, paste("Dye", i, sep=""))
  
  # bind pmt data and flags to data
  for(i in 1:madata$n.array) {
    if(ndyes > 1)
      tmp <- madata$data[,(i*ndyes-1):(i*ndyes)]
    else
      tmp <- matrix(madata$data[,i], ncol=1)
#    if("ArrayName" %in% names(madata))
#      tmpname <- madata$ArrayName[(i*2-1):(i*2)]
#    else
    # User Array1Dye1 ... as column headers for data
    tmpname <- paste("Array", i, dyename, sep="")
    if("flag" %in% names(madata)) {
      tmp <- cbind(tmp, madata$flag[,i])
      tmpname <- c(tmpname, paste("Array", i, "Flag",sep=""))
    }
    # give it column names
    colnames(tmp) <- tmpname
    # bind tmp to data
    data <- cbind(data, tmp)
  }
  
  write.table(data, datafile, sep="\t", row.names=FALSE,
              col.names=TRUE, quote=FALSE)

  #################################
  # write the design file
  #################################
  if("design" %in% names(madata))
    write.table(madata$design, designfile, sep="\t", row.names=FALSE,
              col.names=TRUE, quote=FALSE)
}
  
