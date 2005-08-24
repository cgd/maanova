######################################################################
#
# read.madata.R
#
# copyright (c) 2001-2004, Hao Wu and Gary A. Churchill, The Jackson Lab.
#
# written Nov, 2001
# modified Nov, 2002 for reading design file for mixed model effect
# modified Mar, 2004 for N-dye system
# modified May, 2004 for collapse.madata function
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# This is the function to read in MicroArray experiment data and design
# from two tab delimited files
#
######################################################################

read.madata <-
  function(datafile, designfile="design.txt", header=TRUE, spotflag=TRUE,
           metarow, metacol, row, col, pmt, ...) 
{
  # read in design file first
  design <- read.table(designfile,sep="\t", quote="", header=TRUE)
  # Array, Dye and Sample must be in design
  design.fields <- names(design)
  if( !("Array" %in% design.fields) )
    stop("There is no Array column in design file")
  if( !("Dye" %in% design.fields) )
    stop("There is no Dye column in design file")
  if( !("Sample" %in% design.fields) )
    stop("There is no Sample column in design file")
  if( "Spot" %in% design.fields )
    stop("You cannot have a column called Spot in the design file")
  if("Label" %in% design.fields)
    stop("You cannot have a column called Label in the design file")

  # number of dyes
  n.dye <- length(unique(design$Dye))
  
  # read in data
  if(missing(datafile))
    stop("You must provide a file name for raw data")
  rawdata <- as.matrix(read.table(datafile,sep="\t",quote="",header,
                                  comment.char=""))
  n.row <- nrow(rawdata)
  n.col <- ncol(rawdata)

  # make the output object
  data <- NULL
  data$n.dye <- n.dye
  # number of columns per array
  if(spotflag) ncol.array <- n.dye + 1
  else ncol.array <- n.dye
  n.array <- (n.col-pmt+1)/ncol.array
  if(round(n.array) != n.array)
    stop("Data has wrong number of columns for pmt data")
  data$n.array <- n.array
  
  # pmt data with/without flags
  pmtdata <- matrix(as.numeric(rawdata[,pmt:n.col]),
                     n.row,n.col-pmt+1)
  # If there's no flag info, pmt is the raw data
  # otherwise, we need to seperate pmt and flag
  if(spotflag) {
    idx.flag <- ncol.array * (1:n.array)
    idx.pmt <- setdiff(1:(n.array*ncol.array), idx.flag)
    data$data <- pmtdata[,idx.pmt]
    data$flag <- pmtdata[,idx.flag]
  }
  else{
    data$data <- pmtdata
  }

  # get the CloneID
  #if( !missing(cloneid) )
  #  data$cloneid <- rawdata[,cloneid]
  #else {
  #  warning(paste("Clone ID is not provided. 1 to",n.row,"is used."))
  #  data$cloneid <- 1:n.row
  #}

  # get metarow, metacol, row, col
  if(!missing(metarow))
    data$metarow <- as.integer(rawdata[, metarow])
  else {
    warning("No meta row information, use 1 instead!")
    data$metarow <- rep(1, n.row)
  }
  
  if(!missing(metacol))
    data$metacol <- as.integer(rawdata[, metacol])
  else {
    warning("No meta column information, use 1 instead!")
    data$metacol <- rep(1, n.row)
  }
  
  if(!missing(row))
    data$row <- as.integer(rawdata[, row])
  else {
    warning("No row information.")
  }
  if(!missing(col))
    data$col <- as.integer(rawdata[, col])
  else {
    warning("No column information.")
  }

  # get other data from ...
  args <- list(...)
  nargu <- length(args)
  # number of fileds in data
  n.field <- length(names(data))
  
  if(nargu) { # if there's any additional arguments
    for(i in 1:nargu) { #
      argname <- names(args[i]) # argument name
      argvalue <- args[[i]] # argument value
      if( !is.numeric(argvalue) | argvalue<=0 | argvalue != round(argvalue) )
        stop(paste("The value for",argname,"must be an positive integer"))
      # create new field in data object and read in the data
      data$tmp <- as.character(rawdata[,argvalue])
      # change the field name to argname
      names(data)[n.field+i] <- argname
    }
  }

  # skip reading array name
  # read in array names (if any)
#  if(header == T) {
#    tmp <- strsplit(readLines(datafile,1), "\t")[[1]]
#    tmp <- tmp[pmt:length(tmp)]
#    if(spotflag)
#      data$ArrayName <- tmp[idx.pmt]
#    else
#      data$ArrayName <- tmp
#  }

  # no data transformation at this time
  data$TransformMethod <- "None"

  data$design <- design
  class(data) <- "rawdata"

  data
}
  
