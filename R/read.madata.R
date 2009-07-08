######################################################################
#
# read.madata.R
#
# copyright (c) 2001-2002, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written May, 2002
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
######################################################################

read.madata <- function(datafile=datafile, designfile=designfile, covM = covM,
  arrayType=c("oneColor", "twoColor"),header=TRUE, spotflag=FALSE, n.rep=1, avgreps=0,
  log.trans=FALSE, metarow, metacol, row, col, probeid, intensity, matchDataToDesign=FALSE, ...){

  #================== make the output object 
  data <- NULL
  arrayType <- match.arg(arrayType); 
  if(arrayType  == 'oneColor') 
    cat(paste("Reading one color array.\n Otherwise change arrayType='twoColor' then read the data again\n"))
  else if(arrayType  == 'twoColor')
    cat(paste("Reading two color array.\n Array, Dye and Sample information should be in design file.\n 
    metarow, metacol, row, col, probeid, intensity information should be provided. \n 
    Log transformation ('log.trans=TRUE') is recommended. 
    If you miss some of information, provide them and read the data again.\n"))
  else  stop("Invalid arrayType")

  #============================ read design
  if(missing(designfile))
    stop("You must provide a designfile name or data.frame/matrix object for design file")

  if( is(designfile, 'character') )
    design <- read.table(designfile, sep="\t", quote="", header)
  else if(is(designfile, 'matrix')) design = designfile
  else if(is(designfile, 'data.frame')) design = designfile
  else stop("Invalid design file type")

  # Array (oneCol), Array, Dye and Sample (twoCol) must be in design 
  design.fields <- names(design)
  n.dye = 1
  if( !("Array" %in% design.fields) )
    stop("There is no Array column in design file")
  if(arrayType != 'oneColor'){
    if( !("Dye" %in% design.fields) )
      stop("There is no Dye column in design file")
    if( !("Sample" %in% design.fields) )
      stop("There is no Sample column in design file")
    # number of dyes
    n.dye <- length(unique(design$Dye))
  }
  if( "Spot" %in% design.fields )
      stop("You cannot have a column called 'Spot' in the design file")
  if("Label" %in% design.fields)
      stop("You cannot have a column called 'Label' in the design file")
  if("covM" %in% design.fields)
      stop("covM is reserved for covariate matrix.\n
      You cannot have a column called 'covM' in the design file")

  #============================= read data
  if(missing(datafile))
    stop("You must provide a file name or matrix object for data")
  
  if(is(datafile, 'character'))
  {
    rawdata<- as.matrix(read.table(datafile,sep="\t",
       quote="",header, comment.char=""))
    n.row <- nrow(rawdata)
   
    if(missing(probeid)){
      probeid = 1
      warning(paste("Assume that the first column is probeid. If you have probeid specify it, otherwise set 'probeid=0' then read the data again\n"))
    }
    if( length(probeid) == n.row ) data$probeid = probeid
    else{
      if( probeid > 0) data$probeid <- rawdata[,probeid]
      else{
        warning(paste("No Probe ID is available. 1 to ",n.row," is used."))
        data$probeid <- 1:n.row
      }
    }
    if(missing(intensity)){
      intensity = 2
      warning(paste("Assume that intensity value is saved from the second column. Otherwise provide 'intensity' (first column storing intensity) information, and read the data again\n"))
    }
  }
  else if(is(datafile, 'matrix') || is(datafile, 'ExpressionSet'))
  {
    # if it is an expression set just pull out the expression values
    if(is(datafile, 'ExpressionSet'))
    {
      datafile <- exprs(datafile)
    }
    
    n.row <- nrow(datafile)
    n.col <- ncol(datafile)
    if(missing(probeid)) probeid = -1
    if( length(probeid) == n.row ) data$probeid = probeid
    else{
      if( probeid <1 ){
        rname = rownames(datafile)
        if(length(rname)>0){
          data$probeid=rname; 
          cat(paste("Rowname of matrix is used for Probe ID\n"))
        }
        else{
          warning(paste("Datafile does not have rawname. 1 to",n.row,"is used for probeid."))
          data$probeid <- 1:n.row
        }
      }
      else data$probeid = datafile[,probeid]
    }      
    if(missing(intensity)){
       intensity = 1 
       cat("Assuming that intensity value is saved from the first column. ",
           "Otherwise provide 'intensity' (first column storing intensity) ",
           "information, and read the data again")
    }
    rawdata = datafile
    rm(datafile)
  }
  else
  {
    stop("Invalid data type. Provide either data file name (tab deliminated ",
         "format) or matrix object (class(object)==matrix)")
  }
  
  if(matchDataToDesign)
  {
      # try to match up the design file sample ordering with the expression set
      # ordering (otherwise the user has to do this themselves)
      if(spotflag)
      {
        arrayColSpan <- n.dye + 1
      }
      else
      {
        arrayColSpan <- n.dye
      }
      
      dataArrayNames <- colnames(rawdata)[seq(intensity, ncol(rawdata), by=arrayColSpan)]
      designArrayNames <- rle(as.character(design[, "Array"]))$values
      
      if(is.null(dataArrayNames) )
      {
          stop("matchDataToDesign should be set to FALSE and ",
               "expression set values should be manually ordered because ",
               "there are no array names in the data matrix")
      }
      else if(any(duplicated(dataArrayNames)))
      {
          stop("matchDataToDesign should be set to FALSE and ",
               "expression set values should be manually ordered because ",
               "the array names contain duplicates")
      }
      else if(any(duplicated(designArrayNames)))
      {
          stop("bad design file format. all array duplicates should be in ",
               "consecutive rows")
      }
      else if(length(dataArrayNames) != length(designArrayNames))
      {
          stop("expect the number of arrays in the design and input data to ",
               "match")
      }
      else
      {
          # use match for the initial reordering values
          newDataOrdering <- match(designArrayNames, dataArrayNames)
          if(any(is.na(newDataOrdering)))
          {
              stop("failed to match up the following design array names ",
                   "with data array names: ",
                   paste(designArrayNames[is.na(newDataOrdering)], sep=" ,"))
          }
          
          # account for the gaps in array column spans
          if(arrayColSpan > 1)
          {
              # we need to swith to 0 based indices for a sec
              # for this to work well
              newDataOrdering <- (newDataOrdering - 1) * arrayColSpan
              tmpOrdering <- c()
              for(i in 1:length(newDataOrdering))
              {
                  for(j in 1:arrayColSpan)
                  {
                      tmpOrdering <- c(tmpOrdering, newDataOrdering[i] + j)
                  }
              }
              
              newDataOrdering <- tmpOrdering
          }
          
          # insert any pre-intensity cols without any reordering
          if(intensity > 1)
          {
             # first shift the other indices to the right to make room
             newDataOrdering <- newDataOrdering + (intensity - 1)
             newDataOrdering <- c(1:(intensity - 1), newDataOrdering)
          }
          
          # do a final validity check
          if((length(newDataOrdering) == ncol(rawdata)) &&
             all(sort(newDataOrdering) == (1:length(newDataOrdering))))
          {
              # finally we can apply the new ordering to the data
              rawdata <- rawdata[, newDataOrdering]
          }
          else
          {
              stop("failed to match samples to design. please report this ",
                   "issue to maanova@jax.org")
          }
      }
  }
  
  n.row <- nrow(rawdata)
  n.col <- ncol(rawdata)
  n.gene <- n.row/n.rep
  if(round(n.gene) != n.gene)
    stop("Number of rows in data file does not match the number of replicates. Check n.rep or data file")

  #============================= read covM
  if( !missing(covM) ){
    if( is(covM, 'character') ){
      covm <- as.matrix(read.table(datafile,sep="\t", quote="",header, comment.char=""))
      if( (nrow(covm) != n.row | ncol(covm) != (n.col-intensity+1) ) ) 
        stop("Dimension of Covariate matrix does not match to the that of data matrix.\n
              Covariate matrix should (not) have header if data does (not)
  have header, and should have row name column.\n")
    }
    else{
      if( is(covM, 'matrix') ){
        covm = covM
        if( (nrow(covm) != n.row | ncol(covm) != (n.col-intensity+1)) ) 
          stop("Dimension of Covariate matrix does not match to the that of data matrix \n")
      }   
      else
        stop("Invalid covariate data type. Provide either file name (tab deliminated format)
             or matrix type R object (i.e., class(covM)==matrix)\n")
    }
    covM = covm[1,]
    design = cbind(design, covM)
  }
  else covm = NULL
  #====================== number of columns per array
  if(spotflag) ncol.array <- n.dye + 1
  else ncol.array <- n.dye
  n.array <- (n.col-intensity+1)/ncol.array
  if(round(n.array) != n.array)
    stop("Data has wrong number of columns")

  #====================== separate PM data and check if it has missing
  intensitydata <- matrix(as.numeric(rawdata[,intensity:n.col]),n.row, (n.col-intensity+1))
  
  #====================== divide genes by flag
  # If there's no flag info, intensity is the raw data
  # otherwise, we need to seperate intensity and flag
  if(spotflag){
    idx.flag <- ncol.array * (1:n.array)
    idx.intensity <- setdiff(1:(n.array*ncol.array), idx.flag)
    data$data <- intensitydata[,idx.intensity]
    data$flag <- intensitydata[,idx.flag]
  }
  else data$data <- intensitydata

  if(any(is.na(data$data))){
    class(data) <- "madata"
    data <- fill.missing(data)
    warning("the data contains NA values which have been filled in using the fill.missing(...) ",
            "function. If there should not have been any NA values please check raw data values ",
            "correct the NA's and re-read the data")
  }

  #== get metarow, metacol, row, col for transformation two col only : this is for transformation
  if(arrayType != 'oneColor'){
    if( !missing(metarow) )
      data$metarow <- as.integer(rawdata[, metarow])
    else {
      warning(paste("You have two color array - No meta row information, use 1 instead!\n"))
      data$metarow <- rep(1, n.row)
    }
  
    if( !missing(metacol) )
      data$metacol <- as.integer(rawdata[, metacol])
    else {
      warning(paste("You have two color array - No meta column information, use 1 instead!\n"))
      data$metacol <- rep(1, n.row)
    }
  
    if( !missing(row) )
      data$row <- as.integer(rawdata[, row])
    else {
      warning(paste("You have two color array - No row information.\n"))
    }
    if( !missing(col) )
      data$col <- as.integer(rawdata[, col])
    else {
      warning(paste("You have two color array - No column information.\n"))
    }
  } 
  data$collapse <- FALSE

  #======================== log transform
  data$TransformMethod = "None"
  if(log.trans){
    if(arrayType  == 'oneColor') warning(paste('You are taking log2 transformation to one color array \n'))
    data$data <- log2(data$data)
    data$TransformMethod <- "log2"
  }
  #================================  collapse if replicates > 1 
  if( (avgreps != 0) & (n.rep != 1) ) {
    # initialize data and flag
    tmpdata <- matrix(0, n.gene, n.array*n.dye)
    # n.rep became one (because of collapsing)
    
    # row indices
    idx <- seq(1,dim(data$data)[[1]],n.rep)
    data$probeid = data$probeid[idx] 
 
    if(avgreps == 1) { # take mean of the replicates
      for(i in 0:(n.rep-1)) 
          tmpdata <- tmpdata + data$data[idx+i,]
      data$data <- tmpdata/n.rep
    }
    if(avgreps == 2) { # take median of the replicates
      for(i in 1:(n.array*n.dye)) {
        for(j in 1:n.gene) {
          tmp <- data$data[(n.rep*(j-1)+1):(n.rep*j),i]
          data$data[j,i] <- median(tmp)
        }
      }
    }
    if( !is.null(data$flag) ) {
      tmpflag <- matrix(0, n.gene, n.array)
      for(i in 0:(n.rep-1))
        data$flag <- tmpflag + data$flag[idx+i,]
    }
    n.rep <- 1
    data$collapse <- TRUE
    # drop grid location fields to save some memory,
    # they are no longer useful anyway
    data$metarow <- NULL
    data$metacol <- NULL
    data$row <- NULL
    data$col <- NULL
  }
  #================== make the output object 
  data$n.dye <- n.dye
  data$n.array <- n.array
  data$design <- design
  data$n.gene <- n.gene
  data$n.spot <- n.rep*n.array
  data$colmeans <- apply(data$data, 2, mean)
  data$n.rep <- n.rep
  data$covm = covm
  #=============== other information : DO WE NEED THIS??
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

