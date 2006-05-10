######################################################################
#
# util.R
#
# copyright (c) 2001, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written Nov, 2001
# 
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# This file contains the following functions:
#
# matsort, zeros, ones, repmat, blkdiag, num2yn, matrank, mixed, etc.
#
######################################################################

matsort <- function(mat, index=1)
{
  result <- mat
  n <- dim(mat)
  if(index==1) {
    for(i in 1:n[2])
      result[,i] <- sort(mat[,i])
  }
  else if(index==2) {
    for(i in 1:n[1])
      result[i,] <- sort(mat[i,])
  }

  result
}

    
zeros <- function(dim)
{
  n.dim <- length(dim)
  
  dim.total <- 1
  for(i in 1:n.dim)
    dim.total <- dim.total*dim[i]
  
  result <- array(rep(0,dim.total),dim=dim)
  result
}


ones <- function(dim)
{
  n.dim <- length(dim)
  
  dim.total <- 1
  for(i in 1:n.dim)
    dim.total <- dim.total*dim[i]
  
  result <- array(rep(1,dim.total),dim=dim)
  result
}


repmat <- function(mat, n.row, n.col, ...)
{
  tmp <- NULL
  result <- NULL
  for (i in 1:n.col)
    tmp <- cbind(tmp, mat)
  for(i in 1:n.row)
    result <- rbind(result,tmp)
  
  result

}


blkdiag <- function(...)
{
  # get the variables
  argus <- list(...)
  nargu <- length(argus)
  
  nrow <- NULL
  ncol <- NULL
  result.nrow <- 0
  result.ncol <- 0
  
  # go thru the variables
  for (i in 1:nargu) {
    mat <- argus[[i]]
    if(is.null(mat)) {
      nrow[i] <- 0
      ncol[i] <- 0
    }
    else if( is.vector(mat) ) {
      nrow[i] <- 1
      ncol[i] <- length(mat)
    }
    else if( is.matrix(mat) ) {
      nrow[i] <- dim(mat)[1]
      ncol[i] <- dim(mat)[2]
    }
    else
      stop(paste("Wrong data type for the number",i,"input variable"))
    
    result.nrow <- result.nrow + nrow[i]
    result.ncol <- result.ncol + ncol[i]
  }

  # construct the result matrix
  result <- matrix( rep(0, result.nrow*result.ncol),
                   nrow=result.nrow, ncol=result.ncol )

  row.start <- 1
  col.start <- 1
  for(i in 1:nargu) {
    if(!is.null(argus[[i]])) {
      row.end <- row.start + nrow[i] - 1
      col.end <- col.start + ncol[i] - 1
      result[row.start:row.end, col.start:col.end] <- argus[[i]]
      row.start <- row.end + 1
      col.start <- col.end + 1
    }
  }
  
  result
}



# function to turn the logical value to "Yes" or "No"
# this is used for several summary functions
num2yn <- function(x, ...)
{
  if(x) result <- "Yes"
  else result <- "No"

  result
}


# function to find the maximum/minumum value for row/columns
# input variable is a matrix (2-dimensional array)
rowmax <- function(x)
{
  k <- dim(x)[1]
  result <- rep(0, k)
  for(i in 1:k)
    result[i] <- max(x[i,])

  result
}

rowmin <- function(x)
{
  k <- dim(x)[1]
  result <- rep(0, k)
  for(i in 1:k)
    result[i] <- min(x[i,])
  result
}

colmax <- function(x)
{
  k <- dim(x)[2]
  result <- rep(0, k)
  for(i in 1:k)
    result[i] <- max(x[,i])
  result
}

colmin <- function(x)
{
  k <- dim(x)[2]
  result <- rep(0, k)
  for(i in 1:k)
    result[i] <- min(x[,i])
  result
}


# calculate the sum of rows for a given matrix
sumrow <- function(x)
{
  k <- dim(x)[1]
  result <- rep(0, k)
  for(i in 1:k)
    result[i] <- sum(x[i,])

  result
}


# get the matrix rank (this is no this function in R???
# I just couldn't find it)
matrank <- function(X)
{
  if(is.vector(X))
    # rank is 1 for vectors
    return(1)
  else if(is.matrix(X)) {
    s <- La.svd(X,0,0)
    tol <- max(dim(X)) * s$d[1] * .Machine$double.eps
    r <- sum(s$d>tol)
    return(r)
  }
  else{
    return(0)
  }
}


# calculate matrix or vector norm, working only for vector now
norm <- function(X)
{
  result <- sqrt(sum(X*X))
  result
}


###############################################################
# function to parse the fixed and random formula
###############################################################
parseformula <- function(formula, random, covariate)
{
  formula.terms <- terms(formula)
  random.terms <- terms(random)
  cov.terms <- terms(covariate)
  # get the labels
  formula.labels <- attr(formula.terms, "term.labels")
  random.labels <- attr(random.terms, "term.labels")
  cov.labels <- attr(cov.terms, "term.labels")
  
  # make output object
  result <- NULL
  result$labels <- formula.labels
  result$factors <- attr(formula.terms, "factors")
  result$order <- attr(formula.terms, "order")
  result$random <- rep(0, length(formula.labels))
  result$covariate <- rep(0, length(formula.labels))
  
  # find random terms
  if(length(random.labels) != 0) {
    for(i in 1:length(random.labels)) {
      idx <- which(random.labels[i]==formula.labels)
      # check data, e.g., all labels in random have to be in formula
      if( length(idx) == 0 )
        stop(paste("Random term", random.labels[i], "is not in formula."))
      else
        result$random[idx] <- 1
    }
  }
  # for higher order terms, if any main factor is random,
  # the interaction is random
  if(any(result$order>1)) {
    idx <- which(result$order>1)
    for(i in 1:length(idx)) {
      # find the indices for main terms
      idx.mainterm <- which(result$factors[,idx[i]]==1)
      if(any(result$random[idx.mainterm]))
        result$random[idx[i]] <- 1
    }
  }

  # find covariates
  if(length(cov.labels) != 0) {
    for(i in 1:length(cov.labels)) {
      idx <- which(cov.labels[i]==formula.labels)
      if(length(idx) == 0)
        stop(paste("Covariate", cov.labels[i], "is not in formula."))
      else
        result$covariate[idx] <- 1
    }
  }

  result
}

###############################################################
# function to make the design matrix for interaction terms
# it's working for two way interaction only
###############################################################
intprod <- function(terms, intterm)
{
  mat1 <- terms[[intterm[1]]]
  mat2 <- terms[[intterm[2]]]

  n1 <- dim(mat1)[2]
  n2 <- dim(mat2)[2]

  # init result
  result <- matrix(0, dim(mat1)[1], n1*n2)

  for(i in 1:n1) {
    for(j in 1:n2) {
      result[,(i-1)*n2+j] <- mat1[,i] * mat2[,j]
    }
  }

  result
}

  
###############################################################
# function to make comparison matrix given number of levels
###############################################################
makeCompMat <- function(n)
{
  if(n==1) return(1)
  # if the number of levels is n,
  # the result matrix has dimension (n-1) x n
  C <- matrix(0, n-1, n)
  C[,1] <- 1
  for(i in 1:(n-1))
    C[i, i+1] <- -1

  C
}

###############################################################
# function to find subgroups in unconnected design
# this code is messy. I may rewrite it later
###############################################################
findgroup <- function(varid, ndye)
{
  # local variables
  vargroup <- list(NULL)
  ngroups <- 1
  nvars <- length(varid)
  narrays <- nvars / ndye
  
  # varids on the first array belong to the first group
  vargroup[[ngroups]] <- varid[1:ndye]
  
  # loop for the rest arrays
  if(narrays > 1) {
    for(i in 2:narrays) {
      newgroup <- 1
      varid.thisarray <- varid[((i-1)*ndye+1):(i*ndye)]
      for(j in 1:ngroups) {
        if( any(varid.thisarray %in% vargroup[[j]]) ) {
          # varid.thisarray belong to vargroup j
          vargroup[[j]] <- c(vargroup[[j]], varid.thisarray)
          newgroup <- 0
          break;
        }
      }
      # if varid.thisarray doesn't belong to any existing group
      # make a new group
      if(newgroup) {
        ngroups <- ngroups + 1
        vargroup[[ngroups]] <- varid.thisarray
      }
    }
  }

  if(length(vargroup) == 1)
    finalgroup <- vargroup
  else {
    # start to do merge
    finalgroup <- list(NULL)
    for(i in 1:(length(vargroup)-1)) {
      for(j in (i+1):length(vargroup)) {
        if(length(intersect(vargroup[[i]], vargroup[[j]])) > 0) {
          # merge this two
          vargroup[[i]] <- c(vargroup[[i]], vargroup[[j]])
          vargroup[[j]] <- numeric(0)
        }
      }
    }
    ngroups <- 1
    for(i in 1:length(vargroup)) {
      if(length(vargroup[[i]]) > 0) {
        finalgroup[[ngroups]] <- vargroup[[i]]
        ngroups <- ngroups + 1
      }
    }
  }
      
  # arrange finalgroup - not include zeros in the group
  for(i in 1:length(finalgroup)) 
    finalgroup[[i]] <- setdiff(sort(unique(finalgroup[[i]])), 0)

  # return
  finalgroup
    
}

###############################################################
#
# function to calculate the pseudo-inverse of a singular matrix.
# Note that I was using ginv function in MASS but it is not
# robust, e.g., sometimes have no result. That's because
# the engine function dsvdc set the maximum number of iteration
# to be 30, which is not enough in some case.
# I use La.svd instead of svd in my function.
# I don't want to spend time on it so it doesn't support complex number
#
# Note that method "dgesdd" in La.svd gave me troubles sometimes.
# I couldn't figure out why. I'm using "dgesvd" instead. Seems it generates
# the exact same result as in Matlab. But What's the disadvantage of it?
#
# From R-2.3.0 La.svd(X,method="dgesvd") is deprecated, so change 
# it to "dgesdd". Hopefully it will not have problem like before.
#
# From R-2.4.0 La.svd(X, method="dgesdd" or "dgesvd") is deprecated, so change
# it to La.svd(X). method="dgesdd" is the default.
#
###############################################################
  
pinv <- function(X, tol)
{
  if ( length(dim(X)) > 2 || !(is.numeric(X)) ) 
    stop("X must be a numeric matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xsvd <- La.svd(X)
  
  # find the tolerance
  if(missing(tol))
    tol <- max(dim(X)) * Xsvd$d[1] * .Machine$double.eps
  
  Positive <- Xsvd$d > tol
  if (all(Positive)) 
    t(Xsvd$vt) %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2:1])
  else t(Xsvd$v)[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                            t(Xsvd$u[, Positive, drop = FALSE]))
}

###############################################################
#
# fdr - function to calculate the adjusted P values
# for FDR control
#
###############################################################
fdr <- function(p, method=c("stepup", "adaptive", "stepdown"))
{
  method <- match.arg(method)

  m <- length(p)
  tmp <- sort(p, index.return=TRUE)
  sortp <- tmp$x
  idx <- tmp$ix

  if(method == "stepdown") {
    d <- m:1
    sortp <- (1-(1-sortp)^d) * d/m

    for(i in 1:(m-1)) {
      if(sortp[i+1] < sortp[i])
        sortp[i+1] <- sortp[i]
    }
  }
  else {
    if(method == "stepup")
      m0 <- m
    else if(method == "adaptive") {
      s <- sort(1 - sortp)/(1:m)
      # calculate m0
      m0raw <- m
      i <- m
      while(i > 1 && s[i] <= s[i - 1]) i <- i - 1
      if(i > 1)
      m0raw <- 1/s[i - 1]
      else m0raw <- 1/s[1]
      m0 <- min(floor(1 + m0raw), m)
    }


    # calculate sortp
    sortp <- sortp * m0 / (1:m)
    for(i in (m-1):1) {
      if(sortp[i] > sortp[i+1])
        sortp[i] <- sortp[i+1]
    }
  }

  # return variable
  result <- NULL
  result[idx] <- sortp

  result
}

###############################################################
# meanvarlog - function to generate mean and var for a logrithm
# chi2 distribution
# This is used by JSshrinker
###############################################################

meanvarlog <- function(df)
{
  meanlog <- rep(0,length(df))
  varlog <- rep(0, length(df))
  for(i in 1:length(df)) {
    Chis <- rchisq(100000, df[i])
    # mean
    meanlog[i] = mean(log(Chis/df[i]))
    # variance
    varlog[i] = var(log(Chis))
  }

  # output
  result <- NULL
  result$meanlog <- meanlog
  result$varlog <- varlog

  result
}

###############################################################
# engine function for JSshrinker
# The input need to be a column vector
###############################################################
JS <- function(X, var)
{
  m=mean(X)
  out=sum((X-m)*(X-m))
  p=length(X)
  out=1-(p-3)*var/out
  out=max(out,0)
  out=m+out*(X-m)

  out
}


###############################################################
# James-Stein shrinkage estimator
###############################################################
JSshrinker <- function(X, df, meanlog, varlog)
{
  # X could be a vector or matrix, convert it to matrix if it's not
  if(!is.matrix(X))
    X <- matrix(X)

  # check the length of df 
  if(dim(X)[2] != length(df))
    stop("Degree of freedom vector has wrong length")
  
  if(missing(meanlog)|missing(varlog)) {
    tmp <- meanvarlog(df)
    meanlog <- tmp$meanlog
    varlog <- tmp$varlog
  }
  
  #initialize result
  result <- X
  # loop for columns
  for (i in 1:length(df)) {
    DATA <- X[,i]
    DATA <- log(DATA)
    DATA <- DATA - meanlog[i]
    mans <- JS(DATA,varlog[i])
    result[,i] <- exp(mans);
  }

  result
}

 
#################################################
# the following functions are useful for Jmaanova
#################################################

######################################################################
#
# make.ratio.R
# Calculate the logratio for two dye arrays
#
######################################################################
make.ratio <- function(object, norm.std=TRUE)
{
  # local variables
  n.array <- object$n.array
  
  if( class(object) == "madata" ) {
    x <- object$data
    colmeans <- object$colmeans
  }
  else if(class(object)  == "rawdata") {
    x <- log2(object$data)
    colmeans <- apply(x, 2, mean)
  }
  else
    stop("The first input variable is not an object of class madata or rawdata!")
  
  # calculate column means
  #colmeans <- apply(x, 2, mean)
  
  if(object$n.dye != 2)
    stop("make.ratio only works for 2-dye arrays")
  
  n.row <- dim(x)[[1]]
  n.col <- dim(x)[[2]]

  # calculate the standard deviation for each column
  std.array <- NULL
  for(i in 1:(n.array*2))
      std.array[i] <- sd(x[,i])

  # call C function to normalize the data
  z <- .C("makeratio",
          as.double(x), # adjusted data
          as.double(colmeans), # column means
          as.double(std.array), # column standard deviation
          as.integer(norm.std), # flag to divided by std or not
          as.integer(n.row), # number of rows in the data
          as.integer(n.col), # number of columns in the data
          result=as.double(rep(0, n.row*n.col/2)), # return variable
          PACKAGE="maanova"
          )
  
  result <- matrix(z$result, n.row, n.col/2)

  result
}

