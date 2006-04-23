######################################################################
#
# makeModel.R
#
# copyright (c) 2001, Hao Wu and Gary A. Churchill, The Jackson Lab.
#
# written Nov, 2001
# Modified Dec, 2002 for mixed model
# Modified Mar 2004 for N-dye system
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# This is the function to make model object for microarray experiment
#
######################################################################


makeModel <-
  function(data, design, formula, random=~1, covariate=~1)
{
  # check input variables
  if (class(data) != "madata")
    stop("You have to provide an object of class madata.")
  
  # create local variables
  ndyes <- data$n.dye
  nreps <- data$n.rep
  narrays <- data$n.array
  ngenes <- data$n.gene
  nspots <- data$n.spot
  # parse the input fixed and random formula
  parsed.formula <- parseformula(formula, random, covariate)
  nlabels <- length(parsed.formula$labels)
  
  # by default design will be in data
  if( missing(design) )
    design <- data$design
  
  # check design
  for(i in 1:length(design)) {
    if(length(design[[i]]) != narrays*ndyes)
      stop(paste(names(design)[i], "in design object has wrong length"))
  }

  # add new column(s) to design for interaction terms if it's not there
  if(nlabels > 1) {
    for(i in 1:nlabels) {
      if(parsed.formula$order[i] == 2) {
        if( !(parsed.formula$labels[i] %in% names(design)) ) {
          intterm <- which(parsed.formula$factors[,i]==1)
          newterm <- paste(design[[parsed.formula$labels[intterm[1]]]],
                           design[[parsed.formula$labels[intterm[2]]]], sep=":")
          design <- cbind(design, newterm)
          names(design)[length(design)] <- parsed.formula$labels[i]
        }
      }
    }
  }
  
  # cannot fit SG if there's no replicates or one dye
  if("Spot" %in% parsed.formula$labels) {
    if(nreps==1)
      stop("There is not technical replicates. You cannot fit Spot effect. Modify formula and try again.")
    if(ndyes==1)
      stop("You cannot fit Spot effect for one dye system. Modify formula and try again.")
  }
  # cannot fit labeling if there's no replicates or one dye
  if("Label" %in% parsed.formula$labels) {
    if(nreps==1)
      stop("There is not technical replicates. You cannot fit Labelling effect. Modify formula and try again.")
    if(ndyes==1)
      stop("You cannot fit Labelling effect for one dye system. Modify formula and try again.")
  }
  # make designobj, e.g., convert the fields of input design object to vectors
  designobj <- makeDesign(design)

  # check confounding problem, implement later

  # check if this is fixed or mixed model
  if( any(parsed.formula$random) )
    mixed <- 1
  else
    mixed <- 0

  ######################################################
  # make design matrices, e.g., X, Z,
  ######################################################
  # find the index for reference sample
  refid <- which(designobj$Sample==0)

  # initialization
  # first column in X are all ones, for mu
  X <- matrix(rep(1, ndyes*narrays*nreps),ncol=1)
  Z <- NULL
  dimZ <- NULL
  dimX <- NULL
  nesting <- NULL
  # make a column for reference samples (if any) in X
  if(length(refid) != 0) {
    tmp <- matrix(0, nreps, ndyes*narrays)
    tmp[,refid] <- 1
    X <- cbind(X, as.vector(tmp))
  }
  
  df <- NULL # degree of freedom of each term
  if(nlabels >= 1) {
    # it's not a completely null model, e.g., y~1
    # go throught all labels in formula to make X and Z matrices
    termX <- list(NULL)
    # loop for all terms in the formula
    for(i in 1:nlabels) {
      termX[[i]] <- numeric(0)
      l <- parsed.formula$labels[i]
      # if this is a first order term and
      # it's not in design and it's not Spot or Label, report error
      if( (parsed.formula$order[i] == 1) &
         !(l %in% c("Spot", "Label", names(designobj))) )
        stop(paste("Unrecognized term", l, "in formula"))
      
      # if this term is Spot
      if(l == "Spot") {
        replicate <- rep(1:nreps, ndyes*narrays)
        array <- rep(1:narrays, each=ndyes*nreps)
        spot <- nreps*(array-1)+replicate
        for(j in 1:nspots)
          termX[[i]] <- cbind(termX[[i]], spot==j)
      }
      # if this term is Label
      if(l == "Label") {
        label <- rep(1:(ndyes*narrays), each=nreps)
        for(j in 1:(ndyes*narrays))
          termX[[i]] <- cbind(termX[[i]], label==j)
      }
      # if this term is a factor in design
      if(l %in% names(designobj)) {
        # if l is covariate
        if(parsed.formula$covariate[i] == 1) {
          # note that covariate can be fixed or random
          # covariate must be numeric
          if(mode(as.vector(design[[l]])) != "numeric")
            stop(paste("Covariate", l, "is not numeric."))
          tmp <- rep(design[[l]], each=nreps)
          termX[[i]] <- matrix(tmp, ncol=1)
        }
        else {
          # l is not covariate
          mvarid <- designobj[[l]]
          variety <- as.vector(repmat(t(mvarid), nreps, 1))
          for(j in 1:max(mvarid))
            termX[[i]] <- cbind(termX[[i]], variety==j)
        }
      }
      # if this term is an interaction (for 2 interaction only now)
      if(parsed.formula$order[i] == 2) {
        intterm <- which(parsed.formula$factors[,i]==1)
      #  termX[[i]] <- intprod(termX, intterm)
      }
      
      ########## finish making design matrix for this term
      # make X and Z matrices
      if(parsed.formula$random[i] == 0) { # this term is fixed
        X <- cbind(X, termX[[i]])
        dimX <- c(dimX, dim(termX[[i]])[2])
      }
      else { # this term is random
        Z <- cbind(Z, termX[[i]])
        dimZ <- c(dimZ, dim(termX[[i]])[2])
      }
    }

    ##################################
    #calculate df for each term
    ##################################
    # find nesting relationships
    nesting <- matrix(0, nlabels, nlabels)
    rownames(nesting) <- parsed.formula$labels
    colnames(nesting) <- parsed.formula$labels
    if(nlabels > 1) {
      for(i in 1:(nlabels-1)) {
        nesting[i,i] <- 1
        for(j in (i+1):nlabels) {
          X.combine <- cbind(termX[[i]], termX[[j]])
          r <- matrank(X.combine)
          r1 <- matrank(termX[[i]])
          r2 <- matrank(termX[[j]])
          if(r == r1) { # term i nestes within term j
            nesting[j,i] <- 1
            nesting[i,j] <- -1
          }
          else if(r == r2) { # term j nestes within term i
            nesting[i,j] <- 1
            nesting[j,i] <- -1
          }
        }
      }
    }
    nesting[nlabels, nlabels] <- 1

    # start to calculate df for all terms
    df <- NULL
    # rank for X
    rank.X <- matrank(X)
    # rank for X when all terms are fixed
    rank.full <- matrank(cbind(X,Z))
    # index for fixed terms
    idx.fixed <- which(parsed.formula$random==0)
    # loop for all terms in the formula
    for(i in 1:nlabels) {
      X.drop <- matrix(rep(1, ndyes*narrays*nreps),ncol=1)
      # make a column for reference samples (if any) in X.drop
      if(length(refid) != 0) {
        tmp <- matrix(0, nreps, ndyes*narrays)
        tmp[,refid] <- 1
        X.drop <- cbind(X.drop, as.vector(tmp))
      }

      # find all terms nested with in the ith one
      idx.nest <- which(nesting[i,]==1)
      # combine all other termX together except the nesting ones
      for(j in setdiff(1:nlabels, idx.nest))
        X.drop <- cbind(X.drop, termX[[j]])
      df[i] <- matrank(cbind(X.drop, termX[[i]])) - matrank(X.drop)
    }
  
      
#      if(parsed.formula$random[i] ==0) { # for fixed terms
        # combine all other termX for fixed terms together
        # except the current one
#        for(j in setdiff(idx.fixed, i))
#          X.drop <- cbind(X.drop, termX[[j]])
#        df[i] <- rank.X - matrank(X.drop)
#      }
#      else { # for random terms
#        df[i] <- matrank(cbind(X, termX[[i]])) - rank.X
#      }
  }
  
  # calculate error df
  df <- c(df, dim(X)[1]-matrank(cbind(X,Z)))
  names(df) <- c(parsed.formula$labels, "Error")


  # make return object
  result <- NULL
  result$X <- X
  result$Z <- Z
  result$dimX <- dimX
  result$dimZ <- dimZ
  result$df <- df
  result$mixed <- mixed
  result$design <- design
  result$formula <- formula
  result$random <- random
  result$covariate <- covariate
  result$parsed.formula <- parsed.formula
  result$nesting <- nesting
  class(result) <- "mamodel"

  invisible(result)
}

summary.mamodel <- function(object, ...)
{
  if(class(object) != "mamodel")
    stop("The input variable is not an object of class mamodel!")

  result <- NULL

  # is it mixed or fixed model
  if(object$mixed == 1)
    result$mixed <- "mixed"
  else
    result$mixed <- "fixed"
  
  # model information
  result$formula <- as.character(object$formula)[2]
  # random
  r <- as.character(object$random)[2]
  if(r == "1")
    result$random <- "None"
  else
    result$random <- r
  # covariate
  cov <- as.character(object$covariate)[2]
  if(cov=="1")
    result$covariate <- "None"
  else
    result$covariate <- cov
  
  # class level information
  ref <- NULL; nref <- NULL; effref <- NULL
  if(any(object$design$Sample==0)) {
    ref <- "Reference"
    nref <- 1
    effref <- "fixed"
  }
  Class <- c(ref, object$parsed.formula$labels[object$parsed.formula$random==0],
             object$parsed.formula$labels[object$parsed.formula$random==1])
  Levels <- c(nref, object$dimX, object$dimZ)
  Effect <- c(effref,rep("fixed", sum(object$parsed.formula$random==0)),
              rep("random",sum(object$parsed.formula$random==1)))
  for(i in 1:length(Class)) {
    if(Class[i] != "Reference") {
      idx <- which(object$parsed.formula$labels==Class[i])
      if(object$parsed.formula$covariate[idx] == 1)
        # this term is a covariate
        Effect[i] <- "Covariate"
    }
  }
  clevel <- as.data.frame(cbind(Class, Levels, Effect))
  result$clevel <- clevel
  
  # dimension info
  diminfo <- NULL
  diminfo$nobs <- dim(object$X)[1]
  diminfo$ncolX <- dim(object$X)[2]
  diminfo$ncolZ <- dim(object$Z)[2]
  result$diminfo <- diminfo
  result$df <- as.data.frame(object$df)
  colnames(result$df) <- "df"
  class(result) <- "summary.mamodel"
  result
}

print.summary.mamodel <- function (x, ...)
{
  cat("\n\tModel Summary\n\n")
  cat("This is a", x$mixed, "effect model", "\n\n")
  cat("Gene-specific ANOVA model:\t", x$formula, "\n")
  cat("Gene-specific Random terms:\t", x$random, "\n\n")
  cat("Gene-specific covariate:\t", x$covariate, "\n\n")
  cat("\tClass Level Information\n\n")
  print(x$clevel)
  cat("\n")
  cat("\tDimensions\n\n")
  cat("Observations(per gene):\t", x$diminfo$nobs, "\n")
  cat("Columns in X:\t", x$diminfo$ncolX, "\n")
  if(x$mixed == "mixed")
    cat("Columns in Z:\t", x$diminfo$ncolZ,"\n")
#  cat("\n")
#  cat("\tDegree of freedoms:\n\n")
#  print(x$df)
  cat("\n\n")
}
                 

###############################################################
# function to make a integer list from input design object
###############################################################
makeDesign <- function(design)
{
  # find the idx for non-references
  idx.noref <- which(design$Sample!=0)
  
  # output
  result <- list(NULL)
  
  for(i in 1:length(design)) {
    field <- design[[i]]
    level <- as.vector(sort(unique(field[idx.noref])))
    mvarid <- rep(0, length(field))
    for(j in 1:length(level))
      mvarid[field==level[j]] <- j
    result[[i]] <- mvarid
  }

  names(result) <- names(design)
  
  result
}


