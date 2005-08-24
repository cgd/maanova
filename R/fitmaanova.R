######################################################################
#
# fitmaanova.R
#
# copyright (c) 2001, Hao Wu and Gary A. Churchill, The Jackson Lab.
#
# written Nov, 2001
# Modified Dec, 2002 for mixed effect model
# Modified Mar, 2004 for N-dye system
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
# This is the function to fit ANOVA model
#
######################################################################

fitmaanova <-
  function(madata, mamodel, inits20,
           method=c("REML","ML","MINQE-I","MINQE-UI", "noest"),
           verbose=TRUE)

{
  if (class(madata) != "madata")
    stop("The first input variable is not an object of class madata.")
  if (class(mamodel) != "mamodel")
    stop("The second input variable is not an object of class mamodel.")
  method <- match.arg(method)

  # create local variables
  data <- madata$data
  ndyes <- madata$n.dye
  narrays <- madata$n.array 
  ngenes <- madata$n.gene
  nspots <- madata$n.spot 
  nreps <- madata$n.rep 
  #colmeans <- apply(madata$data, 2, mean)
  # take colmeans from madata
  colmeans <- madata$colmeans
  X <- mamodel$X
  Z <- mamodel$Z
  dimX <- mamodel$dimX
  dimZ <- mamodel$dimZ
  design <- mamodel$design
  formula <- mamodel$formula
  cov <- mamodel$covariate
  parsed.formula <- mamodel$parsed.formula

  # check if we have enough df to do the fit
  nobs <- ndyes*narrays*nreps
  nrankX <- matrank(X)
  if(nrankX >= nobs)
    stop("Not enough degree of freedom to do regression. You need to remove someterms in your model.")
  # find reference id
  refid <- which(design$Sample==0)
  
  # create the output object
  anova <- NULL
#  anova$rss <- rep(0, ngenes)
  anova$yhat <- matrix(0, ngenes*nreps, ndyes*narrays)
  anova$S2 <- matrix(0, ngenes, length(dimZ)+1)
  # error df
  error.df <- dim(mamodel$X)[1] - matrank(mamodel$X)
  # for column-center the data
  colmeans.gene <- repmat(t(colmeans), nreps, 1)
  # fit the model
  if(mamodel$mixed == 0) {
    # fixed model
    Allb <- matrix(0, ngenes, dim(X)[2])
    invX <- pinv(X)
    for(i in 1:ngenes) {
      # take the data for this gene
      y <- data[(nreps*(i-1)+1):(nreps*i),] - colmeans.gene
      y <- as.vector(y)
      b <- invX %*% y
      Allb[i,] <- b
      yfit <- X %*% b
      # residual sum of squares
      rss <- sum((y-yfit)*(y-yfit))
      anova$S2[i] <- rss/error.df
      # fitted values
      anova$yhat[(nreps*(i-1)+1):(nreps*i),] <- matrix(yfit, nrow=nreps) +
        colmeans.gene
    }
  }
  else {
    # this is mixed model
    ######################################################################
    # fit fixed model first then to calculate the variance components.
    # The variance is used as the initial value for mixed effect model
    # fitting. This can speed up the convergence of ML/REML
    # if there's no enough df in fixed model, skip that
    #
    # This part will be skipped if the initial variance components
    # were provided
    ######################################################################
    if(missing(inits20)) {
      idx.random <- which(parsed.formula$random==1)
      nrandom <- length(idx.random)
      if( length(refid)==0 ) # no reference sample
        bstart <- 1
      else # have reference sample
        bstart <- 2
    
      # make a fixed model and fit it
      model.fixed <- makeModel(madata, design, formula, covariate=cov)
      nrankX <- matrank(model.fixed$X)
      if(nrankX >= nobs) {
        # no enough df to fit fixed model, skip variance calculation
        # and use a default one
        inits20 <- matrix(0.01, ngenes, nrandom+1)
      }
      else {
        if(verbose==TRUE)
          cat("Calculating variance components for fixed model...\n")
        invX <- pinv(model.fixed$X)
        inits20 <- matrix(0, ngenes, nrandom+1)
        for(i in 1:ngenes) {
          y <- data[(nreps*(i-1)+1):(nreps*i),] - colmeans.gene
          y <- as.vector(y)
          tmpb <- invX %*% y
          yfit <- model.fixed$X %*% tmpb
          residual <- y - yfit
          # calculate variance for random factors
          for(j in 1:nrandom) {
            if(idx.random[j]>1)
              boffset <- sum(model.fixed$dimX[1:(idx.random[j]-1)])
            else
              boffset <- 0
            lo <- bstart+boffset+1
            nn <- model.fixed$dimX[idx.random[j]]
            tmp <- tmpb[lo:(lo+nn-1)]
            inits20[i,j] <- var(tmp)
          }
          # calculate error variance
          inits20[i,nrandom+1] <- var(residual)
        }
        # sometimes there are numeric problems if the initial s2
        # are too small. I convert the small numbers to 1e-5 here
        inits20[inits20<1e-5] <- 1e-5
      }
    }
    
    if(verbose==TRUE)
      cat("Fitting mixed effect model...\n")
    # initialize
    Allb <- matrix(0, ngenes, dim(X)[2])
    Allu <- matrix(0, ngenes, dim(Z)[2])
#    anova$loglik <- rep(0, ngenes)
    anova$loops <- rep(0, ngenes)
    anova$S2 <- matrix(0, ngenes, length(dimZ)+1)
    anova$S2.level <- parsed.formula$labels[parsed.formula$random==1]
    # prepare some variables for mixed model fitting
    # these values will be constant for all genes so I precompute
    # them to save some time
    XX <- t(X) %*% X
    XZ <- t(X) %*% Z
    ZZ <- t(Z) %*% Z
    Zi <- makeZiZi(Z, dimZ)
    method <- match.arg(method)
    # loop thru all genes
    for(i in 1:ngenes) {
      if(verbose == TRUE)
        if(round(i/100) == i/100)
          cat("Finish gene number", i, "...\n")
      y <- data[(nreps*(i-1)+1):(nreps*i),] - colmeans.gene
      y <- as.vector(y)
      generesult <- mixed(y, X, Z, XX, XZ, ZZ, Zi$Zi, Zi$ZiZi,
                           dimZ, inits20[i,], method)
      yfit <- X%*%generesult$b + Z%*%generesult$u
      # residual sum of squares
#      anova$rss[i] <- sum((y-yfit)*(y-yfit))
      # fitted values
      anova$yhat[(nreps*(i-1)+1):(nreps*i),] <- matrix(yfit, nrow=nreps) +
        colmeans.gene
      # other parameters
      anova$S2[i,] <- generesult$s2
#      anova$loglik[i] <- generesult$loglik
      anova$loops[i] <- generesult$loops
      Allb[i,] <- generesult$b
      Allu[i,] <- generesult$u
    }
  }
      

  # get the estimates from Allb
  anova$G <- Allb[,1]
  next.fix <- 2
  if(length(refid) != 0) {
    # there's references
    anova$reference <- Allb[,2]
    next.fix <- next.fix + 1
  }
  # get the rest of terms from Allb and Allu according to formula
  nfix <- 0
  nrandom <- 0
  next.random <- 1

  if(length(parsed.formula$labels) >= 1) {
    for(i in 1:length(parsed.formula$labels)) {
      anovalength <- length(anova)
      l <- parsed.formula$labels[i]
      if(parsed.formula$random[i] == 0) {
        # this term is fixed
        ncols <- dimX[nfix+1]
        anova[[anovalength+1]] <- Allb[,next.fix:(next.fix+ncols-1)]
        names(anova)[anovalength+1] <- l
        next.fix <- next.fix + ncols
        nfix <- nfix + 1
      }
      else {
        # this is a random term
        ncols <- dimZ[nrandom+1]
        anova[[anovalength+1]] <- Allu[,next.random:(next.random+ncols-1)]
        names(anova)[anovalength+1] <- l
        next.random <- next.random + ncols
        nrandom <- nrandom + 1
      }
      # make level for this term - this telles you the columns
      # of result represent which level of the term
      # skip making level for Spot and Label
      if( !(l %in% c("Spot", "Label")) ) {
        anovalength <- length(anova)
        anova[[anovalength+1]] <- makelevel(mamodel, l)
        names(anova)[anovalength+1] <- paste(l, "level", sep=".")
      }
    }
  }
  
  # get the flag for genes (if any)
  if(!is.null(madata$flag)) {
    flag <- rep(0, ngenes)
    for(i in 1:ngenes) {
      tmp <- madata$flag[(nreps*(i-1)+1):(nreps*i),]
      if(any(tmp))
        flag[i] <- 1
    }
    anova$flag <- flag
  }
    
  # put model into the object
  anova$model <- mamodel
  class(anova) <- "maanova"
  anova
}


###########################################
# function to make level for a term
###########################################
makelevel <- function(model, term)
{
  # local variables
  pf <- model$parsed.formula
  # output
  l <- NULL
  # find out the indices for non-reference sample
  idx.ref <- which(model$design$Sample==0)
  idx.noref <- setdiff(1:length(model$design$Sample), idx.ref)

  l <- as.vector(sort(unique(model$design[[term]][idx.noref])))

  # if it's a first order term
#  if( pf$order[pf$labels==term] == 1)
#    l <- as.vector(unique(model$design[[term]][idx.noref]))
#  else {
#    intterm <- which(pf$factors[,term]==1)
    # note that only twoway interaction are considered here
    # get the levels for two main terms
#    term1 <- pf$labels[intterm[1]]
#    level1 <- as.vector(unique(model$design[[term1]][idx.noref]))
#    term2 <- pf$labels[intterm[2]]
#    level2 <- as.vector(unique(model$design[[term2]][idx.noref]))
#    for(i in 1:length(level1)) 
#      for(j in 1:length(level2)) 
#        l[(i-1)*length(level2)+j] <- paste(level1[i],level2[j],sep=":")
#  }

  l
}


# summary function
print.maanova <- function(x,...)
{
}

    
