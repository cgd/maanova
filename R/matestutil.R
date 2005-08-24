######################################################################
#
# matestutil.R
#
# The utility functions for matest
#
# copyright (c) 2001-2004, Hao Wu and Gary A. Churchill, The Jackson Lab.
# Written Apr, 2004
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
#
######################################################################

###################################################
# locateTerm - locate the terms in given lables
###################################################
locateTerm <- function(labels, term)
{
  termidx <- NULL
  for(i in 1:length(term)) {
    idx <- which(term[i]==labels)
    if(length(idx)==0)
      warning(paste(term[i], "is not in formula, ignored."))
    else
      termidx <- c(termidx, idx)
  }
  termidx
}


###################################################
#
# caldf - the function to calculate the denominator's
# degree of freedom from model and the term to test
#
###################################################
caldf <- function(model, term)
{
  # local variables
  dimX <- model$dimX
  dimZ <- model$dimZ
  labels <- model$parsed.formula$labels
  random <- model$parsed.formula$random
  
  # check the nesting between the term to be tested and
  # all other terms in formula except itself
  confound <- NULL; nlevel <- NULL; tmpterm <- NULL
  # loop for all other terms
  for( i in 1:length(labels) ) {
    if( (labels[i] != term) & (random[i] != 0) ) {
      # if this is a random term
      tmp <- check.confounding(model, term, labels[i])
      confound <- c(confound, tmp$confounding)
      nlevel <- c(nlevel, tmp$dimX2)
      tmpterm <- c(tmpterm, labels[i])
    }
  }
  # note that confound shouldn't be null because
  # the model is mixed and the term is fixed
  
  if(all(confound!=2)) {
    # there's no nesting, df is the df for residuals
    result <- model$df[[length(model$df)]]
  }
  else {
    # this's nesting, calculate df from a lower level factor
    # sort the confound according to nlevels (number of columns)
    idx <- order(nlevel)
    confound <- confound[idx]
    nlevel <- nlevel[idx]
    tmpterm <- tmpterm[idx]
    idx <- min(which(confound==2))
    lowfactor <- tmpterm[idx]
    tmp <-  check.confounding(model, term, lowfactor)
    if(tmp$dimX1 > tmp$dimX2) {
      msg <- ("Model has problem:")
      msg <- paste(msg, "if", lowfactor, "is random,")
      msg <- paste(msg, term, "must be random as well.")
      stop(msg, call.=FALSE)
    }
    # df for lowfactor is the denominator's df for the test
    result <- model$df[[lowfactor]]
   }

  if(result < 0)
    stop("There's no enough degree of freedom to do test", call.=FALSE)
  
  result
}

        
############################################################
# check.confounding - to check the confounding for two terms
############################################################
check.confounding <- function(model, term1, term2)
{
  # local variables
  labels <- model$parsed.formula$labels
  random <- model$parsed.formula$random
  dimX <- model$dimX
  dimZ <- model$dimZ
  
  # get the design matrix for term 1
  idx.term1 <- locateTerm(labels, term1)
  if(random[idx.term1]==0) {
    # this term is fixed
    idx.term1 <- locateTerm(labels[random==0], term1)
    # starting column number in X matrix
    start <- 2
    if(any(model$design$Sample==0))
      start <- 3
    if(idx.term1 > 1)
      start <- start + sum(dimX[1:(idx.term1-1)])
    X.term1 <- model$X[,start:(start+dimX[idx.term1]-1)]
  }
  else {
    # this term is random
    idx.term1 <- locateTerm(labels[random==1], term1)
    start <- 1
    if(idx.term1>1)
      start <- start + sum(dimZ[1:(idx.term1-1)])
    X.term1 <- model$Z[,start:(start+dimZ[idx.term1]-1)]
  }

  # get the design matrix for term 2
  idx.term2 <- locateTerm(labels, term2)
  if(random[idx.term2]==0) {
    # this term is fixed
    idx.term2 <- locateTerm(labels[random==0], term2)
    # starting column number in X matrix
    start <- 2
    if(any(model$design$Sample==0))
      start <- 3
    if(idx.term2 > 1)
      start <- start + sum(dimX[1:(idx.term2-1)])
    X.term2 <- model$X[,start:(start+dimX[idx.term2]-1)]
  }
  else {
    # this term is random
    idx.term2 <- locateTerm(labels[random==1], term2)
    start <- 1
    if(idx.term2>1)
      start <- start + sum(dimZ[1:(idx.term2-1)])
    X.term2 <- model$Z[,start:(start+dimZ[idx.term2]-1)]
  }

  # combine two design matrices and check if it's singular
  X.combine <- cbind(X.term1, X.term2)
  r <- matrank(X.combine)
  if(r == max(dim(X.term1)[2], dim(X.term2)[2]))
    confounding <- 2 # nesting
  else if(r < dim(X.combine)[2]-1 )
    confounding <- 1 # partially confounding
  else
    confounding <- 0 # orthogonal

  if(is.matrix(X.term1)) dimX1 <- dim(X.term1)[2]
  else dimX1 <- 1
  if(is.matrix(X.term2)) dimX2 <- dim(X.term2)[2]
  else dimX2 <- 1

  result <- NULL
  result$confounding <- confounding
  result$X1 <- X.term1
  result$X2 <- X.term2
  result$dimX1 <- dimX1
  result$dimX2 <- dimX2

  result
}


##################################################
# function to calculate permutation P values
##################################################
calPval <- function(fstar, fobs, pool)
{
   
  # local variables
  ngenes <- dim(fstar)[1]
  n.perm <- dim(fstar)[3]
  nContrast <- dim(fstar)[2]
  npts <- length(fstar[,1,]) - ngenes

  # output
  F <- NULL
  F$Pvalpg <- matrix(0, ngenes, nContrast)
  F$Pvalmax <- matrix(0, ngenes, nContrast)

  # loop thru contrasts
  for(icon in 1:nContrast) {
    fstar.con <- fstar[,icon,]
    fobs.con <- fobs[,icon]
    
    # start calculation
    Fstarmax <- colmax(fstar.con)
    # calculate pvalmax (FWER adjusted p values)
    # and pvalpg (permutation p values)
    # changed by Hao Wu May 20, 2004 to remove
    # Fobs in Fstar
    if(pool) { # use pooled p-value
      # pvalmax
      for(i in 1:ngenes) {
        # nominal P value
        F$Pvalmax[i,icon] <- 1 - mean(fobs.con[i]>Fstarmax[-1])
      }
      # pooled P values
      # rank fstar
      Fstar.rank <- rank(fstar.con)
      # rank Fobs
      Fobs.rank <- rank(fobs.con)
      F$Pvalpg[,icon] <- 1 - (Fstar.rank[1:ngenes]-Fobs.rank)/npts
    }
    else { # non-pooled p-value
      for(i in 1:ngenes) {
        # nominal P value
        F$Pvalmax[i,icon] <- 1 - mean(fobs.con[i]>Fstarmax[-1])
        # permutation P values without pooling 
        F$Pvalpg[i,icon] <- 1 - mean(fobs.con[i]>fstar.con[i,-1])
      }
    }
  }
  # return
  F
}


########################################################
#
# Engine function to do test
#
########################################################
matest.engine <- function(anovaobj, term, mv, test.method, Contrast,
                          is.ftest, partC, verbose=FALSE)
{
  
  # local variables
  model <- anovaobj$model
  # parse the formula
  parsed.formula <- model$parsed.formula

  ################################################################
  # find the indices for these terms in information matrix (C)
  # relocate the term(s), e.g., the indices of term in all fixed terms
  fixed.term <- parsed.formula$labels[parsed.formula$random==0]
  termidx <- locateTerm(fixed.term, term)
  # find reference sample id
  refid <- which(model$design$Sample == 0)
  if(length(refid)==0)
    # no reference sample
    bstart <- 1
  else
    bstart <- 2
  colidx <- NULL
  for(i in 1:length(termidx)) {
    if(termidx[i] == 1)
      colstart <- bstart + 1
    else
      colstart <- sum(model$dimX[1:(termidx[i]-1)]) + bstart + 1
    colend <- sum(model$dimX[1:termidx[i]]) + bstart;
    colidx <- c(colidx, colstart:colend)
  }

  if(verbose==TRUE)
    cat("Start doing F test...\n")

  # get the estimates for the tested terms
  b <- NULL
  for(i in 1:length(termidx)) {
    tmpterm <- term[i]
    # get the estimates
    tmpb <- anovaobj[[tmpterm]]
    b <- cbind(b, tmpb)
  }
 
  # Finish making comparison matrix
  # calculate denominator's df from model and the term to test 
  dfFtest <- caldf(model, term[1])

  # setup some parameters
  if(is.ftest) { # this is F-test
    # df for numerator - fixed model only
    dfnu <- dim(Contrast)[1]
    # number of test is 1 for F test
    n.test <- 1
  }
  else {
    # do T-test
    # number of tests 
    n.test <- dim(Contrast)[1]
    dfnu <- 1 #df for T-test is 1
  }
  
                  
  ##################################
  # start to do tests
  ##################################
  # dfs for random terms
  idx.random <- c(which(model$parsed.formula$random==1), length(model$df))
  df.random <- model$df[idx.random]

  # initialize output object
  ngenes <- length(anovaobj$G)
  ftest <- NULL
  if(test.method[1] == 1)
    ftest$F1 <- matrix(0, ngenes, n.test)

  if(test.method[2] == 1)
    ftest$F2 <- matrix(0, ngenes, n.test)

  if(test.method[3] == 1)
    ftest$F3 <- matrix(0, ngenes, n.test)
  
  if(test.method[4] == 1)
    ftest$Fs <- matrix(0, ngenes, n.test)
  
  # for fixed model
  if(model$mixed == 0) {
    if(missing(partC)) { # calculate partC if not given
      X <- model$X
      C <- pinv(t(X) %*% X)
      partC <- C[colidx,colidx]
      ftest$partC <- partC
    }
    
    if(!is.ftest) {
      # do t-test
      ttest <- function(L) {
        tmp <- as.numeric(solve(t(L) %*% partC %*% L))
        fm <- b %*% L
        Fval <- fm^2 *tmp
        Fval
      }
      Fval <- apply(Contrast, 1, ttest)
    }
    else { # do F-test
      tmp <- solve(Contrast%*%partC%*%t(Contrast))
      fm <- b %*% t(Contrast)
      Fval <- apply(fm, 1, function(x) as.numeric((t(x) %*% tmp %*% x) / dfnu))
      Fval <- matrix(Fval,ncol=n.test)
    }

    #F1
    if(test.method[1] == 1) {
      ftest$F1 <- apply(Fval, 2, function(x) x/anovaobj$S2)
    }
    #F2
    if(test.method[2] == 1) {
      meanS2 <- mean(anovaobj$S2)
      tmps2 <- (meanS2+anovaobj$S2) / 2
      ftest$F2 <- apply(Fval, 2, function(x) x/tmps2)
    }
    # F3
    if(test.method[3] == 1) {
      meanS2 <- mean(anovaobj$S2)
      ftest$F3 <- apply(Fval, 2, function(x) x/meanS2)
    }
    #Fs
    if(test.method[4] == 1) {
      if(missing(mv)) {
        mv <- meanvarlog(df.random)
      }
      ftest$mv <- mv
      # JS shrunk variances
      JSs2 <- JSshrinker(anovaobj$S2, df.random, mv$meanlog, mv$varlog);
      ftest$Fs <- apply(Fval, 2, function(x) x/JSs2)
    }
  }
  
  if(model$mixed == 1) { # for mixed model
    X <- model$X
    XX <- t(X) %*% X
    Z <- model$Z
    XZ <- t(X) %*% Z
    ZZ <- t(Z) %*% Z
    dimZ <- model$dimZ
    k <- dim(XZ)[1]; m <- dim(XZ)[2]
    Im <- diag(1, m)
    # number of random variables (excluding residual)
    nrandom <- length(dimZ)

    if(test.method[3] == 1) {
      # make C matrix for F3 test. Note that C3 matrix is the same for all genes
      # so I pull it out of the gene loop
      # mean s2 for all random terms (including residual)
      means2 <- NULL
      for(i in 1:(nrandom+1))
        means2[i] <- mean(anovaobj$S2[,i])
      s0  <- means2[nrandom+1]
      d <- NULL
      for(i in 1:nrandom)
        d <- c(d, means2[i]*rep(1,dimZ[i]))
      D <- diag(d)
      V <- s0*Im + ZZ%*%D
      A <- rbind( cbind(XX, XZ%*%D), cbind(t(XZ), V) )
      A <- pinv(A)
      C3 <- s0 * rbind(A[1:k,], D%*%A[(k+1):(k+m),])
    }
    # finish build C3 matrix

    if(test.method[4] == 1) {
      # calcualate the shrinkage estimator for variance components.
      # This is needed for Fs test
      if(missing(mv)) {
        mv <- meanvarlog(df.random)
      }
      ftest$mv <- mv
      # JS shrunk variances
      JSs2 <- JSshrinker(anovaobj$S2, df.random, mv$meanlog, mv$varlog);
    }

    # loop for all genes
    for(i in 1:ngenes) {
      # make C matrices
      #C1
      if(test.method[1] == 1) {
        s2 <- anovaobj$S2[i,]
        s0 <- s2[nrandom+1]
        d <- NULL
        for(j in 1:nrandom)
          d <- c(d, s2[j]*rep(1,dimZ[j]))
        D <- diag(d)
        V <- s0*Im + ZZ%*%D
        A <- rbind( cbind(XX, XZ%*%D), cbind(t(XZ), V) )
        A <- pinv(A)
        C1 <- s0 * rbind(A[1:k,], D%*%A[(k+1):(k+m),])
      }
    
      # C2
      if(test.method[2] == 1) {
        s2 <- (anovaobj$S2[i,] + means2) / 2
        s0 <- s2[nrandom+1]
        d <- NULL
        for(j in 1:nrandom)
          d <- c(d, s2[j]*rep(1,dimZ[j]))
        D <- diag(d)
        V <- s0*Im + ZZ%*%D
        A <- rbind( cbind(XX, XZ%*%D), cbind(t(XZ), V) )
        A <- pinv(A)
        C2 <- s0 * rbind(A[1:k,], D%*%A[(k+1):(k+m),])
      }
    
      # Cs
      if(test.method[4] == 1) {
        s2 <- JSs2[i,]
        s0 <- s2[nrandom+1]
        d <- NULL
        for(j in 1:nrandom)
          d <- c(d, s2[j]*rep(1,dimZ[j]))
        D <- diag(d)
        V <- s0*Im + ZZ%*%D
        A <- rbind( cbind(XX, XZ%*%D), cbind(t(XZ), V) )
        A <- pinv(A)
        Cs <- s0 * rbind(A[1:k,], D%*%A[(k+1):(k+m),])
      }

      # make F or T test
      bgene <- b[i,]
      if(!is.ftest) { # this is a T-test
        ttest <- function(L) {
          tmp <- as.numeric(solve(t(L) %*% partC %*% L))
          fm <- bgene %*% L
          Fval <- fm^2 *tmp
          Fval
        }
      }
      else
        fm <- Contrast %*% bgene
      
      # F1
      if(test.method[1] == 1) {
#        browser()
        partC <- C1[colidx, colidx]
        if(!is.ftest)
          ftest$F1[i,] <- apply(Contrast, 1, ttest)
        else
          ftest$F1[i] <- as.numeric((t(fm) %*% solve(Contrast%*%partC%*%t(Contrast)) %*% fm)
                                    / dfnu)
      }
      # F2
      if(test.method[2] == 1) {
        partC <- C2[colidx, colidx]
        if(!is.ftest)
          ftest$F2[i,] <- apply(Contrast, 1, ttest)
        else
          ftest$F2[i] <- as.numeric((t(fm) %*% solve(Contrast%*%partC%*%t(Contrast)) %*% fm)
                                    / dfnu)
      }
      # F3
      if(test.method[3] == 1) {
        partC <- C3[colidx, colidx]
        if(!is.ftest)
          ftest$F3[i,] <- apply(Contrast, 1, ttest)
        else
          ftest$F3[i] <- as.numeric((t(fm) %*% solve(Contrast%*%partC%*%t(Contrast)) %*% fm)
                                    / dfnu)
      }
      
      # Fs
      if(test.method[4] == 1) {
        partC <- Cs[colidx, colidx]
        if(!is.ftest)
          ftest$Fs[i,] <- apply(Contrast, 1, ttest)
        else
          ftest$Fs[i] <- as.numeric((t(fm) %*% solve(Contrast%*%partC%*%t(Contrast)) %*% fm)
                                    / dfnu)
      }
      
    }
    # finish gene loop
  }
  
  ftest$dfFtest <- dfFtest
  ftest$dfnu <- dfnu
  ftest$Contrast <- Contrast
  ftest$partC <- partC
  # return
  ftest
}


######################################################################
#
# function to perform permutation F test
#
#######################################################################
matest.perm <- function(n.perm, FobsObj, data, model, term, Contrast, inits20, mv,
                        is.ftest, partC, MME.method, test.method,
                        shuffle.method, pool.pval)
{
  # local variables
  ngenes <- data$n.gene

  # number of contrasts
  if(is.ftest) # this is F-test
    nContrast <- 1
  else # this is T-test
    nContrast <- dim(Contrast)[1]
  
  # initialize result
  # the result is the number of permutation F values bigger than the observed
  # F value for each gene. Plus the maximum permutation F value for each permutaion
  Pval <- NULL
  if(test.method[1] == 1) {
    Pval$F1$Pmax <- array(0, c(ngenes, nContrast))
    Pval$F1$Pperm <- array(0, c(ngenes, nContrast))
  }
  if(test.method[2] == 1) {
    Pval$F2$Pmax <- array(0, c(ngenes, nContrast))
    Pval$F2$Pperm <- array(0, c(ngenes, nContrast))
  }
  if(test.method[3] == 1) {
    Pval$F3$Pmax <- array(0, c(ngenes, nContrast))
    Pval$F3$Pperm <- array(0, c(ngenes, nContrast))
  }
  if(test.method[4] == 1) {
    Pval$Fs$Pmax <- array(0, c(ngenes, nContrast))
    Pval$Fs$Pperm <- array(0, c(ngenes, nContrast))
  }

#  fstar <- NULL
#  if(test.method[1] == 1)
#    fstar$Fstar1 <- array(0, c(ngenes, nContrast, n.perm))
#  if(test.method[2] == 1)
#    fstar$Fstar2 <- array(0, c(ngenes, nContrast, n.perm))
#  if(test.method[3] == 1)
#    fstar$Fstar3 <- array(0, c(ngenes, nContrast, n.perm))
#  if(test.method[4] == 1)
#    fstar$Fstars <- array(0, c(ngenes, nContrast, n.perm))

  # residual shuffle is only for fixed model
  if( (model$mixed==1) & (shuffle.method=="resid") ) {
    warning("You can only do sample shuffling for mixed model test")
    shuffle.method <- "sample"
  }

  # prepare permutation 
  if(shuffle.method == "sample")  {
    # base used in shuffling
    shuffle.base <- NULL
    # if this is fixed model, shuffle tested term randomly
    # if this is mixed model, find nesting terms
    if(model$mixed) {
      nesting <- model$nesting
      # exclude the term itself
      #diag(nesting) <- NA
      # find the nesting terms
      idx.nesting <- which(nesting[term[1],]==1)
      # for multiple terms
      if(length(term) > 1) {
        for(iterm in 2:length(term)) {
          tmpidx <- which(nesting[term[iterm],]==1)
          idx.nesting <- intersect(idx.nesting, tmpidx)
        }
      }
      # only the random terms can be shuffle base
      idx.nesting <- intersect(idx.nesting,
                               which(model$parsed.formula$random==1))
      # if there's no base
      if(length(idx.nesting) == 0)
        shuffle.base <- NULL
      if(length(idx.nesting) > 1) {
        # if we have multiple base here, find the one on the lowest level
        # the term with the least number of columns should be the one
        #index for random terms in model
        idx.random <- which(model$parsed.formula$random==1)
        # index for nesting terms in random terms
        tmpidx <- NULL
        for(itmp in 1:length(idx.random))
          tmpidx <- c(tmpidx, which(idx.nesting[itmp]==idx.random))
        #tmpidx <- which(idx.nesting==idx.random)
        # number of columns for these nesting random terms
        dimZ.nesting <- model$dimZ[tmpidx]
        # which one is the smallest?
        idx.smallest <- which(dimZ.nesting == min(dimZ.nesting))
        # this one should be idx.nesting
        idx.nesting <- idx.nesting[idx.smallest]
      }

      # make shuffle base
      if(length(idx.nesting) == 0)
        shuffle.base <- NULL
      else {
        base.term <- model$parsed.formula$labels[idx.nesting]
        # make shuffle base
        designObj <- makeDesign(model$design)
        shuffle.base <- designObj[[base.term]]
      }
    }
  }      
        
  else if(shuffle.method=="resid") {
    # residual shuffling - for fixed model only
    # this is to shuffle the null model residuals
    # make a null model with the tested term excluded
    # terms in the null model. note that this is only for fixed model
    terms.null <- setdiff(model$parsed.formula$labels, term)
    if(length(terms.null) == 0) # null model has no terms
      nullformula <- "~1"
    else {
      nullformula <- paste(terms.null, collapse="+")
      nullformula <- paste("~", nullformula, sep="")
    }
    nullmodel <- makeModel(data, formula=as.formula(nullformula))
    # fit anova model on null model
    anova0 <- fitmaanova(data, nullmodel, verbose=FALSE)
    resid0 <- data$data - anova0$yhat
    nresid <- length(resid0)
    data.perm <- data
  }

  ##############################
  # permutation loop
  ##############################
  for(b in 1:n.perm) {
    if(round((b+1)/100) == (b+1)/100)
      cat(paste("Finish permutation # ", b+1, "\n"))
    
    design.perm <- model$design
    # for sample shuffling
    if(shuffle.method == "sample") {
      if(is.null(shuffle.base)) {
        # shuffle the design - this could be complicated
        design.perm <- shuffle.maanova(data, model, term)
      }
      else {
        # shuffle tested terms according to shuffle base
        # number of samples to be shuffled
        tmpSample <- shuffle.base[model$design$Sample != 0]
        tmpSample.unique <- unique(tmpSample)
        nsample <- length(tmpSample.unique)
        # index for non-reference sample
        idx.noref <- model$design$Sample!=0 
        # permutation index
        idx.perm <- sample(nsample)
        # shuffle the design and remake a model object
        design.perm <- model$design
        for(i in 1:nsample) {
          # index to be replaced
          # it should be non-ref sample
          idx <- which( (shuffle.base == tmpSample.unique[i]) & (idx.noref) )
          # index used to take value
          newidx <- which( (shuffle.base==tmpSample.unique[idx.perm[i]])
                          & (idx.noref) )
          newvalue <- model$design[[term]][newidx[1]]
          design.perm[[term]][idx] <- newvalue
        }
      }
             
      # remake model
      model.perm <- makeModel(data, design.perm, model$formula,
                              model$random, model$covariate)
      # fit anova model
      anovaobj.perm <- fitmaanova(data, model.perm, inits20=inits20,
                                method=MME.method, verbose=FALSE)
    }
    else if(shuffle.method=="resid") {
      # residual shuffling - this is for fixed model
      # global shuffle residual without replacement
      idx <- sample(nresid, replace=FALSE)
      # remake data
      data.perm$data <- anova0$yhat + resid0[idx]
      # fit anova model
      anovaobj.perm <- fitmaanova(data.perm, model, verbose=FALSE)
    }

    # start to do F-test
    # provide partC for fixed model test
    if(model$mixed == 0)
      # if this is fixed model, pass in partC
      # so we can save some calculation time
      ftest.perm <- matest.engine(anovaobj.perm, term, mv, test.method,
                                Contrast, is.ftest, partC, verbose=FALSE)
    else
      ftest.perm <- matest.engine(anovaobj.perm, term, mv, test.method,
                                  Contrast, is.ftest, verbose=FALSE)

    # update the result
    ffields <- c("F1","F2","F3","Fs")
    for(icon in 1:nContrast) { # loop for multiple contrasts
      for(i in 1:4) { # loop for different F test methods
        if(test.method[i] == 1) { # if this F test is requested
          fobs <- FobsObj[[ffields[i]]]$Fobs[,icon]
          fstar <- ftest.perm[[ffields[i]]][,icon]
          fstar.max <- max(fstar)
          # For Fstarmax
          if(test.method[i] == 1)
            Pval[[ffields[i]]]$Pmax[,icon] <-
              Pval[[ffields[i]]]$Pmax[,icon] + (fstar.max>fobs)
        
          # for P value count
          if(pool.pval) { # use pooled p value
            Fobs.rank <- rank(fobs)
            Fstar.rank <- rank(cbind(fobs, fstar))
            Pval[[ffields[i]]]$Pperm[,icon] <- Pval[[ffields[i]]]$Pperm[,icon] +
              ngenes - (Fstar.rank[1:ngenes]-Fobs.rank)
          }
          else { # not pool
            tmp <- fstar > fobs
            Pval[[ffields[i]]]$Pperm[,icon] <- Pval[[ffields[i]]]$Pperm[,icon] + tmp
          }
        }
      }
    }
  }
  
  # return
  Pval
}
  

################################################
# checkContrast - to check if the given contrast
# matrix is valid or not
################################################
checkContrast <- function(model, term, Contrast)
{
  # local variables
  fullC <- rep(0, dim(model$X)[2])
  parsed.formula <- model$parsed.formula
  fixed.term <- parsed.formula$labels[parsed.formula$random==0]
  termidx <- locateTerm(fixed.term, term)
  nContrast <- dim(Contrast)[1] # number of contrast
  rankX <- matrank(model$X)
  # starting column number for terms
  # 1 if no reference sample
  # 2 if has reference sample 
  refid <- which(model$design$Sample == 0)
  if(length(refid)==0)
    # no reference sample
    bstart <- 1
  else
    bstart <- 2

  # loop for all terms
  colidx <- NULL
  for(j in 1:length(termidx)) {
    if(termidx[j] == 1)
      colstart <- bstart + 1
    else
      colstart <- sum(model$dimX[1:(termidx[j]-1)]) + bstart + 1
    colend <- colstart + model$dimX[termidx[j]] - 1;
    colidx <- c(colidx, colstart:colend)
  }
  
  # loop for contrasts
  for(i in 1:nContrast) {
    # current contrast
    thisContrast <- Contrast[i,]
    fullC[colidx] <- thisContrast
    # see if fullC is a linear combination of all rows in X
    rank.tmp <- matrank(rbind(fullC, model$X))
    if(rank.tmp > rankX)
      stop(paste("The number", i, "test is not estimable"), call.=FALSE)
  }

}


#################################################
# makeContrast - make contrast matrix for F test
# if it's not given
#################################################
makeContrast <- function(model, term)
{
  # local variables
  parsed.formula <- model$parsed.formula
  
  # starting column number for terms
  # 1 if no reference sample
  # 2 if has reference sample 
  refid <- which(model$design$Sample == 0)
  if(length(refid)==0)
    # no reference sample
    bstart <- 1
  else
    bstart <- 2

  # get the column number in design matrix X for tested term
  fixed.term <- parsed.formula$labels[parsed.formula$random==0]
  termidx <- locateTerm(fixed.term, term)

  # loop for all terms
  colidx <- NULL
  for(j in 1:length(termidx)) {
    if(termidx[j] == 1)
      colstart <- bstart + 1
    else
      colstart <- sum(model$dimX[1:(termidx[j]-1)]) + bstart + 1
    colend <- colstart + model$dimX[termidx[j]] - 1;
    colidx <- c(colidx, colstart:colend)
  }
   
  # make X0 (the design matrix for null model)
  # and X.term (the design matrix for tested terms)
  ncol.total <- dim(model$X)[2]
  col.keep <- setdiff(1:ncol.total, colidx)
  X0 <- model$X[,col.keep]
  X.term <- model$X[,colidx]
  if(matrank(X0) == matrank(model$X))
    stop("No degree of freedom to do the test", call.=FALSE)
  
  # make X.new, which spans the same space as X
  # but have redundent columns removed for the tested terms
  X.new <- X0
  idx.keptcol <- NULL
  for(i in 1:length(colidx)) {
    X.tmp <- cbind(X.new, model$X[,colidx[i]])
    if(matrank(X.tmp) > matrank(X.new)) {
      # add this column
      X.new <- X.tmp
      col.keep <- c(col.keep, colidx[i])
    }
  }
  col.keep <- sort(col.keep)
  X.new <- model$X[,col.keep]
  XX.new <- t(X.new) %*% X.new
  invXX.new <- pinv(XX.new)
  # restore to full scale
  ncol.X <- dim(model$X)[2]
  ncol.X0 <- dim(X0)[2]
  invXX <- matrix(0, ncol.X, ncol.X)
  invXX[col.keep, col.keep] <- invXX.new

  # use invXX %*% XX to generate the contrast matrix
  XX <- t(model$X) %*% model$X
  tmp <- invXX %*% XX
  tmpL <- tmp[colidx, colidx]
  # erase the calculation residuals
  tmpL[abs(tmpL) < 1e-10] <- 0
  idx.nonzero <- which(apply(tmpL,1, any))

  # return as a matrix
  matrix(tmpL[idx.nonzero,], nrow=length(idx.nonzero))
  
}


###############################################################
#
# function to shuffle the design for microarray experiment.
# This could be pretty complicated. Basically, Array, Dye,
# And shuffle base need to be kept intact. The corelation
# structure cannot be broken.
#  
# This function is not COMPLETE! It is only working for non-shuffle base cases,
# e.g., shuffle arrays. For shuffle bases in mixed models,
# I still need to think about it.
# 
###############################################################
shuffle.maanova <- function(data, model, term)
{
  # is there Dye effect in the model?
  has.dye <- "Dye" %in% model$parsed.formula$labels
  # is there Array effect in the model?
  has.array <- "Array" %in% model$parsed.formula$labels
  # is there reference sample in the model?
  has.ref <- any(model$design$Sample==0)

  # some local variables
  ndye <- data$n.dye
  idx.dye <- 1:ndye
  narray <- data$n.array
  idx.array <- 1:narray
  # original index - make it a matrix of ndye x narray
  idx.obs <- matrix(1:(ndye*narray), nrow=ndye)

  # result
  design.perm <- model$design

  # make the shuffled index for the tested term
  if(!has.ref) { # if there's no reference samples - this is easier
    if(has.array) { # has Array effect in the model
      # shuffled unit is array
      idx.array <- sample(narray)
      idx.obs <- idx.obs[,idx.array]
      if(!has.dye) {
        # no Dye effect, dye within each array can be shuffled
        for(i in 1:narray) { # shuffle dye within each array
          idx.dye <- sample(ndye)
          idx.obs[,i] <- idx.obs[idx.dye,i]
        }
      }
      idx.perm <- as.vector(idx.obs)
    }

    else { # NO array effect, the tie among arrays can be broken now
      if(has.dye) { 
        # has dye effect, samples on a specific dye need to
        # be on the same dye after shuffling
        # so we are shuffling each row seperately for idx.obs matrix
        for(i in 1:ndye) {
          idx.array <- sample(narray)
          idx.obs[i,] <- idx.obs[i,idx.array]
        }
        idx.perm <- as.vector(idx.obs)
      }
      else { # has no dye effect, shuffle everything freely
        idx.perm <- as.vector(idx.obs[sample(ndye*narray)])
      }
    }
  }

  else  { # there IS reference samples - this is troubler
    # rearrange sample ID to a matrix of ndye x nsample
    sample.mtx <- matrix(model$design$Sample, nrow=ndye)
    if(has.array) {
      # has array effect, must divide the arrays into exchangable groups
      # The arrays in each group are exchangable
      # There'll be ndye+1 groups, representing the arrays with reference sample
      # at ith (i=1:ndye) dye and arrays without reference sample
      # (in the cases of reference and non-reference experiments put together).
      groups <- makeShuffleGroup(sample.mtx, ndye, narray)
      # shuffle arrays within each group
      for(i in 1:length(groups)) {
        narray.group <- length(groups[[i]])
        idx.group.perm <- groups[[i]][sample(narray.group)]
        # shuffle the columns of idx.obs
        idx.obs[,groups[[i]]] <- idx.obs[,idx.group.perm]
      }
      if(!has.dye) {
        # has array, no dye, shuffle the arrays within each group,
        # and shuffle (non-reference) dyes within each array
        for(i in 1:narray) {
          # non-refernce samples on this array
          idx.noref <- which(sample.mtx[,i]!=0)
          n.noref <- length(idx.noref)
          # shuffle
          if(n.noref > 1) {
            idx.dye.perm <- sample(n.noref)
            idx.obs[idx.noref,i] <- idx.obs[idx.noref[idx.dye.perm],i]
          }
        }
      }
      idx.perm <- as.vector(idx.obs)
    }
    
    else{ # no array effect
      if(has.dye) {
        # no array has dye, freely shuffle the non-reference samples
        # in each row of idx.obs
        for(i in 1:ndye) {
          idx.noref <- which(sample.mtx[i,]!=0)
          n.noref <- length(idx.noref)
          # shuffle this row in idx.obs
          idx.obs[i,idx.noref] <- idx.obs[i, idx.noref[sample(n.noref)]]
        }
        idx.perm <- as.vector(idx.obs)
      }
      else{ # no array no dye, shuffle non-reference samples freely
        # non-reference samples
        idx.noref <- which(model$design$Sample!=0)
        n.noref <- length(idx.noref)
        idx.perm <- 1:(ndye*narray)
        idx.perm[idx.noref] <- idx.noref[sample(n.noref)]
      }
    }
    
  }

  # shuffle the tested term
  for(iterm in 1:length(term)) {
    design.perm[[term[iterm]]] <- model$design[[term[iterm]]][idx.perm]
  }

  # return
  design.perm
}

##########################################################
# function to divide the arrays into exchangable groups
# this is called by shuffle.maanova
##########################################################
makeShuffleGroup <- function(sample.mtx, ndye, narray) {
  idx.array <- 1:narray
  groups <- NULL
  for(i in 1:ndye) {
    groups[[i]] <- which(sample.mtx[i,]==0)
    # the arrays left
    idx.array <- setdiff(idx.array, groups[[i]])
  }
  # the last one is for the arrays without reference sample
  if(length(idx.array) > 0)
    groups[[ndye+1]] <- idx.array

  groups
}

