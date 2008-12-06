#####################################################################
#
# matest.R
# copyright (c) 2001-2004, Hao Wu, Hyuna Yang and Gary A. Churchill, 
# The Jackson Lab.
# Written Apr, 2004
# Modified Apr, 2006
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# This is the function to do F permutation test 
#
######################################################################

matest <- function(data, anovaobj, term, Contrast, n.perm=1000, nnodes=1,
                   critical=.9, test.type=c("ttest","ftest"),
                   shuffle.method=c("sample", "resid"),
                   MME.method=c("REML","noest","ML"),
                   test.method=c(1,1),pval.pool=TRUE, verbose=TRUE){
  if( class(data) != "madata" )
    stop("data is not an object of class madata.")
  if( class(anovaobj) != "maanova" )
    stop("anovaobj is not an object of class maanova.")
  if( critical <= 0 || critical >1 ) 
    stop(" Critical value is between 0 to 1: Recommended value =.9")
  model = anovaobj$model 

  #if(model$mixed != 0 & test.method[3]==1){
  #  cat("Fss is not availble in mixed effect model\n")
  #  test.method[3]=0
  #}
  if( test.method[1]+test.method[2]==0 ){
    stop("F1 or Fs test should be included in the test statistics\n")
  }
  
  # get the methods
  shuffle.method <- match.arg(shuffle.method)
  MME.method <- match.arg(MME.method)
  
  # local variables
  nreps <- data$n.rep
  ngenes <- data$n.gene
  narrays <- data$n.array

  ########## initialize cluster ############
  if( (n.perm>1) & (nnodes>1) ) {
    if(!require(snow))
      stop("You have to install SNOW package to use cluster")
    # make cluster object
    cl <- makeMPIcluster(nnodes)
    # correct the possible correlation problem
    clusterApply(cl, runif(length(cl), max=1000000000), set.seed)
  }

  ################################################################
  # if any tested term is random, issue a warning and rebuild model
  ################################################################
  # take the parsed formula from model
  parsed.formula <- model$parsed.formula
  # identify the terms
  termidx <- locateTerm(parsed.formula$labels, term)

  if( is.null(termidx) )
    stop("The term(s) to be tested is not in formula.")
  
  if(any(parsed.formula$random[termidx] == 1)) {
    warning("Testing random terms. Random term will be tested as fixed terms.")
    # rebuild random formula and model object
    # random terms
    random.term <- parsed.formula$labels[parsed.formula$random==1]
    # random term to be tested
    random.term.test <- parsed.formula$labels[termidx]
    random.term.test <- random.term.test[parsed.formula$random[termidx]==1]
    # make the random term to be tested fixed,
    # remake the formula for random
    random.term <- setdiff(random.term, random.term.test)
    if(length(random.term) == 0) {
      random <- ~1
      stop("Cannot test the only random term in the model.")
    }
    else{
      random <- paste(random.term, collapse="+")
      random <- paste("~", random, sep="")
    }
    # remake covariate
    cov.term <- parsed.formula$labels[parsed.formula$covariate==1]
    if(length(cov.term) == 0)
      cov <- ~1
    else {
      cov <- paste(cov.term, collapse="+")
      cov <- paste("~", cov, sep="")
    }
    # remake model
    
    anovaobj <- fitmaanova(data, formula=model$formula, 
      random=as.formula(random),covariate=as.formula(cov),
      method=MME.method, verbose=FALSE, subCol=anovaobj$subCol)
  }
  #########################################################
  # for contrast matrix - make one if not given,
  # check the validity if given
  #########################################################
  # for backward compatibility, if user provide Contrast matrix
  # this will be a T-test by default, unless user specify the test type
  if(missing(Contrast)){# no Contrast matrix, make it
    Contrast <- makeContrast(model, term)
    nContrast <- 1
    # this must be a F-test
    if(missing(test.type))
      is.ftest <- TRUE
    else{# test.type is given
      test.type <- match.arg(test.type)
      if(test.type=="ttest")
        is.ftest <- FALSE
      else
        is.ftest <- TRUE
      if(test.type=="ttest") # cannot be T-test
        stop("You must specify a Contrast matrix for doing T-test")
    }
  }
  else{ # given contrast matrix.
    # check if the contrast matrix is valid
    checkContrast(model, term, Contrast)
    # test type is T-test by default for backward compatibility
    test.type <- match.arg(test.type)
    if(test.type=="ftest") {
      is.ftest <- TRUE
      nContrast <- 1
    }
    else { # this is T-test, but will return F values 
      is.ftest <- FALSE
      # number of tests requested
      nContrast <- dim(Contrast)[1] # number of t-tests (comparisons)
    }
  }
  
  # do F test on the observed data
  if(verbose) cat("Doing F-test on observed data ...\n")
  ftest.obs <- matest.engine(anovaobj, term, test.method=test.method, 
           Contrast=Contrast, is.ftest=is.ftest, verbose=FALSE)
  # get the results
  mv <- ftest.obs$mv
  dfnu <- ftest.obs$dfnu; dfFtest <- ftest.obs$dfFtest
  partC <- ftest.obs$partC
  #mean_est = ftest.obs$mean_est
  #tau_est = ftest.obs$tau_est
  # initialize output object
  ftest <- NULL
  # general info in the output object
  ftest$critical = critical 
  ftest$model <- model
  ftest$term <- term
  ftest$dfde <- ftest.obs$dfFtest
  ftest$dfnu <- dfnu
  ftest$obsAnova <- anovaobj
  ftest$Contrast <- Contrast
  #ftest$mean_est = ftest.obs$mean_est
  #ftest$tau_est = ftest.obs$tau_est
  ftest$probeid = data$probeid
  if(!is.ftest)
    class(ftest) <- c("matest", "ttest")
  else
    class(ftest) <- c("matest", "ftest")

  # calculate P values
  # F1
  if(test.method[1] == 1) {
    ftest$F1 <- NULL
    ftest$F1$Fobs <- ftest.obs$F1
    ftest$F1$Ptab <- 1 - pf(ftest$F1$Fobs, dfnu, dfFtest)
    if(n.perm > 1) {
      ftest$F1$Pvalperm <- array(0, c(ngenes, nContrast))
      ftest$F1$Pvalmax <- array(0, c(ngenes, nContrast))
    }
  }
  #Fs
  if(test.method[2] == 1) {
    ftest$Fs <- NULL
    ftest$Fs$Fobs <- ftest.obs$Fs
    ftest$Fs$Ptab <- 1 - pf(ftest$Fs$Fobs, dfnu, dfFtest)
    if(n.perm > 1) {
      ftest$Fs$Pvalperm <- array(0, c(ngenes, nContrast))
      ftest$Fs$Pvalmax <- array(0, c(ngenes, nContrast))
    }
  }

  #Fss
  #if(test.method[3] == 1) {
  #  ftest$Fss$Fobs <- ftest.obs$Fss
  #  ftest$Fss$Ptab <- 1 - pf(ftest$Fss$Fobs, dfnu, dfFtest)
  #  if(n.perm > 1) {
  #    ftest$Fss$Pvalperm <- array(0, c(ngenes, nContrast))
  #    ftest$Fss$Pvalmax <- array(0, c(ngenes, nContrast))
  #  }
  #}


  # return if no permutation test
  if(n.perm == 1) {
    warning("You are not doing permutation test. Only observed values are calculated.")
    return(ftest)
  }

  ########################################
  # start to do permutation test
  ########################################
  if(verbose)
    cat("Doing permutation. This may take a long time ... \n")
  sdata=data;  sS2=anovaobj$S2;

  if(critical <1){
    if(nContrast==1){
      hsidx = ftest$F1$Fobs <= qf(critical,ftest$dfnu ,ftest$dfde)
    }
    else{
      tmp=as.matrix(apply(ftest$F1$Fobs, 1, min))
      hsidx= tmp <= qf(critical,ftest$dfnu ,ftest$dfde)
    }
    if(ncol(sS2)==1) sS2 = matrix(sS2[hsidx==TRUE,], ncol=1)
    else sS2=sS2[hsidx==TRUE,]
    sdata$n.gene = sum(hsidx)
    sdata$data = sdata$data[hsidx==TRUE,];
    sdata$metarow = sdata$metarow[hsidx==TRUE]
    sdata$metacol = sdata$metacol[hsidx==TRUE]
    sdata$probeid = sdata$probeid[hsidx==TRUE]
    sdata$colmeans = apply(sdata$data, 2, mean)
  }
  # permutation - use MPI cluster if specified
  if(nnodes > 1) { # use MPI cluster
    # use cluster call to do permutation
    # calculate the number of permutation needed in each node
    nperm.cluster <- rep(floor((n.perm-1)/nnodes), nnodes)
    # maybe some leftovers
    leftover <- n.perm - 1 - sum(nperm.cluster)
    if(leftover > 0)
      nperm.cluster[1:leftover] <- nperm.cluster[1:leftover] + 1
    
    # load library on all nodes
    clusterEvalQ(cl, library(maanova))
    # use clusterApply to run permutation on all nodes
    # note that the only different parameter passed to function
    # is ftest.mixed.perm. So in ftest.mixed.perm I put nperm 
    # as the first argument so I can use clusterApply
    cat(paste("Doing permutation on", nnodes, "cluster nodes ... \n"))
   
    pstar.nodes <- clusterApply(cl, nperm.cluster, matest.perm,ftest,sdata,
        model, term,Contrast, mv, is.ftest, partC, 
    	MME.method, test.method, shuffle.method, pval.pool, ngenes)
    ffields <- c("F1","Fs")
    for(i in 1:nnodes) {
      if(nperm.cluster[i] > 0) {
        pstar <- pstar.nodes[[i]]
        for(itest in 1:2) {
          if(test.method[itest] == 1){
            ftest[[ffields[itest]]]$Pvalperm <-
              ftest[[ffields[itest]]]$Pvalperm + pstar[[ffields[itest]]]$Pperm
            ftest[[ffields[itest]]]$Pvalmax <-
              ftest[[ffields[itest]]]$Pvalmax + pstar[[ffields[itest]]]$Pmax
          }
        }
      }
    }
    # clear cluster results to save some memory
    rm(pstar.nodes);

    # stop the cluster
    stopCluster(cl)
  }
  else { # no cluster, do it on single node
    pstar <- matest.perm(n.perm, ftest, sdata, model, term,
      Contrast, mv, is.ftest, partC, MME.method, test.method,
      shuffle.method,pval.pool,ngenes)
    
    # update the pvalues
    ffields <- c("F1","Fs")
    for(itest in 1:2) {
      if(test.method[itest] == 1) {
        ftest[[ffields[itest]]]$Pvalperm <-
          ftest[[ffields[itest]]]$Pvalperm + pstar[[ffields[itest]]]$Pperm
        ftest[[ffields[itest]]]$Pvalmax <-
          ftest[[ffields[itest]]]$Pvalmax + pstar[[ffields[itest]]]$Pmax
      }
    }
  }    
  #finish permutation loop
  
  # calculate the pvalues. Note that the current object contains the "counts"
  ffields <- c("F1","Fs")
  for(itest in 1:2) {
    if(test.method[itest] == 1) {
      # Pvalperm
      if(pval.pool){
        ftest[[ffields[itest]]]$Pvalperm <-
          ftest[[ffields[itest]]]$Pvalperm / (n.perm*sdata$n.gene)
      }
      else
        ftest[[ffields[itest]]]$Pvalperm <-
          ftest[[ffields[itest]]]$Pvalperm / n.perm
      # for Pvalmax
      ftest[[ffields[itest]]]$Pvalmax <-
        ftest[[ffields[itest]]]$Pvalmax / n.perm
    }
  }
  #if(test.method[3] == 1) {
  #  itest = 3
  #  # Pvalperm
  #  if(pval.pool){
  #    ftest[[ffields[itest]]]$Pvalperm <-
  #      ftest[[ffields[itest]]]$Pvalperm / (n.perm*data$n.gene)
  #  }
  #  else
  #    ftest[[ffields[itest]]]$Pvalperm <-
  #      ftest[[ffields[itest]]]$Pvalperm / n.perm
  #  # for Pvalmax
  #  ftest[[ffields[itest]]]$Pvalmax <-
  #    ftest[[ffields[itest]]]$Pvalmax / n.perm
  #}

  # add some other info to the return object
  ftest$n.perm <- n.perm
  ftest$shuffle <- shuffle.method
  ftest$pval.pool <- pval.pool
  
  # return
  ftest
}

