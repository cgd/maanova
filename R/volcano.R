######################################################################
#
# volcano.R
#
# copyright (c) 2001, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written Nov, 2001
# Modified Apr, 2004
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
#
######################################################################


volcano <-
  function(matestobj, threshold=c(0.001,0.05),
           method=c("unadj","unadj"), title="Volcano Plot",
           highlight.flag=TRUE, onScreen=TRUE)
{
  if(class(matestobj)[1] != "matest")
    stop("The first input variable is not an object of class matest.")

  # must have F1
  if( !("F1" %in% names(matestobj))  )
    stop("No Fs test result in input object. Cannot do volcano plot")

  if(class(matestobj)[2] == "ttest")
    # this is a T-test result
    result <- volcano.ttest(matestobj, threshold, method, title,
                            highlight.flag, onScreen)
  else
    # this is a F test volcano plot
    result <- volcano.ftest(matestobj, threshold, method, title,
                            highlight.flag)
}


######################################################################
# function to do f-test volcano plot
######################################################################
volcano.ftest <- function(matestobj, threshold, method, title,
                          highlight.flag){
  anovaobj = matestobj$obsAnova
  probeid = matestobj$probeid

  if(length(method) == 1) method <- rep(method, 2)
  
  # local variables
  th.f1 <- threshold[1]
  th.fs <- threshold[2]
  #th.fss <- threshold[3]
  ngenes <- length(anovaobj$G)

  # get P values for all F tests according to method
  # for F1
  p1 <- getPval.volcano(matestobj, method, 1)
  idx1 <- p1 < th.f1 ; 
  # for Fs
  ps <- getPval.volcano(matestobj, method, 2)
  idxs <- ps < th.fs ; 
  # for Fss
  #pss <- getPval.volcano(matestobj, method, 3)
  #idxss <- pss < th.fss; ; 
  
  # calculate the x axis value
  # the x-axis value should be the numerator of the test, e.g.
  # (Lb)'*pinv((L*pinv(X'X)*L')*(Lb)

  xvalue <- calVolcanoXval(matestobj)
  # y-axis value
  yvalue <- -log10(p1)

  # plot the figure
  plot(xvalue, yvalue, xlab="Log2(FoldChange)", ylab="-log10(Pvalue)",
       col="blue", pch=4, cex=0.5, main=title)

  # draw a reference line based on the threshold of f1
  abline(h=-log10(th.f1))
  
  # replot the significant genes from Fs (if there) in red
  if("Fs" %in% names(matestobj))
    points(xvalue[idxs], yvalue[idxs], col="red", pch=4, cex=0.5)
 
  # plot the significant genes from Fss (if there) in orange
  #if("Fss" %in% names(matestobj))
  #  points(xvalue[idxss], yvalue[idxss], col="orange", pch=4, cex=0.5)

  # circle the flaged genes (if any)
  if(highlight.flag) {
    if(!is.null(anovaobj$flag)) {
      idx.flag <- which(anovaobj$flag==1)
      points(xvalue[idx.flag], yvalue[idx.flag])
    }
  }

  # find the significant genes from all three F tests
  result <- NULL
  idx.F1 <- which(idx1); idx <- idx.F1
  names(idx.F1) = probeid[idx.F1]
  result$idx.F1 =idx.F1

  if("Fs" %in% names(matestobj)) {
    idx.Fs <- which(idxs); idx <- intersect(idx, idx.Fs)
    names(idx.Fs) = probeid[idx.Fs]
    result$idx.Fs <- idx.Fs 
  }
  #if("Fss" %in% names(matestobj)) {
  #  idx.Fss <- which(idxss); idx <- intersect(idx, idx.Fss)
  #  names(idx.Fss) = probeid[idx.Fss]
  #  result$idx.Fss <- idx.Fss 
  #}
  idx.all = idx; names(idx.all) = probeid[idx.all]; result$idx.all <- idx.all
  result
}



######################################################################
# function to do T-test volcano plot
# this will generate multiple plots
######################################################################
volcano.ttest <- function(matestobj, threshold, method, title,
                          highlight.flag, onScreen)
{
  anovaobj = matestobj$obsAnova
  probeid = matestobj$probeid
  
  if(length(method) == 1)
    method <- rep(method, 2)
  
  # local variables
  th.f1 <- threshold[1]
  th.fs <- threshold[2]
  #th.fss <- threshold[3]
  ngenes <- length(anovaobj$G)
  
  # get P values for all F tests according to method
  # for F1
  p1 <- getPval.volcano(matestobj, method, 1)
  idx1 <- p1 < th.f1
  # for Fs
  ps <- getPval.volcano(matestobj, method, 2)
  idxs <- ps < th.fs
  # for Fss
  #pss <- getPval.volcano(matestobj, method, 3)
  #idxss <- pss < th.fss
 
  #############################
  # calculate the x axis value
  #############################
  # the term(s) tested
  diff.terms <- matestobj$term
  Contrast <- matestobj$Contrast
  # number of plots to be generated
  nplots <- dim(Contrast)[1]

  # init result
  result <- NULL

  # calculate x-axis values
  xvalue.all <- calVolcanoXval(matestobj)

  # loop for contrasts
  levels <- NULL
  for(i in 1:length(diff.terms)) {
    # level names
    levels <- c(levels, anovaobj[[paste(diff.terms[i], ".level", sep="")]])
  }
  for(icon in 1:nplots) {
    # xvalue for this contrast
    xvalue <- xvalue.all[,icon]
    
    ######### make x-axis labels
    xlabel = levels[Contrast[icon, ] !=0]
    xlabel <- paste(paste(xlabel, sep="*"), collapse="-")
    xlabel <- paste('Log2(FoldChange) of ', xlabel, sep="", collapse=" ")

    # figure title
    title <- paste("comparison", icon)

    # start to plot
    if(onScreen) {
      # open a window on screen
      dev.new()
    }
    # y-axis value
    yvalue <- -log10(p1[,icon])
    # plot the figure
    plot(xvalue, yvalue, xlab=xlabel, ylab="-log10(Pvalue)",
         col="blue", pch=4, cex=0.5, main=title)
  
    # draw a reference line based on the threshold of f1
    abline(h = -log10(th.f1))
    idx1 <- p1[,icon] < th.f1
 
    # replot the significant genes from Fs (if there) in green
    if("Fs" %in% names(matestobj)) {
      idxs <- ps[,icon] < th.fs
      points(xvalue[idxs], yvalue[idxs], col="red", pch=4, cex=0.5)
    }
    # plot the significant genes from Fss (if there) in orange
    #if("Fss" %in% names(matestobj)) {
    #  idxss <- pss[,icon] < th.fss
    #  points(xvalue[idxss], yvalue[idxss], col="orange", pch=4, cex=0.5)
    #}

    # draw the reference line based on F3 test (if there)
    # note that this is alway vshape
#    vshape <- 0
#    if("F3" %in% names(matestobj)) {
#      idx3 <- p3[,icon] < th.f3
#      if(sum(idx3) != 0) {
        # find the location of the reference line
#        l <- min(abs(xvalue[idx3]))
#        if(vshape == 1) {
#          abline(v=l)
#          abline(v=-l)
#        }
#        else {
#          abline(v=l)
#        }
#      }
#      else
#        warning("There is no significant genes from F3 test.")
#    }
  
    # circle the flaged genes (if any)
    if(highlight.flag) {
      if(!is.null(anovaobj$flag)) {
        idx.flag <- which(anovaobj$flag==1)
        points(xvalue[idx.flag], yvalue[idx.flag])
      }
    }
  
    # find the significant genes from all three F tests

    result.tmp <- NULL
    idx.F1 <- which(idx1); idx <- idx.F1
    names(idx.F1) = probeid[idx.F1]
    result.tmp$idx.F1 <- idx.F1

    #if("Fss" %in% names(matestobj)) {
    #  idx.Fss <- which(idxss); idx <- intersect(idx, idx.Fss)
    # names(idx.Fss) = probeid[idx.Fss]
    # result.tmp$idx.Fss <- idx.Fss 
    #}
    if("Fs" %in% names(matestobj)) {
      idx.Fs <- which(idxs); idx <- intersect(idx, idx.Fs)
      names(idx.Fs) = probeid[idx.Fs]
      result.tmp$idx.Fs <- idx.Fs 
    }
    idx.all <- idx
    names(idx.all) = probeid[idx.all]; result.tmp$idx.all <- idx.all
    result.name <- paste("comparison", icon, sep="")
    result[[result.name]] <- result.tmp
  }
  result
}



###########################################
# function to get P values from matest object
###########################################
getPval.volcano <- function(matestobj, method, idx)
{
  # field name
  if(idx == 2)
    whichF <- "Fs"
  if(idx == 1)
    whichF <- "F1"
  #if(idx == 3)
  #  whichF <- "Fss"

  if( !(whichF %in% names(matestobj)) )
    return(NULL)

  # get the P values
  Fobj <- matestobj[[whichF]]
  method <- method[idx]
  
  if(method == "unadj")
    p <- Fobj$Ptab
  else if(method == "nominal")
    p <- Fobj$Pvalperm
  else if(method[1] == "fwer")
    p <- Fobj$Pvalmax
  else if(method[1] == "fdr")
    p <- Fobj$adjPtab
  else if(method[1] == "fdrperm")
    p <- Fobj$adjPvalperm
  else
    stop(paste("Unrecognized method for F test P value,", method))

  if(is.null(p))
    stop(paste(method, "P value is not available for", whichF))

  # change zeros to 1e-17
  p[p==0] <- 1e-17
  # return
  p
}


##############################################
# function to calculate the x-axis value
# on the volcano plot
# basically it's from Lb
# If L is one row,it's just Lb
# IF L has multiple rows, use sqrt((Lb)' * Lb))
##############################################
calVolcanoXval <- function(matestobj)
{
  model <- matestobj$model
  anovaobj = matestobj$obsAnova

  # get the estimates from ANOVA object on observed data
  term <- matestobj$term
  parsed.formula <- model$parsed.formula
  fixed.term <- parsed.formula$labels[parsed.formula$random==0]
  termidx <- locateTerm(fixed.term, term)

  b <- NULL
  for(i in 1:length(termidx)) {
    tmpterm <- term[i]
    # get the estimates
    tmpb <- anovaobj[[tmpterm]]
    b <- cbind(b, tmpb)
  }

  # Contrast
  L <- matestobj$Contrast
  
  # calculate x values
  if(class(matestobj)[2] == "ftest") {
    # this is a F-test object
    Lb <- L %*% t(b)
   # browser()
    if(nrow(L) == 1)
      xval <- Lb
    else
      xval <- apply(Lb, 2, function(x) sqrt(sum(x^2)))
  }
  else if(class(matestobj)[2] == "ttest") {
    # this is a T-test object
    # The volcano will always be 2-sided
    xval <- apply(L, 1, function(x) x%*%t(b))
  }        
  xval
}
