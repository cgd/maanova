######################################################################
#
# summarytable.R
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
######################################################################

summarytable <-
  function(matestobj, method=c("Fold.change","Pvalperm","adjPvalperm"), 
         test =c("F1","Fs"),whichTest=c("F1.Pvalperm","F1.adjPvalperm", 
           "Fs.Pvalperm","Fs.adjPvalperm"), 
            threshold, outfile="summarytable.csv")
{
 
  if(class(matestobj)[1] != "matest")
    stop("The First input variable is not an object of class matest.")
  fold = FALSE
  if(missing(test)){
    alltest =c("F1", "Fs")
    test = alltest[alltest %in% names(matestobj)]
  }
  else{
    alltest = c("F1", "Fs")
    stest = test[test %in% names(matestobj)]
    notest =setdiff(test, stest) 
    if(length(notest) > 0) 
      warning(paste(notest, ' is not available. ', sep="", collapse=""))
    test = stest
  }

  if(length(test) <1 ) 
    stop("There is no available test statistics.")
  
  if(missing(method)){
    allmeth = c("Pvalperm", "adjPvalperm")
    smethod = allmeth[allmeth %in% names(matestobj[[test[1]]])]
    fold = TRUE
  }
  else{
    smethod = method[method %in% names(matestobj[[test[1]]])]
    nomethod= setdiff(setdiff(method, smethod), 'Fold.change') 

    if( "Fold.change" %in% method) fold = TRUE
    if(length(nomethod) > 0) 
      warning(paste(nomethod, ' is not available. ', sep="", collapse=""))
  }
  if(missing(whichTest)){
    whichTest = NULL
  }

  if(class(matestobj)[2] == "ttest")
    # this is a T-test result
    result <- summarytable.ttest(matestobj, smethod, test, whichTest, threshold, fold)
  else
    # this is a F test 
    result <- summarytable.ftest(matestobj, smethod, test, whichTest, threshold, fold)
  write.csv(result, file=outfile)
  result = result
}


######################################################################
# function to summarize f-test
######################################################################
summarytable.ftest <- function(matestobj, method, test, whichTest,
  threshold, fold){
  xvalue = NULL
  if(length(method)>0){
    pval <- getPval.table(matestobj, method, test,1)
    pvalname= pval$pvalname; pval = pval$pval; 
  }
  else{pvalname= NULL; pval = NULL}
  if(fold == TRUE){
    Fold.change <- matrix(calVolcanoXval(matestobj), ncol=1)
    pval = cbind(Fold.change, pval)
    pvalname=c('Fold.change',pvalname)
  }
  nc=ncol(pval)
  id = matestobj$probeid
  if(length(whichTest)>0){
    if(missing(threshold)){
      warning('No threshold infomation. Save all'); threshold = 1
    }
    whicht= which(pvalname %in% whichTest); lw = length(whicht)
    if(lw < 1) warning('whichTest is not a vaild. Save all results.\n')
    else{
      if(lw==1) this = pval[,whicht] < threshold
      else{
        warning('More than one whichTest. Maximum set is used.\n')
        this = apply(pval[,whicht] < threshold, 1, max) 
      }
      pval <- pval[this==1, ]
      if(nc == 1) pval = matrix(pval, ncol=1)
      id = id[this]
    }
  }
  else warning('No whichTest infomation. Save all');
  
  # check to make sure probe count is divisible by num pvalues
  if((length(id) %% nrow(pval)) != 0)
  {
	  stop("bad matest object. probe id count is not divisible by f-test pvalue count");
  }
  
  # if there are replicates we end up with repeat id's. by taking the
  # modulus we remove the repeats
  if(length(id) > nrow(pval))
  {
	  num_repeats <- length(id) %/% nrow(pval);
      id <- id[seq(from=1, to=length(id), by=num_repeats)];
  }
  
  rownames(pval)=id;
  colnames(pval)=pvalname; 
  pval
}



######################################################################
# function for T-test 
######################################################################
summarytable.ttest <- function(matestobj,method, test, whichTest, 
  threshold, fold) 
{
  anovaobj <- matestobj$obsAnova

  diff.terms <- matestobj$term
  Contrast <- matestobj$Contrast
  # number of plots to be generated
  nplots <- dim(Contrast)[1]
  id = matestobj$probeid

  if(fold == TRUE){
    xvalue.all <- calVolcanoXval(matestobj)
    if(nplots==1) xvalue.all = matrix(xvalue.all, ncol=1)
  }
  # loop for contrasts

  levels <- NULL
  for(i in 1:length(diff.terms)) {
    # level names
    levels <- c(levels, anovaobj[[paste(diff.terms[i], ".level", sep="")]])
  }
  name = NULL; pvalall = NULL; tname=NULL; lm = length(method)
  for(icon in 1:nplots) {
    label = levels[Contrast[icon, ] !=0]
    label <- paste(paste(label, sep="*"), collapse="-")

    if(lm>0){
      pval <- getPval.table(matestobj, method, test, icon)
      pvalname= pval$pvalname; pval = pval$pval;
    } 
    if(fold == TRUE){
      pval <-cbind(xvalue.all[,icon], pval)
      pvalname=c("Fold.change", pvalname)
    }
    pvalall = cbind(pvalall, pval)
    name=c(name, paste(paste(label, sep="*", collapse=""),pvalname))
    tname=c(tname,pvalname)
  }
  nc = ncol(pvalall)
  if(length(whichTest)>0){
    if(missing(threshold)){
     warning('No threshold infomation. Save all'); threshold = 1
    }
    whicht= which(tname %in% whichTest); lw = length(whicht)
    if(lw<1) warning('whichTest is not a vaild. Save all results.\n')
    else{
      if(lw>1) warning('More than one whichTest. Maximium set is saved.')
      if(lw==1) this = pvalall[,whicht] < threshold
      else this = apply(pvalall[,whicht] < threshold, 1, max) 
      pvalall <- pvalall[this==1, ]
      if(nc==1) pvalall = matrix(pvalall, ncol=1)
      id = id[this]
    }
  }
  else warning('No whichTest infomation. Save all');
  
  # check to make sure probe count is divisible by num pvalues
  if((length(id) %% nrow(pvalall)) != 0)
  {
	  stop("bad matest object. probe id count is not divisible by t-test pvalue count");
  }
  
  # if there are replicates we end up with repeat id's. by taking the
  # modulus we remove the repeats
  if(length(id) > nrow(pvalall))
  {
	  num_repeats <- length(id) %/% nrow(pvalall);
	  id <- id[seq(from=1, to=length(id), by=num_repeats)];
  }
  
  rownames(pvalall) = id;
  colnames(pvalall) = name;
  pvalall
}



###########################################
# function to get P values from matest object
###########################################
getPval.table <- function(matestobj, method, test, ii){
  pval = NULL; pvalname= NULL
  for(i in 1:length(test)){
    for(j in 1:length(method)){
      Fobj = matestobj[[test[i]]]
      Fobj = Fobj[[method[j]]]  
      pval = cbind(pval, Fobj[,ii])
      thisname = paste(test[i], '.', method[j], sep="", collapse="")
      pvalname = c(pvalname, thisname)
    }
  }
  pval = list(pval=pval, pvalname=pvalname)
  pval 
}
