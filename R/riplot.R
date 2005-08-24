
######################################################################
#
# riplot.R
#
# copyright (c) 2001, Hao Wu and Gary A. Churchill, The Jackson Lab.
#
# written Nov, 2001
# Modified Nov, 2002 
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# 
######################################################################

riplot <-
  function(object, title, array, color="blue", highlight.flag=TRUE,
           flag.color="Red", idx.highlight, highlight.color="Green",
           rep.connect=FALSE, onScreen=TRUE)
{
  # it works for 2-dye system only
  if(object$n.dye != 2)
    stop("riplot works for 2-dye array only")

  if( !(class(object) %in% c("madata", "rawdata")) )
    stop("The input variable must be an object of madata or rawdata")
  
  # local variables
  if( class(object) == "madata" ) # if this is madata
    x <- object$data
  else if( class(object) == "rawdata" ) # if this is rawdata
    x <- log2(object$data)
  else
    stop("The input variable must be an object of rawdata or madata")

  x.dim <- dim(x)

  # flag
  flag <- object$flag
  
  # calculate R (ratio) and I (intensity)
  R <- x[,seq(1,x.dim[2],2)] - x[,seq(2,x.dim[2],2)]
  I <- x[,seq(1,x.dim[2],2)] + x[,seq(2,x.dim[2],2)]

  if(!is.matrix(R)) # if R is a vector, not a matrix
    R <- matrix(R, length(R),1) # convert R to a matrix
  # do same thing to I
  if(!is.matrix(I))
    I <- matrix(I, length(I),1)
  
  tmp <- max(abs(R));

  if(missing(title)){
    title <- NULL
    for(i in 1:object$n.array) 
      title[i] <- paste("RI plot for array number",i)
  }

  # now draw the figures
  if(missing(array))
    array <- 1:object$n.array
  for (i in array) {
    if(onScreen) {
      get(getOption("device"))()
      # open a window on screen
#      if(.Platform$GUI == "AQUA")
#        quartz()
#      else
#        x11()
    }
    plot( I[,i], R[,i], xlim=c(min(I),max(I)), ylim=c(-tmp,tmp),
         xlab=expression(log[2](R*G)), ylab=expression(log[2](R/G)),
         col=color, pch=4, cex=0.5,
         main=title[i] )
    # if highlight is given, redraw those points
    if(!missing(idx.highlight)) {
      if(class(object) == "madata") {
        idx.gene <- matrix(1:(object$n.gene*object$n.rep),
                     object$n.rep, object$n.gene)
        high <- as.vector(idx.gene[,idx.highlight])
      }
      else
        high <- idx.highlight
      points(I[high,i], R[high,i], col=highlight.color, pch=4, cex=0.5)
    }
    
    # highlight the flagged spot (if any)
    if(highlight.flag & !is.null(flag)) {
      high <- flag[,i]!=0
      points(I[high,i], R[high,i], col=flag.color, pch=4, cex=0.5)
    }
        
    if( rep.connect & (class(object)=="madata") ) {
      if(object$n.rep!=1) {
        # connect the dots between replicates (if any)
        idx.rep <- as.vector(repmat(t(1:object$n.gene), object$n.rep,1))
        for(j in 1:object$n.gene) {
          x <- I[idx.rep==j,i]
          y <- R[idx.rep==j,i]
          lines(x,y, type="l", col="grey")
        }
      }
    }
  }
  # finish loop for arrays
}



