######################################################################
#
# transform.madata.R
# This is the function to do data transformation
# This function was called smooth before.
#
# copyright (c) 2001, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written Nov, 2001
# Modified May 2004 for misc changes
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
#
######################################################################

require("stats")

# for rawdata
transform.rawdata <- function(`_data`, ...)
{
  transform.madata(`_data`, ...)
}

# for madata
transform.madata <-
  function(`_data`,
           method=c("shift","glowess","rlowess","linlog","linlogshift"),
           lolim, uplim,
           f=0.1, iter=3, degree=1,
           cg=0.3, cr=0.3, n.bin=10,
           draw=c("screen", "dev", "off"), ...)
  
{
  object <- `_data`
  if( !(class(object) %in% c("madata", "rawdata")) )
    stop("The input variable must be an object of madata or rawdata")

  # if there's smooth method applied to this object, warning
  if( !(object$TransformMethod %in% c("None","log2")) )
    warning("You have already smoothed the data... I will do it again anyway! ")

  # set the data to be transformed
  # if it's madata, the data field should be log2 based,
  # restore back to raw scale
  if( class(object)=="madata" ) 
    if(object$TransformMethod !="log2") 
      warning(" I assume data is log2 transformed. Otherwise, read madata using
      read.madata(... log.transform=T). \n")   
  object$data <- 2^(object$data)

  # transformation method
  method <- match.arg(method)

  # where to draw the plots
  draw <- match.arg(draw)
  
  # save old par parameters
  if (draw=="screen") {
    old.mar <- par("mar")
    old.las <- par("las")
    on.exit(par(las=old.las,mar=old.mar))
    par(las=1)
  }
  
  # if method is shift
  if(method %in% c("shift", "linlogshift")) { 
    if( missing(lolim) ) # calculate the default lolim
      lolim <- -min(object$data)
    if( missing(uplim) ) # calculate the default uplim
      uplim <- min(object$data)
  
    if(uplim <= lolim)
      stop("Offset upper limit must be greater than lower limit!")
  }

  if(method == "shift")
    object <- shift(object, lolim, uplim, draw)

  # if method is global lowess or loess
  if(method == "glowess")
    object <- glowess(object, method, f, iter, degree, draw)

  # if method is regional lowess or loess
  if(method == "rlowess") {
    # are the grid information there?
    if( sum(c("metarow","metacol","row" ,"col") %in% names(object)) != 4) {
      # object contains the grid location information
      msg <- "There's no grid location information in the input data object."
      msg <- paste(msg, method, "method failed!")
      stop(msg)
    }
    # calculate the grand row and column
    n.row <- max(object$row)
    n.col <- max(object$col)
    grow <- object$row + n.row*(object$metarow-1)
    gcol <- object$col + n.col*(object$metacol-1)
    object <- rlowess(object, method, grow, gcol, f, iter, degree, draw)
  }
  
  if(method == "linlog")
    object <- linlog(object, cg, cr, draw)

  if(method == "linlogshift")
    object <- linlogshift(object, lolim, uplim, cg, cr, n.bin, draw)

  # set the smooth method
  object$TransformMethod <- method

  # recalculate column means for madata
  if(class(object) == "madata")
    object$colmeans <- apply(object$data, 2, mean)

  object
}

  
shift <- function(object, lolim, uplim, draw)
{
  # get the data to be transformed
  data <- object$data
  
  # initialize variables
  sad <- rep(0,100)
  c <- c( seq(lolim, uplim, (uplim-lolim)/99) )
  offset <- rep(0,object$n.array*2)

  for (k in 1:object$n.array) { # loop thru the arrays
    cat("Smoothing array", k, "...\n")
    r <- data[,2*k-1]
    g <- data[,2*k]

    # find the offset
    for (i in 1:100) {
      r.shift <- r-c[i]
      r.shift[r.shift<1] <- 1
      g.shift <- g+c[i]
      g.shift[g.shift<1] <- 1
      x <- log2(r.shift) - log2(g.shift)
      sad[i] <- sum(abs(x-median(x)))
    }

    # find the offset value for this array
    offset[2*k-1] <- -c[sad==min(sad)]
    offset[2*k] <- -offset[2*k-1]
    
    # check if there's negative value after shift
    minr <- min(r)+offset[2*k-1]
    ming <- min(g)+offset[2*k]

    # if there are any negative value after shift,
    # adjust the offset value for this array
    if (minr <= 0) { # r is negative after shift
      offset[2*k-1] <- 1 - min(r);
      offset[2*k] <- -offset[2*k-1];
    }
    else if (ming <= 0) {# g is negative after shift
      offset[2*k] <- 1-min(g);
      offset[2*k-1] <- -offset[2*k];
    }
                
    if( draw != "off") { # draw the offset plot
      if( (draw=="screen") & (k!=1) ) {
        # open a window on screen
        dev.new()
#        if(.Platform$GUI == "AQUA")
#          quartz()
#        else
#          x11()
      }
      # setup the layout of the figures
      layout(matrix(1:4,2,2) )
      
      # plot RI plot before offset at figure 1
      par(mar=c(5.1,4.1,4.1,2.1))
      plot( log2(r)+log2(g), log2(r)-log2(g),
           xlab="log2(R*G)",ylab="log2(R/G)",
           col="blue", pch=4, cex=0.5, 
           main=paste("Array",k,"no offset") )
      
      # plot RI plot after offset at figure 2
      par(mar=c(5.1,4.1,4.1,2.1))
      plot( log2(r+offset[2*k-1])+log2(g+offset[2*k]),
           log2(r+offset[2*k-1])-log2(g+offset[2*k]),
           xlab="log2(R*G)",ylab="log2(R/G)",
           col="blue", pch=4, cex=0.5, 
           main=paste("Array",k," offset=",offset[2*k]) )

      # plot SAD curve at figure 3
      par(mar=c(5.1,4.1,4.1,2.1))
      plot(c, sad, xlab="Offset", ylab="SAD",
           col="blue", type="l", 
           main="SAD vs. Offset")
    }
  }
  # update the fields of the input variable
  if(is.null(object$n.rep))
    nreps <- 1
  else
    nreps <- object$n.rep

  if(is.null(object$n.gene))
    ngenes <- dim(object$data)[1]
  else
    ngenes <- object$n.gene
  offsetmat <- matrix(rep(offset,ngenes*nreps),
                      ngenes*nreps, object$n.array*2, byrow=TRUE)
  
  # set the output object
  if(class(object) == "madata")
    # for madata, the data should be log2 based
    object$data <- log2(data+offsetmat)
  if(class(object) == "rawdata") # for rawdata, data is in raw scale
    object$data <- object$data + offsetmat


  object
}


glowess <- function(object, method, f, iter, degree, draw)
{
  # init adjusted data
  adjdata <- array(0, dim(object$data))

  # loop thru arrays
  for (k in 1:object$n.array) { 
    cat("Smoothing array", k, "...\n")
    r <- object$data[,2*k-1]
    g <- object$data[,2*k]
    ratio <- log2(r/g)
    intensity <- log2(r*g)
    idx <- order(intensity)
    # call loess function in modreg to do loess
    l <- loess(ratio~intensity, span=f, degree=degree,
               control=loess.control(trace.hat="approximate"))
    adjdata[,2*k-1] <- log2(r) - l$fitted/2
    adjdata[,2*k] <- log2(g) + l$fitted/2

    # draw figure if requested
    if(draw != "off") {      
      if( (draw=="screen") & (k!=1) ) {
        # open a window on screen
        dev.new()
#        if(.Platform$GUI == "AQUA")
#          quartz()
#        else
#          x11()
      }
      # setup the layout of the figures
      layout(matrix(1:2,2,1) )
      # plot RI plot before offset at figure 1
      par(mar=c(5.1,4.1,4.1,2.1))
      plot(intensity, ratio,
           xlab="log2(R*G)",ylab="log2(R/G)",
           col="blue", pch=4, cex=0.5, 
           main=paste("Array",k,"before",method, sep=" ") )
      # plot the lowess line
      lines(intensity[idx], l$fitted[idx], col="red")
      # RI plot after lowess
      # recalculate ratio and intensity
      ratio <- adjdata[,2*k-1] - adjdata[,2*k]
      intensity <- adjdata[,2*k-1] + adjdata[,2*k]
      par(mar=c(5.1,4.1,4.1,2.1))
      plot(intensity, ratio,
           xlab="log2(R*G)",ylab="log2(R/G)",
           col="blue", pch=4, cex=0.5, 
           main=paste("Array",k,"after", method, sep=" ") )
    } # end if(draw)
  }
  
  # set the output object
  # not that now adjdata is log2 based
  if(class(object) == "madata") {
    object$data <- adjdata
  }
  if(class(object) == "rawdata") {
    # restore back to raw scale for rawdata
    object$data <- 2^adjdata
  }
  
  object
}


rlowess <- function(object, method, grow, gcol, f, iter, degree, draw)
{
  # init adjusted data
  adjdata <- array(0, dim(object$data))
  
  # control parameters for loess function
  loess.control(surface="interpolate", statistics="approximate",
               tract.hat="approximate", cell=0.2, iterations=4)
  
  for (k in 1:object$n.array) { # loop thru the arrays
    cat("Smoothing array", k, "...\n")
    r <- object$data[,2*k-1]
    g <- object$data[,2*k]
    ratio <- log2(r/g)
    intensity <- log2(r*g)
    idx <- order(intensity)
    # call loess(modreg) to do loess
    l <- loess(ratio~intensity+grow+gcol, span=f, degree=degree,
               control=loess.control(trace.hat="approximate"))
    adjdata[,2*k-1] <- log2(r) - l$fitted/2
    adjdata[,2*k] <- log2(g) + l$fitted/2

    # draw figure if requested
    if(draw!="off") {
      if( (draw=="screen") & (k!=1) ) {
        # open a window on screen
        dev.new()
#    if(.Platform$GUI == "AQUA")
#          quartz()
#        else
#          x11()
      }
      layout(matrix(1:2,2,1) )
      # plot RI plot before offset at figure 1
      par(mar=c(5.1,4.1,4.1,2.1))
      plot(intensity, ratio,
           xlab="log2(R*G)",ylab="log2(R/G)",
           col="blue", pch=4, cex=0.5,
           main=paste("Array",k,"before",method, sep=" ") )
      # plot the lowess line
      lines(intensity[idx], l$fitted[idx], col="red")
      # RI plot after lowess
      # recalculate ratio and intensity
      ratio <- adjdata[,2*k-1] - adjdata[,2*k]
      intensity <- adjdata[,2*k-1] + adjdata[,2*k]
      par(mar=c(5.1,4.1,4.1,2.1))
      plot(intensity, ratio,
           xlab="log2(R*G)",ylab="log2(R/G)",
           col="blue", pch=4, cex=0.5,
           main=paste("Array",k,"after", method, sep=" ") )
    } # end if(draw)
  } # end array loop
  
  # set the output object
  # now adjdata is log2 based
  if(class(object) == "madata") {
    object$data <- adjdata
  }
  if(class(object) == "rawdata") {
    object$data <- 2^adjdata
  }
  
  object

}



linlog <- function(object, cg, cr, draw)
{
  # init adjdata
  adjdata <- array(0, dim(object$data))
  
  # loop thru all arrays
  for(i in 1:object$n.array) {
    cat("Smoothing array", i, "...\n")
    adjdata[,i*2-1] <- linlog.engine(object$data[,i*2-1],cr)
    adjdata[,i*2] <- linlog.engine(object$data[,i*2],cg)
  }

  # do RI plot for all arrays before and after linear log
  if(draw != "off") {
    for(i in 1:object$n.array) {
      if( (draw=="screen") & (i!=1)) {
        # open a window on screen
        dev.new()
#        if(.Platform$GUI == "AQUA")
#          quartz()
#        else
#          x11()
      }
      layout(c(1,2))
      # before linlog
      r <- object$data[,i*2-1]
      g <- object$data[,i*2]
      plot( log2(r)+log2(g), log2(r)-log2(g),
           xlab="log2(R*G)",ylab="log2(R/G)",
           col="blue", pch=4, cex=0.5, 
           main=paste("Array",i,"before linear log") )
      # after linlog
      r <- adjdata[,i*2-1]
      g <- adjdata[,i*2]
      plot( r+g, r-g,
           xlab="log2(R*G)",ylab="log2(R/G)",
           col="blue", pch=4, cex=0.5, 
           main=paste("Array",i,"after linear log") )
    }
  }

  # set output object
  if(class(object) == "madata") {
    object$data <- adjdata
  }
  if(class(object) == "rawdata") {
    object$data <- 2^adjdata
  }
  
  object
}


linlogshift <- function(object, lolim, uplim, cg, cr, n.bin, draw)
{
  # init variables
  adjdata <- array(0, dim(object$data))
  offset <- rep(0,object$n.array*2)
  sad <- rep(0,100)
  c <- c( seq(lolim, uplim, (uplim-lolim)/99) )

  for (k in 1:object$n.array) { # loop thru the arrays
    cat("Smoothing array", k, "...\n")
    r <- object$data[,2*k-1]
    g <- object$data[,2*k]
    
    # find the offset
    for (i in 1:100) {
      r.shift <- r-c[i]
      r.shift[r.shift<1] <- 1
      g.shift <- g+c[i]
      g.shift[g.shift<1] <- 1
      r.linlog <- linlog.engine(r.shift, cr)
      g.linlog <- linlog.engine(g.shift, cg)
      x <- r.linlog - g.linlog
      sad[i] <- sum(abs(x-median(x)))
    }
    # find the offset value for this array
    offset[2*k-1] <- -c[sad==min(sad)]
    offset[2*k] <- -offset[2*k-1]

    # check if there's negative value after shift
    minr <- min(r)+offset[2*k-1]
    ming <- min(g)+offset[2*k]

    # if there are any negative value after shift,
    # adjust the offset value for this array
    if (minr <= 0) { # r is negative after shift
      offset[2*k-1] <- 1 - min(r);
      offset[2*k] <- -offset[2*k-1];
    }
    else if (ming <= 0) {# g is negative after shift
      offset[2*k] <- 1-min(g);
      offset[2*k-1] <- -offset[2*k];
    }

    # do linear-log on shifted data
    r.shift <- r + offset[2*k-1]
    g.shift <- g + offset[2*k]
    # update adjdata
    adjdata[,2*k-1] <- linlog.engine(r.shift, cr)
    adjdata[,2*k] <- linlog.engine(g.shift, cg)

    # draw the figure
    if(draw != "off") {
      if( (draw=="screen") & (k!=1)) {
        # open a window on screen
        dev.new()
#        if(.Platform$GUI == "AQUA")
#          quartz()
#        else
#          x11()
      }
      # setup the layout of the figures
      layout(matrix(1:4,2,2) )

      # RI plot before transformation
      plot( log2(r)+log2(g), log2(r)-log2(g),
           xlab="log2(R*G)",ylab="log2(R/G)",
           col="blue", pch=4, cex=0.5,
           main=paste("Array",k,"before linlog shift") )

      # RI plot after transformation
      r <- adjdata[,k*2-1]
      g <- adjdata[,k*2]
      plot( r+g, r-g,
           xlab="log2(R*G)",ylab="log2(R/G)",
           col="blue", pch=4, cex=0.5, 
           main=paste("Array",k,"after linlog shift") )

      # SAD curve
      plot(c, sad, xlab="Offset", ylab="SAD",
           col="blue", type="l",
           main="SAD vs. Offset")

      # variance plot
      logsum <- (adjdata[,2*k-1] + adjdata[,2*k])/2
      logdiff <- adjdata[,2*k-1] - adjdata[,2*k]
      ratioVarplot(logsum, logdiff, n.bin)
      title("Variance plot")
    } # end of drawing figure
  } # end array loop

  # set output object
  if(class(object) == "madata") {
    object$data <- adjdata
  }
  if(class(object) == "rawdata") {
    object$data <- 2^adjdata
  }

  object
}

ratioVarplot <- function(logsum, logdiff, n)
{
  xlim <- range(logsum)
  for(i in 1:n) {
    x <- which( (logsum>(min(logsum)+i*(max(logsum)-min(logsum))/n)) &
                (logsum<(min(logsum)+(i+1)*(max(logsum)-min(logsum))/n)) )
    if(length(x)==0 | length(x)==1) {
      y <- 0
      My <- 18
    }
    else {
      y <- log(var(logdiff[x]))
      My <- mean(logsum[x])
    }
    if(i==1)
      plot(My, y, xlab="(linlog(R)+linlog(G))/2",
           ylab="variance of linlog(R)-linlog(G)",
           col="blue", pch=4, cex=0.5,
           xlim=xlim,ylim=c(-10,5))
    else
      points(My, y, col="blue", pch=4, cex=0.5)
  }
}
      
                 


linlog.engine <- function(data, cutoff)
{
  n <- length(data)
  transdata <- rep(0,n)
  # index
  idx <- ceiling(n*cutoff)
  # cutoff value
  cutoff.value <- sort(data)[idx]

  # linear-log transform data
  transdata[data>cutoff.value] <- log2(data[data>cutoff.value])
  transdata[data<=cutoff.value] <- data[data<=cutoff.value]/
    cutoff.value + log2(cutoff.value) - 1

  transdata
}
