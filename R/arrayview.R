######################################################################
#
# arrayview.R
#
# copyright (c) 2001, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written Nov, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# View the spatial patter of the logratios for a 2-dye array
#
######################################################################

arrayview <-
  function(object, ratio, array, colormap, onScreen=TRUE, ...)
{
  # local variables
  ndye <- object$n.dye
  
  # stop if there's no grid information
  if(sum(c("metarow","metacol","row" ,"col") %in% names(object)) != 4)
    stop("There's no grid location information in input data.")
  
  # stop if the input object is not madata or rawdata
  if(!(class(object) %in% c("madata","rawdata")) )
    stop("The input variable must be an object of rawdata or madata")

  # if object is madata and reps were collapsed, cannot do arrayview
  if(class(object)=="madata")
    if (object$collapse==TRUE) 
      stop("The replicates were collapsed. Arrayview is unavailable.")
  
  # list of arrays
  if(missing(array))
    array <- 1:object$n.array
  
  if(missing(ratio)) { # if ratio is not provided
    if(ndye == 1) {
      ratio <- object$data
    }
    else if(ndye == 2) {
      cat("No input ratio, make.ratio is being called.\n")
      ratio <- make.ratio(object)
      total.spot <- dim(ratio)[[1]]
      n.array <- dim(ratio)[[2]]
    }
    else {
      stop("arrayview only works for 1 or 2 dye array at this time")
    }
  }
  else{ # ratio is provided, check if it's valid
    if(is.vector(ratio)) { # ratio is a vector
      total.spot <- length(ratio)
      n.array <- 1
    }
    else if(is.matrix(ratio)) { # ratio is a matrix
      total.spot <- dim(ratio)[[1]]
      n.array <- dim(ratio)[[2]]
    }
    # list of arrays
    array <- 1:n.array
    if(total.spot != length(object$row))
      stop("Number of elements in ratio doesn't match row record.")
  }

  # construct data.grid
  n.row <- max(object$row)
  n.col <- max(object$col)
  grow <- object$row + n.row*(object$metarow-1)
  gcol <- object$col + n.col*(object$metacol-1)
  n.grow <- max(grow)
  n.gcol <- max(gcol)
  data.grid <- matrix(rep(0,n.grow*n.gcol), n.grow, n.gcol)

  # assign default color map
  if(missing(colormap)) {
    r <- c(31:0,rep(0,32))/31
    g <- c(rep(0,32),0:31)/31
    b <- rep(0,64)
    colormap <- NULL
    for(i in 1:64)
      colormap[i] <- rgb(r=r[i],g=g[i],b=b[i])
  }

  for(i in array) { # loop for all arrays
    for(j in 1:total.spot) { #loop over all spots
      if(is.matrix(ratio)) # if ratio is a matrix
        data.grid[grow[j], gcol[j]] <- ratio[j,i]
      else # ratio is a vector
        data.grid[grow[j], gcol[j]] <- ratio[j]
    }
    if(onScreen) {
      # open a window on screen
      get(getOption("device"))()
#      if(.Platform$GUI == "AQUA")
#        quartz()
#      else
#        x11()
    }
    image(1:n.gcol, 1:n.grow, t(data.grid), ylab="Row",
          xlab="Column", col=colormap, main=paste("Array", i),... )
  }

}
       
  
