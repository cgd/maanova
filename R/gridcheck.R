######################################################################
#
# gridcheck.R
#
# copyright (c) 2002, Hao Wu and Gary A. Churchill, The Jackson Lab.
#
# written Nov, 2002
#
# Modified Dec, 2002 for mixed effect model
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
#
######################################################################

gridcheck <- function(rawdata, array1, array2, highlight.flag=TRUE,
                      flag.color="Red", margin=c(3.1,3.1,3.1,1.1))
{
  if(class(rawdata) != "rawdata")
    stop("The first input variable is not an object of class rawdata.")

  if(rawdata$n.dye != 2)
    stop("gridcheck only works for 2-dye arrays")

  # if there's no grid info, stop
  if( sum(c("metarow","metacol","row" ,"col") %in% names(rawdata)) != 4) {
    # rawdata contains the grid location information
    msg <- "There's no grid location information in the input data object."
    msg <- paste(msg, "You cannot do grid checking!")
    stop(msg)
  }

  # get metarow and metacol and flag
  mrow <- rawdata$metarow
  mcol <- rawdata$metacol
  n.mrow <- max(mrow)
  n.mcol <- max(mcol)
  flag <- rawdata$flag

  # save old par parameters and setup new layout
  old.mar <- par("mar")
  old.las <- par("las")
  on.exit(par(las=old.las,mar=old.mar))
  par(las=1)

  if(missing(array2)) {
    # no array 2, plot self comparison for array 1
    if(missing(array1))
      # no array 1, plot for all arrays
      array1 <- 1:rawdata$n.array
    # loop for all arrays
    for (i.array in array1) {
      if(i.array != array1[1]) {
        # open a new window
        get(getOption("device"))()
#        if(.Platform$GUI == "AQUA")
#          quartz()
#        else
#          x11()
      }
      # setup the layout and margin
      layout(matrix(1:(n.mrow*n.mcol), n.mrow, n.mcol, byrow=TRUE))
      par(mar=margin)
      for(i in 1:n.mrow) {
        for(j in 1:n.mcol) {
          idx <- which((mrow==i) & (mcol==j))
          plot(log2(rawdata$data[idx,i.array*2-1]),
               log2(rawdata$data[idx,i.array*2]),
               col="blue", pch=4, cex=0.5, xlab="", ylab="")
          # highlight the flagged spot (if any)
          if(highlight.flag & !is.null(flag)) {
            high <- flag[idx,i.array]!=0
            points(log2(rawdata$data[idx[high],i.array*2-1]),
                   log2(rawdata$data[idx[high],i.array*2]),
                   col=flag.color, pch=4, cex=0.5)
          }
        }
      }
    }
  }
  else {
    # have array 2, compare the same sample for array 1 and array 2
    if(missing(array1))
      stop("Miss the first array number")
    if((length(array1)!=1) | (length(array2)!=1) )
      stop("Both array1 and array2 must be an integer")
    
    # get the sample ids for array 1 and array 2 from design
    if(is.null(rawdata$design))
      stop("No experimental design information in rawdata. Cannot do grid check on two arrays.")
    sample1 <- rawdata$design$Sample[c(array1*2-1, array1*2)]
    sample2 <- rawdata$design$Sample[c(array2*2-1, array2*2)]
    if(length(intersect(sample1, sample2)) == 0)
      stop(paste("No common sample in array", array1, "and array", array2,
                 "Cannot do grid check"))
    # start plot
    nplot <- 0
    # get the data for two arrays
    data1 <- rawdata$data[,c(array1*2-1, array1*2)]
    data2 <- rawdata$data[,c(array2*2-1, array2*2)]
    for(i.array1 in 1:2) {
      for(i.array2 in 1:2) {
        if(sample1[i.array1] == sample2[i.array2]) {
          nplot <- nplot + 1
          if(nplot!=1) {
            # open a window on screen
            get(getOption("device"))()
#            if(.Platform$GUI == "AQUA")
#              quartz()
#            else
#              x11()
          }

          # setup the layout and margin
          layout(matrix(1:(n.mrow*n.mcol), n.mrow, n.mcol, byrow=TRUE))
          par(mar=margin)
          for(i in 1:n.mrow) {
            for(j in 1:n.mcol) {
              idx <- which((mrow==i) & (mcol==j))
              plot(log2(data1[idx, i.array1]), log2(data2[idx, i.array2]),
                   col="blue", pch=4, cex=0.5, xlab="", ylab="")
              # highlight the flagged spot (if any)
              if(highlight.flag & !is.null(flag)) {
                high <- sumrow(flag[idx, c(array1, array2)])!=0
                points(log2(data1[idx[high], i.array1]), log2(data2[idx[high], i.array2]),
                       col=flag.color, pch=4, cex=0.5)
              }
            }
          }
        }
      }
    }
  }

}
      
  
  
