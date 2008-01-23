######################################################################
#
# fitmaanova.R
#
# copyright (c) 2001-2004, Hao Wu and Gary A. Churchill, The Jackson Lab.
#
# written Dec, 2003
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
#
######################################################################


fom <- function(anovaobj, idx.gene, term, ngroups)
{
  # check input data
  if (class(anovaobj) != "maanova")
    stop("The first input variable is not an object of class maanova.")
  if( !(term %in% names(anovaobj)) )
    stop(paste(term, "is not a field in anovaobj"))
         
  # create local variables
  VG <- anovaobj[[term]][idx.gene,]
  nvars <- dim(VG)[2]

  # max number of groups
  maxgrp <- min(30, length(idx.gene)-1)
  if(missing(ngroups)) # calculate FOM for only 1 group to 30 groups
    ngroups <- 2:maxgrp
  else {
    if(ngroups < 2)
      stop("Number of groups must be greater than or equal to 2")
    if(length(ngroups) == 1) {
      # provide ngroups as integer
      if(ngroups>maxgrp)
        stop("Group number cannot be greater than number of genes!")
      else
        ngroups <- 2:ngroups
    }
  }

  # sample index
  sampleidx <- 1:nvars

  # normalize VG, subtract the column mean and divided by column sd
  m <- apply(VG, 2, mean)
  sd <- apply(VG, 2, sd)
  normVG <- t(apply(VG, 1, function(x) (x-m)/sd))

  # output variable
  f <- rep(0, length(ngroups))

  # calculate FOM
  i <- 1
  for (ngrp in ngroups) { #loop thru different number of clusters
    for (samplej in sampleidx) { # loop thru all samples
      # leave out sampleidx(samplej) and do cluster
      sampleset <- setdiff(sampleidx, samplej)
      # K-mean cluster the VG without for sampleidx(samplej)
      class <- kmeans(normVG[,sampleset], ngrp)
      # note that class is a vector with length(geneidx) to indicate which
      # group each genes belong to
      sigma <- 0 # weighted variance for the leave-out condition
      for(k in 1:ngrp) { # loop thru the groups
        idx <- which(class$cluster == k) # find the genes in group k
        # calculate the within group variance for VGstar
        if(length(idx) > 1) {# if there are any genes in this group
          # proportion of genes in group k
          alpha <- length(idx)/length(idx.gene)
          sigma <- sigma + alpha*var(normVG[idx,samplej])
        }
      }
      # calculate adjusted FOM for ngrp groups
      f[i] <- f[i] + sqrt(sigma)
    }
    # adjust FOM
    f[i] <- f[i] / sqrt((length(idx.gene)-ngrp)/length(idx.gene))
    i <- i+1
  }

  # make the FOM plot
  plot(ngroups, f, col="blue", main="FOM plot", type="l",
       xlab="Number of groups", ylab="Figure of Merit value")
  points(ngroups, f, col="blue")

  invisible(f)
  
}


