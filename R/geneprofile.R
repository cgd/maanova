######################################################################
#
# geneprofile.R
#
# copyright (c) 2001-2003, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written Feb, 2003
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
######################################################################

geneprofile <- function(anovaobj, term, geneidx,
                        col="blue", type="b", ylim, xlab, ylab, ...)
{

  if(class(anovaobj) != "maanova")
    stop("The first input variable is not an object of class maanova.")

  # gene and sample index
  if(missing(geneidx))
    geneidx <- 1:length(anovaobj$G)

  sampleidx <- 1:dim(anovaobj[[term]])[2]

  # data to be plotted
  data.plot <- anovaobj[[term]][geneidx,sampleidx]
  # create an empty plot
  if(missing(ylim))
    ylim <- range(data.plot)
  if(missing(xlab))
    xlab <- ""
  if(missing(ylab))
    ylab <- ""
  plot(1:length(sampleidx), ylim=ylim, type="n",
       xaxt="n", xlab=xlab, ylab=ylab, ...)
  
  # draw the lines
  apply(data.plot, 1, function(x) lines(x, col=col, type=type))
  
  # draw the axis
  axis(1, at=1:length(sampleidx), labels=anovaobj[[paste(term,".level",sep="")]])       
}
