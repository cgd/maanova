######################################################################
#
# varplot.R
#
# copyright (c) 2001-2004, Hao Wu and Gary A. Churchill, The Jackson Lab.
# Written Sep 2003
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
######################################################################

varplot <- function(anovaobj, xlab, ylab, main){

 
  if (class(anovaobj) != "maanova")
    stop("The first input variable is not an object of class maanova.")

  s2 <- anovaobj$S2

  line.color <- c("black", "blue", "red", "green", "yellow", "cyan")
  
  if(is.null(s2))
    stop("No variance component in input object")

  # calculate number of breaks
  npts <- dim(s2)[1]
  nlevel <- dim(s2)[2]
  nbreaks <- npts/10
  x <- NULL
  y <- NULL

  for(i in 1:nlevel) {
    varcom <- s2[,i]
    tmp <- density(sqrt(s2[,i]))
    x <- cbind(x, tmp$x)
    y <- cbind(y, tmp$y)
  }
  if(missing(xlab))
    xlab <- "Sigma"
  if(missing(ylab))
    ylab <- "Density"
  if(missing(main))
    main = "Density plot for sqrt of variance"
  plot(x, y, type="n", xlab=xlab, ylab=ylab, main=main)
  
  for(i in 1:nlevel)
    lines(x[,i], y[,i], col=line.color[i])

  # make legend
  s2level <- c(anovaobj$S2.level, "error")
  xlim <- par("xaxp")
  ylim <- par("yaxp")
  xpos <- 0.3*xlim[1] + 0.7*xlim[2]
  ypos <- 0.1*ylim[1] + 0.9*ylim[2]
  legend(xpos, ypos, s2level, col=line.color[1:nlevel], lty=1)
}
  
