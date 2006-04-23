######################################################################
#
# resiplot.R
#
# copyright (c) 2001, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written Nov, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# 
######################################################################

resiplot <- function(madata, anovaobj, header)
  
{
  # check input variables
  if(class(madata) != "madata")
    stop("This first input variable is not an object of class madata.")
  if(class(anovaobj) != "maanova")
    stop("This second input variable is not an object of class maanova.")

  # local variables
  y <- madata$data
  yhat <- anovaobj$yhat
  ndyes <- madata$n.dye
  color <- c("red","green", "blue", "black")
  
  # calculate the residual
  resi <- y - yhat

  # create a dye row vector
  dye <- rep(1:ndyes, madata$n.array)

  # figure title
  if(missing(header))
    header <- "Residual vs. Yhat plot"
  
  # save old par parameters
  old.mar <- par("mar")
  old.las <- par("las")
  old.mfrow <- par("mfrow")
  old.mfcol <- par("mfcol")
  on.exit(par(las=old.las,mar=old.mar,mfrow=old.mfrow,mfcol=old.mfcol))
  
  par(las=1)
  layout(matrix(1:ndyes,ndyes))

  # xlim and ylim
  xlim <- range(yhat)
  ylim <- range(resi)
  
  # plot the residual
  for(i in 1:ndyes) 
    plot(x=yhat[,dye==i], y=resi[,dye==i], xlim=xlim, ylim=ylim,
         xlab="Yhat", ylab="Residual", col=color[i],
         pch=4, cex=0.5, main=paste(header,"Dye", i))
}
