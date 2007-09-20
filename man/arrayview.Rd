\author{Hao Wu}
\name{arrayview}
\alias{arrayview}
\title{View the layout of input data}
\description{
  This function reconstructs the input data according to the Microarray
  grid location structure and plots the data according to the user
  specified color map.

  By default, it will plot the log ratios for 2-dye array and raw
  intensity for 1-dye array. It does not work for N-dye (N>2) array at
  this time.
  
  Note that if user collapsed the replicates by using 'avgrep' in
  \code{\link[maanova]{read.madata}}, and then 
  arrayview will not be available. 
}

\usage{
arrayview(object, ratio, array, colormap, onScreen=TRUE, ...)
}

\arguments{
  \item{object}{An object of class \code{madata}.}
  \item{ratio}{The data to be plotted. The length of it must be
    equal to the length of the grid locations, .e.g, madata$row and
    madata$col. If ratio is a vector, there will be one plot. If ratio
    is a matrix, there will be one plot for each column. If ratio is not
    provided, \code{link[maanova]{make.ratio}} will be called to
    calculate the ratios from the original data.}
  \item{array}{A list of arrays to be plotted. This variable is only
    valid when \code{ratio} is not provided. Whenever \code{ratio} is
    provided, all columns in ratio will be plotted.}
  \item{colormap}{User specified color map. See
    \code{\link[base]{colors}} for more detail.}
  \item{\dots}{Other parameters to be passed to
    \code{\link[base]{image}}.}
  \item{onScreen}{A logical value to represent whether to display the
    plots on screen or not. If TRUE, x11() (in Unix/Windows) or
    macintosh (in Mac) will be called inside the function. Otherwise,
    it will plot the figure on the current device. Default is TRUE.}
}

\examples{
\dontrun{
data(kidney)
# arrayview data on screen
arrayview(kidney.raw, array=1)
graphics.off()
# arrayview raw data array 1 and 3 and output to postscript file
postscript(file="kidneyArrayview.ps")
arrayview(kidney.raw, array=c(1,3), onScreen=FALSE)
}
}

\keyword{hplot}
