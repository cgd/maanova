\author{Hao Wu}

\name{transform.madata}
\alias{transform.madata}
\alias{transform.rawdata}

\title{Micro Array experiment data transformation}

\description{
  This is the function to transform the Micro Array experiment data
  based on the given method. 
}

\details{
  The smoothing methods include:
  
  shift -- the calculation of offset is based on the minimum
  sum of absolute deviations (SAD). Each array will have its own offset
  value. The data after shift cannot be smaller than 1.

  glowess -- Intensity-based lowess. This method is to smooth the scatter
  plot of Ratio (R/G) versus Intensity (R*G). The formula in the fitting
  is ratio ~ intensity. 

  rlowess -- Joint lowess. This method is to smooth the scatter plot
  of Ratio versus Intensity and grid locations. It is the joint of
  intensity-based lowess and spatial loess. You have to have the grid
  location for every spot in order to use this method. The formula in
  fitting is ratio ~ intesity + row + col. 

  linlog -- Linear-log transformation.

  linlogshift -- Linear-log shift transformation.
  
  Previously, intensity lowess was called global lowess and joint lowess
  was called regional lowess. So I use "glowess" and "rlowess" in the
  method. Although the method names doesn't make too much sense, I will
  keep them for the reason of backward compatibility.
  
  If you have replicated spots and want to collapse them in
  \code{\link[maanova]{createData}} by providing avgreps=1 or 2, you will
  lose grid information and joint lowess will be unavailable.

  This function works on both \code{rawdata} and \code{madata} objects.
  For \code{rawdata}, the result transformed data will be on the raw
  scale. For \code{madata}, the result transformed data will be on log2
  based scale.

  Note that this function is only working for two-dye array at this
  time. 
}

\usage{
transform.madata(`_data`,method=c("shift","glowess","rlowess","linlog","linlogshift"),
       lolim, uplim, f=0.1, iter=3, degree=1, cg=0.3, cr=0.3, n.bin=10,
       draw=c("screen", "dev", "off"), \dots)
}

\arguments{
  \item{_data}{An object of class \code{madata} or \code{rawdata}.}
  \item{method}{The smoothing method.}
  \item{lolim}{Low shift limit. If this argument is missing, the negative
    of the minimum element in the pmt data is used.}
  \item{uplim}{High shift limit. If this argument is missing, the minimum
    element in the pmt data is used. lolim and uplim are applicable only if
    the method is "shift" or "linlogshift".}
  \item{f}{The smoother span. This gives the proportion of points in the
    plot which influence the smooth at each value. Larger values give
    more smoothness. It is equivalent to the "span" parameter in
    \code{\link[modreg]{loess}}.}
  \item{iter}{The number of robustifying iterations which should be
    performed. Using smaller values of iter will make lowess run
    faster.}
  \item{degree}{The degree of the polynomials to be used in loess, up to
    2. This is used when method is "glowess" or "rlowess".}
  \item{cg}{Percentage of genes to be transformed linearly for the green
    channel.}
  \item{cr}{Percentage of genes to be transformed linearly for the red
    channel.}
  \item{n.bin}{Number of bins for calculating the variance after
    linlogShift.}
  \item{draw}{Where to plot the transformation plots. "off" means
    no plot. "screen" means to display the plots on screen then
    x11() (in Unix/Windows) or
    macintosh() (in Mac) will be called inside the function.
    "dev" means to output the plots to the current device. User can use
    this option to output the plot to a file. Default option is
    "screen". }
  \item{\dots}{Ignored at this point.}
}

\section{Value}{
  The return value is an object of class \code{madata} or
  \code{rawdata}. Compared with the input object, the following fields
  are changed:
  \itemize{
    \item Field \code{data} is the transformed data.
    \item Field \code{TransformMethod} will be the transformation method
    applied.
  }
}

\examples{
# load in data
data(kidney)
# do regional loess on raw data
\dontrun{
raw.lowess <- transform.madata(kidney.raw, method="rlowess")
graphics.off()

# make data object and do normalization on it
# create data object with collapsing replicates
data1<-createData(kidney.raw)
# do shift without displaying the plot
data1.shift <- transform.madata(data1, method="shift", lolim=-50, uplim=50,
         draw="off")

# do global lowess and output the plots to a postscript file
postscript(file="glowess.ps")
data1.glowess <- transform.madata(data1, method="glowess", draw="dev")
graphics.off()

# do linear-log
data1.linlog <- transform.madata(data1, method="linlog")
graphics.off()

# do linear-log shift
data1.linlogshift <- transform.madata(data1, method="linlogshift", lolim=-50,
  uplim=50)
graphics.off()
}
}

\references{
  Kerr and Churchill(2001), Statistical design and the analysis of gene
  expression microarrays, \emph{Genetical Research}, \bold{77}:123-128.
  
  Kerr, Martin and Churchill(2000), Analysis of variance for gene expression
  microarray data, \emph{Journal of Computational Biology},
  \bold{7}:819-837.

  Cui, Kerr and Churchill(2002), Data transformations for cDNA
  Microarray data, submitted, find manuscript in www.jax.org/research/churchill.
}


\seealso{
  \code{\link[modreg]{loess}}
}


\keyword{smooth}

