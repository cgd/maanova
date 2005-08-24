\author{Hao Wu \email{hao@jax.org} }

\name{subset.madata}
\alias{subset.madata}

\title{Subsetting Micro Array data objects}

\usage{
subset.madata(x, arrays, genes, \dots)
}

\description{
  Return subsets of an an object of class \code{madata} meeting given
  conditions.
}

\arguments{
  \item{x}{An object of class \code{madata}. Read
    \code{\link[maanova]{createData}} for details.}
  \item{arrays}{A vector specifying which arrays to keep or discard.}
  \item{genes}{A vector specifying which genes to keep or discard.}
  \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{madata} with specified arrays and genes.
}

\examples{
data(kidney)
# create data object with replicate
data <- createData(kidney.raw)
# take out array 1 and 2
smalldata <- subset(data, arrays=c(1,2))
# take out the all arrays except array 1
idx.array <- 1:data$n.array
smalldata <- subset(data,arrays=(idx.array[-1]))
# take out gene number 1 to 20
smalldata <- subset(data,genes=1:20)
}

\keyword{utilities}
