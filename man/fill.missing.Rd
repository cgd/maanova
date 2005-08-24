\author{Hao Wu \email{hao@jax.org} }

\name{fill.missing}
\alias{fill.missing}

\title{Fill in missing data}

\description{
  This is the function to do missing data imputation.
}

\details{
  This function will take an object of class \code{rawdata} and fill in
  the missing data. Currently only KNN (K nearest neighbour) algorithm
  is implemented. The memory usage is quadratic in the number of
  genes.
}

\usage{
fill.missing(rawdata, method="knn", k=20, dist.method="euclidean")
}

\arguments{
  \item{rawdata}{An object of class \code{rawdata}, which should be the result
    from \code{\link[maanova]{read.madata}}.}
  \item{method}{The method to do missing data imputation. Currently only
    "knn" (K nearest neighbour) is implemented.}
  \item{k}{Number of neighbours used in imputation. Default is 20.}
  \item{dist.method}{The distance measure to be used. See
    \code{\link[stats]{dist}} for detail.}
}

\value{
  An object of class \code{rawdata} with missing data filled in.
}

\examples{
data(abf1)
# randomly generate some missing data 
rawdata <- abf1.raw
ndata <- length(abf1.raw$data)
pct.missing <- 0.05 # 5% missing
idx.missing <- sample(ndata, floor(ndata*pct.missing))
rawdata$data[idx.missing] <- NA
rawdata <- fill.missing(rawdata)
# plot impute data versus original data
plot(rawdata$data[idx.missing], abf1.raw$data[idx.missing])
abline(0,1)
}

\references{
  \item O.Troyanskaya, M. Cantor, G. Sherlock, P. Brown, T. Hastie,
  R. Tibshirani, D. Botstein, & R. B. Altman. Missing Value estimation
  methods for DNA microarrays. Bioinformatics 17(6):520-525, 2001.
}

\keyword{utilities}

