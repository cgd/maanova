\author{Hao Wu}

\name{fill.missing}
\alias{fill.missing}

\title{Fill in missing data}

\description{
  This is the function to do missing data imputation.
}

\details{
  This function will take an object of class \code{madata} and fill in
  the missing data. Currently only KNN (K nearest neighbour) algorithm
  is implemented. The memory usage is quadratic in the number of genes.
}

\usage{
fill.missing(data, method="knn", k=20, dist.method="euclidean")
}

\arguments{
  \item{data}{An object of class \code{madata}, which should be the result
    from \code{\link[maanova]{read.madata}}.}
  \item{method}{The method to do missing data imputation. Currently only
    "knn" (K nearest neighbour) is implemented.}
  \item{k}{Number of neighbours used in imputation. Default is 20.}
  \item{dist.method}{The distance measure to be used. See
    \code{\link[stats]{dist}} for detail.}
}

\value{
  An object of class \code{madata} with missing data filled in.
}

\examples{
data(abf1)
# randomly generate some missing data 
rawdata <- abf1
ndata <- length(abf1$data)
pct.missing <- 0.05 # 5% missing
idx.missing <- sample(ndata, floor(ndata*pct.missing))
rawdata$data[idx.missing] <- NA
rawdata <- fill.missing(rawdata)
# plot impute data versus original data
plot(rawdata$data[idx.missing], abf1$data[idx.missing])
abline(0,1)
}

\references{
  O.Troyanskaya, M. Cantor, G. Sherlock, P. Brown, T. Hastie,
  R. Tibshirani, D. Botstein, & R. B. Altman. Missing Value estimation
  methods for DNA microarrays. Bioinformatics 17(6):520-525, 2001.
}
\keyword{utilities}

