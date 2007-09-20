\author{Hao Wu}

\name{dyeswapfilter}
\alias{dyeswapfilter}

\title{Gene filter for dye-swap experiment}

\description{
  This function is used to flag the questionable spot in any kind of
  dye-swap experiment.

  This function only works for 2-dye arrays.
}

\details{
  For each pair of dye-swap, the difference in log ratios (d) are
  computed. Then compute the IQR (interquartile range) of d and convert
  that to Standard Deviation by SD = IQR/1.35. Any gene with d larger
  than r times SD will be flagged.

  Note that I assume in the input data object, the adjacent arrays is a
  dye-swap pair. 
}
  
\usage{
dyeswapfilter(dataobj, r=4)
}

\arguments{
  \item{dataobj}{An object of class \code{madata}.}
  \item{r}{A cut-off value for bad spot. The genes with log-ratio
    difference larger than r times standard deviation will be flagged.}
}

\value{
  An object of class \code{rawdata} or \code{madata} with the
  \code{flag} field created or updated.
}

\examples{
\dontrun{
data(kidney)
# riplot before filtering
riplot(kidney.raw, array=1)
# filter the gene
rawdata <- dyeswapfilter(kidney.raw)
# riplot again - some genes are highlighted
riplot(rawdata, array=1)}
}

\keyword{models}
