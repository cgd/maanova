\author{Hao Wu \email{hao@jax.org} }

\name{createData}
\alias{createData}

\title{Calculate and create a data object for Micro Array experiment}

\description{
  This is the function to create a \code{madata} object based on the
  given \code{rawdata} and some parameters. 
}

\details{
  The data integrity is checked before making the object. The number of
  rows for the data must be consistent with the number of replicates; the
  number of columns for the data must be consistent with the number of
  dyes, etc.

  Users have the option to collapse the replicated spots. Collapsing
  will be done by taking mean or median of the intensity values for the
  replicated spots. This function assumes the input data is on raw scale.
  So if your data is pre-transformed and on a log2 based,
  collapsing could be wrong.  Also note that once
  the replicated spots are collapsed, you will lose the grid location
  and spatial loess ("rlowess" option in
  \code{\link[maanova]{transform.madata}} will be unapplicable for it.
}

\usage{
createData(rawdata, n.rep=1, avgreps=0, log.trans=TRUE)
}

\arguments{
  \item{rawdata}{An object of class \code{rawdata}, which should be the result
    from \code{\link[maanova]{read.madata}}.}
  \item{n.rep}{An integer to represent the number of replicates.}
  \item{avgreps}{An integer to indicate whether to average the
    replicates or not. 0 means no average; 1 means to take the mean of the
    replicates; 2 means to take the median of the replicates.}
  \item{log.trans}{A logical value to indicate whether to take log2
    transformation on the raw data or not. It is TRUE by default. But
    in the case that your data is pre-transformed, you need to set it to
    FALSE. If this is TRUE, \code{TransformMethod} field will be set to
    "log2".} 
}

\value{
  An object of class \code{madata}, which is a list of following components:
  \item{n.gene}{Total number of genes in the experiment.}
  \item{n.rep}{Number of replicates in the experiment.}
  \item{n.spot}{Number of spots for each gene.}
  \item{data}{data field. It is either the log2 transformed
    rawdata (if log.trans=TRUE), or just the rawdata (if
    log.trans=FALSE).}
  \item{Others}{All other fields in the input object of class \code{rawdata}.}
}

\examples{
#################
# 2-dye arrays
#################
data(paigen)
# create data object with replicate
data2 <- createData(paigen.raw, n.rep=2)
# summarize the data object
summary(data2)
# create data with averaging of replicates
data1 <- createData(paigen.raw, n.rep=2, avgreps=1)
summary(data1)

####################################################################
# affy array - data is pre-transformed so log2 is skipped
####################################################################
data(abf1)
abf1 <- createData(abf1.raw, n.rep=1, log.trans=FALSE)
summary(abf1)
}

\keyword{utilities}

