\author{Hao Wu \email{hao@jax.org} }

\name{adjPval}
\alias{adjPval}

\title{Generate FDR adjusted P values for F test result.}

\description{
  This function takes a result object from \code{\link[maanova]{matest}}
  and calcualte the FDR adjusted P values. The new P values will be
  appended to the input object as additional fields.
}

\usage{
adjPval(matestobj, method=c("stepup","adaptive", "stepdown"))
}

\arguments{
  \item{matestobj}{An object of class \code{matest}, which is the result
    from \code{\link[maanova]{matest}}.}
  \item{method}{The method for FDR control.}
}

\value{
  An object of class \code{matest} with the following fields added for
  each F test:
  \item{adjPtab}{FDR adjusted tabulated P values.}
  \item{adjPvalperm}{FDR adjusted permutation P values.}
}

\examples{
data(paigen)
paigen <- createData(paigen.raw, n.rep=2)
model.noint.fix <- makeModel(data=paigen, formula=~Array+Dye+Spot+Strain+Diet)
# F-test strain effect
\dontrun{test.strain.fix <- matest(paigen, model.noint.fix, term="Strain", n.perm=100,
        shuffle.method="resid", test.method=rep(1,4))
# make FDR adjusted P values
test.strain.fix <- adjPval(test.strain.fix)}
# there will be new fields in test.strain.fix after this
}

\keyword{utilities}
