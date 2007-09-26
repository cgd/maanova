\author{Hao Wu}
\name{adjPval}
\alias{adjPval}
\title{Generate FDR adjusted P values for F test result.}
\description{
  This function takes a result object from \code{\link[maanova]{matest}}
  and calculate the FDR adjusted P values. The new P values will be
  appended to the input object as additional fields. It has four
  options;  "stepup" (Hochberg and Benjamini, 1990), "adaptive" (Benjamini and Hochberg, 2000), "stepdown" (Westfall and Young, 1993) and "jsFDR" (Storey,
  2002). "jsFDR" option uses 'qvalue' package by John Storey and user suppose
  to install 'qvalue' package before using this option. There is no default
  option, thus you need to specify one option. 
}
\usage{
adjPval(matestobj, method=c("stepup","adaptive", "stepdown", "jsFDR"))
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
data(abf1)
\dontrun{ 
fit.full.simple = fitmaanova(abf1,formula = ~Strain)
# F-test strain effect
ftest.all = matest(abf1, fit.full.simple, term="Strain",n.perm= 1000)
# make FDR adjusted P values
ftest.all = adjPval(ftest.all, 'jsFDR')
# there will be new fields in test.strain.fix after this
}
}
\keyword{utilities}
