\author{Hao Wu}

\name{varplot}
\alias{varplot}

\title{Variance component plot}

\description{
  This function plots the density curve of each variance component of a
  result from \code{\link[maanova]{fitmaanova}}.

  If the input is from fixed model ANOVA, it will plot one curve for
  error variance component. If the input is from mixed model ANOVA, it
  will plot multiple curves, one for a random term (including error). 
}

\usage{
varplot(anovaobj, xlab, ylab, main)
}

\arguments{
  \item{anovaobj}{An object of class \code{maanova}.}
  \item{xlab}{Figure xlab. Default is "Sigma".}
  \item{ylab}{Figure ylab. Default is "Density".}
  \item{main}{Figure title. Default is "Density plot for sqrt of variance".}
}

\examples{
# load abf1 data
data(abf1)
\dontrun{
fit.full.mix <- fitmaanova(abf1, formula = ~Strain+Sample, 
   random = ~Sample)
varplot(fit.full.mix)}
}

\seealso{
  \code{\link[maanova]{fitmaanova}},
  \code{\link[base]{density}}
}

\keyword{dplot}
