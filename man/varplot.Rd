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
varplot(anovaobj)
}

\arguments{
  \item{anovaobj}{An object of class \code{maanova}.}
}

\examples{
# load in Paigen's data
data(paigen)
# make data object with rep 2
paigen <- createData(paigen.raw, 2)
\dontrun{model.full.mix <- makeModel(data=paigen,
      formula=~Dye+Array+Spot+Strain+Diet+Strain:Diet+Sample,
      random=~Array+Spot+Sample)
anova.full.mix <- fitmaanova(paigen, model.full.mix, method="REML")
varplot(anova.full.mix)}
}

\seealso{
  \code{\link[maanova]{fitmaanova}},
  \code{\link[base]{density}}
}

\keyword{dplot}
