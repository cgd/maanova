\author{Hao Wu}

\name{geneprofile}
\alias{geneprofile}

\title{Expression plot for selected genes}

\description{
  This function generate a plot with many lines. Each line represents
  a gene. The y-axis is the estimated expression level for the given
  factor from ANOVA model. The x-axis is for the levels of the
  give factor, e.g., different strains. 
}

\usage{
geneprofile(anovaobj, term, geneidx,
            col="blue", type="b", ylim, xlab, ylab, ...)
}

\arguments{
  \item{anovaobj}{An object of class \code{maanova}. It should be the
    result from \code{\link[maanova]{fitmaanova}}.}
  \item{term}{The terms to be plotted.}
  \item{geneidx}{The index of genes to be plotted.}
  \item{col}{The color to be used in plot.}
  \item{type}{The line type.}
  \item{ylim}{Y-axis limit.}
  \item{xlab}{X-axis label.}
  \item{ylab}{Y-axis label.}
  \item{\dots}{Other parameters to be passed to
    \code{\link[base]{plot}}.}
}

\examples{
# load in data
data(abf1)
# fit the anova model
\dontrun{
fit.fix = fitmaanova(abf1,formula = ~Strain)
# test Strain effect 
test.fix = matest(abf1, fit.fix, term="Strain",n.perm= 1000)
# pick significant genes - pick the genes selected by Fs test
idx <- volcano(test.fix)$idx.Fs
geneprofile(fit.fix, "Strain", idx)}
}

\keyword{hplot}
