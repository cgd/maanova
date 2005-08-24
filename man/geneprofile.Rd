\author{Hao Wu \email{hao@jax.org} }

\name{geneprofile}
\alias{geneprofile}

\title{Expression plot for selected genes}

\description{
  This function generate a plot with many lines. Each line represents
  a gene. The y-axis is the estimated expression level for the given
  factor from ANOVA model. The x-axis is for the levels of the
  give factor, e.g.,  different strains. 
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
data(paigen)
# make data object with rep 2
paigen <- createData(paigen.raw, 2)
# make an additive model
model.add.fix <- makeModel(data=paigen, formula=~Dye+Array+Strain+Diet)
# fit ANOVA model
anova.add <- fitmaanova(paigen, model.add.fix)
# test strain effect
\dontrun{test.Strain.fix <- matest(paigen, model.add.fix, term="Strain", n.perm=100)
# volcano plot
idx <- volcano(test.Strain.fix)

# do gene profile for the selected genes
geneprofile(anova.add, "Strain", idx$idx.all)}
}

\keyword{hplot}
