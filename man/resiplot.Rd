\author{Hao Wu \email{hao@jax.org} }

\name{resiplot}
\alias{resiplot}

\title{Residual plot for Micro Array Experiment}

\description{
  This is the function to plot the residuals versus fitted value figure. Two
  channels, e.g., red and green, are drawn in seperate figures.
}

\usage{
resiplot(madata, anovaobj, header)
}

\arguments{
  \item{madata}{An object of class \code{madata}.}
  \item{anovaobj}{An object of class \code{maanova}, which is the output
    from \code{\link[maanova]{fitmaanova}}.}
  \item{header}{Optional. The title of the figure. The default figure
    title will be "Residual vs. Yhat plot".}
}

\examples{
# load in Paigen's data
data(paigen)
# make data object with rep 2
paigen <- createData(paigen.raw, 2)
# Note that the data is normalized so normalization is skipped
# full model
model.full.fix <- makeModel(data=paigen,
      formula=~Dye+Array+Spot+Strain+Diet+Strain:Diet)
anova.full.fix <- fitmaanova(paigen, model.full.fix)
# residual plot
resiplot(paigen, anova.full.fix)
}

\keyword{hplot}
