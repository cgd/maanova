\author{Hao Wu}

\name{resiplot}
\alias{resiplot}

\title{Residual plot for Microarray Experiment}

\description{
  This is the function to plot the residuals versus fitted value figure. Two
  channels, e.g., red and green, are drawn in seperate figures.
}

\usage{
resiplot(madata, anovaobj, header, xlab, ylab)
}

\arguments{
  \item{madata}{An object of class \code{madata}.}
  \item{anovaobj}{An object of class \code{maanova}, which is the output
    from \code{\link[maanova]{fitmaanova}}.}
  \item{header}{Optional. The title of the figure. The default figure
    title will be "Residual vs. Yhat plot".}
  \item{xlab}{Optional. The xlab of the figure. The default will be 
     "Yhat".}
  \item{ylab}{Optional. The ylab of the figure. The default will be 
     "Residual".}
}

\examples{
# load abf1 data
data(abf1)
\dontrun{
fit.full.mix <- fitmaanova(abf1, formula = ~Strain+Sample, 
   random = ~Sample)
resiplot(abf1, fit.full.mix)
}}

\keyword{hplot}
