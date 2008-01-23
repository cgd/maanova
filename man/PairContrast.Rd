\author{Hyuna Yang}

\name{PairContrast}
\alias{PairContrast}

\title{Pairwise comparison matrix}

\description{
 It returns all possible pairwise comparison. 
}

\usage{
PairContrast(n)
}

\arguments{
  \item{n}{Number of levels in test term.}
}

\value{
  \item{res}{All pairwise comparison matrix.}
}
\examples{
# load in abf1 data
data(abf1)
\dontrun{
fit.full.mix <- fitmaanova(abf1, formula = ~Strain+Sample, 
   random = ~Sample)
ftest.all = matest(abf1, fit.full.mix, test.method=c(1,1),
    shuffle.method="sample", term="Strain", n.perm= 100)
C = PairContrast(3)
ftest.pair = matest(abf1, fit.full.mix, Contrast = C, 
   term="Strain", n.perm=100)}
}
\keyword{models}
