\author{Hao Wu}

\name{fom}
\alias{fom}

\title{Figure of Merit}

\description{
  K-means clustering needs a given number of groups, which is difficult to
  guess in most of the cases. This function calculates the Figure of
  Merit values for different number of groups and generates the FOM
  plot (FOM value versus number of groups). Lower FOM value means
  better grouping. User can decide the number of groups in kmeans
  cluster based on that result.
}

\usage{
fom(anovaobj, idx.gene, term, ngroups)
}

\arguments{
  \item{anovaobj}{An object of class \code{maanova}.}
  \item{idx.gene}{The index of genes to be clustered.}
  \item{term}{The factor (in formula) used in clustering. The expression
    level for this term will be used in clustering. This term has to 
    correspond to the gene list, e.g, idx.gene in this function. The
    gene list should be the significant hits in testing this term.}
  \item{ngroups}{The number of groups for K-means cluster. This could be
    a vector or an integer.}
}

\value{
  A vector of FOM values for the given number of groups
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
# generate FOM
m <- fom(fit.fix, idx, "Strain", 10)}
}

\references{
  Yeung, K.Y., D.R. Haynor, and W.L.Ruzzo (2001). Validating clustering
  for gene expression data. \emph{Bioinformatics}, \bold{17:309-318}.
}

\seealso{
  \code{\link[maanova]{macluster}},
  \code{\link[maanova]{consensus}},
  \code{\link[stats]{kmeans}}
}

\keyword{models}

