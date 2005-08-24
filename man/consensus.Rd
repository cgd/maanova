\author{Hao Wu \email{hao@jax.org} }

\name{consensus}
\alias{consensus}

\title{Build consensus tree out of bootstrap cluster result}

\description{
  This is the function to build the consensus tree from the bootstrap
  clustering analysis. If the clustering algorithm is hierarchical
  clustering, the majority rule consensus tree will be built based on
  the given significance level. If the clustering algorithm is K-means, a
  consensus K-means group will be built. 
}

\usage{
consensus(macluster, level = 0.8, draw=TRUE)
}

\arguments{
  \item{macluster}{An object of class \code{macluster}, which is the
    output of \code{\link[maanova]{macluster}}}.
  \item{level}{The significance level for the consensus tree. This is a
    numeric number between 0.5 and 1.}
  \item{draw}{A logical value to indicate whether to draw the consensus
    tree on screen or not.}
}

\value{
  An object of class \code{consensus.hc} or \code{consensus.kmean}
  according to the clustering method.
}

\examples{
# load in data
data(paigen)
# make data object with rep 2
paigen <- createData(paigen.raw, 2)
# make interactive model
model.int.fix <- makeModel(data=paigen,
      formula=~Dye+Array+Strain+Diet+Strain:Diet)
# fit ANOVA model
anova.int <- fitmaanova(paigen, model.int.fix)
# test interaction effect 
\dontrun{test.int.fix <- matest(paigen, model.int.fix, term="Strain:Diet", n.perm=100)
# pick significant genes - pick the genes selected by Fs test
idx <- volcano(test.int.fix)$idx.Fs

# do k-means cluster on genes
gene.cluster <- macluster(anova.int, "Strain:Diet", idx, "gene",
    "kmean", kmean.ngroups=5)
# get the consensus group
consensus(gene.cluster, 0.5)

# HC cluster on samples
sample.cluster <- macluster(anova.int, "Strain:Diet", idx, "sample","hc")
# get the consensus group
consensus(sample.cluster, 0.5)}
}

\seealso{
  \code{\link[maanova]{macluster}}
}

\keyword{cluster}

