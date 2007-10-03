\author{Hao Wu}
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
# load data
data(abf1)
\dontrun{
# fit the anova model
fit.fix = fitmaanova(abf1,formula = ~Strain)
# test Strain effect 
test.fix = matest(abf1, fit.fix, term="Strain",n.perm= 1000)
# pick significant genes - pick the genes selected by Fs test
idx <- volcano(test.fix)$idx.Fs
# do k-means cluster on genes
gene.cluster <- macluster(fit.fix, term="Strain", idx, what="gene", 
   method="kmean", kmean.ngroups=5, n.perm=100)
# get the consensus group
genegroup = consensus(gene.cluster, 0.5)
# get the gene names belonging to each group
genegroupname = genegroup$groupname

# HC cluster on samples
sample.cluster <- macluster(fit.fix, term="Strain", idx, what="sample",method="hc")
# get the consensus group
consensus(sample.cluster, 0.5)
}}
\seealso{
  \code{\link[maanova]{macluster}}
}
\keyword{cluster}

