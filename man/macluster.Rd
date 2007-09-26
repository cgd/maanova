\author{Hao Wu}

\name{macluster}
\alias{macluster}

\title{Clustering analysis for Microarray experiment }

\description{
  This function bootstraps K-means or hierarchical clusters and builds a 
  consensus tree (consensus group for K-means) from the bootstrap
  result.
}

\details{
  Normally after the F test, user can select a list of differentially expressed
  genes. The next step is to investigate the relationship among these
  genes. Using the expression levels of these genes, the user can cluster the
  genes or the samples using either hierarchical or K-means clustering
  algorithm. In order to evaluate the stability of the relationship,
  this function bootstraps the data, re-fits the model and recluster the
  genes/samples. Then for a certain number of bootstrap iterations, say,
  1000, we have 1000 cluster results. We can use
  \code{\link[maanova]{consensus}} to  build the consensus tree from
  these 1000 trees. 

  Note that if you have a large number (say, more than 100) of
  genes/samples to cluster, hierarchical clustering could be very
  unstable. A slight change in the data can result in a big change in the
  tree structure. In that case, K-means will give better results.
}

  
\usage{
macluster(anovaobj, term, idx.gene, what = c("gene", "sample"), 
    method = c("hc", "kmean"), dist.method = "correlation",
    hc.method = "ward", kmean.ngroups, n.perm = 100)
}

\arguments{
  \item{anovaobj}{The result object for fitting ANOVA model.}
  \item{term}{The factor (in formula) used in clustering. The expression
    level for this term will be used in clustering. This term has to 
    correspond to the gene list, e.g, idx.gene in this function. The
    gene list should be the significant hits in testing this term.}
  \item{idx.gene}{A vector indicating the list of differentially expressed
    genes. The expression level of these genes will be used to construct
    the cluster.}
  \item{what}{What to be clustered, either gene or sample.}
  \item{method}{The clustering method. Right now hierarchical clustering
    ("hc") and K-means ("kmean") are available.}
  \item{dist.method}{Distance measure to be used in hierarchical
    clustering. Besides the methods listed in \code{\link[mva]{dist}},
    there is a new method "correlation" (default). The "correlation"
    distance equals to (1 - $r^2$), where r is the sample correlation
    between observations. }
  \item{hc.method}{The agglomeration method to be used in hierarchical
    clustering. See \code{\link[mva]{hclust}} for detail.}
  \item{kmean.ngroups}{The number of groups for K-means cluster.}
  \item{n.perm}{Number of bootstraps. If it is 1, this function will
    cluster the observed data. If it is bigger than 1, a bootstrap
    will be performed.}
}

\value{
  An object of class \code{macluster}.
}

\seealso{
  \code{\link[mva]{hclust}}, \code{\link[mva]{kmeans}},
  \code{\link[maanova]{consensus}}
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
# do k-means cluster on genes
gene.cluster <- macluster(fit.fix, term="Strain", idx, what="gene", 
   method="kmean", kmean.ngroups=5, n.perm=100)
# get the consensus group
consensus(gene.cluster, 0.5)

# HC cluster on samples
sample.cluster <- macluster(fit.fix, term="Strain", idx, what="sample",method="hc")
# get the consensus group
consensus(sample.cluster, 0.5)}
}

\keyword{cluster}
  
