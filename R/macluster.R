######################################################################
#
# macluster.R
#
# copyright (c) 2001-2002, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written May, 2002
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
#
######################################################################

require("stats")

macluster <-
  function(anovaobj, term, idx.gene, what=c("gene","sample"),
           method=c("hc", "kmean"), dist.method="correlation",
           hc.method="ward", kmean.ngroups, n.perm=100)
{
  what <- match.arg(what)
  method <- match.arg(method)
  
  # input data checking
  if(n.perm<1)
    stop("The number of permutations has to be an integer greater than or equal to 1")
  
  if( ! (term %in% names(anovaobj)) )
    stop(paste(term, "is not in input anova result"))

  # use all genes to do clustering by default
  if(missing(idx.gene))
    idx.gene <- 1:length(anovaobj$G)
  
  # get the data to be clustered
  if(what=="gene")
    data.cluster <- anovaobj[[term]][idx.gene,]
  else
    data.cluster <- t(anovaobj[[term]][idx.gene,])

  # number of levels for the term
  nlevel <- dim(data.cluster)[2]
  
  # initialize output variable
  result <- NULL
  result$what <- what
  result$method <- method
  result$n.perm <- n.perm
  result$term <- term
  result$cluster.obs <- NULL
  result$idx.gene <- idx.gene
  result$VG <- data.cluster
  # make leave names and condition names
  # leave names is for HC and condition name is for kmeans
  if(what=="sample") { # cluster sample
    result$leave.names <- anovaobj[[paste(term,"level",sep=".")]]
    result$condition.names <- paste("gene", idx.gene)
  }
  else { #cluster genes
    result$leave.names <- paste("gene", idx.gene)
    result$condition.names <- anovaobj[[paste(term,"level",sep=".")]]
  }
  
  # do cluster on observed data
  if(method=="hc") { # hierarchical clustering
    if(dist.method == "correlation")
      d <- dist.cor(data.cluster)
    else
      d <- dist(data.cluster, dist.method)
    result$cluster.obs <- hclust(d, hc.method)
  }
  else if(method=="kmean") { # Kmean clustering
    if(kmean.ngroups>dim(data.cluster)[1])
      stop("Number of groups is more than the data points")
    tmp <- kmeans(data.cluster, centers=kmean.ngroups,
                  iter.max=50)
    result$cluster.obs <- tmp$cluster
    result$kmean.ngroups <- kmean.ngroups
    center.obs <- tmp$centers
  }
     
  # set up bootstrap if n.perm is bigger than 1
  if(n.perm > 1) {
    # initialize cluster permutation result
    result$cluster.perm <- list(NULL)
    
    VG0 <- data.cluster

    # bootstrap iteration
    for(i in 2:n.perm) {
      cat(paste("Permutation number", i,"\n"))

      ################ shuffle data ###################
      # shuffle the observed VG matrix columns
      idx.shuffle <- sample(nlevel, replace=TRUE)
      data.cluster <- VG0[, idx.shuffle]

      # also shuffle the columns of center if method is kmean
      if(method == "kmean")
        center.perm <- center.obs[,idx.shuffle]

      ######## start clustering ###############
      # do cluster on shuffled data
      if(method=="hc") { # hierarchical clustering
        if(dist.method == "correlation")
          d <- dist.cor(data.cluster)
        else
          d <- dist(data.cluster, dist.method)
        result$cluster.perm[[i]] <- hclust(d, hc.method)$merge
      }
      else if(method=="kmean") { # Kmean clustering
        tmp <- kmeans(data.cluster, centers=center.perm,
                      iter.max=50)
        result$cluster.perm[[i]] <- tmp$cluster
      }

    }
    #### finish permutation ################

  } ## for(i in 1:n.perm) {

  class(result) <- "macluster"
  invisible(result)
}


dist.cor <- function(x) {
  # number of rows of the input
  N <- nrow(x)
  # take a transpose because cor operates on columns
  # but dist operates on rows
  x <- t(as.matrix(x))
  # compute the correlations
  r <- cor(x)
  # take the lower triangle
  r <- r[lower.tri(r)]
  d <- 1 - r^2

  # return variable
  attr(d, "Size") <- N
  attr(d, "Labels") <- dimnames(x)[[1]]
  attr(d, "Diag") <- FALSE
  attr(d, "Upper") <- FALSE
  attr(d, "method") <- "correlation"
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  return(d)
}

  
  
  
