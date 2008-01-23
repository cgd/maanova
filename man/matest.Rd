\author{Hao Wu}

\name{matest}
\alias{matest}

\title{Statistical test for Microarray experiment}

\description{
  This is the function to perform F or T test on one or multiple
  experimental factor(s). Permutation test will be carried upon request.
}

\details{
  If user provide a comparison matrix, this function will perform T-test
  on the given comparison(s). Otherwise, this function will perform
  F-test for the given term.
  
  There are three types of tests available. All three tests are based on the
  gene-specific ANOVA model.  F1 is the usual F statistic, Fs is based on the
  James-Stein shrinkage estimates of the error variance.
  
  Permutation tests can run on MPI cluster. 
  This feature is only available for Unix/Linux system. Several
  other R packages (such like SNOW, Rmpi, etc.) are needed for using
  cluster. You may need help from your system administrator to setup LAM-MPI
  correctly. For detailed information on LAM-MPI cluster setup and the
  cluster usage in R, read "MPI\_README" distributed with the package.
}

\usage{
matest(data, anovaobj, term, Contrast, n.perm=1000, nnodes=1,
       critical=.9, test.type = c("ttest", "ftest"),
       shuffle.method=c("sample", "resid"),
       MME.method=c("REML","noest","ML"),
       test.method=c(1,1),pval.pool=TRUE, verbose=TRUE)
}

\arguments{
  \item{data}{An object of class \code{madata}.}
  \item{anovaobj}{An object of class \code{fitmaanova}.}
  \item{term}{The term(s) to be tested. It can be multiple terms. Note
    that the tested term must be fixed. If the term to be tested is a
    random term, it will be converted to a fixed term than do test.}
  \item{Contrast}{The contrast matrix for the term. The number of
    columns equals to the number of levels in the term. The number of
    rows is the number of T-test you want to carry. Note that it must
    be a matrix. Use \code{\link[maanova]{PairContrast}} to make all possible
  pairwise comparison or \code{\link[base]{matrix}} command to make it
  manually. Note that the the hypothesis test can be formulated as H0: Lb=0
  versus alternative. This contrast matrix is L. For testing a covariate,
    use a one by one contrast matrix of 1.}
  \item{n.perm}{An integer for number of permutations.}
  \item{nnodes}{Number of nodes in the MPI cluster. If 1, the
    permutation test will be running on the local computer.}
  \item{critical}{percentile of F-distribution used to get a subset to 
    calculate p-value. Default is 90th percentile of F-distribution, and 
    permutation analysis is conducted based on genes whose test statistics 
    are smaller than 90th percentile of the F-distribution.}  
  \item{test.type}{Test type. It could be F-test or T-test. If the
    Contrast matrix is missing, this should be a "ftest" and the 
    contrast matrix is generated automatically to cover the whole
    linear space except for testing covariates. If the Contrast matrix 
    is given, this could be "ftest" or "ttest". The default is "ttest" 
    (for backward compatability). For T-test, the code will do a series 
    of T-test, where each T-test corresponds to a row in the contrast matrix.}
  \item{shuffle.method}{Data shuffling method. "sample" for sample
    shuffling and "resid" for residual shuffling. Read "Data Shuffling"
    section for detail.}
  \item{MME.method}{The method used to solve the Mixed Model Equations. See
    \code{\link[maanova]{fitmaanova}} for detail. This parameter only
    applies for mixed effects model permutation test. Default method is
    "REML". The variance components for observed data will be used
    for permuted data. It will greatly increase the speed but you may
    lose power in statistical test in some cases.}
  \item{test.method}{An integer vector of two elements to indicate
    which F test to carry. Default is c(1,1), which means do F1 and Fs test.}
  \item{pval.pool}{A logical value to indicate whether to use pooled
    permutation F values to calculate the P values.}
  \item{verbose}{A logical value to indicate whether to display some
    message for calculation progress.}
}

\section{Data Shuffling}{
  Data shuffling method is a crucial part in the permutation
  test. Currently there are two shuffling method available, residual
  shuffling and sample shuffling. 

  Residual shuffling is to shuffle the null model residuals within gene
  without replacement. 

  Sample shuffling is to shuffle the samples based on the nesting
  relationship among the experimental factors in the model. For sample
  shuffling, you need to make sure you have a 
  good sample size. Otherwise the result may be biased.
}


\value{
  An object of class \code{matest}, which is a list of the following
  components:
  \item{model}{Input model object.}
  \item{term}{The input term(s) to be tested.}
  \item{dfde}{Denominator's degree of freedom for the test.}
  \item{dfnu}{Numerator's degree of freedom for the test. Note that this
    is always 1 for T-test.}
  \item{obsAnova}{An object of \code{fitmaanova}, which is the ANOVA model
    fitting result on the original data.}
  \item{Contrast}{The contrast matrix used in the test.}
  \item{n.perm}{Number of permutations.}
  \item{shuffle}{Shuffle style}
  \item{pval.pool}{Use pooled P value or not.}
  \item{F1, Fs}{Objects of four different F tests results. All or
    any of them could be there according to the requested F test
    method. Each of them contains the following fields:
    \itemize{
      \item{Fobs}{F value for the observed data.}
      \item{Ptab}{Tabulated P values for the observed data.}
      \item{Pvalperm}{Nominal permutation P values for each gene. This
	field will be unavailable if user do not do permutation test.}
      \item{Pvalmax}{FWER one-step adjusted P values from the permutation
	test.}
    }
    All the F values and P values are matrices. The number of rows in
    the matrices equals to the number of genes. For F-test, the number
    of columns will be one. For T-test, the number of columns equals
    to the number of tests carried.
  }
}
\examples{
# load in abf1 data
data(abf1)
\dontrun{
fit.full.mix <- fitmaanova(abf1, formula = ~Strain+Sample, 
   random = ~Sample)
ftest.all = matest(abf1, fit.full.mix, test.method=c(1,1),
    shuffle.method="sample", term="Strain", n.perm= 100)
C = matrix(c(1,-1,0,1,0,-1), ncol=3, byrow=T)
ftest.pair = matest(abf1, fit.full.mix, Contrast = C, 
   term="Strain", n.perm=100)}
}
\keyword{models}

\references{
  Cui, X. and Churchill,GA (2003), Statistical tests for differential
  expression in cDNA Microarray experiments, \emph{Genome Biology} 4:210. 

  Cui, X., Hwang, J.T.G., Blades N., Qiu J. and Churchill GA (2003),
  Improved statistical tests for differential gene expression by
  shrinking variance components, to be submitted.
}

\seealso{
  \code{\link[maanova]{makeModel}}, \code{\link[maanova]{fitmaanova}}
}
