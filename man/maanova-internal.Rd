\author{Hao Wu;  Hyuna Yang, \email{hyuna.yang@jax.org}}

\name{maanova-internal}

\alias{JS}
\alias{JSshrinker}
\alias{buildtree}
\alias{calPval}
\alias{calVolcanoXval}
\alias{caldf}
\alias{check.confounding}
\alias{checkContrast}
\alias{cluster2num}
\alias{consensus.hc}
\alias{consensus.kmean}
\alias{dist.cor}
\alias{findgroup}
\alias{getPval.volcano}
\alias{glowess}
\alias{intprod}
\alias{linlog}
\alias{linlog.engine}
\alias{linlogshift}
\alias{locateTerm}
\alias{make.ratio}
\alias{makeAB}
\alias{makeCompMat}
\alias{makeD}
\alias{makeDesign}
\alias{makeHq}
\alias{makeShuffleGroup}
\alias{makeZiZi}
\alias{makelevel}
\alias{matest.engine}
\alias{matest.perm}
\alias{meanvarlog}
\alias{plot.consensus.hc}
\alias{plot.consensus.kmean}
\alias{print.maanova}
\alias{print.madata}
\alias{print.summary.madata}
\alias{print.summary.mamodel}
\alias{ratioVarplot}
\alias{rlowess}
\alias{shift}
\alias{shuffle.maanova}
\alias{solveMME}
\alias{summary.madata}
\alias{summary.mamodel}
\alias{volcano.ftest}
\alias{volcano.ttest}
\alias{blkdiag}
\alias{colmax}
\alias{colmin}
\alias{makeContrast}
\alias{matrank}
\alias{matsort}
\alias{mixed}
\alias{norm}
\alias{num2yn}
\alias{ones}
\alias{parseformula}
\alias{pinv}
\alias{repmat}
\alias{rowmax}
\alias{rowmin}
\alias{sumrow}
\alias{zeros}
\alias{fdr}

\usage{
JS(X, var)
JSshrinker(X, df, meanlog, varlog)
buildtree(ct, binstr, depth, parent, idx.node, idx.leave)
calPval(fstar, fobs, pool)
calVolcanoXval(matestobj)
caldf(model, term)
check.confounding(model, term1, term2)
checkContrast(model, term, Contrast)
cluster2num(clust)
consensus.hc(macluster, level, draw)
consensus.kmean(macluster, level, draw)
dist.cor(x)
findgroup(varid, ndye)
getPval.volcano(matestobj, method, idx)
glowess(object, method, f, iter, degree, draw)
intprod(terms, intterm)
linlog(object, cg, cr, draw)
linlog.engine(data, cutoff)
linlogshift(object, lolim, uplim, cg, cr, n.bin, draw)
locateTerm(labels, term)
make.ratio(object, norm.std=TRUE)
makeAB(ct, coord, treeidx, startx, maxdepth)
makeCompMat(n)
makeD(s20, dimZ)
makeDesign(design)
makeHq(s20, y, X, Z, Zi, ZiZi, dim, b, method)
makeShuffleGroup(sample.mtx, ndye, narray)
makeZiZi(Z, dimZ)
makelevel(model, term)
matest.engine(anovaobj, term, mv, test.method, Contrast,
              is.ftest, partC, verbose=FALSE)
matest.perm(n.perm, FobsObj, data, model, term, Contrast, inits20,
            mv, is.ftest, partC, MME.method, test.method,
            shuffle.method, pool.pval, ngenes)
meanvarlog(df)
plot.consensus.hc(x, title, \dots)
plot.consensus.kmean(x, \dots)
\method{print}{maanova}(x, \dots)
\method{print}{madata}(x, \dots)
\method{print}{summary.mamodel}(x, \dots)
ratioVarplot(logsum, logdiff, n)
rlowess(object, method, grow, gcol, f, iter, degree, draw)
shift(object, lolim, uplim, draw)
shuffle.maanova(data, model, term)
solveMME(s20, dim, XX, XZ, ZZ, a)
\method{summary}{madata}(object, \dots)
\method{summary}{mamodel}(object, \dots)
volcano.ftest(matestobj, threshold, method, title,highlight.flag)
volcano.ttest(matestobj, threshold, method, title,highlight.flag,
   onScreen)
matsort(mat, index=1)
repmat(mat, n.row, n.col, \dots)
zeros(dim)
ones(dim)
blkdiag(...)
rowmax(x)
rowmin(x)
colmax(x)
colmin(x)
sumrow(x)
matrank(X)
norm(X)
mixed(y, X, Z, XX, XZ, ZZ, Zi, ZiZi, dimZ, s20, method =
      c("noest", "MINQE-I", "MINQE-UI", "ML", "REML"),
      maxiter = 100)
parseformula(formula, random, covariate)
makeContrast(model, term)
pinv(X, tol)
fdr(p, method = c("stepup", "adaptive", "stepdown"))
}

\title{Internal maanova functions}

\description{
  Internal maanova functions. These are generally not to be called by
  the user. 
}

\details{
  Some funtion descriptions are:
    \itemize{
    \item matsort: Sort matrix in ascending order along specified dimension
    \item repmat: Replicate and tile an array
    \item zeros: Create an array with all zeros
    \item ones: Create an array with all ones
    \item blkdiag: Block diagonal concatenation of input arguments
    \item num2yn: convert a logical value to string "Yes" or "No"
    \item rowmax, rowmin, colmax, colmin: find the maximum/minimum value
    for row/columns
    \item sumrow: calculate the sum of rows for a given matrix
    \item matrank: calculate the rank of a matrix
    \item norm: calculate matrix or vector norm, working only for vector
    now
    \item mixed: engine function to solve Mixed Model Equations using EM
    algorithm
    \item parseformula: parse input formula. This is used for mixed
    effect model
    \item makeDesign: function to make a integer list from input design object
    \item intpord: function to make the design matrix for interaction terms
    it's working for two way interaction only
    \item makeContrast: function to make the contrast matrix given model
    and the term to be tested
    number of levels
    \item pinv: calculate the pseudo inverse for a singular matrix.
    Note that I was using ginv function in MASS but it is not
    robust, e.g., sometimes have no result. That's because
    the engine function dsvdc set the maximum number of iteration
    to be 30, which is not enough in some case.
    I use La.svd instead of svd in my function.
    I don't want to spend time on it so it doesn't support complex
    number
    \item fdr: function to calculate the adjusted P values for FDR control.
  }

}

\examples{
# for matsort
a<-matrix(c(1,6,4,3,5,2),2,3)
matsort(a,1)
matsort(a,2)

# for ones and zeros
ones(c(2,2))
zeros(c(2,3,2))

# for repmat
a<-c(1,2)
repmat(a,2,1)
a<-matrix(1:4,2,2)
repmat(a,1,2)

# for blkdiag
a<-matrix(1:4,2,2)
b<-matrix(3:6,2,2)
blkdiag(a,b)
blkdiag(a,b,c(1,2))

# others examples are omitted

}

\keyword{internal}
