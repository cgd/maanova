\author{Hao Wu}
\name{fitmaanova}
\alias{fitmaanova}
\title{Fit ANOVA model for Micro Array experiment}
\description{
  This is the function to fit the ANOVA model for Microarray
  experiment. Given the data and formula, this function fits the
  regression model for each gene and calculates the ANOVA estimates, variance components
  for random terms, fitted values, etc.
  For a mixed effect models, the output estimates will be BLUE and BLUP.

  All terms used in the formula should be corresponding to the factor names in designfile except
  "Spot" and "Label". "Spot" represents the spotting effect and "Label"
  represents the labeling effects. They are from the within slide
  technical replicates. If there is no replicated spots, These two terms
  cannot be fitted. Also these two terms cannot be fitted for one-dye
  system (e.g., Affymetrix arrays). (Note that Dye effect should not be
  fitted in one-dye system). 

  A typical formula will be like "~Array+Dye+Sample", which means you
  want to fit Array, Dye and Sample effect in the ANOVA model. In this case, you need to have Array, Dye and Sample columns in
  your input design file. Make sure you have enough degree of
  freedom when making a model. Also you need to be careful about
  confounding problem.
  
  If you have multiple factors in your experiment, you can specify the main
  and interaction effect in the formula. At this time, only two-way
  interactions are allowed.

  When you have random or covariate effect they should be specified in the 'random' and 'covariate', and
  also in the formula.

  For most mixed effect models, Array should be treated as random
  factor. Sample should be treated as random if you have biological
  replicates. Note that the reference sample (0's in Sample) will
  always be treated as fixed even if you specify Sample as random.

  Note that the calculation could be very slow for mixed effect
  models. The computational time depends on the number of genes, number
  of arrays and the size of the random variables (dimension of Z
  matrix). 

  Array specific covariate should be included in the design matrix, and gene specific covariate 
  should be read by 'covM' in \code{read.madata()}, and need to be specified in covariate term. 
}

\usage{
fitmaanova(madata, formula, random= ~1, covariate = ~1, mamodel,
           inits20,method=c("REML","ML","MINQE-I","MINQE-UI", "noest"),
           verbose=TRUE, subCol=FALSE)
}

\arguments{
  \item{madata}{An object of class \code{madata}.}
  \item{formula}{The ANOVA model formula.}
  \item{random}{The formula for random terms. ~1 means only the residual
    is random (fixed model). Note that all random terms should be in the
    ANOVA model formula.}
  \item{covariate}{The formula for covariates. ~1 means no covariates.
    The array specific covariates should be numeric values in the design matrix, 
    and the gene specific covariates should be read by \code{covM} in  \code{\link[base]{read.madata}}}.
  \item{mamodel}{Inside arguments to save the calculation time.}
  \item{inits20}{The initial value for variance components. This should
    be a matrix with number of rows equals to the number of genes and
    number of columns equals to the number of random terms in the model.
    Good initial values will greatly speed up the calculation. If it is not
    given, it will be calculated based on the corresponding fixed model.}
  \item{method}{The method used to solve the Mixed Model
    Equation. Available options includes: "ML" for maximum liklihood;
    "REML" for restricted maximum liklihood; "MINQE-I" and "MINQE-UI"
    are for minimum norm and "noest" for no estimate
    for variance component (use the initial value). Both "ML" and "REML"
    use method of scoring algorithm to solve MME iteratively. "noest" skips the
    iteration and will be significantly faster (but accurate). Default method is "REML". For details about fitting mixed
    effects models, read the "Fitting mixed Effects model" section.
  } 
  \item{verbose}{A logical value to indicate whether to display some
    message for calculation progress.}
  \item{subCol}{A logical value to indicate whether subtracting column mean from the raw data or not. Default is not subtracting column mean but for two color array it automatically subtracts the column mean.}
}
\section{Fitting mixed Effects model}{
  Fitting mixed effects models needs a lot of computation. A good
  starting value for the variances is very important. This function
  first treats all random factors as fixed and fits a fixed effects
  model. Then variances for random factors are calculated and used as
  the initial values for mixed effects model fitting. 

  There are several methods available for fitting the mixed effects
  model. "noest" does not really fit the mixed effects model. It takes
  the initial variance and solve mixed model equations to get the
  estimates (BLUE and BLUP). "MINQE-I" and "MINQE-UI" are based on
  minimum norm unbiased estimators. It is can be thought as a first
  iterate solution of "ML" and "REML", respectively. "ML" and "REML" are
  based on maximum likelihood and restricted maximum likelihood. Both of
  them need to be solved iteratively so they are very slow to
  compute. For "ML" and "REML", a MINQUE estimates is used as
  the starting value.  "Method of scoring" is used as the iteratively
  algorithm to solve ML and REML. "Method of scoring" algorithm is
  similar to New-Raphson method except that it uses the expected value
  of Hessian (second derivative matrix of the objective function)
  instead of Hessian itself. Method of scoring is more robust to poor
  starting values and the Hessian is easier to calculate than
  Newton-Raphson.

  For more mathematical details please read Searle et al. 
}

\value{
  It returns \code{anova} and \code{anova.subcol}. Depending on 'subCol'
  option, one field may not contain any information. Still it needs two fields
  to calculate Fss test statistics. \code{anova} and \code{anova.subcol} contains the same following fields.
  \item{yhat}{Fitted intensity value which has the same dimension as the input
    intensity data}
  \item{S2}{Variance components for the random terms. It is a matrix
    with number of rows equals to the number of genes and number of
    columns equals to the number of random terms. Note that for fixed
    effect model, S2 is a one column vector for error's variance.}
  \item{G}{Gene effects. A vector with the same length as the number
    of genes.}
  \item{reference}{The estimates for reference sample. If there is no
    reference sample specified in the design, this field will be absent
    in the output object.}
  \item{S2.level}{A list of strings to indicate the order of the S2
    field. Note that the last column of S2 is always the error's
    variance. S2.level is only for the non-error terms.
    For example, if there are three columns in S2 and S2.level is
    c("Strain", "Diet"), then the three columns of S2
    correspond to the variances of Strain, Diet and error
    respectively for each gene. }
  \item{Others}{Estimates (or BLUE/BLUP for mixed effect model) for
    the terms in model. There will be XXX.level field for each term
    representing the order of the estimates (similar to S2.level).}
  \item{flag}{A vector to indicate whether there is bad spot for this 
    gene. 0 means no bad spot and 1 means has bad spot. If there is no
    flag information in input data, this field will not be available.}
  \item{model}{The model object used for this fitting.}
}

\examples{
###################################
# fixed model fitting
###################################
# load in abf1 data
data(abf1)
\dontrun{

# fit model with random effect
fit.full.mix <- fitmaanova(abf1, formula = ~Strain+Sample, 
   random = ~Sample)

# this is to explain the usage of including covariate variable.
# .CEL file is not included in the package, thus use can not use this. 
# array specific covariate : add it to the design matrix 
beforeRma <- ReadAffy() # suppose there are 18 arrays.
rmaData <- rma(beforeRma)
datafile <- exprs(rmaData)
design.table=data.frame(Array=row.names(pData(beforeRma)))
Strain = rep(c('Aj', 'B6', 'B6xAJ'), each=6) 
Sample = rep(c(1:9), each=2) 
Cov1 = sample(1:100,18) # this is artificial example 
designfile.cov1 = cbind(design.table, Strain, Sample,Cov1) 
data.cov1=read.madata(datafile, designfile=designfile.cov1) 
fit.cov1 = fitmaanova(data.cov1,formula = ~Strain+Sample+Cov1, covariate = ~ Cov1) 

# gene specific covariate - make artificial 'covM' matrix 
covm = matrix(rnorm(length(datafile)), nrow=nrow(datafile)) 
designfile.cov2 = cbind(design.table, Strain, Sample) 
data.cov2=read.madata(datafile, designfile=designfile.cov2, covM=covm) 
fit.cov2 = fitmaanova(data.cov2,formula = ~Strain+Sample+covM, covariate = ~ covM) }
}

\references{
  Kerr and Churchill(2001), Statistical design and the analysis of gene
  expression microarrays, \emph{Genetical Research}, \bold{77}:123-128.
  
  Kerr, Martin and Churchill(2000), Analysis of variance for gene expression
  microarray data, \emph{Journal of Computational Biology},
  \bold{7}:819-837.

  Searle, Casella and McCulloch, \emph{Variance Components}, John Wiley
  and sons, Inc. 
}

\seealso{
  \code{\link[maanova]{makeModel}}, \code{\link[maanova]{matest}}
}

\keyword{models}

