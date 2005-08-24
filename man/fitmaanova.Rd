\author{Hao Wu \email{hao@jax.org} }

\name{fitmaanova}
\alias{fitmaanova}

\title{Fit ANOVA model for Micro Array experiment}

\description{
  This is the function to fit the ANOVA model for Micro Array
  experiment. Given the data and model object, this function will fit the
  regression gene by gene and output the estimates, variance components
  for random terms, fitted values, etc.
  For a mixed effect models, the output estimates will be BLUE and BLUP.

  Note that the calculation could be very slow for mixed effect
  models. The computational time depends on the number of genes, number
  of arrays and the size of the random variables (dimension of Z
  matrix). 
}

\usage{
fitmaanova(madata, mamodel, inits20,
           method=c("REML","ML","MINQE-I","MINQE-UI", "noest"),
           verbose=TRUE)
}

\arguments{
  \item{madata}{An object of class \code{madata}.}
  \item{mamodel}{An object of class \code{mamodel}.}
  \item{inits20}{The initial value for variance components. This should
    be a matrix with number of rows equals to the number of genes and
    number of columns equals to the number of random terms in the model.
    Good initial values will greatly speed up the calculation. If not
    given, it will be calculated based on the corresponding fixed
    model.}
  \item{method}{The method used to solve the Mixed Model
    Equation. Available options includes: "ML" for maximum liklihood;
    "REML" for restricted maximum liklihood; "MINQE-I" and "MINQE-UI"
    are for minimum norm and "noest" for no estimate
    for variance component (use the initial value). Both "ML" and "REML"
    use method of scoring algorithm to solve MME iteratively. "noest" skips the
    iteration and will be significantly faster (but 
    accurate). Default method is "REML". For details about fitting mixed
    effects models, read the "Fitting mixed Effects model" section.
  } 
  \item{verbose}{A logical value to indicate whether to display some
    message for calculation progress.}
}

\section{Fitting mixed Effects model}{
  Fitting mixed effects models needs a lot of computation. A good
  starting value for the variances is very important. This function
  first treats all random factors as fixed and fits a fixed effects
  model. Then variances for random factors are calculated and used as
  the initial values for mixed effects model fitting. 

  There are several methods available for fitting the mixed effects
  model. "noest" does not really fit the mixed effects model. It takes
  the initial variance and solve mixed model equatinos to get the
  estimates (BLUE and BLUP). "MINQE-I" and "MINQE-UI" are based on
  minimum norm unbiased estimators. It is can be thought as a first
  iterate solution of "ML" and "REML", respectively. "ML" and "REML" are
  based on maximum likelihood and restricted maximum likelihood. Both of
  them need to be solved iteratively so they are verly slow to
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
  An object of class \code{maanova}, which is a list of following
  components:
  \item{yhat}{Fitted pmt value which has the same dimension as the input
    pmt data}
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
# load in Paigen's data
data(paigen)
# make data object with rep 2
paigen <- createData(paigen.raw, 2)
# Note that the data is normalized so normalization is skipped
# full model
model.full.fix <- makeModel(data=paigen,
      formula=~Dye+Array+Spot+Strain+Diet+Strain:Diet)
anova.full.fix <- fitmaanova(paigen, model.full.fix)
# residual plot
resiplot(paigen, anova.full.fix)

#######################################
# mixed model fitting -
# Array, Spot and biological
# replicates are random effects. 
# This may take a while to finish
#######################################
\dontrun{model.full.mix <- makeModel(data=paigen,
      formula=~Dye+Array+Spot+Strain+Diet+Strain:Diet+Sample,
      random=~Array+Spot+Sample)
anova.full.mix <- fitmaanova(paigen, model.full.mix, method="REML")
# residual plot
resiplot(paigen, anova.full.mix)
# variance component plot
varplot(anova.full.mix)}
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

