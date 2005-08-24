\author{Hao Wu \email{hao@jax.org} }

\name{makeModel}
\alias{makeModel}

\title{Make model object for N-dye Micro Array experiment}

\description{
  This is the function to make an object of class \code{mamodel} 
  for a Micro Array experiment. 
}

\details{
  The user needs to specify the ANOVA model by formula. It can be a fixed or
  mixed effect model. This function will check the validity of the data,
  calculate some parameters, construct and design matrices and wrap up
  everything together to create an output object.

  The model formula is for a gene-specific model. All terms in the
  formula should be corresponding to the factor names in design except
  "Spot" and "Label". "Spot" represents the spotting effect and "Label"
  represents the labelling effects. They are from the within slide
  technical replicates. If there is no replicated spots, These two terms
  cannot be fitted. Also these two terms cannot be fitted for one-dye
  system (e.g., affymetric arrays). (Note that Dye effect should not be
  fitted in one-dye system).

  A typical formual will be like "~Array+Dye+Sample", which means you
  want to fit Array effect, Dye effect and Sample effect in the ANOVA
  model. In this case, you need to have Array, Dye and Sample columns in
  your input design file. Make sure you have enough degree of
  freedom when making a model. Also you need to be careful about
  confounding problem.
  
  If you have multiple factors in your experiment, you can specify the main
  and interaction effect in the formula. At this time, only two-way
  interactions are allowed.

  For most mixed effect models, Array should be treated as random
  factor. Sample should be treated as random if you have biological
  replicates. Note that the reference sample (0's in Sample) will
  always be treated as fixed even if you specify Sample as random.
}

\usage{
makeModel(data, design, formula, random=~1, covariate=~1)
}

\arguments{
  \item{data}{An object of class \code{madata}.}
  \item{design}{A data frame representing the experimental design. By
    default, it is a field in madata. But you can always make a data
    frame and pass it to the function.}
  \item{formula}{The ANOVA model formula.}
  \item{random}{The formula for random terms. ~1 means only the residual
    is random (fixed model). Note that all random terms should be in the
    ANOVA model formula.}
  \item{covariate}{The formula for covariates. ~1 means no covariates.
    The covariates will be continuous values in the design matrix.}
}

\value{
  An object of class \code{mamodel} with the following fields:
  \item{X}{Design matrix for fixed terms.}
  \item{dimX}{Number of columns in X for each fixed term.}
  \item{Z}{Design matrix for random terms. This will be absent for fixed
  model.}
  \item{dimZ}{Number of columns in Z for each random term.This will be
    absent for fixed model.}
  \item{df}{The degree of freedom for each term in the model.}
  \item{mixed}{An integer to indicate whether this is a fixed or mixed
    effect model. 0 means fixed and 1 means mixed.}
  \item{design}{The input experimental design as a data frame.}
  \item{formula}{The input model formula.}
  \item{random}{The input formula for random terms.}
  \item{covariate}{The input formula for covariates.}
}

\examples{
#load in data
data(paigen)
# make data object with rep 2
paigen <- createData(paigen.raw, 2)
# make full model for fixed effect model
model.full.fix <- makeModel(data=paigen,
      formula=~Dye+Array+Spot+Strain+Diet+Strain:Diet)
summary(model.full.fix)
# make full model for mixed effect model
model.full.mix <- makeModel(data=paigen,
      formula=~Dye+Array+Spot+Strain+Diet+Strain:Diet+Sample,
      random=~Array+Spot+Sample)
summary(model.full.mix)
}

\references{
  Kerr and Churchill(2001), Statistical design and the analysis of gene
  expression microarrays, \emph{Genetical Research}, \bold{77}:123-128.
  
  Kerr, Martin and Churchill(2000), Analysis of variance for gene expression
  microarray data, \emph{Journal of Computational Biology},
  \bold{7}:819-837.
}

\keyword{models}
