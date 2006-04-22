\author{Hao Wu}

\name{exprSet2Rawdata}
\alias{exprSet2Rawdata}

\title{Convert an object of exprSet to an object of Rawdata}

\description{
  This function converts an object of exprSet class, which is the main
  class for microarray data into an object of class "rawdata". It serves
  as a bridge between BioConductor and R/maanova. 
}

\usage{
exprSet2Rawdata(exprdata, ndye=1, trans.method="None") 
}

\arguments{
  \item{exprdata}{An object of class \code{\link[Biobase]{exprSet}}.}
  \item{ndye}{Number of dyes used in the experiment. Default is 1.}
  \item{trans.method}{A string for data transformation method
    used. Default is "None". This field is used only for display
    purpose. I suggest user to put something here as a reminder to
    yourself. For example if you run justGCRMA to get the exprdata, put
    "GCRMA" here.} 
}

\value{
  An object of class \code{rawdata}.
}


\keyword{IO}
