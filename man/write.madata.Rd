\author{Hao Wu \email{hao@jax.org} }

\name{write.madata}
\alias{write.madata}

\title{Write Micro Array data to a TAB delimited simple text file}

\description{
  This function is used to write the contents of an object of class
  \code{madata} or \code{rawdata} to a TAB delimited simple text file.
}

\usage{
write.madata(madata, datafile="madata.txt", designfile="design.txt") 
}

\arguments{
  \item{madata}{The object to be output. It must be an object of class
    \code{madata} or \code{rawdata}.}
  \item{datafile}{The output file name for the data.}
  \item{designfile}{The output file name for the design file.}
}

\examples{
# load in data
data(paigen)
# make data object
paigen <- createData(paigen.raw, 2)
# take out first 6 arrays
smallpaigen <- subset(paigen, array=1:6)
# write to file
\dontrun{write.madata(smallpaigen, datafile="smallpaigen.txt",
      designfile="smallpaigendesign.txt")}
}

\keyword{IO}




