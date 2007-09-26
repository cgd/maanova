\author{Hao Wu}

\name{write.madata}
\alias{write.madata}

\title{Write Micro Array data to a TAB delimited simple text file}

\description{
  This function is used to write the contents of an object of class
  \code{madata} to a TAB delimited simple text file.
}

\usage{
write.madata(madata, datafile="madata.txt", designfile="design.txt") 
}

\arguments{
  \item{madata}{The object to be output. It must be an object of class
    \code{madata}.}
  \item{datafile}{The output file name for the data.}
  \item{designfile}{The output file name for the design file.}
}

\examples{
# load abf1 data
data(abf1)
# take out first 6 arrays
\dontrun{
smalldata <- subset(abf1, array=1:6)
# write to file
write.madata(smalldata, datafile="smallabf1.txt",
      designfile="smallabf1design.txt")}
}
\keyword{IO}