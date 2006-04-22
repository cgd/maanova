\author{Hao Wu}

\name{gridcheck}
\alias{gridcheck}

\title{Plot grid-by-grid data comparison for arrays}

\description{
  This function is used to check microarray data quality. It can
  check the data within the same array or cross different arrays.

  Normally, on one array, the pmt data for both channels (CY5 and CY3)
  should be highly correlated (also apparent on the RI plot). The pmt data
  for the same sample on different arrays should be highly correlated
  too. Normally if an error happened in gridding, only a few blocks
  will be misgridded. This function does the scatter plot on a grid basis
  to check the quality of hybridazition and gridding. 

  If you only provide array1 (either an integer or a vector), it will do
  grid check within the same array, that is, for each slide, there will
  be one scatter plot for log2(Red) versus log2(Green) for each grid.
  If you provide array1 and array2 (both need to be one integer), it
  will check the data for the same sample (sample ID information is in
  experimental design) for these two arrays. If there's no common sample
  on these two arrays, the function will report an error. 
  
  In either case, you should see a nearly linear curve in all
  plots. If there were errors in hybridization and/or gridding, some of
  the plots will look messy. Then you have to check if something wrong
  happened, e.g., miss labeling, wrong gridding, etc.

  If you don't have grid information for the data, this function will be
  unavailable.

  Note that this function only works for 2-dye array.
}

\usage{
gridcheck(rawdata, array1, array2, highlight.flag = TRUE, flag.color = "Red",
          margin = c(3.1, 3.1, 3.1, 1.1)) 
}

\arguments{
  \item{rawdata}{An object of class \code{rawdata}.}
  \item{array1}{A list of array numbers for which you want to do grid
    checking. All arrays will be checked by default. If you want to
    compare the same sample across arrays, this parameter must be an
    integer to indicate the first array number.}
  \item{array2}{The second array number if you want to do cross array
    comparisons.} 
  \item{highlight.flag}{A logical parameter to indicate whether to
    highlight the bad spot or not.}
  \item{flag.color}{The color for bad spot; default is red.}
  \item{margin}{A numerical vector of the form c(bottom, left, top,
    right) which gives the lines of the margin to be specified on the four
    sides of the plot. Read \code{\link[base]{par}} for details.}
}

\note{
  This function will plot one figure for each array. So if you have many
  arrays, there will be many figures generated.
}

\examples{
\dontrun{
# load in data
data(kidney)
# grid check on the first arrays
gridcheck(kidney.raw, array1=1, margin=c(1,1,1,1))
graphics.off()
# grid check array 1 versus array 2
gridcheck(rawdata, array1=1, array2=2)
graphics.off()
}
}

\keyword{hplot}


