\author{Hao Wu \email{hao@jax.org} }

\name{volcano}
\alias{volcano}

\title{Volcano plot for F test results}

\description{
  This function generates a volcano-like plot given the F test results. 
}

\details{
  This function allows one to visualize the results from the F or T 
  tests. The figure looks like an erupting volcano. There will be one
  plot For F-test result and multiple plots for T-test result, each plot
  crresponds to one T-test. You must have F1 test result in the input
  object in order to do volcano plot. 

  On the plot, the y-axis value is -log10(P-value) for the F1 test.
  The x-axis value is propotional to the fold changes. 
  A horizontal line represents the significance threshold of
  the F1 test. The red dots are the genes selected by the F2
  test (if there's F2 test result). 
  The green dots are the genes selected by the F3
  test (if there's F3 test result).
  The orange dots are the genes
  selected by the Fs test (if there's Fs test result). If there
  is flag information in the data and the user wants to highlight the
  flagged genes, the genes with any bad spots will be circled by a black
  circle.  
}
  
\usage{
volcano(matestobj, threshold=c(0.001,0.05,0.05,0.05),
        method=c("unadj","unadj","unadj","unadj"), title="Volcano Plot",
        highlight.flag=TRUE, onScreen=TRUE)
}

\arguments{
  \item{matestobj}{An object of class \code{matest}.}
  \item{threshold}{A vector of four double values to indicate the
    thresholds for four F tests. The values should be between 0 and 1.
    Note that you need to put four values here even if you don't have
    all four F tests in matestobj.}
  \item{method}{A flag indicating to use which P values to generate the
    plot and select genes. This is a vector with four elements, which
    corresponds to four F tests. Each element should be one of the
    following five selections:
    \itemize{
      \item{"unadj"}{Unadjusted tabulated P values.}
      \item{"nominal"}{Nominal permutation P values.}
      \item{"fwer"}{FWER one-step adjusted P values.}
      \item{"fdr"}{FDR adjusted tabulated P values.}
      \item{"fdrperm"}{FDR adjusted nominal permutation P values.}
    }
    Default value is c("unadj", "unadj", "unadj", "unadj") which means to use
    tabulated P values for all four tests.

    Note that you need to put four values here even if you don't have
    all four F tests in matestobj.}
  \item{title}{Figure title. Default is "Volcano Plot".}
  \item{highlight.flag}{A logical value to indicate whether to highlight
    the genes with bad spots or not.}
  \item{onScreen}{A logical value to represent whether to display the
    plots on screen or not. If TRUE, the figure will be plotted on the
    screen. Otherwise, it will plot the figure on the current device.
    Default is TRUE.}
}

\value{
  For F-test volcano plot, it returns an object which is a list of the
  following four fields: 
  \item{idx.F1}{The significant genes selected by F1 test.}
  \item{idx.F2}{The significant genes selected by F2 test.}
  \item{idx.F3}{The significant genes selected by F3 test.}
  \item{idx.Fs}{The significant genes selected by Fs test.}
  \item{idx.all}{The significant genes selected by all four F tests.}

  For T-test volcano plot, it returns an array of the above object. Each
  element in the array corresponds to one T-test.
}


\examples{
data(paigen)
paigen <- createData(paigen.raw, n.rep=2)

# make model without interaction
model.noint.fix <- makeModel(data=paigen, formula=~Array+Dye+Spot+Strain+Diet)

# F-test strain effect
\dontrun{
test.strain.fix <- matest(paigen, model.noint.fix, term="Strain", n.perm=500,
        shuffle.method="resid", test.method=rep(1,4))
# volcano plot
idx.strain.fix <- volcano(test.strain.fix, title="Strain test - fixed model")

# T-test all pairwise comparison on strain
C <- matrix(c(1,-1,0,1,0,-1, 0,1,-1), nrow=3, byrow=TRUE)
ttest.strain.fix <- matest(paigen, model.noint.fix, term="Strain",
         Contrast=C, n.perm=500, test.method=rep(1,4))
# volcano plot
volcano(ttest.strain.fix)
}
}

\keyword{hplot}
