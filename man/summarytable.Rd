\author{Hyuna Yang}
\name{summarytable}
\alias{summarytable}
\title{Summarize the matest result.}
\description{
  This function returns list of probe ID with p-value, FDR adjusted p-value or
  fold change selected by given threshold.
}
\details{
This function use 

}
\usage{
summarytable(matestobj, method=c("Fold.change","Pvalperm","adjPvalperm"), 
     test =c("F1","Fs"),whichTest=c("F1.Pvalperm","F1.adjPvalperm", 
         "Fs.Pvalperm","Fs.adjPvalperm"), 
         threshold, outfile="summarytable.csv")
}
\arguments{
  \item{matestobj}{An object of class \code{matest}, which is the result
    from \code{\link[maanova]{matest}}.}
  \item{method}{Results that you want to include in the output among
  Fold.change, Pvalperm and adjPvalperm. Default is possible every field.}
  \item{test}{Test that you want to include in the output among F1, Fs and
  Fss. Default is possible every test.}
  \item{whichTest}{Test result used to get the subset. If whichTest is not
  provided, save all result. }
  \item{threshold}{Threshold that you want to get the subset. If threshold is
  not provided, save all result.}
  \item{outfile}{File name that you want to save the result. Default name is 'summarytable.csv'}

}
\value{
  \item{result}{It returns result fields (among Fold.change, P-value or adjust
  P-value if you selected) from test statistics (among F1, Fs or Fss
  statistics if you selected), subsetted by given statistics ('whichtest') is
  smaller than 'threshold'.}
  \item{outfile}{.csv file saved under the working directory.}
}

\examples{
data(abf1)
\dontrun{ 
fit.full.simple = fitmaanova(abf1,formula = ~ Strain)
# F-test strain effect
ftest.all = matest(abf1, fit.full.simple, term="Strain", n.perm= 1000)
# make FDR adjusted P values
ftest.all = adjPval(ftest.all, 'jsFDR')
summarytable(ftest.all, outfile='all.csv')
smallset = summarytable(ftest.all, method ='Pvalperm',
test=c('F1','Fs'), whichTest='Fs.Pvalperm', threshold = 0.1)
}}
\keyword{utilities}
