#####################################################################
#
# adjPval.R
#
# copyright (c) 2001-2004, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written Apr, 2004
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
# This is the function to adjust F test object
#
######################################################################
adjPval <- function(matestobj, method=c("stepup","adaptive", "stepdown", 
  "jsFDR"))
{
  if(class(matestobj)[1] != "matest")
    stop("First input object is not an object of class matest")

  fdr.method <- function(x) {fdr(x, method)}
  # start to adjust
  # for F1
  if( !(is.null(matestobj$F1)) ) {
    matestobj$F1$adjPtab <- apply(matestobj$F1$Ptab, 2, fdr.method)
    if( !(is.null(matestobj$F1$Pvalperm)) )
      matestobj$F1$adjPvalperm <- apply(matestobj$F1$Pvalperm, 2, fdr.method)
  }
  
  # for Fs
  if( !(is.null(matestobj$Fs)) ) {
    matestobj$Fs$adjPtab <- apply(matestobj$Fs$Ptab, 2, fdr.method)
    if( !(is.null(matestobj$Fs$Pvalperm)) )
      matestobj$Fs$adjPvalperm <- apply(matestobj$Fs$Pvalperm, 2, fdr.method)
  }
  
  # for Fss
  #if( !(is.null(matestobj$Fss)) ) {
  #  matestobj$Fss$adjPtab <- apply(matestobj$Fss$Ptab, 2, fdr.method)
  #  if( !(is.null(matestobj$Fss$Pvalperm)) )
  #    matestobj$Fss$adjPvalperm <- apply(matestobj$Fss$Pvalperm, 2, fdr.method)
  #}
  # return
  invisible(matestobj)
}
