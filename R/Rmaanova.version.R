#####################################################################
#
# Rmaanova.version.R
#
# copyright (c) 2001-2004, Lei Wu and Gary A. Churchill, The Jackson Lab.
# Written July, 2006
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
# This is the function to display the version of the current package
#
######################################################################

Rmaanova.version <- function() {
    cat("\nCurrent R/maanova version is ")
    u <- strsplit(library(help = maanova)[[3]][[1]][2], " ")[[1]]
    cat(u[length(u)])
    cat("\n\n")
}

    
