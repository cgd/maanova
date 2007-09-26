################################################
#
# scrpit for testing CAMDA (Critical Assessment 
# of Microarray Data Analysis) kidney data
# The website for CAMDA is:
# http://www.camda.duke.edu/
# 
############################################

# clear the work space
rm(list=ls())
# close all figures
graphics.off()

#load in library
library(maanova)

# read in data
# kidney <- read.madata(datafile="kidney.txt",designfile="kidneydesign.txt",
#   metarow=1, metacol=2, col=3, row=4, probeid=5, intensity=7,spotflag=TRUE)

data(kidney)

############################################
# STEP I: Data Quality Check on Raw Data
############################################
# grid check
gridcheck(kidney.raw)
# close figures
graphics.off()

# RI plot on raw data
riplot(kidney.raw)
# close figures
graphics.off()

# Array View on raw data
arrayview(kidney.raw,zlim=c(-7,7))
# close figures
graphics.off()


##############################################
# STEP II: Fixed model ANOVA, permutation test
# and clustering analysis
##############################################
# smooth the data using intensity lowess
kidney <- transform.madata(kidney.raw, method="glowess")
graphics.off()

# fit ANOVA model and do residual plot
anova.fix <- fitmaanova(kidney, formula= ~Dye+Array+Sample)
resiplot(kidney, anova.fix)

test.fix <- matest(kidney, anova.fix, term="Sample", n.perm=100, shuffle="sample")

# calculate FDR adjusted F values
test.fix <- adjPval(test.fix, 'jsFDR')

# volcano plot
idx.fix <- volcano(test.fix, title="Volcano plot for fixed effect model")

# clustering and consensus tree
# kmean the genes
cluster.kmean <- macluster(anova.fix,term="Sample",
	idx.gene=idx.fix$idx.all,what="gene", method="kmean", kmean.ngroups=5,
	n.perm=100)
con.kmean <- consensus(cluster.kmean, 0.7)

# HC the samples
cluster.hc <- macluster(anova.fix,term="Sample", 
         idx.gene=idx.fix$idx.all,what="sample", method="hc",
         n.perm=100)
con.hc <- consensus(cluster.hc)
