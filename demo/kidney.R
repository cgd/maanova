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
#kidney.raw <- read.madata("kidney.txt", designfile="kidneydesign.txt", 
#	metarow=1, metacol=2, col=3, row=4, probeid=6,
#	intensity=7, arrayType='twoColor',log.trans=T, spotflag=T)

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
# smooth the data using lowess
kidney <- transform.madata(kidney.raw, method="rlowess")
graphics.off()

# fit ANOVA model and do residual plot
anova.fix <- fitmaanova(kidney, formula= ~Dye+Array+Sample)
resiplot(kidney, anova.fix)

test.fix <- matest(kidney, anova.fix, term="Sample", n.perm=500, shuffle="sample")

# calculate FDR adjusted F values
library(qvalue)
test.fix <- adjPval(test.fix, 'jsFDR')

# get summarytable based on Fs q-value 0.05
summarytable(test.fix, method=c("Fold.change","Pvalperm","adjPvalperm"), test = 
c("F1","Fs"),whichTest="Fs.adjPvalperm", threshold=0.05, outfile="summarytable.csv")

# volcano plot
idx.fix <- volcano(test.fix, title="Volcano plot for fixed effect model")

# clustering and consensus tree
# kmean the genes
cluster.kmean <- macluster(anova.fix,term="Sample",
	idx.gene=idx.fix$idx.all,what="gene", method="kmean", kmean.ngroups=5,
	n.perm=100)
X11();con.kmean <- consensus(cluster.kmean, 0.7)

# to see which genes belongs to which group
con.kmean$groupname
# HC the samples
cluster.hc <- macluster(anova.fix,term="Sample", 
         idx.gene=idx.fix$idx.all,what="sample", method="hc",
         n.perm=100)
X11();con.hc <- consensus(cluster.hc)

#########################################
# STEP III: Mixed model ANOVA and F test
# Mixed model analysis is skipped for this demo
#########################################
# make model object for mixed model
anova.mix <- fitmaanova(madata=kidney, formula=~Dye+Array+Sample, random=~Array)
varplot(anova.mix)
#save(anova.mix, file="anovamix.RData")
#load("anovamix.RData")
# residual plot
resiplot(kidney, anova.mix)
# test Sample
ftest.mix <- matest(data=kidney, anovaobj=anova.mix, term="Sample", n.perm=100)
#save(ftest.mix, file="ftestmix.RData")
#load("ftestmix.RData")

# volcano plot
idx.mix <- volcano(ftest.mix)

# what's the difference for mixed and fixed model?
setdiff(idx.fix$idx.all, idx.mix$idx.all)
