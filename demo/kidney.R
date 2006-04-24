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
#	metarow=1, metacol=2, col=3, row=4, Name=5, ID=6,
#	pmt=7, spotflag=T)

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
# make data object with rep 1
kidney <- createData(kidney.raw);
summary(kidney)
# smooth the data using intensity lowess
kidney <- transform.madata(kidney, method="glowess")
graphics.off()

# make model object for fixed model
# note that because there's no replicates, Spot and Label cannot be fitted
model.fix <- makeModel(data=kidney, formula=~Dye+Array+Sample)
summary(model.fix)
# fit ANOVA model and do residual plot
anova.fix <- fitmaanova(kidney, model.fix)
resiplot(kidney, anova.fix)

test.fix <- matest(kidney, model=model.fix, term="Sample", n.perm=100, shuffle="sample")
# calculate FDR adjusted F values
test.fix <- adjPval(test.fix)
#save(ftest.fix, file="ftestfix.RData")
#load("ftestfix.RData")

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


#########################################
# STEP III: Mixed model ANOVA and F test
# Mixed model analysis is skipped for this demo
#########################################
# make model object for mixed model
model.mix <- makeModel(data=kidney, formula=~Dye+Array+Sample, random=~Array)
summary(model.mix)
# the following line takes long time to finish,
anova.mix <- fitmaanova(kidney, model.mix)
varplot(anova.mix)
#save(anova.mix, file="anovamix.RData")
#load("anovamix.RData")
# residual plot
resiplot(kidney, anova.mix)
# test Sample
ftest.mix <- matest(data=kidney, model=model.mix, term="Sample")
#save(ftest.mix, file="ftestmix.RData")
#load("ftestmix.RData")

# volcano plot
idx.mix <- volcano(anova.mix, ftest.mix, 
	method=c("unadj","unadj","unadj"),
	title="Volcano plot for mixed effect model")

