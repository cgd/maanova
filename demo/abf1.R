 ################################################################
#
# This is the script to demonstrate the data analysis on a
# multiple factor experiment using affymetric arrays.
#
################################################################

# clear the work space
rm(list=ls())
# close all figures
graphics.off()

# load library and read in data
library(maanova)

# read in data
# library(affy)
# beforeRma <- ReadAffy()
# rmaData <- rma(beforeRma)
# datafile <- exprs(rmaData)
# design.table <- data.frame(Array=row.names(pData(beforeRma)));
# Strain <- rep(c('Aj', 'B6', 'B6xAJ'), each=6)
# Sample <- rep(c(1:9), each=2)
# designfile <- cbind(design.table, Strain, Sample)
# abf1 <- read.madata(datafile, designfile=designfile)
#####################################################
# another way to read data
# abf1 <- read.madata("AffyData.txt", designfile="AffyDesign.txt", 
#	probeid=1, intensity=2)

data(abf1)
##############################
# STEP I: fixed model analysis - skipped
##############################
# fit ANOVA model
# full model

anova.full.fix <- fitmaanova(madata=abf1, formula=~Strain+Sample)

# F-test on Sample
test.sample.fix <- matest(data=abf1, anovaobj=anova.full.fix, term="Sample", n.perm=500,
        shuffle.method="sample", test.method=c(1,1))
idx.fix <- volcano(test.sample.fix, title="Volcano Plot - fixed model", threshold=rep(0.05,2))

##############################
# STEP II: mixed model analysis
##############################
# fit mixed effect model object, treat Sample as random factor
anova.full.mix <- fitmaanova(madata = abf1,  formula=~Strain+Sample, random=~Sample)
varplot(anova.full.mix)

# F-test on Strain
test.Strain.mix <- matest(abf1, anova.full.mix, term="Strain", n.perm=nperm,
        test.method=c(1,1))
# FDR adjustment
library(qvalue)
test.Strain.mix <- adjPval(test.Strain.mix, method = 'jsFDR')
summarytable(test.Strain.mix, outfile='all.csv')

# volcano plot
x11();idx.mix <- volcano(test.Strain.mix, 
	title="Volcano Plot - Testing Strain - mixed model",
	threshold=rep(0.05,2))

# T-test on Strain
C <- matrix(c(1,-1,0, 1,0,-1, 0,1,-1),3, byrow=T)
ttest.strain.mix <- matest(abf1, anova.full.mix, term="Strain", Contrast=C, n.p=nperm, 
	test.method=rep(1,2))
#save(ttest.strain.mix, file="ttestmix.RData")
#load("ttestmix.RData")
volcano(ttest.strain.mix, threshold=rep(0.05,2))

