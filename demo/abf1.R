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
#abf1.raw <- read.madata("AffyData.txt", designfile="AffyDesign.txt", 
#	probeID=1, pmt=2, spotflag=F)

# load in data
data(abf1)
# make data object with rep 1
# note that this data was pre-transformed so log.trans=F
abf1 <- createData(abf1.raw, 1, log.trans=F)

##############################
# STEP I: fixed model analysis - skipped
##############################
# make several model objects and fit ANOVA model
# full model
model.full.fix <- makeModel(data=abf1, formula=~Strain+Sample)
anova.full.fix <- fitmaanova(abf1, model.full.fix)

# F-test on Sample
test.sample.fix <- matest(abf1, model.full.fix, term="Sample", n.perm=500,
        shuffle.method="resid", test.method=rep(1,4))
idx.fix <- volcano(test.sample.fix, title="Volcano Plot - fixed model")

##############################
# STEP II: mixed model analysis
##############################
# make mixed effect model object, treat Sample as random factor
model.full.mix <- makeModel(data=abf1, formula=~Strain+Sample, random=~Sample)
anova.full.mix <- fitmaanova(abf1, model.full.mix)
varplot(anova.full.mix)

# F-test on Strain
test.Strain.mix <- matest(abf1, model.full.mix, term="Strain", n.perm=500,
        test.method=rep(1,4))
# FDR adjustment
test.Strain.mix <- adjPval(test.Strain.mix)
#save(test.Strain.mix, file="testStrainMix.RData")
#load("testStrainMix.RData")
# volcano plot
x11();idx.mix <- volcano(test.Strain.mix, 
	title="Volcano Plot - Testing Strain - mixed model",
	threshold=rep(0.05,4))


