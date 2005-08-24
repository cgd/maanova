################################################################
#
# This is the script to demonstrate the data analysis on a
# multiple factor experiment.
#
# The expeirment was done in Bev Paigen's lab in The Jackson Lab
# It's a big 28 array experiment. I hand pick 300 genes from it
# to speed up the calculation. The data is spatial lowessed
# so we don't need to do data normalization here
#
################################################################

# clear the work space
rm(list=ls())
# close all figures
graphics.off()

# load library and read in data
library(maanova)

# read data
#paigen.raw <- read.madata("paigen.txt", designfile="design.txt",cloneid=1, 
#        metarow=2, metacol=3, row=4, col=5, pmt=6, spotflag=F)

# load in data
data(paigen)
# because this is not a completed data so data quality check
# is skipped

# make data object with rep 2
paigen <- createData(paigen.raw, 2)

# Note that the data is normalized so normalization is skipped

##############################
# STEP I: fixed model analysis
##############################
# make several model objects and fit ANOVA model
# full model
model.full.fix <- makeModel(data=paigen, formula=~Array+Dye+Spot+Strain+Diet+Strain:Diet+Sample)
anova.full.fix <- fitmaanova(paigen, model.full.fix)
# model without Sample
model.noSample.fix <- makeModel(data=paigen, formula=~Array+Dye+Spot+Strain+Diet+Strain:Diet)
anova.noSample.fix <- fitmaanova(paigen, model.noSample.fix)
# model without interaction
model.noint.fix <- makeModel(data=paigen, formula=~Array+Dye+Spot+Strain+Diet)
anova.noint.fix <- fitmaanova(paigen, model.noint.fix)
# model without strain effect
model.nostrain.fix <- makeModel(data=paigen, formula=~Array+Dye+Spot+Diet)
anova.nostrain.fix <- fitmaanova(paigen, model.nostrain.fix)
# model without diet effect
model.nodiet.fix <- makeModel(data=paigen, formul=~Array+Dye+Spot+Strain)
anova.nodiet.fix <- fitmaanova(paigen, model.nodiet.fix)

# permutation tests - all test use non-restricted residual shuffling
# test for Sample
test.Sample.fix <- matest(paigen, model.full.fix, term="Sample", 
	n.perm=500, shuffle.method="resid",  test.method=rep(1,4))
#save(test.Sample.fix, file="testsamplefix.RData")
#load("testsamplefix.RData")

# test for interaction effect
test.int.fix <- matest(paigen, model.noSample.fix, term="Strain:Diet", 
	n.perm=500, shuffle.method="resid",  test.method=rep(1,4))
#save(test.int.fix, file="testintfix.RData")
#load("testintfix.RData")
 
# test for strain effect
test.strain.fix <- matest(paigen, model.noint.fix, term="Strain", n.perm=500,
        shuffle.method="resid", test.method=rep(1,4))
#save(test.strain.fix, file="teststrainfix.RData")
#load("teststrainfix.RData")

# test for diet effect
test.diet.fix <- matest(paigen, model.noint.fix, term="Diet", n.perm=500,
        shuffle.method="resid", test.method=rep(1,4))
#save(test.diet.fix, file="testdietfix.RData")
#load("testdietfix.RData")

# calculate the FDR adjusted P values
test.Sample.fix <- adjPval(test.Sample.fix)
test.int.fix <- adjPval(test.int.fix)
test.strain.fix <- adjPval(test.strain.fix)
test.diet.fix <- adjPval(test.diet.fix)

# volcano plot
idx.sample.fix <- volcano(test.Sample.fix, title="Test for Sample - fixed model")
x11();idx.int.fix <- volcano(test.int.fix, title="Interaction test - fixed model")
x11(); idx.strain.fix <- volcano(test.strain.fix, title="Strain test - fixed model")
x11(); idx.diet.fix <- volcano(test.diet.fix, title="Diet test - fixed model")

# T-test pairwise comparison of Strain
C <- matrix(c(1,-1,0,1,0,-1, 0,1,-1), nrow=3, byrow=T)
ttest.strain.fix <- matest(paigen, model.noint.fix, term="Strain", Contrast=C, 
	n.perm=500, test.method=rep(1,4))
#save(ttest.strain.fix, file="tteststrainfix.RData")
#load("tteststrainfix.RData")
# volcano plot
volcano(ttest.strain.fix)

#################################
# STEP II: Mixed model analysis
#################################
# make severl model objects and fit ANOVA
# full model
model.full.mix <- makeModel(data=paigen, 
	formula=~Array+Spot+Dye+Strain+Diet+Strain:Diet+Sample,
	random=~Array+Spot+Sample)
# fit anova
anova.full.mix <- fitmaanova(paigen, model.full.mix, method="REML")
save(anova.full.mix, file="anovafullmix.RData")
#load("anovafullmix.RData")
varplot(anova.full.mix)

# test interaction effect
ftest.int.mix <- matest(data=paigen, model=model.full.mix, term="Strain:Diet",n.perm=100,
	test.method=rep(1,4))
save(ftest.int.mix, file="ftestintmix.RData")
#load("ftestintmix.RData")

# volcano plot
x11(); idx.int.mix <- volcano(ftest.int.mix, title="Interaction test - mixed model")

# model without interaction
model.noint.mix <- makeModel(data=paigen, formula=~Dye+Array+Spot+Strain+Diet+Sample, 
        random=~Array+Spot+Sample)
anova.noint.mix <- fitmaanova(paigen, model.noint.mix, method="REML")
save(anova.noint.mix, file="anovanointmix.RData")
#load("anovanointmix.RData")
varplot(anova.noint.mix)

# use model without interaction to test Strain and Diet effect
ftest.strain.mix <- matest(data=paigen, model=model.noint.mix, term="Strain", 
	n.perm=100, test.method=rep(1,4))
save(ftest.strain.mix, file="fteststrainmix.RData")
#load("fteststrainmix.RData")

ftest.diet.mix <- matest(data=paigen, model=model.noint.mix, term="Diet", 
	n.perm=100, test.method=rep(1,4))
save(ftest.diet.mix, file="ftestdietmix.RData")
#load("ftestdietmix.RData")

# volcano plot
x11();idx.strain.mix <- volcano(ftest.strain.mix, title="Strain test - mixed model")
x11();idx.diet.mix <- volcano(ftest.diet.mix, title="Diet test - mixed model")

