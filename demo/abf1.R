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

data(abf1)
##############################
# STEP I: fixed model analysis - skipped
##############################
# make several model objects and fit ANOVA model
# full model

fit.full.fix <- fitmaanova(abf1, formula = ~Strain)

# F-test on Sample
ftest.fix = matest(abf1, fit.full.fix, test.method=c(1,1,0),
    shuffle.method="sample", term="Strain", n.perm= 100)
idx.fix <- volcano(ftest.fix, title="Volcano Plot - fixed model")

##############################
# STEP II: mixed model analysis
##############################
# make mixed effect model object, treat Sample as random factor
fit.full.mix <- fitmaanova(abf1, formula = ~Strain+Sample, 
    random = ~Sample)
varplot(fit.full.mix)

# F-test on Strain
ftest.mix <- matest(abf1, fit.full.mix, term="Strain", n.perm=500)

# FDR adjustment
ftest.all = adjPval(ftest.mix, 'jsFDR')

#summary table 
summarytable(ftest.all, outfile='all.csv')

# volcano plot
x11();idx.mix <- volcano(ftest.all,
	title="Volcano Plot - Testing Strain - mixed model",
	threshold=rep(0.05,3))


