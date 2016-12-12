# Installation and loading of the required packages for the analysis

'''
#Uncomment this in case the local machine doesnt contain the required packages. 

source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("ballgown", "genefilter", "dplyr", "devtools"))

'''

library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
devtools::install_github('alyssafrazee/RSkittleBrewer')
library(RSkittleBrewer)

setwd("/local/data/BIF30806_2016_2/project/groups/Pyomics")
pheno_data = read.csv("ballgown/phenodata.csv")
bg_bis = ballgown(dataDir = "ballgown", samplePattern = "bg_", pData = pheno_data)
bg_filtered = bg_filtered = subset(bg_obj, "rowVars(texpr(bg_obj)) >1", genomesubset=TRUE)
results_transcripts = stattest(bg_filtered, feature="gene", covariate = "group", getFC=TRUE, meas= "FPKM")