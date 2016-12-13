# Installation and loading of the required packages for the analysis


library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
devtools::install_github('alyssafrazee/RSkittleBrewer')
library(RSkittleBrewer)

# Set working directory to the main group folder.
setwd("/local/data/BIF30806_2016_2/project/groups/Pyomics")

# Specify phenodata.csv as the phenodata file to be used by ballgown
pheno_data = read.csv("ballgown/phenodata.csv")

# Read in the expression data that was calculated by StringTie
bg_bis = ballgown(dataDir = "ballgown", samplePattern = "bg_", pData = pheno_data)

# Filter out those genes whose count variance among all samples is less than 1.
bg_filtered = subset(bg_bis, "rowVars(texpr(bg_bis)) >1", genomesubset=TRUE)

# Make a statistical test for all the transcripts to show statistically significant differences between groups
results_transcripts = stattest(bg_filtered, feature="transcript", covariate = "group", getFC=TRUE, meas= "FPKM")

# Make a statistical test for all the genes to show statistically significant differences between groups
results_genes = stattest(bg_filtered, feature="gene", covariate = "group", getFC=TRUE, meas= "FPKM")

# Add gene names and gene IDs to the results_genes dataframe:

results_transcripts = data.frame(geneNames = ballgown::geneNames(bg_filtered),
                                 geneIDs = ballgown::geneIDs(bg_filtered),
                                 results_transcripts)

# Sort both results dataframes from smallest to largest p-value:
results_transcripts = arrange(results_transcripts, pval)
results_genes = arrange(results_genes, pval)

# Write results to csv file
write.csv(results_transcripts, "ballgown/bis_transcript_results.csv", row.names = FALSE)
write.csv(results_transcripts, "ballgown/bis_gene_results.csv", row.names = FALSE)

# Identify transcripts and genes with a qvalue < 0.05
sig_transcripts = subset(results_transcripts, results_transcripts$qval<0.05)
sig_genes = subset(results_genes, results_genes$qval<0.05)

write.csv(sig_transcripts, "ballgown/sig_transcripts.csv")
write.csv(sig_genes, "ballgown/sig_genes.csv")



