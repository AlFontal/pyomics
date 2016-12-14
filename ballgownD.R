# Installation and loading of the required packages for the analysis

lib_loc = "/home/kersj001/R_libs"

library(ballgown, lib.loc = lib_loc)
library(genefilter, lib.loc = lib_loc)
library(dplyr, lib.loc = lib_loc)
#library(devtools, lib.loc = lib_loc)
#devtools::install_github('alyssafrazee/RSkittleBrewer')
#library(RSkittleBrewer, lib.loc = lib_loc)

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

head(results_transcripts)
head(results_genes)

# Write results to csv file
write.csv(results_transcripts, "ballgown/bis_transcript_results.csv", row.names = FALSE)
write.csv(results_transcripts, "ballgown/bis_gene_results.csv", row.names = FALSE)

# Identify transcripts and genes with a pvalue < 0.05 and a fold change > 2
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05 & results_transcripts$fc>2)
sig_genes = subset(results_genes, results_genes$pval<0.05 & results_transcripts$fc>2)

sig_transcripts
sig_genes

write.csv(sig_transcripts, "ballgown/sig_transcripts.csv")
write.csv(sig_genes, "ballgown/sig_genes.csv")

# Show the distribution of gene abundances (measured as FPKM values) across samples
fpkm = texpr(bg_bis, meas="FPKM")

# Transform the FPKM data with a log2 transformation to visualize the plots more easily
fpkm = log2(fpkm+1)

# Creating boxplot
boxplot(fpkm, col=as.numeric(pheno_data$samples), las=2, ylab='log2(FPKM+1)')

# Show the name of the transcript and the name of the gene
ballgown::transcriptNames(bg_bis)[1]
ballgown::geneNames(bg_bis)[1] 
  
# Creating a plot for the 1st transcript in the data set
plot(fpkm[1,] ~ pheno_data$samples, border=c(1,2), 
		main=paste(ballgown::geneNames(bg_bis)[1], ' : ',
		ballgown::transcriptNames(bg_bis)[1]), pch=19, 
		xlab="Samples", ylab='log2(FPKM+1)')
points(fpkm[1,] ~ jitter(as.numeric(pheno_data$samples)), 
			col=as.numeric(pheno_data$samples))
