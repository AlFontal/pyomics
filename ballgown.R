# Installation and loading of the required packages for the analysis


library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
devtools::install_github('alyssafrazee/RSkittleBrewer')
library(RSkittleBrewer)

# Set working directory to the main group folder.
setwd("/home/alejandro/Documents/MsC Bioinformatics/Advanced Bioinformatics/")

# Specify phenodata.csv as the phenodata file to be used by ballgown
pheno_data = read.csv("Ballgown/phenodata.csv")

# Read in the expression data that was calculated by StringTie
bg_bis = ballgown(dataDir = "Ballgown", samplePattern = "bg_", pData = pheno_data)

bg_bis = bg_bis[ballgown::geneNames(bg_bis) != "."]

# Filter out those genes whose count variance among all samples is less than 1.
bg_filtered = subset(bg_bis, "rowVars(texpr(bg_bis)) >1", genomesubset=TRUE)
bg_filtered = subset(bg_filtered, "geneNames != '.'", genomesubset=TRUE)

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

# Add log fold change column (Will be useful to subset later)
results_genes$logfc = log2(results_genes$fc)
results_transcripts$logfc = log2(results_transcripts$fc)

# Write results to csv file
write.csv(results_transcripts, "ballgown/bis_transcript_results.csv", row.names = FALSE)
write.csv(results_transcripts, "ballgown/bis_gene_results.csv", row.names = FALSE)

# Identify transcripts and genes with a p-value < 0.05 and a logFC bigger than 1 or smaller than -1
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05 & abs(results_transcripts$logfc) > 1)

# Keep only those which we have a name for
de_genes_df = sig_transcripts[sig_transcripts$geneNames != ".", ]
de_genes_list = de_genes_df[1]

write.csv(de_genes_df, "Ballgown/DE_genes.csv")

tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow", "darkred")
palette(tropical)

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
    plot(fpkm[i,] ~ pheno_data$samples, border=c(1,2), 
		    main=paste(ballgown::geneNames(bg_bis)[i], ' : ',
		    ballgown::transcriptNames(bg_bis)[i]), pch=19, 
	    	xlab="Samples", ylab='log2(FPKM+1)')
    points(fpkm[i,] ~ jitter(as.numeric(pheno_data$samples)), 
			col=as.numeric(pheno_data$samples))


