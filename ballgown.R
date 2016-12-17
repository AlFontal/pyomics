# Installation and loading of the required packages for the analysis
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
devtools::install_github('alyssafrazee/RSkittleBrewer')
library(RSkittleBrewer)
library(calibrate)

# Set plot features to fit presentation
par(bg = "#2e3037", fg = "#39c0ba", col.axis = "#39c0ba", col.main = "#39c0ba", col.lab = "#39c0ba")
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow", "darkred")
palette(tropical)



# Set working directory to the main group folder.
setwd("/home/alejandro/Documents/MsC Bioinformatics/Advanced Bioinformatics/")

# Specify phenodata.csv as the phenodata file to be used by ballgown
pheno_data = read.csv("Ballgown/phenodata.csv")

# Read in the expression data that was calculated by StringTie
bg_bis = ballgown(dataDir = "Ballgown", samplePattern = "bg_", pData = pheno_data)
bg_bis = subset(bg_bis, 'ballgown::geneNames(bg_bis) != "."')

# Filter out those genes whose count variance among all samples is less than 1.
bg_filtered = subset(bg_bis, "rowVars(texpr(bg_bis)) >1", genomesubset=TRUE)

# Make a statistical test for all the transcripts to show statistically significant differences between groups
results_transcripts = stattest(bg_filtered, feature="transcript", covariate = "group", getFC=TRUE, meas= "FPKM")

# Add gene names and gene IDs to the results dataframe:

results_transcripts = data.frame(geneNames = ballgown::geneNames(bg_filtered),
                                 geneIDs = ballgown::geneIDs(bg_filtered),
                                 results_transcripts)

# Sort results dataframes from smallest to largest p-value:
results_transcripts = arrange(results_transcripts, pval)

# Add log fold change column (Will be useful to subset later)
results_transcripts$logfc = log2(results_transcripts$fc)



# Identify transcripts and genes with a p-value < 0.05 and a logFC bigger than 1 or smaller than -1
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05 & abs(results_transcripts$logfc) > 1)

# Make histogram of distribution of Log2 Fold Change and P-values among mapped transcripts. 

hist(results_transcripts$logfc, breaks = 200, xlim = c(-3,3), col = "#2e3037", main = "Histogram of logFC", 
     xlab = "Log2 Fold-Change")

hist(results_transcripts$pval, breaks = 100, xlim = c(0, 1), col = "#2e3037", main = "Histogram of p-values", 
     xlab = "p-value")


# Read annotations and add them to the dataframe 
annots = read.delim("/home/alejandro/Documents/GitHub/pyomics/cro_functional_annotation.final.txt", header = FALSE)
colnames(annots) = c("geneNames", "geneFunction")
annots$geneNames = paste(substr(annots$geneNames, 1, 4), substr(annots$geneNames, 6, 11), sep = "")
sig_transcripts = merge(sig_transcripts, annots, by= "geneNames")
results_transcripts = merge(results_transcripts, annots, by= "geneNames")

# Write results to csv file
write.csv(results_transcripts, "ballgown/bis_transcript_results.csv", row.names = FALSE)
write.csv(sig_transcripts, "ballgown/sig_genes.csv", row.names = FALSE)

# Show the distribution of gene abundances (measured as FPKM values) across samples
fpkm = texpr(bg_bis, meas="FPKM")

# Transform the FPKM data with a log2 transformation to visualize the plots more easily
fpkm = log2(fpkm+1)

# Creating boxplot of FPKM among samples
boxplot(fpkm, col=as.numeric(pheno_data$samples), las = 2, ylab='log2(FPKM+1)',
        names = c("bis1", "bis2", "bis3","ctrl1", "ctrl2", "ctrl3"), ylim = c(0,5) )




# Creating a plot for the expression of BIS1 among samples
bis = which(geneNames(bg_bis) == "CRO_003206")
plot(fpkm[bis,] ~ pheno_data$samples, border=c(1,2), main="BIS1", pch=19, 
	    xlab="Samples", ylab='log2(FPKM+1)')


setwd("/home/alejandro/Documents/GitHub/pyomics/")
# Read overlapping gene list
overlapping = read.csv("overlapping_genes.txt", header = FALSE)
overlapping = overlapping[, "V1"]

# Make volcano plot
with(results_transcripts, plot(logfc, -log10(pval), pch=20, main="Volcano plot", xlim=c(-5,5),  col = "#617189"))

# Make genes with qval < 0.5 visible 
with(subset(results_transcripts, qval < .05), points(logfc, -log10(pval), pch=20, col="#f35b69"))

# Make genes with pval <0.5 and |log2FC| > 1 visible
with(subset(results_transcripts, pval < .05 & abs(logfc) > 1), points(logfc, -log10(pval), pch=20, col="#39c0ba"))
with(subset(results_transcripts, qval < .05), points(logfc, -log10(pval), pch=20, col="#f35b69"))

# Make overlapping genes visible
with(subset(results_transcripts, geneNames %in% as.vector(overlapping)), points(logfc, -log10(pval), pch=20, col="orange"))

# Make interesting points visible
with(subset(results_transcripts, geneNames == "CRO_003206"), points(logfc, -log10(pval), pch=20, col="#b910bc"))
with(subset(results_transcripts, geneNames == "CRO_003206"), textxy(logfc, -log10(pval), labs="BIS1", cex=.8, col = "white"))
with(subset(results_transcripts, geneNames == "CRO_005949"), points(logfc, -log10(pval), pch=20, col="#b910bc"))
with(subset(results_transcripts, geneNames == "CRO_005949"), textxy(logfc, -log10(pval), labs="SecS", cex=.8, offset = -0.8, col = "white"))
with(subset(results_transcripts, geneNames == "CRO_021502"), points(logfc, -log10(pval), pch=20, col="#b910bc"))
with(subset(results_transcripts, geneNames == "CRO_021502"), textxy(logfc, -log10(pval), labs="ProgBR", cex=.8, col = "white"))
with(subset(results_transcripts, geneNames == "CRO_033829"), points(logfc, -log10(pval), pch=20, col="#b910bc"))
with(subset(results_transcripts, geneNames == "CRO_033829"), textxy(logfc, -log10(pval), labs="StricS", cex=.8,col = "white"))
with(subset(results_transcripts, geneNames == "CRO_008172"), points(logfc, -log10(pval), pch=20, col="#00ff00"))
with(subset(results_transcripts, geneNames == "CRO_032641"), points(logfc, -log10(pval), pch=20, col="#00ff00"))
with(subset(results_transcripts, geneNames == "CRO_023569"), points(logfc, -log10(pval), pch=20, col="#00ff00"))
with(subset(results_transcripts, geneNames == "CRO_023903"), points(logfc, -log10(pval), pch=20, col="#00ff00"))
with(subset(results_transcripts, geneNames == "CRO_005949"), points(logfc, -log10(pval), pch=20, col="#00ff00"))



