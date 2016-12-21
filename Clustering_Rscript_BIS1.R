#!/usr/bin/env R
#setwd("//SCOMP0856/veen119$/My Documents/Minor 2016/BIF 30806 Advanced Bioinformatics/Project/ballgown")
args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  infile = "gene_expression_table.csv"
  g.o.i = "CRO_003206"
  method = "complete" # can also change method to average.
  parameters = c(6,1500)
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

bis1 = read.delim(file = "./ballgown/bg_bis1/t_data.ctab", header = TRUE)
bis2 = read.delim(file = "./ballgown/bg_bis2/t_data.ctab", header = TRUE)
bis3 = read.delim(file = "./ballgown/bg_bis3/t_data.ctab", header = TRUE)
control1 = read.delim(file = "./ballgown/bg_control1/t_data.ctab", header = TRUE)
control2 = read.delim(file = "./ballgown/bg_control2/t_data.ctab", header = TRUE)
control3 = read.delim(file = "./ballgown/bg_control3/t_data.ctab", header = TRUE)

genenames = bis1$gene_name[bis1$gene_name != "."]
bis1FPKM = bis1$FPKM[bis1$gene_name != "."]
bis2FPKM = bis2$FPKM[bis2$gene_name != "."]
bis3FPKM = bis3$FPKM[bis3$gene_name != "."]
control1FPKM = control1$FPKM[control1$gene_name != "."]
control2FPKM = control2$FPKM[control2$gene_name != "."]
control3FPKM = control3$FPKM[control3$gene_name != "."]

exp_data = data.frame(bis1FPKM, bis2FPKM, bis3FPKM, control1FPKM, control2FPKM, control3FPKM, row.names = genenames)
sds = apply(exp_data, 1, sd)
idxs = which(sds == 0)
exp_data = exp_data[-idxs, ]
dim(exp_data)
head(exp_data)

write.table(exp_data, file=infile, sep=",", row.names=TRUE)

test_data=read.table(
  infile, 
  row.names=1, header=TRUE, sep =",")

# check dimensions before normalization
dim(test_data)

# removes last column and calculates the sd of samples in rows
sds <- apply(test_data[,-7], 1, sd)

###___normalization method 1, cutt-off sd___###

# looking at the plot for further normalization
#plot(sds)
#plot(sds[sds<50], main = "plot of stdev below 500", xlab = "samples", ylab = "stdev")
#abline(8, 0, col = 'red', lty = 1)
#minsds <- 8
#length(which(sds>minsds))
#subdata <- test_data[which(sds>minsds),]
#dim(subdata)

###___normalization method 2, highest 2000 sd's___###

# orders the samples and returns indices of first the sample with higest sd
ind <- order(sds,decreasing=TRUE)
#takes 2000 samples
subnumber <- parameters[2]
#takes only first 2000 with highest sd, will differ most in gene expression
subdata <- test_data[ind[1:subnumber],]
# check dimensions after normalization
dim(subdata)

# transpose the data and rename columns
tsubdata <- data.frame(t(subdata[,-7]))

# Check whether gene(s) of interest is still contained in dataset.
which(rownames(subdata) == g.o.i)

#set the number of clusters we want
numberclusters <- parameters[1]

# calculate the correlation distances
tsubdata_eucl_cor <- as.dist(1-cor(tsubdata))
# Do a hierarchical correlation clustering
tsubdata_cluster_cor = hclust(tsubdata_eucl_cor, method=method)
# split the data into 6 clusters
tsubdata_tree_cluster_cor = cutree(tsubdata_cluster_cor, k=numberclusters)

table(tsubdata_tree_cluster_cor)

# a function for drawing the correlation cluster graph
draw_graph_cluster <- function(clustnumber, colornumber) {
  n=which(tsubdata_tree_cluster_cor==clustnumber)[1]
  plot(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l",col=colornumber, 
       ylab="gene expression", xlab="samples", main=paste("Cluster ",clustnumber, ", gene expression profile"))
  for (j in 2:length(which(tsubdata_tree_cluster_cor==clustnumber))) {
    n=which(tsubdata_tree_cluster_cor==clustnumber)[j]
    points(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l",col=colornumber)
  }
}

# still need to write the plots to a file

# for every cluster draw the hierarchical correlation graph
#pdf(file="cluster_graphs.pdf")
#par(mfrow = c(numberclusters/3, 3))
#for (i in 1:numberclusters) {
#  draw_graph_cluster(i, i)
#}
#dev.off()
#par(mfrow=c(1,1))
# Zoom in on the graph with the gene of interest

target <- tsubdata_tree_cluster_cor[g.o.i]
names(target)<-'Cluster containing the gene of interest'
cl.sel <- which(tsubdata_tree_cluster_cor == target)
#cl.sel      # this will show you all the genes in the gene cluster you selected
cl.sel.df <- as.data.frame(cl.sel)

# write all the interesting genes, that are contained inside the cluster, to a file
write.table(cl.sel.df, file="Interesting_genes.csv", sep=",", row.names=TRUE)

target      # this will show you which cluster you have extracted
#pdf(file="interesting_genes_graph.pdf")
#draw_graph_cluster(target,target)
#dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Some experimenting with EdgeR
"
#### Do this only once per PC! it takes some time
#install the edgeR  library
"
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
"
##########

#load the library this has to be done each time. 
library(edgeR)

length(tsubdata)
DGlist <- DGEList(counts=tsubdata[,1:length(tsubdata)])
##checking that the lib.size contains the sum of the counts for each sample
# These numbers are the same
DGlist$samples$lib.size
colSums(DGlist$counts)

##TMM (Trimmed mean of M-values) normalization 
# this normalization method finds a set of scaling factors for
# the libary sizes that minimize the log-fold changes between the
# samples for most genes
# http://evowiki.haifa.ac.il/index.php?title=The_trimmed_mean_of_M-values_normalization_method_(TMM)
DGlist.norm <- calcNormFactors(DGlist)
DGlist.norm$samples
# lib.size is the number of RNAseq counts
DGlist.norm$samples$lib.size
DGlist.norm$samples$norm.factors

##Data exploration
# multidim scaling plot of distances between gene expression profiles
# distances in the plot approximate the log2 fold changes between the sapmles
plotMDS(DGlist.norm)

##similar plot plot with distances defined in terms of shrunk fold
logCPM <- predFC(DGlist.norm, prior.count=2*ncol(DGlist))
"

#plotMDS(logCPM, main="logFC distance")


