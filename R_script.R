#!/usr/bin/env R
#setwd("//SCOMP0856/veen119$/My Documents/Minor 2016/BIF 30806 Advanced Bioinformatics/Project")
args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  infile = "http://www.bioinformatics.nl/courses/BIF-30806/maize_e3.table"
  g.o.i = "GRMZM2G054123"
  method = "complete"
  parameters = c(6,2000)
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#library(DESeq2)
#source("http://www.bioinformatics.nl/courses/BIF-30806/DEseq2Exercise.R")

test_data=read.table(
  infile, 
  row.names=1, header=TRUE, sep ="\t")

subnumber <- parameters[2]
sds <- apply(test_data[,-7], 1, sd)
ind <- order(sds,decreasing=TRUE)
subdata <- test_data[ind[1:subnumber],]

tsubdata <- data.frame(t(subdata[,-7]))
rownames(tsubdata) = c("first1", "first2", "first3", "second1", "second2", "second3")

numberclusters <- parameters[1]

tsubdata_eucl_cor <- as.dist(1-cor(tsubdata))
tsubdata_cluster_cor = hclust(tsubdata_eucl_cor, method=method)
tsubdata_tree_cluster_cor = cutree(tsubdata_cluster_cor, k=numberclusters)

table(tsubdata_tree_cluster_cor)

draw_graph_cluster <- function(clustnumber, colornumber) {
  n=which(tsubdata_tree_cluster_cor==clustnumber)[1]
  plot(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l",col=colornumber)
  for (j in 2:length(which(tsubdata_tree_cluster_cor==clustnumber))) {
    n=which(tsubdata_tree_cluster_cor==clustnumber)[j]
    points(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l",col=colornumber)
  }
}

par(mfrow = c(2, 3))
for (i in 1:numberclusters) {
  draw_graph_cluster(i, i)
}
par(mfrow=c(1,1))

target <- tsubdata_tree_cluster_cor[g.o.i]
names(target)<-'Cluster containing the gene of interest'
cl.sel <- which(tsubdata_tree_cluster_cor == target)
cl.sel      # this will show you all the genes in the gene cluster you selected
target      # this will show you which cluster you have extracted

draw_graph_cluster(target,target)
