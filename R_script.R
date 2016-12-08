#!/usr/bin/env R
#setwd("//scomp0854/bosch098$/My Documents/Adv. Bioinformatics/project")
args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  infile = "http://www.bioinformatics.nl/courses/BIF-30806/maize_e3.table"
  g.o.i = "GRMZM2G054123"
  method = "complete" # can also change method to average.
  parameters = c(6,2000)
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


test_data=read.table(
  infile, 
  row.names=1, header=TRUE, sep ="\t")


#takes 2000 samples
subnumber <- parameters[2]
# removes last column and calculates the sd of samples in rows
sds <- apply(test_data[,-7], 1, sd)
# orders the samples and returns indices of first the sample with higest sd
ind <- order(sds,decreasing=TRUE)

#takes only first 2000 with highest sd, will differ most in gene expression
subdata <- test_data[ind[1:subnumber],]

tsubdata <- data.frame(t(subdata[,-7]))
rownames(tsubdata) = c("first1", "first2", "first3", "second1", "second2", "second3")

numberclusters <- parameters[1]

tsubdata_eucl_cor <- as.dist(1-cor(tsubdata))
# 
tsubdata_cluster_cor = hclust(tsubdata_eucl_cor, method=method)
# split the data into 6 clusters
tsubdata_tree_cluster_cor = cutree(tsubdata_cluster_cor, k=numberclusters)

table(tsubdata_tree_cluster_cor)

draw_graph_cluster <- function(clustnumber, colornumber) {
  n=which(tsubdata_tree_cluster_cor==clustnumber)[1]
  plot(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l",col=colornumber, 
       ylab="gene expression", xlab="samples", main="Clustering gene expression profiles")
  for (j in 2:length(which(tsubdata_tree_cluster_cor==clustnumber))) {
    n=which(tsubdata_tree_cluster_cor==clustnumber)[j]
    points(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l",col=colornumber)
  }
}

# still need to write the plots to a file
par(mfrow = c(2, 3))
for (i in 1:numberclusters) {
  draw_graph_cluster(i, i)
}
par(mfrow=c(1,1))

target <- tsubdata_tree_cluster_cor[g.o.i]
names(target)<-'Cluster containing the gene of interest'
cl.sel <- which(tsubdata_tree_cluster_cor == target)
cl.sel      # this will show you all the genes in the gene cluster you selected
cl.sel.df <- as.data.frame(cl.sel)

# write all the interesting genes to a file
write.table(cl.sel.df, file="Interesting_genes.csv", sep=",", row.names=TRUE)

target      # this will show you which cluster you have extracted

draw_graph_cluster(target,target)
