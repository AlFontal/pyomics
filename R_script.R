source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
#source("http://www.bioinformatics.nl/courses/BIF-30806/DEseq2Exercise.R")

test_data=read.table(
  "http://www.bioinformatics.nl/courses/BIF-30806/maize_e3.table", 
  row.names=1, header=TRUE, sep ="\t")

n <- 2000
sds <- apply(test_data[,-7], 1, sd)
ind <- order(sds,decreasing=TRUE)
subdata <- test_data[ind[1:n],]
#subdata$tissue <- test_data$tissue

tsubdata <- data.frame(t(subdata[,-7]))
rownames(tsubdata) = c("first1", "first2", "first3", "second1", "second2", "second3")
dim(cor(tsubdata))

tsubdata_eucl_cor <- as.dist(1-cor(tsubdata))
tsubdata_complete_cor = hclust(tsubdata_eucl_cor, method="complete")
tsubdata_average_cor = hclust(tsubdata_eucl_cor, method="average")
tsubdata_tree_complete_cor = cutree(tsubdata_complete_cor, k=6)
tsubdata_tree_average_cor = cutree(tsubdata_average_cor, k=6)

table(tsubdata_tree_complete_cor)
table(tsubdata_tree_average_cor)

plot(tsubdata_complete_cor)
plot(tsubdata_average_cor)

par(mfrow = c(2, 3))

draw_graph_average <- function(clustnumber, colornumber) {
  n=which(tsubdata_tree_average_cor==clustnumber)[1]
  plot(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l",col=colornumber)
  for (j in 2:length(which(tsubdata_tree_average_cor==clustnumber))) {
    n=which(tsubdata_tree_average_cor==clustnumber)[j]
    points(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l",col=colornumber)
  }
}

draw_graph_complete <- function(clustnumber, colornumber) {
  n=which(tsubdata_tree_complete_cor==clustnumber)[1]
  plot(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l",col=colornumber)
  for (j in 2:length(which(tsubdata_tree_complete_cor==clustnumber))) {
    n=which(tsubdata_tree_complete_cor==clustnumber)[j]
    points(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l",col=colornumber)
  }
}

#for (i in 1:6) {
#  draw_graph_average(i, i)
#}
for (i in 1:6) {
  draw_graph_complete(i, i)
}

par(mfrow=c(1,1))

target <- tsubdata_tree_complete_cor["GRMZM2G054123"]
names(target)<-'Cluster containing GRMZM2G054123'
cl.sel <- which(tsubdata_tree_complete_cor == target)
cl.sel      # this will show you all the yeast genes in the gene cluster you selected
target      # this will show you which cluster you have extracted

n = which(tsubdata_tree_complete_cor == target)[1]
plot(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n], type="l")
for (i in 2:length(which(tsubdata_tree_complete_cor==target))) {
  n=which(tsubdata_tree_complete_cor==target)[i]
  points(as.matrix(scale(tsubdata, center=TRUE, scale=TRUE))[,n],type="l")
}


#source("http://bioconductor.org/biocLite.R")
#biocLite(c("yeastCC","GOstats","Category","org.Sc.sgd.db","Mfuzz","biomaRt","r
#GADEM", "rtracklayer", "seqLogo"))
#biocLite(c("Rgraphviz"))

#library(yeastCC)  # expression data of yeast cell cycle
#library(GOstats)  # GO testing tool package
#library(Category)
#library(org.Sc.sgd.db) # yeast gene annotation package
#library(Mfuzz)    # clustering tool package
#library(biomaRt)  # query Ensembl genome database via Biomart

