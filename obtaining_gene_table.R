#!/usr/bin/env R
bis1 = read.delim(file = "t_data_bis1.ctab", header = TRUE)
bis2 = read.delim(file = "t_data_bis2.ctab", header = TRUE)
bis3 = read.delim(file = "t_data_bis3.ctab", header = TRUE)
control1 = read.delim(file = "t_data_control1.ctab", header = TRUE)
control2 = read.delim(file = "t_data_control2.ctab", header = TRUE)
control3 = read.delim(file = "t_data_control3.ctab", header = TRUE)

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

write.table(exp_data, file="gene_expression_table.csv", sep=",", row.names=TRUE)
